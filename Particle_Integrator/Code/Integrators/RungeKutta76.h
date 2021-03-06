/* 7th and 6th order Runge-Kutta-Nystrom method as described in:
   J.R. Dormand & P.J. Prince (1978): New Runge-Kutta algorithms for numerical simulation in dynamical astronomy. Celestial Mechanics, Vol. 18, p. 223-232. 
   
   Steps exceeding the maximum allowed error (e_target) will be repeated. */

int RungeKutta76(configuration_values *config_data, SpiceDouble *nstate, FILE *statefile)
{
	// Create some variables
	int stepcount = 0, substepcount = 0, j = 0, k = 0, m = 0, interp_ret = 0;
	SpiceDouble h = 10000.0				// [s] (initial) step size
		, floating_stepcount = 0.		// not strictly the step counter
		, hp2							// [s^2] h squared
		, tEps = config_data->e_target	// [km] error allowed per step
		, tEps_p = 0.0;					// [km] temporary storage of partial error per space dimension

	config_data->saving = 0;

	// Create body arrays and set initial body positions
	SpiceDouble **(body[9]); // body[0] is t = time[1] - h, body[1] is t = time[1], ..., body[8] is t = time[1] + h

	for (k = 0; k < 9; k++)
	{
		body[k] = (SpiceDouble **)malloc(config_data->N_bodys * sizeof(SpiceDouble *));
		if (body[k] == NULL)
		{
			printf("\n\nerror: could not allocate body state array (OOM)");
			return 1;
		}
	}

	for (j = 0; j < config_data->N_bodys; j++)
	{
		for (k = 0; k < 9; k++)
		{
			body[k][j] = (SpiceDouble *)malloc(6 * sizeof(SpiceDouble));

			if (body[k][j] == NULL)
			{
				printf("\n\nerror: could not allocate body state array (OOM)");
				return 1;
			}
		}
#pragma omp critical(SPICE)
		{
			// Critical section is only executed on one thread at a time (spice is not threadsafe)
			get_body_state(config_data, j, &nstate[6], &body[8]);
		}
	}

	// Create some more variables
	SpiceDouble initPos[3], initVel[3], dir_SSB[3], time[9], dtime[9];
	SpiceDouble f[9][3], sun_dist;

	// Body coefficient variable
	SpiceDouble **body_c[6]; // positions 4-6 may not be used for order = 2, or none for order = 0
	// Allocate memory for coefficients
	if (interp_body_states_malloc(config_data, &body_c))
	{
		return 1; // OOM
	}

	dtimepowers dtp;

	dtime[0] = 0;

#ifdef __WSTEPINFO
	SpiceDouble hmin = config_data->final_time - nstate[6], hmax = 0.0, maxEps = 0;
#endif

	// Integrate
	while (nstate[6] < config_data->final_time)
	{
		// Set initial state for this step
		initPos[0] = nstate[0];
		initPos[1] = nstate[1];
		initPos[2] = nstate[2];
		initVel[0] = nstate[3];
		initVel[1] = nstate[4];
		initVel[2] = nstate[5];

		time[1] = nstate[6];

		for (j = 0; j < config_data->N_bodys; j++)
		{
			for (k = 0; k < 6; k++)
			{
				// Save previous body state (body[1] is only a useful value from the second step onwards)
				body[0][j][k] = body[1][j][k];
				// Set new initial body state (body[8] is always a useful value)
				body[1][j][k] = body[8][j][k];
			}
		}

		// F1
		dir_SSB[0] = -(initPos[0]);
		dir_SSB[1] = -(initPos[1]);
		dir_SSB[2] = -(initPos[2]);
		calc_accel(config_data, dir_SSB, &body[1], f[0], initVel, 0.);

		if ((nstate[6] + h) > config_data->start_time_save)
		{
			if (config_data->saving != 1)
			{
				config_data->saving = 1;
			}

#ifdef __SaveRateOpt
		// save rate optimization for close encounters with un-sunny bodies.
			if (calc_save_factor(config_data, dir_SSB, &body[1], f[0], initVel, 0.))
			{
				printf("\n\nerror: Sun missing.");
				return 1;
			}
#endif
		}

		// dtime: time difference compared to time[0]
		dtime[1] = time[1] - time[0];

		do
		{
			// Set dynamic step size
			if (tEps > (config_data->e_target / 22.)) // do not increase the step size by more than a factor of 1.4 at a time
			{
				h = 0.9 * h * pow(config_data->e_target / tEps, 1. / 7);
			}
			else
			{
				h = h * 1.4;
			}

			// End integration on time
			if (config_data->endontime)
			{
				if (time[1] + h > config_data->final_time)
				{
					h = config_data->final_time - time[1];
				}
			}

			// Calculate times
			dtime[2] = dtime[1] + h / 10.;
			dtime[3] = dtime[1] + h / 5.;
			dtime[4] = dtime[1] + 3. / 8. * h;
			dtime[5] = dtime[1] + h / 2.;
			dtime[6] = dtime[1] + h * ((7. - sqrt(21)) / 14.);
			dtime[7] = dtime[1] + h * ((7. + sqrt(21)) / 14.);
			dtime[8] = dtime[1] + h;

			precompute_dtime_powers(config_data, &dtp, dtime);

			for (j = 2; j < 9; j++)
			{
				time[j] = time[0] + dtime[j];
			}

			// Get body positions
			// Instead of body[9], body[8] is just used below. They would always be identical. The same applies for dtime[8/9] and time[8/9]
			
			for (j = 0; j < config_data->N_bodys; j++)
			{
				if (stepcount > 1) // not during the first two steps
				{
					// Get new end state of body
#pragma omp critical(SPICE)
					{
						get_body_state(config_data, j, &time[8], &body[8]);
					}

					// Interpolate body states
					interp_ret = interp_body_states(config_data, &body, &body_c, dtime, &dtp, h, j);
					if (interp_ret == 0) // all is well
					{
						;
					}
					else if (interp_ret == 2) // cspice everything
					{
#pragma omp critical(SPICE) // Critical section is only executed on one thread at a time (spice is not threadsafe)
						{
							for (m = 2; m < 8; m++) // body[8][j] already done above
							{
								get_body_state(config_data, j, &time[m], &body[m]);
							}
						}
					}
				}
				else // first step, no previous position available for interpolation -> get all body positions with SPICE
				{
#pragma omp critical(SPICE) // Critical section is only executed on one thread at a time (spice is not threadsafe)
					{
						//printf("\nbefore: h = %.8le, time[1] = %.8le, tEpsMax = %.8le, tEps = %.8le ",h,time[1],tEpsMax,tEps);
						for (m = 2; m < 9; m++)
						{
							get_body_state(config_data, j, &time[m], &body[m]);
						}
					}
				}
			}

			hp2 = h * h;

			// F2
			dir_SSB[0] = -(initPos[0] + h / 10 * initVel[0] + 1. / 200 * hp2 * f[0][0]);
			dir_SSB[1] = -(initPos[1] + h / 10 * initVel[1] + 1. / 200 * hp2 * f[0][1]);
			dir_SSB[2] = -(initPos[2] + h / 10 * initVel[2] + 1. / 200 * hp2 * f[0][2]);
			calc_accel(config_data, dir_SSB, &body[2], f[1], initVel, dtime[2] - dtime[1]);

			// F3
			dir_SSB[0] = -(initPos[0] + h / 5 * initVel[0] + hp2 / 150 * (f[0][0] + 2 * f[1][0]));
			dir_SSB[1] = -(initPos[1] + h / 5 * initVel[1] + hp2 / 150 * (f[0][1] + 2 * f[1][1]));
			dir_SSB[2] = -(initPos[2] + h / 5 * initVel[2] + hp2 / 150 * (f[0][2] + 2 * f[1][2]));
			calc_accel(config_data, dir_SSB, &body[3], f[2], initVel, dtime[3] - dtime[1]);

			// F4
			dir_SSB[0] = -(initPos[0] + 3. / 8. * h * initVel[0] + hp2 * (171. / 8192 * f[0][0] + 45. / 4096 * f[1][0] + 315. / 8192 * f[2][0]));
			dir_SSB[1] = -(initPos[1] + 3. / 8. * h * initVel[1] + hp2 * (171. / 8192 * f[0][1] + 45. / 4096 * f[1][1] + 315. / 8192 * f[2][1]));
			dir_SSB[2] = -(initPos[2] + 3. / 8. * h * initVel[2] + hp2 * (171. / 8192 * f[0][2] + 45. / 4096 * f[1][2] + 315. / 8192 * f[2][2]));
			calc_accel(config_data, dir_SSB, &body[4], f[3], initVel, dtime[4] - dtime[1]);

			// F5
			dir_SSB[0] = -(initPos[0] + h / 2 * initVel[0] + hp2 * (5. / 288 * f[0][0] + 25. / 528 * f[1][0] + 25. / 672 * f[2][0] + 16. / 693 * f[3][0]));
			dir_SSB[1] = -(initPos[1] + h / 2 * initVel[1] + hp2 * (5. / 288 * f[0][1] + 25. / 528 * f[1][1] + 25. / 672 * f[2][1] + 16. / 693 * f[3][1]));
			dir_SSB[2] = -(initPos[2] + h / 2 * initVel[2] + hp2 * (5. / 288 * f[0][2] + 25. / 528 * f[1][2] + 25. / 672 * f[2][2] + 16. / 693 * f[3][2]));
			calc_accel(config_data, dir_SSB, &body[5], f[4], initVel, dtime[5] - dtime[1]);

			// F6
			dir_SSB[0] = -(initPos[0] + (7. - sqrt(21)) * h / 14 * initVel[0] + hp2 * ((1003. - 205. * sqrt(21)) / 12348 * f[0][0] - 25. * (751. - 173. * sqrt(21)) / 90552 * f[1][0]
					+ 25. * (624. - 137. * sqrt(21)) / 43218 * f[2][0] - 128. * (361. - 79. * sqrt(21)) / 237699 * f[3][0] + (3411. - 745. * sqrt(21)) / 24696 * f[4][0]));
			dir_SSB[1] = -(initPos[1] + (7. - sqrt(21)) * h / 14 * initVel[1] + hp2 * ((1003. - 205. * sqrt(21)) / 12348 * f[0][1] - 25. * (751. - 173. * sqrt(21)) / 90552 * f[1][1]
					+ 25. * (624. - 137. * sqrt(21)) / 43218 * f[2][1] - 128. * (361. - 79. * sqrt(21)) / 237699 * f[3][1] + (3411. - 745. * sqrt(21)) / 24696 * f[4][1]));
			dir_SSB[2] = -(initPos[2] + (7. - sqrt(21)) * h / 14 * initVel[2] + hp2 * ((1003. - 205. * sqrt(21)) / 12348 * f[0][2] - 25. * (751. - 173. * sqrt(21)) / 90552 * f[1][2]
					+ 25. * (624. - 137. * sqrt(21)) / 43218 * f[2][2] - 128. * (361. - 79. * sqrt(21)) / 237699 * f[3][2] + (3411. - 745. * sqrt(21)) / 24696 * f[4][2]));
			calc_accel(config_data, dir_SSB, &body[6], f[5], initVel, dtime[6] - dtime[1]);

			// F7
			dir_SSB[0] = -(initPos[0] + (7. + sqrt(21)) * h / 14 * initVel[0] + hp2 * ((793. + 187. * sqrt(21)) / 12348 * f[0][0] - 25. * (331. + 113. * sqrt(21)) / 90552 * f[1][0]
					+ 25. * (1044. + 247. * sqrt(21)) / 43218 * f[2][0] - 128. * (14885. + 3779. * sqrt(21)) / 9745659 * f[3][0] + (3327. + 797. * sqrt(21)) / 24696 * f[4][0]
					- (581. + 127. * sqrt(21)) / 1722 * f[5][0]));
			dir_SSB[1] = -(initPos[1] + (7. + sqrt(21)) * h / 14 * initVel[1] + hp2 * ((793. + 187. * sqrt(21)) / 12348 * f[0][1] - 25. * (331. + 113. * sqrt(21)) / 90552 * f[1][1]
					+ 25. * (1044. + 247. * sqrt(21)) / 43218 * f[2][1] - 128. * (14885. + 3779. * sqrt(21)) / 9745659 * f[3][1] + (3327. + 797. * sqrt(21)) / 24696 * f[4][1]
					- (581. + 127. * sqrt(21)) / 1722 * f[5][1]));
			dir_SSB[2] = -(initPos[2] + (7. + sqrt(21)) * h / 14 * initVel[2] + hp2 * ((793. + 187. * sqrt(21)) / 12348 * f[0][2] - 25. * (331. + 113. * sqrt(21)) / 90552 * f[1][2]
					+ 25. * (1044. + 247. * sqrt(21)) / 43218 * f[2][2] - 128. * (14885. + 3779. * sqrt(21)) / 9745659 * f[3][2] + (3327. + 797. * sqrt(21)) / 24696 * f[4][2]
					- (581. + 127. * sqrt(21)) / 1722 * f[5][2]));
			calc_accel(config_data, dir_SSB, &body[7], f[6], initVel, dtime[7] - dtime[1]);

			// F8
			dir_SSB[0] = -(initPos[0] + h * initVel[0] + hp2 * ((-1.) * (157. - 3. * sqrt(21)) / 378 * f[0][0] + 25. * (143. - 10. * sqrt(21)) / 2772 * f[1][0]
					- 25. * (876. + 55. * sqrt(21)) / 3969 * f[2][0] + 1280. * (913. + 18. * sqrt(21)) / 596673 * f[3][0] - (1353. + 26. * sqrt(21)) / 2268 * f[4][0]
					+ 7. * (1777. + 377. * sqrt(21)) / 4428 * f[5][0] + 7. * (5. - sqrt(21)) / 36 * f[6][0]));
			dir_SSB[1] = -(initPos[1] + h * initVel[1] + hp2 * ((-1.) * (157. - 3. * sqrt(21)) / 378 * f[0][1] + 25. * (143. - 10. * sqrt(21)) / 2772 * f[1][1]
					- 25. * (876. + 55. * sqrt(21)) / 3969 * f[2][1] + 1280. * (913. + 18. * sqrt(21)) / 596673 * f[3][1] - (1353. + 26. * sqrt(21)) / 2268 * f[4][1]
					+ 7. * (1777. + 377. * sqrt(21)) / 4428 * f[5][1] + 7. * (5. - sqrt(21)) / 36 * f[6][1]));
			dir_SSB[2] = -(initPos[2] + h * initVel[2] + hp2 * ((-1.) * (157. - 3. * sqrt(21)) / 378 * f[0][2] + 25. * (143. - 10. * sqrt(21)) / 2772 * f[1][2]
					- 25. * (876. + 55. * sqrt(21)) / 3969 * f[2][2] + 1280. * (913. + 18. * sqrt(21)) / 596673 * f[3][2] - (1353. + 26. * sqrt(21)) / 2268 * f[4][2]
					+ 7. * (1777. + 377. * sqrt(21)) / 4428 * f[5][2] + 7. * (5. - sqrt(21)) / 36 * f[6][2]));
			calc_accel(config_data, dir_SSB, &body[8], f[7], initVel, dtime[8] - dtime[1]);

			// F9 - only used for error calculation
			dir_SSB[0] = -(initPos[0] + h * initVel[0] + hp2 * (1. / 20 * f[0][0] + 8. / 45 * f[4][0] + 7. * (7. + sqrt(21)) / 360 * f[5][0] + 7. * (7. - sqrt(21)) / 360 * f[6][0]));
			dir_SSB[1] = -(initPos[1] + h * initVel[1] + hp2 * (1. / 20 * f[0][1] + 8. / 45 * f[4][1] + 7. * (7. + sqrt(21)) / 360 * f[5][1] + 7. * (7. - sqrt(21)) / 360 * f[6][1]));
			dir_SSB[2] = -(initPos[2] + h * initVel[2] + hp2 * (1. / 20 * f[0][2] + 8. / 45 * f[4][2] + 7. * (7. + sqrt(21)) / 360 * f[5][2] + 7. * (7. - sqrt(21)) / 360 * f[6][2]));
			calc_accel(config_data, dir_SSB, &body[8], f[8], initVel, dtime[8] - dtime[1]);

			// Absolute error (2-norm)
			tEps_p = (f[7][0] - f[8][0]);
			tEps = tEps_p * tEps_p;
			tEps_p = (f[7][1] - f[8][1]);
			tEps += tEps_p * tEps_p;
			tEps_p = (f[7][2] - f[8][2]);
			tEps += tEps_p * tEps_p;
			tEps = hp2 / 20 * sqrt(tEps);

			substepcount++;
		} while (tEps > config_data->e_target); // while the error is too big.

#ifdef __WSTEPINFO
		if (tEps > maxEps)
		{
			maxEps = tEps;
		}
#endif

		// Update solution
		// x7th_i+1
		nstate[0] = initPos[0] + h * initVel[0] + hp2 * (18. * f[0][0] + 64. * f[4][0] + (7. * (7. + sqrt(21))) * f[5][0] + (7. * (7. - sqrt(21))) * f[6][0]) / 360;
		nstate[1] = initPos[1] + h * initVel[1] + hp2 * (18. * f[0][1] + 64. * f[4][1] + (7. * (7. + sqrt(21))) * f[5][1] + (7. * (7. - sqrt(21))) * f[6][1]) / 360;
		nstate[2] = initPos[2] + h * initVel[2] + hp2 * (18. * f[0][2] + 64. * f[4][2] + (7. * (7. + sqrt(21))) * f[5][2] + (7. * (7. - sqrt(21))) * f[6][2]) / 360;
		// v7th_i+1
		nstate[3] = initVel[0] + h * (18. * (f[0][0] + f[7][0]) + 128. * f[4][0] + 98. * f[5][0] + 98. * f[6][0]) / 360;
		nstate[4] = initVel[1] + h * (18. * (f[0][1] + f[7][1]) + 128. * f[4][1] + 98. * f[5][1] + 98. * f[6][1]) / 360;
		nstate[5] = initVel[2] + h * (18. * (f[0][2] + f[7][2]) + 128. * f[4][2] + 98. * f[5][2] + 98. * f[6][2]) / 360;

		/* x6th_i+1, not necessary to compute
		nstate[3] = initVel[0] + h * (18 * f[0][0] + 64 * f[4][0] + (7 * (7 + sqrt(21))) * f[5][0] + (7 * (7 - sqrt(21))) * f[6][0] - 18 * f[7][0] + 18 * f[8][0]) / 360;
		nstate[4] = initVel[1] + h * (18 * f[0][1] + 64 * f[4][1] + (7 * (7 + sqrt(21))) * f[5][1] + (7 * (7 - sqrt(21))) * f[6][1] - 18 * f[7][1] + 18 * f[8][1]) / 360;
		nstate[5] = initVel[2] + h * (18 * f[0][2] + 64 * f[4][2] + (7 * (7 + sqrt(21))) * f[5][2] + (7 * (7 - sqrt(21))) * f[6][2] - 18 * f[7][2] + 18 * f[8][2]) / 360; */

		// Update time
		nstate[6] = time[1] + h;
		// Save previous time for body location interpolation
		time[0] = time[1];

#ifdef __WSTEPINFO
		if (h < hmin) // calculate smallest time step
		{
			hmin = h;
		}
		else if (h > hmax) // calculate largest time step
		{
			hmax = h;
		}
#endif

		stepcount++;

		// Increase StepCount
#ifdef __SaveRateOpt
		floating_stepcount += config_data->step_multiplier;
#else
		floating_stepcount += 1.;
#endif // __SaveRateOpt

		// If distance to the Sun is less than 10 sun radii terminate particle by setting coordinates to 99 and breaking integration loop
		sun_dist = sqrt(nstate[0] * nstate[0] + nstate[1] * nstate[1] + nstate[2] * nstate[2]);
		if (sun_dist < 7000000)
		{
			nstate[0] = 99;
			nstate[1] = 99;
			nstate[2] = 99;
			nstate[3] = 99;
			nstate[4] = 99;
			nstate[5] = 99;
			break;
		}

		// Save nth state
		if (config_data->saving == 1)
		{
			if (floating_stepcount >= config_data->n)
			{
				printpdata(statefile, nstate);
				floating_stepcount = 0.;
			}
		}
	}

	// Print last state to file and close file
	if (floating_stepcount != 0.)
	{
		printpdata(statefile, nstate);
	}

	// Deallocate body array
	for (k = 0; k < 9; k++)
	{
		for (j = 0; j < config_data->N_bodys; j++)
		{
			free(body[k][j]);
		}
		free(body[k]);
	}

	// Free body state coefficient memory
	interp_body_states_free(config_data, &body_c);

#ifdef __WSTEPINFO
	printf("\n  Smallest time step: %.4le s", hmin);
	printf("  -  Largest time step: %.4le s", hmax);
	printf("  -  Total number of steps: %d", stepcount);
	printf("  -  Total number of substeps: %d", substepcount);
	printf("  -  Largest error per step: %.4le km", maxEps);
#endif

	return 0;
}