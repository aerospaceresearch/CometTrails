/* 7th and 6th order Runge-Kutta method as described in:
   J.R. Dormand & P.J. Prince (1978): New Runge-Kutta algorithms for numerical simulation in dynamical astronomy. Celestial Mechanics, Vol. 18, p. 223-232. 
   
   Steps exceeding the maximum allowed error (e_target) will be repeated. */

int RungeKutta67(configuration_values *config_out, SpiceDouble *nstate, FILE *statefile)
{
	// Select body position function to use ((*bodyPosFP) for spice, return_SSB for (0,0,0))
	void (*bodyPosFP)(SpiceInt, SpiceDouble, ConstSpiceChar *, ConstSpiceChar *, SpiceInt, SpiceDouble[3], SpiceDouble *);
	if (config_out->ssb_centered == 1)
	{
		bodyPosFP = &return_SSB;
	}
	else
	{
		bodyPosFP = &spkezp_c;
	}

	// Create some variables
	int stepcount = 0, substepcount = 0, j = 0, k = 0, m = 0, zeroEps = 0;
	SpiceDouble lt, h = 3000.0, hp2, tEps = config_out->e_target, tEps_p = 0.0;

	// Create body arrays and set initial body positions
	SpiceDouble **(body[10]); // body[0] is t = time[1] - h, body[1] is t = time[1], body[9] is t = time[1] + h

	for (k = 0; k < 10; k++)
	{
		body[k] = (SpiceDouble **)malloc(config_out->N_bodys * sizeof(SpiceDouble *));
		if (body[k] == NULL)
		{
			printf("\nerror: could not allocate body state array (OOM)");
			return 1;
		}
	}

	for (j = 0; j < config_out->N_bodys; j++)
	{
		for (k = 0; k < 10; k++)
		{
			body[k][j] = (SpiceDouble *)malloc(3 * sizeof(SpiceDouble));

			if (body[k][j] == NULL)
			{
				printf("\nerror: could not allocate body state array (OOM)");
				return 1;
			}
		}
#pragma omp critical(SPICE)
		{
			// Critical section is only executed on one thread at a time (spice is not threadsafe)
			(*bodyPosFP)(config_out->body_int[j], nstate[6], "ECLIPJ2000", "NONE", 0, body[8][j], &lt);
		}
	}

	// Create body interpolation coefficients
	SpiceDouble **bod_a, **bod_b, **bod_c;
	bod_a = malloc(config_out->N_bodys * sizeof(SpiceDouble *));
	bod_b = malloc(config_out->N_bodys * sizeof(SpiceDouble *));
	bod_c = malloc(config_out->N_bodys * sizeof(SpiceDouble *));
	if (bod_a == NULL || bod_b == NULL || bod_c == NULL)
	{
		printf("\nerror: could not allocate body coefficient array (OOM)");
		return 1;
	}
	for (j = 0; j < config_out->N_bodys; j++)
	{
		bod_a[j] = malloc(3 * sizeof(SpiceDouble));
		bod_b[j] = malloc(3 * sizeof(SpiceDouble));
		bod_c[j] = malloc(3 * sizeof(SpiceDouble));
		if (bod_a[j] == NULL || bod_b[j] == NULL || bod_c[j] == NULL)
		{
			printf("\nerror: could not allocate body coefficient array (OOM)");
			return 1;
		}
	}

	// Create some more variables
	SpiceDouble initPos[3], initVel[3], dir_SSB[3], time[9], dtime[9];
	SpiceDouble f[9][3];

	// Constant part of PRD
	SpiceDouble PRDconst = calc_prdc(config_out);

	dtime[0] = 0;

#ifdef __WTIMESTEP
	SpiceDouble hmin = config_out->final_time - nstate[6], hmax = 0.0, maxEps = 0;
#endif

	// Integrate
	while (nstate[6] < config_out->final_time)
	{
		// Set initial state for this step
		initPos[0] = nstate[0];
		initPos[1] = nstate[1];
		initPos[2] = nstate[2];
		initVel[0] = nstate[3];
		initVel[1] = nstate[4];
		initVel[2] = nstate[5];

		time[1] = nstate[6];

		for (j = 0; j < config_out->N_bodys; j++)
		{
			// Save previous body state
			body[0][j][0] = body[1][j][0];
			body[0][j][1] = body[1][j][1];
			body[0][j][2] = body[1][j][2];

			// Set new initial body state
			body[1][j][0] = body[8][j][0];
			body[1][j][1] = body[8][j][1];
			body[1][j][2] = body[8][j][2];
		}

		// F1
		dir_SSB[0] = -(initPos[0]);
		dir_SSB[1] = -(initPos[1]);
		dir_SSB[2] = -(initPos[2]);
		calc_accel(config_out, dir_SSB, &body[1], f[0], initVel, PRDconst);

		// dtime: time difference compared to time[0]
		dtime[1] = time[1] - time[0];

		do
		{
			// Set dynamic step size
			if (tEps > (config_out->e_target / 22.)) // do not increase the step size by more than a factor of 1.4 at a time
			{
				h = 0.9 * h * pow(config_out->e_target / tEps, 1. / 7);
			}
			else
			{
				h = h * 1.4;
				zeroEps = 1;
			}

			// calculate times
			dtime[2] = dtime[1] + h / 10;
			dtime[3] = dtime[1] + h / 5;
			dtime[4] = dtime[1] + 3. * h / 8;
			dtime[5] = dtime[1] + h / 2;
			dtime[6] = dtime[1] + h * ((7 - sqrt(21)) / 14);
			dtime[7] = dtime[1] + h * ((7 + sqrt(21)) / 14);
			dtime[8] = dtime[1] + h;

			for (j = 2; j < 9; j++)
			{
				time[j] = time[0] + dtime[j];
			}

			// Get body positions
			for (j = 0; j < config_out->N_bodys; j++)
			{
				if (stepcount > 1) // not during the first two steps
				{
#pragma omp critical(SPICE)
					{
						(*bodyPosFP)(config_out->body_int[j], time[8], "ECLIPJ2000", "NONE", 0, body[8][j], &lt);
					}

					// solving x = a + bt + ct^2 for quadratic interpolation of body positions
					for (k = 0; k < 3; k++) // loop x,y,z
					{
						bod_a[j][k] = body[0][j][k];

						bod_c[j][k] = (((body[8][j][k] - bod_a[j][k]) / dtime[8]) - ((body[1][j][k] - bod_a[j][k]) / dtime[1])) / h;

						bod_b[j][k] = (body[1][j][k] - bod_a[j][k]) / dtime[1] - bod_c[j][k] * dtime[1];

						// interpolate body states
						for (m = 2; m < 8; m++)
						{
							body[m][j][k] = bod_a[j][k] + (bod_b[j][k] + bod_c[j][k]) * dtime[m];
						}
					}
				}
				else // first step, no previous position available for interpolation -> get all body positions with SPICE
				{
#pragma omp critical(SPICE)
					{
						//printf("\nbefore: h = %.8le, time[1] = %.8le, tEpsMax = %.8le, tEps = %.8le ",h,time[1],tEpsMax,tEps);
						// Critical section is only executed on one thread at a time (spice is not threadsafe)
						for (m = 2; m < 8; m++)
						{
							(*bodyPosFP)(config_out->body_int[j], time[m], "ECLIPJ2000", "NONE", 0, body[m][j], &lt);
						}
					}
				}
				body[9][j][0] = body[8][j][0]; // same time
				body[9][j][1] = body[8][j][1];
				body[9][j][2] = body[8][j][2];
			}

			hp2 = h * h;

			// F2
			dir_SSB[0] = -(initPos[0] + h / 10 * initVel[0] + 1. / 200 * hp2 * f[0][0]);
			dir_SSB[1] = -(initPos[1] + h / 10 * initVel[1] + 1. / 200 * hp2 * f[0][1]);
			dir_SSB[2] = -(initPos[2] + h / 10 * initVel[2] + 1. / 200 * hp2 * f[0][2]);
			calc_accel(config_out, dir_SSB, &body[2], f[1], initVel, PRDconst);

			// F3
			dir_SSB[0] = -(initPos[0] + h / 5 * initVel[0] + hp2 / 150 * (f[0][0] + 2 * f[1][0]));
			dir_SSB[1] = -(initPos[1] + h / 5 * initVel[1] + hp2 / 150 * (f[0][1] + 2 * f[1][1]));
			dir_SSB[2] = -(initPos[2] + h / 5 * initVel[2] + hp2 / 150 * (f[0][2] + 2 * f[1][2]));
			calc_accel(config_out, dir_SSB, &body[3], f[2], initVel, PRDconst);

			// F4
			dir_SSB[0] = -(initPos[0] + 3. * h / 8 * initVel[0] + hp2 * (171. / 8192 * f[0][0] + 45. / 4096 * f[1][0] + 315. / 8192 * f[2][0]));
			dir_SSB[1] = -(initPos[1] + 3. * h / 8 * initVel[1] + hp2 * (171. / 8192 * f[0][1] + 45. / 4096 * f[1][1] + 315. / 8192 * f[2][1]));
			dir_SSB[2] = -(initPos[2] + 3. * h / 8 * initVel[2] + hp2 * (171. / 8192 * f[0][2] + 45. / 4096 * f[1][2] + 315. / 8192 * f[2][2]));
			calc_accel(config_out, dir_SSB, &body[4], f[3], initVel, PRDconst);

			// F5
			dir_SSB[0] = -(initPos[0] + h / 2 * initVel[0] + hp2 * (5. / 288 * f[0][0] + 25. / 528 * f[1][0] + 25. / 672 * f[2][0] + 16. / 693 * f[3][0]));
			dir_SSB[1] = -(initPos[1] + h / 2 * initVel[1] + hp2 * (5. / 288 * f[0][1] + 25. / 528 * f[1][1] + 25. / 672 * f[2][1] + 16. / 693 * f[3][1]));
			dir_SSB[2] = -(initPos[2] + h / 2 * initVel[2] + hp2 * (5. / 288 * f[0][2] + 25. / 528 * f[1][2] + 25. / 672 * f[2][2] + 16. / 693 * f[3][2]));
			calc_accel(config_out, dir_SSB, &body[5], f[4], initVel, PRDconst);

			// F6
			dir_SSB[0] = -(initPos[0] + (7. - sqrt(21)) * h / 14 * initVel[0] + hp2 * ((1003. - 205. * sqrt(21)) / 12348 * f[0][0] - 25. * (751. - 173. * sqrt(21)) / 90552 * f[1][0]
					+ 25. * (624. - 137. * sqrt(21)) / 43218 * f[2][0] - 128. * (361. - 79. * sqrt(21)) / 237699 * f[3][0] + (3411. - 745. * sqrt(21)) / 24696 * f[4][0]));
			dir_SSB[1] = -(initPos[1] + (7. - sqrt(21)) * h / 14 * initVel[1] + hp2 * ((1003. - 205. * sqrt(21)) / 12348 * f[0][1] - 25. * (751. - 173. * sqrt(21)) / 90552 * f[1][1]
					+ 25. * (624. - 137. * sqrt(21)) / 43218 * f[2][1] - 128. * (361. - 79. * sqrt(21)) / 237699 * f[3][1] + (3411. - 745. * sqrt(21)) / 24696 * f[4][1]));
			dir_SSB[2] = -(initPos[2] + (7. - sqrt(21)) * h / 14 * initVel[2] + hp2 * ((1003. - 205. * sqrt(21)) / 12348 * f[0][2] - 25. * (751. - 173. * sqrt(21)) / 90552 * f[1][2]
					+ 25. * (624. - 137. * sqrt(21)) / 43218 * f[2][2] - 128. * (361. - 79. * sqrt(21)) / 237699 * f[3][2] + (3411. - 745. * sqrt(21)) / 24696 * f[4][2]));
			calc_accel(config_out, dir_SSB, &body[6], f[5], initVel, PRDconst);

			// F7
			dir_SSB[0] = -(initPos[0] + (7. + sqrt(21)) * h / 14 * initVel[0] + hp2 * ((793. + 187. * sqrt(21)) / 12348 * f[0][0] - 25. * (331. + 113. * sqrt(21)) / 90552 * f[1][0]
					+ 25. * (1044. + 247. * sqrt(21)) / 43218 * f[2][0] - 128. * (14885. + 3779. * sqrt(21)) / 9745659 * f[3][0] + (3327. + 797. * sqrt(21)) / 24696 * f[4][0]
					+ (581. + 127. * sqrt(21)) / 1722 * f[5][0]));
			dir_SSB[1] = -(initPos[1] + (7. + sqrt(21)) * h / 14 * initVel[1] + hp2 * ((793. + 187. * sqrt(21)) / 12348 * f[0][1] - 25. * (331. + 113. * sqrt(21)) / 90552 * f[1][1]
					+ 25. * (1044. + 247. * sqrt(21)) / 43218 * f[2][1] - 128. * (14885. + 3779. * sqrt(21)) / 9745659 * f[3][1] + (3327. + 797. * sqrt(21)) / 24696 * f[4][1]
					+ (581. + 127. * sqrt(21)) / 1722 * f[5][1]));
			dir_SSB[2] = -(initPos[2] + (7. + sqrt(21)) * h / 14 * initVel[2] + hp2 * ((793. + 187. * sqrt(21)) / 12348 * f[0][2] - 25. * (331. + 113. * sqrt(21)) / 90552 * f[1][2]
					+ 25. * (1044. + 247. * sqrt(21)) / 43218 * f[2][2] - 128. * (14885. + 3779. * sqrt(21)) / 9745659 * f[3][2] + (3327. + 797. * sqrt(21)) / 24696 * f[4][2]
					+ (581. + 127. * sqrt(21)) / 1722 * f[5][2]));
			calc_accel(config_out, dir_SSB, &body[7], f[6], initVel, PRDconst);

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
			calc_accel(config_out, dir_SSB, &body[8], f[7], initVel, PRDconst);
			//printf("\nf8 - dir_SSB[0]: %.16le", dir_SSB[0]);

			// F9 - only used for error calculation
			dir_SSB[0] = -(initPos[0] + h * initVel[0] + hp2 * (1. / 20 * f[0][0] + 8. / 45 * f[4][0] + 7. * (7. + sqrt(21)) / 360 * f[5][0] + 7. * (7. - sqrt(21)) / 360 * f[6][0]));
			dir_SSB[1] = -(initPos[1] + h * initVel[1] + hp2 * (1. / 20 * f[0][1] + 8. / 45 * f[4][1] + 7. * (7. + sqrt(21)) / 360 * f[5][1] + 7. * (7. - sqrt(21)) / 360 * f[6][1]));
			dir_SSB[2] = -(initPos[2] + h * initVel[2] + hp2 * (1. / 20 * f[0][2] + 8. / 45 * f[4][2] + 7. * (7. + sqrt(21)) / 360 * f[5][2] + 7. * (7. - sqrt(21)) / 360 * f[6][2]));
			calc_accel(config_out, dir_SSB, &body[9], f[8], initVel, PRDconst);
			//printf("\nf9 - dir_SSB[0]: %.16le", dir_SSB[0]);

			// Absolute error (2-norm)
			tEps_p = (f[7][0] - f[8][0]);
			//printf("\ntEps_p = %.8le, f[7][0] = %.8le, f[8][0] = %.8le", tEps_p, f[7][0], f[8][0]);
			tEps = tEps_p * tEps_p;
			tEps_p = (f[7][1] - f[8][1]);
			tEps += tEps_p * tEps_p;
			tEps_p = (f[7][2] - f[8][2]);
			tEps += tEps_p * tEps_p;
			tEps = hp2 / 20 * sqrt(tEps);

			substepcount++;
		} while (tEps > config_out->e_target); // while the error is too big.

#ifdef __WTIMESTEP
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
		// Save previous time for body interpolation
		time[0] = time[1];

#ifdef __WTIMESTEP
		if (h < hmin) // calculate smallest time step
		{
			hmin = h;
		}
		else if (h > hmax) // calculate larges time step
		{
			hmax = h;
		}
#endif

		// Increase StepCount
		stepcount++;

		// Save nth state
		if ((stepcount % config_out->n) == 0)
		{
			//printf("\n tEps = %.12le", tEps);
			if (nstate[6] > config_out->start_time_save)
			{
				printpdata(statefile, nstate);
			}
		}
	}

	// Print last state to file and close file
	if (stepcount /*!= 0*/)
	{
		printpdata(statefile, nstate);
	}

	// Deallocate body array
	for (k = 0; k < 10; k++)
	{
		for (j = 0; j < config_out->N_bodys; j++)
		{
			free(body[k][j]);
		}
		free(body[k]);
	}

#ifdef __WTIMESTEP
	if (zeroEps)
	{
		printf("\n            	 Warning: Epsilon was limited.");
	}
	printf("\n            	 Smallest time step: %.4le s", hmin);
	printf("  -  Largest time step: %.4le s", hmax);
	printf("  -  Total number of steps: %d", stepcount);
	printf("  -  Total number of substeps: %d", substepcount);
	printf("  -  Largest error per step: %.4le km", maxEps);
#endif

	return 0;
}