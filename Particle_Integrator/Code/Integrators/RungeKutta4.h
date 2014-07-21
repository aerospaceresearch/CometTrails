/* 4th order Runge-Kutta method */

int RungeKutta4(configuration_values *config_data, SpiceDouble *nstate, FILE *statefile)
{
	//Create some variables
	int j;
	SpiceDouble dt					// [s] time step
		, dt2						// [s] dt/2
		, floating_stepcount = 0.;	// not strictly the step counter

	config_data->saving = 0;
	
	//Create body arrays and set initial body positions
	SpiceDouble **body_pre, **body_mid, **body_end;
	body_pre = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
	body_mid = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
	body_end = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
	if (body_pre == NULL || body_mid == NULL || body_end == NULL)
	{
		printf("\n\nerror: could not allocate body state array (OOM)");
		return 1;
	}
	for (j = 0; j < config_data->N_bodys; j++)
	{
		body_pre[j] = malloc(3 * sizeof(SpiceDouble));
		body_mid[j] = malloc(3 * sizeof(SpiceDouble));
		body_end[j] = malloc(3 * sizeof(SpiceDouble));
		if (body_pre[j] == NULL || body_mid[j] == NULL || body_end[j] == NULL)
		{
			printf("\n\nerror: could not allocate body state array (OOM)");
			return 1;
		}
#pragma omp critical(SPICE)
		{
			//Critical section is only executed on one thread at a time (spice is not threadsafe)
			get_body_state(config_data, j, &nstate[6], &body_end);
		}
	}

	//Create some more variables
	SpiceDouble initPos[3], initVel[3], dir_SSB[3], initTime, nextInitTime, abs_acc;
	SpiceDouble k_acc_1[3], k_acc_2[3], k_acc_3[3], k_acc_4[3];
	SpiceDouble k_vel_1[3], k_vel_2[3], k_vel_3[3], k_vel_4[3];

#ifdef __WSTEPINFO
	SpiceDouble dtmin = config_data->final_time - nstate[6], dtmax = 0.0;
	int stepcount = 0;
#endif

	//Integrate
	while (nstate[6] < config_data->final_time)
	{

		//Set initial state for this step
		initPos[0] = nstate[0];
		initPos[1] = nstate[1];
		initPos[2] = nstate[2];
		initVel[0] = nstate[3];
		initVel[1] = nstate[4];
		initVel[2] = nstate[5];
		initTime = nstate[6];
		for (j = 0; j < config_data->N_bodys; j++)
		{
			body_pre[j][0] = body_end[j][0];
			body_pre[j][1] = body_end[j][1];
			body_pre[j][2] = body_end[j][2];
		}

		//Step 1
		dir_SSB[0] = -initPos[0];
		dir_SSB[1] = -initPos[1];
		dir_SSB[2] = -initPos[2];
		calc_accel(config_data, dir_SSB, &body_pre, k_acc_1, initVel, 0.0);
		k_vel_1[0] = initVel[0];
		k_vel_1[1] = initVel[1];
		k_vel_1[2] = initVel[2];

		//Set dynamic step size
		abs_acc = sqrt(k_acc_1[0] * k_acc_1[0] + k_acc_1[1] * k_acc_1[1] + k_acc_1[2] * k_acc_1[2]);
		dt = (config_data->dv_step / abs_acc);

		if (nstate[6] + dt > config_data->start_time_save)
		{
			if (config_data->saving != 1)
			{
				config_data->saving = 1;
			}

#ifdef __SaveRateOpt
		// save rate optimization for close encounters with un-sunny bodies.
			if (calc_save_factor(config_data, dir_SSB, &body_pre, k_acc_1, initVel, 0.0))
			{
				printf("\n\nerror: Sun missing.");
				return 1;
			}
#endif
		}

#ifdef __WSTEPINFO
		if (dt < dtmin) // calculate smallest time step
		{
			dtmin = dt;
		}
		else if (dt > dtmax) // calculate larges time step
		{
			dtmax = dt;
		}
		stepcount++;
#endif // __WSTEPINFO

		// End integration on time
		if (config_data->endontime)
		{
			if ((initTime + dt) > config_data->final_time)
			{
				dt = config_data->final_time - initTime;
			}
		}

		dt2 = dt / 2;
		nextInitTime = initTime + dt;

		//Get body positions with SPICE
		for (j = 0; j < config_data->N_bodys; j++)
		{
#pragma omp critical(SPICE)
			{
				//Critical section is only executed on one thread at a time (spice is not threadsafe)
				get_body_state(config_data, j, &nextInitTime, &body_end);
			} // ~94% of all computing time is spent here, mostly in spkgps
			body_mid[j][0] = (body_pre[j][0] + body_end[j][0]) / 2;
			body_mid[j][1] = (body_pre[j][1] + body_end[j][1]) / 2;
			body_mid[j][2] = (body_pre[j][2] + body_end[j][2]) / 2;
		}

		//Step 2
		dir_SSB[0] = -(initPos[0] + k_vel_1[0] * dt2);
		dir_SSB[1] = -(initPos[1] + k_vel_1[1] * dt2);
		dir_SSB[2] = -(initPos[2] + k_vel_1[2] * dt2);
		calc_accel(config_data, dir_SSB, &body_mid, k_acc_2, initVel, dt2);
		k_vel_2[0] = initVel[0] + k_acc_1[0] * dt2;
		k_vel_2[1] = initVel[1] + k_acc_1[1] * dt2;
		k_vel_2[2] = initVel[2] + k_acc_1[2] * dt2;

		//Step 3
		dir_SSB[0] = -(initPos[0] + k_vel_2[0] * dt2);
		dir_SSB[1] = -(initPos[1] + k_vel_2[1] * dt2);
		dir_SSB[2] = -(initPos[2] + k_vel_2[2] * dt2);
		calc_accel(config_data, dir_SSB, &body_mid, k_acc_3, initVel, dt2);
		k_vel_3[0] = initVel[0] + k_acc_2[0] * dt2;
		k_vel_3[1] = initVel[1] + k_acc_2[1] * dt2;
		k_vel_3[2] = initVel[2] + k_acc_2[2] * dt2;

		//Step 4
		dir_SSB[0] = -(initPos[0] + k_vel_3[0] * dt);
		dir_SSB[1] = -(initPos[1] + k_vel_3[1] * dt);
		dir_SSB[2] = -(initPos[2] + k_vel_3[2] * dt);
		calc_accel(config_data, dir_SSB, &body_end, k_acc_4, initVel, dt);
		k_vel_4[0] = initVel[0] + k_acc_3[0] * dt;
		k_vel_4[1] = initVel[1] + k_acc_3[1] * dt;
		k_vel_4[2] = initVel[2] + k_acc_3[2] * dt;

		//Update solution
		nstate[0] = initPos[0] + dt * (k_vel_1[0] + 2 * k_vel_2[0] + 2 * k_vel_3[0] + k_vel_4[0]) / 6;
		nstate[1] = initPos[1] + dt * (k_vel_1[1] + 2 * k_vel_2[1] + 2 * k_vel_3[1] + k_vel_4[1]) / 6;
		nstate[2] = initPos[2] + dt * (k_vel_1[2] + 2 * k_vel_2[2] + 2 * k_vel_3[2] + k_vel_4[2]) / 6;
		nstate[3] = initVel[0] + dt * (k_acc_1[0] + 2 * k_acc_2[0] + 2 * k_acc_3[0] + k_acc_4[0]) / 6;
		nstate[4] = initVel[1] + dt * (k_acc_1[1] + 2 * k_acc_2[1] + 2 * k_acc_3[1] + k_acc_4[1]) / 6;
		nstate[5] = initVel[2] + dt * (k_acc_1[2] + 2 * k_acc_2[2] + 2 * k_acc_3[2] + k_acc_4[2]) / 6;
		nstate[6] = nextInitTime;

		// Increase StepCount
#ifdef __SaveRateOpt
		floating_stepcount += config_data->step_multiplier;
#else
		floating_stepcount += 1.;
#endif // __SaveRateOpt

		// Save nth state
		if (floating_stepcount >= config_data->n)
		{
			if (config_data->saving == 1)
			{
				printpdata(statefile, nstate);
			}
			floating_stepcount = 0.;
		}
	}

	//Print last state to file and close file
	if (floating_stepcount != 0.)
	{
		printpdata(statefile, nstate);
	}

	//Deallocate body array
	for (j = 0; j < config_data->N_bodys; j++)
	{
		free(body_pre[j]);
		free(body_mid[j]);
		free(body_end[j]);
	}
	free(body_pre);
	free(body_mid);
	free(body_end);

#ifdef __WSTEPINFO
	printf("\n  Smallest time step: %.6le s", dtmin);
	printf("  -  Largest time step: %.6le s", dtmax);
	printf("  -  Total number of steps: %d", stepcount);
#endif

	return 0;
}