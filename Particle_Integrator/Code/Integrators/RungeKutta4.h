/* 4th order Runge-Kutta method */

int RungeKutta4(configuration_values *config_data, SpiceDouble *nstate, FILE *statefile)
{
	// Select body position function to use ((*bodyPosFP) for spice, return_SSB for (0,0,0))
	void(*bodyPosFP)(SpiceInt, SpiceDouble, ConstSpiceChar *, ConstSpiceChar *, SpiceInt, SpiceDouble[3], SpiceDouble *);
	if (config_data->ssb_centered == 1)
	{
		bodyPosFP = &return_SSB;
	}
	else
	{
		bodyPosFP = &spkezp_c;
	}

	//Create some variables
	int j, i = 0;
	SpiceDouble lt, dt, dt2;

	//Create body arrays and set initial body positions
	SpiceDouble **body_pre, **body_mid, **body_end;
	body_pre = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
	body_mid = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
	body_end = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
	if (body_pre == NULL || body_mid == NULL || body_end == NULL)
	{
		printf("\nerror: could not allocate body state array (OOM)");
		return 1;
	}
	for (j = 0; j < config_data->N_bodys; j++)
	{
		body_pre[j] = malloc(3 * sizeof(SpiceDouble));
		body_mid[j] = malloc(3 * sizeof(SpiceDouble));
		body_end[j] = malloc(3 * sizeof(SpiceDouble));
		if (body_pre[j] == NULL || body_mid[j] == NULL || body_end[j] == NULL)
		{
			printf("\nerror: could not allocate body state array (OOM)");
			return 1;
		}
#pragma omp critical(SPICE)
		{
			//Critical section is only executed on one thread at a time (spice is not threadsafe)
			(*bodyPosFP)(config_data->body_int[j], nstate[6], "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
		}
	}

	//Create some more variables
	SpiceDouble initPos[3], initVel[3], dir_SSB[3], initTime, abs_acc;
	SpiceDouble k_acc_1[3], k_acc_2[3], k_acc_3[3], k_acc_4[3];
	SpiceDouble k_vel_1[3], k_vel_2[3], k_vel_3[3], k_vel_4[3];

#ifdef __WTIMESTEP
	SpiceDouble dtmin = config_data->final_time - nstate[6], dtmax = 0.0;
	int numsteps = 0;
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
		dt2 = dt / 2;

#ifdef __WTIMESTEP
		if (dt < dtmin) // calculate smallest time step
		{
			dtmin = dt;
		}
		else if (dt > dtmax) // calculate larges time step
		{
			dtmax = dt;
		}
		numsteps++;
#endif // __WTIMESTEP

#ifdef __ENDONTIME
		if (initTime + dt > config_data->final_time)
		{
			dt = config_data->final_time - initTime;
		}
#endif // __ENDONTIME

		//Get body positions with SPICE
		for (j = 0; j < config_data->N_bodys; j++)
		{
#pragma omp critical(SPICE)
			{
				//Critical section is only executed on one thread at a time (spice is not threadsafe)
				(*bodyPosFP)(config_data->body_int[j], initTime + dt, "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
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
		nstate[6] = initTime + dt;

		//Increase StepCount
		i++;

		//Save nth state
		if (i == config_data->n)
		{
			if (nstate[6] > config_data->start_time_save)
			{
				printpdata(statefile, nstate);
			}
			i = 0;
		}
	}

	//Print last state to file and close file
	if (i /*!= 0*/)
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

#ifdef __WTIMESTEP
	printf("\n   Smallest time step: %.6le s", dtmin);
	printf("  -  Largest time step: %.6le s", dtmax);
	printf("  -  Total number of steps: %d k", numsteps / 1000);
#endif

	return 0;
}