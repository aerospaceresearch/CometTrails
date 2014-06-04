// Functions called by integration algorithms

// Solar wind drag/Poynting-Robertson drag factor
#ifdef __PRD
	#ifdef __SWD
		#define SWDF 1.35	// PRD and SWD: (1 + sw) with sw = 0.35
	#else
		#define SWDF 1.		// PRD, no SWD
	#endif // SWD
#else // no PRD
	#ifdef __SWD
		#define SWDF 0.35	// SWD, no PRD
	#else
		#define SWDF 0.		// neither SWD nor PRD
	#endif // SWD
#endif // PRD

/* Calculate the acceleration of a particle based on the position of the body relative to the SSB */
void calc_accel(configuration_values *config_data, SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel, SpiceDouble *Vel, SpiceDouble dt)
{
	SpiceDouble r_body[3]	// [km]
		, absr = 0.			// [km]
		, r3				// [km^3]
		, GMr3				// [1/s^2]
		, iVel[3];			// Intermediate Speed [km/s]
	int b; // body

	accel[0] = 0;
	accel[1] = 0;
	accel[2] = 0;

	for (b = 0; b < config_data->N_bodys; b++)
	{
		r_body[0] = (*body_state)[b][0] + dir_SSB[0];
		r_body[1] = (*body_state)[b][1] + dir_SSB[1];
		r_body[2] = (*body_state)[b][2] + dir_SSB[2];

		// Calculate GM*r^3
		absr = sqrt(r_body[0] * r_body[0] + r_body[1] * r_body[1] + r_body[2] * r_body[2]); // abs(distance)
		r3 = absr*absr*absr; // ~ten times faster than pow((r1^2 + r2^2 + r3^2),1.5)
		GMr3 = config_data->GM[b] / r3;

		accel[2] += GMr3 * r_body[2];
		accel[0] += GMr3 * r_body[0];
		accel[1] += GMr3 * r_body[1];
		// printf("\n GMr3 * r_body[0]: %.16le", GMr3 * r_body[0]);

		if (config_data->body_int[b] == 10)
		{
			SpiceDouble c = 299792.458	// [km/s]
				, Sn[3]					// [km]
				, rp					// [km/s]
				, GM_r2;				// [km/s^2]

			// Calculate velocity intermediate value
			iVel[0] = Vel[0] + dt * accel[0];
			iVel[1] = Vel[1] + dt * accel[1];
			iVel[2] = Vel[2] + dt * accel[2];

			// normalized body direction value
			Sn[0] = -r_body[0] / absr;
			Sn[1] = -r_body[1] / absr;
			Sn[2] = -r_body[2] / absr;

			rp = iVel[0] * Sn[0] + iVel[1] * Sn[1] + iVel[2] * Sn[2]; // absolute change of the radius between body and particle [km/s]

			// add radiation pressure acceleration and Poynting-Robertson drag acceleration
			GM_r2 = config_data->betaGM / (absr * absr);
			
			accel[0] += GM_r2 * ((1. - SWDF * rp / c) * Sn[0] - SWDF * iVel[0] / c);
			accel[1] += GM_r2 * ((1. - SWDF * rp / c) * Sn[1] - SWDF * iVel[1] / c);
			accel[2] += GM_r2 * ((1. - SWDF * rp / c) * Sn[2] - SWDF * iVel[2] / c);
		}
	}
}



#ifdef __SaveRateOpt
/* Calculate a factor for more saved steps when planets are significantly influencing the acceleration. Only sun > save_factor = ~1, Only planets: save_factor = 0 */
int calc_save_factor(configuration_values *config_data, SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel)
{
	SpiceDouble r_body[3]		// [km]
		, absr = 0.				// [km]
		, r3					// [km^3]
		, GMr3;					// [1/s^2]

	SpiceDouble solAccel[3]		// [km/s^2]
		, save_factor;			// [1]

	int b, noSun = 1; // body

	for (b = 0; b < config_data->N_bodys; b++)
	{
		if (config_data->body_int[b] == 10) // only for the sun
		{
			r_body[0] = (*body_state)[b][0] + dir_SSB[0];
			r_body[1] = (*body_state)[b][1] + dir_SSB[1];
			r_body[2] = (*body_state)[b][2] + dir_SSB[2];

			// Calculate GM*r^3
			absr = sqrt(r_body[0] * r_body[0] + r_body[1] * r_body[1] + r_body[2] * r_body[2]); // abs(distance)
			r3 = absr*absr*absr;
			GMr3 = config_data->GM[b] / r3;

			solAccel[0] = GMr3 * r_body[0];
			solAccel[1] = GMr3 * r_body[1];
			solAccel[2] = GMr3 * r_body[2];

			noSun = 0;
		}
	}

	if (noSun == 1)
	{
		return 1;
	}

	// Save factor: (solar acceleration / total acceleration)^(0.5)
	save_factor = (solAccel[0] * solAccel[0] + solAccel[1] * solAccel[1] + solAccel[2] * solAccel[2]) / (accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]);

	if (save_factor < 0.8)
	{
		//printf("\n save_factor: %.6le, n_opt: %d", save_factor, config_data->n_opt);
		config_data->n_opt = (int)(sqrt(sqrt(save_factor)) * config_data->n);
		if (config_data->n_opt < 1)
		{
			config_data->n_opt = 1;
		}
	}
	else
	{
		config_data->n_opt = config_data->n;
	}

	return 0;
}
#endif // __SaveRateOpt



/* Print the state vector (x,y,z,vx,vy,vz,t) to the given file */
int printpdata(FILE *statefile, SpiceDouble *nstate)
{
	if (isnan(nstate[0]) || isnan(nstate[1]) || isnan(nstate[2]) || isnan(nstate[3]) || isnan(nstate[4]) || isnan(nstate[5]) || isnan(nstate[6]))
	{
		printf("\n\nerror: printing state to file failed, result is not a number. This value will be skipped.");
		return 1;
	}
	else
	{
		fprintf(statefile, "%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\n", (nstate)[0], (nstate)[1], (nstate)[2], (nstate)[3], (nstate)[4], (nstate)[5], (nstate)[6]);
	}

	return 0;
}



/* Calculate particle-constant part of PRD */
int calc_pInfo(configuration_values *config_data)
{
	int b;

	if (config_data->particle_mass > 0.)
	{
		// Manipulate sun mass to simulate solar pressure
		SpiceDouble PI = 3.1415926535897932		// [1]
			, c = 299792458.;					// [m/s]
		config_data->solar_lum = 3.846e26;		// [W]

		config_data->particle_radius = pow((config_data->particle_mass) / ((1.3333) * PI * (config_data->particle_density)), 0.3333);

		for (b = 0; b < config_data->N_bodys; b++)
		{
			if (config_data->body_int[b] == 10)
			{
				config_data->beta = 3. * (config_data->solar_lum) * config_data->q_pr / (16. * PI * c * (config_data->GM[b] * 1.e9) * config_data->particle_density * config_data->particle_radius);

				config_data->betaGM = config_data->beta * config_data->GM[b];
				break;
			}
		}
	}

	return 0;
}



/* Function imitating spkezp_c but always returning the SSB as the body poisition.
   Only used when SSB_CENTERED	=1
   */
void return_SSB(	SpiceInt            targ,
					SpiceDouble         et,
					ConstSpiceChar    * ref,
					ConstSpiceChar    * abcorr,
					SpiceInt            obs,
					SpiceDouble         ptarg[3],
					SpiceDouble       * lt) // Most parameters are unreferenced but necessary for identical function calls to spkezp_c
{
	// not using: targ, et, ref, abcorr, obs, lt
	(void)targ;
	(void)et;
	(void)ref;
	(void)abcorr;
	(void)obs;
	(void)lt;

	int j;
	for (j = 0; j < 3; j++)
	{
		ptarg[j] = (SpiceDouble)0.0;
	}
}



/* Function imitating spkezr_c but always returning the SSB as the body poisition and zero as the body speed.
   Only used when SSB_CENTERED	=1
   */
void return_SSBr(	ConstSpiceChar    * targ,
					SpiceDouble         et,
					ConstSpiceChar    * ref,
					ConstSpiceChar    * abcorr,
					ConstSpiceChar    * obs,
					SpiceDouble         ptarg[6],
					SpiceDouble       * lt) // Most parameters are unreferenced but necessary for identical function calls to spkezp_c
{
	// not using: targ, et, ref, abcorr, obs, lt
	(void)targ;
	(void)et;
	(void)ref;
	(void)abcorr;
	(void)obs;
	(void)lt;

	int j;
	for (j = 0; j < 6; j++)
	{
		ptarg[j] = (SpiceDouble)0.0;
	}
}



/* Allocate memory for variables in interp_body_states() 
   */
int interp_body_states_malloc(configuration_values *config_data, SpiceDouble **(*body_c)[6])
{
	int i, n;

	int vector_size = 0;

	if (config_data->interp_order == 2)
	{
		vector_size = 3;
	}
	else if (config_data->interp_order == 5)
	{
		vector_size = 6;
	}
	else if (config_data->interp_order == 0)
	{
		;
	}
	else
	{
		;
	}

	if (vector_size != 0)
	{
		// Allocate memory for coefficients
		for (n = 0; n < vector_size; n++){
			(*body_c)[n] = malloc(config_data->N_bodys * sizeof(SpiceDouble *));
			if ((*body_c)[n] == NULL)
			{
				printf("\n\nerror: could not allocate body coefficient array (OOM)");
				return 1;
			}
		}
		for (i = 0; i < config_data->N_bodys; i++)
		{
			for (n = 0; n < vector_size; n++){
				(*body_c)[n][i] = malloc(3 * sizeof(SpiceDouble));
				if ((*body_c)[n][i] == NULL)
				{
					printf("\n\nerror: could not allocate body coefficient array (OOM)");
					return 1;
				}
			}
		}
	}

	return 0;
}


/* Free memory for variables in interp_body_states()
*/
int interp_body_states_free(configuration_values *config_data, SpiceDouble **(*body_c)[6])
{
	int i, n;

	int vector_size = 0;

	if (config_data->interp_order == 2)
	{
		vector_size = 3;
	}
	else if (config_data->interp_order == 5)
	{
		vector_size = 6;
	}
	else if (config_data->interp_order == 0)
	{
		;
	}
	else
	{
		;
	}

	if (vector_size != 0)
	{
		// Free memory of coefficients
		for (i = 0; i < config_data->N_bodys; i++)
		{
			for (n = 0; n < vector_size; n++){
				free((*body_c)[n][i]);
			}
		}
		for (n = 0; n < vector_size; n++){
			free((*body_c)[n]);
		}
	}

	return 0;
}



/* Interpolate body states, either 2nd order or 5th order (coming soon), 
   returns 1 if OOM and 2 if interpolation fails. 
   */
int interp_body_states(configuration_values *config_data, SpiceDouble **(*body)[9], SpiceDouble **(*body_c)[6], SpiceDouble dtime[9], SpiceDouble h, int j)
{
	int k, m;

	if (config_data->interp_order == 2)
	{
		// Solving x = a + bt + ct^2 for quadratic interpolation of body positions
		for (k = 0; k < 3; k++) // loop x,y,z
		{
			(*body_c)[0][j][k] = (*body)[0][j][k];

			(*body_c)[2][j][k] = ((((*body)[8][j][k] - (*body_c)[0][j][k]) / dtime[8]) - (((*body)[1][j][k] - (*body_c)[0][j][k]) / dtime[1])) / h;

			(*body_c)[1][j][k] = ((*body)[1][j][k] - (*body_c)[0][j][k]) / dtime[1] - (*body_c)[2][j][k] * dtime[1];

			// Interpolate body states
			for (m = 2; m < 8; m++)
			{
				(*body)[m][j][k] = (*body_c)[0][j][k] + ((*body_c)[1][j][k] + (*body_c)[2][j][k] * dtime[m]) * dtime[m];
			}
		}
	}
	else if (config_data->interp_order == 5)
	{
		// Pre-compute power dtimes
		SpiceDouble dtime1p2 = dtime[1] * dtime[1];
		SpiceDouble dtime1p3 = dtime1p2 * dtime[1];
		SpiceDouble dtime1p4 = dtime1p3 * dtime[1];
		SpiceDouble dtime1p5 = dtime1p4 * dtime[1];
		SpiceDouble dtime1p6 = dtime1p5 * dtime[1];

		SpiceDouble dtime8p2 = dtime[8] * dtime[8];
		SpiceDouble dtime8p3 = dtime8p2 * dtime[8];
		SpiceDouble dtime8p4 = dtime8p3 * dtime[8];
		SpiceDouble dtime8p5 = dtime8p4 * dtime[8];
		SpiceDouble dtime8p6 = dtime8p5 * dtime[8];

		SpiceDouble dtime81 = (dtime[1] - dtime[8]);
		SpiceDouble dtime81p2 = dtime81 * dtime81;

		SpiceDouble sig[11]; // Precomputed variables needed in more than one formula
		sig[0] = dtime[1] * dtime8p3*(dtime1p2*dtime[8] - dtime1p3)*dtime81p2;

		// Solving x = a + bt + ct^2 + dt^3 + et^4 + ft^5 for 5th order interpolation of body positions
		for (k = 0; k < 3; k++) // loop x,y,z
		{
			// Pre-compute a few things that are needed twice
			sig[1] = 4 * dtime1p4*dtime8p2 * (*body)[0][j][k + 3];
			sig[2] = 4 * dtime1p2*dtime8p4 * (*body)[0][j][k + 3];
			sig[3] = dtime1p4*dtime8p2 * (*body)[8][j][k + 3];
			sig[4] = dtime1p2*dtime8p4 * (*body)[1][j][k + 3];
			sig[5] = 2 * dtime1p5*dtime[8] * (*body)[0][j][k + 3];
			sig[6] = 2 * dtime[1] * dtime8p5 * (*body)[0][j][k + 3];
			sig[7] = 5 * dtime1p4*dtime[8] * (*body)[8][j][k];
			sig[8] = 5 * dtime[1] * dtime8p4 * (*body)[1][j][k];
			sig[9] = 5 * dtime1p4*dtime[8] * (*body)[0][j][k];
			sig[10] = 5 * dtime[1] * dtime8p4 * (*body)[0][j][k];

			(*body_c)[0][j][k] = (*body)[0][j][k];

			(*body_c)[1][j][k] = (*body)[0][j][k + 3]; // k+3 gives the speed vector instead of the position vector

			(*body_c)[2][j][k] = -(3 * dtime1p5 * (*body)[0][j][k] - 3 * dtime8p5 * (*body)[0][j][k] - 3 * dtime1p5 * (*body)[8][j][k] + 3 * dtime8p5 * (*body)[1][j][k] - sig[6] + sig[5]
				- dtime[1] * dtime8p5 * (*body)[1][j][k + 3] + dtime1p5*dtime[8] * (*body)[8][j][k + 3] + sig[2] - sig[1] + sig[4] - sig[3]
				+ sig[10] - sig[9] - sig[8] + sig[7]) / (dtime1p2*dtime8p2*dtime81p2*dtime81);

			(*body_c)[3][j][k] = -(2 * dtime1p6 * (*body)[0][j][k] - 2 * dtime8p6 * (*body)[0][j][k] - 2 * dtime1p6 * (*body)[8][j][k] + 2 * dtime8p6 * (*body)[1][j][k] - dtime[1] * dtime8p6 * (*body)[0][j][k + 3] + dtime1p6*dtime[8] * (*body)[0][j][k + 3] - dtime[1] * dtime8p6 * (*body)[1][j][k + 3]
				+ dtime1p6*dtime[8] * (*body)[8][j][k + 3] + 10 * dtime1p2*dtime8p4 * (*body)[0][j][k] - 10 * dtime1p4*dtime8p2 * (*body)[0][j][k] - 10 * dtime1p2*dtime8p4 * (*body)[1][j][k] + 10 * dtime1p4*dtime8p2 * (*body)[8][j][k] - dtime1p2*dtime8p5 * (*body)[0][j][k + 3]
				+ 8 * dtime1p3*dtime8p4 * (*body)[0][j][k + 3] - 8 * dtime1p4*dtime8p3 * (*body)[0][j][k + 3] + dtime1p5*dtime8p2 * (*body)[0][j][k + 3] - dtime1p2*dtime8p5 * (*body)[1][j][k + 3] + 2 * dtime1p3*dtime8p4 * (*body)[1][j][k + 3] - 2 * dtime1p4*dtime8p3 * (*body)[8][j][k + 3]
				+ dtime1p5*dtime8p2 * (*body)[8][j][k + 3] - 2 * dtime[1] * dtime8p5 * (*body)[0][j][k] + 2 * dtime1p5*dtime[8] * (*body)[0][j][k] + 2 * dtime[1] * dtime8p5 * (*body)[1][j][k] - 2 * dtime1p5*dtime[8] * (*body)[8][j][k]) / (sig[0]);

			(*body_c)[4][j][k] = (4 * dtime1p5 * (*body)[0][j][k] - 4 * dtime8p5 * (*body)[0][j][k] - 4 * dtime1p5 * (*body)[8][j][k] + 4 * dtime8p5 * (*body)[1][j][k] - sig[6] + sig[5] - 2 * dtime[1] * dtime8p5 * (*body)[1][j][k + 3]
				+ 2 * dtime1p5*dtime[8] * (*body)[8][j][k + 3] + 5 * dtime1p2*dtime8p3 * (*body)[0][j][k] - 5 * dtime1p3*dtime8p2 * (*body)[0][j][k] - 5 * dtime1p2*dtime8p3 * (*body)[1][j][k] + 5 * dtime1p3*dtime8p2 * (*body)[8][j][k] + sig[2]
				- sig[1] + sig[4] + dtime1p3*dtime8p3 * (*body)[1][j][k + 3] - dtime1p3*dtime8p3 * (*body)[8][j][k + 3] - sig[3] + sig[10] - sig[9]
				- sig[8] + sig[7]) / (sig[0]);

			(*body_c)[5][j][k] = (2 * dtime1p4 * (*body)[0][j][k] - 2 * dtime8p4 * (*body)[0][j][k] - 2 * dtime1p4 * (*body)[8][j][k]
				+ 2 * dtime8p4 * (*body)[1][j][k] - dtime[1] * dtime8p4 * (*body)[0][j][k + 3] + dtime1p4*dtime[8] * (*body)[0][j][k + 3] - dtime[1] * dtime8p4 * (*body)[1][j][k + 3] + dtime1p4*dtime[8] * (*body)[8][j][k + 3] + 3 * dtime1p2*dtime8p3 * (*body)[0][j][k + 3]
				- 3 * dtime1p3*dtime8p2 * (*body)[0][j][k + 3] + dtime1p2*dtime8p3 * (*body)[1][j][k + 3] - dtime1p3*dtime8p2 * (*body)[8][j][k + 3] + 4 * dtime[1] * dtime8p3 * (*body)[0][j][k] - 4 * dtime1p3*dtime[8] * (*body)[0][j][k]
				- 4 * dtime[1] * dtime8p3 * (*body)[1][j][k] + 4 * dtime1p3*dtime[8] * (*body)[8][j][k]) / (dtime1p3*dtime8p3*dtime81*(dtime1p2 - 2 * dtime[1] * dtime[8] + dtime8p2));

			// Interpolate body states
			for (m = 2; m < 8; m++)
			{
				(*body)[m][j][k] = (*body_c)[0][j][k] + ((*body_c)[1][j][k] + ((*body_c)[2][j][k] + ((*body_c)[3][j][k] + ((*body_c)[4][j][k] + (*body_c)[5][j][k] * dtime[m]) * dtime[m]) * dtime[m]) * dtime[m]) * dtime[m];
			}
		}
	}
	else if (config_data->interp_order == 0)
	{
		return 2;
	}
	else
	{
		printf("\n\nwarning: interpolation order not supported: %d", config_data->interp_order);
		return 2;
	}

	return 0;
}