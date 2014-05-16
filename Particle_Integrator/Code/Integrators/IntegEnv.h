// Functions called by integration algorithms

/* Calculate the acceleration of a particle based on the position of the body relative to the SSB */
void calc_accel(configuration_values *config_data, SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel, SpiceDouble *Vel, SpiceDouble dt)
{
	SpiceDouble r_body[3]	// [km]
		, absr = 0.			// [km]
		, r3				// [km^3]
		, GMr3				// [1/s^2]
		, iVel[3]			// Intermediate Speed [km/s]
		, absiV;			// Intermediate Speed [km/s]
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
				, aPRDv[3]				// [km/s^2]
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

			accel[0] += GM_r2 * ((1. - rp / c) * Sn[0] - iVel[0] / c);
			accel[1] += GM_r2 * ((1. - rp / c) * Sn[1] - iVel[1] / c);
			accel[2] += GM_r2 * ((1. - rp / c) * Sn[2] - iVel[2] / c);
		}
	}
}



/* Print the state vector (x,y,z,vx,vy,vz,t) to the given file */
int printpdata(FILE *statefile, SpiceDouble *nstate)
{
	fprintf(statefile, "%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\n", (nstate)[0], (nstate)[1], (nstate)[2], (nstate)[3], (nstate)[4], (nstate)[5], (nstate)[6]);

	return 0;
}



/* Calculate particle-constant part of PRD */
SpiceDouble calc_pInfo(configuration_values *config_data)
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
void return_SSB(SpiceInt            targ,
				SpiceDouble         et,
				ConstSpiceChar    * ref,
				ConstSpiceChar    * abcorr,
				SpiceInt            obs,
				SpiceDouble         ptarg[3],
				SpiceDouble       * lt) // Most parameters are unreferenced but necessary for identical function calls to spkezp_c
{
	// not using: targ, et, ref, abcorr, obs, lt
	int j;
	for (j = 0; j < 3; j++)
	{
		ptarg[j] = (SpiceDouble)0.0;
	}
}