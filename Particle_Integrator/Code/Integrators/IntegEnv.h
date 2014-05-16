// Functions called by integration algorithms

/* Calculate the acceleration of a particle based on the position of the body relative to the SSB */
void calc_accel(configuration_values *config_out, SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel, SpiceDouble *Vel, SpiceDouble PRDconst, SpiceDouble dt)
{
	/* Units:	lSol	[kg*km^2/s^3] = [W*1e6]
	 *			GM		[km^3/s^2]
	 *			dir_SSB	[km]
	 *			p_rad.	[m]
	 *			c2		[km/s]
	 *			absr	[km]
	 *			ansV	[km/s]
	 */

	SpiceDouble direct_body[3], r3, GMr3, absr = 0.;
	SpiceDouble aPRD; // absolute PRD-based acceleration
	SpiceDouble iVel[3], absiV; // Intermediate Speed
	int b; // body

	accel[0] = 0;
	accel[1] = 0;
	accel[2] = 0;

	for (b = 0; b < config_out->N_bodys; b++)
	{
		direct_body[0] = (*body_state)[b][0] + dir_SSB[0];
		direct_body[1] = (*body_state)[b][1] + dir_SSB[1];
		direct_body[2] = (*body_state)[b][2] + dir_SSB[2];

		// Calculate GM*r^3
		absr = sqrt(direct_body[0] * direct_body[0] + direct_body[1] * direct_body[1] + direct_body[2] * direct_body[2]); // abs(distance)
		r3 = absr*absr*absr; // ~ten times faster than pow((r1^2 + r2^2 + r3^2),1.5)
		GMr3 = config_out->GM[b] / r3;

		accel[2] += GMr3 * direct_body[2];
		accel[0] += GMr3 * direct_body[0];
		accel[1] += GMr3 * direct_body[1];
		// printf("\n GMr3 * direct_body[0]: %.16le", GMr3 * direct_body[0]);
	}

#ifdef __PRD
	for (b = 0; b < config_out->N_bodys; b++)
	{
		// Sun: Calculate PRD
		if (config_out->body_int[b] == 10)
		{
			/* 
			Calculate absolute acceleration due to Poynting-Robertson effect
			
			PRDconst = lSol * (particle_radius / 1000 * particle_radius / 1000) / (4. * c2) * sqrt(GM[b]) / particle_mass;
			*/
			aPRD = PRDconst / sqrt(absr*absr*absr*absr*absr);

			/* Calculate velocity intermediate value */
			iVel[0] = Vel[0] + dt * accel[0];
			iVel[1] = Vel[1] + dt * accel[1];
			iVel[2] = Vel[2] + dt * accel[2];

			absiV = sqrt(iVel[0] * iVel[0] + iVel[1] * iVel[1] + iVel[2] * iVel[2]);

			accel[0] += aPRD * -iVel[0] / absiV;
			accel[1] += aPRD * -iVel[1] / absiV;
			accel[2] += aPRD * -iVel[2] / absiV;
			// printf("\n apr * -Vel[0] / absv: %.16le", aPRD * -Vel[0] / absv);
		}
	}
#endif // __PRD
}



/* Print the state vector (x,y,z,vx,vy,vz,t) to the given file */
int printpdata(FILE *statefile, SpiceDouble *nstate)
{
	fprintf(statefile, "%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\n", (nstate)[0], (nstate)[1], (nstate)[2], (nstate)[3], (nstate)[4], (nstate)[5], (nstate)[6]);

	return 0;
}



/* Calculate particle-constant part of PRD */
SpiceDouble calc_prdc(configuration_values *config_out)
{
#ifdef __PRD
	int j;
	SpiceDouble cp2 = 89875517873.681764; // speed of light squared, Unit: [km^2/s^2]
	SpiceDouble lSol = 3.846e20; // Unit: [kg*km^2/s^3] - equivalent to 3.846e26 [W]
	SpiceDouble PRDconst = 0.; // Unit: [km^3.5/s^2]

	for (j = 0; j < config_out->N_bodys; j++)
	{
		if (config_out->body_int[j] == 10)
		{
			PRDconst = lSol * (config_out->particle_radius / 1000 * config_out->particle_radius / 1000) / (4. * cp2) * sqrt(config_out->GM[j]) / config_out->particle_mass;
			break;
		}
	}

	return PRDconst;
#else // not __PRD
	return 0.0;
#endif // __PRD
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