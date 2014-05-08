// Functions called by integration algorithms

void calc_accel(int N_bodys, SpiceDouble GM[], SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel)
{
	SpiceDouble direct_body[3], r3, GMr3, absr;
	int b; // body

	accel[0] = 0;
	accel[1] = 0;
	accel[2] = 0;

	for (b = 0; b < N_bodys; b++)
	{
		direct_body[0] = (*body_state)[b][0] + dir_SSB[0];
		direct_body[1] = (*body_state)[b][1] + dir_SSB[1];
		direct_body[2] = (*body_state)[b][2] + dir_SSB[2];
		absr = sqrt(direct_body[0] * direct_body[0] + direct_body[1] * direct_body[1] + direct_body[2] * direct_body[2]); // abs(distance)
		r3 = absr*absr*absr; // ~ten times faster than pow((r1^2 + r2^2 + r3^2),1.5)
		GMr3 = GM[b] / r3;
		accel[2] += GMr3 * direct_body[2];
		accel[0] += GMr3 * direct_body[0];
		accel[1] += GMr3 * direct_body[1];
	}
}



int printpdata(FILE *statefile, SpiceDouble *nstate)
{
	fprintf(statefile, "%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\t%.16le\n", (nstate)[0], (nstate)[1], (nstate)[2], (nstate)[3], (nstate)[4], (nstate)[5], (nstate)[6]);

	return 0;
}