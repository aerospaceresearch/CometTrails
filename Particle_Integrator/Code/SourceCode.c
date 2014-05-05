#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
//#include <time.h>
#include <direct.h>
#include <omp.h>
#include <SpiceUsr.h>


void calc_accel(int N_bodys, SpiceDouble GM[], SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel);
int RungeKutta4(int N_bodys, int body_int[], SpiceDouble GM[], SpiceDouble final_time, SpiceDouble start_time, SpiceDouble start_time_save, SpiceDouble dv_step, SpiceDouble *nstate, char outputpath[], int n, int particles_count);
bool particle_already_processed(int p, char already_done_path[]);
bool particle_incomplete(char outputpath[], SpiceDouble *nstate);
int read_configuration(char *inputpath, char *outputpath, int *number_of_threads, SpiceDouble *final_time,SpiceDouble *start_time_final, int *N_bodys, int *body_int, SpiceDouble *GM, SpiceDouble *dv_step, int *n, int *first_particle_number, SpiceDouble *particle_mass, SpiceDouble *particle_density);


//Main Programm
int main(void)
{
	//Print version
	char  * versn;
	versn = "0.12";
	printf("ParticleIntegrator version %s\n", versn);

	//Create some variables
	int j, e, p, g, c, error_code = 0, particles_count = 0, particles_done = 0, number_of_threads = 0, first_particle_number = 0, n, N_bodys, body_int[10];
	SpiceDouble final_time = 0, start_time_save = 0, dv_step = 0, particle_mass = 0, particle_density = 0, GM[10];
	char temp[200], *next_token = NULL, inputpath[200] = "INPUT\\", outputpath[200] = "OUTPUT\\", already_done_path[200] = "INPUT\\processed_particles.txt";
	bool commentLine = false;

	//Load Spice kernels
	printf("\nLoading kernels...		");
	furnsh_c("kernels_generic.txt");
	printf("...done.");

	//READ CONFIG FILE
	printf("\nLoading configuration...	");
	if (read_configuration(inputpath, outputpath, &number_of_threads, &final_time, &start_time_save, &N_bodys, body_int, GM, &dv_step, &n, &first_particle_number, &particle_mass, &particle_density) != 0)
	{
		printf("\n\nerror:	could not read configuration.\n");
		//Sleep(4000);
		return 1;
	}
	printf("...done.");

	//LOAD PARTICLES
	printf("\nLoading particles...		");
	FILE *particles_start_file;
	SpiceDouble **particles_start;
	fopen_s(&particles_start_file, inputpath, "r");
	if (particles_start_file == NULL)
	{
		printf("\n\nerror:	could not load particles.\n");
		//Sleep(4000);
		return 1;
	}
	while ((c = fgetc(particles_start_file)) != EOF) //#tmp# this crashes if the last line is empty! the implementation relies on the number of lines being number of particles + 1 (+ 5 ocmment lines). > not nice
	{
		if (c == '%')
		{
			commentLine = true;
		}
		if (c == '\n')
		{
			if (commentLine == false)
			{
				particles_count++;
			}
			else
			{
				commentLine = false;
			}
		}
	}
	particles_start = malloc((particles_count + 1) * sizeof(int *));
	for (j = 0; j < particles_count; j++)
	{
		particles_start[j] = malloc(8 * sizeof(SpiceDouble));
	}
	fclose(particles_start_file);
	fopen_s(&particles_start_file, inputpath, "r");
	j = 0;
	while (fgets(temp, sizeof(temp), particles_start_file) != NULL)
	{
		if (j > 4)
		{
			char* cval = strtok_s(temp, "\t", &next_token);
			for (g = 0; g < 6; g++)
			{
				sscanf_s(cval, "%lf", &particles_start[j - 5][g]);
				cval = strtok_s(NULL, "\t", &next_token);
			}
			sscanf_s(cval, "%lf", &particles_start[j - 5][6]);
			cval = strtok_s(NULL, "\n", &next_token);
			sscanf_s(cval, "%lf", &particles_start[j - 5][7]);
		}
		j++;
	}
	int last_particle_number = first_particle_number + particles_count - 1;
	fclose(particles_start_file);
	printf("...done. %d particles loaded.\n", particles_count);


	//Print config
	if (number_of_threads > 1)
		printf("\n number of threads	= %d", number_of_threads);
	printf("\n final_time		= %le", final_time);
	if (start_time_save > (double) -3.155e+10)
		printf("\n start_time_save	= %le", start_time_save);
	printf("\n bodys_ID		=");
	for (j = 0; j < N_bodys; j++)
		printf(" %d", body_int[j]);
	printf("\n dv_step		= %le", dv_step);
	if (particle_mass > 0)
		printf("\n particle_mass		= %le", particle_mass);
	printf("\n particle_density	= %le", particle_density);
	if (first_particle_number != 1)
		printf("\n first_particle_number	= %d", first_particle_number);
	printf("\n save_nth		= %d", n);
	
	//Check for progress.txt
	FILE *progress, *already_done;
	fopen_s(&progress, "progress.txt", "r+");
	if (progress == NULL)
	{
		for (e = 0; e < 3; e++)
		{
			fopen_s(&progress, "progress.txt", "w");
			if (progress == NULL)
			{
				Sleep(100);
				if (e == 2)
				{
					printf("\n\nerror: could not create progress.txt");
					return 2;
				}
			}
			else
			{
				fprintf(progress, "0.0");
				break;
			}
		}
	}
	fclose(progress);
	//Check for processed_particles.txt and count particles already processed
	fopen_s(&already_done, already_done_path, "r+");
	if (already_done == NULL)
	{
		for (e = 0; e < 3; e++)
		{
			fopen_s(&already_done, already_done_path, "w");
			if (already_done == NULL)
			{
				Sleep(100);
				if (e == 2)
				{
					printf("\n\nerror: could not create processed_particles.txt");
					return 2;
				}
			}
			else
			{
				break;
			}
		}
	}
	else
	{
		//Count the number of particles already done
		while ((c = fgetc(already_done)) != EOF)
		{
			if (c == '\n')
			{
				particles_done++;
			}
		}
	}
	fclose(already_done);

	
	//Start (parallel) computing
	//clock_t start = clock();
	printf("\n\n Numerical approximation started...");

	/*--------------------------------------------------------------------------------------------------------------------------*/
#pragma omp parallel private(j, e) num_threads(number_of_threads)
	{
		int th_id = omp_get_thread_num();
		
		//Allocate nstate
		SpiceDouble nstate[7];

		//Loop over particles
#pragma omp for
		for (p = first_particle_number; p <= last_particle_number; p++)
		{
			int err = 0;

			//Check if particle has already been processed completely
			if (particle_already_processed(p, already_done_path))
			{
				printf("\n particle #%d	has already been processed", p);
				continue;
			}

			//Set path where to save the particle
			char particle_num[9];
			char particle_path[200] = "";
			sprintf_s(particle_num, 9, "_#%d", p);
			strcat_s(particle_path, 200, outputpath);
			strcat_s(particle_path, 200, particle_num);
			strcat_s(particle_path, 200, ".txt");

			//Set particle start_time
			SpiceDouble start_time = particles_start[p - first_particle_number][7];

			//Check if particle has been processed but is incomplete
			if (particle_incomplete(particle_path, nstate))
			{
				printf("\n particle #%d	will be continued", p);
			}
			else
			{
				//Set initial nstate
				for (j = 0; j < 6; j++)
					nstate[j] = particles_start[p - first_particle_number][j];
				nstate[6] = start_time;
				//Create  initial file
				FILE* init;
				for (e = 0; e < 3; e++)
				{
					fopen_s(&init, particle_path, "w");
					if (init == NULL)
					{
						Sleep(100);
						if (e == 2)
						{
							err = 1;
							printf("\nerror: could not create initial output file");
							break;
						}
					}
					else
					{
						char temp[200] = "", tempval[30] = "";
						for (j = 0; j < 6; j++)
						{
							sprintf_s(tempval, 30, "%.16le\t", (nstate)[j]);
							strcat_s(temp, 200, tempval);
						}
						sprintf_s(tempval, 30, "%.16le\n", (nstate)[6]);
						strcat_s(temp, 200, tempval);
						fprintf(init, "%s", temp);
						fclose(init);
						break;
					}
				}
			}


			//Integrate particle
			if (err == 0)
			{
				err = RungeKutta4(N_bodys, body_int, GM, final_time, start_time, start_time_save, dv_step, nstate, particle_path, n, particles_count);
			}
			

			//Write the particle number to the already-done file and update progress.txt
			FILE* done;
			FILE* progress;
			double fraction;
#pragma omp critical(ALREADYDONE)
			{
				particles_done++;
				for (e = 0; e < 3; e++)
				{
					fopen_s(&done, already_done_path, "a+");
					if (done == NULL)
					{
						Sleep(100);
						if (e == 2)
						{
							err = 1;
							printf("\nerror: could not write to processed_particles.txt");
							break;
						}
					}
					else
					{
						fprintf(done, "%d\n", p);
						fclose(done);
						break;
					}
				}
			}
#pragma omp critical(PROCESS)
			{
				for (e = 0; e < 3; e++)
				{
					fopen_s(&progress, "progress.txt", "w");
					if (progress == NULL)
					{
						Sleep(100);
						if (e == 2)
						{
							//err = 1;
							printf("\nerror: could not write to progress.txt	(non-relevant)");
							break;
						}
					}
					else
					{
						fraction = (float) particles_done / particles_count;
						fprintf(progress, "%f", fraction);
						fclose(progress);
						break;
					}
				}
			}
			if (err != 0)
			{
				error_code += err;
				printf("\n particle #%d	was not succesfully completed", p);
			}
			else
			{
				printf("\n particle #%d	done on thread %d", p, th_id);
			}
		}

		
		//Print elapsed time
		/*
#pragma omp barrier
#pragma omp master
		{
			//remove(already_done_path);
			//Print time
			//clock_t end = clock();
			//double elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;
			//printf("\n\n Elapsed time: %1.3f s", elapsed_time);
		}*/
	}

	
	//Deallocate arrays
	for (j = 0; j < particles_count; j++)
		free(particles_start[j]);
	free(particles_start);
	
	
	if (error_code == 0)
	{
		printf("\n\nAll particles are done. Everything is OK!\n");
		return 0;
	}
	else
	{
		printf("\n\nWarning: %d	particles may have been skipped!\n", error_code);
		return 2;
	}
	//Sleep(3000);
}





//Functions

int RungeKutta4(int N_bodys, int body_int[], SpiceDouble GM[], SpiceDouble final_time, SpiceDouble start_time, SpiceDouble start_time_save, SpiceDouble dv_step, SpiceDouble *nstate, char particle_path[], int n, int particles_count)
{
	//Create some variables
	int j, i = 0, err = 0;
	SpiceDouble lt, dt, dt2;
	
	//Create File for output
	FILE *statefile;
	fopen_s(&statefile, particle_path, "a");
	if (statefile == NULL)
	{
		printf("\nerror: could not write to output file");
		return 1;
	}

	//Create body arrays and set initial body positions
	SpiceDouble **body_pre, **body_mid, **body_end;
	body_pre = malloc(N_bodys * sizeof(int *));
	body_mid = malloc(N_bodys * sizeof(int *));
	body_end = malloc(N_bodys * sizeof(int *));
	if (body_pre == NULL || body_mid == NULL || body_end == NULL)
	{
		printf("\nerror: could not allocate body state array (OOM)");
		return 1;
	}
	for (j = 0; j < N_bodys; j++)
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
			spkezp_c(body_int[j], nstate[6], "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
		}
	}

	//Create some more variables
	SpiceDouble initPos[3], initVel[3], dir_SSB[3], initTime, abs_acc;
	SpiceDouble k_acc_1[3], k_acc_2[3], k_acc_3[3], k_acc_4[3];
	SpiceDouble k_vel_1[3], k_vel_2[3], k_vel_3[3], k_vel_4[3];
	
	//Integrate
	while (nstate[6] < final_time)
	{

		//Set initial state for this step
		initPos[0] = nstate[0];
		initPos[1] = nstate[1];
		initPos[2] = nstate[2];
		initVel[0] = nstate[3];
		initVel[1] = nstate[4];
		initVel[2] = nstate[5];
		initTime   = nstate[6];
		for (j = 0; j < N_bodys; j++)
		{
			body_pre[j][0] = body_end[j][0];
			body_pre[j][1] = body_end[j][1];
			body_pre[j][2] = body_end[j][2];
		}

		//Step 1
		dir_SSB[0] = -initPos[0];
		dir_SSB[1] = -initPos[1];
		dir_SSB[2] = -initPos[2];
		calc_accel(N_bodys, GM, dir_SSB, &body_pre, k_acc_1);
		k_vel_1[0] = initVel[0];
		k_vel_1[1] = initVel[1];
		k_vel_1[2] = initVel[2];

		//Set dynamic step size
		abs_acc = sqrt(k_acc_1[0] * k_acc_1[0] + k_acc_1[1] * k_acc_1[1] + k_acc_1[2] * k_acc_1[2]);
		dt = (dv_step / abs_acc);
		dt2 = dt / 2;

		
		/*if ((*nstate)[i][6] - dt < final_time)
		{
			dt = fabs(final_time - (*nstate)[i][6]);
		}*/
		
		//Get body positions with SPICE
		for (j = 0; j < N_bodys; j++)
		{
#pragma omp critical(SPICE)
			{
				//Critical section is only executed on one thread at a time (spice is not threadsafe)
				spkezp_c(body_int[j], initTime + dt, "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
			}
			body_mid[j][0] = (body_pre[j][0] + body_end[j][0]) / 2;
			body_mid[j][1] = (body_pre[j][1] + body_end[j][1]) / 2;
			body_mid[j][2] = (body_pre[j][2] + body_end[j][2]) / 2;
		}
		
		//Step 2
		dir_SSB[0] = -(initPos[0] + k_vel_1[0] * dt2);
		dir_SSB[1] = -(initPos[1] + k_vel_1[1] * dt2);
		dir_SSB[2] = -(initPos[2] + k_vel_1[2] * dt2);
		calc_accel(N_bodys, GM, dir_SSB, &body_mid, k_acc_2);
		k_vel_2[0] = initVel[0] + k_acc_1[0] * dt2;
		k_vel_2[1] = initVel[1] + k_acc_1[1] * dt2;
		k_vel_2[2] = initVel[2] + k_acc_1[2] * dt2;

		//Step 3
		dir_SSB[0] = -(initPos[0] + k_vel_2[0] * dt2);
		dir_SSB[1] = -(initPos[1] + k_vel_2[1] * dt2);
		dir_SSB[2] = -(initPos[2] + k_vel_2[2] * dt2);
		calc_accel(N_bodys, GM, dir_SSB, &body_mid, k_acc_3);
		k_vel_3[0] = initVel[0] + k_acc_2[0] * dt2;
		k_vel_3[1] = initVel[1] + k_acc_2[1] * dt2;
		k_vel_3[2] = initVel[2] + k_acc_2[2] * dt2;

		//Step 4
		dir_SSB[0] = -(initPos[0] + k_vel_3[0] * dt);
		dir_SSB[1] = -(initPos[1] + k_vel_3[1] * dt);
		dir_SSB[2] = -(initPos[2] + k_vel_3[2] * dt);
		calc_accel(N_bodys, GM, dir_SSB, &body_end, k_acc_4);
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
		if (i == n)
		{
			if (nstate[6] > start_time_save)
			{
				char temp[200] = "", tempval[30] = "";
				for (j = 0; j < 6; j++)
				{
					sprintf_s(tempval, 30, "%.16le\t", nstate[j]);
					strcat_s(temp, 200, tempval);
				}
				sprintf_s(tempval, 30, "%.16le\n", nstate[6]);
				strcat_s(temp, 200, tempval);
				fprintf(statefile, "%s", temp);
			}
			i = 0;
		}
	}

	//Print last state to file and close file
	if (i != 0)
	{
		char temp[200] = "", tempval[30] = "";
		for (j = 0; j < 6; j++)
		{
			sprintf_s(tempval, 30, "%.16le\t", nstate[j]);
			strcat_s(temp, 200, tempval);
		}
		sprintf_s(tempval, 30, "%.16le\n", nstate[6]);
		strcat_s(temp, 200, tempval);
		fprintf(statefile, "%s", temp);
	}
	fclose(statefile);

	
	//Deallocate body array
	for (j = 0; j < N_bodys; j++)
	{
		free(body_pre[j]);
		free(body_mid[j]);
		free(body_end[j]);
	}
	free(body_pre);
	free(body_mid);
	free(body_end);

	return 0;
}



void calc_accel(int N_bodys, SpiceDouble GM[], SpiceDouble dir_SSB[], SpiceDouble **body_state[], SpiceDouble *accel)
{
	SpiceDouble direct_body[3], distance_pow3, GMr3;
	accel[0] = 0;
	accel[1] = 0;
	accel[2] = 0;
	
	int b;
	for (b = 0; b < N_bodys; b++)
	{
		direct_body[0] = (*body_state)[b][0] + dir_SSB[0];
		direct_body[1] = (*body_state)[b][1] + dir_SSB[1];
		direct_body[2] = (*body_state)[b][2] + dir_SSB[2];
		distance_pow3 = pow(direct_body[0] * direct_body[0] + direct_body[1] * direct_body[1] + direct_body[2] * direct_body[2], 1.5);
		GMr3 = GM[b] / distance_pow3;
		accel[0] += GMr3 * direct_body[0];
		accel[1] += GMr3 * direct_body[1];
		accel[2] += GMr3 * direct_body[2];
	}
}



bool particle_already_processed(int p, char already_done_path[])
{
	FILE* check;
	char temp[6] = "";
	int particle_ID;
	bool answer = false;
#pragma omp critical(ALREADYDONE)
	{
		fopen_s(&check, already_done_path, "r");
		if (check != NULL)
		{
			while (fgets(temp, 6, check) != NULL)
			{
				sscanf_s(temp, "%d", &particle_ID);
				if (particle_ID == p)
				{
					answer = true;
				}
			}
			fclose(check);
		}
	}

	return answer;
}



bool particle_incomplete(char particle_path[], SpiceDouble *nstate)
{
	int err = 0;
	FILE *check, *tempfile;
	char temp1[200], temp2[200], temp3[200], *next_token = NULL;
	bool answer = false;
	fopen_s(&check, particle_path, "r");
	if (check != NULL)
	{
#pragma omp critical(TEMP)
		{
			fopen_s(&tempfile, "temp.txt", "w");
			if (tempfile == NULL)
			{
				printf("\nerror: could not create tempfile; particle restarted");
				err = 1;
			}
			else
			{
				int c = 0;
				while (fgets(temp1, sizeof(temp1), check) != NULL)
				{
					if (c > 0)
					{
						strcpy_s(temp3, 200, temp2);
					}
					strcpy_s(temp2, 200, temp1);
					if (c > 0)
					{
						fprintf(tempfile, "%s", temp3);
					}
					c++;
				}
				fclose(tempfile);
				fclose(check);

				if (c > 1)
				{
					char* cval = strtok_s(temp3, "\t", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[0]);
					cval = strtok_s(NULL, "\t", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[1]);
					cval = strtok_s(NULL, "\t", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[2]);
					cval = strtok_s(NULL, "\t", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[3]);
					cval = strtok_s(NULL, "\t", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[4]);
					cval = strtok_s(NULL, "\t", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[5]);
					cval = strtok_s(NULL, "\n", &next_token);
					sscanf_s(cval, "%lf", &(nstate)[6]);
					answer = true;

					//Write the output file again without the last line
					fopen_s(&tempfile, "temp.txt", "r");
					if (tempfile == NULL)
					{
						printf("\nerror: could not create tempfile; particle restarted");
						err = 1;
					}
					else
					{
						fopen_s(&check, particle_path, "w");
						if (check == NULL)
						{
							printf("\nerror: could not read incomplete outputfile; particle restarted");
							err = 1;
						}
						else
						{
							while (fgets(temp1, sizeof(temp1), tempfile) != NULL)
							{
								fprintf(check, "%s", temp1);
							}
							fclose(check);
						}
						fclose(tempfile);
					}
				}
			}
			remove("temp.txt");
		}
		if (err != 0)
		{
			return false;
		}
	}
	return answer;
}

int read_configuration(char *inputpath, char *outputpath, int *number_of_threads, SpiceDouble *final_time, SpiceDouble *start_time_save, int *N_bodys, int *body_int, SpiceDouble *GM, SpiceDouble *dv_step, int *n, int *first_particle_number, SpiceDouble *particle_mass, SpiceDouble *particle_density)
{
	char temp[200], *token, *next_token = NULL, configpath[200] = "";
	int dim, j;
	strcat_s(configpath, 200, inputpath);
	strcat_s(configpath, 200, "configuration.ini");
	
	//Set number of threads (Default is 1)
    *number_of_threads = GetPrivateProfileIntA("simulation", "NUMBER_OF_THREADS", 1, configpath);

	//Set final date of the simulation
    GetPrivateProfileStringA("simulation", "FINAL_TIME", "", temp, sizeof(temp), configpath);
	if (strcmp(temp, "") == 0)
	{
		printf("\n\nerror:	FINAL_TIME not set");
		Sleep(1000);
		return 1;
	}
	str2et_c(temp, final_time);

	//Set start date for saving
    GetPrivateProfileStringA("simulation", "START_TIME_SAVE", "", temp, sizeof(temp), configpath);
	if (strcmp(temp, "") == 0)
	{
		sprintf_s(temp, 200, "1 JAN 1000");
	}
	str2et_c(temp, start_time_save);
	
	//Set bodys
    *N_bodys = GetPrivateProfileIntA("simulation", "N_BODYS", 0, configpath);
	if (*N_bodys == 0)
	{
		printf("\n\nerror:	N_BODYS not set");
		Sleep(1000);
		return 1;
	}
	
    GetPrivateProfileStringA("simulation", "BODYS_ID", "", temp, sizeof(temp), configpath);
	token = strtok_s(temp, " ", &next_token);
	
	for (j = 0; j < *N_bodys; j++)
	{
		if (token == NULL)
		{
			printf("\n\nerror:	BODYS_ID not set or not enough arguments");
			Sleep(1000);
			return 1;
		}
		sscanf_s(token, "%d", &body_int[j]);
		token = strtok_s(NULL, " ", &next_token);
		bodvar_c(body_int[j], "GM", &dim, &GM[j]);
	}

	//Set step size control
    GetPrivateProfileStringA("simulation", "DV_STEP", "10e-5", temp, sizeof(temp), configpath);
	sscanf_s(temp, "%lf", dv_step);
	
	//Set which nth state is saved to disc
    int nf = GetPrivateProfileIntA("simulation", "SAVE_NTH_MULTIPLIER", 20, configpath);
	if (nf == 0)
	{
		//Save every 10nth state. This produces high density of states in the output file and is intended to be used when testing the integrator.
		*n = 10;
	}
	else
	{
		*n = (int)(nf * pow(*dv_step, -1) + 0.5);
	}

	//Set which particle to start and end with (particle number, from 1 to the number of particles in the input file)
    *first_particle_number = GetPrivateProfileIntA("particles", "FIRST_PARTICLE_NUMBER", 1, configpath);
	
	//Set name of the input/output file
    GetPrivateProfileStringA("particles", "PARTICLE_INPUT_FILE_NAME", "", temp, sizeof(temp), configpath);
	if (strcmp(temp, "") == 0)
	{
		printf("\n\nerror:	PARTICLE_INPUT_FILE_NAME not set");
		Sleep(1000);
		return 1;
	}
	strcat_s(inputpath, 200, temp);
	strcat_s(inputpath, 200, ".txt");
    GetPrivateProfileStringA("particles", "PARTICLE_OUTPUT_FILE_NAME", temp, temp, sizeof(temp), configpath);
	_mkdir(outputpath);
	strcat_s(outputpath, 200, temp);

	//Set mass of particles
    GetPrivateProfileStringA("particles", "PARTICLE_MASS", "0", temp, sizeof(temp), configpath);
	sscanf_s(temp, "%lf", particle_mass);
	if (*particle_mass > 0)
	{
		//Set density of particles
        GetPrivateProfileStringA("particles", "PARTICLE_DENSITY", "1000", temp, sizeof(temp), configpath);
		sscanf_s(temp, "%lf", particle_density);

		//Manipulate sun mass to simulate solar pressure
		SpiceDouble Qpr, PI, particle_radius, beta;
		Qpr = 1.0;
		PI = 3.14159;
		particle_radius = pow( (*particle_mass) / ((1.3333) * PI * (*particle_density)), 0.3333 );
		beta = 5.7e-4 * (Qpr / ( (*particle_density) * particle_radius) );
		for (j = 0; j < *N_bodys; j++)
		{
			if (body_int[j] == 10)
			{
				GM[j] = GM[j] * (1 - beta);
				break;
			}
		}
	}

	return 0;
}
