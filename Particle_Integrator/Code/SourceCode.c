#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <SpiceUsr.h>
#include <ini.h>
#include <PIConfig.h>
#ifdef _WIN32
	#include <windows.h> // only needed for Sleep()
	#include <direct.h> // only needed for _mkdir()
#else
	#include <unistd.h> // only needed for usleep()
	#include <sys/stat.h> // only needed for mkdir()
#endif

// activate timing: Add preprocessor definition "__WTIMING"
#ifdef __WTIMING
	#include <time.h>
#endif

// Function cross-platform compatibility
#ifdef _WIN32
	#define SLEEP( a1 ) Sleep( a1 )
	#define mkdir( a1, a2 ) _mkdir( a1 )
#else
	#define SLEEP( a1 ) usleep( a1 * 1000 )
#endif

#ifdef _WIN32 // Avoids MSVC level 3 warning C4996
	#define strdup _strdup
	#define strcpy strcpy_s
	#define sscanf sscanf_s
	#define strtok_r strtok_s
#else
	#define strcpy( a1, a2, a3 ) strcpy( a1, a3 )
	#define fopen_s( a1, a2, a3 ) *a1 = fopen( a2, a3 )
	#define sprintf_s snprintf
#endif

#ifdef _WIN32 // platform-specific separators in paths
	#define OS_SEP "\\"
#else
	#define OS_SEP "/"
#endif
// END Function cross-platform compatibility


#include <RungeKutta4.h>
bool particle_already_processed(int p, char already_done_path[]);
bool particle_incomplete(char outputpath[], SpiceDouble *nstate);
int read_configuration(char *inputfpath, char *outputpath, int *number_of_threads, SpiceDouble *final_time, SpiceDouble *start_time_final, int *N_bodys, int *body_int, SpiceDouble *GM, SpiceDouble *dv_step, int *n, int *first_particle_number, SpiceDouble *particle_mass, SpiceDouble *particle_density, int *save_as_binary);


//Main Program
int main(void)
{
	//Print version
	printf("ParticleIntegrator version " PI_VERSION_MAJOR "." PI_VERSION_MINOR "\n");

	//Create some variables
	int j, e, p, g, c, error_code = 0, particles_count = 0, particles_done = 0, number_of_threads = 0, first_particle_number = 0, n, N_bodys, body_int[10], nCommentLines = 0, save_as_binary = 0;
	SpiceDouble final_time = 0, start_time_save = 0, dv_step = 0, particle_mass = 0, particle_density = 0, GM[10];
	char temp[260], *next_token = NULL, inputfpath[260] = "", outputpath[260] = ("OUTPUT" OS_SEP), already_done_path[260] = "INPUT" OS_SEP "processed_particles.txt";
	bool commentLine = false;

	//Load Spice kernels
	printf("\nLoading kernels...		");
	furnsh_c("kernels_generic.txt");
	printf("...done.");

	//READ CONFIG FILE
	printf("\nLoading configuration...	");
	if (read_configuration(inputfpath, outputpath, &number_of_threads, &final_time, &start_time_save, &N_bodys, body_int, GM, &dv_step, &n, &first_particle_number, &particle_mass, &particle_density, &save_as_binary) != 0)
	{
		printf("\n\nerror:	could not read configuration.\n");
		//SLEEP(4000);
		return 1;
	}
	printf("...done.");

	//LOAD PARTICLES
	printf("\nLoading particles...		");
	FILE *particles_start_file;
	SpiceDouble **particles_start;
	fopen_s(&particles_start_file, inputfpath, "r");
	if (particles_start_file == NULL)
	{
		printf("\n\nerror:	could not load particles.\n");
		//SLEEP(4000);
		return 1;
	}
	j = 0;
	while ((c = fgetc(particles_start_file)) != EOF) //#tmp# requires newline before eof
	{
		if (c == '%')
		{
			commentLine = true;
		}
		else if (c == '\n')
		{
			if (commentLine == false)
			{
				particles_count++;
			}
			else
			{
				commentLine = false;
				nCommentLines++;
			}
		}
	}
	particles_start = malloc((particles_count + 1) * sizeof(SpiceDouble *));
	for (j = 0; j < particles_count; j++)
	{
		particles_start[j] = malloc(8 * sizeof(SpiceDouble));
	}
	fclose(particles_start_file);
	fopen_s(&particles_start_file, inputfpath, "r");
	j = -nCommentLines;
	while (fgets(temp, sizeof(temp), particles_start_file) != NULL)
	{
		if (j >= 0)
		{
			char* cval = strtok_r(temp, "\t", &next_token);
			for (g = 0; g < 6; g++)
			{
				sscanf(cval, "%lf", &particles_start[j][g]);
				cval = strtok_r(NULL, "\t", &next_token);
			}
			sscanf(cval, "%lf", &particles_start[j][6]);
			cval = strtok_r(NULL, "\n", &next_token);
			sscanf(cval, "%lf", &particles_start[j][7]);
		}
		j++;
	}
	int last_particle_number = first_particle_number + particles_count - 1;
	fclose(particles_start_file);
	printf("...done. %d particles loaded.\n", particles_count);

	//Print config
	if (number_of_threads > 1)
		printf("\n number of threads	= %d", number_of_threads);
	printf("\n save as binary   	= %d", save_as_binary);
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
			if (progress == NULL) // If opening file failed, wait 100 ms and try again.
			{
				perror("The following error occurred");
				SLEEP(100);
				if (e == 2) // After 3 failed attempts, abort.
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
				perror("The following error occurred");
				SLEEP(100);
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

#ifdef __WTIMING
	clock_t start = clock(); // Start clock
#endif
	
	//Start (parallel) computing
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
			char particle_path[260] = "";
			sprintf_s(particle_path, 260, "%s_#%d%s", outputpath, p, ".txt");

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
						SLEEP(100);
						if (e == 2)
						{
							err = 1;
							printf("\nerror: could not create initial output file");
							break;
						}
					}
					else
					{
						printpdata(init, nstate);
						fclose(init);
						break;
					}
				}
			}


			//Create File for output
			FILE *statefile;
			fopen_s(&statefile, particle_path, "a");
			if (statefile == NULL)
			{
				printf("\nerror: could not write to output file");
				err = 1;
			}

			//Integrate particle
			if (err == 0)
			{
				err = RungeKutta4(N_bodys, body_int, GM, final_time, start_time_save, dv_step, nstate, statefile, n, particles_count);
			}
			
			fclose(statefile);


			//Write the particle number to the already-done file and update progress.txt
			FILE* done;
			double fraction;
#pragma omp critical(ALREADYDONE)
			{
				particles_done++;
				for (e = 0; e < 3; e++)
				{
					fopen_s(&done, already_done_path, "a+");
					if (done == NULL)
					{
						SLEEP(100);
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
						SLEEP(100);
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
				printf("\n particle #%d	was not successfully completed", p);
			}
			else
			{
				printf("\n particle #%d	done on thread %d", p, th_id);
			}
		}

#ifdef __WTIMING
		//Print elapsed time
#pragma omp barrier
#pragma omp master
		{
			remove(already_done_path);
			//Print time
			clock_t end = clock();
			double elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;
			printf("\n\n Elapsed time: %1.3f s", elapsed_time);
		}
#endif
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
	//SLEEP(3000);
}



//Functions

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
				sscanf(temp, "%d", &particle_ID);
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
	char temp1[260], temp2[260], temp3[260], *next_token = NULL;
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
						strcpy(temp3, 260, temp2);
					}
					strcpy(temp2, 260, temp1);
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
					char* cval = strtok_r(temp3, "\t", &next_token);
					sscanf(cval, "%lf", &(nstate)[0]);
					cval = strtok_r(NULL, "\t", &next_token);
					sscanf(cval, "%lf", &(nstate)[1]);
					cval = strtok_r(NULL, "\t", &next_token);
					sscanf(cval, "%lf", &(nstate)[2]);
					cval = strtok_r(NULL, "\t", &next_token);
					sscanf(cval, "%lf", &(nstate)[3]);
					cval = strtok_r(NULL, "\t", &next_token);
					sscanf(cval, "%lf", &(nstate)[4]);
					cval = strtok_r(NULL, "\t", &next_token);
					sscanf(cval, "%lf", &(nstate)[5]);
					cval = strtok_r(NULL, "\n", &next_token);
					sscanf(cval, "%lf", &(nstate)[6]);
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



/* define the struct type for the config file readout */
typedef struct
{
	// Simulation
	const char* finaltime;
	const char* starttimes;
	int nbodys;
	const char* bodysid;
	const char* dvstep;
	int mult;
	int nthreads;
	int saveas_binary;
	// Particles
	const char* inputfn;
	const char* outputfn;
	const char* pmass;
	int pdensity;
	int fpnum;
} configuration;

static int handler(void* user, const char* section, const char* name, const char* value)
{
	configuration* pconfig = (configuration*)user;

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	if (MATCH("simulation", "FINAL_TIME")) {
		pconfig->finaltime = strdup(value);
	}
	else if (MATCH("simulation", "START_TIME_SAVE")) {
		pconfig->starttimes = strdup(value);
	}
	else if (MATCH("simulation", "N_BODYS")) {
		pconfig->nbodys = atoi(value);
	}
	else if (MATCH("simulation", "BODYS_ID")) {
		pconfig->bodysid = strdup(value);
	}
	else if (MATCH("simulation", "DV_STEP")) {
		pconfig->dvstep = strdup(value);
	}
	else if (MATCH("simulation", "SAVE_NTH_MULTIPLIER")) {
		pconfig->mult = atoi(value);
	}
	else if (MATCH("simulation", "NUMBER_OF_THREADS")) {
		pconfig->nthreads = atoi(value);
	}
	else if (MATCH("simulation", "SAVE_AS_BINARY")) {
		pconfig->saveas_binary = atoi(value);
	}
	else if (MATCH("particles", "PARTICLE_INPUT_FILE_NAME")) {
		pconfig->inputfn = strdup(value);
	}
	else if (MATCH("particles", "PARTICLE_OUTPUT_FILE_NAME")) {
		pconfig->outputfn = strdup(value);
	}
	else if (MATCH("particles", "PARTICLE_MASS")) {
		pconfig->pmass = strdup(value);
	}
	else if (MATCH("particles", "PARTICLE_DENSITY")) {
		pconfig->pdensity = atoi(value);
	}
	else if (MATCH("particles", "FIRST_PARTICLE_NUMBER")) {
		pconfig->fpnum = atoi(value);
	}
	else {
		return 0;  /* unknown section/name, error */
	}
	return 1;
}

int read_configuration(char *inputfpath, char *outputpath, int *number_of_threads, SpiceDouble *final_time, SpiceDouble *start_time_save, int *N_bodys, int *body_int, SpiceDouble *GM, SpiceDouble *dv_step, int *n, int *first_particle_number, SpiceDouble *particle_mass, SpiceDouble *particle_density, int *save_as_binary)
{
	char temp[260], *token, *next_token = NULL, inputpath[260] = ("INPUT" OS_SEP), configpath[260] = "";
    SpiceInt dim, j;
	configuration config;

	sprintf_s(configpath, 260, "%s%s", inputpath, "configuration.ini");
	
	// Set default values
	config.finaltime = "";
	config.starttimes = "1 JAN 1000";
	config.nbodys = 0;
	config.bodysid = "";
	config.dvstep = "10e-5";
	config.mult = 20;
	config.nthreads = 1;
	config.saveas_binary = 0;
	config.inputfn = "";
	config.outputfn = "default";
	config.pmass = 0;
	config.pdensity = 1000;
	config.fpnum = 1;

	// Parse configuration file
	if (ini_parse(configpath, handler, &config) < 0) {
		printf("Can't load 'configuration.ini'\n");
		return 2;
	}

	//Set number of threads
	*number_of_threads = config.nthreads;

	*save_as_binary = config.saveas_binary;

	//Set final date of the simulation
	if (strcmp(config.finaltime, "") == 0)
	{
		printf("\n\nerror:	FINAL_TIME not set");
		SLEEP(1000);
		return 1;
	}
	str2et_c(config.finaltime, final_time);

	//Set start date for saving
	str2et_c(config.starttimes, start_time_save);
	
	//Set bodys
	if (config.nbodys == 0)
	{
		printf("\n\nerror:	N_BODYS not set");
		SLEEP(1000);
		return 1;
	}
	*N_bodys = config.nbodys;
	
	strcpy(temp, sizeof(temp), config.bodysid);
	token = strtok_r(temp, " ", &next_token);
	for (j = 0; j < *N_bodys; j++)
	{
		if (token == NULL)
		{
			printf("\n\nerror:	BODYS_ID not set or not enough arguments");
			SLEEP(1000);
			return 1;
		}
		sscanf(token, "%d", &body_int[j]);
		token = strtok_r(NULL, " ", &next_token);
		bodvcd_c(body_int[j], "GM", *N_bodys, &dim, &GM[j]); // Get standard gravitational parameter of each body (GM)
	}

	//Set step size control
	sscanf(config.dvstep, "%lf", dv_step);
	
	//Set which nth state is saved to disc
	if (config.mult == 0)
	{
		//Save every 10nth state. This produces high density of states in the output file and is intended to be used when testing the integrator.
		*n = 10;
	}
	else
	{
		*n = (int)(config.mult * pow(*dv_step, -1) + 0.5);
	}

	//Set which particle to start and end with (particle number, from 1 to the number of particles in the input file)
	*first_particle_number = config.fpnum;
	
	//Set name of the input/output file
	strcpy(temp, sizeof(temp), config.inputfn);
	if (strcmp(temp, "") == 0)
	{
		printf("\n\nerror:	PARTICLE_INPUT_FILE_NAME not set");
		SLEEP(1000);
		return 1;
	}
	sprintf_s(inputfpath, 260, "%s%s%s", inputpath, temp, ".txt");
	if (strcmp(config.outputfn, "default"))
	{
		strcpy(temp, sizeof(temp), config.outputfn);
	}
	mkdir(outputpath, 0777);
	char outputfile[260] = "";
	sprintf_s(outputfile, 260, "%s%s", outputpath, temp);
	strcpy(outputpath, 260, outputfile);
	//Set mass of particles
	//strcpy(temp, sizeof(temp), config.pmass);
	sscanf(config.pmass, "%lf", particle_mass);
	if (*particle_mass > 0)
	{
		//Set density of particles
		*particle_density = (SpiceDouble)config.pdensity;

		//Manipulate sun mass to simulate solar pressure
		SpiceDouble Qpr, PI, particle_radius, beta;
		Qpr = 1.0;
		PI = 3.1416;
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
