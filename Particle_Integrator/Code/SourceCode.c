#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
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

// Activate timing: Add preprocessor definition "__WTIMING"
#ifdef __WTIMING
	#include <time.h>
#endif



//BEGIN	Function cross-platform compatibility
#ifdef _WIN32
	#define SLEEP( a1 ) Sleep( a1 )
	#define mkdir( a1, a2 ) _mkdir( a1 )
#else
	#define SLEEP( a1 ) usleep( a1 * 1000 )
#endif

// Avoid MSVC level 3 warning C4996
#ifdef _WIN32
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
//END  	Function cross-platform compatibility


// Custom header files
#include <PI_types.h>
#include <IntegEnv.h>
#include <RungeKutta4.h>
#include <RungeKutta67.h>

bool particle_already_processed(int p, char already_done_path[]);
bool particle_incomplete(char outputpath[], SpiceDouble *nstate);
int read_configuration(configuration_values *config_out);
int convert_results_into_binary(configuration_values config_out, int particles_count, double *multiplication_factor);


//Main Program
int main(void)
{
	//Print version
	printf("ParticleIntegrator version " PI_VERSION_MAJOR "." PI_VERSION_MINOR "\n");

	//Create some variables
	int j, e, p, g, c, error_code = 0, particles_count = 0, particles_done = 0, nCommentLines = 0;
	char temp[260], *next_token = NULL, already_done_path[260] = "INPUT" OS_SEP "processed_particles.txt";
	bool commentLine = false;
	configuration_values config_out;

	// Initialize
	config_out.algorithm = 0;
	sprintf_s(config_out.inputfpath, 260, "");
	sprintf_s(config_out.outputpath, 260, "OUTPUT" OS_SEP);
	config_out.number_of_threads = 0;
	config_out.final_time = 0;
	config_out.start_time_save = 0;
	config_out.dv_step = 0;
	config_out.e_target = 0;
	config_out.first_particle_number = 0;
	config_out.particle_mass = 0;
	config_out.particle_density = 0;
	config_out.particle_radius = 0;
	config_out.save_as_binary = 0;
	
	//Load Spice kernels
	printf("\nLoading kernels...		");
	furnsh_c("kernels_generic.txt");
	printf("...done.");

	//READ CONFIG FILE
	printf("\nLoading configuration...	");
	if (read_configuration(&config_out) != 0)
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
	fopen_s(&particles_start_file, config_out.inputfpath, "r");
	if (particles_start_file == NULL)
	{
		printf("\n\nerror:	could not load particles.\n");
		//SLEEP(4000);
		return 1;
	}
	j = 0;
	while ((c = fgetc(particles_start_file)) != EOF) // requires newline before eof
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
	fopen_s(&particles_start_file, config_out.inputfpath, "r");
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
	int last_particle_number = config_out.first_particle_number + particles_count - 1;
	fclose(particles_start_file);
	printf("...done. %d particles loaded.\n", particles_count);

	//Print config
	if (config_out.algorithm == 1)
		printf("\n algorithm		= RK4");
	else if (config_out.algorithm == 2)
		printf("\n algorithm		= RK67");
	else
		printf("\n algorithm unknown.");
	if (config_out.number_of_threads > 1)
		printf("\n number of threads	= %d", config_out.number_of_threads);
	printf("\n final_time		= %le", config_out.final_time);
	if (config_out.start_time_save > (double)-3.155e+10)
		printf("\n start_time_save	= %le", config_out.start_time_save);
	if (config_out.save_as_binary){
		printf("\n saving output in binary format.");
	}
	else {
		printf("\n saving output in text format.");
	}
	printf("\n bodys_ID		=");
	for (j = 0; j < config_out.N_bodys; j++)
		printf(" %d", config_out.body_int[j]);
	if (config_out.ssb_centered == 1)
		printf("\n ssb-centered		= %d", config_out.ssb_centered);
	if (config_out.algorithm == 1)
		printf("\n dv_step		= %le", config_out.dv_step);
	else if (config_out.algorithm == 2)
		printf("\n e_target		= %le", config_out.e_target);
	if (config_out.particle_mass > 0)
		printf("\n particle_mass		= %le", config_out.particle_mass);
	printf("\n particle_density	= %le", config_out.particle_density);
	printf("\n particle_radius	= %le", config_out.particle_radius);
	if (config_out.first_particle_number != 1)
		printf("\n first_particle_number	= %d", config_out.first_particle_number);
	printf("\n save_nth		= %d", config_out.n);

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
#pragma omp parallel private(j, e) num_threads(config_out.number_of_threads)
	{
		int th_id = omp_get_thread_num();

		//Allocate nstate
		SpiceDouble nstate[7];

		//Loop over particles
#pragma omp for
		for (p = config_out.first_particle_number; p <= last_particle_number; p++)
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
			sprintf_s(particle_path, 260, "%s_#%d%s", config_out.outputpath, p, ".txt");

			//Set particle start_time
			SpiceDouble start_time = particles_start[p - config_out.first_particle_number][7];

			//Check if particle has been processed but is incomplete
			if (particle_incomplete(particle_path, nstate))
			{
				printf("\n particle #%d	will be continued", p);
			}
			else
			{
				//Set initial nstate
				for (j = 0; j < 6; j++)
					nstate[j] = particles_start[p - config_out.first_particle_number][j];
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
				switch (config_out.algorithm)
				{
				case 1:
					err = RungeKutta4(&config_out, nstate, statefile);
					break;
				case 2:
					err = RungeKutta67(&config_out, nstate, statefile);
					break;
				default:
					err = 1;
					printf("\nerror: unknown integration algorithm: %d", config_out.algorithm);
				}
				fclose(statefile);
			}
			
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
						fraction = (float)particles_done / particles_count;
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
#endif // __WTIMING
	}

	//Convert .txt output into binary
	if (config_out.save_as_binary == 1)
	{
		double *multiplication_factor;
		multiplication_factor = malloc(particles_count * sizeof(double));
		if (multiplication_factor == NULL)
		{
			printf("\nerror: could not allocate multiplication_factor array (OOM)");
			return 2;
		}
		for (j = 0; j < particles_count; j++)
		{
			multiplication_factor[j] = particles_start[j][6];
		}
		if (convert_results_into_binary(config_out, particles_count, multiplication_factor) != 0)
		{
			printf("\nerror: could not convert to binary");
			return 2;	
		}
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



static int handler(void* user, const char* section, const char* name, const char* value)
{
	configuration_readout* pconfig = (configuration_readout*)user;

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

	if (MATCH("simulation", "ALGORITHM")) {
		pconfig->algo = strdup(value);
	}
	if (MATCH("simulation", "SSB_CENTERED")) {
		pconfig->ssbc = atoi(value);
	}
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
	else if (MATCH("rk4", "DV_STEP")) {
		pconfig->dvstep = strdup(value);
	}
	else if (MATCH("rk67", "E_TARGET")) {
		pconfig->etarget = strdup(value);
	}
	else if (MATCH("simulation", "SAVE_NTH_MULTIPLIER")) {
		pconfig->mult = strdup(value);
	}
	else if (MATCH("simulation", "NUMBER_OF_THREADS")) {
		pconfig->nthreads = atoi(value);
	}
	else if (MATCH("simulation", "SAVE_AS_BINARY")) {
		pconfig->savebin = atoi(value);
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

int read_configuration(configuration_values *config_out)
{
	char temp[260], *token, *next_token = NULL, inputpath[260] = ("INPUT" OS_SEP), configpath[260] = "";
	SpiceInt dim, j;
	configuration_readout config;
	SpiceDouble mult = 0.0;

	sprintf_s(configpath, 260, "%s%s", inputpath, "configuration.ini");

	// Set default values
	config.algo = "RK4";
	config.ssbc = 0;
	config.finaltime = "";
	config.starttimes = "1 JAN 1000";
	config.nbodys = 0;
	config.bodysid = "";
	config.dvstep = "10e-3";
	config.etarget = "10e-15";
	config.mult = "20.";
	config.nthreads = 1;
	config.savebin = 0;
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

	//Set algorithm
	if (config.algo == 0)
	{
		printf("\n\nerror:	ALGORITHM not set");
		SLEEP(1000);
		return 1;
	}
	if (strcmp(config.algo, "RK4") == 0)
	{
		config_out->algorithm = 1;
	}
	else if (strcmp(config.algo, "RK67") == 0)
	{
		config_out->algorithm = 2;
	}
	else
	{
		config_out->algorithm = 0; // invalid input
	}

	//Center bodies at SSB?
	config_out->ssb_centered = config.ssbc;

	//Set number of threads
	config_out->number_of_threads = config.nthreads;

	//Save output as binary?
	config_out->save_as_binary = config.savebin;

	//Set final date of the simulation
	if (strcmp(config.finaltime, "") == 0)
	{
		printf("\n\nerror:	FINAL_TIME not set");
		SLEEP(1000);
		return 1;
	}
	str2et_c(config.finaltime, &config_out->final_time);

	//Set start date for saving
	if ((strcmp(config.starttimes, "0") == 0) || (strcmp(config.starttimes, "") == 0))
	{
		str2et_c("1 JAN 1000", &config_out->start_time_save);
	}
	else 
	{
		str2et_c(config.starttimes, &config_out->start_time_save);
	}


	//Set bodies
	if (config.nbodys == 0)
	{
		printf("\n\nerror:	N_BODYS not set");
		SLEEP(1000);
		return 1;
	}
	config_out->N_bodys = config.nbodys;

	strcpy(temp, sizeof(temp), config.bodysid);
	token = strtok_r(temp, " ", &next_token);
	for (j = 0; j < config_out->N_bodys; j++)
	{
		if (token == NULL)
		{
			printf("\n\nerror:	BODYS_ID not set or not enough arguments");
			SLEEP(1000);
			return 1;
		}
		sscanf(token, "%d", &config_out->body_int[j]);
		token = strtok_r(NULL, " ", &next_token);
		bodvcd_c(config_out->body_int[j], "GM", config_out->N_bodys, &dim, &config_out->GM[j]); // Get standard gravitational parameter of each body (GM)
	}

	//Set step size control (rk4)
	sscanf(config.dvstep, "%lf", &config_out->dv_step);

	//Set target error per step
	sscanf(config.etarget, "%lf", &config_out->e_target);

	//Set which nth state is saved to disc
	sscanf(config.mult, "%lf", &mult);
	if (mult < 0.0000000000001)
	{
		//Save every 10nth state. This produces high density of states in the output file and is intended to be used when testing the integrator.
		config_out->n = 10;
	}
	else
	{
		config_out->n = (int)(mult / config_out->dv_step + 0.5);
	}

	//Set which particle to start and end with (particle number, from 1 to the number of particles in the input file)
	config_out->first_particle_number = config.fpnum;

	//Set name of the input/output file
	strcpy(temp, sizeof(temp), config.inputfn);
	if (strcmp(temp, "") == 0)
	{
		printf("\n\nerror:	PARTICLE_INPUT_FILE_NAME not set");
		SLEEP(1000);
		return 1;
	}
	sprintf_s(config_out->inputfpath, 260, "%s%s%s", inputpath, temp, ".txt");
	if (strcmp(config.outputfn, "default"))
	{
		strcpy(temp, sizeof(temp), config.outputfn);
	}

	if (mkdir(config_out->outputpath, 0777))
	{
		printf("\n ...skip mkdir... ");
	}

	char outputfile[260] = "";
	sprintf_s(outputfile, 260, "%s%s", config_out->outputpath, temp);
	strcpy(config_out->outputpath, 260, outputfile);
	//Set mass of particles
	//strcpy(temp, sizeof(temp), config.pmass);
	sscanf(config.pmass, "%lf", &config_out->particle_mass);
	if (config_out->particle_mass > 0)
	{
		//Set density of particles
		config_out->particle_density = (SpiceDouble)config.pdensity;

		//Manipulate sun mass to simulate solar pressure
		SpiceDouble Qpr, PI, beta;
		Qpr = 1.0;
		PI = 3.1416;
		config_out->particle_radius = pow((config_out->particle_mass) / ((1.3333) * PI * (config_out->particle_density)), 0.3333); // Unit: [m]!
		beta = 5.7e-4 * (Qpr / ((config_out->particle_density) * config_out->particle_radius));
		for (j = 0; j < config_out->N_bodys; j++)
		{
			if (config_out->body_int[j] == 10)
			{
				config_out->GM[j] = config_out->GM[j] * (1 - beta);
				break;
			}
		}
	}

	return 0;
}

int convert_results_into_binary(configuration_values config_out, int particles_count, double *multiplication_factor)
{
	printf("\n Converting text output into binary...	");
	//Create some variables
	int j, e, c, h, state_count, i, result_array_length=1, particle_header_row, l, g;
	FILE *output_file;
	char *next_token = NULL, temp[260];
	double tempdouble;

	//Allocate beginning of result_array
	float **result_array;
	result_array = malloc(1 * sizeof(float *));
	if (result_array == NULL)
	{
		printf("\nerror: could not allocate result_array (OOM)");
		return 2;
	}
	result_array[0] = malloc(7 * sizeof(float));
	if (result_array[0] == NULL)
	{
		printf("\nerror: could not allocate result_array (OOM)");
		return 2;
	}
	//Set file header
	result_array[0][0] = config_out.first_particle_number;
	result_array[0][1] = (config_out.first_particle_number + particles_count);
	result_array[0][2] = (float)config_out.particle_mass;
	result_array[0][3] = (float)config_out.particle_density;
	result_array[0][4] = 0;
	result_array[0][5] = (float)config_out.start_time_save;
	result_array[0][6] = (float)config_out.final_time;

	//Read in all the particles and save them in result_array
	for (j = 0; j < particles_count; j++)
	{
		state_count = 0;
		char particle_path[260] = "";
		sprintf_s(particle_path, 260, "%s_#%d%s", config_out.outputpath, (j + config_out.first_particle_number), ".txt");
		for (e = 0; e < 3; e++)
		{
			fopen_s(&output_file, particle_path, "r");
			if (output_file == NULL)
			{
				SLEEP(100);
				if (e == 2)
				{
					printf("\nerror: could not open .txt file for conversion");
					return 2;
				}
			}
		}
		while ((c = fgetc(output_file)) != EOF)
		{
			if (c == '\n')
			{
				state_count++;
			}
		}
		particle_header_row = result_array_length;
		result_array_length += state_count + 1;
		result_array = realloc(result_array, (result_array_length)*sizeof(float *));
		if (result_array == NULL)
		{
			printf("\nerror: could not allocate result_array (OOM)");
			return 2;
		}
		for (i = particle_header_row; i < result_array_length; i++)
		{
			result_array[i] = malloc(7 * sizeof(float));
			if (result_array[i] == NULL)
			{
				printf("\nerror: could not allocate result_array (OOM)");
				return 2;
			}
		}
		for (i = 0; i < 5; i++)
		{
			result_array[particle_header_row][i] = 0;
		}
		result_array[particle_header_row][5] = (j + config_out.first_particle_number);
		result_array[particle_header_row][6] = (float)(multiplication_factor[j]);
		fclose(output_file);
		for (e = 0; e < 3; e++)
		{
			fopen_s(&output_file, particle_path, "r");
			if (output_file == NULL)
			{
				SLEEP(100);
				if (e == 2)
				{
					printf("\nerror: could not open .txt file for conversion");
					return 2;
				}
			}
		}
		l = particle_header_row;
		while (fgets(temp, sizeof(temp), output_file) != NULL)
		{
			l++;
			char* cval = strtok_r(temp, "\t", &next_token);
			for (g = 0; g < 5; g++)
			{
				sscanf(cval, "%lf", &tempdouble);
				result_array[l][g] = (float)tempdouble;
				cval = strtok_r(NULL, "\t", &next_token);
			}
			sscanf(cval, "%lf", &tempdouble);
			result_array[l][5] = (float)tempdouble;
			cval = strtok_r(NULL, "\n", &next_token);
			sscanf(cval, "%lf", &tempdouble);
			result_array[l][6] = (float)tempdouble;
		}
		fclose(output_file);
	}
	//Save result_array as binary file and delete text files
	FILE *binout;
	char binary_path[260] = "";
	sprintf_s(binary_path, 260, "%s.ctwu", config_out.outputpath);
	fopen_s(&binout, binary_path, "wb");
	for (h = 0; h < result_array_length; h++)
	{
		fwrite(result_array[h], sizeof(float), 7, binout);
	}
	fclose(binout);
	free(result_array);
	printf("...done.");
	fcloseall();
	for (j = 0; j < particles_count; j++)
	{
		char particle_path[260] = "";
		sprintf_s(particle_path, 260, "%s_#%d%s", config_out.outputpath, (j + config_out.first_particle_number), ".txt");
		if (remove(particle_path) != 0)
		{
			printf("\nerror: could not delete .txt file after conversion");
			return 2;
		}
	}
	return 0;
}