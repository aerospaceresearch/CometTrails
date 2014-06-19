#ifdef _WIN32
	#include <stdio.h>
#else
	#define _GNU_SOURCE // for fcloseall() on linux
	#include <stdio.h>
#endif
#include <stdbool.h>
#include <string.h>
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

// Check for leaked memory at the end of the run (MSVC only): Add preprocessor definition "__CHKMEMLEAK"
#ifdef __CHKMEMLEAK
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>
#else // Always include
	#include <stdlib.h>
#endif // __CHKMEMLEAK



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
	#define fcloseall _fcloseall
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
#include <PI_types.h> // configuration_values, configuration_readout
#include <IntegEnv.h> // calc_accel(), printpdata(), calc_pInfo(), return_SSB()
#include <RungeKutta4.h>
#include <RungeKutta76.h>

bool particle_already_processed(int p, char already_done_path[]);
bool particle_incomplete(char outputpath[], SpiceDouble *nstate);
int read_configuration(configuration_values *config_data);
int convert_results_into_binary(configuration_values config_data, int particles_count, double *multiplication_factor, char already_done_path[]);
void printinfo();

//Main Program
int main(int argc, char *argv[])
{
	// Version and build info output
	printinfo();

	//Create some variables
	int j, k, e, p, g, c, error_code = 0, particles_count = 0, particles_done = 0, nCommentLines = 0;
	char temp[260], *next_token = NULL, already_done_path[260] = "INPUT" OS_SEP "processed_particles.txt";
	bool commentLine = false;

	// Check for command-line argument -t: save runtime.txt file, for example for MATLAB processing
	FILE* runtimefn;
	bool printtimefile = 0;
	if (argc >= 2)
	{
		for (j = 0; j < argc; j++)
		{
			if (j)
			{
				printf("\nInput argument set: ");
				for (k = 0; k < (int)strlen(argv[j]); k++)
				{
					printf("%c", argv[j][k]);
				}
			}
			if (strcmp(argv[j], "-t") == 0)
			{
				printf("\n Will save runtime.txt.");
				fopen_s(&runtimefn, "runtime.txt", "w");
				printtimefile = 1;
				break;
			}
		}
	}

	// Initialize config_data
	configuration_values config_data = 
	{
		.algorithm = 0,
		.outputpath = "OUTPUT" OS_SEP,
		.number_of_threads = 0,
		.final_time = 0.,
		.start_time_save = 0.,
		.dv_step = 0.,
		.e_target = 0.,
		.first_particle_number = 0,
		.particle_mass = 0.,
		.particle_density = 0.,
		.particle_radius = 0.,
		.save_as_binary = 0
	};
	// Allocate memory below body char pointers
	for (j = 0; j < 10; j++)
	{
		config_data.body_char[j] = malloc(3 * sizeof(ConstSpiceChar));
	}

	//Load Spice kernels
	printf("\nLoading kernels...		");
	furnsh_c("kernels_generic.txt");
	printf("...done.");

	// Read configuration file
	printf("\nLoading configuration...	");
	if (read_configuration(&config_data) != 0)
	{
		printf("\n\nerror:	could not read configuration.\n");
		return 1;
	}
	printf("...done.");

	printf("\nProcessing configuration...	");
	if (calc_pInfo(&config_data) != 0)
	{
		printf("\n\nerror:	could not process configuration values.\n");
		return 1;
	}
	printf("...done.");

	//LOAD PARTICLES
	printf("\nLoading particles...		");
	FILE *particles_start_file;
	SpiceDouble **particles_start;
	fopen_s(&particles_start_file, config_data.inputfpath, "r");
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
	fopen_s(&particles_start_file, config_data.inputfpath, "r");
	j = -nCommentLines;
	while (fgets(temp, sizeof(temp), particles_start_file) != NULL)
	{
		if (j >= 0) // Warning: sscanf functions will cause a crash if parameters are missing in the input file
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
	int last_particle_number = config_data.first_particle_number + particles_count - 1;
	fclose(particles_start_file);
	printf("...done. %d particles loaded.\n", particles_count);

	//Print config
	if (config_data.algorithm == 1)
		printf("\n algorithm		= RK4");
	else if (config_data.algorithm == 2)
	{
		printf("\n algorithm		= RK76");
		printf("\n interpolation order	= %d", config_data.interp_order);
	}
	else
		printf("\n algorithm unknown.");
	if (config_data.number_of_threads > 1)
		printf("\n number of threads	= %d", config_data.number_of_threads);
	printf("\n final_time		= %.16le", config_data.final_time);
	if (config_data.start_time_save > (double)-3.155e+10)
		printf("\n start_time_save	= %.6le", config_data.start_time_save);
	if (config_data.save_as_binary){
		printf("\n saving output as	  binary (.ctwu)");
	}
	else {
		printf("\n saving output as	  text (.txt)");
	}
	printf("\n e_save_slope		= %d", config_data.e_save_slope);
	printf("\n e_save_max  		= %d", config_data.e_save_max);
	if (config_data.endontime){
		printf("\n end on time		  yes");
	}
	else {
		printf("\n end on time		  no");
	}
	printf("\n bodys_ID		=");
	for (j = 0; j < config_data.N_bodys; j++)
		printf(" %d", config_data.body_int[j]);
	if (config_data.ssb_centered == 1)
		printf("\n ssb-centered		= %d", config_data.ssb_centered);
	if (config_data.algorithm == 1)
		printf("\n dv_step		= %le", config_data.dv_step);
	else if (config_data.algorithm == 2)
		printf("\n e_target		= %le", config_data.e_target);
	if (config_data.particle_mass > 0.)
	{
		printf("\n particle_mass   	= %le", config_data.particle_mass);
		printf("\n particle_density	= %le", config_data.particle_density);
		printf("\n particle_radius	= %.12le", config_data.particle_radius);
		printf("\n beta           	= %.12le", config_data.beta);
	}
	for (j = 0; j < config_data.N_bodys; j++)
	{
		if (config_data.body_int[j] == 10)
			printf("\n beta-corr. solar GM	= %.12le", config_data.GM[j] * (1. - config_data.beta));
	}
	if (config_data.first_particle_number != 1)
		printf("\n first_particle_number	= %d", config_data.first_particle_number);
	printf("\n save_nth		= %d", config_data.n);

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
		printf("\n\n First run. Creating processed_particles.txt");
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
#pragma omp parallel private(j, e) num_threads(config_data.number_of_threads)
	{
		int th_id = omp_get_thread_num();

		//Allocate nstate
		SpiceDouble nstate[7];

		//Loop over particles
#pragma omp for
		for (p = config_data.first_particle_number; p <= last_particle_number; p++)
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
			sprintf_s(particle_path, 260, "%s_#%d%s", config_data.outputpath, p, ".txt");

			//Set particle start_time
			SpiceDouble start_time = particles_start[p - config_data.first_particle_number][7];

			//Check if particle has been processed but is incomplete
			if (particle_incomplete(particle_path, nstate))
			{
				printf("\n particle #%d	will be continued", p);
			}
			else
			{
				//Set initial nstate
				for (j = 0; j < 6; j++)
					nstate[j] = particles_start[p - config_data.first_particle_number][j];
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
							printf("\n\nerror: could not create initial output file");
							break;
						}
					}
					else
					{
						if (printpdata(init, nstate))
						{
							err = 1;
							printf("\n\nerror: initial output is NAN.");
							break;
						}
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
				printf("\n\nerror: could not write to output file");
				err = 1;
			}

			//Integrate particle
			if (err == 0)
			{
				switch (config_data.algorithm)
				{
				case 1:
					err = RungeKutta4(&config_data, nstate, statefile);
					break;
				case 2:
					err = RungeKutta76(&config_data, nstate, statefile);
					break;
				default:
					err = 1;
					printf("\n\nerror: unknown integration algorithm: %d", config_data.algorithm);
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
							printf("\n\nerror: could not write to processed_particles.txt");
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
							printf("\n\nwarning: could not write to progress.txt	(non-relevant)");
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

			if (printtimefile == 1)
			{
				if (runtimefn == NULL)
				{
					printf("\n\nwarning: could not write to runtime.txt	(non-relevant)");
				}
				else
				{
					fprintf(runtimefn, "%.8le", elapsed_time);
				}
			}
		}
#endif // __WTIMING
	}

	//Convert .txt output into binary
	if (config_data.save_as_binary == 1)
	{
		double *multiplication_factor;
		multiplication_factor = malloc(particles_count * sizeof(double));
		if (multiplication_factor == NULL)
		{
			printf("\n\nerror: could not allocate multiplication_factor array (OOM)");
			return 2;
		}
		for (j = 0; j < particles_count; j++)
		{
			multiplication_factor[j] = particles_start[j][6];
		}
		if (convert_results_into_binary(config_data, particles_count, multiplication_factor, already_done_path) != 0)
		{
			printf("\n\nerror: could not convert to binary");
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

#ifdef __CHKMEMLEAK
	_CrtDumpMemoryLeaks();
#endif
}



//Functions

bool particle_already_processed(int p, char already_done_path[])
{
	FILE* check;
	char temp[16] = "";
	int particle_ID;
	bool answer = false;
#pragma omp critical(ALREADYDONE)
	{
		fopen_s(&check, already_done_path, "r");
		if (check != NULL)
		{
			while (fgets(temp, 16, check) != NULL)
			{
				sscanf(temp, "%d", &particle_ID);
				if (particle_ID == p)
				{
					answer = true;
				}
			}
			fclose(check);
		}
		else
		{
			printf("\n\nwarning: could not access processed_particles.txt");
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
				printf("\n\nerror: could not create tempfile; particle restarted");
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
						printf("\n\nerror: could not create tempfile; particle restarted");
						err = 1;
					}
					else
					{
						fopen_s(&check, particle_path, "w");
						if (check == NULL)
						{
							printf("\n\nerror: could not read incomplete outputfile; particle restarted");
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
		free(pconfig->algo);
		pconfig->algo = strdup(value);
	}
	else if (MATCH("simulation", "SSB_CENTERED")) {
		pconfig->ssbc = atoi(value);
	}
	else if (MATCH("simulation", "FINAL_TIME")) {
		free(pconfig->finaltime);
		pconfig->finaltime = strdup(value);
	}
	else if (MATCH("simulation", "START_TIME_SAVE")) {
		free(pconfig->starttimes);
		pconfig->starttimes = strdup(value);
	}
	else if (MATCH("simulation", "N_BODYS")) {
		pconfig->nbodys = atoi(value);
	}
	else if (MATCH("simulation", "BODYS_ID")) {
		free(pconfig->bodysid);
		pconfig->bodysid = strdup(value);
	}
	else if (MATCH("simulation", "ENDONTIME")) {
		pconfig->endontime = atoi(value);
	}
	else if (MATCH("rk4", "DV_STEP")) {
		free(pconfig->dvstep);
		pconfig->dvstep = strdup(value);
	}
	else if (MATCH("rk76", "E_TARGET")) {
		free(pconfig->etarget);
		pconfig->etarget = strdup(value);
	}
	else if (MATCH("rk76", "IORDER")) {
		pconfig->iorder = atoi(value);
	}
	else if (MATCH("simulation", "SAVE_NTH_MULTIPLIER")) { // backward compatibility, now in [saving]
		free(pconfig->mult);
		pconfig->mult = strdup(value);
	}
	else if (MATCH("simulation", "NUMBER_OF_THREADS")) {
		pconfig->nthreads = atoi(value);
	}
	else if (MATCH("simulation", "SAVE_AS_BINARY")) { // backward compatibility, now in [saving]
		pconfig->savebin = atoi(value);
	}
	else if (MATCH("particles", "PARTICLE_INPUT_FILE_NAME")) {
		free(pconfig->inputfn);
		pconfig->inputfn = strdup(value);
	}
	else if (MATCH("particles", "PARTICLE_OUTPUT_FILE_NAME")) {
		free(pconfig->outputfn);
		pconfig->outputfn = strdup(value);
	}
	else if (MATCH("particles", "PARTICLE_MASS")) {
		free(pconfig->pmass);
		pconfig->pmass = strdup(value);
	}
	else if (MATCH("particles", "PARTICLE_DENSITY")) {
		free(pconfig->pdensity);
		pconfig->pdensity = strdup(value);
	}
	else if (MATCH("particles", "Q_PR")) {
		free(pconfig->q_pr);
		pconfig->q_pr = strdup(value);
	}
	else if (MATCH("particles", "FIRST_PARTICLE_NUMBER")) {
		pconfig->fpnum = atoi(value);
	}
	else if (MATCH("saving", "SAVE_NTH_MULTIPLIER")) {
		free(pconfig->mult);
		pconfig->mult = strdup(value);
	}
	else if (MATCH("saving", "SAVE_AS_BINARY")) {
		pconfig->savebin = atoi(value);
	}
	else if (MATCH("saving", "ENCOUNTER_SLOPE")) {
		pconfig->e_slope = atoi(value);
	}
	else if (MATCH("saving", "ENCOUNTER_MAX")) {
		pconfig->e_max = atoi(value);
	}
	else {
		return 0;  /* unknown section/name, error */
	}
	return 1;
}

int read_configuration(configuration_values *config_data)
{
	char temp[260] = "", *token, *next_token = NULL, inputpath[260] = ("INPUT" OS_SEP), configpath[260] = "";
	SpiceInt dim, j;
	SpiceDouble mult = 0.0;
	
	sprintf_s(configpath, 260, "%s%s", inputpath, "configuration.ini");

	// Set default values: initialize config struct
	configuration_readout config =
	{
		/* [particles] */
		.inputfn = (char *)malloc(261), // not C90 compatible
		.outputfn = (char *)malloc(261),
		.pmass = (char *)malloc(31),
		.q_pr = (char *)malloc(31),
		.pdensity = (char *)malloc(31),
		.fpnum = 1,

		/* [simulation] */
		.algo = (char *)malloc(11),
		.ssbc = 0,
		.finaltime = (char *)malloc(101),
		.starttimes = (char *)malloc(101),
		.nbodys = 0,
		.bodysid = (char *)malloc(41),
		.nthreads = 1,
		.endontime = 0,

		/* [saving] */
		.mult = (char *)malloc(31),
		.savebin = 1,
		.e_slope = 4,
		.e_max = 40,

		/* Algorithm-specific */
		/* [rk4] */
		.dvstep = (char *)malloc(31),
		/* [rk76] */
		.etarget = (char *)malloc(31),
		.iorder = 5
	};

	if (config.algo == NULL 
		|| config.finaltime == NULL 
		|| config.starttimes == NULL 
		|| config.bodysid == NULL
		|| config.mult == NULL
		|| config.inputfn == NULL
		|| config.outputfn == NULL
		|| config.pmass == NULL
		|| config.q_pr == NULL
		|| config.pdensity == NULL
		|| config.dvstep == NULL
		|| config.etarget == NULL) // At least one alloc ran OOM
	{
		printf("\n\nerror: could not allocate memory for config char*s (OOM)");
		return 1;
	}
	
	// set default values for char*s
	strcpy(config.algo, 10, "RK76");
	strcpy(config.finaltime, 100, "");
	strcpy(config.starttimes, 100, "1 JAN 1000");
	strcpy(config.bodysid, 40, "10");
	strcpy(config.mult, 30, "20.");
	strcpy(config.inputfn, 260, "");
	strcpy(config.outputfn, 260, "default");
	strcpy(config.pmass, 30, "0.");
	strcpy(config.q_pr, 30, "1.");
	strcpy(config.pdensity, 30, "1000.");
	strcpy(config.dvstep, 30, "10e-3");
	strcpy(config.etarget, 30, "10e-18");

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
		config_data->algorithm = 1;
	}
	else if (strcmp(config.algo, "RK76") == 0)
	{
		config_data->algorithm = 2;
	}
	else
	{
		config_data->algorithm = 0; // invalid input
	}

	//Center bodies at SSB?
	config_data->ssb_centered = (bool)config.ssbc;

	//Set number of threads
	config_data->number_of_threads = config.nthreads;

	//Save output as binary?
	config_data->save_as_binary = (bool)config.savebin;

	//Saving increase slope
	config_data->e_save_slope = config.e_slope;

	//Maximum increase in save rate
	config_data->e_save_max = config.e_max;

	//End on time?
	config_data->endontime = (bool)config.endontime;

	//Set final date of the simulation
	if (strcmp(config.finaltime, "") == 0)
	{
		printf("\n\nerror:	FINAL_TIME not set");
		SLEEP(1000);
		return 1;
	}
	str2et_c(config.finaltime, &config_data->final_time);

	//Set start date for saving
	if ((strcmp(config.starttimes, "0") == 0) || (strcmp(config.starttimes, "") == 0))
	{
		config_data->start_time_save = -1e+11; // set start time to a long time ago.
	}
	else 
	{
		str2et_c(config.starttimes, &config_data->start_time_save);
	}


	//Set bodies
	if (config.nbodys == 0)
	{
		printf("\n\nerror:	N_BODYS not set");
		SLEEP(1000);
		return 1;
	}
	config_data->N_bodys = config.nbodys;

	strcpy(temp, sizeof(temp), config.bodysid);
	token = strtok_r(temp, " ", &next_token);
	for (j = 0; j < config_data->N_bodys; j++)
	{
		if (token == NULL)
		{
			printf("\n\nerror:	BODYS_ID not set or not enough arguments");
			SLEEP(1000);
			return 1;
		}
		sscanf(token, "%d", &config_data->body_int[j]);
		token = strtok_r(NULL, " ", &next_token);
		sprintf_s((char *)config_data->body_char[j], 3, "%d", config_data->body_int[j]);
		bodvcd_c(config_data->body_int[j], "GM", config_data->N_bodys, &dim, &config_data->GM[j]); // Get standard gravitational parameter of each body (GM)
	}

	//Set step size control (rk4)
	sscanf(config.dvstep, "%lf", &config_data->dv_step);

	//Set target error per step
	sscanf(config.etarget, "%lf", &config_data->e_target);

	//Set which nth state is saved to disc
	sscanf(config.mult, "%lf", &mult);
	if (mult < 0.00000000001)
	{
		//Save every 10nth state. This produces high density of states in the output file and is intended to be used when testing the integrator.
		config_data->n = 10;
	}
	else
	{
		if (config_data->algorithm == 1) // RK4
		{
			config_data->n = (int)(mult / config_data->dv_step + 0.5);
		}
		else if (config_data->algorithm == 2) // RK76
		{
			// Close to constant number of total steps saved across e_target values. 1.4e3 is a factor imitating the number of steps that would be saved with RK4.
			config_data->n = (int)(mult / 1.4e3 * pow(10,(3.6072 - 0.0746 * log10( config_data->e_target ))) + 0.5);

			// Interpolation order?
			config_data->interp_order = config.iorder;
		}
		else // Unknown algorithm
		{
			config_data->n = 10;
		}
	}
#ifdef __SaveRateOpt
	// Initialize step multiplier value with 1
	config_data->step_multiplier = 1.;
#endif

	//Set which particle to start and end with (particle number, from 1 to the number of particles in the input file)
	config_data->first_particle_number = config.fpnum;

	//Set name of the input/output file
	strcpy(temp, sizeof(temp), config.inputfn);
	if (strcmp(temp, "") == 0)
	{
		printf("\n\nerror:	PARTICLE_INPUT_FILE_NAME not set");
		SLEEP(1000);
		return 1;
	}
	sprintf_s(config_data->inputfpath, 260, "%s%s%s", inputpath, temp, ".txt");
	if (strcmp(config.outputfn, "default"))
	{
		strcpy(temp, sizeof(temp), config.outputfn);
	}

	if (mkdir(config_data->outputpath, 0777))
	{
		printf("...skip mkdir... ");
	}

	char outputfile[260] = "";
	sprintf_s(outputfile, 260, "%s%s", config_data->outputpath, temp);
	strcpy(config_data->outputpath, 260, outputfile);

	// Set Mie scattering coefficient
	sscanf(config.q_pr, "%lf", &config_data->q_pr);

	//Set mass of particles
	sscanf(config.pmass, "%lf", &config_data->particle_mass);

	if (config_data->particle_mass > 0.)
	{
		//Set density of particles
		sscanf(config.pdensity, "%lf", &config_data->particle_density);
	}

	// Free memory allocated for config char*s
	free(config.algo);
	free(config.finaltime);
	free(config.starttimes);
	free(config.bodysid);
	free(config.mult);
	free(config.inputfn);
	free(config.outputfn);
	free(config.pmass);
	free(config.q_pr);
	free(config.pdensity);
	free(config.dvstep);
	free(config.etarget);

	return 0;
}



int convert_results_into_binary(configuration_values config_data, int particles_count, double *multiplication_factor, char already_done_path[])
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
		printf("\n\nerror: could not allocate result_array (OOM)");
		return 2;
	}
	result_array[0] = malloc(7 * sizeof(float));
	if (result_array[0] == NULL)
	{
		printf("\n\nerror: could not allocate result_array (OOM)");
		return 2;
	}
	//Set file header
	result_array[0][0] = (float)config_data.first_particle_number;
	result_array[0][1] = (float)(config_data.first_particle_number + particles_count - 1);
	result_array[0][2] = (float)config_data.particle_mass;
	result_array[0][3] = (float)config_data.particle_density;
	result_array[0][4] = (float)config_data.beta;
	result_array[0][5] = 0;
	result_array[0][6] = 0;

	//Read in all the particles and save them in result_array
	for (j = 0; j < particles_count; j++)
	{
		state_count = 0;
		char particle_path[260] = "";
		sprintf_s(particle_path, 260, "%s_#%d%s", config_data.outputpath, (j + config_data.first_particle_number), ".txt");
		for (e = 0; e < 3; e++)
		{
			fopen_s(&output_file, particle_path, "r");
			if (output_file == NULL)
			{
				SLEEP(100);
				if (e == 2)
				{
					printf("\n\nerror: could not open .txt file particle #%d for conversion", j + config_data.first_particle_number);
					return 2;
				}
			}
			else
			{
				break;
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
			printf("\n\nerror: could not allocate result_array (OOM)");
			return 2;
		}
		for (i = particle_header_row; i < result_array_length; i++)
		{
			result_array[i] = malloc(7 * sizeof(float));
			if (result_array[i] == NULL)
			{
				printf("\n\nerror: could not allocate result_array (OOM)");
				return 2;
			}
		}
		for (i = 0; i < 4; i++)
		{
			result_array[particle_header_row][i] = 0;
		}
		result_array[particle_header_row][4] = (float)(j + config_data.first_particle_number);
		result_array[particle_header_row][5] = (float)(multiplication_factor[j]);
		result_array[particle_header_row][6] = 0;
		rewind(output_file);
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
		if (j == 0)
		{
			result_array[0][5] = result_array[particle_header_row + 1][6];
		}
		if (j == particles_count - 1)
		{
			result_array[0][6] = result_array[particle_header_row + 1][6];
		}
	}
	//Save result_array as binary file and delete text files
	FILE *binout;
	char binary_path[260] = "OUTPUT" OS_SEP "binary_output.ctwu";
	fopen_s(&binout, binary_path, "wb");
	if (binout == NULL)
	{
		printf("\n\nerror:	could not create binary output file.\n");
		return 1;
	}
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
		remove(already_done_path);
		char particle_path[260] = "";
		sprintf_s(particle_path, 260, "%s_#%d%s", config_data.outputpath, (j + config_data.first_particle_number), ".txt");
		if (remove(particle_path) != 0)
		{
			printf("\n\nerror: could not delete .txt file after conversion");
			return 2;
		}
	}
	return 0;
}


void printinfo()
{
	// Print version
	printf("ParticleIntegrator version " PI_VERSION_MAJOR "." PI_VERSION_MINOR "\n");
	// Print build type
#ifdef RELTYPEDEB
	printf("\n Debug build");
#endif
#ifdef RELTYPERWDI
	printf("\n Release build with debug symbols");
#endif
#ifdef RELTYPEREL
	printf("\n Release build");
#endif

	// Print active options in debug builds
#if defined(RELTYPERWDI) || defined(RELTYPEDEB)
	printf(", "__DATE__ " " __TIME__ "\n Options: ");

	#ifdef __WTIMING
		printf("TIMING ");
	#endif // __WTIMING

	#ifdef __WTIMESTEP
		printf("WTIMESTEP ");
	#endif // __WSTEPINFO

	#ifdef __PRD
		printf("PRD ");
	#endif // __PRD

	#ifdef __SWD
		printf("SWD ");
	#endif // __SWD

	#ifdef __SaveRateOpt
		printf("SaveRateOpt ");
	#endif // __SaveRateOpt

	#ifdef __CHKMEMLEAK
		printf("CHKMEMLEAK ");
	#endif // __CHKMEMLEAK

	printf("\n");
#endif // RELTYPERWDI || RELTYPEDEB
}