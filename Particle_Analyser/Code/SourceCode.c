#ifdef _WIN32
#include <stdio.h>
#else
#define _GNU_SOURCE // for fcloseall() on linux
#include <stdio.h>
#endif
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <SpiceUsr.h>
#include <PAConfig.h>
#ifdef _WIN32
	#include <windows.h> // needed for search wu directory
	#include <direct.h> // only needed for _mkdir()
#else
	#include <unistd.h> // only needed for usleep()
	#include <dirent.h> // needed for search wu directory
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

#ifdef _WIN32
	// Avoid MSVC level 3 warning C4996
	#define strdup _strdup
	#define strcpy strcpy_s
	#define sscanf sscanf_s
	#define strtok_r strtok_s
	#define strcat strcat_s
	#define fcloseall _fcloseall
#else
	#define strcat( a1, a2, a3 ) strcat(a1, a3)
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
#define PI_half 3.14159265358979323846/2

//Global variables
int velflag = 0;
int suncflag = 0;
int nodeflag = 0;
int tminflag = 0;
int tmaxflag = 0;
int mminflag = 0;
int mmaxflag = 0;
int objectflag = 0;
int intgflag = 0;
int pc2flag = 0;
int WUcheckflag = 0;


int wu_count = 0;
int all_wu_count = 0;
int particles_count = 0;
int all_particles_count = 0;
char *cometwu_path, *output_path;
double timespec[3];
float massspec[2];
char *object_id;
double object_state[6];
double node_orbitplane_normvec[3];
double Sun_GM, Sun_state[6];
int N_bodys = 0;
int body_ID[10];
double body_GM[10];
double dv_step;
double max_distance;


//Functions
int parse_input(int argc, char* argv[]);
char **get_WU_paths(void);
char *search_WUsummary(void);
char **select_WUs(char *WUsummary_name);
char **search_WUs(void);
float *get_particles(char **wu_paths);
int get_nearest_state(int z, float *fwu_array);
int get_state_nearest_to_orbital_plane(int z, float *fwu_array);
float *get_all_target_states(float *nearest_states);
void calculate_target_state(double *nstate, double GM);
void integrate_target_state(double *nstate, double GM);
void calc_accel(int N_bodys, double Sun_GM, double pos[], double **body_state[], double *accel);
float *filter_particles_out_of_range(float *target_states);
int save_target_states(float *target_states);
void WUcheck(char *WUsummary_name);


//Main Program
int main(int argc, char* argv[])
{
#ifdef __WTIMING
	clock_t start = clock(); // Start clock
#endif
	//Print version
	printf("ParticleAnalyser version " PA_VERSION_MAJOR "." PA_VERSION_MINOR "\n");

	//Load Spice kernels
	printf("Loading kernels...		");
	furnsh_c("kernels_generic.txt");
	printf("...done.");
	
	//Parse input arguments
	if (parse_input(argc, argv) != 0){
		return 1;
	}

	//Get paths to all work units that contain relevant particles
	char **wu_paths;
	wu_paths = get_WU_paths();
	if (wu_paths == NULL){
		if (WUcheckflag == 1){
			return 0;
		}
		return 1;
	}	
	//If there are no relevant WUs, but no error has occurred, return success
	if (wu_count == 0){
		free(wu_paths);
		return 0;
	}
	
	//Load all relevant particles and get the states that are closest to target time / orbital plane
	float *nearest_states;
	nearest_states = get_particles(wu_paths);
	if (nearest_states == NULL){
		return 1;
	}
	//If there are no relevant particles, but no error has occurred, return success
	if (particles_count == 0){
		free(nearest_states);
		return 0;
	}
	
	//Calculate target states
	float *target_states;
	target_states = get_all_target_states(nearest_states);
	if (target_states == NULL){
		return 1;
	}

	//Filter particles
	if (objectflag == 1 && nodeflag == 0){
		target_states = filter_particles_out_of_range(target_states);
		if (target_states == NULL){
			return 1;
		}
		if (particles_count == 0){
			free(target_states);
			return 0;
		}
	}

	//Write output file
	if (save_target_states(target_states) != 0){
		return 1;
	}


#ifdef __WTIMING 
	//Print elapsed time
	clock_t end = clock();
	double elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;
	printf("\n\n Elapsed time: %1.3f s\n", elapsed_time);
#endif
	return 0;
}


//Functions
int parse_input(int argc, char *argv[])
{
	int i = 1, x = 3, valid;
	while (i < argc)									/* Scan through args. */
	{
		if (strncmp(argv[i], "--", 2) == 0)				/* Check for option character. */
		{
			if (strspn("input", argv[i]) == 5){			//This is the directory containing the WUs
				cometwu_path = argv[i + 1];
				i = i + 2;
				x--;
			}
			else if (strspn("output", argv[i]) == 6){	//This is the file to write the results to
				output_path = argv[i + 1];
				i = i + 2;
				x--;
			}
			else if (strspn("WUcheck", argv[i]) == 7){
				WUcheckflag = 1;
				i = i + 1;
				x--;
			}
			else if (strspn("time", argv[i]) == 4){		//This is the target time
				timespec[0] = atof(argv[i + 1]);
				i = i + 2;
				x--;
			}
			else if (strspn("tmin", argv[i]) == 4){		//Particles must have formed after this date
				timespec[1] = atof(argv[i + 1]);
				tminflag = 1;
				i = i + 2;
			}
			else if (strspn("tmax", argv[i]) == 4){		//Particles must have formed before this date
				timespec[2] = atof(argv[i + 1]);
				tmaxflag = 1;
				i = i + 2;
			}
			else if (strspn("mmin", argv[i]) == 4){		//Particles must have a mass higher than mmin
				massspec[0] = (float)atof(argv[i + 1]);
				mminflag = 1;
				i = i + 2;
			}
			else if (strspn("mmax", argv[i]) == 4){		//Particles must have a mass lower than mmax
				massspec[1] = (float)atof(argv[i + 1]);
				mmaxflag = 1;
				i = i + 2;
			}
			else if (strspn("object", argv[i]) == 6){	//Specify a SPICE conform body ID and a distance
				object_id = argv[i + 1];				// only particles within that distance will be written to output
				max_distance = atof(argv[i + 2]);
				objectflag = 1;
				i = i + 3;
			}
			else if (strspn("intg", argv[i]) == 4){		//Specifies wether or not to integrate particle target states
				dv_step = atof(argv[i + 1]);					// first parameter is dv_step
				N_bodys = atoi(argv[i + 2]) + 1;				// second parameter is number of planets! Following parameters are NAIF IDs of the planets
				int j;
				SpiceInt dim;
				body_ID[0] = 10;								// Sun is automatically added
				for (j = 1; j < N_bodys; j++){
					body_ID[j] = atoi(argv[i + 2 + j]);
					bodvcd_c(body_ID[j], "GM", N_bodys, &dim, &body_GM[j]); // Get standard gravitational parameter of each body (GM)
				}
				intgflag = 1;
				i = i + 3 + N_bodys - 1;
			}
			else{
				printf("\n\nerror: input argument %d is invalid", i);
				return 1;
			}
		}
		else if (strncmp(argv[i], "-", 1) == 0)
		{
			valid = 0;
			if (strspn("s", argv[i]) != 0){				//Specifies wether or not to convert to Sun-centered coordinates
				suncflag = 1;							// before calculating orbital elements
				valid = 1;
			}
			if (strspn("v", argv[i]) != 0){		//Specifies wether or not velocities will be written to output
				velflag = 1;
				valid = 1;
			}
			if (strspn("n", argv[i]) != 0){		//Specifies wether or not to compute particle nodes
				nodeflag = 1;
				valid = 1;
			}
			if (strspn("p", argv[i]) != 0){		//Specifies wether or not to write pointcache output format
				pc2flag = 1;
				valid = 1;
			}
			if (valid == 1){
				i++;
			}
			else {					//No valid character found
				printf("\n\nerror: input argument %d is invalid", i);
				return 1;
			}
		}
		else{
			printf("\n\nerror: input argument %d is invalid", i);
			return 1;
		}
	}

	//Check overall validity of input
	if (x > 0){
		printf("\n\nerror: not enough input arguments\n arguments --input , --output and --time must always be set\n");
		return 1;
	}

	if (nodeflag == 1 && objectflag == 0){
		printf("\n\nerror: specify an object in order to compute nodes (-n requires --object to be set)\n");
		return 1;
	}

	//Load additional kernels
	if (objectflag == 1 || suncflag == 1){
		furnsh_c("kernels_spk.txt");	//object kernel directory must be specified in "kernels_spk.txt"
	}

	//Set up object location for relative coordinates computation
	double lt;
	if (objectflag == 1){
		if (suncflag == 1){
			spkezr_c(object_id, timespec[0], "ECLIPJ2000", "NONE", "SUN", object_state, &lt);
		}
		else{
			spkezr_c(object_id, timespec[0], "ECLIPJ2000", "NONE", "SSB", object_state, &lt);
		}

		//Set up orbital plane for node computation
		if (nodeflag == 1){
			double length;
			node_orbitplane_normvec[0] = object_state[1] * object_state[5] - object_state[2] * object_state[4];
			node_orbitplane_normvec[1] = object_state[2] * object_state[3] - object_state[0] * object_state[5];
			node_orbitplane_normvec[2] = object_state[0] * object_state[4] - object_state[1] * object_state[3];
			length = sqrt(node_orbitplane_normvec[0] * node_orbitplane_normvec[0] + node_orbitplane_normvec[1] * node_orbitplane_normvec[1] + node_orbitplane_normvec[2] * node_orbitplane_normvec[2]);
			if (length != 0){
				node_orbitplane_normvec[0] = node_orbitplane_normvec[0] / length;
				node_orbitplane_normvec[1] = node_orbitplane_normvec[1] / length;
				node_orbitplane_normvec[2] = node_orbitplane_normvec[2] / length;
			}
			else{
				printf("\n\nwarning: length of orbit-plane normal vector is zero. Body may have coordinates 0/0/0\n");
				printf("		...computing nodes in J2000 eclpiptic\n\n");
				node_orbitplane_normvec[0] = 0;
				node_orbitplane_normvec[1] = 0;
				node_orbitplane_normvec[2] = 1;
			}
		}
	}

	//Get Sun's position relative to Barycentric center (little to no effect)
	if (suncflag == 1){
		spkezr_c("SUN", timespec[0], "ECLIPJ2000", "NONE", "SSB", Sun_state, &lt);
	}
	
	return 0;
}

char **get_WU_paths(void)
{
	char *WUsummary_name, **wu_paths_buffer, **wu_paths;
	WUsummary_name = search_WUsummary();
	if (WUsummary_name == NULL){
		if (WUcheckflag == 1){
			printf("\n Failed to check for missing WUs because no WU summary file was found");
			return NULL;
		}
		wu_paths_buffer = search_WUs();
	}
	else{
		if (WUcheckflag == 1){
			WUcheck(WUsummary_name);
			return NULL;
		}
		wu_paths_buffer = select_WUs(WUsummary_name);
	}
	if (wu_paths_buffer == NULL){
		return NULL;
	}

	int i;
	wu_paths = malloc((wu_count+1) * sizeof(char*));
	if (wu_paths == NULL){
		printf("...failed.");
		printf("\n\nerror: could not allocate wu_paths array (OOM)");
		return NULL;
	}
	for (i = 0; i < wu_count; i++){
		wu_paths[i] = strdup(wu_paths_buffer[i]);
		free(wu_paths_buffer[i]);
	}
	free(wu_paths_buffer);
	return wu_paths;
}

char *search_WUsummary(void){

	//Looking for WU summary file
	printf("\nLooking for WU summary file...	");
	char *WUsummary_name;

#ifdef _WIN32					// Get path WUsummary under windows systems
	HANDLE hFind;
	WIN32_FIND_DATA FindData;
	char findpath[256];
	strcpy(findpath, 256, cometwu_path);
	strcat(findpath, 256, OS_SEP);
	strcat(findpath, 256, "*.txt");
	hFind = FindFirstFile(findpath, &FindData);
	if (hFind == INVALID_HANDLE_VALUE){
		printf("...not found.");
		return NULL;
	}
	WUsummary_name = strdup(FindData.cFileName);
	FindClose(hFind);
	printf("...found.");
	return WUsummary_name;
#else							// Get path WUsummary of WUs under linux
	DIR *dir;
	char *ext;
	struct dirent *ent;
	dir = opendir(cometwu_path);
	if (dir){
		while ((ent = readdir(dir)) != NULL){
			ext = strrchr(ent->d_name, '.');
			if (!ext){
				continue;
			}
			else if (strncmp((ext + 1), "txt", 3) == 0){		// makes sure you only try to read .txt files
				WUsummary_name = strdup(ent->d_name);
				printf("...found.");
				return WUsummary_name;
			}
		}
	}
	printf("...not found.");
	return NULL;
#endif

}

char **select_WUs(char *WUsummary_name)
{
	printf("\nSelecting work units...		");

	char **all_wu_paths, WUsummary_path[256];
	all_wu_paths = malloc(100000 * sizeof(char*));		// Max number of WUs is 100,000
	if (all_wu_paths == NULL){
		printf("...failed.");
		printf("\n\nerror: could not allocate all_wu_paths array (OOM)");
		return NULL;
	}

	FILE *WUsummary_file;
	strcpy(WUsummary_path, 256, cometwu_path);
	strcat(WUsummary_path, 256, OS_SEP);
	strcat(WUsummary_path, 256, WUsummary_name);
	fopen_s(&WUsummary_file, WUsummary_path, "r");
	if (WUsummary_file == NULL){
		printf("...failed.");
		printf("\n\nerror:  could not open WU summary file; may be open in another program");
		return NULL;
	}
	else{
		char temp[512], *next_token = NULL;
		float WU_mass = 0;
		double WU_times[2];
		int wu_missing_count = 0, WU_number_of_particles = 0, first_line = 1;
		while ((fgets(temp, sizeof(temp), WUsummary_file)) != NULL){
			if (first_line == 1){
				first_line = 0;						//Skips first line in the WUsummary
				continue;
			}
			all_wu_count++;
			char* cval = strtok_r(temp, "\t", &next_token);
			char* WU_name = cval;
			cval = strtok_r(NULL, "\t", &next_token);
			sscanf(cval, "%d", &WU_number_of_particles);
			cval = strtok_r(NULL, "\t", &next_token);
			sscanf(cval, "%f", &WU_mass);
			cval = strtok_r(NULL, "\t", &next_token);
			sscanf(cval, "%lf", &WU_times[0]);
			cval = strtok_r(NULL, "\n", &next_token);
			sscanf(cval, "%lf", &WU_times[1]);
			if (WU_times[0] > timespec[0]){			// Particles must have been created before target time
				continue;
			}
			if (tminflag == 1){
				if (WU_times[1] < timespec[1]){		// Particles must have been created after tmin
					continue;
				}
			}
			if (tmaxflag == 1){
				if (WU_times[0] > timespec[2]){		// Particles must have been created before tmax
					continue;
				}
			}
			if (mminflag == 1){
				if (WU_mass < massspec[0]){		// Particles must be more massive than mmin
					continue;
				}
			}
			if (mmaxflag == 1){
				if (WU_mass > massspec[1]){		// Particles must be less massive than mmax
					continue;
				}
			}	
			FILE *ftest;
			char full_wu_path[256];
			strcpy(full_wu_path, 256, cometwu_path);
			strcat(full_wu_path, 256, OS_SEP);
			strcat(full_wu_path, 256, WU_name);
			strcat(full_wu_path, 256, ".ctwu");
			fopen_s(&ftest, full_wu_path, "rb");
			if (ftest == NULL){
				wu_missing_count++;
				continue;
			}
			fclose(ftest);
			all_wu_paths[wu_count] = strdup(full_wu_path);
			wu_count++;
			all_particles_count += WU_number_of_particles;
		}
		fclose(WUsummary_file);
		printf("...done.\n	The summary lists %d WUs, %d of which are relevant.", all_wu_count, (wu_count + wu_missing_count));
		if (wu_missing_count != 0){
			printf("\n	%d of the selected WUs are missing.", wu_missing_count);
		}
		if (wu_count == 0){
			printf("	No WUs to process; quiting.");
			return all_wu_paths;
		}
		return all_wu_paths;
	}
}

char **search_WUs(void)
{
	//This function searches the INPUT path for all files named *.ctwu
	// It determines what WUs are relevant and saves their names	
	printf("\nSearching for work units...		");

	char **all_wu_paths;
	all_wu_paths = malloc(100000 * sizeof(char*));		// Max number of WUs is 100,000
	if (all_wu_paths == NULL){
		printf("...failed.");
		printf("\n\nerror: could not allocate all_wu_paths array (OOM)");
		return NULL;
	}


#ifdef _WIN32					// Get paths of WUs under windows systems
	HANDLE hFind;
	WIN32_FIND_DATA FindData;
	char findpath[256];
	strcpy(findpath, 256, cometwu_path);
	strcat(findpath, 256, OS_SEP);
	strcat(findpath, 256, "*.ctwu");		// makes sure you only try to read .ctwu files
	hFind = FindFirstFile(findpath, &FindData);
	if (hFind == INVALID_HANDLE_VALUE){
		printf("...failed.");
		printf("\n\nerror: could not find work units in %s", cometwu_path);
		free(all_wu_paths);
		return NULL;
	}
	all_wu_paths[0] = strdup(FindData.cFileName);
	all_wu_count++;
	while (FindNextFile(hFind, &FindData)){
		all_wu_paths[all_wu_count] = strdup(FindData.cFileName);
		all_wu_count++;
	}
	FindClose(hFind);
#else							// Get paths of WUs under linux
	DIR *dir;
	char *ext;
	struct dirent *ent;
	dir = opendir(cometwu_path);
	if (dir){
		while ((ent = readdir(dir)) != NULL){
			ext = strrchr(ent->d_name, '.');
			if (!ext){
				continue;
			}
			else if (strncmp((ext + 1), "ctw", 3) == 0){		// makes sure you only try to read .ctwu files
				all_wu_paths[all_wu_count] = strdup(ent->d_name);
				all_wu_count++;
			}
		}
	}
	else{
		printf("...failed.");
		printf("\n\nerror: could not find directory %s", cometwu_path);
		free(all_wu_paths);
		return NULL;
	}

	if (all_wu_count == 0){
		printf("...failed.");
		printf("\n\nerror: could not find work units in %s", cometwu_path);
		free(all_wu_paths);
		return NULL;
	}
	closedir(dir);
#endif	

	int i, wu_fails = 0;
	char **wu_paths;
	wu_paths = malloc(all_wu_count * sizeof(char*));
	if (wu_paths == NULL){
		printf("...failed.");
		printf("\n\nerror: could not allocate wu_paths array (OOM)");
		for (i = 0; i < all_wu_count; i++){
			free(all_wu_paths[i]);
		}
		free(all_wu_paths);
		return NULL;
	}
	//Check WUs' relevance in parallel
#pragma omp parallel
	{
		FILE *finput;
		float wu_header[7];

#pragma omp for
		for (i = 0; i < all_wu_count; i++){

			//Open WU and read header line
			char full_wu_path[256];
			strcpy(full_wu_path, 256, cometwu_path);
			strcat(full_wu_path, 256, OS_SEP);
			strcat(full_wu_path, 256, all_wu_paths[i]);
			fopen_s(&finput, full_wu_path, "rb");
			if (finput == NULL){
#pragma omp atomic
				wu_fails++;
				continue;
			}
			fread(wu_header, sizeof(float), 7, finput);
			fclose(finput);

			//Check relevance
			if (wu_header[5] > timespec[0]){			// Particles must have been created before target time
				continue;
			}
			if (tminflag == 1){
				if (wu_header[6] < timespec[1]){		// Particles must have been created after tmin
					continue;
				}
			}
			if (tmaxflag == 1){
				if (wu_header[5] > timespec[2]){		// Particles must have been created before tmax
					continue;
				}
			}
			if (mminflag == 1){
				if (wu_header[2] < massspec[0]){		// Particles must be more massive than mmin
					continue;
				}
			}
			if (mmaxflag == 1){
				if (wu_header[2] > massspec[1]){		// Particles must be less massive than mmax
					continue;
				}
			}
#pragma omp critical(WRITE_WU_PATHS)
			{
				wu_paths[wu_count] = strdup(full_wu_path);
				wu_count++;
			}
#pragma omp atomic
			all_particles_count += (int)(wu_header[1] - wu_header[0] + 1.5);
		}
	}
	for (i = 0; i < all_wu_count; i++){
		free(all_wu_paths[i]);
	}
	free(all_wu_paths);
	if (wu_fails != 0){
		printf("...failed.");
		printf("\n\nerror: %d work units could not be read", wu_fails);
		for (i = 0; i < wu_count; i++){
			free(wu_paths[i]);
		}
		free(wu_paths);
		return NULL;
	}

	//fcloseall();		//Causes a Spice bug when calling kernels. Possibly closes kernel files without Spice knowing?
	printf("...done.\n	%d WUs were found, %d of which are relevant.", all_wu_count, wu_count);
	if (wu_count == 0){
		printf("	No WUs to process; quiting.");
		return wu_paths;
	}

	return wu_paths;
}

float *get_particles(char **wu_paths)
{
	//This funtcion looks up relevant particles within the WUs found by **get_WU_paths,
	// gets the properties of these particles and the nearest ouput state

	printf("\n\nLoading particles...		");
	float *nearest_states;
	nearest_states = malloc(all_particles_count * 12 * sizeof(float));
	if (nearest_states == NULL){
		printf("...failed to load particles.\n\nerror: could not allocate nearest_states array (OOM)");
		return nearest_states;
	}
	SpiceInt dim;
	int i, wu_fails = 0;
	bodvcd_c(10, "GM", 1, &dim, &Sun_GM);				//Get standard gravitational parameter of the Sun

	
#pragma omp parallel num_threads(1)
	{
		FILE *finput;
		int wu_rows, fsize, corrupt_count;

		//Loop over all relevant WUs in parallel
#pragma omp for
		for (i = 0; i < wu_count; i++)
		{
			corrupt_count = 0;
			fopen_s(&finput, wu_paths[i], "rb");		//Open WU
			if (finput == NULL){
#pragma omp atomic
				wu_fails++;
				continue;
			}
			
			fseek(finput, 0, SEEK_END);					//Check WU number of rows to determine required array size
			fsize = ftell(finput);
			wu_rows = fsize / (7 * sizeof(float));
			rewind(finput);
			float *fwu_array;
			fwu_array = calloc((wu_rows + 1) * 7, sizeof(float));
			if (fwu_array == NULL){
#pragma omp atomic
				wu_fails++;
				continue;
			}

			fread(fwu_array, sizeof(float), wu_rows * 7, finput);		//Load WU into array
			fclose(finput);
			int j = 1, z, k, wu_particles_total, wu_particle_count = 0;
			wu_particles_total = (int)(fwu_array[1] - fwu_array[0] + 1.5);	//Check number of particles in WU
			float *wu_nearest_states;
			wu_nearest_states = malloc(wu_particles_total * 12 * sizeof(float));	//Create array to save the nearest particle states in
			if (wu_nearest_states == NULL){
#pragma omp atomic
				wu_fails++;
				continue;
			}

			//Find particles in ouput file and get nearest state
			while (j < wu_rows){
				if (fwu_array[j * 7] == 0){
					z = j + 1;
					if (tminflag == 1){
						if (fwu_array[z * 7 + 6] < timespec[1]) {
							j++;
							continue;
						}
					}
					if (tmaxflag == 1){
						if (fwu_array[z * 7 + 6] > timespec[2]) {
							j++;
							continue;
						}
					}
					if (fwu_array[z * 7 + 6] > timespec[0]){
						j++;
						continue;
					}
					z = get_nearest_state(z, fwu_array);
					if (z == j){
						j++;
						continue;
					}
					else if (z == -1){		// particle has 0/0/0 coordinates -> corrupt
						j++;
						corrupt_count++;
						continue;
					}
					else if (z == -2){		// particle fell into the sun and will be ignored
						j++;
						continue;
					}
					if (nodeflag == 1){
						z = get_state_nearest_to_orbital_plane(z, fwu_array);
					}

					//Write particle to private array with following formatting:
					// 0-2 Position  ;  3-5 Velocity  ;  6 time  ;  7 particle number  ;  8 Multiplication factor  ;  9 mass  ;  10 beta  ;  11 time or origin    
					for (k = 0; k < 7; k++){
						wu_nearest_states[wu_particle_count * 12 + k] = fwu_array[z * 7 + k];
					}
					wu_nearest_states[wu_particle_count * 12 + 7] = fwu_array[j * 7 + 4];
					wu_nearest_states[wu_particle_count * 12 + 8] = fwu_array[j * 7 + 5];
					wu_nearest_states[wu_particle_count * 12 + 9] = fwu_array[2];
					wu_nearest_states[wu_particle_count * 12 + 10] = fwu_array[4];
					wu_nearest_states[wu_particle_count * 12 + 11] = fwu_array[(j+1) * 7 + 6];
					wu_particle_count++;
					j = z;
				}
				j++;
			}


			//Alert if WU might be corrupt
			if (corrupt_count != 0){
				printf("\n %d particles left out because they have 0/0/0 coordinates. WU name:\n %s\n			", corrupt_count, wu_paths[i]);
			}

			//Write particles of current WU to shared array
#pragma omp critical(WRITE_NEAREST)
			{
				for (j = 0; j < wu_particle_count; j++){
					for (k = 0; k < 12; k++){
						nearest_states[(particles_count + j) * 12 + k] = wu_nearest_states[j * 12 + k];
					}
				}
				particles_count += wu_particle_count;
			}
			free(wu_nearest_states);
			free(fwu_array);
		}
	}

	for (i = 0; i < wu_count; i++){
		free(wu_paths[i]);
	}
	free(wu_paths);

	if (wu_fails != 0){
		printf("...failed to load particles.");
		printf("\n\nerror: %d work units failed to be processed", wu_fails);
		return NULL;
	}
	else{
		if (particles_count != all_particles_count){
			if (particles_count == 0){				//This can occur when no particle is within specified distance
				printf("...no relevant particles found (%d).\n", all_particles_count);
				return nearest_states;
			}
			float *temp_nearest_states;
			temp_nearest_states = realloc(nearest_states, particles_count * 12 * sizeof(float));
			if (temp_nearest_states == NULL){
				printf("...failed to load particles.");
				printf("\n\nerror: could not reallocate nearest_states array (OOM)");
				free(nearest_states);
				return NULL;
			}
			else{
				nearest_states = temp_nearest_states;
			}
		}
		printf("...done. %d relevant particles loaded.", particles_count);
	}
	//fcloseall();		//Causes a Spice bug when calling kernels. Possibly closes kernel files without Spice knowing?

	return nearest_states;
}

int get_nearest_state(int z, float *fwu_array)
{
	// This function searches the particle output for the state closest to the target time

	int i = z;
	while (fwu_array[(i + 1) * 7 + 6] < timespec[0]){
		if (fwu_array[(i + 1) * 7 + 6] == 0){
			break;
		}
		i++;
	}

	if (fwu_array[i * 7 + 1] == 99){
		return -2;							// This particle was terminated because it fell into the sun. Particle will be ignored.
	}

	if (fwu_array[i * 7] == 0){
		if (fwu_array[i * 7 + 1] == 0){
			if (fwu_array[i * 7 + 2] == 0){
				return -1;						// SPICE can't process particles with coordinates 0/0/0, WU might be corrupt
			}
		}
	}
	
	// Particles too distant to body are filtered out (only for nodes)
	if (nodeflag == 1 || (intgflag == 1 && objectflag == 1)){
		int k;
		double orb_elts[8], state[6], etime, beta, Sun_GM_beta, rel_pos[3], distance;
		beta = (double)fwu_array[4];
		Sun_GM_beta = Sun_GM * (1 - beta);
		for (k = 0; k < 6; k++){
			state[k] = (double)fwu_array[i * 7 + k];
		}
		etime = (double)fwu_array[i * 7 + 6];
		if (suncflag == 1){
			for (k = 0; k < 6; k++){
				state[k] -= Sun_state[k];
			}
		}
#pragma omp critical(SPICE)
		{
			oscelt_c(state, etime, Sun_GM_beta, orb_elts);
			conics_c(orb_elts, timespec[0], state);
		}
		rel_pos[0] = state[0] - object_state[0];
		rel_pos[1] = state[1] - object_state[1];
		rel_pos[2] = state[2] - object_state[2];
		if (suncflag == 1){
			for (k = 0; k < 6; k++){
				state[k] += Sun_state[k];
			}
		}
		distance = sqrt(rel_pos[0] * rel_pos[0] + rel_pos[1] * rel_pos[1] + rel_pos[2] * rel_pos[2]);
		if (distance > max_distance){
			return (z - 1);
		}
	}
	return i;
}

int get_state_nearest_to_orbital_plane(int z, float *fwu_array)
{
	//This funtion searches the output for the nearest state of a particle to the orbital plane
	// This state is later used to compute the node

	int i;
	double plane_distance[3];

	for (i = 0; i < 100; i++){
		if (i == 100){
			printf("\nwarning: could not find closest state to orbital plane with 100 iterations.\n		Particle inclination may be very low. Consider lowering max distance\n");
		}
		plane_distance[0] = fabs(fwu_array[z * 7 + 0] * node_orbitplane_normvec[0] + fwu_array[z * 7 + 1] * node_orbitplane_normvec[1] + fwu_array[z * 7 + 2] * node_orbitplane_normvec[2]);
		if (fwu_array[(z - 1) * 7 + 6] != 0){
			plane_distance[1] = fabs(fwu_array[(z - 1) * 7 + 0] * node_orbitplane_normvec[0] + fwu_array[(z - 1) * 7 + 1] * node_orbitplane_normvec[1] + fwu_array[(z - 1) * 7 + 2] * node_orbitplane_normvec[2]);
			if (plane_distance[1] < plane_distance[0]){
				z--;
				continue;
			}
		}
		if (fwu_array[(z + 1) * 7 + 6] != 0){
			plane_distance[2] = fabs(fwu_array[(z + 1) * 7 + 0] * node_orbitplane_normvec[0] + fwu_array[(z + 1) * 7 + 1] * node_orbitplane_normvec[1] + fwu_array[(z + 1) * 7 + 2] * node_orbitplane_normvec[2]);
			if (plane_distance[2] < plane_distance[0]){
				z++;
				continue;
			}
		}
		break;
	}
	
	return z;
}

float *get_all_target_states(float *nearest_states)
{	
	// This function administrates the calculation of all relevant particles in parallel

	if (intgflag == 0){
		printf("\nCalculating target states...	");
	}
	else{
		printf("\nIntegrating target states...	");
	}

	int i;
#pragma omp parallel
	{
		double nstate[7], beta, Sun_GM_beta;
		int row, k;

#pragma omp for
		for (i = 0; i < particles_count; i++){
			row = i * 12;
			beta = (double) nearest_states[row + 10];
			Sun_GM_beta = Sun_GM*(1 - beta);
			for (k = 0; k < 7; k++){
				nstate[k] = (double) nearest_states[row + k];
			}
			if (intgflag == 0){
				calculate_target_state(nstate, Sun_GM_beta);
			}
			else{
				integrate_target_state(nstate, Sun_GM_beta);
			}		
			for (k = 0; k < 7; k++){
				nearest_states[row + k] = (float) nstate[k];
			}
		}

	}

	printf("...done.");
	return nearest_states;
}

void calculate_target_state(double *nstate, double Sun_GM_beta)
{
	//This function calculates the particle state at the target time or the partcile node
	// (depending on what has been set) based on the nearest output state

	int k;
	double etime, orb_elts[8], state[6];
	for (k = 0; k < 6; k++){
		state[k] = nstate[k];
	}
	etime = nstate[6];

	if (suncflag == 1){
		for (k = 0; k < 6; k++){
			state[k] -= Sun_state[k];
		}
	}

#pragma omp critical(SPICE)
	{
		oscelt_c(state, etime, Sun_GM_beta, orb_elts);		//Convert state to orbital elements
	}

	// If nodeflag is not set, calculate particle positions at target time using KEPLER orbits computed with SPICE
	if (nodeflag == 0){
			
#pragma omp critical(SPICE)
		{
			conics_c(orb_elts, timespec[0], state);			//Convert orbital elements to state at target time
		}

	}
	else{		// If nodeflag is set, calculate the point at which the particle intersects the orbital plane (numerically)
				// -> converge towards orbital plane the distance is less than 100km. Then save the state.
	
		double plane_distance, particle_speed, angle_vel_norm, dt;

		plane_distance = state[0] * node_orbitplane_normvec[0] + state[1] * node_orbitplane_normvec[1] + state[2] * node_orbitplane_normvec[2];
		while ( fabs(plane_distance) > 100){
			plane_distance = state[0] * node_orbitplane_normvec[0] + state[1] * node_orbitplane_normvec[1] + state[2] * node_orbitplane_normvec[2];
			particle_speed = sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
			dt = fabs(plane_distance) / particle_speed;
			angle_vel_norm = acos((state[3] * node_orbitplane_normvec[0] + state[4] * node_orbitplane_normvec[1] + state[5] * node_orbitplane_normvec[2]) / particle_speed);

			//Decide wether the particle at the current state is moving towards or away from the orbital plane.
			//If it is moving away -> choose negative time step
			if ((plane_distance > 0 && angle_vel_norm < PI_half) || (plane_distance < 0 && angle_vel_norm > PI_half)){
				dt = -dt;
			}
			etime += dt;
#pragma omp critical(SPICE)
			{
				conics_c(orb_elts, etime, state);
			}
		}	

		nstate[6] = etime;
	}

	if (suncflag == 1){
		for (k = 0; k < 6; k++){
			state[k] += Sun_state[k];
		}
	}
	for (k = 0; k < 6; k++){
		nstate[k] = state[k];
	}
}

void integrate_target_state(double *nstate, double Sun_GM_beta)
{
	//This function integrates the particle state at the target time
	// (depending on what has been set) based on the nearest output state

	int j;
	double lt, dt, dt2, abs_acc, pos[3];
	double **body_pre, **body_mid, **body_end;
	double k_acc_1[3], k_acc_2[3], k_acc_3[3], k_acc_4[3];
	double k_vel_1[3], k_vel_2[3], k_vel_3[3], k_vel_4[3];

	body_pre = malloc(N_bodys * sizeof(int *));
	body_mid = malloc(N_bodys * sizeof(int *));
	body_end = malloc(N_bodys * sizeof(int *));
	if (body_pre == NULL || body_mid == NULL || body_end == NULL)
	{
		printf("\nerror: could not allocate body state array (OOM)");
		//Abort...
	}
	for (j = 0; j < N_bodys; j++)
	{
		body_pre[j] = malloc(3 * sizeof(double));
		body_mid[j] = malloc(3 * sizeof(double));
		body_end[j] = malloc(3 * sizeof(double));
		if (body_pre[j] == NULL || body_mid[j] == NULL || body_end[j] == NULL)
		{
			printf("\n\nerror: could not allocate body state array (OOM)");
			//Abort...
		}
#pragma omp critical(SPICE)
		{
			//Critical section is only executed on one thread at a time (spice is not threadsafe)
			spkezp_c(body_ID[j], nstate[6], "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
		}
	}

	//Integrate
	if (nodeflag == 0){
		while (nstate[6] < timespec[0])
		{
			for (j = 0; j < N_bodys; j++)
			{
				body_pre[j][0] = body_end[j][0];
				body_pre[j][1] = body_end[j][1];
				body_pre[j][2] = body_end[j][2];
			}

			//Step 1
			pos[0] = nstate[0];
			pos[1] = nstate[1];
			pos[2] = nstate[2];
			calc_accel(N_bodys, Sun_GM_beta, nstate, &body_pre, k_acc_1);
			k_vel_1[0] = nstate[3];
			k_vel_1[1] = nstate[4];
			k_vel_1[2] = nstate[5];

			//Set dynamic step size
			abs_acc = sqrt(k_acc_1[0] * k_acc_1[0] + k_acc_1[1] * k_acc_1[1] + k_acc_1[2] * k_acc_1[2]);
			dt = (dv_step / abs_acc);
			if (nstate[6] + dt >  timespec[0])
			{
				dt = fabs(timespec[0] - nstate[6]);			//End on final_time exactly
			}
			dt2 = dt / 2;

			//Get body positions with SPICE
			for (j = 0; j < N_bodys; j++)
			{
#pragma omp critical(SPICE)
				{
					//Critical section is only executed on one thread at a time (not thread-safe)
					spkezp_c(body_ID[j], nstate[6] + dt, "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
				}
				body_mid[j][0] = (body_pre[j][0] + body_end[j][0]) / 2;
				body_mid[j][1] = (body_pre[j][1] + body_end[j][1]) / 2;
				body_mid[j][2] = (body_pre[j][2] + body_end[j][2]) / 2;
			}

			//Step 2
			pos[0] = nstate[0] + k_vel_1[0] * dt2;
			pos[1] = nstate[1] + k_vel_1[1] * dt2;
			pos[2] = nstate[2] + k_vel_1[2] * dt2;
			calc_accel(N_bodys, Sun_GM_beta, pos, &body_mid, k_acc_2);
			k_vel_2[0] = nstate[3] + k_acc_1[0] * dt2;
			k_vel_2[1] = nstate[4] + k_acc_1[1] * dt2;
			k_vel_2[2] = nstate[5] + k_acc_1[2] * dt2;

			//Step 3
			pos[0] = nstate[0] + k_vel_2[0] * dt2;
			pos[1] = nstate[1] + k_vel_2[1] * dt2;
			pos[2] = nstate[2] + k_vel_2[2] * dt2;
			calc_accel(N_bodys, Sun_GM_beta, pos, &body_mid, k_acc_3);
			k_vel_3[0] = nstate[3] + k_acc_2[0] * dt2;
			k_vel_3[1] = nstate[4] + k_acc_2[1] * dt2;
			k_vel_3[2] = nstate[5] + k_acc_2[2] * dt2;

			//Step 4
			pos[0] = nstate[0] + k_vel_3[0] * dt;
			pos[1] = nstate[1] + k_vel_3[1] * dt;
			pos[2] = nstate[2] + k_vel_3[2] * dt;
			calc_accel(N_bodys, Sun_GM_beta, pos, &body_end, k_acc_4);
			k_vel_4[0] = nstate[3] + k_acc_3[0] * dt;
			k_vel_4[1] = nstate[4] + k_acc_3[1] * dt;
			k_vel_4[2] = nstate[5] + k_acc_3[2] * dt;

			//Update solution
			nstate[0] = nstate[0] + dt*(k_vel_1[0] + 2 * (k_vel_2[0] + k_vel_3[0]) + k_vel_4[0]) / 6;
			nstate[1] = nstate[1] + dt*(k_vel_1[1] + 2 * (k_vel_2[1] + k_vel_3[1]) + k_vel_4[1]) / 6;
			nstate[2] = nstate[2] + dt*(k_vel_1[2] + 2 * (k_vel_2[2] + k_vel_3[2]) + k_vel_4[2]) / 6;
			nstate[3] = nstate[3] + dt*(k_acc_1[0] + 2 * (k_acc_2[0] + k_acc_3[0]) + k_acc_4[0]) / 6;
			nstate[4] = nstate[4] + dt*(k_acc_1[1] + 2 * (k_acc_2[1] + k_acc_3[1]) + k_acc_4[1]) / 6;
			nstate[5] = nstate[5] + dt*(k_acc_1[2] + 2 * (k_acc_2[2] + k_acc_3[2]) + k_acc_4[2]) / 6;
			nstate[6] = nstate[6] + dt;
		}
	}
	else{
		double plane_distance, particle_speed, particle_speed_perpendicular, angle_vel_norm, dt_max;
		plane_distance = nstate[0] * node_orbitplane_normvec[0] + nstate[1] * node_orbitplane_normvec[1] + nstate[2] * node_orbitplane_normvec[2];
		while (fabs(plane_distance) > 100)
		{
			for (j = 0; j < N_bodys; j++)
			{
				body_pre[j][0] = body_end[j][0];
				body_pre[j][1] = body_end[j][1];
				body_pre[j][2] = body_end[j][2];
			}

			//Step 1
			pos[0] = nstate[0];
			pos[1] = nstate[1];
			pos[2] = nstate[2];
			calc_accel(N_bodys, Sun_GM_beta, nstate, &body_pre, k_acc_1);
			k_vel_1[0] = nstate[3];
			k_vel_1[1] = nstate[4];
			k_vel_1[2] = nstate[5];

			//Set dynamic step size
			abs_acc = sqrt(k_acc_1[0] * k_acc_1[0] + k_acc_1[1] * k_acc_1[1] + k_acc_1[2] * k_acc_1[2]);
			dt = (dv_step / abs_acc);

			particle_speed = sqrt(nstate[3] * nstate[3] + nstate[4] * nstate[4] + nstate[5] * nstate[5]);
			angle_vel_norm = acos((nstate[3] * node_orbitplane_normvec[0] + nstate[4] * node_orbitplane_normvec[1] + nstate[5] * node_orbitplane_normvec[2]) / particle_speed);
			if ((plane_distance > 0 && angle_vel_norm < PI_half) || (plane_distance < 0 && angle_vel_norm > PI_half)){
				break;
			}
			particle_speed_perpendicular = fabs(cos(angle_vel_norm)) * particle_speed;
			dt_max = 0.9 * plane_distance / particle_speed_perpendicular;
			if (dt >  dt_max)
			{
				dt = dt_max;			//Don't overshoot orbital plane
			}
			dt2 = dt / 2;

			//Get body positions with SPICE
			for (j = 0; j < N_bodys; j++)
			{
#pragma omp critical(SPICE)
				{
					//Critical section is only executed on one thread at a time (not thread-safe)
					spkezp_c(body_ID[j], nstate[6] + dt, "ECLIPJ2000", "NONE", 0, body_end[j], &lt);
				}
				body_mid[j][0] = (body_pre[j][0] + body_end[j][0]) / 2;
				body_mid[j][1] = (body_pre[j][1] + body_end[j][1]) / 2;
				body_mid[j][2] = (body_pre[j][2] + body_end[j][2]) / 2;
			}

			//Step 2
			pos[0] = nstate[0] + k_vel_1[0] * dt2;
			pos[1] = nstate[1] + k_vel_1[1] * dt2;
			pos[2] = nstate[2] + k_vel_1[2] * dt2;
			calc_accel(N_bodys, Sun_GM_beta, pos, &body_mid, k_acc_2);
			k_vel_2[0] = nstate[3] + k_acc_1[0] * dt2;
			k_vel_2[1] = nstate[4] + k_acc_1[1] * dt2;
			k_vel_2[2] = nstate[5] + k_acc_1[2] * dt2;

			//Step 3
			pos[0] = nstate[0] + k_vel_2[0] * dt2;
			pos[1] = nstate[1] + k_vel_2[1] * dt2;
			pos[2] = nstate[2] + k_vel_2[2] * dt2;
			calc_accel(N_bodys, Sun_GM_beta, pos, &body_mid, k_acc_3);
			k_vel_3[0] = nstate[3] + k_acc_2[0] * dt2;
			k_vel_3[1] = nstate[4] + k_acc_2[1] * dt2;
			k_vel_3[2] = nstate[5] + k_acc_2[2] * dt2;

			//Step 4
			pos[0] = nstate[0] + k_vel_3[0] * dt;
			pos[1] = nstate[1] + k_vel_3[1] * dt;
			pos[2] = nstate[2] + k_vel_3[2] * dt;
			calc_accel(N_bodys, Sun_GM_beta, pos, &body_end, k_acc_4);
			k_vel_4[0] = nstate[3] + k_acc_3[0] * dt;
			k_vel_4[1] = nstate[4] + k_acc_3[1] * dt;
			k_vel_4[2] = nstate[5] + k_acc_3[2] * dt;

			//Update solution
			nstate[0] = nstate[0] + dt*(k_vel_1[0] + 2 * (k_vel_2[0] + k_vel_3[0]) + k_vel_4[0]) / 6;
			nstate[1] = nstate[1] + dt*(k_vel_1[1] + 2 * (k_vel_2[1] + k_vel_3[1]) + k_vel_4[1]) / 6;
			nstate[2] = nstate[2] + dt*(k_vel_1[2] + 2 * (k_vel_2[2] + k_vel_3[2]) + k_vel_4[2]) / 6;
			nstate[3] = nstate[3] + dt*(k_acc_1[0] + 2 * (k_acc_2[0] + k_acc_3[0]) + k_acc_4[0]) / 6;
			nstate[4] = nstate[4] + dt*(k_acc_1[1] + 2 * (k_acc_2[1] + k_acc_3[1]) + k_acc_4[1]) / 6;
			nstate[5] = nstate[5] + dt*(k_acc_1[2] + 2 * (k_acc_2[2] + k_acc_3[2]) + k_acc_4[2]) / 6;
			nstate[6] = nstate[6] + dt;

			plane_distance = nstate[0] * node_orbitplane_normvec[0] + nstate[1] * node_orbitplane_normvec[1] + nstate[2] * node_orbitplane_normvec[2];
		}
	}

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
}

void calc_accel(int N_bodys, double Sun_GM_beta, double pos[], double **body_state[], double *accel)
{
	double direct_body[3], distance_pow3, GMr3;
	int b;

	//Sun
	direct_body[0] = (*body_state)[0][0] - pos[0];		
	direct_body[1] = (*body_state)[0][1] - pos[1];
	direct_body[2] = (*body_state)[0][2] - pos[2];
	distance_pow3 = pow(direct_body[0] * direct_body[0] + direct_body[1] * direct_body[1] + direct_body[2] * direct_body[2], 1.5);
	GMr3 = Sun_GM_beta / distance_pow3;
	accel[0] = GMr3 * direct_body[0];
	accel[1] = GMr3 * direct_body[1];
	accel[2] = GMr3 * direct_body[2];

	//Other bodys
	for (b = 1; b < N_bodys; b++)		
	{
		direct_body[0] = (*body_state)[b][0] - pos[0];
		direct_body[1] = (*body_state)[b][1] - pos[1];
		direct_body[2] = (*body_state)[b][2] - pos[2];
		distance_pow3 = pow(direct_body[0] * direct_body[0] + direct_body[1] * direct_body[1] + direct_body[2] * direct_body[2], 1.5);
		GMr3 = body_GM[b] / distance_pow3;
		accel[0] += GMr3 * direct_body[0];
		accel[1] += GMr3 * direct_body[1];
		accel[2] += GMr3 * direct_body[2];
	}
}

float *filter_particles_out_of_range(float *target_states)
{
	// This function filters all the particles that are not within the specified distance to the object

	printf("\nCheck distances to object...	");
	float *filtered_states;
	filtered_states = malloc(particles_count * 12 * sizeof(float));
	if (filtered_states == NULL){
		printf("...failed to filter particles.\n\nerror: could not allocate filtered_states array (OOM)");
		free(target_states);
		return NULL;
	}
	int i, filtered_particles_count = 0;

#pragma omp parallel
	{
		double state[6], rel_state[6], distance;
		int row, k;

#pragma omp for
		for (i = 0; i < particles_count; i++){

			row = i * 12;
			for (k = 0; k < 6; k++){
				state[k] = (double)target_states[row + k];
			}

			for (k = 0; k < 6; k++){								//Calculate coordinates relative to object
				rel_state[k] = state[k] - object_state[k];
			}		

			distance = sqrt(rel_state[0] * rel_state[0] + rel_state[1] * rel_state[1] + rel_state[2] * rel_state[2]);
			if (distance < max_distance){
#pragma omp critical(WRITE_FILTERED)			//If particle is within specified range, save to filtered_states array
				{
					for (k = 0; k < 6; k++){
						filtered_states[filtered_particles_count * 12 + k] = (float)rel_state[k];
					}
					for (k = 6; k < 12; k++){
						filtered_states[filtered_particles_count * 12 + k] = target_states[row + k];
					}
					filtered_particles_count++;
				}
			}
		}
	}

	free(target_states);
	particles_count = filtered_particles_count;

	printf("...done. %d particles within distance.", particles_count);
	return filtered_states;
}

int save_target_states(float *target_states)
{
	//This function saves the results in an ouput file, of which the path has been defined in the application input arguments

	printf("\nWriting output...		");
	FILE *outputfile;
	fopen_s(&outputfile, output_path, "wb");
	if (outputfile == NULL){
		printf("...failed.\n\nerror: could not create output file");
		return 1;
	}

	if (pc2flag == 1){					//Produce .PC2 format output (pointCache)		//NOT TESTED!

		float *output_array;
		int i, k;

		char cacheSignature[12] = "POINTCACHE2";
		int fileVersion = 1;
		int numPoints = particles_count;
		float startFrame = 1;
		float sampleRate = 1;
		int numSamples = 1;

		output_array = malloc(particles_count * 3 * sizeof(float));
		if (output_array == NULL){
			printf("...failed.\n\nerror: could not allocate output_array (OOM)");
			fclose(outputfile);
			return 1;
		}
		for (i = 0; i < particles_count; i++){
			for (k = 0; k < 3; k++){
				output_array[i * 3 + k] = target_states[i * 12 + k];
			}
		}

		fwrite(cacheSignature, sizeof(char), 12, outputfile);
		fwrite(&fileVersion, sizeof(int), 1, outputfile);
		fwrite(&numPoints, sizeof(int), 1, outputfile);
		fwrite(&startFrame, sizeof(float), 1, outputfile);
		fwrite(&sampleRate, sizeof(float), 1, outputfile);
		fwrite(&numSamples, sizeof(int), 1, outputfile);

		fwrite(output_array, sizeof(float), particles_count * 3, outputfile);
		
		free(output_array);
	}
	else if (nodeflag == 1){			//If nodeflag is set, then target_states is already properly formatted and ready for saving
										// 3 Position and 3 Velocity coordinates  ;  time  ;  particle number  ;  multiplication-factor  ;  mass  ;  beta  ;  time of origin 
		fwrite(target_states, sizeof(float), particles_count * 12, outputfile);

	}
	else{								//Else take out time and potentially velocity coordinates

		float *output_array;
		int i, k;

		if (velflag == 1){		//Keep velocities in output file
			output_array = malloc(particles_count * 11 * sizeof(float));
			if (output_array == NULL){
				printf("...failed.\n\nerror: could not allocate output array (OOM)");
				fclose(outputfile);
				return 1;
			}
			for (i = 0; i < particles_count; i++){
				for (k = 0; k < 6; k++){							
										// 3 Position coordinates and 3 Velocity coordinates
					output_array[i*11 + k] = target_states[i*12 + k];
				}
				for (k = 6; k < 11; k++){								
										// particle number  ;  multiplication-factor  ;  mass  ;  beta  ;  time of origin
					output_array[i*11 + k] = target_states[i*12 + k + 1];
				}
			}
			fwrite(output_array, sizeof(float), particles_count * 11, outputfile);		//Write to file
		}
		else{					//Throw out velocities in output file
			output_array = malloc(particles_count * 8 * sizeof(float));
			if (output_array == NULL){
				printf("...failed.\n\nerror: could not allocate output array (OOM)");
				fclose(outputfile);
				return 1;
			}
			for (i = 0; i < particles_count; i++){
				for (k = 0; k < 3; k++){								
										// 3 Position coordinates
					output_array[i*8 + k] = target_states[i*12 + k];
				}
				for (k = 3; k < 8; k++){
										// particle number  ;  multiplication-factor  ;  mass  ;  beta  ;  time of origin
					output_array[i*8 + k] = target_states[i*12 + k + 4];
				}
			}
			fwrite(output_array, sizeof(float), particles_count * 8, outputfile);		//Write to file
		}

		free(output_array);
	}
	
	fclose(outputfile);

	free(target_states);
	
	printf("...done.");
	return 0;
}

void WUcheck(char *WUsummary_name){
	printf("\nChecking work units...		");

	char **missing_wu_names, WUsummary_path[256];
	missing_wu_names = malloc(100000 * sizeof(char*));		// Max number of WUs is 100,000
	if (missing_wu_names == NULL){
		printf("...failed.");
		printf("\n\nerror: could not allocate all_wu_paths array (OOM)");
	}
	else{
		FILE *WUsummary_file;
		strcpy(WUsummary_path, 256, cometwu_path);
		strcat(WUsummary_path, 256, OS_SEP);
		strcat(WUsummary_path, 256, WUsummary_name);
		fopen_s(&WUsummary_file, WUsummary_path, "r");
		if (WUsummary_file == NULL){
			printf("...failed.");
			printf("\n\nerror:  could not open WU summary file; may be open in another program");
		}
		else{
			char temp[512], *next_token = NULL;
			int wu_missing_count = 0, first_line = 1;
			while ((fgets(temp, sizeof(temp), WUsummary_file)) != NULL){
				if (first_line == 1){
					first_line = 0;						//Skips first line in the WUsummary
					continue;
				}
				all_wu_count++;
				char* cval = strtok_r(temp, "\t", &next_token);
				char* WU_name = cval;
				FILE *ftest;
				char full_wu_path[256];
				strcpy(full_wu_path, 256, cometwu_path);
				strcat(full_wu_path, 256, OS_SEP);
				strcat(full_wu_path, 256, WU_name);
				strcat(full_wu_path, 256, ".ctwu");
				fopen_s(&ftest, full_wu_path, "rb");
				if (ftest == NULL){
					strcat(WU_name, 256, "\n");
					missing_wu_names[wu_missing_count] = strdup(WU_name);
					wu_missing_count++;
					continue;
				}
				fclose(ftest);
			}
			fclose(WUsummary_file);
			if (wu_missing_count == 0){
				printf("...done.\n	The summary lists %d WUs, none of which are missing.\n", all_wu_count);
			}
			else{
				printf("...done.\n	The summary lists %d WUs, %d of which are missing.", all_wu_count, wu_missing_count);
				FILE *outputfile;
				fopen_s(&outputfile, output_path, "w");
				if (outputfile == NULL){
					printf("\n\nerror: could not create output file");
				}
				else{
					int i;
					for (i = 0; i < wu_missing_count; i++){
						fputs(missing_wu_names[i], outputfile);
					}
					fclose(outputfile);
					printf("\n	Names of the missing WUs have been written to:\n	%s\n", output_path);
				}
			}
		}
	}
}
