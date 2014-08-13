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
int pc2flag = 0;


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
double Sun_GM, Sun_pos[3];
double max_distance;


//Functions
int parse_input(int argc, char* argv[]);
char **get_WU_paths(void);
float *get_particles(char **wu_paths);
int get_nearest_state(int z, float *fwu_array);
int get_state_nearest_to_orbital_plane(int z, float *fwu_array);
float *get_all_target_states(float *nearest_states);
void calculate_target_state(double *nstate, double GM);
float *filter_particles_out_of_range(float *target_states);
int save_target_states(float *target_states);


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
	int i = 1, x = 3;
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
			else{
				printf("\n\nerror: input argument %d is invalid", i);
				return 1;
			}
		}
		else if (strncmp(argv[i], "-", 1) == 0)
		{
			if (strspn("s", argv[i]) != 0){				//Specifies wether or not to convert to Sun-centered coordinates
				suncflag = 1;							// before calculating orbital elements
				i++;
			}
			else if (strspn("v", argv[i]) != 0){		//Specifies wether or not velocities will be written to output
				velflag = 1;
				i++;
			}
			else if (strspn("n", argv[i]) != 0){		//Specifies wether or not to compute particle nodes
				nodeflag = 1;
				i++;
			}
			else if (strspn("p", argv[i]) != 0){		//Specifies wether or not to compute particle nodes
				pc2flag = 1;
				i++;
			}
			else{
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

	//Set up object location for relative coordinates computation
	if (objectflag == 1){
		double lt;
		furnsh_c("kernels_spk.txt");	//object kernel directory must be specified in "kernels_spk.txt"
		spkezr_c(object_id, timespec[0], "ECLIPJ2000", "NONE", "SSB", object_state, &lt);

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
		double lt;
		if (objectflag == 0){
			furnsh_c("kernels_spk.txt");	//de430.bsp kernel directory must be specified in "kernels_spk.txt"
		}
		spkezp_c(10, timespec[0], "ECLIPJ2000", "NONE", 0, Sun_pos, &lt);
	}
	
	return 0;
}

char **get_WU_paths(void)
{
	//This function searches the INPUT path for all files named *.ctwu
	// It determines what WUs are relevant and saves their names

	printf("\nSeeking work units...		");
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
			ext = strrchr(ent->d_name,'.');
			if (!ext){
				continue;
			}
			else if (strncmp((ext+1), "ctw", 3) == 0){		// makes sure you only try to read .ctwu files
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
	if (wu_count == 0){
		printf("...no relevant WUs found (%d).", all_wu_count);
		return wu_paths;
	}
	
	fcloseall();
	printf("...done. %d (%d) work units found.", wu_count, all_wu_count);

	return wu_paths;
}

float *get_particles(char **wu_paths)
{
	//This funtcion looks up relevant particles within the WUs found by **get_WU_paths,
	// gets the properties of these particles and the nearest ouput state

	printf("\nLoading particles...		");
	float *nearest_states;
	nearest_states = malloc(all_particles_count * 12 * sizeof(float));
	if (nearest_states == NULL){
		printf("...failed to load particles.\n\nerror: could not allocate nearest_states array (OOM)");
		return nearest_states;
	}
	SpiceInt dim;
	int i, wu_fails = 0;
	bodvcd_c(10, "GM", 1, &dim, &Sun_GM);				//Get standard gravitational parameter of the Sun

	
#pragma omp parallel
	{
		FILE *finput;
		int wu_rows, fsize;

		//Loop over all relevant WUs in parallel
#pragma omp for
		for (i = 0; i < wu_count; i++)
		{
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
			fwu_array = malloc(fsize + 1);
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
		printf("...done. %d (%d) particles loaded.", particles_count, all_particles_count);
	}
	fcloseall();
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
	
	// Particles too distant to body are filtered out (only for nodes)
	if (nodeflag == 1){
		int k;
		double orb_elts[8], state[6], etime, beta, bGM, rel_pos[3], distance;
		beta = (double)fwu_array[4];
		bGM = Sun_GM * (1 - beta);
		for (k = 0; k < 6; k++){
			state[k] = (double)fwu_array[i * 7 + k];
		}
		etime = (double)fwu_array[i * 7 + 6];
#pragma omp critical(SPICE)
		{
			oscelt_c(state, etime, bGM, orb_elts);
			conics_c(orb_elts, timespec[0], state);
		}
		rel_pos[0] = state[0] - object_state[0];
		rel_pos[1] = state[1] - object_state[1];
		rel_pos[2] = state[2] - object_state[2];
		
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

	printf("\nCalculating target states...	");
	int i;

#pragma omp parallel
	{
		double nstate[7], beta, bGM;
		int row, k;

#pragma omp for
		for (i = 0; i < particles_count; i++){
			row = i * 12;
			beta = (double) nearest_states[row + 10];
			bGM = Sun_GM*(1 - beta);
			for (k = 0; k < 7; k++){
				nstate[k] = (double) nearest_states[row + k];
			}
			calculate_target_state(nstate, bGM);	
			for (k = 0; k < 7; k++){
				nearest_states[row + k] = (float) nstate[k];
			}
		}

	}

	printf("...done.");
	return nearest_states;
}

void calculate_target_state(double *nstate, double bGM)
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
		for (k = 0; k < 3; k++){
			state[k] -= Sun_pos[k];
		}
	}

	// If nodeflag is not set, calculate particle positions at target time using KEPLER orbits computed with SPICE
	if (nodeflag == 0){
			
#pragma omp critical(SPICE)
		{
			oscelt_c(state, etime, bGM, orb_elts);		//Convert state to orbital elements
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
#pragma omp critical(SPICE)
			{
				oscelt_c(state, etime, bGM, orb_elts);
				etime += dt;
				conics_c(orb_elts, etime, state);
			}
		}	

		nstate[6] = etime;
	}

	if (suncflag == 1){
		for (k = 0; k < 3; k++){
			state[k] += Sun_pos[k];
		}
	}
	for (k = 0; k < 6; k++){
		nstate[k] = state[k];
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

	printf("...done. %d particles left.", particles_count);
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
