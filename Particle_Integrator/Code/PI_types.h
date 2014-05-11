/* PI_types.h - defines types used across functions */

/* struct type to store the processed values from the config file */
typedef struct
{
	int algorithm;
	char inputfpath[260];
	char outputpath[260];
	int number_of_threads;
	int save_as_binary;
	SpiceDouble final_time;
	SpiceDouble start_time_save;
	int N_bodys;
	int body_int[10];
	SpiceDouble GM[10];
	SpiceDouble dv_step;
	SpiceDouble e_target;
	int n;
	int first_particle_number;
	SpiceDouble particle_mass;
	SpiceDouble particle_density;
	SpiceDouble particle_radius;
} configuration_values;

/* struct type for the config file readout */
typedef struct
{
	// Simulation
	const char* algo;
	const char* finaltime;
	const char* starttimes;
	int nbodys;
	const char* bodysid;
	const char* mult;
	int nthreads;
	int savebin;
	// Particles
	const char* inputfn;
	const char* outputfn;
	const char* pmass;
	int pdensity;
	int fpnum;
	// Algorithm-specific
	const char* dvstep;
	const char* etarget;
} configuration_readout;