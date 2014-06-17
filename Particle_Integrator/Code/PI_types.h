/* PI_types.h - defines types used across functions */

/* struct type to store the processed values from the config file */
typedef struct
{
	int algorithm;
	bool ssb_centered;
	char inputfpath[260];
	char outputpath[260];
	int number_of_threads;
	bool save_as_binary;
	bool endontime;
	SpiceDouble final_time;			// [s]
	SpiceDouble start_time_save;	// [s]
	int N_bodys;
	int body_int[10];
	ConstSpiceChar *body_char[10];
	SpiceDouble GM[12];				// [km^3/s^2]
	int n;
	SpiceDouble step_multiplier;	// [1]
	int first_particle_number;
	SpiceDouble particle_mass;		// [kg]
	SpiceDouble particle_density;	// [kg/m^-3]
	SpiceDouble particle_radius;	// [m]
	SpiceDouble q_pr;				// [-]
	SpiceDouble solar_lum;			// [W]
	SpiceDouble beta;				// [-]
	SpiceDouble betaGM;				// [km^3/s^2]
	// Algorithm-specific
	SpiceDouble dv_step;			// [km/s^2]
	SpiceDouble e_target;			// [km]
	int interp_order;
} configuration_values;

/* struct type for the config file readout */
typedef struct
{
	// Simulation
	char *algo;
	int ssbc;
	char *finaltime;
	char *starttimes;
	int nbodys;
	char *bodysid;
	char *mult;
	int nthreads;
	int savebin;
	int endontime;
	// Particles
	char *inputfn;
	char *outputfn;
	char *pmass;
	char *q_pr;
	char *pdensity;
	int fpnum;
	// Algorithm-specific
	char *dvstep;
	char *etarget;
	int iorder;
} configuration_readout;

/* struct type for precomputed dtime powers for interpolation */
typedef struct
{
	SpiceDouble dtime1p2;
	SpiceDouble dtime1p3;
	SpiceDouble dtime1p4;
	SpiceDouble dtime1p5;
	SpiceDouble dtime1p6;

	SpiceDouble dtime8p2;
	SpiceDouble dtime8p3;
	SpiceDouble dtime8p4;
	SpiceDouble dtime8p5;
	SpiceDouble dtime8p6;

	SpiceDouble dtime81;
	SpiceDouble dtime81p2;
} dtimepowers;