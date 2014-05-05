ParticleIntegrator
Version 0.12

-	put your particle file in \\INPUT

-	to change configuration edit configuration.ini in \\INPUT

-	particles are saved to \\OUTPUT

-	executable returns:
	0:	the app runs successfully to its end
	1: 	something is wrong with the input
	2: 	the app wasn't able to access required files or ran out of memory
	

Regards, Max


---Changelog---
0.13	-14.04.2014-:
	-	switched to cmake
	-	linux compatibility
	-	changed library for .ini-file parsing
	-	easily switch timing printouts

0.12	-14.03.2014-:
	-	fixed a bug that led to a wrong evaluation of solar radiation pressure
	
0.11	-06.03.2014-:
	-	introduced the PARTICLE_DENSITY parameter, which sets the density of the material the
		particle is made of (this was always 1000 kg/m3 before this update)

0.10	-28.02.2014-:
	-	several tweaks that should improve performance, especially if run with multithreading
	
0.09	-24.02.2014-:
	-	implemented the START_TIME_SAVE parameter, which sets the point in time at which the app
		starts to save particle states
		
0.08	-03.02.2014-:
	-	added much more error checking that should prevent the app from crashes which are caused by
		the app failing to read from / write to files. Now the app tries to access the file once
		again after 100ms. If this doesn't work even for second time, the executable exits and 
		returns 2
	-	progress.txt is now updated everytime a particle is finished and not on every checkpoint,
		like before (should safe resources and simplify the code)
	
0.07	-27.01.2014-:
	-	fixed bug where the app would crash, if it reads a corrupted incomplete particle output file
	-	turned on optimisation options for the compiler
	
0.06	-20.01.2014-:
	-	fixed low accuracy bug
	-	simplified several sections in the code
	-	changed the SAVE_NTH parameter to a new SAVE_NTH_MULTIPLIER parameter
	
0.05	-07.01.2014-:
	-	solar pressure is now accounted for in the simulation
	
0.04	-19.12.2013-:
	-	application does now support Windows XP OS
	-	added a version for 32bit systems
	-	application does now skip the first five lines in a particle input file,
		which contain some header information about the particle
	-	the output values are now printed with scientific notation
		
0.03	-17.12.2013-:
	-	application does now resume to process incomplete particles
		and starts at the point whereever the particle integration was
		interrupted in the previous run
	-	on startup, if not already existent, progress.txt is created
		where the progress of the current task is saved (from 0.0 to 1.0)
	-	number of digits in the output increased
	
0.02	-13.12.2013-:
	- 	stepsize is now recalculated every step
	- 	added particle mass input (solar pressure is not accounted for as yet)
	- 	added defaults for some input parameters
	-	some minor changes
	
0.01	-05.12.2013-:
	- 	initial alpha release
	- 	use of more than one processor core (/thread) is possible (has to be
		specified in the configfile if wanted)
	- 	if the application is interrupted and started again, it will skip
		particles that have already been completed
		