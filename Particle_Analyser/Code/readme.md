# ParticleAnalyser
**Version 0.01**

The ParticleAnalyser is a tool to read a database created by the ParticleIntegrator.
It retrieves the states of all particles of a comet stream at any point in time.
Furthermore, you can specify particle properties like age or mass.
In normal mode it returns all particles of a stream in barycentric coordinates, but
can also be used to compute nodes in an orbital plane or coordinates relative to an 
object.

##Instructions

-	input arguments:

		-	required:
			-	"--input COMETPATH" 	where COMETPATH is the path to the folder containing all a comet's WUs.
			-	"--output OUTPUTPATH" 	where OUTPUTPATH is the path/file where the output particle file is stored.
			- 	"--time TIME"			where TIME is the requested time at which particle positions will be computed.
										format is ephermeris time [s] (before/past J2000).
										
		-	optional:
			-	"--tmin TMIN"			Particles must have formed after TMIN.
			-	"--tmax TMAX"			Particles must have formed before TMAX.
			-	"--mmin TMAX"			Particles must have a mass higher than MMIN.
			-	"--mmax TMAX"			Particles must have a mass lower than MMAX.
			
			-	"--object OBJECT_ID MAX_DISTANCE"		Specify an object with a Spice-conform ID, e.g. 'EARTH' ('3') or 'ROSETTA'.
														Specify a distance to the object within that particles must be, e.g. '150000000' (~1AU).
														The output will only contain particles within MAX_DISTANCE and if -n is not set
														position and velocity coordinates will be relative to that of the object.
			
			-	"-s"  		sun_centered: 	particle positions will be converted to Sun centered coordinates before calculating Kepler orbits
			-	"-v"		velocities:		output in normal mode will contain velocity coordinates
			-	"-n"		nodes:			compute particle nodes in the orbital plane of OBJECT_ID (requires --object to be set)
			-	"-p"		pc2-format:		the output file will have .pc2 file format.
			
	If "--object" is set it is necessary to specify the directory of the SPK kernel ,corresponding to OBJECT_ID, in "kernels_spk.txt".
	Also if "-s" is set it is necessary to specify the directory of a Sun SPK kernel (e.g. de430.bsp) in "kernels_spk.txt".
		
			
			
-	executable returns:

    >   0:	the app runs successfully to its end  
    >   1: 	no success  (see console)
	
	
Regards, Max


##Changelog

__0.01__	-29.06.2014-:
-	initial alpha release