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
										format is ephemeris time [s] (before/past J2000).
										
		-	optional:
			-	"--tmin TMIN"			Particles must have formed after TMIN.
			-	"--tmax TMAX"			Particles must have formed before TMAX.
			-	"--mmin MMIN"			Particles must have a mass higher than MMIN.
			-	"--mmax MMAX"			Particles must have a mass lower than MMAX.
			
			-	"--object OBJECT_ID MAX_DISTANCE"		Specify an object with a Spice-conform ID, e.g. 'EARTH' ('3') or 'ROSETTA'.
														Specify a distance [km] to the object within that particles must be, e.g. '150000000' (~1AU).
														The output will only contain particles within MAX_DISTANCE and if -n is not set
														position and velocity coordinates will be relative to that of the object.
			
			-	"-s"  		sun_centered: 	particle positions will be converted to Sun centred coordinates before calculating Kepler orbits
			-	"-v"		velocities:		output in normal mode will contain velocity coordinates
			-	"-n"		nodes:			compute particle nodes in the orbital plane of OBJECT_ID (requires --object to be set)
			-	"-p"		pc2-format:		the output file will have .pc2 file format.
			
	If "--object" is set it is necessary to specify the directory of the SPK kernel ,corresponding to OBJECT_ID, in "kernels_spk.txt".
	Also if "-s" is set it is necessary to specify the directory of a Sun SPK kernel (e.g. de430.bsp) in "kernels_spk.txt".
		
			
-	output format (all values are stored as floats)

		-	normal mode: 	dimension  8*particle_count 
		
						1		2		3			4					5				6		7		 8
						X		Y		Z	 particle_number  MultiplicationFactor	  mass	  beta	creatio_time
						
						
		-	"-v: 			dimension 11*particle_count 
		
						1	2	3	4	5	6		    7				   8				9	   10		 11
						X	Y	Z	Vx	Vy	Vz	 particle_number  MultiplicationFactor	  mass	  beta	creatio_time
						
						
		-	"-n: 			dimension 12*particle_count 
		
						1	2	3	4	5	6		 7	  	  	   		8				  9				   10	   11		 12
						X	Y	Z	Vx	Vy	Vz	 node_time		particle_number  MultiplicationFactor	  mass	  beta	creation_time
						
						
		-	"-p": 	.pc2 (PointCache 2) format
				
					for details see http://www.footools.com/plugins/docs/PointCache2_README.html#pc2format
					
					
					
-	example input (values can be double precision):

		panalyser --input D:\DATABASE\1000095 --output \storeOutput\leonids1999 --time -3.8447e6 --tmin -5.3647e9 --mmax 1e-6 --object 3 5e7 -sn
		
		Using the above input arguments PAnalyser will do the following:
		- Search for work units in D:\DATABASE\1000095 	(exemplary directory of TempelTuttle results).
		- It looks for all particles created since TMIN (-5.3647e9 ~ 1 JAN 1830) and before TIME (-3.8447e6 ~ 18 NOV 1999) and have a mass of 1e-6 kg and less.
		- It computes the nodes of all those particles, which are within a distance of 5e7 km of Earth at TIME, in the Earth's orbital plane using Sun-centred Kepler orbits.
		- It saves these nodes in a binary file called "leonids1999" in a folder ".\storeOutput" using the format "-n" described above.
						
						
-	executable returns:

    >   0:	the app runs successfully to its end  
    >   1: 	no success  (see console)
	
	
Regards, Max


##Changelog

__0.01__	-29.06.2014-:
-	initial alpha release

