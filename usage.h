#include <iostream>

#ifndef USAGE_H
#define USAGE_H

// Print correct usage to command line
static void usage()
{
	printf("Usage: mpirun -np [number of processors(1 if not MPI)] tsp <options> [inputfile.txt]\n");
	printf("OPTIONS: \n");
	printf("1. 'serial DP' \n");
	printf("2. 'serial EMST' \n");
	printf("3. 'parallel DP' \n");
	printf("4. 'parallel EMST' \n");
	printf("3. 'MPI DP' \n");
	printf("4. 'MPI EMST' \n");

	printf("FILENAME IS ONLY NEEDED FOR SERIAL VERSIONS \n");
	printf("where inputfile is in the format:\n");
	printf("   0 x-coord y-coord\n");
	printf("   1 x-coord y-coord\n");
	printf("   .\n   .\n   .\n");
	printf("   n-1 x-coord y-coord\n");
	printf("and n is the number of cities\n");
	printf("\n");

	exit(1);
}

#endif
