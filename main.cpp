#include <iostream>
#include <climits>
#include "tsp.h"
#include "usage.h"
#include "MyThread.h"
#define CPMS CLOCKS_PER_SEC * 1000

#define NUM_THREADS 16

#define TOTAL_CITY_COUNT (NUM_THREADS * 256)

#define BLOCK_SIDE_LENGTH 100

int main(int argc, char **argv)
{
	// Check argument count
	if (argc < 2)
		usage();

	// Flags to check if mode and algorithm are set correctly
	int mode_flag = 0;
	int algo_flag = 0;

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int my_world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

	char *mode = argv[1];
	char *algo = argv[2];

	// Print total number of cities and number of threads
	if (my_world_rank == 0)
	{
		if (strcmp(mode, "serial") != 0)
		{
			cout << "Total city count: " << TOTAL_CITY_COUNT << endl;
			cout << "Number of threads: " << NUM_THREADS << endl;
		}
	}

	clock_t start;
	clock_t end;

	if (strcmp(mode, "MPI") == 0)
	{
		mode_flag = 1;
		TSP tsp(TOTAL_CITY_COUNT, NUM_THREADS, BLOCK_SIDE_LENGTH);
		if (strcmp(algo, "DP") == 0)
		{
			algo_flag = 1;
			start = clock();
			tsp.MPI_solver(TSP::Algorithm::DP, my_world_rank, world_size);
			end = clock();
			if (my_world_rank == 0)
				cout << "MPI DP Time: " << ((float)(end - start)) / CPMS << " ms\n\n";
		}
		else if (strcmp(algo, "EMST") == 0)
		{
			algo_flag = 1;
			start = clock();
			tsp.MPI_solver(TSP::Algorithm::EMST, my_world_rank, world_size);
			end = clock();
			if (my_world_rank == 0)
				cout << "MPI EMST Time: " << ((float)(end - start)) / CPMS << " ms\n\n";
		}
		if (my_world_rank == 0)
		{
			tsp.fix_inversions();
			tsp.printResults();
		}
	}

	if (my_world_rank == 0)
	{
		if (strcmp(mode, "serial") == 0)
		{
			mode_flag = 1;
			if (argc < 4)
				usage();
			// Read file name
			string f, o;
			f = o = argv[3];
			o.append(".tour");

			TSP tsp(f, o);
			tsp.readCities();
			if (strcmp(algo, "DP") == 0)
			{
				algo_flag = 1;
				start = clock();
				tsp.calculate_distances();
				tsp.serial_DP();
				end = clock();
				cout << "Serial DP Time: " << ((float)(end - start)) / CPMS << " ms\n\n";
			}
			else if (strcmp(algo, "EMST") == 0)
			{
				algo_flag = 1;
				start = clock();
				tsp.calculate_distances();
				tsp.serial_EMST();
				end = clock();
				cout << "Serial EMST Time: " << ((float)(end - start)) / CPMS << " ms\n\n";
			}
			tsp.printResults();
		}
		else if (strcmp(mode, "parallel") == 0)
		{
			mode_flag = 1;
			TSP tsp(TOTAL_CITY_COUNT, NUM_THREADS, BLOCK_SIDE_LENGTH);
			if (strcmp(algo, "DP") == 0)
			{
				algo_flag = 1;
				start = clock();
				tsp.parallel_solver(NUM_THREADS, TSP::Algorithm::DP);
				end = clock();
				cout << "Parallel DP Time: " << ((float)(end - start)) / CPMS << " ms\n\n";
			}
			else if (strcmp(algo, "EMST") == 0)
			{
				algo_flag = 1;
				start = clock();
				tsp.parallel_solver(NUM_THREADS, TSP::Algorithm::EMST);
				end = clock();
				cout << "Parallel EMST Time: " << ((float)(end - start)) / CPMS << " ms\n\n";
			}
			// tsp.dump_all_cities();
			// tsp.fix_inversions();
			tsp.printResults();
		}

		if (mode_flag == 0 || algo_flag == 0)
			usage();
	}

	return 0;
}
