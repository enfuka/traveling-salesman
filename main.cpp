#include <iostream>
#include <climits>
#include "tsp.h"
#include "usage.h"
#include "MyThread.h"
#define CPMS CLOCKS_PER_SEC * 1000

#define NUM_THREADS 4

#define TOTAL_CITY_COUNT (NUM_THREADS * 4096)

#define BLOCK_SIDE_LENGTH 100

int main(int argc, char **argv)
{
	// Check argument count
	if (argc < 2)
		usage();

	// Flags to check if mode and algorithm are set correctly
	int mode_flag = 0;
	int algo_flag = 0;

	///// MPI SETUP /////
	int ndims = 2;
	int dims[ndims], my_coord[ndims];
	int wrap_around[ndims], reorder;
	int my_cart_rank;
	int nrows, ncols;
	int cnt, rlen;
	int row_val, row_sum = 0, row_root = -999;
	int col_val, col_sum = 0, col_root = -999;

	char name[200], nameout[200], rname[100];

	MPI_Comm comm2D;	 // Cartesian communicator
	MPI_Comm comm1D_row; // Communicator for row split
	MPI_Comm comm1D_col; // Communicator for col split

	// For split
	int my_row_rank;
	int my_col_rank;
	int free_coords[ndims];
	int my_row_coord[ndims];
	int my_col_coord[ndims];

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int my_world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

	// /* get the name of communicator */
	// nameout[0] = 0;
	// MPI_Comm_get_name(MPI_COMM_WORLD, nameout, &rlen);
	// if (my_world_rank == 0)
	// 	printf("Name of comm world is: %s\n", nameout);
	// fflush(stdout);

	// Create a Cartesian communicator
	nrows = ncols = sqrt(world_size);
	dims[0] = nrows;					 // number of rows
	dims[1] = ncols;					 // number of columns
	wrap_around[0] = wrap_around[1] = 1; // periodic shift is .true.
	reorder = 1;						 // reorder is .true.
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm2D);

	/* set the name of cartesian communicator */
	strcpy(name, "comm2D");
	MPI_Comm_set_name(comm2D, name);
	nameout[0] = 0;

	// /* get the name of communicator */
	// MPI_Comm_get_name(comm2D, nameout, &rlen);
	// if (my_world_rank == 0)
	// 	printf("Name of comm world is: %s\n", nameout);
	// fflush(stdout);

	/* find my coordinates in the cartesian communicator group */
	MPI_Cart_coords(comm2D, my_world_rank, ndims, my_coord);

	/* use my cartesian coordinates to find my rank in cartesian group*/
	MPI_Cart_rank(comm2D, my_coord, &my_cart_rank);

	// Split comm2D into a 1D row-based communicator: comm1D_row
	free_coords[0] = 0; /* rows */
	free_coords[1] = 1; /* cols */
	MPI_Cart_sub(comm2D, free_coords, &comm1D_row);

	/* get my_row_rank in my comm1D_row group */
	MPI_Comm_rank(comm1D_row, &my_row_rank);

	/* use my_row_rank to find my coordinates in my comm1D_row group */
	MPI_Cart_coords(comm1D_row, my_row_rank, 1, my_row_coord);

	/* set the name of cartesian communicator */
	if (my_row_rank == 0)
	{
		row_root = my_world_rank;
	}

	MPI_Bcast(&row_root, 1, MPI_INT, 0, comm1D_row);

	// Split comm2D into a 1D col-based communicator: comm1D_col
	free_coords[0] = 1; /* rows */
	free_coords[1] = 0; /* cols */
	MPI_Cart_sub(comm2D, free_coords, &comm1D_col);

	/* get my_col_rank in my comm1D_col group */
	MPI_Comm_rank(comm1D_col, &my_col_rank);

	/* use my_col_rank to find my coordinates in my comm1D_col group */
	MPI_Cart_coords(comm1D_col, my_col_rank, 1, my_col_coord);

	/* set the name of cartesian communicator */
	if (my_col_rank == 0)
	{
		col_root = my_world_rank;
	}
	MPI_Bcast(&col_root, 1, MPI_INT, 0, comm1D_col);

	/* test row com: row_sum processor ranks across rows */
	if (my_coord[1] == 0)
		row_val = my_coord[0];
	else
		row_val = -1;
	MPI_Bcast(&row_val, 1, MPI_INT, 0, comm1D_row);
	MPI_Reduce(&my_world_rank, &row_sum, 1, MPI_INT, MPI_SUM, 0, comm1D_row);
	MPI_Bcast(&row_sum, 1, MPI_INT, 0, comm1D_row);

	printf("WRank: %d, RowRank: %d, ColRank:%d, Coord: (%d, %d), RowRoot: %d, ColRoot: %d, RowSum: %d\n",
		   my_world_rank, my_row_rank, my_col_rank, my_coord[0], my_coord[1], row_root, col_root, row_sum);

	MPI_Barrier(MPI_COMM_WORLD);

	///// MPI SETUP /////

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
			tsp.dump_all_cities();
			// tsp.fix_inversions();
			tsp.printResults();
		}

		if (mode_flag == 0 || algo_flag == 0)
			usage();
	}

	// Finalize the MPI environment
	MPI_Comm_free(&comm2D);
	MPI_Comm_free(&comm1D_row);
	MPI_Comm_free(&comm1D_col);
	MPI_Finalize();

	return 0;
}
