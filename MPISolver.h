#include "tsp.h"

using namespace std;

// MPI implementation
void TSP::MPI_solver(Algorithm algorithm, int my_world_rank, int world_size)
{
    ///// MPI SETUP /////
    int ndims = 2;
    int dims[ndims], my_coord[ndims];
    int wrap_around[ndims], reorder;
    int my_cart_rank;
    int nrows, ncols;
    int local_val, row_sum = 0;
    int col_val, col_sum = 0;

    MPI_Comm comm2D;     // Cartesian communicator
    MPI_Comm comm1D_row; // Communicator for row split
    MPI_Comm comm1D_col; // Communicator for col split

    // For split
    int my_row_rank;
    int my_col_rank;
    int free_coords[ndims];
    int my_row_coord[ndims];
    int my_col_coord[ndims];

    // Create a Cartesian communicator
    nrows = ncols = sqrt(world_size);
    dims[0] = nrows;                     // number of rows
    dims[1] = ncols;                     // number of columns
    wrap_around[0] = wrap_around[1] = 1; // periodic shift is .true.
    reorder = 1;                         // reorder is .true.
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm2D);

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

    // Split comm2D into a 1D col-based communicator: comm1D_col
    free_coords[0] = 1; /* rows */
    free_coords[1] = 0; /* cols */
    MPI_Cart_sub(comm2D, free_coords, &comm1D_col);

    /* get my_col_rank in my comm1D_col group */
    MPI_Comm_rank(comm1D_col, &my_col_rank);

    /* use my_col_rank to find my coordinates in my comm1D_col group */
    MPI_Cart_coords(comm1D_col, my_col_rank, 1, my_col_coord);
    ///// MPI SETUP /////

    //~~~~~~~~ Parallel Work (begin) ~~~~~~~~//

    // Regular array to store the results
    City *full_raw_path;
    int total_path_length;

    // Create a new array to store the path
    full_raw_path = new City[world_size * cities_per_block];

    // Create a new vector of cities for the thread
    vector<City> cities;
    std::vector<City> tour_with_coordinates;
    int **emst_distances = new int *[cities_per_block];
    std::vector<std::vector<double>> dp_distances = std::vector<std::vector<double>>(cities_per_block, std::vector<double>(cities_per_block));

    thread_populate_block(my_world_rank, cities);

    if (DEBUG)
        print_thread_cities(my_world_rank, cities);

    if (algorithm == TSP::Algorithm::DP)
    {
        // Calculate the distance between each pair of cities
        for (int i = 0; i < cities_per_block; i++)
        {
            for (int j = 0; j < cities_per_block; j++)
            {
                dp_distances[i][j] = calculate_distance(cities[i], cities[j]);
            }
        }

        if (DEBUG)
            print_thread_distances_DP(my_world_rank, dp_distances);

        // Find the path and length for the thread
        std::list<int> tour;
        DPSolver dp(0, dp_distances, 0);
        tour = dp.getTour();
        pathLength = dp.getTourCost();

        // Map tour to city objects

        for (std::list<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
        {
            tour_with_coordinates.push_back(cities[*it]);
        }
    }
    else
    {
        // Calculate the distance between each pair of cities
        for (int i = 0; i < cities_per_block; i++)
        {
            emst_distances[i] = new int[cities_per_block];
            for (int j = 0; j < cities_per_block; j++)
            {
                emst_distances[i][j] = calculate_distance(cities[i], cities[j]);
            }
        }

        if (DEBUG)
            print_thread_distances_EMST(my_world_rank, emst_distances);

        // Find the path and length for the thread
        std::vector<int> tour;
        EMSTSolver emst(cities_per_block, emst_distances, 0);
        tour = emst.find_path();
        pathLength = emst.get_cost();

        // Map tour to city objects, skip the last city
        for (std::vector<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
        {
            tour_with_coordinates.push_back(cities[*it]);
        }
    }

    // Convert the tour to a simple array of city ids
    City *block_path = new City[cities_per_block];
    for (int i = 0; i < cities_per_block; i++)
    {
        block_path[i] = tour_with_coordinates[i];
    }

    // MOCK STICHING (linear)
    local_val = my_world_rank;
    int recv_val;

    // Stich together block TSP solutions row-wise, each row separately in parallel
    for (int cur_sticher = 1; cur_sticher < ncols; cur_sticher++)
    {
        // SEND
        if (my_row_rank == cur_sticher - 1) // if I am the sender (one before the processor that stiches)
        {
            int destination = my_row_rank + 1;
            if (my_row_rank == 0)
                MPI_Send(&local_val, 1, MPI_INT, destination, 0, comm1D_row);
            else
                MPI_Send(&row_sum, 1, MPI_INT, destination, 0, comm1D_row);
        }
        // RECEIVE
        if (my_row_rank == cur_sticher) // if I am the sticher
        {
            int source = my_row_rank - 1;
            MPI_Recv(&recv_val, 1, MPI_INT, source, 0, comm1D_row, MPI_STATUS_IGNORE);
            // STICH
            row_sum = local_val + recv_val;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // sum col values across cols one by one
    if (my_row_rank == ncols - 1) // only the last processors in each row will stich
    {
        for (int cur_sticher = 1; cur_sticher < nrows; cur_sticher++)
        {
            // SEND
            if (my_col_rank == cur_sticher - 1) // if I am the sender (one before the processor that stiches)
            {
                int destination = my_col_rank + 1;
                if (my_col_rank == 0)
                    MPI_Send(&row_sum, 1, MPI_INT, destination, 0, comm1D_col);
                else
                    MPI_Send(&col_sum, 1, MPI_INT, destination, 0, comm1D_col);
            }
            // RECEIVE
            if (my_col_rank == cur_sticher) // if I am the sticher
            {
                int source = my_col_rank - 1;
                MPI_Recv(&recv_val, 1, MPI_INT, source, 0, comm1D_col, MPI_STATUS_IGNORE);
                // STICH
                col_sum = row_sum + recv_val;
            }
        }
    }

    if (my_row_rank == ncols - 1)
        printf("P%d: RowSum: %d\n", my_world_rank, row_sum);

    if (my_world_rank == world_size - 1)
        printf("P%d: ColSum: %d\n", my_world_rank, col_sum);

    printf("WRank: %d, RowRank: %d, ColRank:%d, Coord: (%d, %d), RowSum: %d\n",
           my_world_rank, my_row_rank, my_col_rank, my_coord[0], my_coord[1], row_sum);

    //~~~~~~~~ Parallel Work (end) ~~~~~~~~//

    // Gather all the results
    MPI_Gather(block_path, cities_per_block * sizeof(City), MPI_BYTE, full_raw_path, cities_per_block * sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Gather all the path lengths
    MPI_Reduce(&pathLength, &total_path_length, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Stich together the paths from each thread
    if (my_world_rank == 0)
    {
        final_tour.clear();
        // At the end of a row, continue adding the next row in the reverse order
        int grid_side = sqrt(world_size);

        for (int row = 0; row < grid_side; row++)
        {
            if (row % 2 != 0)
            {
                for (int column = grid_side - 1; column >= 0; column--)
                {
                    int process_id = row * grid_side + column;
                    // Insert the block's path into the total path
                    for (int i = 0; i < cities_per_block; i++)
                    {
                        final_tour.push_back(full_raw_path[process_id * cities_per_block + i]);
                    }

                    // Add the distance between the last city of the previous block and the first city of the current block
                    if (column != grid_side - 1) // unless last column
                    {
                        // Add the distance between next block's last city and current block's first city
                        total_path_length += calculate_distance(full_raw_path[((process_id + 1) * cities_per_block) + (cities_per_block - 1)], full_raw_path[process_id * cities_per_block]);
                    }
                    else // if last columns
                    {
                        total_path_length += calculate_distance(full_raw_path[((process_id - grid_side) * cities_per_block) + (cities_per_block - 1)], full_raw_path[process_id * cities_per_block]);
                    }
                }
            }
            else
            {
                for (int column = 0; column < grid_side; column++)
                {
                    int process_id = row * grid_side + column;
                    // Insert the block's path into the total path
                    for (int i = 0; i < cities_per_block; i++)
                    {
                        final_tour.push_back(full_raw_path[process_id * cities_per_block + i]);
                    }
                    // Add the distance between the last city of the previous block and the first city of the current block
                    if (column != 0) // unless first column
                    {
                        total_path_length += calculate_distance(full_raw_path[((process_id - 1) * cities_per_block) + (cities_per_block - 1)], full_raw_path[process_id * cities_per_block]);
                    }
                    else if (row != 0) // unless first row but first column
                    {
                        total_path_length += calculate_distance(full_raw_path[((process_id - grid_side) * cities_per_block) + (cities_per_block - 1)], full_raw_path[process_id * cities_per_block]);
                    }
                }
            }
        }

        // Add the distance between the last city and the first city
        City last_city = final_tour.back();
        total_path_length += calculate_distance(last_city, full_raw_path[0]);

        // Add first city to end of tour
        final_tour.push_back(full_raw_path[0]);

        // Store the final tour
        pathLength = total_path_length;
    }

    // Finalize the MPI environment
    MPI_Comm_free(&comm2D);
    MPI_Comm_free(&comm1D_row);
    MPI_Comm_free(&comm1D_col);
    MPI_Finalize();

    // // Free memory
    // if (Algorithm::DP == algorithm)
    // {
    //     for (int i = 0; i < cities_per_block; i++)
    //     {
    //         delete[] &dp_distances[i];
    //     }
    // }
    // else
    // {
    //     for (int i = 0; i < cities_per_block; i++)
    //     {
    //         delete[] emst_distances[i];
    //     }
    // }
}