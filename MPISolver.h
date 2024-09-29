#include "tsp.h"

using namespace std;

// MPI implementation
void TSP::MPI_solver(Algorithm algorithm, int world_rank, int world_size)
{

    // Regular array to store the results
    City *path;
    int total_path_length;

    // Create a new array to store the path
    path = new City[world_size * cities_per_block];

    // Create a new vector of cities for the thread
    vector<City> cities;
    std::vector<City> tour_with_coordinates;
    int **emst_distances = new int *[cities_per_block];
    std::vector<std::vector<double>> dp_distances = std::vector<std::vector<double>>(cities_per_block, std::vector<double>(cities_per_block));

    thread_populate_block(world_rank, cities);

    if (DEBUG)
        print_thread_cities(world_rank, cities);

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
            print_thread_distances_DP(world_rank, dp_distances);

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
            print_thread_distances_EMST(world_rank, emst_distances);

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
    City *tour_simple = new City[cities_per_block];
    for (int i = 0; i < cities_per_block; i++)
    {
        tour_simple[i] = tour_with_coordinates[i];
    }

    // Sync all the threads
    // MPI_Barrier(MPI_COMM_WORLD);

    // Gather all the results
    MPI_Gather(tour_simple, cities_per_block * sizeof(City), MPI_BYTE, path, cities_per_block * sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Gather all the path lengths
    MPI_Reduce(&pathLength, &total_path_length, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Stich together the paths from each thread
    if (world_rank == 0)
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
                        final_tour.push_back(path[process_id * cities_per_block + i]);
                    }

                    // Add the distance between the last city of the previous block and the first city of the current block
                    if (column != grid_side - 1) // unless last column
                    {
                        // Add the distance between next block's last city and current block's first city
                        total_path_length += calculate_distance(path[((process_id + 1) * cities_per_block) + (cities_per_block - 1)], path[process_id * cities_per_block]);
                    }
                    else // if last columns
                    {
                        total_path_length += calculate_distance(path[((process_id - grid_side) * cities_per_block) + (cities_per_block - 1)], path[process_id * cities_per_block]);
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
                        final_tour.push_back(path[process_id * cities_per_block + i]);
                    }
                    // Add the distance between the last city of the previous block and the first city of the current block
                    if (column != 0) // unless first column
                    {
                        total_path_length += calculate_distance(path[((process_id - 1) * cities_per_block) + (cities_per_block - 1)], path[process_id * cities_per_block]);
                    }
                    else if (row != 0) // unless first row but first column
                    {
                        total_path_length += calculate_distance(path[((process_id - grid_side) * cities_per_block) + (cities_per_block - 1)], path[process_id * cities_per_block]);
                    }
                }
            }
        }

        // Add the distance between the last city and the first city
        City last_city = final_tour.back();
        total_path_length += calculate_distance(last_city, path[0]);

        // Add first city to end of tour
        final_tour.push_back(path[0]);

        // Store the final tour
        pathLength = total_path_length;
    }

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