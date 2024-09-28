#include "tsp.h"
#include "MyThread.h"

using namespace std;

// Constructor for serial TSP
TSP::TSP(string in, string out)

{
	/////////////////////////////////////////////////////
	// Constructor
	/////////////////////////////////////////////////////
	inFname = in;
	outFname = out;

	// set n to number of lines read from input file
	getNodeCount();

	// Allocate memory
	graph = new int *[n];
	for (int i = 0; i < n; i++)
	{
		graph[i] = new int[n];
		for (int j = 0; j < n; j++)
			graph[i][j] = 0;
	}
}

// Constructor for parallel TSP
TSP::TSP(int total_city_count, int thread_count, int block_side_length)

{
	// calculate number of cities per block if running parallel
	cities_per_block = total_city_count / thread_count;

	// set block side length
	this->block_side_length = block_side_length;

	thread_results = std::vector<std::pair<std::vector<TSP::City>, int>>(thread_count);
};

TSP::~TSP()
{
	/////////////////////////////////////////////////////
	// Destructor
	/////////////////////////////////////////////////////
	// Free memory
	if (graph != NULL && n > 0)
	{
		for (int i = 0; i < n; i++)
		{
			delete[] graph[i];
		}
		delete[] graph;
	}
}

void TSP::getNodeCount()
{
	int count = 0;
	ifstream inStream;
	inStream.open(inFname.c_str(), ios::in);

	if (!inStream)
	{
		cerr << "Can't open input file " << inFname << endl;
		exit(1);
	}
	std::string unused;
	while (std::getline(inStream, unused))
		++count;
	n = count;
	inStream.close();
};

void TSP::readCities()
{
	/////////////////////////////////////////////////////
	ifstream inStream;
	inStream.open(inFname.c_str(), ios::in);
	if (!inStream)
	{
		cerr << "Can't open input file " << inFname << endl;
		exit(1);
	}
	int c, x, y;
	int i = 0;
	while (!inStream.eof())
	{
		inStream >> c >> x >> y;
		// Push back new city to vector
		struct City c = {i, x, y};
		cities.push_back(c);
		i++;
	}
	inStream.close();
};

int TSP::get_distance(struct TSP::City c1, struct TSP::City c2)
{
	/////////////////////////////////////////////////////
	// Calculate distance between c1 and c2
	/////////////////////////////////////////////////////
	int dx = pow((float)(c1.x - c2.x), 2);
	int dy = pow((float)(c1.y - c2.y), 2);
	return (floor((float)(sqrt(dx + dy)) + 0.5));
};

void TSP::calculate_distances()
{
	/////////////////////////////////////////////////////
	// Calculate distance between each pair of cities
	/////////////////////////////////////////////////////
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			graph[i][j] = graph[j][i] = get_distance(cities[i], cities[j]);
		}
	}
}

void TSP::serial_DP()
{
	/////////////////////////////////////////////////////
	// Dynamic programming to solve TSP
	/////////////////////////////////////////////////////
	// Create a copy of the graph
	std::vector<std::vector<double>> distance = std::vector<std::vector<double>>(n, std::vector<double>(n));
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
			distance[i][j] = distance[j][i] = graph[i][j];
	}

	DPSolver dp(0, distance, 0);
	cout << "Final path: ";
	std::list<int> tour = dp.getTour();
	pathLength = dp.getTourCost();

	final_tour.clear();
	for (std::list<int>::iterator it = tour.begin(); it != tour.end(); ++it)
		final_tour.push_back(*it);
}

void TSP::serial_EMST()
{
	/////////////////////////////////////////////////////
	// Eucledian minimum spanning tree using Prim's algorithm to solve TSP
	/////////////////////////////////////////////////////

	// Create a copy of the graph

	int **distance = new int *[n];
	for (int i = 0; i < n; i++)
	{
		distance[i] = new int[n];
		for (int j = 0; j < n; j++)
			distance[i][j] = graph[i][j];
	}

	// creating an object of the class
	EMSTSolver emst(n, distance, 0);

	// finding the path
	std::vector<int> tour = emst.find_path();

	// finding the cost of the path
	pathLength = emst.get_cost();

	final_tour.clear();
	for (std::vector<int>::iterator it = tour.begin(); it != tour.end(); ++it)
		final_tour.push_back(*it);
}

void TSP::parallel_solver(int num_threads, Algorithm algorithm)
{
	// Create array of thread objects
	MyThread threads[num_threads];

	// Create a thread for each block
	for (int i = 0; i < num_threads; i++)
	{
		threads[i].mytsp = this;
		threads[i].algorithm = algorithm;
		threads[i].my_id = i;
		threads[i].start();
		// sleep for a bit so we can see the output
		sleep(1);
	}

	// Wait for all the threads
	for (int i = 0; i < num_threads; i++)
	{
		threads[i].join();
	}

	print_thread_paths(num_threads);

	// Stich together the paths from each thread
	// At the end of a row, continue adding the next row in the reverse order
	int grid_side = sqrt(num_threads);
	std::vector<TSP::City> tour;
	int pathLength = 0;
	for (int row = 0; row < grid_side; row++)
	{
		if (row % 2 != 0)
		{
			for (int column = grid_side - 1; column >= 0; column--)
			{
				int thread_id = row * grid_side + column;
				// Insert the block's path into the total path
				tour.insert(tour.end(), thread_results[thread_id].first.begin(), thread_results[thread_id].first.end());
				// Add the block's path length to the total path length
				pathLength += thread_results[thread_id].second;
				// Add the distance between the last city of the previous block and the first city of the current block
				if (column != grid_side - 1) // unless last column
				{
					pathLength += get_distance(thread_results[thread_id + 1].first.back(), thread_results[thread_id].first.front());
				}
				else // if last column
				{
					pathLength += get_distance(thread_results[thread_id - grid_side].first.back(), thread_results[thread_id].first.front());
				}
			}
		}
		else
		{
			for (int column = 0; column < grid_side; column++)
			{
				int thread_id = row * grid_side + column;
				// Insert the block's path into the total path
				tour.insert(tour.end(), thread_results[thread_id].first.begin(), thread_results[thread_id].first.end());
				// Add the block's path length to the total path length
				pathLength += thread_results[thread_id].second;
				// Add the distance between the last city of the previous block and the first city of the current block
				if (column != 0) // unless first column
				{
					pathLength += get_distance(thread_results[thread_id - 1].first.back(), thread_results[thread_id].first.front());
				}
				else if (row != 0) // unless first row but first column
				{
					pathLength += get_distance(thread_results[thread_id - grid_side].first.back(), thread_results[thread_id].first.front());
				}
			}
		}
	}

	// Add the distance between the last city and the first city
	City last_city = (grid_side % 2 != 0) ? thread_results[num_threads - 1].first.back() : thread_results[num_threads - grid_side].first.back();

	pathLength += get_distance(last_city, thread_results[0].first.front());

	// Add first city to end of tour
	tour.push_back(thread_results[0].first.front());

	// Print the final path and length
	cout << "Final path: ";
	for (int i = 0; i < tour.size(); i++)
	{
		cout << tour[i].id << " ";
	}
	cout << endl;
	cout << "Final length: " << pathLength << endl;
}

// Function that each thread will run for the EMST parallel implementation
void TSP::openTSP_EMST(int thread_id)
{
	// Create a new vector of cities for the thread
	vector<City> cities;

	thread_populate_block(thread_id, cities);

	print_thread_cities(thread_id, cities);

	// Calculate the distance between each pair of cities
	int **emst_distances = new int *[cities_per_block];
	for (int i = 0; i < cities_per_block; i++)
	{
		emst_distances[i] = new int[cities_per_block];
		for (int j = 0; j < cities_per_block; j++)
		{
			emst_distances[i][j] = get_distance(cities[i], cities[j]);
		}
	}

	print_thread_distances_EMST(thread_id, emst_distances);

	// Find the path and length for the thread
	std::vector<int> tour;
	int pathLength;
	EMSTSolver emst(cities_per_block, emst_distances, 0);
	tour = emst.find_path();
	pathLength = emst.get_cost();

	// Map tour to city objects, skip the last city
	std::vector<TSP::City> tour_with_coordinates;
	for (std::vector<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
	{
		tour_with_coordinates.push_back(cities[*it]);
	}

	// Record the path and length in the right index in thread_paths
	thread_results[thread_id] = make_pair(tour_with_coordinates, pathLength);
}

void TSP::openTSP_DP(int thread_id)
{
	// Create a new vector of cities for the thread
	vector<City> cities;

	thread_populate_block(thread_id, cities);

	print_thread_cities(thread_id, cities);

	// Calculate the distance between each pair of cities
	std::vector<std::vector<double>> dp_distances = std::vector<std::vector<double>>(cities_per_block, std::vector<double>(cities_per_block));
	for (int i = 0; i < cities_per_block; i++)
	{
		for (int j = 0; j < cities_per_block; j++)
		{
			dp_distances[i][j] = get_distance(cities[i], cities[j]);
		}
	}

	print_thread_distances_DP(thread_id, dp_distances);

	// Find the path and length for the thread
	std::list<int> tour;
	int pathLength;
	DPSolver dp(0, dp_distances, 0);
	tour = dp.getTour();
	pathLength = dp.getTourCost();

	// Print tour
	cout << "Thread " << thread_id << " tour: ";
	for (std::list<int>::iterator it = tour.begin(); it != tour.end(); ++it)
	{
		cout << *it << " ";
	}

	// Map tour to city objects
	std::vector<TSP::City> tour_with_coordinates;
	for (std::list<int>::iterator it = tour.begin(); it != prev(tour.end()); ++it)
	{
		tour_with_coordinates.push_back(cities[*it]);
	}

	// Record the path and length in the right index in thread_paths
	thread_results[thread_id] = make_pair(tour_with_coordinates, pathLength);
}

void TSP::make_shorter()
{
	// Modify final_tour & pathLength
	twoOpt(graph, final_tour, pathLength, n);
}

void TSP::thread_populate_block(int thread_id, std::vector<TSP::City> &cities)
{
	// seed random number generator
	srand(thread_id);

	for (int i = 0; i < cities_per_block; i++)
	{
		// Create a random city in the cell, have relative id's
		int x = rand() % block_side_length + thread_id * block_side_length;
		int y = rand() % block_side_length + thread_id * block_side_length;
		City city = {thread_id * cities_per_block + i, x, y};
		cities.push_back(city);
	}
}

void TSP::dump_all_cities()
{
	// Extract all city coordinates from thread_results and write to file
	ofstream outputStream;
	outputStream.open("all_cities.txt", ios::out);
	int index = 0;
	for (int i = 0; i < thread_results.size(); i++)
	{
		for (int j = 0; j < thread_results[i].first.size(); j++)
		{
			outputStream << index << " " << thread_results[i].first[j].x << " " << thread_results[i].first[j].y << endl;
			index++;
		}
	}
	outputStream.close();
}

void TSP::MPI_solver()
{
	// MPI implementation

	// Regular array to store the results
	City *path;
	int total_path_length;

	// Initialize MPI
	MPI_Init(NULL, NULL);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	// Print off a hello world message
	printf("Hello world from processor %s, rank %d out of %d processors\n",
		   processor_name, world_rank, world_size);

	// // allocate memory for the results
	// process_results = std::vector<std::pair<std::vector<TSP::City>, int>>(world_size);

	path = new City[world_size * cities_per_block];

	// Create a new vector of cities for the thread
	vector<City> cities;

	thread_populate_block(world_rank, cities);

	print_thread_cities(world_rank, cities);

	// Calculate the distance between each pair of cities
	int **emst_distances = new int *[cities_per_block];
	for (int i = 0; i < cities_per_block; i++)
	{
		emst_distances[i] = new int[cities_per_block];
		for (int j = 0; j < cities_per_block; j++)
		{
			emst_distances[i][j] = get_distance(cities[i], cities[j]);
		}
	}

	print_thread_distances_EMST(world_rank, emst_distances);

	// Find the path and length for the thread
	std::vector<int> tour;
	int pathLength;
	EMSTSolver emst(cities_per_block, emst_distances, 0);
	tour = emst.find_path();
	pathLength = emst.get_cost();

	// Map tour to city objects, skip the last city
	std::vector<TSP::City> tour_with_coordinates;
	for (std::vector<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
	{
		tour_with_coordinates.push_back(cities[*it]);
	}

	// Convert the tour to a simple array of city ids
	City *tour_simple = new City[cities_per_block];
	for (int i = 0; i < cities_per_block; i++)
	{
		tour_simple[i] = tour_with_coordinates[i];
	}

	// Sync all the threads
	MPI_Barrier(MPI_COMM_WORLD);

	// Gather all the results
	MPI_Gather(tour_simple, cities_per_block * sizeof(City), MPI_BYTE, path, cities_per_block * sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);

	// Gather all the path lengths
	MPI_Reduce(&pathLength, &total_path_length, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	// Print the results
	if (world_rank == 0)
	{
		// Add the distance between the last city and the first city
		City last_city = path[world_size * cities_per_block - 1];
		total_path_length += get_distance(last_city, path[0]);
		// Print the final path
		cout << "Final path: " << endl;
		for (int i = 0; i < world_size * cities_per_block; i++)
		{
			cout << path[i].id << " " << path[i].x << " " << path[i].y << endl;
		}
		cout << "Final length: " << total_path_length << endl;
	}

	// Finalize the MPI environment
	MPI_Finalize();

	// Free memory
	free(path);
}

//================================ PRINT FUNCTIONS ================================//

void TSP::printResult()
{
	ofstream outputStream;
	outputStream.open(outFname.c_str(), ios::out);
	outputStream << pathLength << endl;
	for (vector<int>::iterator it = final_tour.begin(); it != final_tour.end(); ++it)
	{
		// for (vector<int>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
		outputStream << *it << endl;
	}
	// outputStream << *(circuit.end()-1);
	outputStream.close();
};

void TSP::printPath()
{
	cout << endl;
	for (vector<int>::iterator it = final_tour.begin(); it != final_tour.end() - 1; ++it)
	{
		cout << *it << " to " << *(it + 1) << " ";
		cout << graph[*it][*(it + 1)] << endl;
	}

	cout << "\nLength: " << pathLength << endl
		 << endl;
};

void TSP::printCities()
{
	cout << endl;

	for (vector<City>::iterator it = cities.begin(); it != cities.end(); ++it)
		cout << it->id << ":  " << it->x << " " << it->y << endl;
}

void TSP::printDistanceGraph()
{
	cout << endl;
	for (int i = 0; i < n; i++)
	{
		cout << i << ": ";
		for (int j = 0; j < n; j++)
			cout << graph[i][j] << " ";
		cout << endl;
	}
};

void TSP::print_thread_cities(int thread_id, std::vector<TSP::City> &cities)
{
	// Print cities
	cout << endl;
	cout << "Thread " << thread_id << " cities:\n";
	for (int i = 0; i < cities_per_block; i++)
	{
		cout << cities[i].id << " (" << cities[i].x << ", " << cities[i].y << ")\n";
	}
}

void TSP::print_thread_distances_EMST(int thread_id, int **emst_distances)
{
	// Print distances
	cout << endl;
	cout << "Thread " << thread_id << " distances:\n";
	for (int i = 0; i < cities_per_block; i++)
	{
		for (int j = 0; j < cities_per_block; j++)
		{
			cout << emst_distances[i][j] << " ";
		}
		cout << endl;
	}
}

void TSP::print_thread_distances_DP(int thread_id, std::vector<std::vector<double>> &dp_distances)
{
	// Print distances
	cout << endl;
	cout << "Thread " << thread_id << " distances:\n";
	for (int i = 0; i < cities_per_block; i++)
	{
		for (int j = 0; j < cities_per_block; j++)
		{
			cout << dp_distances[i][j] << " ";
		}
		cout << endl;
	}
}

void TSP::print_thread_paths(int num_threads)
{
	cout << endl;
	cout << "Results:" << endl;
	// Print thread_paths
	for (int i = 0; i < num_threads; i++)
	{
		cout << "Thread " << i << " path: ";
		for (int j = 0; j < thread_results[i].first.size(); j++)
		{
			cout << thread_results[i].first[j].id << " ";
		}
		cout << endl;
		cout << "Thread " << i << " length: " << thread_results[i].second << endl;
	}
}