//==================================================================
// File			: tsp.h
// Author		: rsagalyn
// Date			: Aug 18, 2013
// Description	:
//==================================================================
#ifndef MWM_H_
#define MWM_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <queue>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>

#include <mpi.h>

#include "twoOpt.h"
#include "EMSTSolver.h"
#include "DPSolver.h"

using namespace std;

// Toggle printing debugging info to console
#define DEBUG 0

class TSP
{
private:
	// x and y coords of a node
	struct City
	{
		int id;
		int x;
		int y;
	};

	// Filename supplied on command line to read from
	string inFname;

	// Program-generated filename to output to
	string outFname;

	// Initialization function
	void getNodeCount();

protected:
public:
	// Number of nodes
	int n;

	int thread_count;
	int total_city_count;
	int block_side_length;
	int cities_per_block;

	// Algorithm type
	enum Algorithm
	{
		DP,
		EMST
	};

	// Final tour
	vector<int> final_tour;

	// Store cities and coords read in from file
	vector<City> cities;

	// Full n x n cost matrix of distances between each city
	int **graph;

	// Shortest path length
	int pathLength;

	std::vector<std::pair<std::vector<TSP::City>, int>> thread_results;

	// Constructor for serial TSP
	TSP(string in, string out);

	// Constructor for parallel TSP
	TSP(int total_city_count, int thread_count, int block_side_length);

	// Destructor
	~TSP();

	// Calculate distance
	int get_distance(struct City c1, struct City c2);

	// Initialization functions
	void readCities();

	void calculate_distances();

	// Create tour starting at specified node
	void create_tour(int);

	void serial_DP();

	void serial_EMST();

	void openTSP_DP(int thread_id);

	void openTSP_EMST(int thread_id);

	void parallel_solver(int num_threads, Algorithm algorithm);

	void thread_populate_block(int thread_id, std::vector<TSP::City> &cities);

	void dump_all_cities();

	void print_thread_distances_EMST(int thread_id, int **emst_distances);

	void print_thread_cities(int thread_id, std::vector<TSP::City> &cities);

	void print_thread_distances_DP(int thread_id, std::vector<std::vector<double>> &dp_distances);

	void print_thread_paths(int num_threads);

	// Calls twoOpt function
	void make_shorter();

	// Debugging functions
	void printCities();
	void printDistanceGraph();
	void MPI_solver();
	void printResult();
	void printPath();

	// Get node count
	int get_size() { return n; };
};

#endif /* MWM_H_ */
