#include "tsp.h"
#include "MyThread.h"

using namespace std;

struct thread_data
{
	int tid;
	TSP *tsp;
};
struct thread_data *threadData;

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

	cost = new int *[n];
	for (int i = 0; i < n; i++)
	{
		cost[i] = new int[n];
	}

	path_vals = new int *[n];
	for (int i = 0; i < n; i++)
	{
		path_vals[i] = new int[n];
	}

	// Adjacency lsit
	adjlist = new vector<int>[n];

	distance = std::vector<std::vector<double>>(n, std::vector<double>(n));

	thread_distance = std::vector<std::vector<std::vector<double>>>(THREADS, std::vector<std::vector<double>>(n, std::vector<double>(n)));
	// Vector to hold the mapping of regular city index to thread city index
	thread_city_index = std::vector<std::vector<int>>(THREADS, std::vector<int>(n));

	thread_city_mapping = std::vector<std::vector<City>>(THREADS);

	thread_results = std::vector<std::pair<std::vector<TSP::City>, int>>(THREADS);
};

TSP::~TSP()
{
	/////////////////////////////////////////////////////
	// Destructor
	/////////////////////////////////////////////////////

	for (int i = 0; i < n; i++)
	{
		delete[] graph[i];
		delete[] cost[i];
		delete[] path_vals[i];
	}
	delete[] path_vals;
	delete[] graph;
	delete[] cost;
	delete[] adjlist;
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

void *F(void *args)
{
	struct thread_data *my_data = (struct thread_data *)args;
	int tid = my_data->tid;
	TSP *tsp = my_data->tsp;
	int **graph = tsp->graph;
	int start, end;
	// start = START_AT(tid, THREADS, tsp->n);
	// end = END_AT(tid, THREADS, tsp->n);

	start = tsp->start_idx[tid];
	end = tsp->end_idx[tid];
	// cout << "thread " << setw(4) << left << tid << setw(8) << left << " start: " << setw(5) << left << start;
	// cout << setw(6) << left << " end: " << setw(5) << left << end << " load: " << end- start + 1 << endl;

	// clock_t t = clock();
	//  fill matrix with distances from every city to every other city
	for (int i = start; i <= end; i++)
	{
		for (int j = i; j < tsp->n; j++)
		{
			// Don't delete this line  it's supposed to be there.
			graph[i][j] = graph[j][i] = tsp->get_distance(tsp->cities[i], tsp->cities[j]);
		}
	}
	// t = clock() - t;
	// t = clock();
	// cout << "thread " << tid << " time: " << 1000*(((float)clock())/CLOCKS_PER_SEC) << " s"<< endl;
	pthread_exit(NULL);
}

void TSP::fillMatrix_threads()
{
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	int amount = (n / THREADS) * 0.2;
	int x = (n / THREADS) - amount; // min amount given to threads
	int rem = n - (x * THREADS);
	int half = THREADS / 2 + 1;

	int pos = 0;
	for (int i = 0; i < half; i++)
	{
		start_idx[i] = pos;
		pos += (x - 1);
		end_idx[i] = pos;
		pos++;
	}
	int remainingThreads = THREADS - half + 1;
	int extraToEach = rem / remainingThreads;
	// Divide remainer among second half of threads
	for (int i = half; i < THREADS; i++)
	{
		start_idx[i] = pos;
		pos += (x + extraToEach - 1);
		end_idx[i] = pos;
		pos++;
	}
	end_idx[THREADS - 1] = n - 1;

	int rc;
	void *status;
	threadData = new struct thread_data[n];

	// allocate space for n thread ids
	pthread_t *thread = new pthread_t[n];
	pthread_attr_t attr;

	// Initialize and set thread detached attribute
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (long t = 0; t < THREADS; t++)
	{
		// printf("Creating thread %ld\n", t);
		threadData[t].tid = t;
		threadData[t].tsp = this;
		rc = pthread_create(&thread[t], &attr, F, (void *)&threadData[t]);
		if (rc)
		{
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	// Free attribute and wait for the other threads
	pthread_attr_destroy(&attr);
	for (long t = 0; t < THREADS; t++)
	{
		rc = pthread_join(thread[t], &status);
		if (rc)
		{
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
		// printf("Completed join with thread %ld having a status of %ld\n",t,(long)status);
	}
	delete[] threadData;
}

// Function that takes a vector of cities and returns the cell id of each city
void TSP::create_2D_grid(vector<City> cities, int n, int grid_size, vector<int> &cell_ids)
{
	/////////////////////////////////////////////////////
	// Divide the space that contains all the cities into a square grid of cells
	/////////////////////////////////////////////////////
	// Find min and max x and y values
	int min_x = cities[0].x;
	int max_x = cities[0].x;
	int min_y = cities[0].y;
	int max_y = cities[0].y;
	for (int i = 1; i < n; i++)
	{
		if (cities[i].x < min_x)
			min_x = cities[i].x;
		if (cities[i].x > max_x)
			max_x = cities[i].x;
		if (cities[i].y < min_y)
			min_y = cities[i].y;
		if (cities[i].y > max_y)
			max_y = cities[i].y;
	}

	int x, y;
	int num_rows = grid_size;
	int num_cols = grid_size;

	// Find side length of each cell
	int min_side_length = max(max_x, max_y) - min(min_x, min_y);
	int cell_side_length = ceil((float)min_side_length / (float)grid_size);

	int range_floor = min(min_x, min_y);
	int range_ceil = range_floor + cell_side_length * grid_size;

	cout << "range_floor: " << range_floor << " range_ceil: " << range_ceil << endl;
	cout << "cell_side_length: " << cell_side_length << endl;

	for (int i = 0; i < n; i++)
	{
		x = cities[i].x;
		y = cities[i].y;
		int row = (int)((float)(x - range_floor) / (float)cell_side_length);
		int col = (int)((float)(y - range_floor) / (float)cell_side_length);
		cell_ids.push_back(row * grid_size + col);
	}

	// print out each cell, its boundaries, and the cities in it
	for (int i = 0; i < grid_size; i++)
	{
		for (int j = 0; j < grid_size; j++)
		{
			cout << "Cell " << i * grid_size + j << " (" << i << ", " << j << "): ";
			cout << "x: " << i * cell_side_length + range_floor << " to " << (i + 1) * cell_side_length + range_floor << ", ";
			cout << "y: " << j * cell_side_length + range_floor << " to " << (j + 1) * cell_side_length + range_floor << endl;
			for (int k = 0; k < n; k++)
			{
				if (cell_ids[k] == i * grid_size + j)
				{
					cout << k << " ";
					cout << "(" << cities[k].x << ", " << cities[k].y << ") ";
				}
			}
			cout << endl;
		}
	}
}

//================================ IRRELEVANT ================================//
void TSP::findMST_old()
{
	/////////////////////////////////////////////////////
	// In each iteration, we choose a minimum-weight
	// edge (u, v), connecting a vertex v in the set A to
	// the vertex u outside of set A
	/////////////////////////////////////////////////////
	int key[n];		// Key values used to pick minimum weight edge in cut
	bool in_mst[n]; // To represent set of vertices not yet included in MST
	int parent[n];

	// For each vertex v in V
	for (int v = 0; v < n; v++)
	{
		// Initialize all keys to infinity
		key[v] = std::numeric_limits<int>::max();

		// Mark as not being in mst yet
		in_mst[v] = false;
	}

	// Node 0 is the root node so give it the lowest distance (key)
	key[0] = 0;
	parent[0] = -1; // First node is always root of MST

	for (int i = 0; i < n - 1; i++)
	{
		// Find closest remaining (not in tree) vertex
		// TO DO : This would be better represented by heap/pqueue
		int v = minKey(key, in_mst);

		// Add vertex v to the MST
		in_mst[v] = true;

		// Look at each vertex u adjacent to v that's not yet in mst
		for (int u = 0; u < n; u++)
		{
			if (graph[v][u] && in_mst[u] == false && graph[v][u] < key[u])
			{
				// Update parent index of u
				parent[u] = v;

				// Update the key only if dist is smaller than key[u]
				key[u] = graph[v][u];
			}
		}
	}

	// map relations from parent array onto matrix
	for (int v1 = 0; v1 < n; v1++)
	{
		// there is an edge between v1 and parent[v1]
		int v2 = parent[v1];
		if (v2 != -1)
		{
			adjlist[v1].push_back(v2);
			adjlist[v2].push_back(v1);
		}
	}
};
// findMST helper function
int TSP::minKey(int key[], bool mstSet[])
{
	// Initialize min value
	int min = std::numeric_limits<int>::max();
	int min_index;
	for (int v = 0; v < n; v++)
		if (mstSet[v] == false && key[v] < min)
		{
			min = key[v];
			min_index = v;
		}
	return min_index;
};

void TSP::findOdds()
{
	/////////////////////////////////////////////////////
	// Find nodes with odd degrees in T to get subgraph O
	/////////////////////////////////////////////////////

	// store odds in new vector for now
	for (int r = 0; r < n; r++)
	{
		// cities[r].isOdd = ((adjlist[r].size() % 2) == 0) ? 0 : 1;
		if ((adjlist[r].size() % 2) != 0)
		{
			odds.push_back(r);
		}
	}
}

void TSP::perfect_matching()
{
	/////////////////////////////////////////////////////
	// find a perfect matching M in the subgraph O using greedy algorithm
	// but not minimum
	/////////////////////////////////////////////////////
	int closest, length; // int d;
	std::vector<int>::iterator tmp, first;

	// Find nodes with odd degrees in T to get subgraph O
	findOdds();

	// for each odd node
	while (!odds.empty())
	{
		first = odds.begin();
		vector<int>::iterator it = odds.begin() + 1;
		vector<int>::iterator end = odds.end();
		length = std::numeric_limits<int>::max();
		for (; it != end; ++it)
		{
			// if this node is closer than the current closest, update closest and length
			if (graph[*first][*it] < length)
			{
				length = graph[*first][*it];
				closest = *it;
				tmp = it;
			}
		} // two nodes are matched, end of list reached
		adjlist[*first].push_back(closest);
		adjlist[closest].push_back(*first);
		odds.erase(tmp);
		odds.erase(first);
	}
}

void TSP::euler(int pos, vector<int> &path)
{
	/////////////////////////////////////////////////////////
	// Based on this algorithm:
	//	http://www.graph-magics.com/articles/euler.php
	// we know graph has 0 odd vertices, so start at any vertex
	// O(V+E) complexity
	/////////////////////////////////////////////////////////

	// make copy of original adjlist to use/modify
	vector<int> *temp = new vector<int>[n];
	for (int i = 0; i < n; i++)
	{
		temp[i].resize(adjlist[i].size());
		temp[i] = adjlist[i];
	}

	path.clear();

	// Repeat until the current vertex has no more neighbors and the stack is empty.
	stack<int> stk;
	while (!stk.empty() || temp[pos].size() > 0)
	{
		// If current vertex has no neighbors -
		if (temp[pos].size() == 0)
		{
			// add it to circuit,
			path.push_back(pos);
			// remove the last vertex from the stack and set it as the current one.
			int last = stk.top();
			stk.pop();
			pos = last;
		}
		// Otherwise (in case it has neighbors)
		else
		{
			// add the vertex to the stack,
			stk.push(pos);
			// take any of its neighbors,
			int neighbor = temp[pos].back();
			// remove the edge between selected neighbor and that vertex,
			temp[pos].pop_back();
			for (unsigned int i = 0; i < temp[neighbor].size(); i++)
				if (temp[neighbor][i] == pos)
				{													  // find position of neighbor in list
					temp[neighbor].erase(temp[neighbor].begin() + i); // remove it
					break;
				}
			// and set that neighbor as the current vertex.
			pos = neighbor;
		}
	}
	path.push_back(pos);
}

void TSP::make_hamilton(vector<int> &path, int &path_dist)
{
	// remove visited nodes from Euler tour
	bool visited[n]; // boolean value for each node if it has been visited yet
	memset(visited, 0, n * sizeof(bool));

	path_dist = 0;

	int root = path.front();
	vector<int>::iterator curr = path.begin();
	vector<int>::iterator next = path.begin() + 1;
	visited[root] = true;

	// loop until the end of the circuit list is reached
	while (next != path.end())
	{
		// if we haven't been to the next city yet, go there
		if (!visited[*next])
		{
			path_dist += graph[*curr][*next];
			curr = next;
			visited[*curr] = true;
			next = curr + 1;
		}
		else
		{
			next = path.erase(next); // remove it
		}
	}

	// add the distance back to the root
	path_dist += graph[*curr][*next];
}

void TSP::create_tour(int pos)
{
	// call euler with actual circuit vector
	euler(pos, circuit);

	// make it hamiltonian
	// pass actual vars
	make_hamilton(circuit, pathLength);
}
// Does euler and hamilton but doesn't modify any variables
// Just finds path length from the node specified and returns it
int TSP::find_best_path(int pos)
{

	// create new vector to pass to euler function
	vector<int> path;
	euler(pos, path);

	// make it hamiltonian, pass copy of vars
	int length;
	make_hamilton(path, length);

	// Optimize
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);
	twoOpt(graph, path, length, n);

	return length;
}
//================================ ^^ IRRELEVANT ^^ ================================//

void TSP::DP()
{
	/////////////////////////////////////////////////////
	// Dynamic programming to solve TSP
	/////////////////////////////////////////////////////
	// Create a copy of the graph

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			distance[i][j] = graph[i][j];
	}

	DPSolver dp(0, distance, 0);
	std::list<int> tour = dp.getTour();
	pathLength = dp.getTourCost();
	circuit.clear();
	for (std::list<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
		circuit.push_back(*it);
}

void TSP::parallel_DP(int num_threads)
{
	// Find the number of blocks in the grid
	int num_blocks = num_threads;

	// find grid side length
	int grid_side = sqrt(num_blocks);

	// Map each city to a cell
	for (int i = 0; i < n; i++)
	{
		thread_city_mapping[cell_ids[i]].push_back(cities[i]);
	}

	// Print mapping
	for (int i = 0; i < num_blocks; i++)
	{
		cout << "Block " << i << ":\n";
		for (int j = 0; j < thread_city_mapping[i].size(); j++)
		{
			cout << thread_city_mapping[i][j].id << " ";
		}
		cout << endl;
	}

	// Create array of thread objects
	MyThread threads[num_threads];

	// Create a seperate distance matrix for each thread including only the cities in the thread's cell
	// While recording the mapping of regular city index to thread city index

	for (int i = 0; i < num_threads; i++)
	{
		for (int j = 0; j < thread_city_mapping[i].size(); j++)
		{
			for (int k = 0; k < thread_city_mapping[i].size(); k++)
			{
				thread_distance[i][j][k] = graph[thread_city_mapping[i][j].id][thread_city_mapping[i][k].id];
				thread_city_index[i][thread_city_mapping[i][j].id] = j;
			}
		}
	}

	// Print distance matrix for each thread
	for (int i = 0; i < num_threads; i++)
	{
		cout << "Thread " << i << ":\n";
		for (int j = 0; j < thread_city_mapping[i].size(); j++)
		{
			for (int k = 0; k < thread_city_mapping[i].size(); k++)
			{
				cout << thread_distance[i][j][k] << " ";
			}
			cout << endl;
		}
	}

	// Create a thread for each block
	for (int i = 0; i < num_threads; i++)
	{
		threads[i].mytsp = this;
		threads[i].my_id = i;
		threads[i].distance = thread_distance[i];
		threads[i].start();
	}

	// Wait for all the threads
	for (int i = 0; i < num_threads; i++)
	{
		threads[i].join();
	}
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
		sleep(2);
	}

	// Wait for all the threads
	for (int i = 0; i < num_threads; i++)
	{
		threads[i].join();
	}

	// Print thread_paths_EMST
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
	int cell_side_length = 100;
	int cities_per_cell = 10;
	// seed random number generator

	srand(thread_id);

	// Create a new vector of cities for the thread
	vector<City> cities;
	for (int i = 0; i < cities_per_cell; i++)
	{

		// Create a random city in the cell, have relative id's
		int x = rand() % cell_side_length + thread_id * cell_side_length;
		int y = rand() % cell_side_length + thread_id * cell_side_length;
		City city = {thread_id * cities_per_cell + i, x, y};
		cities.push_back(city);
	}

	// Print cities
	cout << "Thread " << thread_id << " cities:\n";
	for (int i = 0; i < cities_per_cell; i++)
	{
		cout << cities[i].id << " (" << cities[i].x << ", " << cities[i].y << ")\n";
	}

	// Calculate the distance between each pair of cities
	int **emst_distances = new int *[cities_per_cell];
	for (int i = 0; i < cities_per_cell; i++)
	{
		emst_distances[i] = new int[cities_per_cell];
		for (int j = 0; j < cities_per_cell; j++)
		{
			emst_distances[i][j] = get_distance(cities[i], cities[j]);
		}
	}

	// Print distances
	cout << "Thread " << thread_id << " distances:\n";
	for (int i = 0; i < cities_per_cell; i++)
	{
		for (int j = 0; j < cities_per_cell; j++)
		{
			cout << emst_distances[i][j] << " ";
		}
		cout << endl;
	}

	// Find the path and length for the thread
	std::vector<int> tour;
	int pathLength;
	EMSTSolver emst(cities_per_cell, emst_distances, 1);
	tour = emst.find_path();
	pathLength = emst.get_cost();

	// Map tour to city objects, skip the last city
	std::vector<TSP::City> tour_with_coordinates;
	for (std::vector<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
	{
		tour_with_coordinates.push_back(cities[*it]);
	}

	// Record the path and length in the right index in thread_paths_EMST
	thread_results[thread_id] = make_pair(tour_with_coordinates, pathLength);
}

void TSP::openTSP_DP(int thread_id)
{
	// Find how many cities are in the thread's cell
	int num_cities = thread_city_mapping[thread_id].size();

	if (num_cities == 0 || num_cities == 1)
	{
		cout << "Thread " << thread_id << " has no cities\n";
		return;
	}

	// Print size of cell
	cout << "Thread " << thread_id << " num cities: " << num_cities << endl;

	// Create a new distance matrix for the thread with the correct size
	vector<vector<double>> correct_distances(num_cities, vector<double>(num_cities));

	// Fill the distance matrix with the correct distances from thread_distance
	for (int i = 0; i < num_cities; i++)
	{
		for (int j = 0; j < num_cities; j++)
		{
			correct_distances[i][j] = thread_distance[thread_id][i][j];
		}
	}

	// Print correct distances
	cout << "Thread " << thread_id << " correct distances:\n";
	for (int i = 0; i < num_cities; i++)
	{
		for (int j = 0; j < num_cities; j++)
		{
			cout << correct_distances[i][j] << " ";
		}
		cout << endl;
	}

	std::list<int> tour;

	// Run open dp for the thread
	if (num_cities > 2)
	{
		DPSolver dp(0, correct_distances, 1);
		tour = dp.getTour();
		pathLength = dp.getTourCost();
	}
	else if (num_cities == 2)
	{
		pathLength = correct_distances[0][1];
		tour = {0, 1};
	}

	// // Map the tour back to the original city indices
	// std::list<int> correct_tour;
	// for (std::list<int>::iterator it = tour.begin(); it != tour.end(); ++it)
	// {
	// 	correct_tour.push_back(thread_city_mapping[thread_id][*it].id);
	// }

	// // Record the path and length in the right index
	// thread_paths[thread_id] = make_pair(correct_tour, pathLength);
	// // Print the path and length
	// cout << "Thread " << thread_id << " path: ";
	// for (std::list<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
	// 	cout << *it << " ";
	// cout << tour.back() << endl;
	// cout << "Thread " << thread_id << " length: " << pathLength << endl;
}

void TSP::emst()
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
	int cost = emst.get_cost();

	circuit.clear();
	for (std::vector<int>::iterator it = tour.begin(); it != std::prev(tour.end()); ++it)
		circuit.push_back(*it);
}

void TSP::make_shorter()
{
	// Modify circuit & pathLength
	twoOpt(graph, circuit, pathLength, n);
}

//================================ PRINT FUNCTIONS ================================//

void TSP::printResult()
{
	ofstream outputStream;
	outputStream.open(outFname.c_str(), ios::out);
	outputStream << pathLength << endl;
	for (vector<int>::iterator it = circuit.begin(); it != circuit.end(); ++it)
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
	for (vector<int>::iterator it = circuit.begin(); it != circuit.end() - 1; ++it)
	{
		cout << *it << " to " << *(it + 1) << " ";
		cout << graph[*it][*(it + 1)] << endl;
	}
	cout << *(circuit.end() - 1) << " to " << circuit.front() << " ";
	cout << graph[*(circuit.end() - 1)][circuit.front()] << endl;

	cout << "\nLength: " << pathLength << endl
		 << endl;
};

void TSP::printEuler()
{
	for (vector<int>::iterator it = circuit.begin(); it != circuit.end(); ++it)
		cout << *it << endl;
}

void TSP::printAdjList()
{
	for (int i = 0; i < n; i++)
	{
		cout << i << ": "; // print which vertex's edge list follows
		for (int j = 0; j < (int)adjlist[i].size(); j++)
		{
			cout << adjlist[i][j] << " "; // print each item in edge list
		}
		cout << endl;
	}
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
