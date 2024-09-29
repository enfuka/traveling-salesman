#ifndef THREADS_H
#define THREADS_H

#include <cstring>
#include <iostream>
#include <pthread.h>
#include "tsp.h"

using namespace std;

// Abstract Thread class
class Thread
{
private:
	// Internally track thread ids
	pthread_t _id;

	// Serves as global function for pthraed_create
	// Internally calls member function run
	static void *exec(void *thr);

	// Prevent copying
	Thread(const Thread &arg);

	// Prevent assignment
	Thread &operator=(const Thread &rhs);

protected:
	// Pure virtual function, to implement with code thread should run
	virtual void run() = 0;

public:
	// Constructor
	Thread()
	{
		mytsp = NULL;
		_id = 0;
		my_id = start_node = -1;
		algorithm = TSP::Algorithm::DP;
	};

	// Destructor
	virtual ~Thread() {};

	// Starts the internal thread
	void start();

	// Wait for internal thread to exit
	void join();

	// Get thread id of internal thread
	long get_tid() { return (long)_id; }

	// Index of node to start touring at
	int start_node;

	// distance matrix
	std::vector<std::vector<double>> distance;

	// Pointer to TSP object
	TSP *mytsp;

	// aLgorithm to use
	TSP::Algorithm algorithm;

	// Assigned ID (for tracking purposes)
	int my_id;
};

#endif
