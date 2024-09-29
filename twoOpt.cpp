#include "tsp.h"

using namespace std;

// Input: edge 1's v, edge 2's u
// Remove edge 1 and edge 2, reconnect using new path
void TSP::reverse(vector<City> &path, int start, int end, int n)
{
	while (end - start > 0)
	{
		// Start, end is relative value to the start,
		// the index is start|slut % size
		City temp = path[start % n];
		path[start % n] = path[end % n];
		path[end % n] = temp;
		start++;
		end--;
	}
}

int TSP::is_path_shorter(int **graph, City v1, City v2, City v3, City v4, int &total_dist)
{
	if ((get_distance(graph, v1, v3) + get_distance(graph, v2, v4)) < (get_distance(graph, v1, v2) + get_distance(graph, v3, v4)))
	{
		total_dist -= (get_distance(graph, v1, v2) + get_distance(graph, v3, v4) - get_distance(graph, v1, v3) - get_distance(graph, v2, v4));

		return 1;
	}
	return 0;
}

// Non-looping version of two-opt heuristic
int TSP::twoOpt(int **graph, vector<City> &path, int &pathLength, int n)
{
	int counter = 0;
	int term_cond = 5;
	int old_distance = pathLength;
	// int size = path.size();
	int v1, v2, u1, u2;

	// Iterate over each edge
	for (int i = 0; i < n; i++)
	{
		// first edge
		u1 = i;
		v1 = (i + 1) % n; // wrap around to first node if u1 is last node

		// Skip adjacent edges, start with node one past v1
		for (int j = i + 2; (j + 1) % n != i; j++)
		{
			// mod by length to go back to beginning
			u2 = j % n;
			v2 = (j + 1) % n;

			// Check if new edges would shorten the path length
			// If so, decreases pathLength
			if (is_path_shorter(graph, path[u1], path[v1], path[u2],
								path[v2], pathLength))
			{

				// Swap u1--v1 and u2--v2
				reverse(path, i + 1, j, n); // v1, u2

				if (pathLength - old_distance < 5 && counter == term_cond)
					break;

				// reset i
				if (pathLength - old_distance > term_cond)
					i = 0;

				old_distance = pathLength;
				counter++;
			}
		}
	}
	return pathLength;
}

int TSP::get_path_length(int **graph, vector<City> &path, int size)
{
	int length = 0;
	for (int i = 0; i < size - 1; i++)
	{
		length += get_distance(graph, path[i], path[i + 1]);
	}
	length += get_distance(graph, path[size - 1], path[0]);
	return length;
}

int TSP::get_distance(int **graph, City c1, City c2)
{
	/////////////////////////////////////////////////////
	// Calculate distance between c1 and c2
	/////////////////////////////////////////////////////

	// if edge exists in graph, return it
	if (graph[c1.id][c1.id])
		return graph[c1.id][c1.id];
	else
	{

		int dx = pow((float)(c1.x - c2.x), 2);
		int dy = pow((float)(c1.y - c2.y), 2);
		int distance = (floor((float)(sqrt(dx + dy)) + 0.5));
		graph[c1.id][c2.id] = graph[c2.id][c1.id] = distance;
		return distance;
	}
};