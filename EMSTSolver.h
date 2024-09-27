
#include <bits/stdc++.h>
using namespace std;

class EMSTSolver
{
private:
    // Number of vertices in the graph
    int V;
    int **graph;
    int is_open;

public:
    // Dynamic array to store the final answer
    vector<int> final_ans;
    vector<vector<int>> v;

    // Constructor
    EMSTSolver(int V, int **graph, int is_open)
    {
        this->V = V;
        this->graph = graph;
        this->is_open = is_open;

        // if open, set all distances to node 0 to 0
        if (is_open)
        {
            for (int i = 0; i < V; i++)
            {
                graph[i][0] = 0;
            }
        }
    }

    // Destructor
    ~EMSTSolver()
    {
        for (int i = 0; i < V; i++)
        {
            delete[] graph[i];
        }
        delete[] graph;
    }

    // Function to find the vertex with minimum key value
    int minimum_key(int key[], bool mstSet[])
    {
        int min = INT_MAX, min_index;

        for (int v = 0; v < V; v++)
            if (mstSet[v] == false && key[v] < min)
                min = key[v], min_index = v;

        return min_index;
    }

    vector<vector<int>> MST(int parent[])
    {
        vector<vector<int>> v;
        for (int i = 1; i < V; i++)
        {
            vector<int> p;
            p.push_back(parent[i]);
            p.push_back(i);
            v.push_back(p);
            p.clear();
        }
        return v;
    }

    // getting the Minimum Spanning Tree from the given graph
    // using Prim's Algorithm
    vector<vector<int>> primMST()
    {
        int parent[V];
        int key[V];

        // to keep track of vertices already in MST
        bool mstSet[V];

        // initializing key value to INFINITE & false for all mstSet
        for (int i = 0; i < V; i++)
            key[i] = INT_MAX, mstSet[i] = false;

        // picking up the first vertex and assigning it to 0
        key[0] = 0;
        parent[0] = -1;

        // The Loop
        for (int count = 0; count < V - 1; count++)
        {
            // checking and updating values wrt minimum key
            int u = minimum_key(key, mstSet);
            mstSet[u] = true;
            for (int v = 0; v < V; v++)
                if (graph[u][v] && mstSet[v] == false && graph[u][v] < key[v])
                    parent[v] = u, key[v] = graph[u][v];
        }
        vector<vector<int>> v;
        v = MST(parent);
        return v;
    }

    // getting the preorder walk of the MST using DFS
    void DFS(int **edges_list, int num_nodes, int starting_vertex, bool *visited_nodes)
    {
        // adding the node to final answer
        final_ans.push_back(starting_vertex);

        // checking the visited status
        visited_nodes[starting_vertex] = true;

        // using a recursive call
        for (int i = 0; i < num_nodes; i++)
        {
            if (i == starting_vertex)
            {
                continue;
            }
            if (edges_list[starting_vertex][i] == 1)
            {
                if (visited_nodes[i])
                {
                    continue;
                }
                DFS(edges_list, num_nodes, i, visited_nodes);
            }
        }
    }

    int tour_cost(vector<int> path)
    {
        int cost = 0;
        for (int i = 0; i < path.size() - 1; i++)
        {
            cost += graph[path[i]][path[i + 1]];
        }
        return cost;
    }

    vector<int> find_path()
    {
        // getting the output as MST
        v = primMST();

        // creating a dynamic matrix
        int **edges_list = new int *[V];
        for (int i = 0; i < V; i++)
        {
            edges_list[i] = new int[V];
            for (int j = 0; j < V; j++)
            {
                edges_list[i][j] = 0;
            }
        }

        // setting up MST as adjacency matrix
        for (int i = 0; i < v.size(); i++)
        {
            int first_node = v[i][0];
            int second_node = v[i][1];
            edges_list[first_node][second_node] = 1;
            edges_list[second_node][first_node] = 1;
        }

        // a checker function for the DFS
        bool *visited_nodes = new bool[V];
        for (int i = 0; i < V; i++)
        {
            bool visited_node;
            visited_nodes[i] = false;
        }

        // performing DFS
        DFS(edges_list, V, 0, visited_nodes);

        // adding the source node to the path
        final_ans.push_back(final_ans[0]);

        return final_ans;
    }

    // Function to find the cost of the path
    int get_cost()
    {
        return tour_cost(final_ans);
    }
};