
#include <bits/stdc++.h>
#include "par_boruvka_openmp.h"
using namespace std;

typedef std::chrono::milliseconds time_type;

class EMSTSolver
{
private:
    int V; // Number of vertices in the graph
    int **distance_matrix;
    int is_open;

public:
    // Dynamic array to store the final answer
    vector<int> final_ans;
    vector<vector<int>> v;

    // Constructor
    EMSTSolver(int V, int **graph, int is_open)
    {
        this->V = V;
        this->distance_matrix = graph;
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
            delete[] distance_matrix[i];
        }
        delete[] distance_matrix;
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
                if (distance_matrix[u][v] && mstSet[v] == false && distance_matrix[u][v] < key[v])
                    parent[v] = u, key[v] = distance_matrix[u][v];
        }
        vector<vector<int>> v;
        v = MST(parent);
        return v;
    }

    vector<vector<int>> boruvkaMST()
    {
        // Map distance matrix to graph
        graph<edge> edge_vector;
        for (int i = 0; i < V; i++)
        {
            for (int j = i + 1; j < V; j++)
            {
                edge_vector.push_back({i + 1, j + 1, distance_matrix[i][j]});
            }
        }

        // // print the edges
        // for (auto e : edge_vector)
        // {
        //     cout << get<0>(e) << " " << get<1>(e) << " " << get<2>(e) << endl;
        // }

        int NTHREADS = 4;
        auto start = std::chrono::high_resolution_clock::now();
        direct_flat_graph dfg(edge_vector, V, NTHREADS);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<time_type>(stop - start);

        cout << "Time taken to create the graph: " << duration.count() << "ms" << endl;

        // // print the dfg.gr
        // for (int i = 0; i < dfg.gr.size(); i++)
        // {
        //     cout << get<0>(dfg.gr[i]) << " " << get<1>(dfg.gr[i]) << " " << get<2>(dfg.gr[i]) << endl;
        // }

        // // print nums
        // for (int i = 1; i <= V; i++)
        // {
        //     cout << dfg.first_id(i) << " " << dfg.last_id(i) << endl;
        // }

        graph<edge> mst_result = boruvka_mst_par_openmp(dfg, V, NTHREADS);

        // Convert the result to the format expected by the TSP solver
        vector<vector<int>> mst;
        for (auto e : mst_result)
        {
            vector<int> edge;
            edge.push_back(get<0>(e) == 0 ? 0 : get<0>(e) - 1);
            edge.push_back(get<1>(e) == 0 ? 0 : get<1>(e) - 1);
            mst.push_back(edge);
            edge.clear();
        }

        // // Print the MST
        // for (auto e : mst)
        // {
        //     cout << e[0] << " " << e[1] << endl;
        // }

        return mst;
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
            cost += distance_matrix[path[i]][path[i + 1]];
        }
        return cost;
    }

    vector<int> find_path()
    {
        // getting the output as MST
        //  v = primMST();
        v = boruvkaMST();

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