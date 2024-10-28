#pragma once 

#include <vector> 
#include <tuple>
#include <iostream>

template <typename T>
using graph = std::vector<T>;

typedef std::tuple<int,int,int> edge;

template <typename T>
using directed_graph = std::vector<std::vector<T>>;

typedef std::pair<int,int> edge_map;


struct direct_flat_graph
{
    std::vector<int> nums;
    graph<edge> gr;
    
    direct_flat_graph(const graph<edge> &g, int n, int);

    //neighbourhood lists
    inline int first_id(int v) const
    {
        return nums[v-1];
    }

    inline int last_id(int v) const
    {
        return nums[v]-1;
    }
};

