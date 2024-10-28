#include "dsu.h"

std::vector<int> parent,rank;
volatile int numsets;

void make_set(int size) 
{
    parent.resize(size+1);
    rank.resize(size+1);
    for(int i=1; i<=size;++i)
    {
        parent[i] = i;
        rank[i] = 0;
    }
    numsets = size;
}

int find_set(int v) 
{
    if (v != parent[v]) 
    {
        parent[v] = find_set(parent[v]);
    }
    return parent[v];
}

int atomic_find_set(int v)
{
    int oldpar = parent[v];
    if(v == oldpar) return v;
    int par = atomic_find_set(parent[v]);
    __sync_val_compare_and_swap(&parent[v],oldpar,par);
    return parent[v];
}

void union_sets(int a, int b) 
{
    a = find_set(a);
    b = find_set(b);
    if (a != b) 
    {
        if (rank[a] < rank[b])
            std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b])
            rank[a]++;
        --numsets;
    }
}


int components_num()
{
    return numsets;
}