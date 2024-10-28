#include "graph_types.h"
#include <omp.h>

direct_flat_graph::direct_flat_graph(const graph<edge> &g, int n, const int THREADS_NUM)
{
    omp_lock_t locks[n+1];
    for(int i=1;i<=n;++i)
        omp_init_lock(locks+i);    

    nums.assign(n+1,0);
    for(auto [v,w,weight] : g)
    {
        ++nums[v];
        ++nums[w];
    }
    for(int i=1;i<=n;++i)
    {
        nums[i] = nums[i]+nums[i-1];
    }


    std::vector<int> positions = nums;
    gr.resize(nums.back());
    #pragma omp parallel for num_threads(THREADS_NUM)
    for(auto [v,w,weight] : g)
    {
        omp_set_lock(locks+v); 
        gr[positions[v-1]++]={v,w,weight};
        omp_unset_lock(locks+v); 

        omp_set_lock(locks+w); 
        gr[positions[w-1]++]={w,v,weight};
        omp_unset_lock(locks+w); 

    }
    for(int i=1;i<=n;++i)
        omp_destroy_lock(locks+i);
    

}

