#include "par_boruvka_openmp.h"
#include "dsu.h"

#include <omp.h>
#include <chrono>
#include <climits>

graph<edge> boruvka_mst_par_openmp(const direct_flat_graph &dfg, int n, const int OMP_THREADS_NUM)
{
    omp_set_num_threads(OMP_THREADS_NUM);
    omp_lock_t locks[n + 1];
    for (int i = 1; i <= n; ++i)
        omp_init_lock(locks + i);

    graph<edge> mst;
    mst.resize(n - 1);
    make_set(n);

    std::vector<std::pair<int, int>> min(n + 1, {-1, -1});
    std::vector<int> processed(n + 1);
    int iters = 0;

    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();

    auto init = (stop - start);
    auto find = (stop - start);
    auto uni = (stop - start);

    while (components_num() > 1)
    {
        start = std::chrono::high_resolution_clock::now();

        min.assign(n + 1, {-1, -1});
        processed.assign(n + 1, 0);

        stop = std::chrono::high_resolution_clock::now();
        init += (stop - start);
        start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic, CHUNK) shared(dfg)
        for (int i = 1; i <= n; ++i)
        {
            int com = atomic_find_set(i);

            int minimal = INT_MAX;
            int min_id = -1;

            for (int j = dfg.first_id(i); j <= dfg.last_id(i); ++j)
            {
                auto [w, v, weight] = dfg.gr[j];
                v = atomic_find_set(v);

                if (v == com)
                    continue;
                if (weight < minimal)
                {
                    min_id = j;
                    minimal = weight;
                }
            }
            if (min_id == -1)
                continue;

            omp_set_lock(locks + com);

            std::pair<int, int> p = min[com];
            if (p.first == -1 || p.first > minimal)
                min[com] = {minimal, min_id};
            omp_unset_lock(locks + com);
        }

#pragma omp barrier
        stop = std::chrono::high_resolution_clock::now();
        find += (stop - start);
        start = std::chrono::high_resolution_clock::now();

        for (int i = 1; i <= n; ++i)
        {
            if (min[i].second != -1)
            {
                auto [v, w, weight] = dfg.gr[min[i].second];
                int cmp1 = find_set(v);
                int cmp2 = find_set(w);
                if (cmp1 != cmp2)
                {
                    union_sets(cmp1, cmp2);
                    mst.push_back(dfg.gr[min[i].second]);
                }
            }
        }
        stop = std::chrono::high_resolution_clock::now();
        uni += (stop - start);
    }

    // std::cout<< "INIT:"<<  std::chrono::duration_cast<std::chrono::milliseconds>(init).count()<<std::endl;
    // std::cout<< "FIND:" <<std::chrono::duration_cast<std::chrono::milliseconds>(find).count()<<std::endl;
    // std::cout<< "UNION:"<<std::chrono::duration_cast<std::chrono::milliseconds>(uni).count()<<std::endl;

    return mst;
}
