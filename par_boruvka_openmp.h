#pragma once

#include <vector>
#include "graph_types.h"

constexpr int CHUNK = 16;
// returns vector of edges type (from, to, weight)
graph<edge> boruvka_mst_par_openmp(const direct_flat_graph&, int,int);

