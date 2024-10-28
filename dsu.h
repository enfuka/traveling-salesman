#pragma once

#include <vector>

// makes set with index from 1 to size
void make_set(int size);
int find_set(int v);
int atomic_find_set(int v);

void union_sets(int a, int b);
int components_num();
