//
// Created by carlosrey on 28/09/20.
//

#ifndef GREEDYBASED_UTILS_H
#define GREEDYBASED_UTILS_H

#include <iostream>
#include <vector>
#include <algorithm>    // std::sort
#include <map>
#include "OPTUtils.h"
#include <math.h>

using namespace std;
int ij2k(int n, int i, int j);
int ij2k2(int n, int i, int j);
void report(int n, int m, int *p, int *Q, int *w, int *C,vector< unsigned> variables);

#endif //GREEDYBASED_UTILS_H
