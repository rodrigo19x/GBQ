//
// Created by carlosrey on 28/09/20.
//

#ifndef GREEDYBASED_QMKP_01QPCPX2_H
#define GREEDYBASED_QMKP_01QPCPX2_H


#include "OPTUtils.h"
#include <math.h>
#include "utils.h"
//#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ilcplex/ilocplex.h>

// End magic tricks
/* This method blah blah blah ...*/

int solveQMKP_01qp_CPX2( int n, int m, int *p, int *Q, int *w, int *C, vector<vector<double> > x_con,double * objval, int TL, bool intflag,  char *logfile);
//double QKP( int n, int *p, int **Q, int *w, int C);

#endif //GREEDYBASED_QMKP_01QPCPX2_H
