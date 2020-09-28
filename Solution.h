//
// Created by carlosrey on 28/09/20.
//

#ifndef GREEDYBASED_SOLUTION_H
#define GREEDYBASED_SOLUTION_H


#include <iostream>
#include <vector>
#include <fstream>
#include <assert.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>


using namespace std;

class Solution {
public:
    vector<vector<int> >solution;
    int n;
    int m;
    Solution();
    Solution(vector<vector<int> > & solution);
    Solution(vector<vector<double> > & solution);
    double check(int n ,int m,int * w ,int * C,int *p, int * Q);
    void show();
    virtual ~Solution( );


};

#endif //GREEDYBASED_SOLUTION_H
