//
// Created by carlosrey on 28/09/20.
//

#include "Solution.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <limits>
#include <cstddef>
#include "Solution.h"
#include "OPTUtils.h"
#include <math.h>
#include "utils.h"

using namespace std;

Solution::Solution(){}


//solution= matrix dimension m x n
Solution::Solution(vector<vector<int> > & solution){
    for(int i=0;i<solution.size();i++){
        vector<int> aux(solution[i]);
        this->solution.push_back(aux);
    }
    this->n= solution[0].size();
    this->m = solution.size();

}

Solution::Solution(vector<vector<double> > & solution){
    for(int i=0;i<solution.size();i++){
        vector<int> aux;
        for (int j = 0; j < solution[i].size(); j++)
        {
            aux.push_back((int)solution[i][j]);
        }
        this->solution.push_back(aux);
    }
    this->n= solution[0].size();
    this->m = solution.size();

}

Solution::~Solution( )
{

}

void Solution::show(){
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout<<solution[i][j]<<" ";
        }
        cout<<endl;
    }

    cout<<endl;
}
double Solution::check(int n ,int m,int * w ,int * C,int *p, int * Q){
    /**********************************************/
    /**********************************************/
    /**********************************************/
    if(solution.size()!=m){
        cout<<"Error m"<<endl;
        exit(0);
    }else{
        for(int i=0;i<solution.size();i++){
            if(solution[i].size()!=n){
                cout<<"Error n"<<endl;
                exit(0);
            }
        }
    }
    /**********************************************/
    /**********************************************/
    /**********************************************/
    vector<double> c_use;
    vector<double> total_use(n,0.0);
    for(int k=0;k<m;k++){
        double sum_w = 0.0;
        for(int i=0;i<n;i++){
            sum_w = sum_w + (w[i]*solution[k][i]);
            total_use[i] = total_use[i] + solution[k][i];
        }
        c_use.push_back(sum_w);
    }
    cout<<endl;
    cout<<"C ";
    for(int i=0;i<c_use.size();i++) cout<<c_use[i]<<"/"<<C[i]<<"\t";
    cout<<endl;


    for(int i=0;i<c_use.size();i++) {
        if(c_use[i] > C[i]){
            cout<<"Error C"<<endl;
            exit(0);
        }
    }

    cout<<endl;
    cout<<"Sum ";
    for(int i=0;i<total_use.size();i++) cout<<total_use[i]<<" ";
    cout<<endl;

    for(int i=0;i<total_use.size();i++) {
        if(total_use[i] > 1){
            cout<<"Error Number of obj"<<endl;
            exit(0);
        }
    }
    //Calculate the profit
    double acum =0.0;
    for(int k=0;k<m;k++){
        for(int i=0;i<n;i++){
            if(solution[k][i]==1){
                acum = acum + p[i];
                for(int j=i+1;j<n;j++){
                    int index_q =  ij2k(n,i,j);
                    acum = acum + ((solution[k][j])*Q[index_q]);
                }
            }
        }
    }
    cout<<"SumProfit "<<acum<<endl;
    return acum;

    /**********************************************/
    /**********************************************/
    /**********************************************/

}