/*
 * utils.cpp
 *
 *  Created on: January 11, 2017
 * 	Last Update: February 22 2017
 *      Author: laura
 */

/* An IQP instance looks as follows:
 * min x^T Q x + c^T x
 * s.t. x \in {1,...,u}^n
 * */

#include "utils.h"

//#define DEBUG_UTILS

int ij2k(int n, int i, int j) {
/*
 * Returns linear index of upper triangular matrix, using this formula:
 * k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
 * * */


    if (i == j) {
        cout << "\nerror: ij2k wrong index i==j";
        exit(1);
    }

    int k;

    if (j < i) k = (n*(n-1)/2) - (n-j)*((n-j)-1)/2 + i - j - 1;
    else k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;

    return k;
}




int ij2k2(int n, int i, int j){
    if (i == j) {
        cout << "\nerror: ij2k wrong index i==j";
        exit(1);
    }
    int k;
    if (i < j) k = i*(n-1)+j-1;
    else k = i*(n-1)+j;

    return k;
}


void report(int n, int m, int *p, int *Q, int *w, int *C,vector<unsigned> variables){

    std::map<int,vector<int> >::iterator it ;
    std::map<int,vector<int> > var_knapsack;
    std::map<int,vector<int>>  knapsack_var;
    vector<int> real_variables;

    std::vector< std::vector<int> > Q_2;
    for(int i=0;i<n;i++){
        std::vector<int> q_aux(n,-1);
        Q_2.push_back(q_aux);
    }

    int t = 0;
    /*
    for (int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
             Q_2[i][j] = Q[t];
             t++;
        }

    }*/

    for(int i=0;i<m;i++){
        vector<int> v ;
        knapsack_var.insert( std::pair<int,vector<int> >(i,v) );
    }

    for(int i=0;i<variables.size();i++){
        //cout<<variables[i]<<" ";

        int var = variables[i] % n;
        real_variables.push_back(var);
        vector<int> v ;
        //cout<<"["<< variables[i]<<" "<< var<<"]"<<endl;
        var_knapsack.insert ( std::pair<int,vector<int> >(var,v) );
    }
    //cout<<endl;

    //cout<<endl;
    for(int i=0;i<variables.size();i++){
        int var = variables[i] % n;
        int knap = variables[i] / n;
        //cout<<variables[i]<<" "<<var<<" "<<knap<<endl;
        var_knapsack[var].push_back(knap);
        knapsack_var[knap].push_back(var);
    }
    //cout<<endl;
    /*
    std::cout << "var_knapsackContains:\n";
    for (it=var_knapsack.begin(); it!=var_knapsack.end(); ++it){
          //std::cout << it->first << ": [";
          std::cout << it->first << ":";

          vector<int> ::iterator it2 ;
          for (it2=it->second.begin(); it2!=it->second.end(); ++it2){
               cout<<*it2<<" ";
           }
           //cout<<"]";
           cout<<endl;
    }*/
    std::cout << "knapsack_varContains:\n";
    for (it=knapsack_var.begin(); it!=knapsack_var.end(); ++it){
        //std::cout << it->first << ": [";
        std::cout << it->first << ": ";
        vector<int> ::iterator it2 ;
        for (it2=it->second.begin(); it2!=it->second.end(); ++it2){
            cout<<*it2<<" ";
        }
        //cout<<"]";
        cout<<endl;
    }
    cout<<"CheckFeasibleObject "<<endl;
    for (it=var_knapsack.begin(); it!=var_knapsack.end(); ++it){
        if(it->second.size()>1){
            cout<<it->first<<" Object in 2 knapsacks "<<endl;
        }else{
            cout<<it->first<<" Object Ok "<<endl;
        }
    }
    cout<<"CheckCapacity "<<endl;
    int y=0;
    for (it=knapsack_var.begin(); it!=knapsack_var.end(); ++it){
        double sum = 0.0;
        vector<int> ::iterator it2 ;
        for (it2=it->second.begin(); it2!=it->second.end(); ++it2){
            int ind = (*it2);
            sum = sum + double(w[ind]);
        }
        /********************/
        if(sum>C[it->first]){
            cout<<y<<" Exceeded Capacity "<<C[it->first]<<" "<<sum<<endl;
        }else{
            cout<<y<<" Capacity OK "<<C[it->first]<<" "<<sum<<endl;
        }
        /********************/
        y++;
    }
    cout<<"UsePercentage:"<<double(real_variables.size())/double(n) * 100<<endl;

    /*
    sort(real_variables.begin(),real_variables.end());

    cout<<"Get profit from instance"<<endl;
    double linear_profit = 0.0;
    double quadratick = 0.0;

    for(int i=0;i<real_variables.size();i++){
            linear_profit = linear_profit + p[real_variables[i]];
    }
     for(int i = 0;i<n;i++){
         for(int j = 0;j<n;j++){
             cout<<Q_2[i][j]<<" ";
         }
         cout<<endl;
     }
    for(int i=0;i<real_variables.size();i++){
            cout<<real_variables[i]<<" ";
    }
    cout<<endl;
     for(int i=0;i<real_variables.size();i++){
         for(int j=i+1;j<real_variables.size();j++){
             cout<< Q_2[real_variables[i]][real_variables[j]] <<" ";
             cout<< real_variables[i] <<" ";
             cout<< real_variables[j] <<endl;
                quadratick = quadratick + Q_2[real_variables[i]][real_variables[j]];
            }
        }
    cout<<linear_profit<<endl;
    cout<<quadratick<<endl;
    cout<<linear_profit+quadratick<<endl;
     */


}



#ifdef NO

int ij2k(int i, int j, int n)
/* This routine returns the k-th index
 * corresponding to element (i,j) i<j
 * in an uppert triangular matrix:
 * 	1 	2 	3 	4 	5 	...
 * 	x 	6 	7 	8 	9 	...
 * 	x 	x 	10 	11 	12 	...
 * ...
 * ...
 *  x	x	x			...
 * Note that the indexing starts from 0,
 * so for example element (1,2) is in fact (0,1).
 * This routine is used to index z_ij variables.
 * */
{
	if(i <= j)
		return i*n+j - i*(i+1)/2;
	else
		return j*n+i - j*(j+1)/2;
}

#endif

