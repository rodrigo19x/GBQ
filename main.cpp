//============================================================================
// Name        : QMKP.cpp
// Author      : Laura Galli
// Version     :
// Copyright   : Your copyright notice
// Description : algorithm for QMKP, C++, Ansi-style
//============================================================================

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

//Escribiendo un cambio !!!!

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "OPTUtils.h"
//#include "QMKP_01qpCPX.h"
#include "QMKP_01qpCPX2.h"
#include "localsearch.h"
#include "qkp_grdy.h"
//#include "QMKP_01qpGRB.h"
#include "utils.h"
#include <random>
#include <omp.h>

extern "C" int quadknap(int no, int cap, int *ptab, int *wtab, int *xtab, double T_L); // one way

#define NO_debug

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

//#define debug_PRINT

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- GLOBALS ------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------------------------*/
/*------------------------------- TYPES ------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void str2val( const char* const str , T &sthg ) {
    istringstream( str ) >> sthg;
}

int generateQMKP(int n, int *w, int **C, int m)
/*
 * This function creates a QMKP instance with m knapsacks starting from a QKP instance.
 * The (identical) knapsacks capacities are set to 80% of the sum
 * of the instances objects weights divided by the number of knapsacks.
 * C1 = C2 =...=Cm = 0.8 \sum_{i=1...n} w_i/m.
 *
 * */
{
    int cap = 0;

    (*C) = new int[m]; //knapsack capacities vector

    for(int i = 0; i < n; i++)
        cap += w[i];

    cap = 0.8*cap/m;

    for(int i = 0; i < m; i++)
        (*C)[i] = cap;

    return 0;
}

int readBillionnetQKP(char *filename, int *n, int **p, int **Q, int **w, int *c)
/*
 *  This function reads QKP instances according to the Billionnet format described in
 *
 *  http://cedric.cnam.fr/%7Esoutif/QKP/format.html
 *
 *  Each (QKP) instance corresponds to a file containing the following information:
 *  - the reference of the instance (r_10_100_13 in the following example)
 *  - the number of variables (n) (10 in the following example)
 *  - the linear coefficients (c_i) of the objective function (91 78 22 4 48 85 46 81 3 26)
 *  - the quadratic coefficients (c_ij) of the objective function
 *  - a blank line
 *  - 0 if the constraint is of type <= (i.e. always since we are considering (QKP) instances) and 1 if the constraint is an equality constraint
 *  - the capacity of the knapsack (145)
 *  - the coefficients of the capacity constraints (weights, a_i) (34 33 12 3 43 26 10 2 48 39 )
 *  - some comments
 *
 *  Here is an example concerning a 10-variable instance:
 *  r_10_100_13
 *  10
 *  91 78 22 4 48 85 46 81 3 26
 *  55 23 35 44 5 91 95 26 40
 *  92 11 20 43 71 83 27 65
 *  7 57 33 38 57 63 82
 *  100 87 91 83 44 48
 *  69 57 79 89 21
 *  9 40 22 26
 *  50 6 7
 *  71 52
 *  17
 *
 *  0
 *  145
 *  34 33 12 3 43 26 10 2 48 39
 *
 *  Comments
 *  Density : 100.00 %
 *  Seed : 13
 *
 * */
{
    string line;
    int num;

    // open problem file- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ifstream QKPfile(filename);
    if (!QKPfile.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        return (1);
    }

    // read instance- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //skip first line of the file with instance name
    getline(QKPfile, line);

    //read n (number of variables)
    QKPfile >> *n;
    if (QKPfile.fail()) {
        cerr << "Error reading from file " << filename << endl;
        return (1);
    }

    //allocate vectors/matrices
    *p = new int[(*n)]; //linear profit vector
    *Q = new int[(*n)*(*n-1)/2]; //quadratic profit matrix
    *w = new int[(*n)]; //weights of capacity constraint

    //read linear profits
    for (int i = 0; i < (*n); i++) {
        QKPfile >> (*p)[i];
        if (QKPfile.fail()) {
            cerr << "Error reading from file " << filename << endl;
            return (1);
        }
    }

    //read upper triangular profit (quadratic) matrix
    int k = 0;
    for (int i = 0; i < (*n-1); i++) {
        for(int j = 0; j < (*n -1 -i); j++) {
            QKPfile >> (*Q)[k++];
            if (QKPfile.fail()) {
                cerr << "Error reading from file " << filename << endl;
                return (1);
            }
        }
    }

    //skip empty line of the file
    getline(QKPfile, line);

    //skip 0 in the file
    QKPfile >> num;

    //read c (knapsack capacity)
    QKPfile >> *c;
    if (QKPfile.fail()) {
        cerr << "Error reading from file " << filename << endl;
        return (1);
    }

    //read weights
    for (int i = 0; i < (*n); i++) {
        QKPfile >> (*w)[i];
        if (QKPfile.fail()) {
            cerr << "Error reading from file " << filename << endl;
            return (1);
        }
    }

    return 0;
}


int printInstance(int n, int m, int *p, int *Q, int *w, int *C)
{

    cout << "n= " << n << endl;
    cout << "m= " << m << endl;

    cout << "linear profits p: " << endl;
    for(int i = 0; i < n; i++)
        cout << p[i] << " ";
    cout << endl;

    /*
    cout << "\nquadratic profit matrix Q: " << endl;
    int t = 0;
    for (int i = 0; i < n-1; i++) {
        for(int j = 0; j < n -1 -i; j++) {
            cout << Q[t++] << " ";
        }
        cout << endl;
    }
*/
    cout << "\nweights w: " << endl;
    for(int i = 0; i < n; i++)
        cout << w[i] << " ";
    cout << endl;

    cout << "\nCapacities C: " << endl;
    for(int i = 0; i < m; i++)
        cout << C[i] << " ";
    cout << endl;

    return 0;
}

int main(int argc, char **argv) {

    /* data for QKP instance
     * */
    int n; //number of objects
    int *Q; //matrix for quadratic profit term
    int *p; //array for linear profit term
    int *w; //array of weights for capacity constraint
    int c; //capacity (for single knapsack)
    int *C; //array of capacities (for multiple knapsacks)

    /* input data from command line
     * */
    char *fin;//instance file name
    char *fout;//output file name
    char *flog;//log file name
    int alg;//solution algorithm to use
    int m; //number of knapsacks
    int TL; //time limit in seconds
    bool intflag = true;

    /* output data
     * */
    double objval = 0;
    int status = 0;
    OPTtimers *timer = new OPTtimers();  ///< timer
    double duration;

    /*others*/
    int k;

    /* how to call this little program:
     * ./main <fin> <m> <alg> TL <fout> <flog>
     *
     * <fin> name of the instance file to read
     *
     * <m> number of knapsack (m=1 the problem is QKP instead of QMKP)
     *
     * <alg> type of solution algorithm to use:
     * 1 = use CPLEX 01 QP solver
     * 2 = use quadknap program of David Pisinger (it only works for QKP, i.e., m=1)
     *
     * <fout> output file
     *
     * <flog> log file
     * */


    /* read param values from command line
     * */
    fin = argv[1];
    str2val(argv[2], m);
    str2val(argv[3], alg);
    str2val(argv[4], TL);
    fout = argv[5];
    flog = argv[6];

    /* if alg quadknap is chosen,
     * only QKP instances can be solved
     * */
    if(alg == 2 & m > 1)
    {
        cout << "Error: quadknapk cannot solve multiple knapsack!" << endl;
        exit(1);
    }


    /* read Billionnet & Soutif QKP instance
     * NB. These are QKP instances,
     * i.e., with a single knapsack
     * */
    readBillionnetQKP(fin, &n, &p, &Q, &w, &c);

    /* generate QMKP instance with m knapsacks,
     * the (identical) capacities are set to 80% of the sum
     * of the instances objects' weights divided by the number of knapsacks
     * as in Hiley & Julstrom paper.
     * */
    if (m == 1) { C = new int[1]; C[0] = c;}
    else generateQMKP(n, w, &C, m);


#ifdef NO_debug
    /* Print QMKP instance */
    printInstance(n, m, p, Q, w, C);

    cout << "\nQuadratic: " << endl;

    for (int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i!=j){
                k = ij2k(n,i,j);
                cout << Q[k] << " ";
            }else{
                cout << "- ";
            }
        }
        cout << endl;
    }

#endif



    /* switch solution algorithm
     * 1 = CPLEX 01 solver (for QMKP, which includes QKP as a special case with m=1)
     * 2 = David Pisinger "quadknap" algorithm (only for QKP)
     * 3 = GRB
     * */
    status = 0;

    int *x = new int[n];//solution vector for "quadknap"
    int *xQ = new int[n*n];//profit matrix for "quadknap"
    vector<vector<int> > x_con;

    timer->ReSet();
    timer->Start();

    //status = solveQMKP_01qp_CPX2(n, m, p, Q, w, C, x_con, &objval, TL, intflag, flog);
    //cout<<"Greedy one"<<endl;
    //greedy_one(n, m, p, Q, w, C);
    //greedy_one(n, m, p, Q, w, C, x_con);
    double start;
    double end;
    double end2;
    double end3;
    start = omp_get_wtime();
    qkp_grdy(n, m, p, Q, w, C, x_con);

    vector<vector<int> > x_con2;
    for(int i=0;i<x_con.size();i++){
        vector<int> aux(x_con[i]);
        x_con2.push_back(aux);
    }
    end = omp_get_wtime();
    printf("GreedyTime %f\n", end - start);
    /*
    exchange_ls2A(n, m, p, Q, w, C, x_con);
    end2 = omp_get_wtime();
    printf("TotalTimeA %f\n", end2 - start);
    exchange_ls2B(n, m, p, Q, w, C, x_con2);
    end3 = omp_get_wtime();
    printf("TotalTimeB %f\n", (end3 - start)-(end2 - end));
    */
    timer->Stop();
    duration = timer->Read();


    /* print output
     * */
    cout << "\n" << fin << "\t" << n << "\t" << m << "\t" << alg << "\t" << status << "\t" << objval << "\t" << duration << endl;

    ofstream outfile;
    outfile.open(fout, ios::app);

    outfile << "\n" << fin << "\t" << n << "\t" << m << "\t" << alg << "\t" << status << "\t" << objval << "\t" << duration << endl;

    outfile.close();

    /* free memory
     */
    delete [] x;
    delete [] xQ;

    delete [] Q;
    delete [] C;
    delete [] p;
    delete [] w;

    return 0;
}

