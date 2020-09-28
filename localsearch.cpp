//
// Created by carlosrey on 28/09/20.
//

#include "localsearch.h"
#include "Solution.h"
#include <random>
#include <omp.h>

int exchange_ls(int n, int m, int * p, int * Q, int * w, int * C, vector < vector < int > > & x) {

    //Solution * S_var2 = new Solution(x);
    //S_var2->show();
    //S_var2->check(n ,m,w ,C,p,Q);

    vector < double > R(m, 0.0);
    vector < double > K(n, -1.0);

    //vector<int> selected;
    for (int k = 0; k < m; k++) {
        double sum_w = 0.0;
        for (int i = 0; i < n; ++i) {
            if (x[k][i] == 1) {
                K[i] = k;
                //selected.push_back(i);
            }
            sum_w = sum_w + (x[k][i] * w[i]);
        }
        R[k] = C[k] - sum_w;
    }

    /*
    cout<<"R"<<endl;
    for(int k=0;k<m;k++) cout<<R[k]<<" ";
        cout<<endl;
    cout<<"K"<<endl;
    for(int i=0;i<n;i++) cout<<K[i]<<" ";
    cout<<endl;
    //for (std::vector<int>::iterator it = selected.begin() ; it != selected.end(); it++){
    getchar();
    */
    for (int i = 0; i < n - 1; i++) {
        //int i = (*it);
        if (K[i] == -1) {
            continue;
        }
        //comment: compute P (i), the pairwise profit of item i;
        int k_prima = K[i];
        vector < double > P_sup(n, 0.0);

        for (int j = 0; j < n; j++) {
            if (K[j] == k_prima && i != j) P_sup[i] = P_sup[i] + Q[ij2k(n, i, j)];
        }

        int best_l = -1;
        int best_k_hat = -1;
        double current_best_value = -1;

        for (int l = i + 1; l < n; l++) {

            //cout<<"Checking "<<i<<" "<<l<<endl;
            if (K[l] == -1.0 || K[l] == K[i]) {
                continue;
            } else if ((R[K[l]] + w[l] - w[i]) < 0 || (R[k_prima] + w[i] - w[l]) < 0) {
                continue;
            } else {
                //comment: the exchange is feasible, compute the pairwise profit of item ;
                //and the pairwise profits produced by the exchange (P new (i) and P new ());
                int k_hat = K[l];

                vector < double > P(n, 0.0);
                copy(P_sup.begin(), P_sup.end(), P.begin());
                /*
                cout<<"P"<<endl;
                for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                cout<<endl;
                */
                vector < double > P_new(n, 0.0);
                P[l] = 0.0;
                //cout<<"cambiando el objeto "<<l<<" de la mochila "<<k_hat<<endl;

                for (int j = 0; j < n; j++) {
                    if (K[j] == k_hat) {
                        if (l != j) {
                            P[l] = P[l] + Q[ij2k(n, l, j)];
                            //cout<<"Sumando en P "<<Q[ij2k(n,l,j)]<<endl;
                        }
                        if (i != j && j != l) {
                            //if(i!=j){
                            P_new[i] = P_new[i] + Q[ij2k(n, i, j)];
                            //cout<<"Sumando en PNew"<<Q[ij2k(n,i,j)]<<endl;
                        }
                    }
                }
                /*
                cout<<"P"<<endl;
                for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                cout<<endl;
                cout<<"P_new"<<endl;
                for (int i = 0; i < P.size() ; ++i) cout<<P_new[i]<<" ";
                cout<<endl;
                cout<<"--"<<endl;
                */
                for (int j = 0; j < n; j++) {
                    if (K[j] == k_prima && l != j && i != j) {
                        //if(K[j]==k_prima && l!=j){
                        P_new[l] = P_new[l] + Q[ij2k(n, l, j)];
                        //cout<<"Sumando en PNew"<<Q[ij2k(n,l,j)]<<endl;
                    }
                }
                /*
                cout<<"P"<<endl;
                for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                cout<<endl;
                cout<<"P_new"<<endl;
                for (int i = 0; i < P.size() ; ++i) cout<<P_new[i]<<" ";
                cout<<endl;
                */
                //double gain = (P_new[i] + P_new[l]) - (P[i] + P[l]);
                if (P_new[i] + P_new[l] > P[i] + P[l]) {
                    //if(gain>current_best_value && gain > 0 ){

                    //cout<<"-------------------"<<endl;

                    //cout<<"Improve :D"<<endl;
                    //cout<<k_prima<<" "<<k_hat<<"|"<<i<<" "<<l<<"->"<< P_new[i] + P_new[l]<<" | "<< P[i] + P[l] <<endl;
                    x[k_prima][i] = 0;
                    x[k_hat][i] = 1;
                    x[k_hat][l] = 0;
                    x[k_prima][l] = 1;
                    K[i] = k_hat;
                    K[l] = k_prima;

                    for (int k_up = 0; k_up < m; k_up++) {
                        double sum_w = 0.0;
                        for (int i_up = 0; i_up < n; i_up++) {
                            sum_w = sum_w + (x[k_up][i_up] * w[i_up]);
                        }
                        R[k_up] = C[k_up] - sum_w;
                    }
                    Solution * S_var = new Solution(x);
                    S_var -> show();
                    /*
                    cout<<"R"<<endl;
                    for(int k=0;k<m;k++) cout<<R[k]<<" ";
                        cout<<endl;
                    cout<<"K"<<endl;
                    for(int i=0;i<n;i++) cout<<K[i]<<" ";
                    cout<<endl;
                    */
                    S_var -> check(n, m, w, C, p, Q);
                    //cout<<"-------------------"<<endl;
                    break;
                    /*
                    cout<<"Improve "<<l<<" "<<k_hat<<" "<<gain<<endl;
                    best_l = l;
                    best_k_hat = k_hat;
                    current_best_value = gain;
                    */
                }
            }
        }
        /*
        cout<<"Final best "<<best_l<<" "<<best_k_hat<<" "<<current_best_value<<" item: "<<i<<endl;
        if(current_best_value!=-1){
                    cout<<"-------------------"<<endl;
                    cout<<"Improve :D"<<endl;
                    x[k_prima][i] = 0;
                    x[best_k_hat][i] = 1;
                    x[best_k_hat][best_l] = 0;
                    x[k_prima][best_l] = 1;
                    K[i]=best_k_hat;
                    K[best_l]=k_prima;

                    for(int k_up=0;k_up<m;k_up++){
                        double sum_w = 0.0;
                        for (int i_up = 0; i_up < n; i_up++){
                            sum_w = sum_w + (x[k_up][i_up]*w[i_up]);
                        }
                        R[k_up] = C[k_up] - sum_w;
                    }
                    Solution * S_var = new Solution(x);
                    S_var->show();
                    cout<<"R"<<endl;
                    for(int k=0;k<m;k++) cout<<R[k]<<" ";
                        cout<<endl;
                    cout<<"K"<<endl;
                    for(int i=0;i<n;i++) cout<<K[i]<<" ";
                    cout<<endl;
                    S_var->check(n ,m,w ,C,p,Q);
                    cout<<"-------------------"<<endl;

        }*/
    }

    Solution * S_var = new Solution(x);
    S_var -> show();
    S_var -> check(n, m, w, C, p, Q);

    return 0;
}

int exchange_ls2B(int n, int m, int * p, int * Q, int * w, int * C, vector < vector < int > > & x) {

    //Solution * S_var2 = new Solution(x);
    //S_var2->show();
    //S_var2->check(n ,m,w ,C,p,Q);

    vector < double > R(m, 0.0);
    vector < double > K(n, -1.0);

    //vector<int> selected;
    for (int k = 0; k < m; k++) {
        double sum_w = 0.0;
        for (int i = 0; i < n; ++i) {
            if (x[k][i] == 1) {
                K[i] = k;
                //selected.push_back(i);
            }
            sum_w = sum_w + (x[k][i] * w[i]);
        }
        R[k] = C[k] - sum_w;
    }


    bool flag = true;
    while (flag) {
        flag = false;

        int best_l = -1;
        int best_k_hat = -1;
        double current_best_value = -1;
        int k_prima = -1;
        int i = 0;
        int i_aux = -1;
        int k_prima_aux = -1;
        while (i < n - 1) {
            //int i = (*it);
            //cout<<"i "<<i<<endl;
            if (K[i] == -1) {
                i++;
                continue;
            }
            //comment: compute P (i), the pairwise profit of item i;
            k_prima = K[i];
            vector < double > P_sup(n, 0.0);

            for (int j = 0; j < n; j++) {
                if (K[j] == k_prima && i != j) P_sup[i] = P_sup[i] + Q[ij2k(n, i, j)];
            }




            /*****************************************************/
            /************inner loop*******************************/
            /*****************************************************/
            for (int l = i + 1; l < n; l++) {
                //cout<<"L "<<l<<endl;
                //cout<<"Checking "<<i<<" "<<l<<endl;
                if (K[l] == -1.0 || K[l] == K[i]) {
                    continue;
                } else if ((R[K[l]] + w[l] - w[i]) < 0 || (R[k_prima] + w[i] - w[l]) < 0) {
                    continue;
                } else {
                    //comment: the exchange is feasible, compute the pairwise profit of item ;
                    //and the pairwise profits produced by the exchange (P new (i) and P new ());
                    int k_hat = K[l];

                    vector < double > P(n, 0.0);
                    copy(P_sup.begin(), P_sup.end(), P.begin());
                    /*
                    cout<<"P"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                    cout<<endl;
                    */
                    vector < double > P_new(n, 0.0);
                    P[l] = 0.0;
                    //cout<<"cambiando el objeto "<<l<<" de la mochila "<<k_hat<<endl;

                    for (int j = 0; j < n; j++) {
                        if (K[j] == k_hat) {
                            if (l != j) {
                                P[l] = P[l] + Q[ij2k(n, l, j)];
                                //cout<<"Sumando en P "<<Q[ij2k(n,l,j)]<<endl;
                            }
                            if (i != j && j != l) {
                                //if(i!=j){
                                P_new[i] = P_new[i] + Q[ij2k(n, i, j)];
                                //cout<<"Sumando en PNew"<<Q[ij2k(n,i,j)]<<endl;
                            }
                        }
                    }
                    /*
                    cout<<"P"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                    cout<<endl;
                    cout<<"P_new"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P_new[i]<<" ";
                    cout<<endl;
                    cout<<"--"<<endl;
                    */
                    for (int j = 0; j < n; j++) {
                        if (K[j] == k_prima && l != j && i != j) {
                            //if(K[j]==k_prima && l!=j){
                            P_new[l] = P_new[l] + Q[ij2k(n, l, j)];
                            //cout<<"Sumando en PNew"<<Q[ij2k(n,l,j)]<<endl;
                        }
                    }
                    /*
                    cout<<"P"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                    cout<<endl;
                    cout<<"P_new"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P_new[i]<<" ";
                    cout<<endl;
                    */
                    double gain = (P_new[i] + P_new[l]) - (P[i] + P[l]);
                    //if(P_new[i] + P_new[l] > P[i] + P[l] ){
                    if (gain > 0) cout<<"Posible improve "<<i<<" "<<l<<endl;
                    if (gain > current_best_value && gain > 0) {

                        //cout<<"-------------------"<<endl;

                        //cout<<"Improve :D"<<endl;
                        //cout<<k_prima<<" "<<k_hat<<"|"<<i<<" "<<l<<"->"<< P_new[i] + P_new[l]<<" | "<< P[i] + P[l] <<endl;
                        /*
                        x[k_prima][i] = 0;
                        x[k_hat][i] = 1;
                        x[k_hat][l] = 0;
                        x[k_prima][l] = 1;
                        K[i]=k_hat;
                        K[l]=k_prima;

                        for(int k_up=0;k_up<m;k_up++){
                            double sum_w = 0.0;
                            for (int i_up = 0; i_up < n; i_up++){
                                sum_w = sum_w + (x[k_up][i_up]*w[i_up]);
                            }
                            R[k_up] = C[k_up] - sum_w;
                        }
                        Solution * S_var = new Solution(x);
                        S_var->show();

                        cout<<"R"<<endl;
                        for(int k=0;k<m;k++) cout<<R[k]<<" ";
                            cout<<endl;
                        cout<<"K"<<endl;
                        for(int i=0;i<n;i++) cout<<K[i]<<" ";
                        cout<<endl;

                        S_var->check(n ,m,w ,C,p,Q);
                        //cout<<"-------------------"<<endl;
                        break;
                        */

                        cout << "Improve " << l << " " << k_hat << " " << gain << endl;
                        best_l = l;
                        best_k_hat = k_hat;
                        current_best_value = gain;
                        i_aux = i;
                        k_prima_aux = k_prima;
                        flag = true;
                    }
                }
            }
            /*****************************************************/
            /*****************************************************/
            /*****************************************************/
            //cout << "Final best " << best_l << " " << best_k_hat << " " << current_best_value << " item: " << i << endl;
            //cout << i << " " << endl;
            i++;

        }
        if (current_best_value != -1) {
            cout << "-------------------" << endl;
            cout << "Improve :D" << endl;
            x[k_prima_aux][i_aux] = 0;
            x[best_k_hat][i_aux] = 1;
            x[best_k_hat][best_l] = 0;
            x[k_prima_aux][best_l] = 1;
            K[i_aux] = best_k_hat;
            K[best_l] = k_prima_aux;

            for (int k_up = 0; k_up < m; k_up++) {
                double sum_w = 0.0;
                for (int i_up = 0; i_up < n; i_up++) {
                    sum_w = sum_w + (x[k_up][i_up] * w[i_up]);
                }
                R[k_up] = C[k_up] - sum_w;
            }
            Solution * S_var = new Solution(x);
            S_var -> show();
            cout << "R" << endl;
            for (int k = 0; k < m; k++) cout << R[k] << " ";
            cout << endl;
            cout << "K" << endl;
            for (int i = 0; i < n; i++) cout << K[i] << " ";
            cout << endl;
            S_var -> check(n, m, w, C, p, Q);
            cout << "-------------------" << endl;
        }

    }

    Solution * S_var = new Solution(x);
    S_var -> show();
    double f_value = S_var -> check(n, m, w, C, p, Q);
    cout<<"FinalValueB "<<f_value<<endl;

    return 0;
}


int exchange_ls2A(int n, int m, int * p, int * Q, int * w, int * C, vector < vector < int > > & x) {

    //Solution * S_var2 = new Solution(x);
    //S_var2->show();
    //S_var2->check(n ,m,w ,C,p,Q);

    vector < double > R(m, 0.0);
    vector < double > K(n, -1.0);

    //vector<int> selected;
    for (int k = 0; k < m; k++) {
        double sum_w = 0.0;
        for (int i = 0; i < n; ++i) {
            if (x[k][i] == 1) {
                K[i] = k;
                //selected.push_back(i);
            }
            sum_w = sum_w + (x[k][i] * w[i]);
        }
        R[k] = C[k] - sum_w;
    }


    int i = 0;
    bool flag = true;
    while (flag) {
        flag = false;
        while (i < n - 1) {
            //int i = (*it);
            //cout<<"i "<<i<<endl;
            if (K[i] == -1) {
                i++;
                continue;
            }
            //comment: compute P (i), the pairwise profit of item i;
            int k_prima = K[i];
            vector < double > P_sup(n, 0.0);

            for (int j = 0; j < n; j++) {
                if (K[j] == k_prima && i != j) P_sup[i] = P_sup[i] + Q[ij2k(n, i, j)];
            }

            int best_l = -1;
            int best_k_hat = -1;
            double current_best_value = -1;


            /*****************************************************/
            /************inner loop*******************************/
            /*****************************************************/
            for (int l = i + 1; l < n; l++) {
                //cout<<"L "<<l<<endl;
                //cout<<"Checking "<<i<<" "<<l<<endl;
                if (K[l] == -1.0 || K[l] == K[i]) {
                    continue;
                } else if ((R[K[l]] + w[l] - w[i]) < 0 || (R[k_prima] + w[i] - w[l]) < 0) {
                    continue;
                } else {
                    //comment: the exchange is feasible, compute the pairwise profit of item ;
                    //and the pairwise profits produced by the exchange (P new (i) and P new ());
                    int k_hat = K[l];

                    vector < double > P(n, 0.0);
                    copy(P_sup.begin(), P_sup.end(), P.begin());
                    /*
                    cout<<"P"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                    cout<<endl;
                    */
                    vector < double > P_new(n, 0.0);
                    P[l] = 0.0;
                    //cout<<"cambiando el objeto "<<l<<" de la mochila "<<k_hat<<endl;

                    for (int j = 0; j < n; j++) {
                        if (K[j] == k_hat) {
                            if (l != j) {
                                P[l] = P[l] + Q[ij2k(n, l, j)];
                                //cout<<"Sumando en P "<<Q[ij2k(n,l,j)]<<endl;
                            }
                            if (i != j && j != l) {
                                //if(i!=j){
                                P_new[i] = P_new[i] + Q[ij2k(n, i, j)];
                                //cout<<"Sumando en PNew"<<Q[ij2k(n,i,j)]<<endl;
                            }
                        }
                    }
                    /*
                    cout<<"P"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                    cout<<endl;
                    cout<<"P_new"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P_new[i]<<" ";
                    cout<<endl;
                    cout<<"--"<<endl;
                    */
                    for (int j = 0; j < n; j++) {
                        if (K[j] == k_prima && l != j && i != j) {
                            //if(K[j]==k_prima && l!=j){
                            P_new[l] = P_new[l] + Q[ij2k(n, l, j)];
                            //cout<<"Sumando en PNew"<<Q[ij2k(n,l,j)]<<endl;
                        }
                    }
                    /*
                    cout<<"P"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P[i]<<" ";
                    cout<<endl;
                    cout<<"P_new"<<endl;
                    for (int i = 0; i < P.size() ; ++i) cout<<P_new[i]<<" ";
                    cout<<endl;
                    */
                    double gain = (P_new[i] + P_new[l]) - (P[i] + P[l]);
                    //if(P_new[i] + P_new[l] > P[i] + P[l] ){
                    if (gain > current_best_value && gain > 0) {

                        cout << "Improve " << l << " " << k_hat << " " << gain << endl;
                        best_l = l;
                        best_k_hat = k_hat;
                        current_best_value = gain;
                        flag = true;

                    }
                }
            }
            /*****************************************************/
            /*****************************************************/
            /*****************************************************/
            cout << "Final best " << best_l << " " << best_k_hat << " " << current_best_value << " item: " << i << endl;
            if (current_best_value != -1) {
                cout << "-------------------" << endl;
                cout << "Improve :D" << endl;
                x[k_prima][i] = 0;
                x[best_k_hat][i] = 1;
                x[best_k_hat][best_l] = 0;
                x[k_prima][best_l] = 1;
                K[i] = best_k_hat;
                K[best_l] = k_prima;

                for (int k_up = 0; k_up < m; k_up++) {
                    double sum_w = 0.0;
                    for (int i_up = 0; i_up < n; i_up++) {
                        sum_w = sum_w + (x[k_up][i_up] * w[i_up]);
                    }
                    R[k_up] = C[k_up] - sum_w;
                }
                Solution * S_var = new Solution(x);
                S_var -> show();
                cout << "R" << endl;
                for (int k = 0; k < m; k++) cout << R[k] << " ";
                cout << endl;
                cout << "K" << endl;
                for (int i = 0; i < n; i++) cout << K[i] << " ";
                cout << endl;
                S_var -> check(n, m, w, C, p, Q);
                cout << "-------------------" << endl;
            }
            cout << i << " " << endl;
            i++;
        }
    }

    Solution * S_var = new Solution(x);
    S_var -> show();
    double f_value = S_var -> check(n, m, w, C, p, Q);
    cout<<"FinalValueA "<<f_value<<endl;
    return 0;
}