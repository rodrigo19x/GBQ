//
// Created by carlosrey on 28/09/20.
//


#include "QMKP_01qpCPX2.h"
#include <random>



#define Link1
#define Link2
//#define Link3
#define Qcap
//#define Qcard
//#define Qcap2

int solveQMKP_01qp_CPX2( int n, int m, int *p, int *Q, int *w, int *C, vector<vector<double> > x_con,double * objval, int TL, bool intflag,  char *logfile){

    IloEnv env;
    IloModel model(env);
    //GRBQuadExpr obj = 0;
    IloNumVarArray x(env, n*m);
    int size_y = (n*(n-1));
    IloNumVarArray y(env, m*size_y);

    int c=0;
    for(int k=0;k<m;k++){
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i==j)
                    continue;
                string s = "y_" + to_string(i)+"_"+to_string(j)+"_"+to_string(k);
                y[c]  = IloNumVar(env, 0, 1, IloNumVar::Float, s.c_str());
                c++;

            }
        }
    }





    /*

   IloArray<IloArray< IloNumVarArray> > y(env,m);
            for(int i=0;i<m;i++){
      y[i] = IloArray(env, m);
      for(int k=0;k<m;k++){
        y[i][k] = IloNumVarArray(env, n);
        for(int j=0;j<n;j++){
              string s = "y" + to_string();
             y[i][k][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, s.c_str());
        }
      }
    }*/


    /***************************************************************/
    /******************Objetive Function****************************/
    /***************************************************************/
    //Normal profit


    c=0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            string s = "x_" + to_string(i)+"_"+to_string(j);
            x[c]  = IloNumVar(env, 0, 1, IloNumVar::Float, s.c_str());
            c++;
        }
    }

    IloExpr expr(env);
    for(auto i = 0u; i < n*m; ++i) {
        //string s = "x" + to_string(i+1);
        //x[i]  = IloNumVar(env, 0, 1, IloNumVar::Bool, s.c_str());
        int k = i % n;
        expr = expr + (p[k]) * x[i];
    }

    //Quadratc profit
    int count_var=0;
    //cout<<"---"<<endl;
    //cout<<n<<endl;
    //cout<<"---"<<endl;
    int index = 0;
    for(auto k=0;k<m;k++){
        for(auto i = 0; i < n; i++) {
            for(auto j = 0; j < n; j++) {
                //expr = expr + (Q[ij2k(n, i-count_var,j - count_var)])*(x[i] * x[j]);//FIXME: double profit
                if(i==j)
                    continue;

                //int index= k*size_y+(i*(n)+j);

                //cout<<k<<" "<<i<<" "<<j<<" "<<index<<endl;
                //double value1 = Q[ij2k(n, i-count_var,j - count_var)];
                //I value2 = y[index];
                double aux = (double)(Q[ij2k(n,i,j)])/2.0;
                expr = expr + ( aux )*( y[index]);//FIXME: double profit
                index++;

            }
        }
        count_var=count_var+n;
    }

    IloObjective obj(env, expr, IloObjective::Maximize);

    // Add the objective function to the model
    model.add(obj);
    expr.clear();


    /***************************************************************/
    /******************F Constraints********************************/
    /***************************************************************/
    //1. Weight Constrained per knapsack
    int cont_constr = 0;
    IloRangeArray constraint_1array(env, m);

    for(auto i = 0u; i < m; ++i) {

        for(auto j = cont_constr; j < n+cont_constr; ++j) {
            expr = expr + (x[j] * w[j-cont_constr]);
        }

        cont_constr = cont_constr +n;
        string s = "c" + to_string(i+1);
        constraint_1array[i] = IloRange(env, -IloInfinity, expr, C[i], s.c_str());
        expr.clear();

    }
    model.add(constraint_1array);


    //2) One object per knapsack
    auto cont = 0u;
    IloRangeArray constraint_2array(env, n);
    for(auto i = 0u; i < n; ++i) {
        //GRBLinExpr expr = 0;
        auto j = cont;
        while(j < (n*m)) {
            expr = expr + x[j];
            j=j+n;
        }
        string s = "c" + to_string((i+m)+1);
        //model.addConstr(expr, GRB_LESS_EQUAL, 1.0, s);
        constraint_2array[i] = IloRange(env, -IloInfinity, expr, 1, s.c_str());
        expr.clear();
        cont++;
    }
    model.add(constraint_2array);

    /***************************************************************/
    /******************Linking Constraints**************************/
    /***************************************************************/
#ifdef Link1

    IloRangeArray Linking_1(env, m*size_y);

    expr.clear();
    c = 0;
    for(int k=0;k<m;k++){
        for(int i=0;i<n;i++){
            for(int j=0; j<n;j++){
                if(i==j)
                    continue;

                int index_y = k*size_y+ij2k2(n,i,j);
                int index_x = k*n+i;
                expr = y[index_y] - x[index_x];

                string s = "l1_" + to_string(i)+"_"+to_string(j)+"_"+to_string(k);

                Linking_1[c] = IloRange(env, -IloInfinity, expr, 0, s.c_str());
                c++;
                expr.clear();
            }
        }
    }
    model.add(Linking_1);
#endif

    //cout<<"L1 complete"<<endl;
#ifdef Link3
    IloRangeArray Linking_3(env, m*size_y);
          expr.clear();
          c = 0;
          for(int k=0;k<m;k++){
            for(int i=0;i<n;i++){
              for(int j=0;j<n;j++){
              if(i==j)
                  continue;

                int index_y = k*size_y+ij2k2(n,i,j);
                int index_x_1 = k*n+i;
                int index_x_2 = k*n+j;
                expr = y[index_y] - x[index_x_1] - x[index_x_2];
                string s = "l3_" + to_string(i)+"_"+to_string(j)+"_"+to_string(k);
                Linking_3[c] = IloRange(env, -1.0, expr, IloInfinity, s.c_str());
                c++;
                expr.clear();
              }
            }
          }
          model.add(Linking_3);
#endif

#ifdef Link2
    IloRangeArray Linking_2(env, (m*size_y)/2);
    expr.clear();
    c = 0;
    for(int k=0;k<m;k++){
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                int index_y_1 = k*size_y+ij2k2(n,i,j);
                int index_y_2 = k*size_y+ij2k2(n,j,i);

                //cout<<k<<" "<<i<<" "<<j<<" "<<index_y_1<<" "<<index_y_2<<endl;
                expr = y[index_y_1] - y[index_y_2];
                string s = "l2_" + to_string(i)+"_"+to_string(j)+"_"+to_string(k);
                //Linking_2[c] = IloRange(env, -IloInfinity, expr, 0, s.c_str());
                Linking_2[c] = IloRange(env, 0, expr, 0, s.c_str());
                c++;
                expr.clear();
            }
        }
    }
    model.add(Linking_2);
#endif

    /***************************************************************/
    /******************Q Constraints********************************/
    /***************************************************************/

#ifdef Qcap
    //Q Capacity
    IloRangeArray Q_Capacity(env, (n*m));
    expr.clear();
    c = 0;
    for(int i=0;i<n;i++){
        for(int k=0;k<m;k++){


            for (int j = 0; j < n; j++)
            {
                if(i==j)
                    continue;

                int index_y = k*size_y+ij2k2(n,i,j);
                expr = expr + (w[j] * y[index_y]);
            }


            int index_x = k*n+i;

            expr = expr + ((w[i] - C[k])*x[index_x]);

            string s = "qc1_" + to_string(i)+"_"+to_string(k);
            //Linking_2[c] = IloRange(env, -IloInfinity, expr, 0, s.c_str());
            Q_Capacity[c] = IloRange(env, -IloInfinity, expr, 0, s.c_str());
            c++;
            expr.clear();

        }
    }
    model.add(Q_Capacity);
#endif
    //Q Capacity
#ifdef Qcard
    IloRangeArray Q_Cardinality(env, (n*(n-1)*m));
          expr.clear();
          c = 0;
          for(int i=0;i<n;i++){
            for (int j = 0; j < n; j++){
              if(i==j)
                    continue;
                for (int h = 0; h < m; h++){

                    for (int k = 0; k < m; k++){
                      if(h==k)
                        continue;
                      int index_y = k*size_y+ij2k2(n,i,j);
                      expr = expr + (y[index_y]);
                    }

                    int index_x = h*n+i;
                    expr = expr + x[index_x];

                  string s = "qc2_" + to_string(i)+"_"+to_string(j)+"_"+to_string(h);
                  //Linking_2[c] = IloRange(env, -IloInfinity, expr, 0, s.c_str());
                  Q_Cardinality[c] = IloRange(env, -IloInfinity, expr, 1, s.c_str());
                  c++;
                  expr.clear();
                }
            }

          }

          model.add(Q_Cardinality);
#endif

    /***************************************************************/
    /******************Q Constraints2********************************/
    /***************************************************************/

#ifdef Qcap2
    //Q Capacity
          IloRangeArray Q_Capacity2(env, (n*m));
          expr.clear();
          c = 0;
          for(int i=0;i<n;i++){
            for(int k=0;k<m;k++){
                for (int j = 0; j < n; j++)
                {
                  if(i==j)
                    continue;

                  int index_y = k*size_y+ij2k2(n,i,j);
                  expr = expr - (w[j] * y[index_y]);
                  int index_x = k*n+j;
                  expr = expr + (w[j]*x[index_x]);

                }


                int index_x = k*n+i;

                expr = expr + (C[k]*x[index_x]);

                string s = "qc3_" + to_string(i)+"_"+to_string(k);
                //Linking_2[c] = IloRange(env, -IloInfinity, expr, 0, s.c_str());
                Q_Capacity2[c] = IloRange(env, -IloInfinity, expr, C[k], s.c_str());
                c++;
                expr.clear();

            }
          }
          model.add(Q_Capacity2);
#endif
    /***************************************************************/
    /***************************************************************/
    /***************************************************************/
    expr.end();
    // Create the solver object
    IloCplex cplex(model);


    cplex.setParam(IloCplex::Param::Read::DataCheck,1);
    //cplex.setParam(IloCplex::Param::OptimalityTarget,CPX_OPTIMALITYTARGET_OPTIMALGLOBAL);
    cplex.setParam(IloCplex::Param::Threads, 1);
    // cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, 1);
    cplex.setOut(env.getNullStream()) ;
    if(TL > 0){
        cplex.setParam(IloCplex::Param::TimeLimit, TL);
    }

    // Export model to file (useful for debugging!)

    int i=index;
    stringstream ss;
    ss << i;
    //string lp = "model_2.lp";
    //cplex.exportModel(lp.c_str());
    cout<<"Generating.."<<endl;
    cplex.exportModel("model_2.lp");
    //exit(0);



    bool solved = false;
    clock_t before,end;

    try {
        // Try to solve with CPLEX (and hope it does not raise an exception!)
        before = clock();
        solved = cplex.solve();
        end = clock();
    } catch(const IloException& e) {
        std::cerr << "\n\nCPLEX Raised an exception:\n";
        std::cerr << e << "\n";
        env.end();
        throw;
    }

    if(solved) {
        // If CPLEX successfully solved the model, print the results
        std::cout << "\n\nCplex success!\n";
        std::cout << "\tStatus: " << cplex.getStatus() << "\n";
        std::cout << "Objective value: " << cplex.getObjValue()<<"\n";
        std::cout << "Get Nodes  : " << cplex.getNnodes()<<"\n";
        std::cout << "Get Bound  : " << cplex.getBestObjValue()<< "\n";
        std::cout << "Get Cutoff : " << cplex.getCutoff()<<"\n";
        std::cout << "Get Gap    : " << cplex.getMIPRelativeGap()<<"\n";
        std::cout << "Total Time : " << float(end-before)/CLOCKS_PER_SEC <<"\n";
        double objval_aux = cplex.getObjValue();
        *objval = (objval_aux);
        vector<unsigned> variables;
        cout<<"Variables : "<<endl;
        int cont_print = 0;
        for(int i=0;i<n*m;i++){

            cout<<(cplex.getValue(x[i]))<<" ";

            if(cont_print==n){
                cout<<endl;
                cont_print=0;
            }


        }
        cout<<endl;
        /*
        cout<<"Profit: ";
        for(int i=0;i<n*m;i++){
            int k = i % n;
            cout<<p[k]<<" ";
        }
        cout<<endl;
        cout<<"Weight: ";
        for(int i=0;i<n*m;i++){
            int k = i % n;
            cout<<w[k]<<" ";
        }
        cout<<endl;

        report(n,m, p,Q, w, C,variables);
        */
    } else {
        std::cerr << "\n\nCplex error!\n";
        std::cerr << "\tStatus: " << cplex.getStatus() << "\n";
        std::cerr << "\tSolver status: " << cplex.getCplexStatus() << "\n";
    }

    env.end();
    return 1;

}



