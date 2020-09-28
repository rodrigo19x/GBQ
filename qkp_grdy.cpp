//
// Created by carlosrey on 28/09/20.
//

#include "qkp_grdy.h"

#include "qkp_grdy.h"
#include "Solution.h"
#include "QMKP_01qpCPX2.h"
#include <random>
#include <omp.h>
#include "gurobi_c++.h"


extern "C" int quadknap(int no, int cap, int **ptab, int *wtab, int *xtab,double T_L); // one way

int ij2k_main(int n, int i, int j) {



    if (i == j) {
        cout << "\nerror: ij2k wrong index i==j";
        exit(1);
    }

    int k;

    if (j < i) k = (n*(n-1)/2) - (n-j)*((n-j)-1)/2 + i - j - 1;
    else k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;

    return k;
}



double QKP2( int n, int *p, int **Q, int *w, int *x_val,int C){

    GRBEnv *env  = 0;
    GRBVar *Elem = 0;
    double obj_get = -1;
    try{
        int m=1;

        // Create environment
        env = new GRBEnv();

        // Create initial model
        //GRBModel model = GRBModel(*env,"QMKP_CPX.lp");
        GRBModel model = GRBModel(*env);
        //model.read("QMKP_GRB_test.lp");

        //model.optimize();
        //return 0;




        model.set(GRB_StringAttr_ModelName, "QMKP GUROBI");


        GRBQuadExpr obj = 0;
        GRBVar x[n*m];

        //Normal profit
        for(auto i = 0u; i < n*m; ++i) {
            string s = "x" + to_string(i+1);
            //
            //x[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS,s);
            x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,s);
            int k = i % n;
            obj = obj + (p[k]) * x[i];////FIXME: double profit
        }

        //Quadratc profit
        int count_var=0;
        for(auto k=0u;k<m;k++){
            for(auto i = (count_var); i < n+count_var; ++i) {
                for(auto j = i+1; j < n+count_var; ++j) {
                    obj = obj + (Q[i-count_var][j - count_var])*(x[i] * x[j]);//FIXME: double profit
                }
            }
            count_var=count_var+n;
        }


        model.setObjective(obj, GRB_MAXIMIZE);

        //Constraint
        int cont_constr = 0;
        for(auto i = 0u; i < m; ++i) {
            GRBLinExpr expr = 0;
            for(auto j = cont_constr; j < n+cont_constr; ++j) {
                expr = expr + x[j] * w[j-cont_constr];
            }
            cont_constr = cont_constr +n;
            string s = "c" + to_string(i+1);
            model.addConstr(expr, GRB_LESS_EQUAL, C, s);
        }

        auto cont = 0u;
        for(auto i = 0u; i < n; ++i) {
            GRBLinExpr expr = 0;
            auto j = cont;
            while(j < (n*m)) {
                expr = expr + x[j];
                j=j+n;
            }
            string s = "c" + to_string((i+m)+1);
            model.addConstr(expr, GRB_LESS_EQUAL, 1.0, s);
            cont++;
        }
        //Set params
        model.set(GRB_IntParam_PreCrush,0);
        model.set(GRB_IntParam_DualReductions,1);
        model.set(GRB_IntParam_LazyConstraints, 0);
        model.set(GRB_IntParam_Threads,1);

        model.set(GRB_IntParam_PreQLinearize,1);


        /*
        if(root){
          model.set(GRB_DoubleParam_NodeLimit,1);
        }*/

        //LP File
        model.write("QMKP_GRB.lp");

        cout<<"Summary Parameters"<<endl;
        cout<<"GRB_IntParam_PreCrush:\t"<<model.get(GRB_IntParam_PreCrush)<<endl;
        cout<<"GRB_IntParam_DualReductions:\t"<<model.get(GRB_IntParam_DualReductions)<<endl;
        cout<<"GRB_IntParam_LazyConstraints:\t"<<model.get(GRB_IntParam_LazyConstraints)<<endl;
        cout<<"GRB_IntParam_Threads:\t"<<model.get(GRB_IntParam_Threads)<<endl;
        cout<<"GRB_IntParam_PreQLinearize:\t"<<model.get(GRB_IntParam_PreQLinearize)<<endl;
        /*
        if(TL>0){
          model.set(GRB_DoubleParam_TimeLimit,TL);
          cout<<"GRB_DoubleParam_TimeLimit:\t"<<model.get(GRB_DoubleParam_TimeLimit)<<endl;
        }*/
        model.set(GRB_DoubleParam_TimeLimit,1800.0);
        // Optimize model

        model.optimize();
        vector<unsigned> variables;
        //cout<<"------------------"<<endl;
        for(int i=0;i<n*m;i++){
            if(x[i].get(GRB_DoubleAttr_X)>0.5){
                cout<<x[i].get(GRB_StringAttr_VarName)<<"->"<<x[i].get(GRB_DoubleAttr_X)<<endl;
                //cout<<i<<endl;
                variables.push_back(i);
                x_val[i] = 1.0;
            }else{
                x_val[i] = 0.0;
            }
        }
        //cout<<endl;
        //cout<<"------------------"<<endl;
        obj_get =  model.get(GRB_DoubleAttr_ObjVal);
        std::cout << "Obj: " << obj_get<< std::endl;

    }

    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }

    catch (...) {
        cout << "Exception during optimization" << endl;
    }

    // Free environment/vars
    delete[] Elem;
    delete env;
    return obj_get;
}

double QKP( int n, int *p, int **Q, int *w, int *x_val,int C){

    int m=1;
    IloEnv env;

    IloModel model(env);


    //GRBQuadExpr obj = 0;
    IloNumVarArray x(env, n*m);
    IloExpr expr(env);
    /***************************************************************/
    /******************Objetive Function****************************/
    /***************************************************************/
    //Normal profit

    for(auto i = 0u; i < n*m; ++i) {
        string s = "x" + to_string(i+1);
        x[i]  = IloNumVar(env, 0, 1, IloNumVar::Bool, s.c_str());
        int k = i % n;
        expr = expr + (p[k]) * x[i];
        //lin = lin + (2*p[k]) * x[i];////FIXME: double profit
        //lin.addTerm(2*p[k] * x[i]);
    }

    //Quadratc profit
    int count_var=0;

    for(auto k=0u;k<m;k++){
        for(auto i = (count_var); i < n+count_var; ++i) {
            for(auto j = i+1; j < n+count_var; ++j) {
                expr = expr + (Q[i-count_var][j - count_var])*(x[i] * x[j]);//FIXME: double profit
                //quad.addTerm(1, n[j][k]ii,y[i][j]);
            }
        }
        count_var=count_var+n;
    }

    IloObjective obj(env, expr, IloObjective::Maximize);

    // Add the objective function to the model
    model.add(obj);
    expr.clear();


    //Constraint
    int cont_constr = 0;
    IloRangeArray constraint_1array(env, m);  // NORMAÑ Constraints 1)

    for(auto i = 0u; i < m; ++i) {

        for(auto j = cont_constr; j < n+cont_constr; ++j) {
            expr = expr + (x[j] * w[j-cont_constr]);
        }

        cont_constr = cont_constr +n;
        string s = "c" + to_string(i+1);
        constraint_1array[i] = IloRange(env, -IloInfinity, expr, C, s.c_str());
        expr.clear();

    }

    model.add(constraint_1array);

    auto cont = 0u;
    IloRangeArray constraint_2array(env, n);  // Constraints 2)
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
    expr.end();



    // Create the solver object
    IloCplex cplex(model);


    //cplex.setParam(IloCplex::Param::Read::DataCheck,1);
    //cplex.setParam(IloCplex::Param::Preprocessing::Reduce,3);
    //cplex.setParam(IloCplex::Param::OptimalityTarget,CPX_OPTIMALITYTARGET_OPTIMALGLOBAL);
    cplex.setParam(IloCplex::Param::Threads, 1);
    //if(flag_root){
    //    cplex.setParam(IloCplex::Param::MIP::Limits::Nodes,1);
    //}
    //cplex.setOut(env.getNullStream()) ;
    //if(TL > 0){
    cplex.setParam(IloCplex::Param::TimeLimit, 1800.0);
    //}

    // Export model to file (useful for debugging!)
    /*
    int i=index;
    stringstream ss;
    ss << i;
    string lp = "model.lp";
    cplex.exportModel(lp.c_str());
    */
    cplex.exportModel("QMKP_CLPX.lp");
    bool solved = false;
    clock_t before,end;
    double objval_aux = -1;
    before = clock();

    try {
        // Try to solve with CPLEX (and hope it does not raise an exception!)
        solved = cplex.solve();

    } catch(const IloException& e) {
        std::cerr << "\n\nCPLEX Raised an exception:\n";
        std::cerr << e << "\n";
        env.end();
        throw;
    }
    end = clock();
    if(solved) {
        // If CPLEX successfully solved the model, print the results
        std::cout << "\n\nCplex success!\n";
        std::cout << "\tStatus: " << cplex.getStatus() << "\n";
        std::cout << "Obj. Value : " << cplex.getObjValue()<<"\n";
        std::cout << "Get Nodes  : " << cplex.getNnodes()<<"\n";
        //std::cout << "Get Times  : " << cplex.getCplexTime()<<"\n";
        std::cout << "Get Bound  : " << cplex.getBestObjValue()<< "\n";
        std::cout << "Get Cutoff : " << cplex.getCutoff()<<"\n";
        std::cout << "Get Gap    : " << cplex.getMIPRelativeGap()<<"\n";
        std::cout << "Total Time : " << float(end-before)/CLOCKS_PER_SEC <<"\n";
        objval_aux = cplex.getObjValue();
        //*objval = (objval_aux);

        cout<<" ";
        for(int i=0;i< n*m;i++){
            if(cplex.getValue(x[i]) > .5){
                x_val[i] = 1;
            }else{
                x_val[i] = 0;
            }              }
    } else {
        std::cerr << "\n\nCplex error!\n";
        std::cerr << "\tStatus: " << cplex.getStatus() << "\n";
        std::cerr << "\tSolver status: " << cplex.getCplexStatus() << "\n";
    }

    env.end();
    return objval_aux;
}


double normal_qkp(
        int n,
        int *p,
        int *Q,
        int *w,
        int *C,
        vector<double> & ind_remain,
        vector<double> & qkp){

    int size_qkp = ind_remain.size();

    int *x  = (int *)calloc((int)size_qkp,sizeof(int));
    int **xQ = (int **)calloc((int)size_qkp, sizeof(int*));
    for(int i = 0; i < size_qkp; i++) xQ[i] = (int *)calloc((int)size_qkp,sizeof(int));
    //int xQ[MSIZE][MSIZE];

    int *wV  =(int *)calloc((int)size_qkp,sizeof(int));
    int CV = C[0];
    //cout<<endl;

    cout<<"Current Size QKP "<<size_qkp<<endl;
    cout<<"QKPA ";
    for(int av=0;av<ind_remain.size();av++){
        cout<<ind_remain[av]<<" ";
    }
    cout<<endl;

    /*
    for (int i = 0; i < size_qkp; i++) {
      int ind_i = ind_remain[i];
      wV[i] = w[ind_i];
      for(int j = i; j < size_qkp; j++) {
          int ind_j = ind_remain[j];
          if(ind_i == ind_j){
            xQ[i][j] = p[ind_i];
          }else{
                int ind_q =  ij2k_main(n, ind_i, ind_j);
                //xQ[i][j] =  (Q[ind_q])*100;
                xQ[i][j] =  (Q[ind_q]);
                xQ[j][i] =  (Q[ind_q]);
          }
        }
    }*/

    //cout<<endl;
    

    ofstream myfile;

    stringstream ss;
    ss << size_qkp;
    string str = ss.str();
    myfile.open("instance_m_"+str+".txt");
    myfile << "instance_m_"+str+".txt.\n";
    myfile << size_qkp <<"\n";
    for (int i = 0; i < size_qkp; i++) {
        int ind_i = ind_remain[i];
        myfile << p[ind_i] <<" ";
    }
    myfile << "\n";
    for (int i = 0; i < size_qkp; i++) {
        int ind_i = ind_remain[i];
        for(int j = i+1; j < size_qkp; j++) {
            int ind_j = ind_remain[j];
            int ind_q =  ij2k_main(n, ind_i, ind_j);
            myfile << (Q[ind_q]) << " ";
        }
        myfile << "\n";
    }
    myfile << "\n";
    myfile << 0 <<"\n";
    myfile << CV <<"\n";
    for (int i = 0; i < size_qkp; i++) {
        int ind_i = ind_remain[i];
        myfile << w[ind_i] <<" ";
    }
    myfile << "\n";
    myfile.close();

    int *pv  = (int *)calloc((int)size_qkp,sizeof(int));
    for (int i = 0; i < size_qkp; i++) {
        int ind_i = ind_remain[i];
        wV[i] = w[ind_i];
        pv[i] = p[ind_i];
        for(int j = i+1; j < size_qkp; j++) {
            int ind_j = ind_remain[j];
            int ind_q =  ij2k_main(n, ind_i, ind_j);
            xQ[i][j] =  (Q[ind_q]);
            xQ[j][i] =  (Q[ind_q]);

        }
    }

    /*
    for(int i=0;i<size_qkp;i++){
      for(int j=0;j<size_qkp;j++){
        cout<<xQ[i][j]<<" ";
      }
      cout<<endl;
    }
    */
    //double objval_aux =QKP( int n, int *p, int **Q, int *w, int C);
    double objval_aux = QKP2( size_qkp,pv, xQ, wV, x,CV);
    //double objval_aux = (double)quadknap(size_qkp, CV, xQ , wV, x,1000);
    //if(sum_lambda!=0.0){
    //objval_aux=objval_aux/100;
    //cout<<objval_aux<<endl;

    cout<<"QKPR ";
    for(int i=0;i<size_qkp;i++){
        cout<<x[i]<<" ";
    }
    cout<<endl;

    for (int i = 0; i < n; i++) {
        std::vector<double>::iterator it = find (ind_remain.begin(), ind_remain.end(), i);
        if(it != ind_remain.end()){
            int index_i = std::distance(ind_remain.begin(), it);
            qkp.push_back((double)x[index_i]);
        }else{
            qkp.push_back(0.0);
        }
    }
    //cout<<qkp.size()<<endl;
    free(x);
    free(wV);
    for(int i = 0; i < size_qkp; i++) free(xQ[i]);
    free(xQ);
    return objval_aux;
}



int qkp_grdy(int n_input, int m, int *p, int *Q, int *w, int *C, vector<vector<int> > & x){

    double aux = 0.0;
    int n = n_input;
    vector<double> available;
    for(int i=0;i<n_input;i++) available.push_back(i);
    vector<vector<double> > fnl;


    double start;
    double end;
    start = omp_get_wtime();

    for(int k=0;k<m;k++){
        vector<double> qkp;
        /*
        cout<<"Available ";
        for(int av=0;av<available.size();av++){
            cout<<available[av]<<" ";
        }
        cout<<endl;
            */
        aux = aux + normal_qkp(n,p, Q, w, C, available,qkp);

        fnl.push_back(qkp);

        vector<double>  remain_cycle(n,0);
        for(int k=0;k<fnl.size();k++){
            for(int i=0;i<n;i++){
                remain_cycle[i] = remain_cycle[i] +  fnl[k][i];
            }
        }
        available.clear();

        for(int i=0;i<remain_cycle.size();i++){
            if(remain_cycle[i]==0.0){
                available.push_back((double)i);
            }
        }

    }

    cout<<endl;
    for (int k = 0; k < m; k++)
    {
        for (int i = 0; i < n; ++i)
        {
            cout<<fnl[k][i]<<" ";
        }
        cout<<endl;
    }

    cout<<endl;
    double acum =0.0;
    for(int k=0;k<m;k++){
        for(int i=0;i<n;i++){
            if(fnl[k][i]==1){
                acum = acum + p[i];
                for(int j=i+1;j<n;j++){
                    int index_q =  ij2k_main(n,i,j);
                    acum = acum + ((fnl[k][j])*Q[index_q]);
                }
            }
        }
    }

    end = omp_get_wtime();
    cout<<"Value "<<acum<<endl;
    cout<<"Time "<<end - start<<endl;

    Solution * s = new Solution(fnl);
    //s->show();
    s->check(n ,m,w ,C,p,Q);

    for(int k=0;k<m;k++){
        vector<int> aux;
        for(int i=0;i<n;i++){
            aux.push_back(fnl[k][i]);
        }
        x.push_back(aux);
    }


    return 0;

}
