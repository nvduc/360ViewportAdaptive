/* Copyright 2018, Gurobi Optimization, LLC */

/* This example formulates and solves the following simple MIP model:

     maximize    x +   y + 2 z
     subject to  x + 2 y + 3 z <= 4
                 x +   y       >= 1
     x, y, z binary
*/

#include "gurobi_c++.h"
using namespace std;

int
main(int   argc,
     char *argv[])
{
    int No_tile = 64, i, j;
    const int NO_VER = 7;
    int LEN = 200;
    GRBVar** x;
    char*** x_name;
    double** d; // tiles' distortions
    double*w; // tiles' weights (overlaped area/ viewport area)
    double*p;
    double** TILE_SEG_BR; 
    double BW;
    int No_tile_h;
  try {
    GRBEnv env = GRBEnv();

    GRBModel model = GRBModel(env);

    // init
    x = new GRBVar*[No_tile];
    x_name = new char**[No_tile];
    for(i=0; i < No_tile; i++){
        x[i] = new GRBVar[NO_VER];
        x_name[i] = new char*[NO_VER];
        for(j=0; j < NO_VER; j++)
            x_name[i][j] = new char[LEN];
    }

    // Create variables
    for(i=0; i < No_tile; i++){
        for(j=0; j < NO_VER; j++){
            x_name[i][j] = new char[LEN];
            sprintf(x_name[i][j], "x_%d_%d", i, j);
            x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, x_name[i][j]);
        }
    }
    // objective funtion
    GRBQuadExpr obj;
    GRBQuadExpr avgVPQL = 0;
    GRBQuadExpr varVPQL = 0;
    for(i=0; i < No_tile; i++){
        for(j=0; j < NO_VER; j++){
            avgVPQL += p[i]*w[i]*d[i][j]*x[i][j];
        }
    }
    // for(i=0; i < No_tile; i++){
    //     for(j=0; j < NO_VER; j++){
    //         varVPQL += p[i]*w[i]*s[i] * d[i][j] * d[i][j] * x[i][j];
    //         for(ii=0; ii < No_tile; ii++){
    //             for(jj=0; jj < No_tile; jj++){
    //                 varVPQL += p[i]*w[i]*s[i]  * x[i][j] * p[ii] * w[ii] * d[ii][jj] * x[ii]*x[jj];
    //             }
    //         }
    //     }
    // }
    model.setObjective(avgVPQL, GRB_MINIMIZE);
    // set constraints
    // c1: total tiles' bitrates does not exceed bandwidth
    GRBLinExpr constr_1 = 0;
    GRBLinExpr* constr_2 = new GRBLinExpr[No_tile];
    for(i=0; i < No_tile; i++){
        constr_2[i] = 0;
        for(j=0; j < NO_VER; j++){
            constr_1 += TILE_SEG_BR[i][j] * x[i][j];
            constr_2[i] += x[i][j];
        }
        model.addConstr(constr_2[i] <= 1, "constr_2");
    }
    model.addConstr(constr_1 <= BW, "constr_1");
    // c2: at most one version will be selected for a specific tile

    // GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    // GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    // GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    // // Set objective: maximize x + y + 2 z

    // model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

    // // Add constraint: x + 2 y + 3 z <= 4

    // model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

    // // Add constraint: x + y >= 1

    // model.addConstr(x + y >= 1, "c1");

    // Optimize model

    model.optimize();
    for(i=0; i < No_tile; i++){
        for(j=0; j < NO_VER; j++){
            if(x[i][j].get(GRB_DoubleAttr_X) == 1)
                cout << j << " ";
        }
        if((i+1)%No_tile_h == 0)
            cout << endl;
    }

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
