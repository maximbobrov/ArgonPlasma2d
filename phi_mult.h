

#ifndef PHI_MULT_H
#define PHI_MULT_H

//#include "globals.h"

//for multigrid
#define N_X 257
#define N_Y 257
typedef struct
{
    double a,bm,bp,cm,cp;
    int w_bc_type,e_bc_type,n_bc_type,s_bc_type; //0--fixed value, 1--fixed gradient,2 --cyclic; 3 --init
    int w_bc_val,e_bc_val,n_bc_val,s_bc_val;
    double dx,dy;
}INPUT_PARAM;
double jacobi_N(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y],double mask[N_X][N_Y],double maskValue[N_X][N_Y], int itn,int nx,int ny);
double multigrid_N(INPUT_PARAM par, double** field, double** rhs, double** mask, double** maskValue, int itn,int N);
#endif
