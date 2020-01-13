
#include <stdio.h>
#include <stdlib.h>

#include "phi_mult.h"
#include <math.h>

double jacobi_N(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y],double mask[N_X][N_Y],double maskValue[N_X][N_Y], int itn,int nx,int ny)
{
    int i,j,n;
    double res=0;
    double eps =1.0;
    double a,b_p,b_m,c_p,c_m;
    a=par.a;
    b_p=par.bp;
    b_m=par.bm;
    c_p=par.cp;
    c_m=par.cm;
    for(n=0;n<itn;n++)
    {
        for (i=1; i<nx; i++)
        {
            for (j=1; j<ny; j++)
            {
                if (j==1)
                {
                    if (par.s_bc_type==0)//fixed_value
                        field[i][j-1]=par.s_bc_val;
                    if (par.s_bc_type==1)//fixed_gradient
                        field[i][j-1]=field[i][1]-par.dy*par.s_bc_val;
                    if (par.s_bc_type==2)//cyclic
                        field[i][j-1]=field[i][ny-1];
                }

                if (j==ny-1)
                {
                    if (par.n_bc_type==0)//fixed_value
                        field[i][j+1]=par.n_bc_val;
                    if (par.n_bc_type==1)//fixed_gradient
                        field[i][j+1]=field[i][j]+par.dy*par.n_bc_val;
                    if (par.n_bc_type==2)//cyclic
                        field[i][j+1]=field[i][1];
                }

                if (i==1)
                {
                    if (par.w_bc_type==0)//fixed_value
                        field[i-1][j]=par.w_bc_val;
                    if (par.w_bc_type==1)//fixed_gradient
                        field[i-1][j]=field[1][j]-par.dx*par.w_bc_val;
                    if (par.w_bc_type==2)//cyclic
                        field[i-1][j]=field[nx-1][j];
                }

                if (i==nx-1)
                {
                    if (par.e_bc_type==0)//fixed_value
                        field[i+1][j]=par.e_bc_val;
                    if (par.e_bc_type==1)//fixed_gradient
                        field[i+1][j]=field[i][j]+par.dx*par.e_bc_val;
                    if (par.e_bc_type==2)//cyclic
                        field[i+1][j]=field[1][j];
                }

                field[i][j]=((rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))
                        /a)*mask[i][j]+maskValue[i][j];

            }
        }
    }
    return 0;
}
double interp_up(double field[N_X][N_Y], double field2[N_X][N_Y], int nx,int ny)
{
    int i,j,k,l;

    for(k=1; k<(nx)/2; k++)
    {
        i=k*2;
        for (l=1; l<(ny)/2; l++)
        {
            j=l*2;
            field[i][j]+=field2[k][l];
            field[i-1][j]+=0.5*(field2[k][l]+field2[k-1][l]);
            field[i][j-1]+=0.5*(field2[k][l]+field2[k][l-1]);
            field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
        }

        l=(ny)/2; j=l*2;

        field[i][j-1]+=0.5*(field2[k][l]+field2[k][l-1]);
        field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
    }

    k=(nx)/2;
    i=k*2;

    for (l=1;l<(ny)/2;l++)
    {
        j=l*2;
        field[i-1][j]+=0.5*(field2[k][l]+field2[k-1][l]);
        field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
    }

    l=(ny)/2; j=l*2;
    field[i-1][j-1]+=0.25*(field2[k][l]+field2[k][l-1]+field2[k-1][l]+field2[k][l]);
}


double multigrid_N(INPUT_PARAM par, double** field, double** rhs, double** mask, double** maskValue, int itn,int N)
{
    int i,j,l,k,nn;
    double a,b_p,b_m,c_p,c_m;
    double res=0;
    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=par.bp;//-1.0/(dx*dx);
    b_m=par.bm;//-1.0/(dx*dx);
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);
    static double rhs_[10][N_X][N_Y];
    static double field_[10][N_X][N_Y]; 
    static double mask_[10][N_X][N_Y];
    static double maskValue_[10][N_X][N_Y];

    for (i=0; i<N_X; i++)
    {
        for (j=0;j<N_Y;j++)
        {
            field_[0][i][j]=field[i][j];
            rhs_[0][i][j]=rhs[i][j]/*-rho_[i][j]*/;
            mask_[0][i][j]= mask[i][j];
            maskValue_[0][i][j] = maskValue[i][j];
        }
    }

    for (nn=0;nn<N;nn++)
    {
        jacobi_N(par,field_[nn],rhs_[nn],mask_[nn], maskValue_[nn],itn,(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));

        for (i=1; i<(N_X-1)/pow(2,nn+1); i++)
        {
            for (j=1;j<(N_Y-1)/pow(2,nn+1);j++)
            {
                k=i*2; l=j*2;
                rhs_[nn+1][i][j]=rhs_[nn][k][l] - 4.0*(a*field_[nn][k][l]+b_p*field_[nn][k+1][l]+b_m*field_[nn][k-1][l]+c_p*field_[nn][k][l+1]+c_m*field_[nn][k][l-1]);
                mask_[nn+1][i][j]=mask_[nn][k][l];
                maskValue_[nn+1][i][j]=0.0;
                field_[nn+1][i][j]=0.0;
            }
        }
    }

    jacobi_N(par,field_[N],rhs_[N],mask_[N], maskValue_[N],itn,(N_X-1)/pow(2,N),(N_Y-1)/pow(2,N));

    for (nn=N-1;nn>=0;nn--)
    {
        interp_up(field_[nn],field_[nn+1],(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));
        jacobi_N(par,field_[nn],rhs_[nn],mask_[nn], maskValue_[nn],itn,(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));
    }

    for (i=0; i<N_X; i++)
    {
        for (j=0;j<N_Y;j++)
        {
            field[i][j]=field_[0][i][j];
        }
    }
    return res;
}

