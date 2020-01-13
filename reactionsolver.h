#ifndef REACTIONSOLVER_H
#define REACTIONSOLVER_H
#include <vector>
#include <stdio.h>
#include "reaction.h"

class simulationData;
//number of input species

#define I_NUM 3
class reactionSolver
{
public:
    reactionSolver(simulationData *mdata);

   // simulationData* m_data;
   // std::vector<simulationData::simulationField* > m_iFields;
    std::vector<reaction* > m_reactions;

    double n_n,ne_i,nars_i; //initial values
    double ne_o,eps_o,neps_o,nars_o,narp_o; //output values
    double ne_0,nars_0; //values from previous timestep
    double mask;

    double rhs_ne,rhs_neps,rhs_nars,rhs_narp; //right hand side
    double nu_ne,nu_neps,nu_nars,nu_narp; //viscosities


   void solve(int itn); //fully implicit no diffusuion
   bool solve_diffuse(int itn); //fully implicit with diffusuion

private:

    double dValues[I_NUM];
    double oValues[I_NUM];//ne, eps,n_ars

    double f[I_NUM];
    double df[I_NUM][I_NUM];




public:

    double detM_3x3(double mat[3][3])
    {
        double ans;
        ans = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
              - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
              + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
        return ans;
    }
    void solve3x3(double A[3][3],double b[3],double x[3]) //Ax=b
    {
        // Matrix d using coeff as given in cramer's rule
        double D = detM_3x3(A);

        // Matrix d1 using coeff as given in cramer's rule
        static double d1[3][3] = {
            { b[0], A[0][1], A[0][2] },
            { b[1], A[1][1], A[1][2] },
            { b[2], A[2][1], A[2][2] },
        };
        // Matrix d2 using coeff as given in cramer's rule
        static double d2[3][3] = {
            { A[0][0], b[0], A[0][2] },
            { A[1][0], b[1], A[1][2] },
            { A[2][0], b[2], A[2][2] },
        };
        // Matrix d3 using coeff as given in cramer's rule
        static double d3[3][3] = {
            { A[0][0], A[0][1], b[0] },
            { A[1][0], A[1][1], b[1] },
            { A[2][0], A[2][1], b[2] },
        };

        // Calculating Determinant of Matrices d, d1, d2, d3

        double D1 = detM_3x3(d1);
        double D2 = detM_3x3(d2);
        double D3 = detM_3x3(d3);
        /*printf("D is : %lf \n", D);
        printf("D1 is : %lf \n", D1);
        printf("D2 is : %lf \n", D2);
        printf("D3 is : %lf \n", D3);*/

        // Case 1
        //if (D != 0) {
            // Coeff have a unique solution. Apply Cramer's Rule
            x[0] = D1 / D;
            x[1] = D2 / D;
            x[2] = D3 / D; // calculating z using cramer's rule
            printf("Value of x is : %lf\n", x[0]);
            printf("Value of y is : %lf\n", x[1]);
            printf("Value of z is : %lf\n", x[2]);
        //}
        // Case 2
        /*else {
            if (D1 == 0 && D2 == 0 && D3 == 0)
                printf("Infinite solutions\n");
            else if (D1 != 0 || D2 != 0 || D3 != 0)
                printf("No solutions\n");
        }*/
    }

};

#endif // REACTIONSOLVER_H
