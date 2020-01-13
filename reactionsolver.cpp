#include "reactionsolver.h"
#include "simulationdata.h"
#include <qdebug.h>
#include <stdio.h>
#include <math.h>
reactionSolver::reactionSolver(simulationData* mdata)
{

    for (int i=0;i<mdata->m_reactions.size();i++)
        m_reactions.push_back(mdata->m_reactions[i]);
    qDebug()<<m_reactions.size()<< "AA!!";

}

void reactionSolver::solve(int itn)
{
    /*   oValues[0]=ne_i;
    oValues[1]=eps_i;
    oValues[2]=nars_i;
    double k[7];
    double dk[7];
    double dE[7]; //energy gain/loss of reaction

    double x0=ne_0;
    double y0=nars_0;


    double A,D,C,B,E,F;

    for (int i=0;i<7;i++)
    {
        dE[i]=m_reactions[i]->getDe();
    }


    double dt_cur=dt;
    double x,y,eps;
    x=ne_i/x0;y=nars_i/y0;

    eps=eps_i;
   //  qDebug()<<"x0="<<x0<<" y0="<<y0<<" dt="<<dt<<" eps="<<eps_i<<" nn="<<n_n;

    for (int nn=0;nn<itn;nn++)
    {
        for (int i=0;i<7;i++)
        {
            k[i]=m_reactions[i]->getRate(eps*2.0/3.0);
            dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }
        /*        //n_e
                (oValues[0] - ne_i) / dt = oValues[0]*(k[3]*n_n + k[4]*oValues[2]) + k[5]*oValues[2]*nars_i;
                -(oValues[0]-ne_i)/dt+rhs_ne + oValues[0]*(k[3]*n_n +k[4]*oValues[2])+k[5]*oValues[2]*oValues[2];
                //n_ars
                 (oValues[2]-nars_i)/dt=oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2]) -2.0*k[5]*oValues[2]*nars_i;


     (x - x0) / dt - nu * (xr - 2.0 * x + xl) / (dz*dz) - rhs = ...
     x*(1/dt+2.0*nu/(dz*dz)) -x0/dt-nu*(xr+xl)/(dz*dz) - rhs =...




        a=dt*k[3]*n_n;   b=dt*k[1]*n_n; c=dt*nars_i*k[5];
        d=dt*k[4];  e=dt*(k[2]+k[4]);

        //ne
        x - (ne_i+rhs*dt) = x*(a + d*y) + c*y;
        //n_ars
         y-(nars_i+rhs*dt)=x*(b -e*y) -2.0*c*y;
                */

    /*    A=(dt_cur*k[3]*n_n-1.0);
        D=dt_cur*k[4]*y0;
        C=dt_cur*k[5]*y0*(y0/x0);
        //Ax+Dxy+Cy+(1+rhs*dt/x0)=0  ne eqn;  x is x/x0 y is y/y0

        B=dt_cur*k[1]*n_n*(x0/y0);
        E=-dt_cur*(k[2]+k[4])*x0;
        F=-1.0-2.0*dt_cur*k[5]*y0;
        //Bx+Exy+Fy+(1+rhs*dt/y0)=0  nars eqn;  x is x/x0 y is y/y0


       // qDebug()<<"nn="<<nn<<" A="<<A<<" D="<<D<<" C="<<C;
       // qDebug()<<"B="<<B<<" E="<<E<<" F="<<F;



        double f_ne,f_nars, dfx_ne,dfy_ne, dfx_nars,dfy_nars;

        f_ne=A*x+D*x*y+C*y+(1+rhs_ne*dt/x0);
        dfx_ne=A+D*y;
        dfy_ne=C+D*x;

        f_nars=B*x+E*x*y+F*y+(1+rhs_nars*dt/y0);
        dfx_nars=B+E*y;
        dfy_nars=F+E*x;

        double det = dfx_ne*dfy_nars - dfy_ne*dfx_nars;

        double d_x = -f_ne*dfy_nars + f_nars*dfy_ne;
        double d_y = f_ne*dfx_nars - f_nars*dfx_ne;

        x+=d_x;
        y+=d_y;

        //here we've got ne and ars now for two other species (en and arp)

      //  qDebug()<<"f_ne"<<f_ne<<" f_nars="<<f_nars<<" x="<<x<<" y="<<y;


/*
          eps = eps_i +
                        dt*n_n * dE[1]*k[1] + dt*n_n*k3*dE[3] - dt*n_n*k3*eps  +
                        dt*y*y0 *dE[2]*k[2] + dt*y*y0 *k[4]*dE[4] - dt*y*y0 *k[4]*eps   +
                       (dt/ (x*x0))* k[5] * y*y0 * nars_i*dE[5] - (dt/ (x*x0))* k[5] * y*y0 * nars_i*eps  ;
 */

    /*      eps  = (eps_0 + dt*(rhs_eps-rhs_ne)/(x*x0) +
                      dt*n_n * dE[1]*k[1] + dt*n_n*k[3]*dE[3] +
                      dt*y*y0 *dE[2]*k[2] + dt*y*y0 *k[4]*dE[4]    +
                     (dt/ (x*x0))* k[5] * y*y0 * nars_0*dE[5])/(1.0  + dt*(n_n*k[3] + y*y0 *k[4] +  ( k[5] * y*y0 * nars_0/ (x*x0)) )) ;


    }


    ne_o=x*x0;//oValues[0];
    eps_o=eps;//oValues[1];
    nars_o=y*y0;//oValues[2];*/

}

bool reactionSolver::solve_diffuse(int itn)
{

    double k[7];
    //  double dk[7];
    double dE[7]; //energy gain/loss of reaction

    double x0=ne_0;
    double y0=nars_0;


    double A,D,C,B,E,F;

    /*for (int i=0;i<7;i++)
    {
        dE[i]=m_reactions[i]->getDe();
    }*/



    double x,y,z,eps;


    x=ne_i;y=nars_i;

    eps=eps_o;
    for (int nn=0;nn<itn;nn++)
    {

        // lets get neps (z=neps/x0)
        /*
                z*(1/dt+2.0*nu/(dz*dz)) -z0/dt-nu*(zr+zl)/(dz*dz) - rhs =
                                x*x0*(
                                n_n *( dE[1]*k[1] +  dE[3]*k[3]) +
                                y*y0*( dE[2]*k[2] +  dE[4]*k[4])
                                ) +
                                 (dE[5]) * k[5] * y*y0 * nars_i  ;
        */


        for (int i=0;i<7;i++)
        {
          //  k[i]=m_reactions[i]->getRate(fmax(fmin(eps*2.0/3.0,40),0.0001));
            //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }
        /*        //n_e
                (oValues[0] - ne_i) / dt = oValues[0]*(k[3]*n_n + k[4]*oValues[2]) + k[5]*oValues[2]*nars_i;
                -(oValues[0]-ne_i)/dt+rhs_ne + oValues[0]*(k[3]*n_n +k[4]*oValues[2])+k[5]*oValues[2]*oValues[2];
                //n_ars
                 (oValues[2]-nars_i)/dt=oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2]) -2.0*k[5]*oValues[2]*nars_i;


     (x - x0) / dt - nu * (xr - 2.0 * x + xl) / (dz*dz) - rhs = ...
     x*(1/dt+2.0*nu/(dz*dz)) -x0/dt-nu*(xr+xl)/(dz*dz) - rhs =...

       x/dt_x - (x0/dt-nu*(xr+xl)/(dz*dz) + rhs) =

                  x - rhs=r*mask ;

       x/dt_x -rhs_x =


        a=k[3]*n_n;   b=k[1]*n_n; c=nars_i*k[5];
        d=k[4];  e=(k[2]+k[4]);

        //ne
        x - rhs_x = (x*(k[3]*n_n + k[4]*y) + nars_i*k[5]*y)*mask[i][j];
        //n_ars
         y- rhs_y=(x*(k[1]*n_n -(k[2]+k[4])*y) -2.0*nars_i*k[5]*y)*mask[i][j];
                */

        //x=((m_aRHS[i][j]-(bp*m_field->arr[i+1][j]+bm*m_field->arr[i-1][j]+cp*m_field->arr[i][j+1]+cm*m_field->arr[i][j-1]))
        //        /a)*m_mask[i][j]+m_maskValue[i][j] +R*mask/a;
        /*ne=rhs+R*mask/a;
        -ne+rhs+R*mask/a=0;

        nars=rhs+R*mask/a;*/


        /*
        A=-1.0;//(mask*k[3]*n_n-1.0);
        D=0.0;//k[4]*y0*mask;
        C=0.0;//k[5]*y0*(y0/x0)*mask;
        //Ax+Dxy+Cy+(1+rhs*dt/x0)=0  ne eqn;  x is x/x0 y is y/y0

        B=0.0;//k[1]*n_n*(x0/y0)*mask;
        E=0.0;//-(k[2]+k[4])*x0*mask;
        F=-1.0;//-2.0*k[5]*y0*mask;
        //Bx+Exy+Fy+(1+rhs*dt/y0)=0  nars eqn;  x is x/x0 y is y/y0

        double f_ne,f_nars, dfx_ne,dfy_ne, dfx_nars,dfy_nars;

        f_ne=A*x+D*x*y+C*y+(rhs_ne);
        dfx_ne=A+D*y;
        dfy_ne=C+D*x;

        f_nars=B*x+E*x*y+F*y+(rhs_nars);
        dfx_nars=B+E*y;
        dfy_nars=F+E*x;

        double det = dfx_ne*dfy_nars - dfy_ne*dfx_nars;

        double d_x = -f_ne*dfy_nars + f_nars*dfy_ne;
        double d_y = f_ne*dfx_nars - f_nars*dfx_ne;

        x+=d_x;
        y+=d_y;

        if ((x<0.9) || (x>1.1) || (y<0.9) || (y>1.1) )
            return false;//x=0.9;*/




        z=(
                    rhs_neps/* + x*(n_n *( dE[1]*k[1] +  dE[3]*k[3]) + y*y0*( dE[2]*k[2] +  dE[4]*k[4])) + dE[5]*k[5]* y*y0/x0*nars_0*/
                );

        //    z=0.1;
        eps  = z/x;//z/x;


    } //n

    //now lets get arp
    //narp_o*(1/dt+2.0*nu/(dz*dz)) -x0/dt-nu*(xr+xl)/(dz*dz) - rhs = (ne_o)*(k[3]*n_n + k[4]*y*y0) + k[5]*nars_o*nars_i;


    ne_o=x;//*x0;//oValues[0];
    neps_o=z;//*x0;//oValues[1];
    eps_o=z/x;
    nars_o=y;//*y0;//oValues[2];
    narp_o= rhs_narp ;// + ((ne_o)*(k[3]*n_n + k[4]*nars_o) + k[5]*nars_o*nars_0)*mask;
    return true;
} //eof
