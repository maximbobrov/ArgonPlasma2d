#include "simulationsolver.h"
#include "simulationtools.h"
#include "phi_mult.h"
#include "math.h"
#include <QDebug>


simulationSolver::simulationSolver(simulationData* ipData/*= nullptr*/)
{
    m_pData = ipData;


    if (ipData!=nullptr)
    {

        int cellsX =  ipData->getCellsXNumber();
        int cellsY =  ipData->getCellsYNumber();

        m_a = new double*[cellsX];
        m_bm  = new double*[cellsX];
        m_bp  = new double*[cellsX];
        m_cm  = new double*[cellsX];
        m_cp  = new double*[cellsX];

        for (int i = 0; i < cellsX; ++i) {
            m_a[i] = new double [cellsY];
            m_bm[i] = new double [cellsY];
            m_bp[i] = new double [cellsY];
            m_cm[i] = new double [cellsY];
            m_cp[i] = new double [cellsY];
        }

    }
    m_aRHS = nullptr;
    m_field = nullptr;
    m_aNu = nullptr;
    m_aMu = nullptr;
}

simulationSolver::~simulationSolver()
{

}

double simulationSolver::solve(int iNumberIteration)
{
    double res = 0.0;
    if(m_field != nullptr)
    {
        double dx = m_pData->getDx();
        double dy = m_pData->getDy();
        double dt = m_pData->getDt();

        updateMatrix();

        simulationData::simulationParameters* pParams = m_pData->getParameters();
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();
            setBc();

            for (int i = 1; i < m_field->cellsX-1; ++i)
            {
                for (int j = 1; j < m_field->cellsY-1; ++j)
                {
                    //double D = m_aNu[i][j];
                    double a = m_a[i][j];//1.0/dt+D*((2.0)/(dx*dx)+(2.0)/(dy*dy));
                    double bp = m_bp[i][j];//-D/(dx*dx);
                    double bm = m_bm[i][j];
                    double cp = m_cp[i][j];
                    double cm = m_cm[i][j];
                    m_field->arr[i][j]=((m_aRHS[i][j]-(bp*m_field->arr[i+1][j]+bm*m_field->arr[i-1][j]+cp*m_field->arr[i][j+1]+cm*m_field->arr[i][j-1]))
                            /a)*m_mask[i][j]+m_maskValue[i][j];




                    // 0=-m_field->arr[i][j] + m_maskValue[i][j] + ((m_aRHS[i][j]-(bp*m_field->arr[i+1][j]+bm*m_field->arr[i-1][j]+cp*m_field->arr[i][j+1]+cm*m_field->arr[i][j-1]))
                    //         /a)*m_mask[i][j];
                }
            }

        }
    }
}


double simulationSolver::getRhs()
{
    return -1;
}

double simulationSolver::getRhsAt(int i, int j)
{
    return 0;
}

double simulationSolver::getNewtonRhs(int i,int j)
{
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();
    double dt = m_pData->getDt();
    double ** phi  = m_pData->getFieldPhi()->arr;


   /*     double ar= m_charge *  (phi[i+1][j] - phi[i][j]) * m_aMu[i][j];
        double al= m_charge *  (phi[i][j] - phi[i-1][j]) * m_aMu[i][j];
        double au= m_charge *  (phi[i][j+1] - phi[i][j]) * m_aMu[i][j];
        double ad= m_charge *  (phi[i][j] - phi[i][j-1]) * m_aMu[i][j];

      double Grij = m_mask[i+1][j]*(-m_aNu[i][j] - (m_charge>0)*ar) / dx;
      double Glij = m_mask[i-1][j]*( m_aNu[i][j] - (m_charge>0)*al) / dx;
      double Guij = m_mask[i][j+1]*(-m_aNu[i][j] - (m_charge>0)*au) / dy;
      double Gdij = m_mask[i][j-1]*( m_aNu[i][j] - (m_charge>0)*ad) / dy;
      double Gr = m_mask[i+1][j]*(m_field->arr[i+1][j] * m_aNu[i][j]   + (m_charge<=0)*m_field->arr[i][j] * ar) / dx;
      double Gl = - m_mask[i-1][j]*(m_field->arr[i-1][j] * m_aNu[i][j] + (m_charge<=0)*m_field->arr[i][j] * al) / dx;
      double Gu = m_mask[i][j+1]*(m_field->arr[i][j+1] * m_aNu[i][j]   + (m_charge<=0)*m_field->arr[i][j] * au) / dy;
      double Gd = - m_mask[i][j-1]*(m_field->arr[i][j-1] * m_aNu[i][j] + (m_charge<=0)*m_field->arr[i][j] * ad) / dy;*/


      ///////////////



      double eu=/*calcE(j+1);*/phi[i][j+2] - phi[i][j+1];//phi[j+2] - phi[j+1];               //j+1
      double ed=/*calcE(j-1);*/phi[i][j+1] - phi[i][j];

      double er=/*calcE(j+1);*/phi[i+2][j] - phi[i+1][j];//phi[j+2] - phi[j+1];               //j+1
      double el=/*calcE(j-1);*/phi[i+1][j] - phi[i][j];


      double au= m_charge * eu *m_aMu[i][j] ;
      double ad= m_charge * ed *m_aMu[i][j];

      double ar= m_charge * er *m_aMu[i][j] ;
      double al= m_charge * el *m_aMu[i][j];


      double Ga_u = m_mask[i][j+1]*(-m_aNu[i][j] ) / dy;
      double Ga_d = m_mask[i][j-1]*( m_aNu[i][j] ) / dy;
      double Ga_r = m_mask[i+1][j]*(-m_aNu[i][j] ) / dx;
      double Ga_l = m_mask[i-1][j]*( m_aNu[i][j] ) / dx;



      double Gu = m_mask[i][j+1]*(m_field->arr[i][j+1] * m_aNu[i][j]   +2.0*(m_field->arrPrev[i][j+1]  )* au) / dy;
      double Gd = m_mask[i][j-1]*(-m_field->arr[i][j-1] * m_aNu[i][j]  +2.0*(m_field->arrPrev[i][j]) * ad) / dy;

      double Gr = m_mask[i+1][j]*(m_field->arr[i+1][j] * m_aNu[i][j]   +2.0*(m_field->arrPrev[i+1][j]  )* au) / dx;
      double Gl = m_mask[i-1][j]*(-m_field->arr[i-1][j] * m_aNu[i][j]  +2.0*(m_field->arrPrev[i][j]) * ad) / dx;



      ////////////////

    double rhs=getRhsAt(i,j);

    double a = 1.0/dt -  (Ga_u - Ga_d) / dy -  (Ga_r - Ga_l) / dx;

    double res=((rhs + m_field->arrPrev[i][j] / dt + ((Gu - Gd) / dy) + ((Gr - Gl) / dx) ) / a)*m_mask[i][j];
    return res;
}

void simulationSolver::getStepEuler()
{
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        for (int j = 0; j < m_field->cellsY; ++j)
        {
            m_field->arrPrev[i][j] = m_field->arr[i][j];
        }
    }
}

void simulationSolver::setBc()
{
    for (int j = 0; j < m_field->cellsY; ++j)
    {
        m_field->arr[0][j] =0.0;// m_field->arr[m_field->cellsX-2][j];
        m_field->arr[m_field->cellsX-1][j] =0.0;// m_field->arr[1][j];
        m_field->arrPrev[0][j] =0.0;// m_field->arrPrev[m_field->cellsX-2][j];
        m_field->arrPrev[m_field->cellsX-1][j] =0.0;// m_field->arrPrev[1][j];
    }
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        m_field->arr[i][0] =0.0;// m_field->arr[i][m_field->cellsY-2];
        m_field->arr[i][m_field->cellsY-1] = 0.0;//m_field->arr[i][1];
        m_field->arrPrev[i][0] = 0.0;//= m_field->arrPrev[i][m_field->cellsY-2];
        m_field->arrPrev[i][m_field->cellsY-1]= 0.0;// = m_field->arrPrev[i][1];
    }

    simulationData::boundary* b= &(this->m_pData->getParameters()->bound);
    for (int n=0; n<b->num; n++)
    {
        int ii=b->arrI[n];
        int jj=b->arrJ[n];
        int iip=b->arrI_inter[n];
        int jjp=b->arrJ_inter[n];
        m_maskValue[ii][jj]=m_field->arr[iip][jjp];
    }
}

void simulationSolver::updateMatrix()
{
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();
    double dt = m_pData->getDt();

    int cellsx=m_pData->getCellsXNumber();
    int cellsy=m_pData->getCellsYNumber();

    for (int i = 0; i < cellsx; ++i)
    {
        for (int j = 0; j < cellsy; ++j)
        {
            double D = m_aNu[i][j];
            double a = 1.0/dt+D*((2.0)/(dx*dx)+(2.0)/(dy*dy));
            double bp = -D/(dx*dx);
            double bm = bp;
            double cp = -D/(dy*dy);
            double cm = cp;

            m_a[i][j] = a;
            m_bm[i][j] = bm;
            m_bp[i][j] = bp;
            m_cm[i][j] = cp;
            m_cp[i][j] = cm;
        }
    }
}

double simulationSolver::init(double value)
{
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        for (int j = 0; j < m_field->cellsY; ++j)
        {
            m_field->arr[i][j] = value;
            m_field->arrPrev[i][j] = value;
        }
    }
}

solverNe::~solverNe()
{

}

solverNe::solverNe(simulationData* pData):simulationSolver(pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldNe();
    m_aNu = m_pData->getParameters()->arrDe;
    m_aRHS = new double*[m_field ->cellsX];
    m_mask = m_pData->getParameters()->arrMaskNe;
    m_maskValue= m_pData->getParameters()->arrMaskNeValue;
    m_aMu = m_pData->getParameters()->arrMue;
    m_charge =  - 1;

    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }

}

double solverNe::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();


    // m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    double** R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double mult=pParams->p/(pParams->T*8.314);
    double ** phi  = m_pData->getFieldPhi()->arr;
    double dx=m_pData->getDx();
    double dy=m_pData->getDy();
    double dt=m_pData->getDt();
    double mu=1.0;
    for (int i = 1; i < m_field->cellsX-1; ++i)
    {
        for (int j = 1; j < m_field->cellsY-1; ++j)
        {

            int mm=(m_field->arr[i][j]>1e-5);
            m_aRHS[i][j]=// mult*R1_Ar_e[i][j] //*m_field->arr[i][j]+
                    +m_field->arrPrev[i][j]/dt
                    + mm*pParams->arrMue[i][j]*(pParams->arrEx[i][j] * simulationTools::ddxUpWind(m_field->arr, m_field->cellsX, dx, i, j, pParams->arrEx[i][j])
                                                + pParams->arrEy[i][j] * simulationTools::ddyUpWind(m_field->arr, m_field->cellsY, dy, i, j, pParams->arrEy[i][j])
                                                + m_field->arr[i][j] * pParams->arrLaplPhi[i][j])
                    + simulationTools::ddxCentral(m_field->arr,m_field->cellsX, dx, i, j) * simulationTools::ddxCentral(pParams->arrDe,m_field->cellsX, dx, i, j)
                    + simulationTools::ddyCentral(m_field->arr,m_field->cellsY, dy, i, j) * simulationTools::ddyCentral(pParams->arrDe,m_field->cellsY, dy, i, j);
        }
    }
    return 1;
}

double solverNe::getRhsAt(int i, int j)
{


    double ne=m_pData->getFieldNe()->arrPrev[i][j];
    double nars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev[i][j];
    double Te=m_pData->getParameters()->arrTe[i][j];
    static double n_n=m_pData->getParameters()->N;

    static double k[7];


    //for (int i=0;i<7;i++)
    {
        k[3]=m_pData->m_reactions[3]->getRate(fmax(fmin(Te,40),0.0001));
        k[4]=m_pData->m_reactions[4]->getRate(fmax(fmin(Te,40),0.0001));
        k[5]=m_pData->m_reactions[5]->getRate(fmax(fmin(Te,40),0.0001));

        //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
    }

    return ne*(k[3]*n_n + k[4]*nars) + k[5]*nars*nars;
}

void solverNe::setBc()
{
    for (int j = 0; j < m_field->cellsY; ++j)
    {
        m_field->arr[0][j] =0.0;// m_field->arr[m_field->cellsX-2][j];
        m_field->arr[m_field->cellsX-1][j] =0.0;// m_field->arr[1][j];
        m_field->arrPrev[0][j] =0.0;// m_field->arrPrev[m_field->cellsX-2][j];
        m_field->arrPrev[m_field->cellsX-1][j] =0.0;// m_field->arrPrev[1][j];
    }
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        m_field->arr[i][0] =0.0;// m_field->arr[i][m_field->cellsY-2];
        m_field->arr[i][m_field->cellsY-1] = 0.0;//m_field->arr[i][1];
        m_field->arrPrev[i][0] = 0.0;//= m_field->arrPrev[i][m_field->cellsY-2];
        m_field->arrPrev[i][m_field->cellsY-1]= 0.0;// = m_field->arrPrev[i][1];
    }


    //double electronFluxX = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEx[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddxCentral(pNe->arr, m_field->cellsX, dx, i, j));
    //double electronFluxY = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEy[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddyCentral(pNe->arr, m_field->cellsY, dy, i, j));


    //pParams->arrMue[i][j] * pParams->arrEx[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddxCentral(pNe->arr, m_field->cellsX, dx, i, j);

    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double dx=m_pData->getDx();
    double dy=m_pData->getDy();


    simulationData::boundary* b= &(this->m_pData->getParameters()->bound);
    for (int n=0; n<b->num; n++)
    {
        int ii=b->arrI[n];
        int jj=b->arrJ[n];

        int iip=b->arrI_inter[n];
        int jjp=b->arrJ_inter[n];

        //m_maskValue[ii][jj]=m_field->arr[iip][jjp];


        double d=pParams->arrDe[iip][jjp];
        double mu=pParams->arrMue[iip][jjp];
        double Ex=pParams->arrEx[iip][jjp];
        double Ey=pParams->arrEy[iip][jjp];



        double dxy=sqrt(dx*dx*(iip-ii)*(iip-ii) + dy*dy*(jjp-jj)*(jjp-jj));
        double nx=dx*(iip-ii)/dxy;
        double ny=dy*(jjp-jj)/dxy;

        double En=-(Ex*nx+Ey*ny);
        double ep=m_field->arr[iip][jjp];
        double ebound=0;

        //pParams->arrMue[iip][jjp] * pParams->arrEx[iip][jjp] * pNe->arr[iip][jjp] - pParams->arrDe[iip][jjp]*(pNe->arr[iip][jjp]-pNe->arr[ii][jj])/dxy=0;
        //mu*En*ep - d*(ep-ebound)/dxy=0;
        ebound = ep - mu*En*ep*dxy/d;
        // m_maskValue[ii][jj]=m_maskValue[ii][jj]*0.98+0.02*ebound;//m_field->arr[iip][jjp];

    }
}

solverEnergy::~solverEnergy()
{

}

solverEnergy::solverEnergy(simulationData* pData):simulationSolver(pData)
{
    m_pData = pData;

    m_field = m_pData->getFieldEnergy();
    m_aNu = m_pData->getParameters()->arrDeps;
    m_aRHS = new double*[m_field ->cellsX];
    m_mask = m_pData->getParameters()->arrMaskEnergy;
    m_maskValue= m_pData->getParameters()->arrMaskEnergyValue;
     m_charge = - 1;
    m_aMu = m_pData->getParameters()->arrMueps;
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

double solverEnergy::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField2d* pNe = m_pData->getFieldNe();
    //double dz = m_pData->getDz();


    // m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);

    double** R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double de=m_pData->getReactionDe(simulationData::ReactionName::comsol_eAr_2eArp);

    double mult=fabs(m_pData->q*de*pParams->p/(pParams->T*8.314));
    /*for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*simulationTools::ddzCentral(pNe->arr, m_field ->cellsNumber, dz, i);
        //(pNe->arr[i+1]-pNe->arr[i]);
        m_aRHS[i]= mult*R1_Ar_e[i]*pNe->arr[i]
                   + pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                   + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i)
                   - electronFlux * pParams->arrE[i];
    }*/
    double dx=m_pData->getDx();
    double dy=m_pData->getDy();
    double dt=m_pData->getDt();
    double mu=1.0;
    for (int i = 1; i < m_field->cellsX-1; ++i)
    {
        for (int j = 1; j < m_field->cellsY-1; ++j)
        {

            int mm=(m_field->arr[i][j]>1e-5);
            int mm2=(pNe->arr[i][j]>1e-5);
            double electronFluxX = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEx[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddxCentral(pNe->arr, m_field->cellsX, dx, i, j));
            double electronFluxY = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEy[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddyCentral(pNe->arr, m_field->cellsY, dy, i, j));
            m_aRHS[i][j]=// mult*R1_Ar_e[i][j]//*pNe->arr[i][j]+
                   +m_field->arrPrev[i][j]/dt
                    + mm*pParams->arrMueps[i][j]*(pParams->arrEx[i][j] * simulationTools::ddxUpWind(m_field->arr, m_field->cellsX, dx, i, j, pParams->arrEx[i][j])
                                                  + pParams->arrEy[i][j] * simulationTools::ddyUpWind(m_field->arr, m_field->cellsY, dy, i, j, pParams->arrEy[i][j])
                                                  + m_field->arr[i][j] * pParams->arrLaplPhi[i][j])
                    + simulationTools::ddxCentral(m_field->arr,m_field->cellsX, dx, i, j) * simulationTools::ddxCentral(pParams->arrDeps,m_field->cellsX, dx, i, j)
                    + simulationTools::ddyCentral(m_field->arr,m_field->cellsY, dy, i, j) * simulationTools::ddyCentral(pParams->arrDeps,m_field->cellsY, dy, i, j)
                    - electronFluxX * pParams->arrEx[i][j] - electronFluxY * pParams->arrEy[i][j];

        }
    }
    return 1;
}

double solverEnergy::getRhsAt(int i, int j)
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField2d* pNe = m_pData->getFieldNe();

    double dx=m_pData->getDx();
    double dy=m_pData->getDy();
    double dt=m_pData->getDt();
    double mu=1.0;


    int mm2=(pNe->arr[i][j]>1e-5);
    double electronFluxX = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEx[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddxCentral(pNe->arr, m_field->cellsX, dx, i, j));
    double electronFluxY = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEy[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddyCentral(pNe->arr, m_field->cellsY, dy, i, j));

    double ne=m_pData->getFieldNe()->arrPrev[i][j];
    double nars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev[i][j];
    double Te=m_pData->getParameters()->arrTe[i][j];
    static double n_n=m_pData->getParameters()->N;

    static double k[7],dE[7];


    for (int i=1;i<6;i++)
    {
        k[i]=m_pData->m_reactions[i]->getRate(fmax(fmin(Te,40),0.0001));
        dE[i]=m_pData->m_reactions[i]->getDe();
    }


    return   - electronFluxX * pParams->arrEx[i][j] - electronFluxY * pParams->arrEy[i][j] +
            ne*(n_n*(dE[1]*k[1] + dE[3]*k[3])
            +nars*(dE[2]*k[2] + dE[4]*k[4])) + dE[5]*k[5]*nars*nars;

}

void solverEnergy::setBc()
{
    for (int j = 0; j < m_field->cellsY; ++j)
    {
        m_field->arr[0][j] =0.0;// m_field->arr[m_field->cellsX-2][j];
        m_field->arr[m_field->cellsX-1][j] =0.0;// m_field->arr[1][j];
        m_field->arrPrev[0][j] =0.0;// m_field->arrPrev[m_field->cellsX-2][j];
        m_field->arrPrev[m_field->cellsX-1][j] =0.0;// m_field->arrPrev[1][j];
    }
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        m_field->arr[i][0] =0.0;// m_field->arr[i][m_field->cellsY-2];
        m_field->arr[i][m_field->cellsY-1] = 0.0;//m_field->arr[i][1];
        m_field->arrPrev[i][0] = 0.0;//= m_field->arrPrev[i][m_field->cellsY-2];
        m_field->arrPrev[i][m_field->cellsY-1]= 0.0;// = m_field->arrPrev[i][1];
    }

    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double dx=m_pData->getDx();
    double dy=m_pData->getDy();

    simulationData::boundary* b= &(this->m_pData->getParameters()->bound);
    for (int n=0; n<b->num; n++)
    {
        int ii=b->arrI[n];
        int jj=b->arrJ[n];

        int iip=b->arrI_inter[n];
        int jjp=b->arrJ_inter[n];

        //    m_maskValue[ii][jj]=m_field->arr[iip][jjp];


        double d=pParams->arrDeps[iip][jjp];
        double mu=pParams->arrMueps[iip][jjp];
        double Ex=pParams->arrEx[iip][jjp];
        double Ey=pParams->arrEy[iip][jjp];



        double dxy=sqrt(dx*dx*(iip-ii)*(iip-ii) + dy*dy*(jjp-jj)*(jjp-jj));
        double nx=dx*(iip-ii)/dxy;
        double ny=dy*(jjp-jj)/dxy;

        double En=-(Ex*nx+Ey*ny);
        double ep=m_field->arr[iip][jjp];
        double ebound=0;

        //pParams->arrMue[iip][jjp] * pParams->arrEx[iip][jjp] * pNe->arr[iip][jjp] - pParams->arrDe[iip][jjp]*(pNe->arr[iip][jjp]-pNe->arr[ii][jj])/dxy=0;
        //mu*En*ep - d*(ep-ebound)/dxy=0;
        ebound = ep - mu*En*ep*dxy/d;
        //  m_maskValue[ii][jj]=m_maskValue[ii][jj]*0.98+0.02*ebound;//m_field->arr[iip][jjp];


    }
}

solverPhi::~solverPhi()
{

}

solverPhi::solverPhi(simulationData* pData):simulationSolver(pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldPhi();
    m_aRHS = new double*[m_field ->cellsX];
    for (int i = 0; i < m_field ->cellsX; ++i)
    {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

double solverPhi::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double q_e= 1.6*10e-19;
    double eps_0=8.85*10e-12;
    double q_over_eps0=q_e/eps_0;

    double** pNe = m_pData->getFieldNe()->arr;
    double** pArp = m_pData->getFieldHeavySpicies(0)->arr;

    double mult= (pParams->p*6.022e23)/(pParams->T*8.314);
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();

    for (int i=1; i<m_field ->cellsX-1; i++)
    {
        for (int j=1; j<m_field ->cellsY-1; j++)
        {
            m_aRHS[i][j]= /*q_over_eps0*0.5*((pNe[i][j]+pNe[i+1][j])+mult*(pArp[i][j]+pArp[i+1][j]))*/
                    -q_over_eps0*((pNe[i][j])-(pArp[i][j]))
                    - ((m_field->arrPrev[i+1][j]-m_field->arrPrev[i-1][j])*(pParams->arrEps[i+1][j]-pParams->arrEps[i-1][j])/(4*dx*dx)
                    +(m_field->arr[i][j+1]-m_field->arr[i][j])*(pParams->arrEps[i][j+1]-pParams->arrEps[i][j])/(dy*dy))/pParams->arrEps[i][j];
            pParams->arrLaplPhi[i][j] = m_aRHS[i][j];
        }
    }
    return 1;
}

void solverPhi::setBc()
{

}

double solverPhi::solve(int iNumberIteration)
{

    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();

    for (int i = 0; i < iNumberIteration; ++i)
    {

        for (int i=0;i<m_field ->cellsX;i++)
        {
            m_field->arr[i][m_field ->cellsY-1]=0.0;
            m_field->arr[i][0]=0.0;
            //if ((i<20)&&(i>0)) m_field->arr[i][m_field ->cellsY-1]=-1500.0;
            if ((i>=20)&&(i<m_field ->cellsX-1)) m_field->arr[i][0]=+1500;
        }

        getRhs();

        for (int nnn=0;nnn<13;nnn++)
        {
            for (int i=1;i<m_field ->cellsX-1;i++)
            {
                m_field->arr[i][m_field ->cellsY-1]=(m_field->arr[i][m_field ->cellsY-1]+m_field->arr[i-1][m_field ->cellsY-1]+m_field->arr[i+1][m_field ->cellsY-1])/3.0;//-A;//.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
                m_field->arr[i][0]=(m_field->arr[i-1][0]+m_field->arr[i][0]+m_field->arr[i+1][0])/3.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
            }
        }

        INPUT_PARAM par;
        par.a=((2.0)/(dx*dx)+2.0/(dy*dy));
        par.bp=-1.0/(dx*dx);
        par.bm=-1.0/(dx*dx);
        par.cp=-1.0/(dy*dy);
        par.cm=-1.0/(dy*dy);

        par.w_bc_type=1;
        par.e_bc_type=1;
        par.n_bc_type=1;
        par.s_bc_type=3;

        par.w_bc_val=0.0;
        par.e_bc_val=0.0;
        par.n_bc_val=0.0;
        par.s_bc_val=0.0;
        par.dx = dx;
        par.dy = dy;

        for (int i=1;i<10;i++)
        {
            multigrid_N(par, m_field->arr, m_aRHS, pParams->arrMaskPhi, pParams->arrMaskPhiValue, 8, 3);
        }

    }

    return -1;
}

solverHeavySpicies::solverHeavySpicies(simulationData *pData, int num):simulationSolver(pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldHeavySpicies(num);
    m_charge = m_pData->getHeavySpiciesCharge(num);
    m_aRHS = new double*[m_field ->cellsX];
    m_mask = m_pData->getParameters()->arrMaskHeavy;

    m_aNu = m_pData->getParameters()->arrDomega;
    m_aMu = m_pData->getParameters()->arrMuomega;

   m_specie=m_field->m_specie;

    m_maskValue= m_pData->getParameters()->arrMaskHeavyValue;
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

solverHeavySpicies::~solverHeavySpicies()
{

}

double solverHeavySpicies::getRhs()
{

    simulationData::simulationField2d* pNe = m_pData->getFieldNe();
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();
    double dt = m_pData->getDt();

    //  m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    double** R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double mult=pParams->p / ( pParams->T * 8.314 * 6.022e23 * pParams->rho);


    for (int i=1; i<m_field ->cellsX-1; i++)
    {
        for (int j=1; j<m_field ->cellsY-1; j++)
        {
            m_aRHS[i][j]= //mult*R1_Ar_e[i][j]//*pNe->arr[i][j]+
                    + m_field->arrPrev[i][j]/dt +
                    -  m_charge  * pParams->arrMuomega[i][j] * (pParams->arrEx[i][j] * simulationTools::ddxCentral(m_field->arr, m_field ->cellsX, dx, i, j)
                                                                + pParams->arrEy[i][j] * simulationTools::ddyCentral(m_field->arr, m_field ->cellsY, dy, i, j)
                                                                + m_field->arr[i][j] * simulationTools::ddxCentral(pParams->arrEx, m_field ->cellsX, dx, i, j)
                                                                + m_field->arr[i][j] * simulationTools::ddyCentral(pParams->arrEy, m_field ->cellsY, dy, i, j));
        }
    }
    return 1;
}

double solverHeavySpicies::getRhsAt(int i, int j)
{


    double ne=m_pData->getFieldNe()->arrPrev[i][j];
    double nars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev[i][j];
    double Te=m_pData->getParameters()->arrTe[i][j];
    static double n_n=m_pData->getParameters()->N;

    static double k[7];

    if (m_specie==simulationData::SpecieName::Ar_plus)
    {
        //for (int i=0;i<7;i++)
        {
            k[3]=m_pData->m_reactions[3]->getRate(fmax(fmin(Te,40),0.0001));
            k[4]=m_pData->m_reactions[4]->getRate(fmax(fmin(Te,40),0.0001));
            k[5]=m_pData->m_reactions[5]->getRate(fmax(fmin(Te,40),0.0001));

            //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }

        return ne*(k[3]*n_n + k[4]*nars) + k[5]*nars*nars;

    }else if (m_specie==simulationData::SpecieName::Ar_star)
    {

        //        ne*(k[1]*n_n -(k[2]+k[4])*nars) -2.0*k[5]*nars*nars;

        //for (int i=0;i<7;i++)
        {
            k[1]=m_pData->m_reactions[1]->getRate(fmax(fmin(Te,40),0.0001));
            k[2]=m_pData->m_reactions[2]->getRate(fmax(fmin(Te,40),0.0001));
            k[4]=m_pData->m_reactions[4]->getRate(fmax(fmin(Te,40),0.0001));
            k[5]=m_pData->m_reactions[5]->getRate(fmax(fmin(Te,40),0.0001));

            //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }

        return ne*(k[1]*n_n -(k[2]+k[4])*nars) -2.0*k[5]*nars*nars;;
    }
    return 0;
}

void solverHeavySpicies::setBc()
{
    for (int j = 0; j < m_field->cellsY; ++j)
    {
        m_field->arr[0][j] =0.0;// m_field->arr[m_field->cellsX-2][j];
        m_field->arr[m_field->cellsX-1][j] =0.0;// m_field->arr[1][j];
        m_field->arrPrev[0][j] =0.0;// m_field->arrPrev[m_field->cellsX-2][j];
        m_field->arrPrev[m_field->cellsX-1][j] =0.0;// m_field->arrPrev[1][j];
    }
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        m_field->arr[i][0] =0.0;// m_field->arr[i][m_field->cellsY-2];
        m_field->arr[i][m_field->cellsY-1] = 0.0;//m_field->arr[i][1];
        m_field->arrPrev[i][0] = 0.0;//= m_field->arrPrev[i][m_field->cellsY-2];
        m_field->arrPrev[i][m_field->cellsY-1]= 0.0;// = m_field->arrPrev[i][1];
    }
    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double dx=m_pData->getDx();
    double dy=m_pData->getDy();

    simulationData::boundary* b= &(this->m_pData->getParameters()->bound);
    for (int n=0; n<b->num; n++)
    {
        int ii=b->arrI[n];
        int jj=b->arrJ[n];

        int iip=b->arrI_inter[n];
        int jjp=b->arrJ_inter[n];

        //   m_maskValue[ii][jj]=m_field->arr[iip][jjp];

        double d=pParams->arrDomega[iip][jjp];
        double mu=-pParams->arrMuomega[iip][jjp];
        double Ex=pParams->arrEx[iip][jjp];
        double Ey=pParams->arrEy[iip][jjp];



        double dxy=sqrt(dx*dx*(iip-ii)*(iip-ii) + dy*dy*(jjp-jj)*(jjp-jj));
        double nx=dx*(iip-ii)/dxy;
        double ny=dy*(jjp-jj)/dxy;

        double En=-(Ex*nx+Ey*ny);
        double ep=m_field->arr[iip][jjp];
        double ebound=0;

        //pParams->arrMue[iip][jjp] * pParams->arrEx[iip][jjp] * pNe->arr[iip][jjp] - pParams->arrDe[iip][jjp]*(pNe->arr[iip][jjp]-pNe->arr[ii][jj])/dxy=0;
        //mu*En*ep - d*(ep-ebound)/dxy=0;
        ebound = ep - mu*En*ep*dxy/d;
        // m_maskValue[ii][jj]=m_maskValue[ii][jj]*0.98+0.02*ebound;//m_field->arr[iip][jjp];


    }
}


