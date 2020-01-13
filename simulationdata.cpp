#include "simulationdata.h"
#include "simulationtools.h"
#include "string.h"
#include "reaction.h"
#include <QDebug>
#include <math.h>

simulationData::simulationData(int iCellsX,int iCellsY)
{
    m_cellsX = iCellsX;
    m_cellsY = iCellsY;
    /*m_fieldNe = new simulationField2d(iCellsX, iCellsY,"electrons",simulationData::SpecieName::e);
    m_fieldEnergy = new simulationField2d(iCellsX, iCellsY,"energy",simulationData::SpecieName::En);
    m_fieldPhi = new simulationField2d(iCellsX, iCellsY,"potential",simulationData::SpecieName::phi);
    m_fieldsHeavySpecies.push_back(new simulationField2d(iCellsX, iCellsY,"Ar+",simulationData::SpecieName::Ar_plus));
    m_chargeHeavySpecies.push_back(1);
    m_fieldsHeavySpecies.push_back(new simulationField2d(iCellsX, iCellsY,"Ars",simulationData::SpecieName::Ar_star));
    m_chargeHeavySpecies.push_back(0);*/
    m_fieldNe = new simulationField2d(iCellsX, iCellsY,"electrons",simulationData::SpecieName::e);
    m_fieldEnergy = new simulationField2d(iCellsX, iCellsY,"energy",simulationData::SpecieName::En);
    m_fieldPhi = new simulationField2d(iCellsX, iCellsY,"potential",simulationData::SpecieName::phi);
    m_fieldsHeavySpecies.push_back(new simulationField2d(iCellsX, iCellsY,"Ar+",simulationData::SpecieName::Ar_plus));
    m_chargeHeavySpecies.push_back(1);
    m_fieldsHeavySpecies.push_back(new simulationField2d(iCellsX, iCellsY,"Ars",simulationData::SpecieName::Ar_star));
    m_chargeHeavySpecies.push_back(0);
    //  m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar+"));
    // m_chargeHeavySpecies.push_back(1);
    m_numberHeavySpicies = m_fieldsHeavySpecies.size();
    m_params = new simulationData::simulationParameters(iCellsX,iCellsY);

    m_reactions.push_back(new reactionEAr_EAr_comsol(this));
    m_reactions.push_back(new reactionEAr_EArs_comsol(this));
    m_reactions.push_back(new reactionEArs_EAr_comsol(this));
    m_reactions.push_back(new reactionEAr_2EArp_comsol(this));
    m_reactions.push_back(new reactionEArs_2EArp_comsol(this));
   m_reactions.push_back(new reactionArsArs_EArArp_comsol(this));
    m_reactions.push_back(new reactionArsAr_ArAr_comsol(this));
}

simulationData::~simulationData()
{

}

void simulationData::setDt(double idt)
{
    m_dt = idt;
}

void simulationData::setDx(double dx)
{
    m_dx = dx;
}

void simulationData::setDy(double dy)
{
    m_dy = dy;
}

void simulationData::setCellsXNumber(int iCellsXNumber)
{
    m_cellsX = iCellsXNumber;
}

void simulationData::setCellsYNumber(int iCellsYNumber)
{
    m_cellsY = iCellsYNumber;
}

double simulationData::getDt()
{
    return m_dt;
}

double simulationData::getDx()
{
    return m_dx;
}

double simulationData::getDy()
{
    return m_dy;
}

simulationData::simulationField2d* simulationData::getFieldNe()
{
    return m_fieldNe;
}

double **simulationData::getArrTe()
{
    return m_params->arrTe;
}


simulationData::simulationField2d* simulationData::getFieldHeavySpicies(int num)
{
    return m_fieldsHeavySpecies[num];
}

int simulationData::getHeavySpiciesCharge(int num)
{
    return m_chargeHeavySpecies[num];
}

int simulationData::getNumberHeavySpicies()
{
    return m_numberHeavySpicies;
}

double simulationData::getN()
{
    // pv=nkT
    // n/v=p/kT;
    double pres=101505;
    double T=300;
    double k=1.38e-23;


    return pres/(k*T);
}

double** simulationData::getReactionRate(simulationData::ReactionName reactName)
{
    return m_reactions[reactName]->getR();
}

double simulationData::getReactionDe(simulationData::ReactionName reactName)
{
    return m_reactions[reactName]->getDe();
}

void simulationData::calcReaction(simulationData::ReactionName reactName)
{
    m_reactions[reactName]->calc();
}

simulationData::simulationField2d* simulationData::getFieldEnergy()
{
    return m_fieldEnergy;
}

simulationData::simulationField2d* simulationData::getFieldPhi()
{
    return m_fieldPhi;
}

simulationData::simulationParameters *simulationData::getParameters()
{
    return m_params;
}

int simulationData::getCellsXNumber()
{
    return m_cellsX;
}

int simulationData::getCellsYNumber()
{
    return m_cellsY;
}

simulationData::simulationField::simulationField(int iCellsNumber, char* iName)
{
    init(iCellsNumber, iName);
}

void simulationData::simulationField::init(int iCellsNumber, char* iName)
{
    cellsNumber =  iCellsNumber;
    arr = new double[iCellsNumber];
    arrPrev = new double[iCellsNumber];
    int len=strlen(iName);
    name = new char[len];
    strcpy(name, iName);
}

simulationData::simulationField2d::simulationField2d(int iCellsX, int iCellsY, char* iName,SpecieName specie)
{

    init(iCellsX, iCellsY, iName,specie);
}

void simulationData::simulationField2d::init(int iCellsX, int iCellsY, char* iName, SpecieName specie)
{
    cellsX =  iCellsX;
    cellsY =  iCellsY;
    arr = new double*[cellsX];
    arrPrev = new double*[cellsX];
    for (int i = 0; i < cellsX; ++i) {
        arr[i] = new double [cellsY];
        arrPrev[i] = new double [cellsY];
    }
    int len=strlen(iName);
    name = new char[len+1];
    name[len]=0;
    strcpy(name, iName);
    m_specie=specie;
}

void simulationData::simulationParameters::init(int iCellsX, int iCellsY)
{
    cellsX =  iCellsX;
    cellsY =  iCellsY;
    arrDe = new double*[cellsX];
    arrDeps  = new double*[cellsX];
    arrDomega  = new double*[cellsX];
    arrMue  = new double*[cellsX];
    arrMueps  = new double*[cellsX];
    arrMuomega  = new double*[cellsX];
    arrEx = new double*[cellsX];
    arrEy = new double*[cellsX];
    arrLaplPhi = new double*[cellsX];
    arrTe  = new double*[cellsX];
    arrMaskPhi = new double*[cellsX];
    arrMaskPhiValue = new double*[cellsX];
    arrMaskNe = new double*[cellsX];
    arrMaskNeValue = new double*[cellsX];
    arrMaskEnergy = new double*[cellsX];
    arrMaskEnergyValue = new double*[cellsX];
    arrMaskHeavy = new double*[cellsX];
    arrMaskHeavyValue = new double*[cellsX];
    arrEps = new double*[cellsX];
    arrBoundMask = new int*[cellsX];

    for (int i = 0; i < cellsX; ++i) {
        arrDe[i] = new double [cellsY];
        arrDeps[i] = new double [cellsY];
        arrDomega[i] = new double [cellsY];
        arrMue[i] = new double [cellsY];
        arrMueps[i] = new double [cellsY];
        arrMuomega[i] = new double [cellsY];
        arrEx[i] = new double [cellsY];
        arrEy[i] = new double [cellsY];
        arrLaplPhi[i] = new double [cellsY];
        arrTe[i] = new double [cellsY];
        arrMaskPhi[i] = new double[cellsY];
        arrMaskPhiValue[i] = new double[cellsY];
        arrMaskNe[i] = new double[cellsY];
        arrMaskNeValue[i] = new double[cellsY];
        arrMaskEnergy[i] = new double[cellsY];
        arrMaskEnergyValue[i] = new double[cellsY];
        arrMaskHeavy[i] = new double[cellsY];
        arrMaskHeavyValue[i] = new double[cellsY];
        arrEps[i] = new double[cellsY];

        arrBoundMask[i] = new int[cellsY];
    }

    for (int i=0; i<cellsX; i++)
    {
        for (int j=0; j<cellsY; j++)
        {

            if((i >=  (cellsX-1) / 4  && i < 2*(cellsX-1)/4 && j >= (cellsY-1) / 4  && j < 2*(cellsY-1)/4))
            {
                arrMaskPhi[i][j] = 0.0;
                arrMaskPhiValue[i][j] = -1500;
            }
            else
            {
                arrMaskPhi[i][j] = 1.0;
                arrMaskPhiValue[i][j] =0.0;
            }


            if(((i >=  (cellsX-1) / 4  && i < 2*(cellsX-1)/4 && j >= (cellsY-1) / 4  && j < 2*(cellsY-1)/4)) || j <= (cellsY-1) / 4)
            {
               arrMaskNe[i][j] = 0.0;
                arrMaskNeValue[i][j] = 1e5;
                arrMaskEnergy[i][j] = 0.0;
                arrMaskEnergyValue[i][j] = 1e5;
               arrMaskHeavy[i][j] = 0.0;
                arrMaskHeavyValue[i][j] = 1e-25;

               arrBoundMask[i][j] = 0;
            }
            else
            {
                arrMaskNe[i][j] = 1.0;
                arrMaskNeValue[i][j] =0.0;
                arrMaskEnergy[i][j] = 1.0;
                arrMaskEnergyValue[i][j] =0.0;
                arrMaskHeavy[i][j] = 1.0;
                arrMaskHeavyValue[i][j] =0.0;
                arrBoundMask[i][j] = 1;
            }
        }
    }

    bound.init(cellsX,cellsY,arrBoundMask);
    T=400; //K
    p=101325/7600;//101325;//pa
    mAr=39.948/1000.0;//kg/mol
    rho=p*mAr/(8.314*T);//kg/m^3
    double Na=6.022e23; //1/mol
    N=p*Na/(8.314*T);

}

void simulationData::updateParams()
{
    simulationData::simulationParameters* pParams = m_params;
    simulationData::simulationField2d* pEn = m_fieldEnergy;
    simulationData::simulationField2d* pNe = m_fieldNe;
    simulationData::simulationField2d* pPhi = m_fieldPhi;


    for (int i=0; i<pParams->cellsX; i++)
    {
        for (int j=0; j<pParams->cellsY; j++)
        {
            if((i >= 0/* (pParams->cellsX-1) / 4 */ && i < (pParams->cellsX-1)/8 && j >= (pParams->cellsY-1) / 8  && j < 0.75*2*(pParams->cellsY-1)/8))
            {
                pParams->arrMaskPhi[i][j] = 0.0;
                pParams->arrMaskPhiValue[i][j] = -1500;
            }
            else
            {
                pParams->arrMaskPhi[i][j] = 1.0;
                pParams->arrMaskPhiValue[i][j] =0.0;
            }
            pParams->arrEps[i][j] = ( j >= (pParams->cellsY-1) / 8) ? 1.0 : 500.0;

            if(((i >=0 /* (pParams->cellsX-1) / 4*/  && i < (pParams->cellsX-1)/8 && j >= (pParams->cellsY-1) / 8  && j < 0.75*2*(pParams->cellsY-1)/8)) || j <= (pParams->cellsY-1) / 8
                    || (i == 0) || (i == pParams->cellsX - 1) || (j == 0) || (j == pParams->cellsY - 1))
            {
                pParams->arrMaskNe[i][j] = 0.0;
                pParams->arrMaskNeValue[i][j] = 0.0;//1e5;
                pParams->arrMaskEnergy[i][j] = 0.0;
                pParams->arrMaskEnergyValue[i][j] = 0.0;//1e5;
                pParams->arrMaskHeavy[i][j] = 0.0;
                pParams->arrMaskHeavyValue[i][j] = 0.0;//1e-25;

                pParams->arrBoundMask[i][j] = 0;
            }
            else
            {
                pParams->arrMaskNe[i][j] = 1.0;
                pParams->arrMaskNeValue[i][j] =0.0;
                pParams->arrMaskEnergy[i][j] = 1.0;
                pParams->arrMaskEnergyValue[i][j] =0.0;
                pParams->arrMaskHeavy[i][j] = 1.0;
                pParams->arrMaskHeavyValue[i][j] =0.0;
                pParams->arrBoundMask[i][j] = 1;
            }

            pParams->arrMue[i][j] = 4e24/m_params->N; ; //4e4; m^2/(V*s)
            pParams->arrMueps[i][j] =5.0 * pParams->arrMue[i][j] / 3.0;

            pParams->arrTe[i][j]=(2.0/3.0)*fabs(pEn->arr[i][j])/(fabs(pNe->arr[i][j])+1);
            if (pParams->arrTe[i][j]>60.0) pParams->arrTe[i][j]=60.0; //some limiter here


            pParams->arrDe[i][j] = pParams->arrMue[i][j] * pParams->arrTe[i][j];
            pParams->arrDeps[i][j] = pParams->arrMueps[i][j] * pParams->arrTe[i][j];
            pParams->arrEx[i][j] = - simulationTools::ddxCentral(pPhi->arr, pParams->cellsX, getDx(), i, j);
            pParams->arrEy[i][j] = - simulationTools::ddyCentral(pPhi->arr, pParams->cellsY, getDy(), i, j);
            pParams->arrLaplPhi[i][j] = 0.0;
            pParams->arrDomega[i][j] = 0.01;// m^2/s
            pParams->arrMuomega[i][j]=pParams->arrDomega[i][j]*q/(k_B_const*pParams->T);
        }

    }
}



simulationData::simulationParameters::simulationParameters(int iCellsX, int iCellsY)
{
    init(iCellsX,iCellsY);
}

simulationData::boundary::boundary()
{
    arrI=nullptr;
}

void simulationData::boundary::init(int cellsX, int cellsY, int **mask)
{
    if (arrI!=nullptr)
    {
        delete[] arrI;
        delete[] arrJ;
        delete[] arrI_inter;
        delete[] arrJ_inter;
    }

    std::vector<int> ii;
    std::vector<int> jj;

    std::vector<int> ii_intern;
    std::vector<int> jj_intern;


    for (int i=1; i<cellsX-1; i++)
    {
        for (int j=1; j<cellsY-1; j++)
        {
            int grad_x=mask[i+1][j]-mask[i-1][j];
            int grad_y=mask[i][j+1]-mask[i][j-1];

            if ((grad_x*grad_x+grad_y*grad_y>0.1)&&(mask[i][j]==0))
            {
                ii.push_back(i);
                jj.push_back(j);

                ii_intern.push_back(i+grad_x);
                jj_intern.push_back(j+grad_y);
            }

        }
    }

    num=ii.size();

    arrI=new int[num];
    arrJ=new int[num];
    arrI_inter=new int[num];
    arrJ_inter=new int[num];

    for (int i=0;i<num;i++)
    {
        arrI[i]=ii[i];
        arrJ[i]=jj[i];
        arrI_inter[i]=ii_intern[i];
        arrJ_inter[i]=jj_intern[i];
    }

}
