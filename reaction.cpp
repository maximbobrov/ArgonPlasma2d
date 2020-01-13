#include "reaction.h"
#include "simulationdata.h"
#include <QDebug>

static const int comsol_EAr_2EArp_num =101;

static double comsol_EAr_2EArp[comsol_EAr_2EArp_num][2] = {
    #include "comsol/e+Ar_2e+Ar+.txt"
};

static const int comsol_EAr_EAr_num =201;

static double comsol_EAr_EAr[comsol_EAr_EAr_num][2] = {
    #include "comsol/e+Ar_e+Ar.txt"
};

static const int comsol_EAr_EArs_num =201;

static double comsol_EAr_EArs[comsol_EAr_EArs_num][2] = {
    #include "comsol/e+Ar_e+Ars.txt"
};

static const int comsol_EArs_EAr_num =201;

static double comsol_EArs_EAr[comsol_EArs_EAr_num][2] = {
    #include "comsol/e+Ars_e+Ar.txt"
};

static const int comsol_EArs_2EArp_num =201;

static double comsol_EArs_2EArp[comsol_EArs_2EArp_num][2] = {
    #include "comsol/e+Ars_2e+Ar+.txt"
};



reaction::reaction(simulationData* data)
{
    m_pData=data;
    if (data != nullptr)
    {
        int cx=data->getCellsXNumber();
        int cy=data->getCellsYNumber();
        m_R=new double*[cx];
        for (int i = 0 ; i < cx ; i++)
        {
            m_R[i] = new double[cy];
            for (int j = 0; j< cy; j++)
                m_R[i][j] = 0;
        }
    }
}

void reaction::calc()
{
    //fill from spline
    //rate only depends on electron temperature

    int cx=m_pData->getCellsXNumber();
    int cy=m_pData->getCellsYNumber();
    for (int i = 0 ; i < cx ; i++)
    {
        for (int j = 0; j< cy; j++)
          m_R[i][j] = m_spline->getSpline(m_pData->getArrTe()[i][j]);//i*100.0/m_pData->getCellsNumber());
        //m_cs->getSpline(En[i])*N*Ar[i]*Ne[i]*0.0;

        // qDebug()<<"i="<<i<<" R="<<m_R[i];
    }
}

double reaction::getDe()
{
    return 0.0;
}

double **reaction::getR()
{
    return m_R;
}

double reaction::getRate(double eps)
{
    return m_spline->getSpline(eps);
}

double reaction::getDeriv(double eps)
{
    return m_spline->getPrime(eps);
}

reactionEAr_EAr_comsol::reactionEAr_EAr_comsol(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        /*     m_cs=new crossSection(36,g_EAr_EAr_num,12);
        m_cs->fillSigmas2(g_EAr_EAr,g_EAr_EAr_num); //fill from static arrays (safer then messing up with external txt files)*/

        m_spline=new splineInterp(20);

        double Na=6.022e23;
        m_spline->fillData(comsol_EAr_EAr,comsol_EAr_EAr_num,Na);

    }
}



double reactionEAr_EAr_comsol::getDe()
{
    return m_massRatio;//m_energy;
}

reactionEAr_EArs_comsol::reactionEAr_EArs_comsol(simulationData *data):reaction(data)
{
    /* if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EAr_EArs_num,12);
        m_cs->fillSigmas2(g_EAr_EArs,g_EAr_EArs_num);
    }*/

    if (data != nullptr)
    {
        m_spline=new splineInterp(20);

        double Na=6.022e23;
        m_spline->fillData(comsol_EAr_EArs,comsol_EAr_EArs_num,Na);
        //m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}


double reactionEAr_EArs_comsol::getDe()
{
    return m_energy;
}

reactionEArs_EAr_comsol::reactionEArs_EAr_comsol(simulationData *data):reaction(data)
{
    /* if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EAr_EArs_num,12);
        m_cs->fillSigmas2(g_EAr_EArs,g_EAr_EArs_num);
    }*/

    if (data != nullptr)
    {
        m_spline=new splineInterp(20);

        double Na=6.022e23;
        m_spline->fillData(comsol_EArs_EAr,comsol_EArs_EAr_num,Na);
        //m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}


double reactionEArs_EAr_comsol::getDe()
{
    return m_energy;
}



reactionEArs_2EArp_comsol::reactionEArs_2EArp_comsol(simulationData *data):reaction(data)
{
    /* if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EArs_2EArp_num,12);
        m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }*/

    if (data != nullptr)
    {
        m_spline=new splineInterp(20);

        double Na=6.022e23;
        m_spline->fillData(comsol_EArs_2EArp,comsol_EArs_2EArp_num,Na);
        //m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}


double reactionEArs_2EArp_comsol::getDe()
{
    return m_energy;
}

///////////////////////////////good below
reactionEAr_2EArp_comsol::reactionEAr_2EArp_comsol(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        m_spline=new splineInterp(20);

        double Na=6.022e23;
        m_spline->fillData(comsol_EAr_2EArp,comsol_EAr_2EArp_num,Na);
        //m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}


double reactionEAr_2EArp_comsol::getDe()
{
    return m_energy;
}

reactionArsArs_EArArp_comsol::reactionArsArs_EArArp_comsol(simulationData *data):reaction(data)
{

}

void reactionArsArs_EArArp_comsol::calc()
{
    int cx=m_pData->getCellsXNumber();
    int cy=m_pData->getCellsYNumber();
    for (int i = 0 ; i < cx ; i++)
    {
        for (int j = 0; j< cy; j++)
              m_R[i][j] = m_k;
    }
}

double reactionArsArs_EArArp_comsol::getRate(double eps)
{
    return m_k/6.022e23;
}

double reactionArsArs_EArArp_comsol::getDeriv(double eps)
{
    return 0.0;
}

reactionArsAr_ArAr_comsol::reactionArsAr_ArAr_comsol(simulationData *data):reaction(data)
{

}

void reactionArsAr_ArAr_comsol::calc()
{
    int cx=m_pData->getCellsXNumber();
    int cy=m_pData->getCellsYNumber();
    for (int i = 0 ; i < cx ; i++)
    {
        for (int j = 0; j< cy; j++)
              m_R[i][j] = m_k;
    }
}

double reactionArsAr_ArAr_comsol::getRate(double eps)
{
    return m_k/6.022e23;
}

double reactionArsAr_ArAr_comsol::getDeriv(double eps)
{
    return 0.0;
}

