#include "crosssection.h"
#include "math.h"
#include <cstddef>
#include <stdio.h>
#include <vector>
#include <QDebug>


splineInterp::splineInterp(int splN)
{
    if (splN > 0)
    {
        m_splNum = splN;
        m_splX = new double[m_splNum+1];
        m_splA = new double[m_splNum+1];
        m_splB = new double[m_splNum+1];
        m_splC = new double[m_splNum+1];
        m_splD = new double[m_splNum+1];
    }else
    {
        m_splNum = 0;
        m_splX = nullptr;
        m_splA = nullptr;
        m_splB = nullptr;
        m_splC = nullptr;
        m_splD = nullptr;
    }
}



void splineInterp::buildSpline()
{

    int n=m_splNum;

    double  h[n], A[n], l[n + 1],u[n + 1], z[n + 1];

    /** Step 1 */
    int i=0;
    int j=0;
    for ( i = 0; i <= n - 1; ++i) h[i] = m_splX[i + 1] - m_splX[i];


    /** Step 2 */
    for ( i = 1; i <= n - 1; ++i)
        A[i] = 3 * (m_splA[i + 1] - m_splA[i]) / h[i] - 3 * (m_splA[i] - m_splA[i - 1]) / h[i - 1];

    /** Step 3 */
    l[0] = 1;
    u[0] = 0;
    z[0] = 0;

    /** Step 4 */
    for (i = 1; i <= n - 1; ++i) {
        l[i] = 2 * (m_splX[i + 1] - m_splX[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    /** Step 5 */
    l[n] = 1;
    z[n] = 0;
    m_splC[n] = 0;

    /** Step 6 */
    for (j = n - 1; j >= 0; --j) {
        m_splC[j] = z[j] - u[j] * m_splC[j + 1];
        m_splB[j] = (m_splA[j + 1] - m_splA[j]) / h[j] - h[j] * (m_splC[j + 1] + 2 * m_splC[j]) / 3;
        m_splD[j] = (m_splC[j + 1] - m_splC[j]) / (3 * h[j]);
    }

}

void splineInterp::fillData(double esNew[][2], int n,double Na)//filling the data array into a spline
{

    double x0,x1;
    x0=esNew[0][0];
    m_x0=x0;
    x1=esNew[n-1][0];
    qDebug()<<"x0 = "<<x0<<" x1= "<<x1;


    double dx2=pow((x1-x0),m_xPow)/m_splNum;
    m_dXp=dx2;

    m_splX[0]= x0;
    m_splA[0]=esNew[0][1]/Na;
    m_splB[0]=0.0;
    m_splC[0]=0.0;
    m_splD[0]=0.0;

    for (int i=1;i<=m_splNum;i++)
    {
        m_splX[i]= x0+pow(dx2*i,1.0/m_xPow);
        int j=0;
        while  ((m_splX[i]>esNew[j][0])&&(j<n)) j++;
        j--;
        if (j<n-1)
        {
            double alpha=(m_splX[i]-esNew[j][0])/(esNew[j+1][0]-esNew[j][0]);
            qDebug()<<"alpha "<<alpha<<" j= "<<j;
            m_splA[i]=(esNew[j][1]*(1.0-alpha)+esNew[j+1][1]*alpha)/Na;
        }else
            m_splA[i]=esNew[n-1][1]/Na;

        m_splB[i]=0.0;
        m_splC[i]=0.0;
        m_splD[i]=0.0;
        qDebug()<<"i = "<<i<<" xi = "<<m_splX[i]<<" yi ="<<m_splA[i];
    }


    buildSpline();
}

double splineInterp::getSpline(double xp)
{
    int i=0;

   /* while ((m_splX[i]<xp)&&(i<m_splNum))
    {
        i++;
    }



    i--;*/
   // double dxp=pow((m_splX[0]-m_splX[m_splNum]),m_xPow)/m_splNum;

    //pow(m_splX[i]-x0,m_xPow)/dx2= i;
    i=(int)(pow(xp-m_x0,m_xPow)/m_dXp);

    if (i>m_splNum)
        return m_splA[m_splNum];

    double dx=xp-m_splX[i];
    return m_splA[i]+m_splB[i]*dx+m_splC[i]*dx*dx+m_splD[i]*dx*dx*dx;

}

double splineInterp::getPrime(double xp)
{
    int i=(int)(pow(xp-m_x0,m_xPow)/m_dXp);

    if (i>m_splNum)
        return m_splA[m_splNum];

    double dx=xp-m_splX[i];
    return m_splB[i]+2.0*m_splC[i]*dx+3.0*m_splD[i]*dx*dx;
}

splineInterp::~splineInterp()
{
    if (m_splX!=nullptr)
    {
        delete[] m_splX;
        delete[] m_splA;
        delete[] m_splB;
        delete[] m_splC;
        delete[] m_splD;
    }
}



crossSection::crossSection()
{
    m_sigmas = nullptr;
    m_energs = nullptr;
    m_K = nullptr;
    m_phis = nullptr;
    m_splX = nullptr;
    m_splA = nullptr;
    m_splB = nullptr;
    m_splC = nullptr;
    m_splD = nullptr;
    m_kNum = 0;
    m_sigNum = 0;
    m_splNum = 0;
}

crossSection::crossSection(int m_KN, int sigN, int splN)
{
    m_sigmas = nullptr;
    m_energs = nullptr;
    m_K = nullptr;
    m_phis = nullptr;
    m_splX = nullptr;
    m_splA = nullptr;
    m_splB = nullptr;
    m_splC = nullptr;
    m_splD = nullptr;
    m_kNum = 0;
    m_sigNum = 0;
    m_splNum = 0;


    init(m_KN,sigN,splN);

}

void crossSection::init(int m_KN, int sigN, int splN)
{
    m_kNum = m_KN;
    m_sigNum = sigN;
    m_splNum = splN;

    if (m_sigmas != nullptr)
    {
        delete[] m_sigmas;
        delete[] m_energs;
        delete[] m_K;
        delete[] m_phis;
        delete[] m_splX;
        delete[] m_splA;
        delete[] m_splB;
        delete[] m_splC;
        delete[] m_splD;
    }
    m_sigmas = new double[m_sigNum];
    m_energs = new double[m_sigNum];
    m_K = new double[m_kNum];
    m_phis = new double[m_kNum];
    m_splX = new double[m_splNum+1];
    m_splA = new double[m_splNum+1];
    m_splB = new double[m_splNum+1];
    m_splC = new double[m_splNum+1];
    m_splD = new double[m_splNum+1];
}

crossSection::~crossSection()
{
    if (m_sigmas!=nullptr)
    {
        delete[] m_sigmas;
        delete[] m_energs;
        delete[] m_K;
        delete[] m_phis;
        delete[] m_splX;
        delete[] m_splA;
        delete[] m_splB;
        delete[] m_splC;
        delete[] m_splD;
    }
}

void crossSection::fillSigmas(double *sNew, double *eNew, int n)
{
    m_sigNum=n;
    if (m_sigmas!=nullptr)
    {
        delete[] m_sigmas;
        delete[] m_energs;
    }

    m_sigmas = new double[m_sigNum];
    m_energs = new double[m_sigNum];

    for (int i=0;i<m_sigNum;i++)
    {
        m_sigmas[i]=sNew[i];
        m_energs[i]=eNew[i];
    }
    fillK();
}


void crossSection::fillSigmas2(double esNew[][2], int n)
{
    m_sigNum=n;
    if (m_sigmas!=nullptr)
    {
        delete[] m_sigmas;
        delete[] m_energs;
    }

    m_sigmas = new double[m_sigNum];
    m_energs = new double[m_sigNum];


    for (int i=0;i<m_sigNum;i++)
    {
        m_sigmas[i]=esNew[i][1];
        m_energs[i]=esNew[i][0];

    }

    fillK();
}

void crossSection::readFromFile(char *fileName)
{
    FILE* file=fopen(fileName,"r");

    if (file == NULL)
        return;

    char str[4096];
    fgets(str,4096,file);
    fgets(str,4096,file);
    fgets(str,4096,file);
    fgets(str,4096,file);
    fgets(str,4096,file);
    std::vector<double> ens;
    std::vector<double> sgs;
    while (fgets(str,4096,file) != NULL)
    {
        double ee,sg;
        int ret=sscanf(str,"%lf %lf",&ee,&sg);
        if (ret == 2)
        {
            ens.push_back(ee);
            sgs.push_back(sg);
        }
    }

    if (ens.size() == 0)
        return;

    m_sigNum=ens.size();
    if (m_sigmas != nullptr)
    {
        delete[] m_sigmas;
        delete[] m_energs;
    }

    m_sigmas = new double[m_sigNum];
    m_energs = new double[m_sigNum];


    for (int i=0;i<m_sigNum;i++)
    {
        m_sigmas[i]=sgs[i];
        m_energs[i]=ens[i];
    }

    fillK();
}

double crossSection::getKICoarse(double phi)
{
    double beta1=2.07296489683;
    double gamma=5.9312e5;
    double beta2=1.5;

    double A=beta2/fabs(phi);
    double sum=0.0;
    for (int i=1;i<m_sigNum;i++)
    {
        double em,ep;
        em=m_energs[i-1];
        ep=m_energs[i];

        double sm,sp;
        sm=m_sigmas[i-1];
        sp=m_sigmas[i];

        double fm,fp;
        fm=em*exp(-em*A);
        fp=ep*exp(-ep*A);

        sum+=(ep-em)*(fm*sm+fp*sp)*0.5;
    }
    return gamma*beta1*pow(phi,-1.5)*sum;


}

double crossSection::getKIFine(double phi)
{
    double beta1=2.07296489683;
    double gamma=5.9312e5;

    double sum=0.0;
    for (int i=1;i<m_sigNum;i++)
    {
        double em,ep;
        em=m_energs[i-1];
        ep=m_energs[i];

        double sm,sp;
        sm=m_sigmas[i-1];
        sp=m_sigmas[i];

        sum+=integrate(sm,sp,em,ep,phi);
    }
    return gamma*beta1*pow(phi,-1.5)*sum;
}

double crossSection::integrate(double sm, double sp, double em, double ep, double phi)
{
    const double beta2=1.5;
    double de=phi/(30.0*beta2);
    double A=beta2/fabs(phi);
    double sum=0.0;

    if (em>30*phi/(beta2))
    {
        double fm,fp;
        fm=em*exp(-em*A);
        fp=ep*exp(-ep*A);
        sum=(ep-em)*(fm*sm+fp*sp)*0.5;
    }
    else
    {
        int n_sub=(int)((ep-em)/de)+1;
        for (int n=0;n<n_sub;n++)
        {
            double e_curr=em+n*(ep-em)/n_sub;
            double e_curr_p=em+(n+1)*(ep-em)/n_sub;

            double s_curr=sm+n*(sp-sm)/n_sub;
            double s_curr_p=sm+(n+1)*(sp-sm)/n_sub;

            double f_curr=e_curr*exp(-e_curr*A);
            double f_curr_p=e_curr_p*exp(-e_curr_p*A);
            sum+=(e_curr_p-e_curr)*(f_curr*s_curr+f_curr_p*s_curr_p)*0.5;
        }
    }
    return sum;
}

void crossSection::fillK()
{

    double phi_max=100.0; //eV
    double dphi=phi_max/m_kNum;
    for (int i=0;i<m_kNum;i++)
    {
        m_phis[i]=i*dphi+0.0001;
        m_K[i]=getKIFine(m_phis[i]);
    }

    for (int i=0;i<=m_splNum;i++)
    {
        m_splX[i]= m_phis[i*((m_kNum-1)/m_splNum)];
        m_splA[i]= m_K[i*((m_kNum-1)/m_splNum)];
        m_splB[i]=0.0;
        m_splC[i]=0.0;
        m_splD[i]=0.0;
    }

    buildSpline();

}
