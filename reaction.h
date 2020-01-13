#ifndef REACTION_H
#define REACTION_H

//#include "simulationdata.h"
#include "crosssection.h"
class simulationData;

class reaction
{
public:
    reaction(simulationData* data);// = nullptr);
    virtual void calc();
    virtual double getDe();
    double **getR();
    virtual double getRate(double eps);
    virtual double getDeriv(double eps);
protected:
    double **m_R;
    simulationData* m_pData;
    splineInterp * m_spline;
};

class reactionEAr_2EArp_comsol:public reaction //e+Ar=>2e+Ar+  Ionization
{
public:
    reactionEAr_2EArp_comsol(simulationData* data);// = nullptr);
   // virtual void calc();
     virtual double getDe();
private:
    double m_energy =  -15.80;// eV threshold energy

};

class reactionEArs_2EArp_comsol:public reaction //e+Ars=>2e+Ar+  Ionization
{
public:
    reactionEArs_2EArp_comsol(simulationData* data);// = nullptr);
    //virtual void calc();
     virtual double getDe();
private:
    double m_energy =  -4.427;// eV threshold energy
  //  crossSection * m_cs;
};


class reactionEAr_EAr_comsol:public reaction //e+Ar=>e+Ar  ELASTIC
{
public:
    reactionEAr_EAr_comsol(simulationData* data);// = nullptr);
   // virtual void calc();
     virtual double getDe();
private:
    double m_massRatio = 0.136E-04; //m_e/m_Ar
   // crossSection * m_cs;
};


class reactionEAr_EArs_comsol:public reaction //e+Ar=>e+Ars  EXCITATION
{
public:
    reactionEAr_EArs_comsol(simulationData* data);// = nullptr);
   // virtual void calc();
     virtual double getDe();
private:
    double m_energy =  -11.50 ; //excitation energy eV
    double m_wRat = 12;  // statistical weight ratio of initial state to the excited state
    int m_Detailed = 1; //use detailed balance (if 0-- otherwise)
   // crossSection * m_cs;
};


class reactionEArs_EAr_comsol:public reaction //e+Ars=>e+Ar  EXCITATION
{
public:
    reactionEArs_EAr_comsol(simulationData* data);// = nullptr);
   // virtual void calc();
     virtual double getDe();
private:
    double m_energy =  11.50 ; //excitation energy eV
    double m_wRat = 12;  // statistical weight ratio of initial state to the excited state
    int m_Detailed = 1; //use detailed balance (if 0-- otherwise)
   // crossSection * m_cs;
};

class reactionArsArs_EArArp_comsol:public reaction //Ars+Ars=>e+Ar+Ar+
{
public:
    reactionArsArs_EArArp_comsol(simulationData* data);
        virtual void calc();
    virtual double getRate(double eps);
    virtual double getDeriv(double eps);
private:
    double m_k=3.3734e8;
};

class reactionArsAr_ArAr_comsol:public reaction //Ars+Ar=>Ar+Ar
{
public:
    reactionArsAr_ArAr_comsol(simulationData* data);
        virtual void calc();
    virtual double getRate(double eps);
    virtual double getDeriv(double eps);
private:
    double m_k=1807;
};

#endif // REACTION_H
