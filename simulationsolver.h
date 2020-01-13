#ifndef SIMULATIONSOLVER_H
#define SIMULATIONSOLVER_H
#include "simulationdata.h"


class simulationSolver
{
public:

    simulationSolver(simulationData* pData = nullptr);


    virtual double solve(int numberIteration);
    virtual double getRhs();
    virtual double getRhsAt(int i, int j);
    virtual void getStepEuler();
    virtual void setBc();
    virtual void updateMatrix();
    virtual double getNewtonRhs(int i, int j);
    virtual double init(double value);

    virtual ~simulationSolver();

protected:
    simulationData* m_pData;
    simulationData::simulationField2d* m_field;
    double** m_a;
    double** m_bm;
    double** m_bp;
    double** m_cm;
    double** m_cp;

    double** m_aRHS;
    double** m_aNu;
    double** m_aMu;
    double** m_mask;
    double** m_maskValue;
        int m_charge;
};

class solverNe : public simulationSolver
{
public:
    solverNe(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double getRhsAt(int i,int j);
    virtual void setBc();
    virtual ~solverNe();
};


class solverEnergy : public simulationSolver
{
public:
    solverEnergy(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double getRhsAt(int i,int j);
    virtual void setBc();
    virtual ~solverEnergy();
};

class solverPhi : public simulationSolver
{
public:
    solverPhi(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double solve(int numberIteration);
    virtual void setBc();

    virtual ~solverPhi();
private:
    double q_e = 1.6*10e-19;
    double eps_0 = 8.85*10e-12;
};

class solverHeavySpicies : public simulationSolver
{
public:
    solverHeavySpicies(simulationData* pData = nullptr, int num = 0);
    virtual double getRhs();
    virtual double getRhsAt(int i,int j);
    virtual void setBc();

    virtual ~solverHeavySpicies();

    simulationData::SpecieName m_specie;
private:
};

#endif // SIMULATIONSOLVER_H
