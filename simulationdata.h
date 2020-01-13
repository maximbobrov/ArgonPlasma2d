#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H
#include <iostream>
#include <vector>
#include "reaction.h"

class simulationData
{
public:
    simulationData(int iCellsX,int iCellsY);
    ~simulationData();

    enum SpecieName
    {

        Ar_plus,
        Ar_star,

        e,
        En,
        phi,
        Ar
    };

    enum ReactionName
    {
        comsol_eAr_eAr,
        comsol_eAr_eArs,
        comsol_eArs_eAr,
        comsol_eAr_2eArp,
        comsol_eArs_2eArp,
        comsol_ArsArs_EArArp,
        comsol_ArsAr_ArAr
        //comsol_eAr_2eArp
    };

    /*enum SpecieName
    {
       Ar,
       Ar_star,
       Ar_plus
    };

    enum ReactionName
    {
       comsol_eAr_eAr,
       comsol_eAr_eArs,
       comsol_eAr_2eArp,
       comsol_eArs_2eArp
    };
*/

    struct simulationField
    {
        char* name;
        double* arr;
        double* arrPrev;
        int cellsNumber;
        simulationField(int cellsNumber, char* name);
        private: void init(int cellsNumber, char* name);
    };

    struct boundary
    {
        int num;
        int* arrI;
        int* arrJ;
        int* arrI_inter;
        int* arrJ_inter;

        boundary();
        void init(int cellsX, int cellsY, int** mask);

    };


    struct simulationField2d
    {
        char* name;
        double** arr;
        double** arrPrev;
        int cellsX;
        int cellsY;
        simulationData::SpecieName m_specie;
        simulationField2d(int cellsX, int cellsY, char* name,simulationData::SpecieName specie);
        private: void init(int cellsX, int cellsY, char* name,simulationData::SpecieName specie);
    };

    struct simulationParameters
    {
        double** arrDe;
        double** arrDeps;
        double** arrDomega;
        double** arrMue;
        double** arrMueps;
        double** arrMuomega;
        double** arrTe;
        double** arrEx;
        double** arrEy;
        double** arrLaplPhi;
        double** arrMaskPhi;
        double** arrMaskPhiValue;
        double** arrMaskNe;
        double** arrMaskNeValue;
        double** arrMaskEnergy;
        double** arrMaskEnergyValue;
        double** arrMaskHeavy;
        double** arrMaskHeavyValue;
        int ** arrBoundMask;
        double** arrEps;

        boundary bound;
        int cellsX;
        int cellsY;
        double rho; //mixture density
        double p; //pressure
        double T; //temperature
        double mAr; //argon molar mass
        double N; //neutral number denisty
        simulationParameters(int cellsX, int cellsY);
        private: void init(int cellsX, int cellsY);
    };

    const double q=1.6022e-19; //coulumbs elementary charge
    const double k_B_const=1.3806e-23;// boltzmann constant

    void setDt(double dt = 0.000003);
    void setDx(double dx = 0.01);
    void setDy(double dy = 0.01);

    void setCellsXNumber(int cellsXNumber);
    void setCellsYNumber(int cellsYNumber);

    int getCellsXNumber();
    int getCellsYNumber();
    double getDt();
    double getDx();
    double getDy();
    void updateParams();
    double** getArrTe();
    simulationField2d* getFieldNe();
    simulationField2d* getFieldEnergy();
    simulationField2d* getFieldPhi();
    simulationField2d* getFieldHeavySpicies(int num);
    int getHeavySpiciesCharge(int num);
    int getNumberHeavySpicies();
    double getN(); //total number of particles in one m^3
    double** getReactionRate(simulationData::ReactionName reactName);
    double getReactionDe(simulationData::ReactionName reactName);

    void calcReaction(ReactionName reactName);
    simulationParameters *getParameters();
    std::vector<reaction *> m_reactions;
private:
    friend class reactionSolver;
    double m_dt;
    double m_dx;
    double m_dy;
    int m_cellsX;
    int m_cellsY;
    int m_numberHeavySpicies;
    simulationField2d* m_fieldNe;
    simulationField2d* m_fieldEnergy;
    simulationField2d* m_fieldPhi;
    std::vector<simulationField2d* > m_fieldsHeavySpecies;
    std::vector<int> m_chargeHeavySpecies;
    simulationParameters* m_params;

};

#endif // SIMULATIONDATA_H
