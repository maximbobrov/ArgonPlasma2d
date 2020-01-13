#ifndef SIMULATIONTOOLS_H
#define SIMULATIONTOOLS_H

class simulationTools
{
public:
    static double ddx(double **arr, double dx, int i, int j);
    static double ddy(double **arr, double dy, int i, int j);

    static double ddxCentral(double **arr, int NX, double dx,  int i, int j);
    static double ddyCentral(double **arr, int NY, double dy,  int i, int j);

    static double ddxUpWind(double **arr, int NX, double dx, int i, int j, double flux);
    static double ddyUpWind(double **arr, int NY, double dy, int i, int j, double flux);

    static double gauss(double x, double d);
    simulationTools();
};
#endif // SIMULATIONTOOLS_H
