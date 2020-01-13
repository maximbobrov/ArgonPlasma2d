#include "simulationtools.h"
#include "math.h"

simulationTools::simulationTools()
{

}



double simulationTools::ddx(double** arr, double dx, int i, int j)
{
    return i == 0 ? (arr[1][j] - arr[0][j]) / (dx) : (arr[i][j]-arr[i-1][j])/(dx);
}

double simulationTools::ddy(double** arr, double dy, int i, int j)
{
    return j == 0 ? (arr[i][1] - arr[i][0]) / (dy) : (arr[i][j]-arr[i][j-1])/(dy);
}


double simulationTools::ddxCentral(double **arr, int NX, double dx, int i, int j)
{
    return i == 0 ? (-3.0*arr[i][j] + 4.0*arr[i+1][j] - arr[i+2][j])/(2.0*dx) : (i == NX-1 ? (3.0*arr[i][j] - 4.0*arr[i-1][j] + arr[i-2][j])/(2.0*dx) : (arr[i+1][j]-arr[i-1][j])/(2.0*dx));
}

double simulationTools::ddyCentral(double **arr, int NY, double dy,  int i, int j)
{
    return j == 0 ? (-3.0*arr[i][j] + 4.0*arr[i][j+1] - arr[i][j+2])/(2.0*dy) : (j == NY-1 ? (3.0*arr[i][j] - 4.0*arr[i][j-1] + arr[i][j-2])/(2.0*dy) : (arr[i][j+1]-arr[i][j-1])/(2.0*dy));
}

double simulationTools::ddxUpWind(double **arr, int NX, double dx, int i, int j, double flux)
{
    int ip = (flux>0);
    int im = (flux<=0);
    if(i-im < 0)
    {
        im = 0;
        ip = 1;
    }

    if(i+ip > NX-1)
    {
        im = 1;
        ip = 0;
    }
        return (arr[i+ip][j]-arr[i-im][j])/(dx);
}

double simulationTools::ddyUpWind(double **arr, int NY, double dy,  int i, int j, double flux)
{
    int jp = (flux>0);
    int jm = (flux<=0);
    if(j-jm < 0)
    {
        jm = 0;
        jp = 1;
    }

    if(j+jp > NY-1)
    {
        jm = 1;
        jp = 0;
    }
    return (arr[i][j+jp]-arr[i][j-jm])/(dy);
}


double simulationTools::gauss(double x, double d)
{
    return exp(-x*x/(d*d));
}
