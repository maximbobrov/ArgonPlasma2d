#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <QGridLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QScrollBar>
#include <QLabel>
#include <QVector>
#include <QCheckBox>
#include <QProgressBar>
#include "glwidget.h"
#include <QComboBox>
#include <QFrame>

#include "simulationdata.h"
#include "simulationsolver.h"
#include "simulationtools.h"
#define NX 257
#define NY 257

class reactionSolver;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

    struct plotStruct{
        double arr[NX][NY];
        double **arrRef;
        QString name;
        double max;
        double min;
        double scale;
        bool visible;
    };

    struct storeStruct{
        QVector<plotStruct> plots;
        double time;
    };

public slots:
    bool solveNewton();
    void initData();
    void updateData(int);
    void simulateData(bool);
    void stopAnim(bool);
    void replotGraph(int);
    void setFields();
    void setVisualArr(const QString &text);


private:
    QWidget* m_widget;
    GLWidget *glWidget;
    QCustomPlot* m_customPlot;

    QGridLayout* m_grid;

    QGridLayout* m_visualGrid;
    QGridLayout* m_settingGrid;

    QGridLayout* m_simulSettingGrid;
    QGridLayout* m_visualSettingGrid;

    QPushButton* m_simulateButton;
    QPushButton* m_stopButton;
    QPushButton* m_initButton;

    QLineEdit* m_textStartTime;
    QLineEdit* m_textEndTime;
    QLineEdit* m_textDeltaTime;

    QLabel* m_timeLabel;
    QProgressBar* m_progressBar;


    QComboBox* m_fieldsComboBox;
    QScrollBar* m_timeScrollBar;
    QScrollBar* m_scallingBar;
    QScrollBar* m_crossBar;



    simulationData* m_data;
    solverNe* m_sNe;
    solverEnergy* m_sEn;
    solverPhi* m_sPhi;
    QVector<solverHeavySpicies*> m_sHeavy;
    int m_numberHeavySpicies;
    simulationData::simulationField2d* m_fNe;
    simulationData::simulationField2d* m_fEnergy;
    simulationData::simulationField2d* m_fPhi;
    QVector<simulationData::simulationField2d*> m_fHeavy;
    simulationData::simulationParameters* m_pParam;
    reactionSolver * m_rSolver;

    QVector<plotStruct> m_plots;
    QVector<storeStruct> m_storage;
    double m_maxY, m_minY;
    bool m_animStopped;
    double m_startTime,m_endTime;
    double m_time;
    double** m_visualArr;
    QString m_visualArrName;
    bool m_initiolized;

private:
    void saveInStorage();
    void addPlot(double** arr,char* name,double scale = 1.0);
    void addPlotXY(double *arr,double*xx, char *name, int size, double scale = 1.0);
};

#endif // MAINWINDOW_H
