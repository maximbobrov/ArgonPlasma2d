
#include "mainwindow.h"
#include "qcustomplot.h"
#include <QGridLayout>
#include <QPushButton>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <stdio.h>
#include <reactionsolver.h>

#define NX 257
#define NY 257



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    m_widget=new QWidget();

    m_widget->setObjectName("central");
    m_widget->setMinimumSize(1100,700);

    m_grid = new QGridLayout(this);
    m_visualGrid = new QGridLayout(this);
    m_settingGrid = new QGridLayout(this);

    m_simulSettingGrid = new QGridLayout(this);
    m_visualSettingGrid = new QGridLayout(this);

    m_grid->addLayout(m_settingGrid,0,0,3,1);
    m_grid->addLayout(m_visualGrid,0,1,3,3);

    QFrame *lineSim = new QFrame(this);
    lineSim->setFrameShape(QFrame::Box);
    lineSim->setFrameShadow(QFrame::Sunken);

    QFrame *lineVis = new QFrame(this);
    lineVis->setFrameShape(QFrame::Box);
    lineVis->setFrameShadow(QFrame::Sunken);

    QFont font;
    font.setWeight(QFont::Bold);
    font.setPixelSize(12);

    QLabel * simLabel = new QLabel("Simulation settings");
    simLabel->setAlignment(Qt::AlignCenter);
    simLabel->setFont(font);

    QLabel * visLabel = new QLabel("Visualize settings");
    visLabel->setAlignment(Qt::AlignCenter);
    visLabel->setFont(font);

    m_settingGrid->addWidget(lineSim                                    , 0, 0, 6, 20);
    m_settingGrid->addWidget(simLabel                                   , 0, 1, 1, 18);
    m_settingGrid->addLayout(m_simulSettingGrid                         , 1, 1, 5, 18);
    m_settingGrid->addWidget(lineVis                                    , 7, 0, 6, 20);
    m_settingGrid->addWidget(visLabel                                   , 7, 1, 1, 18);
    m_settingGrid->addLayout(m_visualSettingGrid                        , 8, 1, 5, 18);


    m_timeScrollBar = new QScrollBar(Qt::Horizontal,this);
    m_scallingBar = new QScrollBar(Qt::Horizontal,this);
    m_crossBar = new QScrollBar(Qt::Horizontal,this);
    m_fieldsComboBox = new QComboBox(this);


    m_visualSettingGrid->addWidget(new QLabel("Visualized field:")        , 0, 0, 1, 3);
    m_visualSettingGrid->addWidget(m_fieldsComboBox                       , 1, 0, 1, 3);
    m_visualSettingGrid->addWidget(new QLabel("Time:")                    , 2, 0, 1, 3);
    m_visualSettingGrid->addWidget(m_timeScrollBar                        , 3, 0, 1, 3);
    m_visualSettingGrid->addWidget(new QLabel("Scalling:")                , 4, 0, 1, 3);
    m_visualSettingGrid->addWidget(m_scallingBar                          , 5, 0, 1, 3);
    m_visualSettingGrid->addWidget(new QLabel("Cross section coordinate:"), 6, 0, 1, 3);
    m_visualSettingGrid->addWidget(m_crossBar                             , 7, 0, 1, 3);
    m_visualSettingGrid->setAlignment(Qt::AlignTop);


    m_textStartTime = new QLineEdit(this);
    m_textStartTime->setMaximumWidth(65);
    m_textDeltaTime = new QLineEdit(this);
    m_textDeltaTime->setMaximumWidth(65);
    m_textEndTime =   new QLineEdit(this);
    m_textEndTime->setMaximumWidth(65);

    m_simulateButton = new QPushButton("Simulate");
    m_stopButton = new QPushButton("Stop");
    m_initButton = new QPushButton("Init");
    m_progressBar = new QProgressBar(this);
    m_timeLabel = new QLabel("Progress:");
    m_progressBar->setRange(0,100);


    m_simulSettingGrid->addWidget(new QLabel("Start time:")            , 0, 0, 1, 2);
    m_simulSettingGrid->addWidget(m_textStartTime                      , 1, 0, 1, 2);
    m_simulSettingGrid->addWidget(new QLabel("dt:")                    , 0, 2, 1, 2);
    m_simulSettingGrid->addWidget(m_textDeltaTime                      , 1, 2, 1, 2);
    m_simulSettingGrid->addWidget(new QLabel("End time:")              , 0, 4, 1, 2);
    m_simulSettingGrid->addWidget(m_textEndTime                        , 1, 4, 1, 2);
    m_simulSettingGrid->addWidget(m_simulateButton                     , 2, 0, 1, 2);
    m_simulSettingGrid->addWidget(m_stopButton                         , 2, 2, 1, 2);
    m_simulSettingGrid->addWidget(m_initButton                         , 2, 4, 1, 2);
    m_simulSettingGrid->addWidget(new QLabel(" ")                      , 3, 0, 1, 6);
    m_simulSettingGrid->addWidget(m_timeLabel                          , 4, 0, 1, 6);
    m_simulSettingGrid->addWidget(m_progressBar                        , 5, 0, 1, 6);
    m_simulSettingGrid->setAlignment(Qt::AlignTop);


    m_customPlot = new QCustomPlot(this);
    glWidget = new GLWidget;
    m_visualGrid->addWidget(glWidget,0,0,2,2);
    m_visualGrid->addWidget(m_customPlot,2,0,2,2);

    m_data = new simulationData(NX,NY);
    m_visualArr = new double*[NX];
    for (int i = 0; i < NX; ++i)
        m_visualArr[i] = new double [NY];

    m_data->setDx(1.0/NX);
    m_data->setDy(0.3/NY);
    m_data->setCellsXNumber(NX);
    m_data->setCellsYNumber(NY);
    m_data->setDt(1e-6);

    m_sNe = new solverNe(m_data);
    m_sEn = new solverEnergy(m_data);
    m_sPhi = new solverPhi(m_data);

    m_rSolver= new reactionSolver(m_data);

    m_textStartTime->setText(QString().number(0.0));
    m_textDeltaTime->setText(QString().number(m_data->getDt()));
    m_textEndTime->setText(QString().number(1e-4));

    m_scallingBar->setRange(1,100);
    m_scallingBar->setValue(1);

    m_crossBar->setRange(0, NY-1);
    m_crossBar->setValue((NY-1)/2);


    connect (m_simulateButton, SIGNAL(clicked(bool)), this, SLOT(simulateData(bool)));
    connect (m_stopButton, SIGNAL(clicked(bool)), this, SLOT(stopAnim(bool)));
    connect (m_initButton, SIGNAL(clicked(bool)), this, SLOT(initData()));
    connect(m_timeScrollBar, SIGNAL(valueChanged(int)), this, SLOT(replotGraph(int)));
    connect(m_scallingBar, SIGNAL(valueChanged(int)), this, SLOT(setFields()));
    connect(m_crossBar, SIGNAL(valueChanged(int)), this, SLOT(setFields()));
    connect(m_fieldsComboBox, SIGNAL(currentIndexChanged(const QString)), this, SLOT(setVisualArr(const QString)));

    m_widget->setLayout(m_grid);
    setCentralWidget(m_widget);
    setWindowTitle("PlasmaSolver");
    m_animStopped=true;
    initData();
}

MainWindow::~MainWindow()
{

}

void MainWindow::setFields()
{
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}

void MainWindow::replotGraph(int number)
{
    m_timeLabel->setText("Progress(time = " + QString().number(m_storage[number].time)+"):");

    m_maxY=-1e30;
    m_minY=+1e30;
    m_plots.clear();
    m_plots = m_storage[number].plots;
    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel("x");
    m_customPlot->yAxis->setLabel("y");
    QColor colors[6] = {QColor(255,0,0), QColor(255,155,0), QColor(0,255,0),QColor(0,0,255), QColor(255,0,255),QColor(0,255,255)};
    QVector<double> x, y;
    x.resize(NX);
    y.resize(NX);

    for (int j = 0; j < m_plots.size(); ++j)
    {
        m_plots[j].max=-1e30;
        m_plots[j].min=+1e30;
        for (int i = 0; i < NX; ++i)
        {
            m_plots[j].max = m_plots[j].arr[i][(int)(m_crossBar->value())] > m_plots[j].max ? m_plots[j].arr[i][(int)(m_crossBar->value())] : m_plots[j].max;
            m_plots[j].min = m_plots[j].arr[i][(int)(m_crossBar->value())] < m_plots[j].min ? m_plots[j].arr[i][(int)(m_crossBar->value())] : m_plots[j].min;
        }
        m_maxY = m_plots[j].max > m_maxY ? m_plots[j].max : m_maxY;
        m_minY = m_plots[j].min < m_minY ? m_plots[j].min : m_minY;
    }
    double globalMax = m_maxY;
    double globalMin = m_minY;
    for (int j = 0; j < m_plots.size(); ++j)
    {
        m_plots[j].scale = (m_maxY - m_minY) / (m_plots[j].max - m_plots[j].min);
        globalMax = m_plots[j].scale * m_plots[j].max > globalMax ? m_plots[j].scale * m_plots[j].max : globalMax;
        globalMin = m_plots[j].scale * m_plots[j].min < globalMin ? m_plots[j].scale * m_plots[j].min : globalMin;
    }



    for (int j = 0; j < m_plots.size(); ++j)
    {
        if(m_visualArrName ==  m_plots[j].name)
        {
            for (int ii = 0; ii < NX; ++ii)
                for (int jj = 0; jj < NY; ++jj)
                    m_visualArr[ii][jj] = m_plots[j].arr[ii][jj];
            glWidget->setField( m_visualArr , m_pParam->arrEps , m_pParam->arrMaskPhi, NX, NY, m_data->getDx(), m_data->getDy(), m_crossBar->value(), 1.0*m_scallingBar->value());
            glWidget->repaint();
        }

        m_customPlot->addGraph();
        for (int i = 0; i < NX; ++i)
        {
            x[i] = i / (NX * 1.0 / 2)-1;
            y[i] =  m_plots[j].scale * m_plots[j].arr[i][(int)(m_crossBar->value())];
        }

        m_customPlot->graph(j)->setName(m_plots[j].name + "( scale = " + QString().number(m_plots[j].scale) + ")");
        m_customPlot->graph(j)->addToLegend();
        m_customPlot->graph(j)->setData(x, y);
        m_customPlot->graph(j)->setLineStyle((QCPGraph::LineStyle)(1));
        QPen graphPen;
        graphPen.setColor(colors[j]);
        graphPen.setWidthF(2);
        m_customPlot->graph(j)->setPen(graphPen);
        m_customPlot->xAxis->setRange(-1, 1);
        m_customPlot->yAxis->setRange(globalMin, globalMax);
        m_customPlot->legend->setVisible(true);
        m_customPlot->replot();

    }
}

void MainWindow::saveInStorage()
{
    for (int k = 0; k < m_plots.size(); ++k)
    {
        for (int i = 0; i < NX; ++i)
        {
            for (int j = 0; j < NY; ++j)
            {
                m_plots[k].arr[i][j] = m_plots[k].arrRef[i][j];
            }
        }
    }
    storeStruct storeItem;
    storeItem.plots = m_plots;
    storeItem.time = m_time;
    m_storage.push_back(storeItem);
}

void MainWindow::addPlot(double **arr, char *name, double scale)
{
    m_fieldsComboBox->addItem(QString(name));
    plotStruct plot;
    for (int i = 0; i < NX; ++i)
    {
        for (int j = 0; j < NY; ++j)
        {
            plot.arr[i][j] = arr[i][j];
        }
    }
    plot.arrRef = arr;
    plot.name = name;
    plot.scale = scale;
    m_plots.push_back(plot);
}

void MainWindow::initData()
{
    m_initiolized = false;

    //////INIT FIELDS AND VARIABLES
    {

        m_pParam =m_data->getParameters();
        m_fNe = m_data->getFieldNe();
        m_fEnergy = m_data->getFieldEnergy();
        m_fPhi = m_data->getFieldPhi();
        m_numberHeavySpicies = m_data->getNumberHeavySpicies();
        m_time = 0.0;
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_fHeavy.push_back(m_data->getFieldHeavySpicies(j));
            m_sHeavy.push_back(new solverHeavySpicies(m_data, j));
        }

        m_sPhi->init(0.0);
        for (int i = 0; i < NX; ++i)
        {
            for (int j = 0; j < NY; ++j)
            {
                double x=(i-3*NX/4)*m_data->getDx();
                double y=(j-NY/2)*m_data->getDy();
                double r=x*x+y*y;
                m_fNe->arr[i][j] = 1e5 + 1e11*simulationTools::gauss(sqrt(x*x+y*y), NY*m_data->getDy()*0.1);
                m_fNe->arrPrev[i][j] =m_fNe->arr[i][j];
                m_fEnergy->arr[i][j] = 0.01*m_fNe->arr[i][j];
                m_fEnergy->arrPrev[i][j] = m_fEnergy->arr[i][j];
                for (int h = 0; h < m_numberHeavySpicies; ++h)
                {
                    m_fHeavy[h]->arr[i][j] =(m_fNe->arr[i][j] * m_pParam->T*8.314)/(m_pParam->p*6.022e23);
                    m_fHeavy[h]->arrPrev[i][j] = m_fHeavy[h]->arr[i][j];
                }
            }
        }
        m_data->updateParams();
        m_sPhi->solve(20);
    }


    //////INIT GUI
    {
        m_plots.clear();
        m_fieldsComboBox->clear();
        m_storage.clear();
        m_progressBar->setValue(0.0);
        m_visualArrName = m_fNe->name;

        addPlot(m_fNe->arr, m_fNe->name);
        addPlot(m_pParam->arrTe, m_fEnergy->name);
        addPlot(m_fPhi->arr, m_fPhi->name);
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            addPlot(m_fHeavy[j]->arr, m_fHeavy[j]->name);
        }
        saveInStorage();
        replotGraph(0);
        m_timeScrollBar->setRange(0,0);
    }
    m_initiolized = true;
}

void MainWindow::updateData(int numberIt)
{
    double dt= m_data->getDt();


    //dt*= 1.05;
    m_data->setDt(1e-10);//dt);
    //for (int i = 0; i < numberIt; ++i)
    {

        solveNewton();
        m_sNe->getStepEuler();
        m_sEn->getStepEuler();
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_sHeavy[j]->getStepEuler();
        }
    }
}


/*void MainWindow::updateData(int numberIt)
{
    double dt= m_data->getDt();


    dt*= 1.05;
    m_data->setDt(dt);
    for (int i = 0; i < 2; ++i)
    {
        m_data->updateParams();
        m_sNe->solve(5);
        m_sEn->solve(5);
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_sHeavy[j]->solve(5);
        }
        m_sPhi->solve(50);

    }



    m_sNe->getStepEuler();
    m_sEn->getStepEuler();
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_sHeavy[j]->getStepEuler();
    }
}*/

bool MainWindow::solveNewton()
{

    double **ne=m_data->getFieldNe()->arr;
    double **ars=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;
    double **arp=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_plus)->arr;

    double **en=m_data->getFieldEnergy()->arr;

    double **ne0=m_data->getFieldNe()->arrPrev;
    double **ars0=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev;
    double **arp0=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_plus)->arrPrev;
    double **en0=m_data->getFieldEnergy()->arrPrev;

    for (int nn=0;nn<10;nn++)
    {
        m_sPhi->solve(1);
        m_data->updateParams();

        for (int i=1; i <m_data->getCellsXNumber() - 2; i++)
        {
            for (int j=1; j <m_data->getCellsXNumber() - 2; j++)
            {
                ne[i][j]=ne[i][j]*0.5+0.5*m_sNe->getNewtonRhs(i,j)*m_pParam-> arrMaskNe[i][j];
                en[i][j]=en[i][j]*0.5+0.5*m_sEn->getNewtonRhs(i,j)*m_pParam-> arrMaskNe[i][j];
                ars[i][j]=ars[i][j]*0.5+0.5*m_sHeavy[simulationData::SpecieName::Ar_star]->getNewtonRhs(i,j)*m_pParam-> arrMaskNe[i][j];//*m_pParam-> arrMaskNe[i][j];;
                arp[i][j]=arp[i][j]*0.5+0.5*m_sHeavy[simulationData::SpecieName::Ar_plus]->getNewtonRhs(i,j)*m_pParam-> arrMaskNe[i][j];//*m_pParam-> arrMaskNe[i][j];;
            }
        }
    }
    m_sPhi->solve(2);
    return true;

}

void MainWindow::simulateData(bool status)
{
    m_animStopped=false;
    int saveNum = (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()) / m_textDeltaTime->text().toDouble();
    if(saveNum< 20)
        saveNum = 1;
    else
        saveNum = saveNum * 1.0 / 20.0;
    while(!m_animStopped &&  m_time < (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()) - saveNum * m_textDeltaTime->text().toDouble())
    {
        m_time += m_textStartTime->text().toDouble() + saveNum * m_textDeltaTime->text().toDouble();
        updateData(saveNum);

        saveInStorage();
        m_timeScrollBar->setRange(0, m_storage.size() - 1);
        m_timeScrollBar->setValue(m_storage.size() - 1);
        replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
        QCoreApplication::processEvents();

        m_progressBar->setValue(100.0 * (m_time - m_textStartTime->text().toDouble()) / (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()));

    }
    m_progressBar->setValue(100);
}

void MainWindow::stopAnim(bool)
{
    m_animStopped=true;
}

void MainWindow::setVisualArr(const QString &text)
{
    if(m_initiolized)
    {
        m_visualArrName = text;
        replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
    }
}





