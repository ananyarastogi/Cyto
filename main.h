#ifndef i_main
#define i_main
#include <iostream>
//#include <stdlib.h>
#include <math.h>
#include "cell.h"
#include "interact.h"
#include "outsidefunc.h"
#include "output.h"
#include "distribution.h"
#include "parameters.h"
#include <sstream>
#include <vector>

#define HeatMapHyp 1

#define MaxContacts 10000


class SumOfSquare
{
public:
    SumOfSquare();
    float FracDeadCells;
    float Fig4I;
    float DeadFig4E;
    float AliveFig4E;
    float Intact4F;
    float Killed4F;
    float Dead4D;
    float Intact4D;
    float AbsoluteSumOfSq;


    float FracDeadCells_AfterAT;
    float Fig4I_AfterAT;
    float DeadFig4E_AfterAT;
    float AliveFig4E_AfterAT;
    float Intact4F_AfterAT;
    float Killed4F_AfterAT;
    float Dead4D_AfterAT;
    float Intact4D_AfterAT;
    float AbsoluteSumOfSq_AfterAT;



    float MMH_FracDeadCells;
    float MMH_Fig4I;
    float MMH_DeadFig4E;
    float MMH_AliveFig4E;
    float MMH_Intact4F;
    float MMH_Killed4F;
    float MMH_Dead4D;
    float MMH_Intact4D;
    float MMH_AbsoluteSumOfSq;


    float MMH_FracDeadCells_AfterAT;
    float MMH_Fig4I_AfterAT;
    float MMH_DeadFig4E_AfterAT;
    float MMH_AliveFig4E_AfterAT;
    float MMH_Intact4F_AfterAT;
    float MMH_Killed4F_AfterAT;
    float MMH_Dead4D_AfterAT;
    float MMH_Intact4D_AfterAT;
    float MMH_AbsoluteSumOfSq_AfterAT;

    int n_AbsSS;
    int MMH_n_AbsSS;
    int n_AbsSSAfterAT;
    int MMH_n_AbsSSAfterAT;





};



void analysis(sim_parameters &ParValue);


void setSimulationStatus(bool stop);


void RunSimulation(vector<float> &AvgAliveSom, vector<float> &AvgDeadSom, sim_parameters & ParValue,
                   vector<int> &DurationInt, int x, int y, float &AvgDead, string appAddress,
                   vector<float> &AvgDeadSomAfterAT, vector<float> &AvgAliveSomAfterAT, DataOutput &Observation,
                   vector <SumOfSquare> &Results_SumSq, int &MMH_n_AbsSumSq, int &MMH_n_AbsSumSq_AfterAT, int &n_AbsSumSq,
                   int &n_AbsSumSq_AfterAT, vector<vector<float>> &AllAliveData, vector<vector<float>> &AllDeadData,
                   vector<vector<float>> &AllAliveData_AfterAT, vector<vector<float>> &AllDeadData_AfterAT);



void DoSimulate();

void WriteReport(string fn, sim_parameters&  ParValue, DefaultParameterRange ParRange);

void GenerateRScript(string fn, sim_parameters & ParValue);


void ParametersChanging(sim_parameters &ParValue, float TempRatio, float vmin, float vmax);

//bool BullshitDetector(sim_parameters ParValue);


void AnalysisMethod(sim_parameters & ParValue, int x, string appAddress,
                      vector<int> DurationInt, int AfterAT, vector <SumOfSquare> &Results_SumSq, int &MMH_n_Abs, int &n_AbsSumSq, DataOutput &Observation);

void CalcOutputVector(vector<float> &AvgOutputVector, vector<float> &StdDevOutputVector, vector<vector<float>> DataVector, sim_parameters & ParValue);



#endif


