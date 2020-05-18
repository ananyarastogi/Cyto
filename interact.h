#ifndef i_interact
#define i_interact
#include "config.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "cell.h"
#include "output.h"
#include <vector>

using namespace std;
class sim_parameters;
class cell;
class CTLCellType;
class SomaticCellType;
class CTLHistory;



class Interact
{
public:
    Interact();
    int Somatic_Index;
    int CTL_Index;
    float TimeofInteraction; //Duration of ongoing interaction
    /*
    float DurationFlux;
    bool FluxStatus; //There is 60% ish chance that flux will start. So when interaction starts there is a 60% chance of flux
    bool FluxYesorNo;
    double FluxStartAt;
    */
    float TotalInteractTime;
    bool CheckedForDeath;



};


void CheckInteraction (vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, float t,  vector <CTLHistory> &TrumpHistory, vector<CTLHistory> &CompleteHistory,
                       vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory, DataOutput &Observation);

void AlreadyInteracting (vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL, sim_parameters & ParValue, vector<Interact> &InteractList,
                         vector<int> &DurationInt, float t, vector<CTLHistory> &TrumpHistory, vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory,
                         vector<CTLHistory> &ShortSpaceHistory);

void CheckforDeath (vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL, sim_parameters & ParValue,
                    vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt,  vector <CTLHistory> &TrumpHistory,
                    vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory, vector <CTLHistory> &CompleteHistory);

float CalcDist(float a[3], float b[3]);

void InteractionInitiated(vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL, sim_parameters & ParValue,
                          Interact &I, vector<CTLHistory> &TrumpHistory, vector<CTLHistory> &ShortTimeHistory,
                          vector<CTLHistory> &ShortSpaceHistory, vector <CTLHistory> &CompleteHistory, int i, int j, float t, DataOutput &Observation);

#endif
