#ifndef i_hypothesis
#define i_hypothesis
#include "config.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "cell.h"

#include <vector>


//Equal probability
void hypothesis1_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector <CTLHistory> &CompleteHistory, vector <CTLHistory> &ShortTimeHistory, vector <CTLHistory> &ShortSpaceHistory);

//Damage rate based on duration of contact
void hypothesis2_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector <CTLHistory> &CompleteHistory, vector <CTLHistory> &ShortTimeHistory, vector <CTLHistory> &ShortSpaceHistory);

//Damage and repair- both time based,
//Damage is same as hypothesis 2, repair is proportional to damage.
void hypothesis3_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters  & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);



void hypothesis4_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);

void hypothesis5_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);

void hypothesis6_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);

void hypothesis7_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);

void hypothesis8_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);

void hypothesis9_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory);
#endif
