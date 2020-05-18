#ifndef i_dataStorage
#define i_dataStorage
#include "config.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "cell.h"

#include <vector>

using namespace std;
class sim_parameters;
class cell;
class CTLCellType;
class SomaticCellType;


void StoreCTLMoveHistory(vector<CTLCellType> &CTL, sim_parameters & ParValue, float t, vector<CTLHistory> &TrumpHistory,
                         vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory, vector<CTLHistory> &CompleteHistory, vector<Interact> &InteractList,
                         vector<SomaticCellType> &somatic);

void UpdateHistory(vector<CTLCellType> &CTL, sim_parameters & ParValue, vector<CTLHistory> &TrumpHistory,
                   vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory,
                   vector <CTLHistory> &CompleteHistory, int CTLIndex, float t);

#endif

