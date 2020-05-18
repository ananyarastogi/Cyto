#ifndef i_cell
#define i_cell
#include "config.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "main.h"

#include <vector>

using namespace std;
class sim_parameters;


class cell
{
public:
    cell();
    //cell(const cell & ) {cerr << "you copy a cell without realizing" << endl;}
    float Posn[3];
    int Numofcontacts;
    int PriorCellContacts;
    float NotInteractingSince;
    float StoppedInteracting;
};


class CTLCellType: public cell
{
public:
    CTLCellType();
    //CTLCellType(const CTLCellType & ) : cell() {cerr << "you copy a CTLCellType without realizing" << endl;}

    float PrevPosn[3];
    float Time_SameDir;
    float velocity_mag;
    float velocity[3];
    float theta;
    float phi;
    float RemainingLatentTime;
    int Trump_HistoryPosn;
    int LocationStatus;
    int SS_HistoryPosn;
    int ST_HistoryPosn;
    int Complete_HistoryPosn;
    float KillingCapacity;

    vector <float> TimeBetContacts;
};

class SomaticCellType: public cell
{
public:
    SomaticCellType();
    //SomaticCellType(const SomaticCellType & ) : cell() {cerr << "you copy a SomaticCellType without realizing" << endl;}
    float DurationFlux;
    bool ApoptoticState;
    float DyingSince;
    bool Dead;
    float Damage;
    int ContactsAfterAnalyzeTime;
    int ZombieContacts;

    int Contacts0_240;
    int Contacts0_360;
    int Contacts0_480;

    bool CheckedDeath;
    int interactionDuringDeath;
    float GrandTotalInteractTime;
    float GrandTotalInteractTimeAfterAT;
    float Data4I;
    float Data4I_AfterAT;
    float DiedAtTime;
    float TimeTakenToDie;
    vector<float> SingleContactDuration;
    vector<float> SingleContactDuration_AfterAT;



    vector<float> DamageOverTime;
};

class CTLHistory
{
public:
    CTLHistory();
    //CTLHistory(const CTLHistory & toCopy) {cerr << "you copy a history without realizing" << endl;}
    int CTLIndex;
    bool Ongoing;
    int InteractedWith;
    /*
     * Killed saves cells that died during the ongoing interaction
     * */
    int Killed;
    int DyingCellsInteract;


    float TotalContactDuration;
    float OngoingContactTime;
    bool Interacting;


    float HistoryStart;
    float HistoryEnd;
};


enum CTLMovementStatus {MovedOut, MovedIn, StayedOut, StayedIn};

void AddCells(int num_cell, int X, vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL,
              sim_parameters &ParValue, vector <CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortSpaceHistory);

void CTLmove(vector <CTLCellType> &CTL, vector<SomaticCellType> &somatic, sim_parameters & ParValue, float t);

void StoreSomaticHistory(vector <SomaticCellType> &somatic, sim_parameters & ParValue, float t);

void transform(float trans_vector[], float theta, float phi);

void AssignInitialVelocityAndAngle(CTLCellType &CTL, sim_parameters &ParValue, bool randomTimeSameDir);

float DistPoint_Vector(float OriginalPosn[], float NewPosn[], float DistFrom[], float denominator);

int CellMovementStatus(vector<CTLCellType> &CTL, sim_parameters & ParValue, int Update);

bool Inside(float CheckingPosn[], sim_parameters & ParValue);

#endif

