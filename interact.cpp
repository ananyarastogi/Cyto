#include <iostream>
#include <stdlib.h>
#include <math.h>

//#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include "main.h"
#include "random.h"
#include "cell.h"
#include "interact.h"
#include <vector>
#include <algorithm>
#include "hypothesis.h"

using namespace std;

/*
 * This class stores all details about interacting pairs
 * */
Interact::Interact()
{
    Somatic_Index=-1; ///Stores index of infected cell that is involved in interaction
    CTL_Index=-1; ///Stores index of CTL that is involved in interaction
    TimeofInteraction=0.0; ///Stores duration of ongoing interaction
    CheckedForDeath=0; ///During an interaction, in hypothesis 1 and 6, the infected cell should be checked once for death. This keeps count
    TotalInteractTime=0.0; ///At the start of each interaction, the duration of interaction is decided and stored in this field.
}

/*
 * Calculate distance between two points
 * */
float CalcDist(float a[], float b[])
{
    float distance=0, diff=0;

    for (int i=0;i<3;i++)
    {
        diff=a[i]-b[i];
        distance+= pow(diff, 2);
    }
    return pow(distance, 0.5);
}

/*
 * Checks if an infected cell and a CTL can initiate interaction on the basis of the distance between them
 * */
void CheckInteraction (vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL, sim_parameters & ParValue,
                       vector <Interact> &InteractList, float t, vector<CTLHistory> &TrumpHistory, vector<CTLHistory> &CompleteHistory,
                       vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory, DataOutput &Observation)
{
    int num_somatic_temp= somatic.size();
    int num_CTL_temp=CTL.size();



    for (int i=0;i<num_somatic_temp;i++)
    {
        for (int j=0;j<num_CTL_temp;j++)
        {

            bool InteractionStarts=false;

            /*
             * CASE 1: No cell can have more than one contact
             * */
            float ThresholdDist=(ParValue.SomaticRadius+ParValue.CTLRadius)*ParValue.ThresholdInt;
            if(ParValue.MultipleCTLSingleInfected==0)
            {

                if(CalcDist(somatic[i].Posn, CTL[j].Posn)<ThresholdDist && somatic[i].Numofcontacts==0 &&
                    CTL[j].Numofcontacts==0 && CTL[j].RemainingLatentTime>=2.0)
                {
                    if(somatic[i].ApoptoticState==0)
                    {
                        InteractionStarts=true;

                    }
                    else if(somatic[i].ApoptoticState==1 && ParValue.InteractDuringApoptosis==1)
                    {
                        InteractionStarts=true;
                        somatic[i].ZombieContacts++;
                    }
                }
            }

            /*
             * CASE 2: The infected cells can have multiple interactions.
             * */
            else if(ParValue.MultipleCTLSingleInfected==1)
            {

                if (CalcDist(somatic[i].Posn, CTL[j].Posn)<ThresholdDist &&
                    CTL[j].Numofcontacts==0 && CTL[j].RemainingLatentTime>=2.0)
                {
                    if(somatic[i].ApoptoticState==0)
                    {
                        InteractionStarts=true;
                    }
                    else if(somatic[i].ApoptoticState==1 && ParValue.InteractDuringApoptosis==1)
                    {
                        InteractionStarts=true;
                        somatic[i].ZombieContacts++;
                    }
                }
            }
            if(InteractionStarts==true)
            {
                Interact I;
                InteractionInitiated(somatic, CTL, ParValue, I, TrumpHistory, ShortTimeHistory, ShortSpaceHistory,
                                     CompleteHistory, i, j, t, Observation);
                InteractList.push_back(I);
            }

        }
    }
}

/*
 * This function checks all interactions to see if any of them need to be broken
 * */
void AlreadyInteracting(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue, vector<Interact> &InteractList,
                        vector<int> &DurationInt, float t, vector<CTLHistory> &TrumpHistory, vector<CTLHistory> &CompleteHistory,
                        vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt= InteractList.size();
    vector <int> DeleteList;

    for (int i=0;i<NumInt;i++)
    {
        InteractList[i].TimeofInteraction+=ParValue.dT;
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        int SS_Posn= CTL[c].SS_HistoryPosn;
        int ST_Posn= CTL[c].ST_HistoryPosn;
        int Trump_Posn= CTL[c].Trump_HistoryPosn;

        CompleteHistory[c].OngoingContactTime+=ParValue.dT;

        if(t>ParValue.AnalyzeAfter)
        {
            ShortTimeHistory[ST_Posn].OngoingContactTime+=ParValue.dT;
            if(ShortTimeHistory[ST_Posn].Ongoing==0)
                cout<<endl<<"Bullshit to the power bullshit."<<endl;
        }

        if(SS_Posn>=0)
        {
            ShortSpaceHistory[SS_Posn].OngoingContactTime+=ParValue.dT;
            if(ShortSpaceHistory[SS_Posn].Ongoing==0)
                cout<<endl<<"Bullshit to the power bullshit."<<endl;
        }

        if(Trump_Posn>=0 && t>ParValue.AnalyzeAfter)
        {

            TrumpHistory[Trump_Posn].OngoingContactTime+=ParValue.dT;
            if(TrumpHistory[Trump_Posn].Ongoing==0)
                cout<<endl<<"Bullshit to the power bullshit."<<endl;
        }


        somatic[s].GrandTotalInteractTime+=ParValue.dT;
        if(t>ParValue.AnalyzeAfter)
            somatic[s].GrandTotalInteractTimeAfterAT+=ParValue.dT;

        ///Breaking Interaction
        if (InteractList[i].TimeofInteraction>=InteractList[i].TotalInteractTime)
        {


            int SS_Posn= CTL[c].SS_HistoryPosn;
            int ST_Posn= CTL[c].ST_HistoryPosn;
            int Trump_Posn= CTL[c].Trump_HistoryPosn;

            somatic[s].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);

            CompleteHistory[c].OngoingContactTime=0.0;

            if(t>ParValue.AnalyzeAfter)
            {
                //random check
                if(c!= ST_Posn)
                    cout<<"Something is going wrong while storing histories!!!!!!!!!!!!!"<<endl;


                somatic[s].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_Posn].OngoingContactTime);
                ShortTimeHistory[ST_Posn].OngoingContactTime=0.0;
            }

            if(SS_Posn>=0)
            {

                ShortSpaceHistory[SS_Posn].OngoingContactTime=0.0;
            }

            if(Trump_Posn>=0 && t>ParValue.AnalyzeAfter)
            {
                TrumpHistory[Trump_Posn].OngoingContactTime=0.0;
            }


            AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);

            somatic[s].Numofcontacts--;

            CTL[c].Numofcontacts--;

            CTL[c].RemainingLatentTime=0.0;

            int RangePosn=int (InteractList[i].TimeofInteraction);
            DurationInt[RangePosn]++;

            DeleteList.push_back(i);

            CTL[c].StoppedInteracting=t;
            if(somatic[s].Numofcontacts==0)
                somatic[s].StoppedInteracting=t;
        }
        if(somatic[s].ApoptoticState==1 && InteractList[i].TimeofInteraction>=InteractList[i].TotalInteractTime && CTL[c].Numofcontacts!=0)
            cerr<<"This is crap"<<endl;



    }

    int i_max=(DeleteList.size()-1);
    for (int i=i_max;i>=0;i--)
    {
        int d=DeleteList[i];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * Once interaction has initiated, the values of contacts etc for all cells are updated
 * */
void InteractionInitiated(vector <SomaticCellType> &somatic, vector <CTLCellType> &CTL, sim_parameters & ParValue,
                          Interact &I, vector<CTLHistory> &TrumpHistory, vector<CTLHistory> &ShortTimeHistory,
                          vector<CTLHistory> &ShortSpaceHistory, vector <CTLHistory> &CompleteHistory, int i, int j, float t, DataOutput &Observation)
{
    if(CTL[j].PriorCellContacts!=0)
        CTL[j].TimeBetContacts.push_back(CTL[j].NotInteractingSince);
    CTL[j].NotInteractingSince=0.0;


    somatic[i].PriorCellContacts++;

    if (somatic[i].ApoptoticState==1)
    {
        somatic[i].interactionDuringDeath++;
        CompleteHistory[j].DyingCellsInteract++;
    }

    if(somatic[i].Numofcontacts==0 && somatic[i].StoppedInteracting!=0)
    {
        float TimeDiff=t-somatic[i].StoppedInteracting;
        Observation.TimeBetContacts_Infected.push_back(TimeDiff);

    }
    CTL[j].Numofcontacts++;
    CTL[j].PriorCellContacts++;

    somatic[i].Numofcontacts++;
    CompleteHistory[j].InteractedWith++;

    int SS_Posn= CTL[j].SS_HistoryPosn;
    int ST_Posn= CTL[j].ST_HistoryPosn;
    int Trump_Posn= CTL[j].Trump_HistoryPosn;

    if(t>=ParValue.AnalyzeAfter)
    {

        somatic[i].ContactsAfterAnalyzeTime++;
        if(CTL[j].ST_HistoryPosn!=j || ShortTimeHistory[ST_Posn].CTLIndex!=j)
            cout<<"Error in ST_HistoryPosn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl<<endl;

        if(somatic[i].ApoptoticState==0)
            ShortTimeHistory[ST_Posn].InteractedWith++;
        else if (somatic[i].ApoptoticState==1)
            ShortTimeHistory[ST_Posn].DyingCellsInteract++;
    }




    //cout<<"Alert..!! Pay attention here. If you don't change this the results will change"<<endl;
    //cout<<"This was done to change the time window observation. Just follow instructions below to set it right."<<endl;
    if(t<=240)
        somatic[i].Contacts0_240++;
    if(t<=360)
        somatic[i].Contacts0_360++;
    if(t<=480)
        somatic[i].Contacts0_480++;



    if(SS_Posn>=0)
    {
        if(somatic[i].ApoptoticState==0)
            ShortSpaceHistory[SS_Posn].InteractedWith++;
        else if (somatic[i].ApoptoticState==1)
            ShortSpaceHistory[SS_Posn].DyingCellsInteract++;
    }

    if(Trump_Posn>=0 && t>=ParValue.AnalyzeAfter)
    {
        if(somatic[i].ApoptoticState==0)
            TrumpHistory[Trump_Posn].InteractedWith++;
        else if (somatic[i].ApoptoticState==1)
            TrumpHistory[Trump_Posn].DyingCellsInteract++;
    }

    I.Somatic_Index=i;
    I.CTL_Index=j;

    do{
        I.TotalInteractTime= ParValue.DistributionInteractionDuration.getRandValue();
    }while(I.TotalInteractTime<0);



}

/*
 * Based on hypothesis value, the function calls the respective hypothesis functions
 * */
void CheckforDeath(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                   vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                   vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory, vector <CTLHistory> &CompleteHistory)
{

    int HypothesisTemp=ParValue.Hypothesis;
    switch(HypothesisTemp)
    {
    case 1:
        hypothesis1_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;

    case 2:
        hypothesis2_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;

    case 3:
        hypothesis3_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;
    case 4:
        hypothesis4_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;

    case 5:
        hypothesis5_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;
    case 6:
        hypothesis6_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;
    case 7:
        hypothesis7_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;
    case 8:
        hypothesis8_death(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);
        break;
    default:
        cout<<"Bullshit hypothesis!!"<<endl;
        exit(0);
    }
}
