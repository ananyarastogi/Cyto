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
#include "dataStorage.h"
#include <vector>
#include <algorithm>

using namespace std;

/*
 * Equal probability of killing at each interaction (does not depend on previuos interactions had)- null hypothesis
 * The decision to die is independent of duration of interaction
 * At the start of each interaction- decide if cell will die
 * */
void hypothesis1_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{

    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;

    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;

        if(ParValue.DeathDecisionAtStart==1)
        {
            if(InteractList[i].CheckedForDeath==0)
            {
                if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=ParValue.ProbDie_Interact_Hyp1))
                {
                    somatic[s].ApoptoticState=1;
                }
                InteractList[i].CheckedForDeath=1;


            }
        }
        else if(ParValue.DeathDecisionAtStart==0)
        {

            if ((somatic[s].ApoptoticState==0 && (InteractList[i].TotalInteractTime-InteractList[i].TimeofInteraction)<=ParValue.dT)
                    && InteractList[i].CheckedForDeath==0)
            {
                if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=ParValue.ProbDie_Interact_Hyp1))
                {

                    somatic[s].ApoptoticState=1;
                }
                InteractList[i].CheckedForDeath=1;
            }
        }
    }

    /*
     * All infected cells are checked to see if they'll die
     * If time to die for an infected cell is over, then the infected cell disappers
     * The CTl that was interacting is assigned a new velocity
     * This is true for all hypothesis
     * */
    for (int i=0;i<NumSom;i++)
    {

        if ((somatic[i].DyingSince>(somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0)))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;

            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;

            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;


            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {
                        int c=InteractList[j].CTL_Index;


                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);
                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);

                        somatic[i].Numofcontacts--;

                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        //ASK PHILIPPE
                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;
                        DeleteList.push_back(j);
                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}


/*
 * Damage rate based on duration of interaction- exponential damage
 * When damage equals 1, the decision to die is taken
 * */
void hypothesis2_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;

    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        if ((somatic[s].ApoptoticState==0))
        {
            somatic[s].Damage+=((ParValue.AlphaDamage_Hyp2)*ParValue.dT*exp(-1*(InteractList[i].TimeofInteraction)*(ParValue.AlphaDamage_Hyp2)));
            if (somatic[s].Damage>=1.0)
            {
                somatic[s].ApoptoticState=1;
            }
        }
    }

    //Apoptotic state over- cell dead!!!!!!!!!!!!
    for (int i=0;i<NumSom;i++)
    {

        if ((somatic[i].DyingSince>=somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;
            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;

            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {


                        int c=InteractList[j].CTL_Index;

                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);



                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);

                        somatic[i].Numofcontacts--;

                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        //Check???
                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;


                        DeleteList.push_back(j);
                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * Damage and repair: Exponential damage and repair rate is proportional to damage
 * Repair can only take place when infected cell is not interacting
 * */
void hypothesis3_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;



    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        if(ParValue.DeathDecisionAtStart==1)
        {
            if(InteractList[i].CheckedForDeath==0)// && somatic[s].CheckedDeath==0)
            {
                float DeathProb_s=exp(-1*(somatic[s].PriorCellContacts)/(ParValue.AlphaDamage_Hyp3));

                //cout<<somatic[s].PriorCellContacts<<"zz "<<ParValue.AlphaDamage_Hyp3<<" zz"<<DeathProb_s<<endl;
                //float DeathProb_s=(somatic[s].PriorCellContacts*ParValue.GammaFactor_Hyp6);
                if (DeathProb_s<0.0)
                    DeathProb_s=0.0;
                else if (DeathProb_s>1.0)
                    DeathProb_s=1.0;

                if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=DeathProb_s))
                {
                    somatic[s].ApoptoticState=1;
                }
                InteractList[i].CheckedForDeath=1;

                //somatic[s].CheckedDeath=1;
            }
        }
        else if(ParValue.DeathDecisionAtStart==0)
        {
            if ((somatic[s].ApoptoticState==0 && (InteractList[i].TotalInteractTime-InteractList[i].TimeofInteraction)<=ParValue.dT)
                    && InteractList[i].CheckedForDeath==0)
            {
                float DeathProb_s=exp(-1*(somatic[s].PriorCellContacts)/(ParValue.AlphaDamage_Hyp3));

                //cout<<somatic[s].PriorCellContacts<<"zz "<<ParValue.AlphaDamage_Hyp3<<" zz"<<DeathProb_s<<endl;
                //float DeathProb_s=(somatic[s].PriorCellContacts*ParValue.GammaFactor_Hyp6);
                if (DeathProb_s<0.0)
                    DeathProb_s=0.0;
                else if (DeathProb_s>1.0)
                    DeathProb_s=1.0;

                if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=DeathProb_s))
                {
                    somatic[s].ApoptoticState=1;
                }
                InteractList[i].CheckedForDeath=1;
            }
        }
    }

    /*
    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        if ((somatic[s].ApoptoticState==0))
        {
            somatic[s].Damage+=((ParValue.AlphaDamage_Hyp3)*ParValue.dT*exp(-1*(InteractList[i].TimeofInteraction)*(ParValue.AlphaDamage_Hyp3)));

            if (somatic[s].Damage>=1.0)
            {
                somatic[s].Damage=10.0;
                somatic[s].ApoptoticState=1;
            }
        }
    }
    */

    //Repair
    /*
    for(int i=0;i<NumSom;i++)
    {
        if(somatic[i].Numofcontacts==0 && somatic[i].Damage!=0.0 && somatic[i].Dead==0)
        {
            somatic[i].NotInteractingSince+=ParValue.dT;
            if (somatic[i].Damage<1e-6)
            {
                somatic[i].Damage=0.0;
            }
            else
            {
                somatic[i].Damage-=(ParValue.BetaRepair_Hyp3*somatic[i].Damage*ParValue.dT);
            }
        }
    }
    */

    //Apoptotic state over- cell dead!!!!!!!!!!!!
    for (int i=0;i<NumSom;i++)
    {

        if ((somatic[i].DyingSince>=somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;
            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;
            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {
                        // Philippe: this is copied-pasted for all hypothesis, bery high danger of problem. Why not making a function ?

                        int c=InteractList[j].CTL_Index;


                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);
                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);


                        somatic[i].Numofcontacts--;

                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;

                        DeleteList.push_back(j);
                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * CTLs get more lethal with increasing number of interactions
 * */
void hypothesis4_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();
    vector <int> DeleteList;
    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        if ((somatic[s].ApoptoticState==0))
        {
            int PriorContact=CTL[c].PriorCellContacts;
            bool ToDamage=true;
            if(ParValue.LongContactsControlledDamage==1 && InteractList[i].TimeofInteraction>ParValue.ThresholdContactDurDamage)
                ToDamage=false;

            if(ToDamage==true)
            {
                somatic[s].Damage+=float(ParValue.AlphaDamage_Hyp4*PriorContact*ParValue.dT);
            }

            if (somatic[s].Damage>=1.0)
            {
                somatic[s].ApoptoticState=1;
            }
        }
    }

    //Apoptotic state over- cell dead!!!!!!!!!!!!
    for (int i=0;i<NumSom;i++)
    {


        if ((somatic[i].DyingSince>=somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;

            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;
            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {

                        int c=InteractList[j].CTL_Index;

                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);
                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);


                        somatic[i].Numofcontacts--;


                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;

                        DeleteList.push_back(j);
                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * Time dependent linear damage
 * */
void hypothesis5_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;



    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;

        if ((somatic[s].ApoptoticState==0))
        {

            //InteractList[i].TimeofInteraction

            bool ToDamage=true;
            if(ParValue.LongContactsControlledDamage==1 && InteractList[i].TimeofInteraction>ParValue.ThresholdContactDurDamage)
                ToDamage=false;

            if(ToDamage==true)
            {
                //cout<<ParValue.ThresholdContactDurDamage<<"  "<<ParValue.AlphaDamage_Hyp5<<endl;

                somatic[s].Damage+=(ParValue.dT*ParValue.AlphaDamage_Hyp5);
            }
            if (somatic[s].Damage>=1.0)
            {
                somatic[s].ApoptoticState=1;
            }

        }
    }

    for (int i=0;i<NumSom;i++)
    {


        if ((somatic[i].DyingSince>=somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;

            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;
            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {

                        int c=InteractList[j].CTL_Index;

                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);
                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);


                        somatic[i].Numofcontacts--;

                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;

                        DeleteList.push_back(j);
                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * Probability of death at each interaction increases linearly as no. of previous interactions goes up
 * Basically, it is directly proportional to no. of previous contacts
 * */
void hypothesis6_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;

    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        if(ParValue.DeathDecisionAtStart==1)
        {
            if(InteractList[i].CheckedForDeath==0)// && somatic[s].CheckedDeath==0)
            {
                float DeathProb_s=(somatic[s].PriorCellContacts*ParValue.GammaFactor_Hyp6);
                if (DeathProb_s<0.0)
                    DeathProb_s=0.0;
                else if (DeathProb_s>1.0)
                    DeathProb_s=1.0;

                if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=DeathProb_s))
                {
                    somatic[s].ApoptoticState=1;
                }
                InteractList[i].CheckedForDeath=1;

                //somatic[s].CheckedDeath=1;
            }
        }
        else if(ParValue.DeathDecisionAtStart==0)
        {
            if ((somatic[s].ApoptoticState==0 && (InteractList[i].TotalInteractTime-InteractList[i].TimeofInteraction)<=ParValue.dT)
                    && InteractList[i].CheckedForDeath==0)
            {
                float DeathProb_s=(somatic[s].PriorCellContacts*ParValue.GammaFactor_Hyp6);
                if (DeathProb_s<0.0)
                    DeathProb_s=0.0;
                else if (DeathProb_s>1.0)
                    DeathProb_s=1.0;

                if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=DeathProb_s))
                {
                    somatic[s].ApoptoticState=1;
                }
                InteractList[i].CheckedForDeath=1;
            }
        }
    }


    for (int i=0;i<NumSom;i++)
    {


        if ((somatic[i].DyingSince>=(somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0)))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;

            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;

            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {



                        int c=InteractList[j].CTL_Index;

                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);


                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);


                        somatic[i].Numofcontacts--;


                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;

                        DeleteList.push_back(j);

                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * Time dependent linear damage and time dependent linear repair
 * */
void hypothesis7_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;



    for(int i=0;i<NumSom;i++)
    {
        if(somatic[i].Numofcontacts==0 && somatic[i].Damage!=0.0 && somatic[i].Dead==0)
        {
            somatic[i].NotInteractingSince+=ParValue.dT;
            if (somatic[i].Damage<1e-6)
            {
                somatic[i].Damage=0.0;
            }
            else
            {

                somatic[i].Damage=somatic[i].Damage- (ParValue.BetaRepair_Hyp7*ParValue.dT*somatic[i].Damage);
                if(somatic[i].Damage<0)
                    somatic[i].Damage=0;
            }
        }
        else if(somatic[i].Numofcontacts!=0)
        {

            bool ToDamage=true;
            /*
            if(ParValue.LongContactsControlledDamage==1 && InteractList[i].TimeofInteraction>ParValue.ThresholdContactDurDamage)
                ToDamage=false;
                */

            if(ToDamage==true)
            {
                somatic[i].Damage= somatic[i].Damage + (somatic[i].Numofcontacts*ParValue.dT*ParValue.AlphaDamage_Hyp7)
                        -(ParValue.BetaRepair_Hyp7*ParValue.dT*somatic[i].Damage);
            }
            else
            {
                somatic[i].Damage=somatic[i].Damage- (ParValue.BetaRepair_Hyp7*ParValue.dT*somatic[i].Damage);
                if(somatic[i].Damage<0)
                    somatic[i].Damage=0;
            }


            if (somatic[i].Damage>=1.0)
            {

                int c=-1;
                for(int j=0;j<NumInt;j++)
                {
                    if(i==InteractList[j].Somatic_Index)
                    {
                        somatic[i].ApoptoticState=1;
                        c=InteractList[j].CTL_Index;
                    }
                }
            }
        }
    }

    for (int i=0;i<NumSom;i++)
    {


        if ((somatic[i].DyingSince>=somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;

            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;
            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {



                        int c=InteractList[j].CTL_Index;

                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);
                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);


                        somatic[i].Numofcontacts--;

                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        //Check???
                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;

                        DeleteList.push_back(j);
                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }
}

/*
 * Non uniform T cell killing capacity
 * */
void hypothesis8_death(vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL, sim_parameters & ParValue,
                       vector<Interact> &InteractList, int &DeadCells, int &Dead_TimeBin, int &Dead_TimeBinAfterAT, int &Dead_TimeBinSS, float t, vector<int> &DurationInt, vector<CTLHistory> &TrumpHistory,
                       vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    int NumInt=InteractList.size();
    int NumSom=somatic.size();

    vector <int> DeleteList;

    for(int i=0;i<NumInt;i++)
    {
        int s=InteractList[i].Somatic_Index;
        int c=InteractList[i].CTL_Index;
        if(InteractList[i].CheckedForDeath==0)// && somatic[s].CheckedDeath==0)
        {
            /// Philippe: this will not work, to sample always a new number for the T cell
            //float DeathProb_s=min(1.0,max(0., (double) (CTL[c].KillingCapacity*random::normal(ParValue.KillingFactorMean_Hyp8, ParValue.KillingFactorSD_Hyp8))));

            if((somatic[s].ApoptoticState==0) && (random::uniformDouble(0.0, 1.0)<=CTL[c].KillingCapacity))
            {
                somatic[s].ApoptoticState=1;
            }
            InteractList[i].CheckedForDeath=1;
            //somatic[s].CheckedDeath=1;
        }
    }

    for (int i=0;i<NumSom;i++)
    {


        if ((somatic[i].DyingSince>=(somatic[i].TimeTakenToDie) && (somatic[i].ApoptoticState==1) && (somatic[i].Dead==0)))
        {
            somatic[i].Dead=1;
            somatic[i].DiedAtTime=t;
            DeadCells++;
            Dead_TimeBin++;
            if(t>ParValue.AnalyzeAfter)
                Dead_TimeBinAfterAT++;

            if(Inside(somatic[i].Posn, ParValue)==1)
                Dead_TimeBinSS++;

            while(somatic[i].Numofcontacts!=0)
            {
                for(int j=0;j<NumInt;j++)
                {
                    if(InteractList[j].Somatic_Index==i)
                    {



                        int c=InteractList[j].CTL_Index;

                        somatic[i].SingleContactDuration.push_back(CompleteHistory[c].OngoingContactTime);
                        int ST_HistoryPosn= CTL[c].ST_HistoryPosn;
                        if(t>ParValue.AnalyzeAfter)
                            somatic[i].SingleContactDuration_AfterAT.push_back(ShortTimeHistory[ST_HistoryPosn].OngoingContactTime);


                        AssignInitialVelocityAndAngle(CTL[c], ParValue, 0);
                        UpdateHistory(CTL, ParValue, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, c, t);


                        somatic[i].Numofcontacts--;


                        CTL[c].Numofcontacts--;

                        CTL[c].RemainingLatentTime=0;

                        int RangePosn=int (InteractList[j].TimeofInteraction);
                        DurationInt[RangePosn]++;

                        DeleteList.push_back(j);

                    }
                }
            }
        }
        else if (somatic[i].ApoptoticState==1 && somatic[i].DyingSince<=somatic[i].TimeTakenToDie)
        {
            somatic[i].DyingSince+=ParValue.dT;
        }
    }

    std::sort(DeleteList.begin(), DeleteList.end());

    for(int ii=(DeleteList.size()-1);ii>=0;ii--)
    {
        int d=DeleteList[ii];
        InteractList.erase(InteractList.begin()+d);
    }

}


