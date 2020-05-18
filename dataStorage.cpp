#include <iostream>
#include <stdlib.h>
#include <math.h>

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
#include "dataStorage.h"
#include <QString>

using namespace std;

/*
 * Storing history of CTLs after each movement
 * For physical window, CTLs whose LocationStatus is MovedIn are assigned a new history from scratch
 * When LocationStatus is MovedOut, the history ends.
 * The above steps are done to ensure that it mimics the experimental setup as closely as possible
 * */
void StoreCTLMoveHistory(vector<CTLCellType> &CTL, sim_parameters & ParValue, float t, vector<CTLHistory> &TrumpHistory,
                         vector <CTLHistory> &ShortTimeHistory, vector <CTLHistory> &ShortSpaceHistory, vector<CTLHistory> &CompleteHistory, vector<Interact> &InteractList,
                         vector<SomaticCellType> & somatic)
{
    // starting point of analysis
    if((t-ParValue.AnalyzeAfter)>=0 && (t-ParValue.AnalyzeAfter)<ParValue.dT)
    {
        for (int i=0; i < (int) CTL.size(); ++i)
        {

            CTLHistory ShortTimeTemp;
            ShortTimeTemp.CTLIndex=i;
            ShortTimeTemp.InteractedWith=CTL[i].Numofcontacts;
            ShortTimeTemp.Killed=0;
            ShortTimeTemp.Ongoing=1;
            ShortTimeTemp.DyingCellsInteract=0;
            ShortTimeTemp.HistoryStart=t;

            if(CTL[i].Numofcontacts!=0)
            {
                for(int find=0; find< (int) InteractList.size(); ++find)
                {
                    int SomInd=InteractList[find].Somatic_Index;
                    if((InteractList[find].CTL_Index==i) && (somatic[SomInd].ApoptoticState==1))
                    {
                        ShortTimeTemp.DyingCellsInteract=1;
                    }
                }
            }


            CTL[i].ST_HistoryPosn=i;
            ShortTimeHistory.push_back(ShortTimeTemp);

            //ALERT: cout<<"Still have to change this. On reaching analyzeTime, observation starts. What if DyingCellsInteract is not zero?"<<endl;
            bool location= Inside(CTL[i].Posn, ParValue);
            if (location==1)
            {
                CTLHistory TrumpH;
                TrumpH.CTLIndex=i;
                TrumpH.Killed=0;
                TrumpH.Ongoing=1;
                TrumpH.InteractedWith=CTL[i].Numofcontacts;
                TrumpH.DyingCellsInteract=0;
                TrumpH.HistoryStart=t;
                if(CTL[i].Numofcontacts!=0)
                {
                    for(unsigned int find=0;find<InteractList.size();find++)
                    {
                        int SomInd=InteractList[find].Somatic_Index;
                        if(InteractList[find].CTL_Index==i && somatic[SomInd].ApoptoticState==1)
                        {
                            TrumpH.DyingCellsInteract=1;
                        }
                    }
                }


                CTL[i].Trump_HistoryPosn=TrumpHistory.size();
                TrumpHistory.push_back(TrumpH);
            }
        }
    }
    else if(t>ParValue.AnalyzeAfter)
    {
        for (unsigned int i=0;i<CTL.size();i++)
        {
            if (CTL[i].LocationStatus==MovedOut)
            {
                int tempHistoryPosn=CTL[i].Trump_HistoryPosn;

                TrumpHistory[tempHistoryPosn].Ongoing=0;
                TrumpHistory[tempHistoryPosn].HistoryEnd=t;
                CTL[i].Trump_HistoryPosn=-2;
            }
            else if(CTL[i].LocationStatus==MovedIn)
            {
                CTLHistory H;
                H.CTLIndex=i;
                H.InteractedWith=0;
                H.Killed=0;
                H.Ongoing=1;
                H.HistoryStart=t;
                CTL[i].Trump_HistoryPosn=TrumpHistory.size();
                TrumpHistory.push_back(H);
            }
            else if (CTL[i].LocationStatus==StayedOut)
            {
                CTL[i].Trump_HistoryPosn=-3;
            }
        }
    }

    for (unsigned int i=0;i<CTL.size();i++)
    {
        if (CTL[i].LocationStatus==MovedOut)
        {
            int temp_SSHistPosn=CTL[i].SS_HistoryPosn;
            ShortSpaceHistory[temp_SSHistPosn].Ongoing=0;
            ShortSpaceHistory[temp_SSHistPosn].HistoryEnd=t;
            CTL[i].SS_HistoryPosn=-2;
        }
        else if(CTL[i].LocationStatus==MovedIn)
        {
            CTLHistory H;
            H.CTLIndex=i;
            H.InteractedWith=0;
            H.Killed=0;
            H.Ongoing=1;
            H.HistoryStart=t;
            CTL[i].SS_HistoryPosn=ShortSpaceHistory.size();
            ShortSpaceHistory.push_back(H);
        }
        else if (CTL[i].LocationStatus==StayedOut)
        {
            CTL[i].SS_HistoryPosn=-3;
        }
    }

    for(unsigned int i=0;i<ShortSpaceHistory.size();i++)
    {
        int ID=ShortSpaceHistory[i].CTLIndex;
        if(ShortSpaceHistory[i].Ongoing==1 && CTL[ID].Numofcontacts!=0)
            ShortSpaceHistory[i].TotalContactDuration+=ParValue.dT;
    }

    for(unsigned int i=0;i<CompleteHistory.size();i++)
    {
        int ID=CompleteHistory[i].CTLIndex;
        if(CompleteHistory[i].Ongoing==1 && CTL[ID].Numofcontacts!=0)
            CompleteHistory[i].TotalContactDuration+=ParValue.dT;
    }
    if(t>ParValue.AnalyzeAfter)
    {
        for(unsigned int i=0;i<ShortTimeHistory.size();i++)
        {
            int ID=ShortTimeHistory[i].CTLIndex;
            if(ShortTimeHistory[i].Ongoing==1 && CTL[ID].Numofcontacts!=0)
                ShortTimeHistory[i].TotalContactDuration+=ParValue.dT;
        }

        for(unsigned int i=0;i<TrumpHistory.size();i++)
        {
            int ID=TrumpHistory[i].CTLIndex;
            if(TrumpHistory[i].Ongoing==1 && CTL[ID].Numofcontacts!=0)
                TrumpHistory[i].TotalContactDuration+=ParValue.dT;
        }
    }
}



void UpdateHistory(vector<CTLCellType> &CTL, sim_parameters &ParValue, vector<CTLHistory> &TrumpHistory,
                   vector<CTLHistory> &ShortTimeHistory, vector<CTLHistory> &ShortSpaceHistory,
                   vector <CTLHistory> &CompleteHistory, int CTLIndex, float t)
{
    CompleteHistory[CTLIndex].Killed++;
    CompleteHistory[CTLIndex].OngoingContactTime=0.0;

    int SS_Posn= CTL[CTLIndex].SS_HistoryPosn;
    int ST_Posn= CTL[CTLIndex].ST_HistoryPosn;
    int Trump_Posn= CTL[CTLIndex].Trump_HistoryPosn;

    if(t>ParValue.AnalyzeAfter)
    {
        //random check
        if(CTLIndex!= ST_Posn)
            cout<<"Something is going wrong while storing histories!!!!!!!!!!!!!"<<endl;

        ShortTimeHistory[ST_Posn].Killed++;
        ShortTimeHistory[ST_Posn].OngoingContactTime=0.0;
    }

    if(SS_Posn>=0)
    {
        ShortSpaceHistory[SS_Posn].Killed++;
        ShortSpaceHistory[SS_Posn].OngoingContactTime=0.0;
    }

    if(Trump_Posn>=0 && t>ParValue.AnalyzeAfter)
    {
        TrumpHistory[Trump_Posn].Killed++;
        TrumpHistory[Trump_Posn].OngoingContactTime=0.0;
    }
}
