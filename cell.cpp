#include <iostream>
#include <stdlib.h>
#include <math.h>

//#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include <vector>
#include "random.h"
#include "cell.h"
#include <math.h>
#include <signal.h>  // for raise(SIGSEGV);
#include "main.h"
#define PI 3.14159265
using namespace std;


/*
 * Parent cell class
 * */
cell:: cell()
{
    Posn[0]=0; ///x coordinate- Current posn
    Posn[1]=0; ///y coordinate- Current posn
    Posn[2]=0; ///z coordinate- Current posn
    Numofcontacts=0; ///No. of existing contacts
    PriorCellContacts=0; ///Total number of contacts formed over time
    NotInteractingSince=0.0; ///Time since cell has had no interaction
    StoppedInteracting=0.0;
}

/*
 * CTL class
 * */
CTLCellType:: CTLCellType()
{
    PrevPosn[0]=0; ///x coordinate- Current posn
    PrevPosn[1]=0; ///y coordinate- Current posn
    PrevPosn[2]=0; ///z coordinate- Current posn
    Time_SameDir=0.0; ///For what time has the direction been same (Persistent time)
    velocity_mag = 0;
    velocity[0]=0; ///x vector of velocity
    velocity[1]=0; ///y vector of velocity
    velocity[2]=0; ///z vector of velocity
    theta=0; /// Angles to define velocity
    phi=0; /// Angles to define velocity
    RemainingLatentTime=0.0; ///Duration since when CTL hasn't been interacting
    Trump_HistoryPosn=-1; ///Posn in vector TrumpHistory where history for that cell is being stored
    SS_HistoryPosn=-1; ///Posn in vector ShortSpaceHistory where history for that cell is being stored
    ST_HistoryPosn=-1; ///Posn in vector ShortTimeHistory where history for that cell is being stored
    Complete_HistoryPosn=-1; ///Posn in vector CompleteHistory where history for that cell is being stored
    LocationStatus=-1; ///Calculated wrt physical window. Can be: MovedIn, MovedOut, StayedIn, StayedOut
    KillingCapacity=0.0;


}

/*
 * infected cell class
 * */
SomaticCellType::SomaticCellType()
{

    DurationFlux=0.0;
    ApoptoticState=0; ///To store if cell is in apoptotic state
    DyingSince=0.0; ///Time since cell has been in apoptotic state
    Dead=0; ///Alive or dead
    Damage=0.0; ///Damage value of infected cell
    ContactsAfterAnalyzeTime=0; ///Contacts observed during time window
    ZombieContacts=0;

    Contacts0_240=0;
    Contacts0_360=0;
    Contacts0_480=0;


    CheckedDeath=0; ///During each interaction, cell should only be checked for death once
    interactionDuringDeath=0; ///No. of contacts formed during apoptotic state
    GrandTotalInteractTime=0.0; ///Total duration of all contacts that the cell has had
    GrandTotalInteractTimeAfterAT=0.0; ///Total duration of all contacts that the cell has had after analyze time
    DiedAtTime=0.0;
    Data4I=0.0;
    Data4I_AfterAT=0.0;
    TimeTakenToDie=0.0;

}


CTLHistory::CTLHistory()
{
    CTLIndex=-1; ///Posn in CTL vector whose history is being stored
    Ongoing=1; ///If cell moves out of physical observation window, ongoing =0 and storing of history stops
    InteractedWith=0; ///No. on infected cells seen
    Killed=0; ///No. of infected cells that died during this cell was interacting with it
    DyingCellsInteract=0; ///No. of interactions with infected cells in apoptotic state

    OngoingContactTime = 0;
    Interacting=0;
    TotalContactDuration=0.0;
    HistoryStart=0.0;
    HistoryEnd=0.0;

}

enum {CTLTYPE=1, SOMATYPE=2};

/*
 * New random velocities need to be assigned to CTLs at various points-
 * initialization, breaking off interaction or when infected cell dies
 * */
void AssignInitialVelocityAndAngle(CTLCellType &currentCTL, sim_parameters & ParValue, bool randomTimeSameDir)
//        float vel[], sim_parameters ParValue, float &theta, float &phi)
{
    currentCTL.theta =random::uniformDouble(0.0, 360.0)-180;
    currentCTL.phi =random::uniformDouble(0.0, 360.0);

    do
    {
        currentCTL.velocity_mag = ParValue.DistributionSpeed.getRandValue();
    }while(currentCTL.velocity_mag<0);
    currentCTL.velocity[0] = currentCTL.velocity_mag* sin((currentCTL.theta*PI)/180)* cos(currentCTL.phi*PI/180);
    currentCTL.velocity[1] = currentCTL.velocity_mag* sin(currentCTL.theta*PI/180)* sin(currentCTL.phi*PI/180);
    currentCTL.velocity[2] = currentCTL.velocity_mag* cos(currentCTL.theta*PI/180);

    // Philippe: If you don't want the cells to be synchronized at start, would be smart to put a different
    // time before persistent time. But when the persistent time is from a distribution, would be peter to have
    // a field 'time_remaining_sameDir' instead ...
    // so this solution next line is bad, you will need to improve it
    //Ananya: I have kept this for now. I am confused. Should I make it completely random? The persistent time?
    if(randomTimeSameDir)
    {
        currentCTL.Time_SameDir=random::uniformDouble(0, ParValue.DistributionPersistentTime.getRandValue());
    }
    else
        currentCTL.Time_SameDir = 0.;
}

void GetNewVelocityAndAngle(CTLCellType & currentCTL, sim_parameters & ParValue){
    float theta_new;
    do{
        theta_new= (PI/180.) * ParValue.DistributionTurningAngles.getRandValue();
    }while(theta_new<0);

    float phi_new= (PI/180.) * random::uniformDouble(0.0, 360.0);
    float trans_vector[3];
    trans_vector[0]=sin(theta_new)*cos(phi_new);
    trans_vector[1]=sin(theta_new)*sin(phi_new);
    trans_vector[2]=cos(theta_new);

    transform(trans_vector, currentCTL.theta, currentCTL.phi);

    if (trans_vector[0] > 0.)
        currentCTL.phi= atan(trans_vector[1]/trans_vector[0])*180/PI;
    else
        currentCTL.phi= 180 + atan(trans_vector[1]/trans_vector[0])*180/PI;

    currentCTL.theta= acos(trans_vector[2])*180/PI;

    do
    {
        currentCTL.velocity_mag = ParValue.DistributionSpeed.getRandValue();
    }while(currentCTL.velocity_mag<0);

    currentCTL.velocity[0]=currentCTL.velocity_mag * trans_vector[0];
    currentCTL.velocity[1]=currentCTL.velocity_mag * trans_vector[1];
    currentCTL.velocity[2]=currentCTL.velocity_mag * trans_vector[2];

    currentCTL.Time_SameDir=0;
}



/*
 * Function to calculate perpendicular distance between a point and a line
 * */
float DistPoint_Vector(float OriginalPosn[], float NewPosn[], float DistFrom[], float denominator)
{

    float dist_array[3]={0}, var1, var3, DistFromVector=0, dotproduct=0;

    for(int k=0;k<3;k++)
        dotproduct+=(DistFrom[k]-OriginalPosn[k])*(NewPosn[k]-OriginalPosn[k]);

    dotproduct=dotproduct/(denominator*denominator);

    for (int k=0;k<3;k++)
    {
        var1=DistFrom[k]-OriginalPosn[k];

        var3=(NewPosn[k]-OriginalPosn[k]);
        dist_array[k]=var1-(dotproduct*var3);
        DistFromVector+=pow(dist_array[k],2);
    }

    DistFromVector=pow (DistFromVector, 0.5);

    return DistFromVector;
}

/*
 * Function to assign positions to cells. Works for both CTL and somatic. There should be no overlap between cells
 * To decide which vector (somatic or CTL) the cell has to be put in, an identifying nymber X is sent
 * It also stores history of CTLs for CompleteHistory and ShortSpaceHistory. (Both have histores stored right from start of simulation)
 * */
void AddCells(int num_cell, int X, vector<SomaticCellType> &somatic, vector<CTLCellType> &CTL,
              sim_parameters & ParValue, vector<CTLHistory> &CompleteHistory, vector<CTLHistory> &ShortSpaceHistory)
{
    if(X==SOMATYPE)
    {
        for (int i=0;i<num_cell;i++)
        {
            SomaticCellType temp;
            float LowerLimitPosn=ParValue.box_Z*ParValue.factorConfinementInfected;
            int FoundPosn=0;
            while(FoundPosn==0)
            {
                temp.Posn[0]=random::uniformDouble(0.0, ParValue.box_X);
                temp.Posn[1]=random::uniformDouble(0.0, ParValue.box_Y);
                temp.Posn[2]=random::uniformDouble(LowerLimitPosn, ParValue.box_Z);
                int NoProblem=0;
                for (int k=0;k<i;k++)
                {
                    if(CalcDist(temp.Posn, somatic[k].Posn)>=(ParValue.SomaticRadius*2))
                        NoProblem++;
                }

                if (CTL.size()!=0)
                {
                    for (int k=0;k< (int) CTL.size();k++)
                    {
                        if(CalcDist(temp.Posn,CTL[k].Posn)>=(ParValue.CTLRadius+ParValue.SomaticRadius))
                            NoProblem++;
                    }
                }
                if (NoProblem== int(i+CTL.size()))
                    FoundPosn=1;
            }
            float TTD=-1;

            do{
                TTD= ParValue.DistributionTimeToDie.getRandValue(); //TimeTodie
            }while(TTD<0);
            temp.TimeTakenToDie=TTD;

            somatic.push_back(temp);
        }
    }
    else if (X==CTLTYPE)
    {
        for (int i=0;i<num_cell;i++)
        {
            CTLCellType temp;

            temp.KillingCapacity=min(1.0, max(0., random::normal(ParValue.KillingFactorMean_Hyp8, ParValue.KillingFactorSD_Hyp8)));


            AssignInitialVelocityAndAngle(temp, ParValue, 1);

            bool FoundPosn=false;
            while(!FoundPosn)
            {
                temp.Posn[0]=random::uniformDouble(0.0, ParValue.box_X);
                temp.Posn[1]=random::uniformDouble(0.0, ParValue.box_Y);
                temp.Posn[2]=random::uniformDouble(0.0, ParValue.box_Z);

                int NoProblem=0;
                for (int k=0;k<i;k++)
                {
                    if(CalcDist(temp.Posn, CTL[k].Posn)>=(ParValue.CTLRadius*2))
                        NoProblem++;
                }

                if (somatic.size()!=0)
                {
                    for (int k=0;k< (int) somatic.size();k++)
                    {
                        if(CalcDist(temp.Posn,somatic[k].Posn)>=(ParValue.CTLRadius+ParValue.SomaticRadius))
                            NoProblem++;
                    }
                }

                if (NoProblem==(i+ (int) somatic.size()))
                    FoundPosn=true;
            }
            CTL.push_back(temp);
        }

        for (int i=0; i < (int) CTL.size(); ++i)
        {
            CTLHistory TempH, ShortSpaceTempH;
            TempH.CTLIndex=i;
            TempH.InteractedWith=0;
            TempH.Ongoing=1;
            TempH.Killed=0;
            CTL[i].Complete_HistoryPosn=i;
            CompleteHistory.push_back(TempH);

            if(Inside(CTL[i].Posn, ParValue)==1)
            {
                ShortSpaceTempH.CTLIndex=i;
                ShortSpaceTempH.InteractedWith=0;
                ShortSpaceTempH.Ongoing=1;
                ShortSpaceTempH.Killed=0;
                ShortSpaceTempH.HistoryStart=0.0;
                CTL[i].SS_HistoryPosn=ShortSpaceHistory.size();
                ShortSpaceHistory.push_back(ShortSpaceTempH);
            }
        }
    }
}

/*
 * Function to make cells move- The cells have an assigned direction.
 * After persistent time, cells change direction. Turning angle taken from data sent by S. Halle
 * Collsion detection to avoid cells 'moving through each other'
 * The CTLs are moved in random order to avoid bias.
 * */
void CTLmove(vector<CTLCellType> &CTL, vector<SomaticCellType> &somatic, sim_parameters & ParValue, float t)
{
    int num_CTL_temp=CTL.size();
    int num_somatic_temp=somatic.size();

    vector <int> CTL_Order(num_CTL_temp);
    std::iota(CTL_Order.begin(), CTL_Order.end(), 0);   // fills with 0, 1, 2, 3 ...
    random::shuffle(CTL_Order);

    for(int i=0; i <num_CTL_temp; ++i)
    {
        /*
         * If CTL has no contact and has already had at least one contact (To ensure that we are actually calculating time between contacts)
         * */
        if(CTL[i].Numofcontacts==0 && CTL[i].PriorCellContacts!=0)
            CTL[i].NotInteractingSince+=ParValue.dT;
    }


    //vector<int> EmptyCTL(num_CTL_temp,0);
    for (int j=0; j < num_CTL_temp; ++j)
    {
        // Philippe: here, you will update the same cell multiple times
        //Ananya: Reply to above statement. It works already. All cells are updated and done so in a non uniform manner. I like it.
        // int U=random::uniformInteger(0, (CTL_Order.size()-1)); ///Cell to be updated in that cycle
        // int Update=CTL_Order[U];
        int Update=CTL_Order[j];
        //EmptyCTL[Update]++;


        for(int ii=0;ii<3;ii++)
        {
            CTL[Update].PrevPosn[ii]=CTL[Update].Posn[ii];
        }

        if(CTL[Update].Numofcontacts==0)
        {
            CTL[Update].RemainingLatentTime+=ParValue.dT;

            if(CTL[Update].Time_SameDir>=ParValue.PersistentTimeMean)
            {
                GetNewVelocityAndAngle(CTL[Update], ParValue);
            }
            else // means CTL[Update].Numofcontacts == 0 and Time_SameDir < ParValue.PersistentTime
            {
                CTL[Update].Time_SameDir+=ParValue.dT;
            }
        }
        else // means (CTL[Update].Numofcontacts!=0)
        {
            CTL[Update].velocity_mag = 0;
            CTL[Update].velocity[0]= 0.0;
            CTL[Update].velocity[1]= 0.0;
            CTL[Update].velocity[2]= 0.0;
        }

        float TrialPosn[3]={0};
        for (int k=0;k<3;k++)
            TrialPosn[k]=CTL[Update].Posn[k] +(ParValue.dT* CTL[Update].velocity[k]);

        while (TrialPosn[1]>=ParValue.box_Y)
        {
            AssignInitialVelocityAndAngle(CTL[Update], ParValue, 0);
            for (int k=0;k<3;k++)
            {
                TrialPosn[k]=CTL[Update].Posn[k] +(ParValue.dT* CTL[Update].velocity[k]);
            }
        }

        /*
         * Geometry approach - collision detection
         * In the direction of the velocity, the maximum distance less than (velocity*dT) is found
         * This dist is such that the cell can move without intersecting/colliding with any other cell
         * */
        char CellType='X'; ///To know if the closest cell is a CTL or somatic
        int ClosestCell=-1; ///CC (to store ID)

        float time_max=ParValue.dT;

        if(CTL[Update].Numofcontacts==0)
        {
            for (int z=0;z<num_CTL_temp;z++)
            {
                float diagonal=(ParValue.CTLRadius*2)*ParValue.ThresholdMove;
                if (z!=Update)
                {
                    float denominator=CTL[Update].velocity_mag*ParValue.dT;

                    float dotproduct=0;
                    for(int k=0;k<3;k++)
                        dotproduct+=((CTL[z].Posn[k]-CTL[Update].Posn[k])*(TrialPosn[k]-CTL[Update].Posn[k]));

                    if(dotproduct>0.0)
                    {
                        float DistFromVector=DistPoint_Vector(CTL[Update].Posn, TrialPosn, CTL[z].Posn, denominator);
                        if(DistFromVector<(ParValue.CTLRadius*2*ParValue.ThresholdMove))
                        {
                            float Dist_Cells= CalcDist(CTL[z].Posn, CTL[Update].Posn);
                            float DistAlongVector=pow(Dist_Cells,2)-pow(DistFromVector, 2);
                            DistAlongVector=std::abs(pow(DistAlongVector,0.5));

                            if (DistAlongVector<=(CTL[Update].velocity_mag*ParValue.dT))
                            {
                                float var=std::abs(pow((pow(diagonal, 2)- pow(DistFromVector,2)),0.5));

                                float DistAway=DistAlongVector-var;
                                float time=DistAway/CTL[Update].velocity_mag;
                                if (time_max>time)
                                {
                                    time_max=time;
                                    CellType='C';
                                    ClosestCell=z;
                                }
                            }
                        }
                    }
                }
            }

            for (int z=0;z<num_somatic_temp;z++)
            {
                if(somatic[z].Dead ==0)
                {
                    float diagonal=(ParValue.CTLRadius+ParValue.SomaticRadius)*ParValue.ThresholdMove;
                    float denominator=0, Dist_Cells=0, DistAlongVector=0;
                    denominator=CTL[Update].velocity_mag*ParValue.dT;

                    float dotproduct=0;
                    for(int k=0;k<3;k++)
                        dotproduct+=((somatic[z].Posn[k]-CTL[Update].Posn[k])*(TrialPosn[k]-CTL[Update].Posn[k]));

                    if(dotproduct>0)
                    {
                        float DistFromVector=DistPoint_Vector(CTL[Update].Posn, TrialPosn, somatic[z].Posn, denominator);
                        if(DistFromVector<((ParValue.CTLRadius+ParValue.SomaticRadius)*ParValue.ThresholdMove))
                        {
                            Dist_Cells= CalcDist(somatic[z].Posn, CTL[Update].Posn);
                            DistAlongVector=pow(Dist_Cells,2)-pow(DistFromVector, 2);

                            DistAlongVector=std::abs(pow(DistAlongVector,0.5));
                            if (DistAlongVector<=((CTL[Update].velocity_mag*ParValue.dT)+(2*ParValue.CTLRadius)))
                            {
                                float var=std::abs(pow((pow(diagonal, 2)- pow(DistFromVector,2)),0.5));
                                float DistAway=DistAlongVector-var;
                                float time=DistAway/CTL[Update].velocity_mag;
                                if (time_max>time)
                                {
                                    time_max=time;
                                    CellType='S';
                                    ClosestCell=z;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (CTL[Update].Numofcontacts!=0)
            time_max=0.0;

        float dimension[3]={ParValue.box_X, ParValue.box_Y, ParValue.box_Z};
        for (int k=0;k<3;k++)
        {
            TrialPosn[k]=CTL[Update].Posn[k]+(CTL[Update].velocity[k]*time_max);

            if (TrialPosn[k]>dimension[k])
            {
                TrialPosn[k]-=dimension[k];
            }
            else if (TrialPosn[k]<0.0)
            {
                TrialPosn[k]+=dimension[k];
            }
            CTL[Update].Posn[k]=TrialPosn[k];
        }
        //CTL_Order.erase(CTL_Order.begin()+U);
    }

    for (unsigned int i=0;i<CTL.size();i++)
        CTL[i].LocationStatus=CellMovementStatus(CTL, ParValue, i);

    /*
     * Updating damage of infected cell
     * */

    float t_temp=float(int(t));
    if(abs(t_temp-t)<1e-6)
    {
        for(unsigned int i=0;i<somatic.size();i++)
            somatic[i].DamageOverTime.push_back(somatic[i].Damage);
    }
}

void StoreSomaticHistory(vector <SomaticCellType> &somatic, sim_parameters &ParValue, float t)
{

    for(unsigned int i=0;i<somatic.size();i++)
    {
        if(somatic[i].PriorCellContacts>0 && somatic[i].Dead==0)
        {
            somatic[i].Data4I+=ParValue.dT;


        }
        if((t-ParValue.AnalyzeAfter)>=0 && (t-ParValue.AnalyzeAfter)<ParValue.dT)
        {
            if(somatic[i].Numofcontacts!=0)
            {
                somatic[i].ContactsAfterAnalyzeTime+=somatic[i].Numofcontacts;
            }
        }

        if(t>ParValue.AnalyzeAfter)
        {
            if(somatic[i].ContactsAfterAnalyzeTime>0 && somatic[i].Dead==0)
                somatic[i].Data4I_AfterAT+=ParValue.dT;
            //else if(somatic[i].DyingSince>0)
                //somatic[i].Data4I_AfterAT+=ParValue.dT;
        }



        //Analyze after 120


        /*
        if((t-120)>=0 && (t-120)<ParValue.dT)
        {
            if(somatic[i].Numofcontacts!=0)
            {
                somatic[i].Contacts120_360+=somatic[i].Numofcontacts;
            }
        }




        if((t-240)>=0 && (t-240)<ParValue.dT)
        {
            if(somatic[i].Numofcontacts!=0)
            {
                somatic[i].Contacts240_480+=somatic[i].Numofcontacts;
            }
        }
        */

        /*
        if((t-180)>=0 && (t-180)<ParValue.dT)
        {
            if(somatic[i].Numofcontacts!=0)
            {
                somatic[i].Contacts180+=somatic[i].Numofcontacts;
            }
        }
        */



    }
}


/*
 * Got this function from Jaber
 * Sends a unit vector in direction that we want it to be turned in along with old theta and phi
 * Returns new unit vector in new direction
 * */
void transform(float trans_vector[], float theta, float phi)
{
    double	x	= trans_vector[0];
    double	y	= trans_vector[1];
    double	z	= trans_vector[2];

    double	X1	=  x*cos(theta*PI/180)+z*sin(theta*PI/180);
    double	Y1	=  y;
    double	Z1	= -x*sin(theta*PI/180)+z*cos(theta*PI/180);

    double	X2	=  X1*cos(phi*PI/180)-Y1*sin(phi*PI/180);
    double	Y2	=  X1*sin(phi*PI/180)+Y1*cos(phi*PI/180);
    double	Z2	=  Z1;

    trans_vector[0]=X2;
    trans_vector[1]=Y2;
    trans_vector[2]=Z2;
}


/*
 * Function to calculate whether CTL moved in or moved out or stayed in or stayed out of the physical window
 * */
int CellMovementStatus(vector <CTLCellType> & CTL, sim_parameters & ParValue, int i)
{
    int CellLocationStatus = StayedIn; // basal value
    bool InsideBeforeMoving, InsideAfterMoving;
    InsideBeforeMoving= Inside(CTL[i].PrevPosn, ParValue);
    InsideAfterMoving= Inside(CTL[i].Posn, ParValue);

    if(InsideBeforeMoving== 1 && InsideAfterMoving== 1 )
        CellLocationStatus= StayedIn;
    else if (InsideBeforeMoving== 0 && InsideAfterMoving== 0)
        CellLocationStatus=StayedOut;
    else if (InsideBeforeMoving== 1 && InsideAfterMoving== 0)
        CellLocationStatus=MovedOut;
    else if (InsideBeforeMoving== 0 && InsideAfterMoving== 1)
        CellLocationStatus=MovedIn;
    else
        exit(0);

    return CellLocationStatus;
}

/*
 * Function to check if CTL is inside or outside of physical window
 * */
bool Inside(float CheckingPosn[], sim_parameters & ParValue)
{
    float LowerLim[3], UpperLim[3];
    int Counter=0;

    LowerLim[0]= ParValue.X_LowerLim;
    LowerLim[1]= ParValue.Y_LowerLim;
    LowerLim[2]= ParValue.Z_LowerLim;
    UpperLim[0]= LowerLim[0] + ParValue.X_AnalyzeRange;
    UpperLim[1]= LowerLim[1] + ParValue.Y_AnalyzeRange;
    UpperLim[2]= LowerLim[2] + ParValue.Z_AnalyzeRange;

    for(int i=0;i<3;i++)
    {
        if(CheckingPosn[i]>=LowerLim[i] && CheckingPosn[i]<=UpperLim[i])
            Counter++;
    }
    if (Counter==3)
        return 1;
    else
        return 0;
}
