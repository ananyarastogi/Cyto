 #include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
using namespace std;

#include "starter.h"
#include "output.h"
#include "main.h"
#include "random.h"
#include "cell.h"
#include "interact.h"
#include "outsidefunc.h"
#include "dataStorage.h"
#include "parameters.h"


//#include <random>
//#include <stdlib.h>
//#include <chrono>

// comment or uncomment
#define USE_GRAPHICS
#define SAVE_VIDEO false

#ifdef USE_GRAPHICS
#include "plot3d.h"
#endif
#include <QtWidgets/QApplication>

#define Readouts 8

//#define HeatMapHyp 15
//Number 1 change

using namespace std;

/*
 * Types of cells
 * */
enum {CTLTYPE=1, SOMATYPE=2};


SumOfSquare::SumOfSquare()
{
    FracDeadCells=0.0;
    Fig4I=0.0;
    DeadFig4E=0.0;
    AliveFig4E=0.0;
    Intact4F=0.0;
    Killed4F=0.0;
    Dead4D=0.0;
    Intact4D=0.0;
    AbsoluteSumOfSq=0.0;


    FracDeadCells_AfterAT=0.0;
    Fig4I_AfterAT=0.0;
    DeadFig4E_AfterAT=0.0;
    AliveFig4E_AfterAT=0.0;
    Intact4F_AfterAT=0.0;
    Killed4F_AfterAT=0.0;
    Dead4D_AfterAT=0.0;
    Intact4D_AfterAT=0.0;
    AbsoluteSumOfSq_AfterAT=0.0;


    MMH_FracDeadCells=0.0;
    MMH_Fig4I=0.0;
    MMH_DeadFig4E=0.0;
    MMH_AliveFig4E=0.0;
    MMH_Intact4F=0.0;
    MMH_Killed4F=0.0;
    MMH_Dead4D=0.0;
    MMH_Intact4D=0.0;
    MMH_AbsoluteSumOfSq=0.0;


    MMH_FracDeadCells_AfterAT=0.0;
    MMH_Fig4I_AfterAT=0.0;
    MMH_DeadFig4E_AfterAT=0.0;
    MMH_AliveFig4E_AfterAT=0.0;
    MMH_Intact4F_AfterAT=0.0;
    MMH_Killed4F_AfterAT=0.0;
    MMH_Dead4D_AfterAT=0.0;
    MMH_Intact4D_AfterAT=0.0;
    MMH_AbsoluteSumOfSq_AfterAT=0.0;

    n_AbsSS=0;
    MMH_n_AbsSS=0;
    n_AbsSSAfterAT=0;
    MMH_n_AbsSSAfterAT=0;
}



int main(int argc, char** argv)
{
    //cout<<"If you want or don't want videos, comment/uncomment the USE_GRAPHICS."<<endl<<endl;


#ifdef _WIN32
system("copy ..\\Cytokill\\DataInteractionDuration.txt DataInteractionDuration.txt");
system("copy ..\\Cytokill\\DataTurningAngle.txt DataTurningAngle.txt");
system("copy ..\\Cytokill\\SpeedHist.txt SpeedHist.txt");
#endif
#ifdef __linux__
system("cp ../Cytokill/DataInteractionDuration.txt DataInteractionDuration.txt");
system("cp ../Cytokill/DataTurningAngle.txt DataTurningAngle.txt");
system("cp ../Cytokill/SpeedHist.txt SpeedHist.txt");
#endif


    random::initialize(); // put your seed here
    sim_parameters* ParValue = new sim_parameters();

    //To read parameter values from txt file
    //ASK PHILIPPE: This needs to be updated as it does not have parameter ranges.
    // Philippe: I guess you don't need to change parameter ranges ? so just one set of ranges forever ?

    // If the parameter file is given as argument, just run it
    if(argc == 2)
    {
        string fileName = string(argv[1]);
        ParValue->ReadParameters(fileName);

        analysis(*ParValue);
        return 0;  

    }


    #ifdef USE_GRAPHICS
    //char* v = (char*) "abc"; char** argv = &v; // just to create a 'argv' because opengl needs one.
    if(ParValue->video) initOpenGL(argc, argv);
    #endif


    // if no argument is given, then launch graphical interface
    QApplication a(argc, argv);


    Starter* start = new Starter(ParValue);


    start->show();
    a.exec(); // looping
    cerr << "Finished" << endl;
    return 0;
}

bool stopSimu = false;
// function called from the interface
void setSimulationStatus(bool stop)
{
    stopSimu = stop;
}


/*
 * Function is called by main.
 * It calls all other functions that are involved in running of simulations
 * */
void RunSimulation(vector<float> &AvgAliveSom, vector<float> &AvgDeadSom, sim_parameters & ParValue,
                   vector<int> &DurationInt, int x, int y, float &AvgDeadPerSimSet, string appAddress,
                   vector<float> &AvgDeadSomAfterAT, vector<float>&AvgAliveSomAfterAT,
                   DataOutput &Observation,vector<SumOfSquare> &Results_SumSq, int &MMH_n_AbsSumSq,
                   int &MMH_n_AbsSumSq_AfterAT, int &n_AbsSumSq, int &n_AbsSumSq_AfterAT,
                   vector<vector<float>> &AllAliveData, vector<vector<float>> &AllDeadData,
                   vector<vector<float>> &AllAliveData_AfterAT, vector<vector<float>> &AllDeadData_AfterAT)
{
    ///Initialize
    vector <SomaticCellType> somatic;
    vector <CTLCellType> CTL;
    vector <Interact> InteractList;

    /*
     * TrumpHistory--> Stores CTL history in defined space and time window
     * CompleteHistory--> Stores complete CTL history over all of the space and during the whole simulation
     * ShortTimeHistory--> Stores history in complete physical space but only the designated time window
     * ShortSpaceHistory--> Stores history during all time but only the designated space window
     *
     * */
    vector <CTLHistory> TrumpHistory;
    vector <CTLHistory> CompleteHistory;
    vector <CTLHistory> ShortTimeHistory;
    vector <CTLHistory> ShortSpaceHistory;

    int DeadCells=0, Dead_TimeBin=0, Dead_TimeBinAfterAT=0,  Dead_TimeBinSS=0; /// To keep a count of the cells have died
    float t=0;
    bool Complete=0;

    /*
     * AddCells asignes the initial positions of the CTLs and infected cells
     * */
    AddCells(ParValue.num_CTL, CTLTYPE, somatic, CTL, ParValue, CompleteHistory, ShortSpaceHistory);
    AddCells(ParValue.num_somatic, SOMATYPE, somatic, CTL, ParValue, CompleteHistory, ShortSpaceHistory);


    /*
     * At each time step the following things happen:
     * (1) CTLs move
     * (2) After the CTLs finish moving, their histories are stored
     * (3) The CTLs and infected cells are checked for interactions
     * (4) The infected cells are checked for death. CheckforDeath calls another function based on the parameter Hypothesis
     * (5) Cells that are already interacting are checked to see if interaction goes on or it ends
     *
     * */
    for(t=ParValue.dT;t<=(ParValue.SimulationTotalTime+1.0);t+=ParValue.dT)
    {


        #ifdef USE_GRAPHICS
        bool saveThistime = ParValue.video && (fmod(t, ParValue.intervalOutput) < ParValue.dT);
        if(ParValue.video) display(somatic , CTL, t, saveThistime);
        if(saveThistime && ParValue.stacks) makeStacks(ParValue.nStacks, ParValue.Z_AnalyzeRange);
        //glutMainLoop();

        if(stopSimu) return;
        #endif


        CTLmove(CTL, somatic, ParValue, t);

        // Analysis function, records things that happened to each T cells in a list of Histories
        StoreCTLMoveHistory(CTL, ParValue, t, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, InteractList, somatic);

        // Analysis function, records things happening to somatic cells (records inside the somatic cell fields)
        StoreSomaticHistory(somatic, ParValue, t);

        // Initiation of interaction  when CTL and somatic get close (max 1 interaction per CTL, and more than one for somatic can be allowed or not)
        CheckInteraction (somatic, CTL, ParValue, InteractList, t, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory, Observation);


        // Dying decision, will depend on the chosen hypothesis
        CheckforDeath(somatic, CTL, ParValue, InteractList, DeadCells, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, t, DurationInt, TrumpHistory, ShortTimeHistory,
                      ShortSpaceHistory, CompleteHistory);

        // Decision to break an interaction. Note: better to die before removing an interaction, because some death decisions are taken at the end of an interaction. It should still exist
        AlreadyInteracting(somatic, CTL, ParValue, InteractList, DurationInt, t, TrumpHistory, CompleteHistory, ShortTimeHistory, ShortSpaceHistory);




        if((t!=0) &&  (fmod(t, ParValue.TimeBin)<ParValue.dT))//int(int(t)% int(ParValue.TimeBin))==0 && (abs((Factor*ParValue.TimeBin)-t)< ParValue.dT))
        {
            DeadCellsPerTimeBin(t, Dead_TimeBin, Dead_TimeBinAfterAT, Dead_TimeBinSS, ParValue, Observation, x, y, appAddress, CTL);
            Dead_TimeBin=0;
            Dead_TimeBinAfterAT=0;
            Dead_TimeBinSS=0;
        }


        /*
        int TM=int (int(t)/int(60));
        if(float(float(t)-float(TM*60))<ParValue.dT && float(float(t)-float(TM*60))>0)
        {
            cout<<t<<endl;
            int InfContacts[25];
            for(int qq=0;qq<25;qq++)
            {
                InfContacts[qq]=0;
            }
            for(int pp=0;pp<somatic.size();pp++)
            {
                int TempPosnCon=somatic[pp].PriorCellContacts;
                InfContacts[TempPosnCon]++;
            }

            for(int qq=0;qq<25;qq++)
            {
                cout<<qq<<"  "<<InfContacts[qq]<<endl;
            }
            cout<<endl<<endl;
        }
        */

    }



    float AvgZombieCont=0;
    int DeadSomNoZom=0, DeadSomWithZom=0;

    for(int zc=0;zc<somatic.size();zc++)
    {
        if(somatic[zc].Dead==1 && somatic[zc].ZombieContacts==0)
        {
            DeadSomNoZom++;
        }
        else if(somatic[zc].Dead==1 && somatic[zc].ZombieContacts!=0)
        {
            DeadSomWithZom++;
            AvgZombieCont+=somatic[zc].ZombieContacts;
        }

    }


    /*
     * These vectors are for each simulation. At the end of each simulation, the values are deleted
     * Vectors to store number of dead and alive infected cells at exactly x number of contacts with CTLs
     * */
    vector<float> NumCTLContacts_DeadSom(MaxContacts, 0);
    vector<float> NumCTLContacts_AliveSom(MaxContacts, 0);

    /*
     * These vectors are for each simulation. At the end of each simulation, the values are deleted
     * Vectors to store number of dead and alive infected cells at exactly x number of contacts with CTLs
     * But for these the observations are taken after the AnalyzeAfter time defined in sim_parameters
     * */
    vector<float>CTLContacts_DeadSom_AfterAT(MaxContacts, 0); //Cells analyzed after 60 mins ... AT= AnalyzeTime
    vector<float>CTLContacts_AliveSom_AfterAT(MaxContacts, 0);




    /*
     * Storing values in all 4 vectors defined above
     * */
    for (int i=0;i< (int) somatic.size();i++)
    {
        int k=somatic[i].PriorCellContacts;
        if(somatic[i].Dead==1)
            NumCTLContacts_DeadSom[k]++;
        else if(somatic[i].Dead==0)
            NumCTLContacts_AliveSom[k]++;
    }
    for (int i=0;i< (int) somatic.size();i++)
    {
        int k=somatic[i].ContactsAfterAnalyzeTime;
        if(somatic[i].Dead==1)
        {
            if(k!=0)
                CTLContacts_DeadSom_AfterAT[k]++;
        }
        else if(somatic[i].Dead==0)
            CTLContacts_AliveSom_AfterAT[k]++;
    }


    /*
    for (int i=0;i< (int) somatic.size();i++)
    {
        int k=somatic[i].ContactsAfter120;
        if(somatic[i].Dead==1)
        {
            if(k!=0)
                CTLContacts_DeadSom_120[k]++;
        }
        else if(somatic[i].Dead==0)
            CTLContacts_AliveSom_120[k]++;
    }


    for (int i=0;i< (int) somatic.size();i++)
    {
        int k=somatic[i].ContactsAfter180;
        if(somatic[i].Dead==1)
        {
            if(k!=0)
                CTLContacts_DeadSom_180[k]++;
        }
        else if(somatic[i].Dead==0)
            CTLContacts_AliveSom_180[k]++;
    }
    */

    /*
     * Adding values stored in the NumCTLContacts_AliveSom, NumCTLContacts_DeadSom, CTLContacts_AliveSom_AfterAT, CTLContacts_DeadSom_AfterAT
     * to vectors AvgAliveSom, AvgDeadCells, AvgAliveSomAfterAT, AvgDeadSomAfterAT
     * Later used to calculate average over all simulations in each parameter set
     * */
    for(int i=0;i<MaxContacts;i++)
    {
        AvgAliveSom[i]+=NumCTLContacts_AliveSom[i];

        AvgDeadSom[i]+=NumCTLContacts_DeadSom[i];

        AvgAliveSomAfterAT[i]+=CTLContacts_AliveSom_AfterAT[i];

        AvgDeadSomAfterAT[i]+=CTLContacts_DeadSom_AfterAT[i];


    }
    AllAliveData.push_back(NumCTLContacts_AliveSom);
    AllDeadData.push_back(NumCTLContacts_DeadSom);

    AllAliveData_AfterAT.push_back(CTLContacts_AliveSom_AfterAT);
    AllDeadData_AfterAT.push_back(CTLContacts_DeadSom_AfterAT);





    AvgDeadPerSimSet=AvgDeadPerSimSet+((float (DeadCells))/(float (ParValue.SimPerParSet)));

    if((y==(ParValue.SimPerParSet-1)) && (x==(ParValue.SimulationSets-1)))
    {
        Complete=1;
    }
    else
        Complete=0;

    /*
     * When the last simulation for a particular parameter set has finished running
     * */
    if(y==(ParValue.SimPerParSet-1))
    {
        NumDeadCells_vs_Parameter(AvgDeadPerSimSet, ParValue, Complete, appAddress, Observation);
        //TimeBetweenContacts_Infected(x, Observation, ParValue, appAddress);
    }

    CalcIndividualPKCR(x, y, Observation, ParValue, appAddress, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory);

    NumCells_vs_Contacts(ParValue, somatic, appAddress, x, y, Observation, Results_SumSq, MMH_n_AbsSumSq, MMH_n_AbsSumSq_AfterAT,
                         n_AbsSumSq, n_AbsSumSq_AfterAT);


    TotalDurationInfected(ParValue, somatic, appAddress, x, y, Observation, Results_SumSq, MMH_n_AbsSumSq, MMH_n_AbsSumSq_AfterAT,
                          n_AbsSumSq, n_AbsSumSq_AfterAT);

    KilledPerCTL(x, y, Observation, ParValue, appAddress, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory);

    //TotalDurationCTL(x, y, ParValue, appAddress, TrumpHistory, ShortTimeHistory, ShortSpaceHistory, CompleteHistory, Observation);
    IndividualDurationCTL(x, y, ParValue, appAddress, somatic, Observation, Results_SumSq, n_AbsSumSq, n_AbsSumSq_AfterAT, MMH_n_AbsSumSq, MMH_n_AbsSumSq_AfterAT);


    /*
    int Hy=ParValue.Hypothesis;
    if(Hy==2 || Hy==3 || Hy==4 || Hy==5 || Hy==7)
        PrintDamageInfected(x, y, ParValue, appAddress, somatic);
        */

    /*
     * Clearing all vectors at the end of each simulation
     * */
    NumCTLContacts_AliveSom.clear();
    NumCTLContacts_DeadSom.clear();
    somatic.clear();
    CTL.clear();
    InteractList.clear();

    ///Clear History
    TrumpHistory.clear();
    CompleteHistory.clear();
    ShortSpaceHistory.clear();
    ShortTimeHistory.clear();

}


// Performs an analysis with varying a parameter
void analysis(sim_parameters &ParValue)
{

    DataOutput Observation;
    /*
     * BullshitDetector returns 0 whenever some parameter combinations don't add up
     * */

    if((ParValue.SimulationSets > 1) && (ParValue.BullshitDetector()))
    {
        cout<<"Wrong parameter combination!"<<endl<<"Shit has happened."<<endl<<"Save yourself and RUN!!! "<<endl;
        exit(0);
    }

    /*
     * Generation of file name. Includes date and time.
     * */
    std::time_t ct = std::time(NULL);
    char foldername[35];
    std::strftime(foldername,
                  sizeof(foldername),
                  "Results%d-%m-%Yat%H-%M-%S",
                  std::localtime(&ct));
    string fn= foldername;



    if(ParValue.InteractDuringApoptosis==true)
    {
        if(ParValue.Hypothesis==5 && ParValue.LongContactsControlledDamage==true)
            fn="H"+std::to_string(int(ParValue.Hypothesis))+"_Z_LD_"+QString::number(ParValue.ThresholdContactDurDamage).toStdString()+"_"+fn;
        else if (ParValue.Hypothesis==7)
            fn="H"+std::to_string(int(ParValue.Hypothesis))+"_Z_"+QString::number(ParValue.TimetoDieMean).toStdString()+"_"+
                    QString::number(ParValue.AlphaDamage_Hyp7).toStdString()+"_"+fn;
        else
            fn="H"+std::to_string(int(ParValue.Hypothesis))+"_Z_"+fn;

    }
    else if(ParValue.InteractDuringApoptosis==false)
    {
        if(ParValue.Hypothesis==5 && ParValue.LongContactsControlledDamage==true)
            fn="H"+std::to_string(int(ParValue.Hypothesis))+"_NZ_LD_"+QString::number(ParValue.ThresholdContactDurDamage).toStdString()+"_"+fn;
        else if (ParValue.Hypothesis==7)
            fn="H"+std::to_string(int(ParValue.Hypothesis))+"_NZ_"+QString::number(ParValue.TimetoDieMean).toStdString()+"_"+
                    QString::number(ParValue.AlphaDamage_Hyp7).toStdString()+"_"+fn;
        else
            fn="H"+std::to_string(int(ParValue.Hypothesis))+"_NZ_"+fn;
    }


    createFolder(fn);

#ifdef USE_GRAPHICS
 setPrefix3Dfiles(fn);
 setObsWindow(ParValue.X_LowerLim, ParValue.X_LowerLim + ParValue.X_AnalyzeRange, ParValue.Y_LowerLim, ParValue.Y_LowerLim + ParValue.Y_AnalyzeRange, ParValue.Z_LowerLim, ParValue.Z_LowerLim + ParValue.Z_AnalyzeRange);
 set3DBox(0, ParValue.box_X, 0, ParValue.box_Y, 0, ParValue.box_Z);
 setFlatColors(ParValue.flatColors);
 setSizeOneStack(((double) ParValue.Z_AnalyzeRange) / (((double) ParValue.nStacks)));
#endif

    setVideoSaveFolder(fn);
    /*
     * Asks users for parameter values. Has to be updated to be used by S. Halle
     * */
    //ALERT: Check if these values are stored amd reproduced in the report!!!!!!!!!!


    DefaultParameterRange ParRange;
    if(ParValue.VariableParMin==0.0 && ParValue.VariableParMax==0.0 && ParValue.SimulationSets!=1)
    {
        cout<<"Variable parameter has no range!!"<<endl;
        cout<<"The program will crash now, you dumbass."<<endl;
        exit(10);
    }
    ParRange.MinValues[ParValue.VariableParameter]=ParValue.VariableParMin;
    ParRange.MaxValues[ParValue.VariableParameter]=ParValue.VariableParMax;


    WriteReport(fn, ParValue, ParRange);
    GenerateRScript(fn, ParValue);



    vector <SumOfSquare> Results_SumSq;

    std::string appAddress;         /// Address of application (ending with "/")
    appAddress = currentDir();
    appAddress.append("/"+fn+"/");
    cout<<appAddress<<endl;

    float ParHypMin, ParHypMax;


    float ParHypStep=0.0;

    HypDependentParRange(ParValue, ParHypMin, ParHypMax);
    ParHypStep=(ParHypMax-ParHypMin)/float(HeatMapHyp-1);

    ofstream myfileVarPar(string(appAddress+"VarPar")+string(".txt"));
    for(float ParHyp=ParHypMin;ParHyp<=(ParHypMax+(1e-4));ParHyp+=ParHypStep)
    {
        myfileVarPar<<ParHyp<<endl;
    }
    myfileVarPar.close();



    vector <vector<int>> UltimateDeadCellsList;
    vector <int> UltimateNumSomatic;
    vector <int> UltimateNumCTL;
    if(ParValue.velocity_mag_Distribution != FromData) ParValue.DistributionSpeed.set(ParValue.velocity_mag_Distribution, ParValue.velocity_mag_Mean, ParValue.velocity_mag_SD);
    else ParValue.DistributionSpeed.set(string("SpeedHist.txt"));

    if(ParValue.InteractTimeDistribution != FromData) ParValue.DistributionInteractionDuration.set(ParValue.InteractTimeDistribution, ParValue.InteractTimeMean, ParValue.InteractTimeSD);
    else ParValue.DistributionInteractionDuration.set(string("DataInteractionDuration.txt"));

    if(ParValue.PersistentTimeDistribution != FromData) ParValue.DistributionPersistentTime.set(ParValue.PersistentTimeDistribution, ParValue.PersistentTimeMean, ParValue.PersistentTimeSD);
    else ParValue.DistributionPersistentTime.set(string(""));

    if(ParValue.TurningAngleDistribution != FromData) ParValue.DistributionTurningAngles.set(ParValue.TurningAngleDistribution, ParValue.TurningAngleMean, ParValue.TurningAngleSD);
    else ParValue.DistributionTurningAngles.set(string("DataTurningAngle.txt"));

    for(float ParHyp=ParHypMin;ParHyp<=(ParHypMax+(1e-4));ParHyp+=ParHypStep)
    {

        if(HeatMapHyp!=1)
        {
            HypDependentChangePar(ParValue, ParHyp);
        }
        cout<<ParValue.ProbDie_Interact_Hyp1<<"ff"<<endl;
        vector<float> VarParList;

        for(int x=0;x<ParValue.SimulationSets;x++)
        {

            int MMH_n_AbsSumSq=0, MMH_n_AbsSumSq_AfterAT=0;

            int n_AbsSumSq=0, n_AbsSumSq_AfterAT=0;
            cout<<"Value of x: "<<x<<" / "<<ParHyp<<" "<<ParValue.KillingFactorMean_Hyp8<<endl;
            float TempRatio=(float(x)/float(ParValue.SimulationSets-1));

            ParametersChanging(ParValue, TempRatio, ParRange.MinValues[ParValue.VariableParameter], ParRange.MaxValues[ParValue.VariableParameter]);



            VarParList.push_back(ParValue.ValChangingVariable);
            cout<<"Pay Attention here. If you want the duration of simulation to be longer then change the size of these vectors otherwise there will be segmentation error"<<endl<<endl;
            vector<float> AvgDeadSom(MaxContacts, 0), AvgAliveSom(MaxContacts, 0);
            vector<float>AvgDeadSomAfterAT(MaxContacts, 0), AvgAliveSomAfterAT(MaxContacts, 0);

            float AvgDeadPerSimSet=0.0;
            SumOfSquare temp;
            Results_SumSq.push_back(temp);

            if(ParValue.TimetoDieDistribution != FromData) ParValue.DistributionTimeToDie.set(ParValue.TimetoDieDistribution, ParValue.TimetoDieMean, ParValue.TimetoDieSD);
            else ParValue.DistributionTimeToDie.set(string(""));











            //ParValue.num_CTL=20*(x+1);
            UltimateNumSomatic.push_back(ParValue.num_somatic);
            UltimateNumCTL.push_back(ParValue.num_CTL);












            vector<vector<float>> AllAliveData, AllDeadData, AllAliveData_AfterAT, AllDeadData_AfterAT;

            vector<int> SingleSimDeadCell;

            for(int y=0;y<ParValue.SimPerParSet;y++)
            {



                vector<int> DurationInt(ParValue.SimulationTotalTime, 0);

                RunSimulation(AvgAliveSom, AvgDeadSom, ParValue, DurationInt, x, y, AvgDeadPerSimSet, appAddress,
                              AvgDeadSomAfterAT, AvgAliveSomAfterAT, Observation, Results_SumSq, MMH_n_AbsSumSq, MMH_n_AbsSumSq_AfterAT,
                              n_AbsSumSq, n_AbsSumSq_AfterAT, AllAliveData, AllDeadData, AllAliveData_AfterAT, AllDeadData_AfterAT);




                if(ParValue.SimulationSets==1)
                {
                    exit(10);
                }

                int Temp_Dead=0;
                for(uint td=0;td<AllDeadData[y].size();td++)
                {
                    Temp_Dead=Temp_Dead + AllDeadData[y][td];
                }
                SingleSimDeadCell.push_back(Temp_Dead);

                /*
                 * Analysis of the data that has been generated
                 * */
                if(y==(ParValue.SimPerParSet-1))
                {
                    int AfterAT;
                    AfterAT=0;
                    AnalysisMethod(ParValue, x, appAddress, DurationInt, AfterAT, Results_SumSq, MMH_n_AbsSumSq, n_AbsSumSq,
                                   Observation);
                    AfterAT=1;
                    AnalysisMethod(ParValue, x, appAddress, DurationInt, AfterAT, Results_SumSq, MMH_n_AbsSumSq_AfterAT,
                                     n_AbsSumSq_AfterAT, Observation);


                    AfterAT=-1;
                    AnalysisMethod(ParValue, x, appAddress, DurationInt, AfterAT, Results_SumSq, MMH_n_AbsSumSq_AfterAT,
                                     n_AbsSumSq_AfterAT, Observation);

                    AfterAT=-2;
                    AnalysisMethod(ParValue, x, appAddress, DurationInt, AfterAT, Results_SumSq, MMH_n_AbsSumSq_AfterAT,
                                     n_AbsSumSq_AfterAT, Observation);

                    AfterAT=-3;
                    AnalysisMethod(ParValue, x, appAddress, DurationInt, AfterAT, Results_SumSq, MMH_n_AbsSumSq_AfterAT,
                                     n_AbsSumSq_AfterAT, Observation);


                    AvgDeadSom.clear();
                    AvgAliveSom.clear();
                    AvgAliveSomAfterAT.clear();
                    AvgDeadSomAfterAT.clear();
                }
            }
            UltimateDeadCellsList.push_back(SingleSimDeadCell);
            SingleSimDeadCell.clear();


            Results_SumSq[x].n_AbsSS= n_AbsSumSq;
            Results_SumSq[x].MMH_n_AbsSS= MMH_n_AbsSumSq;
            Results_SumSq[x].n_AbsSSAfterAT= n_AbsSumSq_AfterAT;
            Results_SumSq[x].MMH_n_AbsSSAfterAT= MMH_n_AbsSumSq_AfterAT;

        }

        ofstream myfileTotDead(string(appAddress+"UltimateDeadList") + string(".txt"));



        for(uint td=0;td<UltimateDeadCellsList.size();td++)
        {

            if(td==0)
            {
                myfileTotDead<<"\"Number of somatic cells\""<<"\t"<<"\"Number of CTLs\"";
                for(uint sd=0;sd<UltimateDeadCellsList[0].size();sd++)
                {
                    myfileTotDead<<"\t"<<"\"Sim"<<sd<<"\"";
                }
                myfileTotDead<<endl;
            }


            for(uint sd=0;sd<UltimateDeadCellsList[0].size();sd++)
            {
                if(sd==0)
                {
                    myfileTotDead<<UltimateNumSomatic[td]<<"\t"<<UltimateNumCTL[td];
                }
                myfileTotDead<<"\t"<<UltimateDeadCellsList[td][sd];
            }
            myfileTotDead<<endl;
        }
        myfileTotDead.close();

        SumOfSquare AvgSumSq;


        ofstream myfileSumOfSq(string(appAddress+"SumSq") + QString::number(ParHyp).toStdString()+ string(".txt"));

        myfileSumOfSq<<"#The parameter(s) being varied is/are: ";//----------------------------------change

        myfileSumOfSq<<NameParameter(ParValue.VariableParameter)<<endl<<endl;

        myfileSumOfSq<<"#Value of above parameter is "<<ParValue.ValChangingVariable<<endl;

        myfileSumOfSq<<"\"Variable Parameter\""<<'\t'<<"\"Fig7L\""<<'\t';
        myfileSumOfSq<<"\"Fig4I\""<<'\t';
        myfileSumOfSq<<"\"Fig4E Dead\""<<'\t';
        myfileSumOfSq<<"\"Fig4E Intact\""<<'\t';
        myfileSumOfSq<<"\"Fig4F Killed\""<<'\t';
        myfileSumOfSq<<"\"Fig4F Intact\""<<'\t';
        myfileSumOfSq<<"\"Fig 4D Dead\""<<'\t';
        myfileSumOfSq<<"\"Fig 4D Intact\""<<'\t';
        myfileSumOfSq<<"\"Absolute\""<<endl;

        /*
        if(VarParList.size()!=ParValue.SimulationSets)
        {
            cout<<"Crap. Crap. Bullshit";
            exit(100);
        }
        */
        for(uint i=0;i<VarParList.size();i++)
        {
            myfileSumOfSq<<VarParList[i]<<'\t'<<Results_SumSq[i].FracDeadCells;
            AvgSumSq.FracDeadCells+=Results_SumSq[i].FracDeadCells;

            myfileSumOfSq<<'\t'<<Results_SumSq[i].Fig4I<<'\t';
            AvgSumSq.Fig4I+=Results_SumSq[i].Fig4I;

            myfileSumOfSq<<Results_SumSq[i].DeadFig4E<<'\t';
            AvgSumSq.DeadFig4E+=Results_SumSq[i].DeadFig4E;

            myfileSumOfSq<<Results_SumSq[i].AliveFig4E<<'\t';
            AvgSumSq.AliveFig4E+=Results_SumSq[i].AliveFig4E;

            myfileSumOfSq<<Results_SumSq[i].Killed4F<<'\t';
            AvgSumSq.Killed4F+=Results_SumSq[i].Killed4F;

            myfileSumOfSq<<Results_SumSq[i].Intact4F<<'\t';
            AvgSumSq.Intact4F+=Results_SumSq[i].Intact4F;

            myfileSumOfSq<<Results_SumSq[i].Dead4D<<'\t';
            AvgSumSq.Dead4D+=Results_SumSq[i].Dead4D;

            myfileSumOfSq<<Results_SumSq[i].Intact4D<<'\t';
            AvgSumSq.Intact4D+=Results_SumSq[i].Intact4D;

            SumOfSquare S=Results_SumSq[i];
            float tempSumSq=S.FracDeadCells + S.Fig4I + S.DeadFig4E + S.AliveFig4E + S.Killed4F + S.Intact4F + S.Dead4D + S.Intact4D;


            //myfileSumOfSq<<float(Results_SumSq[i].AbsoluteSumOfSq)/float(Results_SumSq[i].n_AbsSS)<<endl;


            myfileSumOfSq<<(tempSumSq/float(Readouts))<<endl;
            AvgSumSq.AbsoluteSumOfSq+=(tempSumSq/float(Readouts));
        }
        float sizeVarPar=float(VarParList.size());
        myfileSumOfSq<<"\"Average\""<<'\t'<<AvgSumSq.FracDeadCells/sizeVarPar<<'\t'<<AvgSumSq.Fig4I/sizeVarPar<<'\t'<<AvgSumSq.DeadFig4E/sizeVarPar<<'\t';
        myfileSumOfSq<<AvgSumSq.AliveFig4E/sizeVarPar<<'\t'<<AvgSumSq.Killed4F/sizeVarPar<<'\t'<<AvgSumSq.Intact4F/sizeVarPar<<'\t'<<AvgSumSq.Dead4D/sizeVarPar<<'\t';
        myfileSumOfSq<<AvgSumSq.Intact4D/sizeVarPar<<'\t'<<AvgSumSq.AbsoluteSumOfSq/sizeVarPar<<endl;


        myfileSumOfSq.close();

        ofstream myfileSumOfSq_AfterAT(string(appAddress+"SumSq_AfterAT")+QString::number(ParHyp).toStdString()+string(".txt"));

        myfileSumOfSq_AfterAT<<"#The parameter(s) being varied is/are: ";//----------------------------------change

        myfileSumOfSq_AfterAT<<NameParameter(ParValue.VariableParameter)<<endl<<endl;

        myfileSumOfSq_AfterAT<<"#Value of above parameter is "<<ParValue.ValChangingVariable<<endl;

        myfileSumOfSq_AfterAT<<"\"Variable Parameter\""<<'\t'<<"\"Fig7L AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig4I AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig4E Dead AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig4E Intact AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig4F Killed AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig4F Intact AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig 4D Dead AfterAT\""<<'\t';
        myfileSumOfSq_AfterAT<<"\"Fig 4D Intact AfterAT\""<<'\t';

        myfileSumOfSq_AfterAT<<"\"Absolute AfterAT\""<<endl;
        for(uint i=0;i<VarParList.size();i++)
        {
            myfileSumOfSq_AfterAT<<VarParList[i]<<'\t'<<Results_SumSq[i].FracDeadCells_AfterAT;
            AvgSumSq.FracDeadCells_AfterAT+=Results_SumSq[i].FracDeadCells_AfterAT;

            myfileSumOfSq_AfterAT<<'\t'<<Results_SumSq[i].Fig4I_AfterAT<<'\t';
            AvgSumSq.Fig4I_AfterAT+=Results_SumSq[i].Fig4I_AfterAT;

            myfileSumOfSq_AfterAT<<Results_SumSq[i].DeadFig4E_AfterAT<<'\t';
            AvgSumSq.DeadFig4E_AfterAT+=Results_SumSq[i].DeadFig4E_AfterAT;

            myfileSumOfSq_AfterAT<<Results_SumSq[i].AliveFig4E_AfterAT<<'\t';
            AvgSumSq.AliveFig4E_AfterAT+=Results_SumSq[i].AliveFig4E_AfterAT;

            myfileSumOfSq_AfterAT<<Results_SumSq[i].Killed4F_AfterAT<<'\t';
            AvgSumSq.Killed4F_AfterAT+=Results_SumSq[i].Killed4F_AfterAT;

            myfileSumOfSq_AfterAT<<Results_SumSq[i].Intact4F_AfterAT<<'\t';
            AvgSumSq.Intact4F_AfterAT+=Results_SumSq[i].Intact4F_AfterAT;

            myfileSumOfSq_AfterAT<<Results_SumSq[i].Dead4D_AfterAT<<'\t';
            AvgSumSq.Dead4D_AfterAT+=Results_SumSq[i].Dead4D_AfterAT;

            myfileSumOfSq_AfterAT<<Results_SumSq[i].Intact4D_AfterAT<<'\t';
            AvgSumSq.Intact4D_AfterAT+=Results_SumSq[i].Intact4D_AfterAT;

            SumOfSquare S=Results_SumSq[i];
            float tempSumSq=S.FracDeadCells_AfterAT + S.Fig4I_AfterAT + S.DeadFig4E_AfterAT + S.AliveFig4E_AfterAT +
                    S.Killed4F_AfterAT + S.Intact4F_AfterAT + S.Dead4D_AfterAT + S.Intact4D_AfterAT;


            myfileSumOfSq_AfterAT<<(tempSumSq/float(Readouts))<<endl;

            AvgSumSq.AbsoluteSumOfSq_AfterAT+=(tempSumSq/float(Readouts));

            //myfileSumOfSq_AfterAT<<float(Results_SumSq[i].AbsoluteSumOfSq_AfterAT)/float(Results_SumSq[i].n_AbsSSAfterAT)<<endl;

        }

        myfileSumOfSq_AfterAT<<"\"Average\""<<'\t'<<AvgSumSq.FracDeadCells_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.Fig4I_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.DeadFig4E_AfterAT/sizeVarPar<<'\t';
        myfileSumOfSq_AfterAT<<AvgSumSq.AliveFig4E_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.Killed4F_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.Intact4F_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.Dead4D_AfterAT/sizeVarPar<<'\t';
        myfileSumOfSq_AfterAT<<AvgSumSq.Intact4D_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.AbsoluteSumOfSq_AfterAT/sizeVarPar<<endl;
        myfileSumOfSq_AfterAT.close();



        ofstream myfileSumOfSqMMH(string(appAddress+"MMHSumSq")+QString::number(ParHyp).toStdString()+string(".txt"));

        myfileSumOfSqMMH<<"#The parameter(s) being varied is/are: ";//----------------------------------change

        myfileSumOfSqMMH<<NameParameter(ParValue.VariableParameter)<<endl<<endl;

        myfileSumOfSqMMH<<"#Value of above parameter is "<<ParValue.ValChangingVariable<<endl;

        myfileSumOfSqMMH<<"\"Variable Parameter\""<<'\t'<<"\"Fig7L\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4I\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4E Dead\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4E Intact\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4F Killed\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4F Intact\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4D Killed\""<<'\t';
        myfileSumOfSqMMH<<"\"Fig4D Intact\""<<'\t';
        myfileSumOfSqMMH<<"\"Absolute\""<<endl;

        for(uint i=0;i<VarParList.size();i++)
        {
            myfileSumOfSqMMH<<VarParList[i]<<'\t'<<Results_SumSq[i].MMH_FracDeadCells;
            AvgSumSq.MMH_FracDeadCells+=Results_SumSq[i].MMH_FracDeadCells;

            myfileSumOfSqMMH<<'\t'<<Results_SumSq[i].MMH_Fig4I<<'\t';
            AvgSumSq.MMH_Fig4I+=Results_SumSq[i].MMH_Fig4I;

            myfileSumOfSqMMH<<Results_SumSq[i].MMH_DeadFig4E<<'\t';
            AvgSumSq.MMH_DeadFig4E+=Results_SumSq[i].MMH_DeadFig4E;

            myfileSumOfSqMMH<<Results_SumSq[i].MMH_AliveFig4E<<'\t';
            AvgSumSq.MMH_AliveFig4E+=Results_SumSq[i].MMH_AliveFig4E;

            myfileSumOfSqMMH<<Results_SumSq[i].MMH_Killed4F<<'\t';
            AvgSumSq.MMH_Killed4F+=Results_SumSq[i].MMH_Killed4F;

            myfileSumOfSqMMH<<Results_SumSq[i].MMH_Intact4F<<'\t';
            AvgSumSq.MMH_Intact4F+=Results_SumSq[i].MMH_Intact4F;

            myfileSumOfSqMMH<<Results_SumSq[i].MMH_Dead4D<<'\t';
            AvgSumSq.MMH_Dead4D+=Results_SumSq[i].MMH_Dead4D;

            myfileSumOfSqMMH<<Results_SumSq[i].MMH_Intact4D<<'\t';
            AvgSumSq.MMH_Intact4D+=Results_SumSq[i].MMH_Intact4D;

            SumOfSquare S=Results_SumSq[i];
            float tempSumSq=S.MMH_FracDeadCells + S.MMH_Fig4I + S.MMH_DeadFig4E + S.MMH_AliveFig4E +
                    S.MMH_Killed4F + S.MMH_Intact4F + S.MMH_Dead4D + S.MMH_Intact4D;


            myfileSumOfSqMMH<<(tempSumSq/float(Readouts))<<endl;

            //myfileSumOfSqMMH<<float(Results_SumSq[i].MMH_AbsoluteSumOfSq)/float(Results_SumSq[i].MMH_n_AbsSS)<<endl;
            AvgSumSq.MMH_AbsoluteSumOfSq+=(tempSumSq/float(Readouts));

        }

        myfileSumOfSqMMH<<"\"Average\""<<'\t'<<AvgSumSq.MMH_FracDeadCells/sizeVarPar<<'\t'<<AvgSumSq.MMH_Fig4I/sizeVarPar<<'\t'<<AvgSumSq.MMH_DeadFig4E/sizeVarPar<<'\t';
        myfileSumOfSqMMH<<AvgSumSq.MMH_AliveFig4E/sizeVarPar<<'\t'<<AvgSumSq.MMH_Killed4F/sizeVarPar<<'\t'<<AvgSumSq.MMH_Intact4F/sizeVarPar<<'\t'<<AvgSumSq.MMH_Dead4D/sizeVarPar<<'\t';
        myfileSumOfSqMMH<<AvgSumSq.MMH_Intact4D/sizeVarPar<<'\t'<<AvgSumSq.MMH_AbsoluteSumOfSq/sizeVarPar<<endl;


        myfileSumOfSqMMH.close();


        ofstream myfileSumOfSqMMH_AfterAT(string(appAddress+"MMHSumSq_AfterAT")+QString::number(ParHyp).toStdString()+string(".txt"));

        myfileSumOfSqMMH_AfterAT<<"#The parameter(s) being varied is/are: ";//----------------------------------change

        myfileSumOfSqMMH_AfterAT<<NameParameter(ParValue.VariableParameter)<<endl<<endl;

        myfileSumOfSqMMH_AfterAT<<"#Value of above parameter is "<<ParValue.ValChangingVariable<<endl;

        myfileSumOfSqMMH_AfterAT<<"\"Variable Parameter\""<<'\t'<<"\"Fig7L AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4I AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4E Dead AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4E Intact AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4F Killed AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4F Intact AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4D Killed AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Fig4D Intact AfterAT\""<<'\t';
        myfileSumOfSqMMH_AfterAT<<"\"Absolute AfterAT\""<<endl;
        for(uint i=0;i<VarParList.size();i++)
        {
            myfileSumOfSqMMH_AfterAT<<VarParList[i]<<'\t'<<Results_SumSq[i].MMH_FracDeadCells_AfterAT;
            AvgSumSq.MMH_FracDeadCells_AfterAT+=Results_SumSq[i].MMH_FracDeadCells_AfterAT;

            myfileSumOfSqMMH_AfterAT<<'\t'<<Results_SumSq[i].MMH_Fig4I_AfterAT<<'\t';
            AvgSumSq.MMH_Fig4I_AfterAT+=Results_SumSq[i].MMH_Fig4I_AfterAT;

            myfileSumOfSqMMH_AfterAT<<Results_SumSq[i].MMH_DeadFig4E_AfterAT<<'\t';
            AvgSumSq.MMH_DeadFig4E_AfterAT+=Results_SumSq[i].MMH_DeadFig4E_AfterAT;

            myfileSumOfSqMMH_AfterAT<<Results_SumSq[i].MMH_AliveFig4E_AfterAT<<'\t';
            AvgSumSq.MMH_AliveFig4E_AfterAT+=Results_SumSq[i].MMH_AliveFig4E_AfterAT;

            myfileSumOfSqMMH_AfterAT<<Results_SumSq[i].MMH_Killed4F_AfterAT<<'\t';
            AvgSumSq.MMH_Killed4F_AfterAT+=Results_SumSq[i].MMH_Killed4F_AfterAT;

            myfileSumOfSqMMH_AfterAT<<Results_SumSq[i].MMH_Intact4F_AfterAT<<'\t';
            AvgSumSq.MMH_Intact4F_AfterAT+=Results_SumSq[i].MMH_Intact4F_AfterAT;

            myfileSumOfSqMMH_AfterAT<<Results_SumSq[i].MMH_Dead4D_AfterAT<<'\t';
            AvgSumSq.MMH_Dead4D_AfterAT+=Results_SumSq[i].MMH_Dead4D_AfterAT;

            myfileSumOfSqMMH_AfterAT<<Results_SumSq[i].MMH_Intact4D_AfterAT<<'\t';
            AvgSumSq.MMH_Intact4D_AfterAT+=Results_SumSq[i].MMH_Intact4D_AfterAT;



            SumOfSquare S=Results_SumSq[i];
            float tempSumSq=S.MMH_FracDeadCells_AfterAT + S.MMH_Fig4I_AfterAT + S.MMH_DeadFig4E_AfterAT + S.MMH_AliveFig4E_AfterAT +
                    S.MMH_Killed4F_AfterAT + S.MMH_Intact4F_AfterAT + S.MMH_Dead4D_AfterAT + S.MMH_Intact4D_AfterAT;


            myfileSumOfSqMMH_AfterAT<<(tempSumSq/float(Readouts))<<endl;
            AvgSumSq.MMH_AbsoluteSumOfSq_AfterAT+=(tempSumSq/float(Readouts));


            //myfileSumOfSqMMH_AfterAT<<float(Results_SumSq[i].MMH_AbsoluteSumOfSq_AfterAT)/float(Results_SumSq[i].MMH_n_AbsSSAfterAT)<<endl;
        }

        myfileSumOfSqMMH_AfterAT<<"\"Average\""<<'\t'<<AvgSumSq.MMH_FracDeadCells_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.MMH_Fig4I_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.MMH_DeadFig4E_AfterAT/sizeVarPar<<'\t';
        myfileSumOfSqMMH_AfterAT<<AvgSumSq.MMH_AliveFig4E_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.MMH_Killed4F_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.MMH_Intact4F_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.MMH_Dead4D_AfterAT/sizeVarPar<<'\t';
        myfileSumOfSqMMH_AfterAT<<AvgSumSq.MMH_Intact4D_AfterAT/sizeVarPar<<'\t'<<AvgSumSq.MMH_AbsoluteSumOfSq_AfterAT/sizeVarPar<<endl;


        myfileSumOfSqMMH_AfterAT.close();


        Results_SumSq.clear();
    }

     ofstream myfileEnd(string(appAddress+"AAAA_OVER")+string(".txt"));
     myfileEnd<<ParValue.TimetoDieMean<<endl;
     myfileEnd<<ParValue.TimetoDieSD<<endl;
     myfileEnd<<"Time to die mean"<<endl;
     myfileEnd<<"Time to die SD"<<endl;
     myfileEnd<<"#Done!!";
     myfileEnd.close();

#ifdef USE_GRAPHICS


    // glutMainLoop(); // this crashes ??
#endif


}




/*
 * fracDead(k)=D(k)/[D(k)+A(k)+D(k+1)+A(k+1)+...+D(n)+A(n)]
 * This is the analysis that would fit according to Bernoulli trial
 * Generating text file to store results
 * */
/*
void AnalysisMethod(vector<float> &AvgAliveSom, vector<float> &AvgDeadSom, sim_parameters & ParValue, int x, string appAddress,
                      vector<int> DurationInt, int AfterAT, vector <SumOfSquare> &Results_SumSq, int &MMH_n_Abs, int &n_AbsSumSq,
                    vector<vector<float>> &AllAliveData, vector<vector<float>> &AllDeadData)
{
    int Contacts[7]={0, 1, 2, 3, 4, 5, 6};

    float ExpFracDead[7]={0.0106, 0.0903, 0.0833, 0.2143, 0.2424, 0.375, 0.5};


    float temp_SumSqFracDead=0.0;
    float MMH_Temp_AbsSumSq=0.0;
    int n=0;

    int MMH_n=0;
    string TextFileName;
    if(AfterAT==0)
        TextFileName="AvgCTLContacts";
    else if(AfterAT==1)
        TextFileName="AvgCTLContactsAfterAT";







    vector<float> FracDeadVector;
    for(unsigned int i=0;i<AvgAliveSom.size();i++)
    {
        float fracAlive=0.0, fracDead=0.0, fracAlive_PaperAnalysis=0.0, fracDead_PaperAnalysis=0.0;

        fracAlive_PaperAnalysis=AvgAliveSom[i]/(AvgAliveSom[i]+AvgDeadSom[i]);
        fracDead_PaperAnalysis=AvgDeadSom[i]/(AvgAliveSom[i]+AvgDeadSom[i]);
        float denominator=0.0;
        for(unsigned int d=i;d<AvgAliveSom.size();d++)
        {
            denominator=denominator+AvgAliveSom[d]+AvgDeadSom[d];
        }
        if(fabs(denominator)<1e-6)
        {
            fracAlive=0.0;
            fracDead=0.0;
        }
        else
        {
            fracDead=AvgDeadSom[i]/denominator;
            fracAlive=1-fracDead;
        }

        FracDeadVector.push_back(fracDead);
        if(i<7 && ExpFracDead[i]!=0)
        {
            temp_SumSqFracDead+= pow((fracDead-(ExpFracDead[i])) ,2);
            n_AbsSumSq++;
            n++;


            if(ExpFracDead[i]>1e-6)
            {
                MMH_Temp_AbsSumSq=MMH_Temp_AbsSumSq+((pow((fracDead-ExpFracDead[i]),2))/pow((ExpFracDead[i]) ,2));
                MMH_n_Abs++;
                MMH_n++;
            }
        }
    }

    if(AfterAT==0)
    {



        Results_SumSq[x].FracDeadCells=float(temp_SumSqFracDead/float(n));
        Results_SumSq[x].MMH_AbsoluteSumOfSq+=MMH_Temp_AbsSumSq;
        Results_SumSq[x].MMH_FracDeadCells=float(MMH_Temp_AbsSumSq/float(MMH_n));
        Results_SumSq[x].AbsoluteSumOfSq+=temp_SumSqFracDead;
    }
    else if(AfterAT==1)
    {


        Results_SumSq[x].FracDeadCells_AfterAT=float(temp_SumSqFracDead/float(n));
        Results_SumSq[x].MMH_AbsoluteSumOfSq_AfterAT+=MMH_Temp_AbsSumSq;
        Results_SumSq[x].MMH_FracDeadCells_AfterAT=float(MMH_Temp_AbsSumSq/float(MMH_n));
        Results_SumSq[x].AbsoluteSumOfSq_AfterAT+=temp_SumSqFracDead;

    }
}

*/








/*
 * fracDead(k)=D(k)/[D(k)+A(k)+D(k+1)+A(k+1)+...+D(n)+A(n)]
 * This is the analysis that would fit according to Bernoulli trial
 * Generating text file to store results
 * */

void AnalysisMethod(sim_parameters & ParValue, int x, string appAddress,
                      vector<int> DurationInt, int AfterAT, vector <SumOfSquare> &Results_SumSq, int &MMH_n_Abs, int &n_AbsSumSq,
                     DataOutput &Observation)
{
    // Philippe: why ?
    //cout<<endl<<"check duration interaction. it doesn't seem to be correct!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

    //int Contacts[7]={0, 1, 2, 3, 4, 5, 6};

    float ExpFracDead[7]={0.0106, 0.0903, 0.0833, 0.2143, 0.2424, 0.375, 0.5};
    if((AfterAT==0))
    {
        ofstream myfileDur(string(appAddress+"DurationInteraction")+QString::number(x).toStdString()+string(".txt"));
        myfileDur<<"#The parameter(s) being varied is/are: ";//----------------------------------change

        myfileDur<<NameParameter(ParValue.VariableParameter)<<endl<<endl;

        myfileDur<<"#Value of above parameter is "<<ParValue.ValChangingVariable<<endl;

        myfileDur<<"\"Duration in minutes\""<<'\t'<<"\"Number of interaction\""<<endl;

        for(int i=0;i<ParValue.SimulationTotalTime;i++)
        {
            if(DurationInt[i]!=0)
            {
                myfileDur<<float(i)+0.5<<'\t'<<DurationInt[i]<<endl;
            }
        }
        myfileDur.close();
    }


    float temp_SumSqFracDead=0.0;
    float MMH_Temp_AbsSumSq=0.0;
    int n=0;

    int MMH_n=0;
    vector<float> Alive;
    vector<float> Dead;
    string TextFileName;
    if(AfterAT==0)
    {
        TextFileName="AvgCTLContacts";
        for(unsigned int i=0;i<Observation.InfectedAlive.size();i++)
        {
            Alive.push_back(Observation.InfectedAlive[i]);
            Dead.push_back(Observation.InfectedDead[i]);
        }

    }
    else if(AfterAT==1)
    {
        TextFileName="AvgCTLContactsAfterAT";
        for(unsigned int i=0;i<Observation.InfectedAliveAfterAT.size();i++)
        {
            Alive.push_back(Observation.InfectedAliveAfterAT[i]);
            Dead.push_back(Observation.InfectedDeadAfterAT[i]);
        }
    }
    else if(AfterAT==-1)
    {
        TextFileName="AvgCTLContactsSet0_240_";
        for(unsigned int i=0;i<Observation.InfectedAlive0_240.size();i++)
        {
            Alive.push_back(Observation.InfectedAlive0_240[i]);
            Dead.push_back(Observation.InfectedDead0_240[i]);
        }

    }
    else if(AfterAT==-2)
    {
        TextFileName="AvgCTLContactsSet0_360_";
        for(unsigned int i=0;i<Observation.InfectedAlive0_360.size();i++)
        {
            Alive.push_back(Observation.InfectedAlive0_360[i]);
            Dead.push_back(Observation.InfectedDead0_360[i]);
        }
    }
    else if(AfterAT==-3)
    {
        TextFileName="AvgCTLContactsSet0_480_";
        for(unsigned int i=0;i<Observation.InfectedAlive0_480.size();i++)
        {
            Alive.push_back(Observation.InfectedAlive0_480[i]);
            Dead.push_back(Observation.InfectedDead0_480[i]);
        }
    }




    ofstream myfile(string(appAddress+TextFileName)+QString::number(x).toStdString()+string(".txt"));
    myfile<<"#The parameter(s) being varied is/are: ";//----------------------------------change

    myfile<<NameParameter(ParValue.VariableParameter)<<endl<<endl;

    myfile<<"#Value of above parameter is "<<ParValue.ValChangingVariable<<endl;



    vector<float> FracDeadVector;
    for(unsigned int i=0;i<Alive.size();i++)
    {


        if(i==0)
        {
            myfile<<"\"Number of CTL contacts\""<<'\t';
            myfile<<"\"Avg Frac Dead\""<<'\t'<<"\"Avg Frac Alive\""<<'\t';
            myfile<<"\"Avg Frac Dead (PaperAnalysis)\""<<'\t'<<"\"Avg Frac Alive (PaperAnalysis)\""<<'\t';
            myfile<<"\"Average Dead Cells\""<<'\t'<<"\"Average Alive Cells\""<<endl;
        }


        float fracAlive=0.0, fracDead=0.0, fracAlive_PaperAnalysis=0.0, fracDead_PaperAnalysis=0.0;

        fracAlive_PaperAnalysis=Alive[i]/(Alive[i]+Dead[i]);
        fracDead_PaperAnalysis=Dead[i]/(Alive[i]+Dead[i]);
        float denominator=0.0;
        for(unsigned int d=i;d<Alive.size();d++)
        {
            denominator=denominator+Alive[d]+Dead[d];
        }
        if(fabs(denominator)<1e-6)
        {
            fracAlive=0.0;
            fracDead=0.0;
        }
        else
        {
            fracDead=Dead[i]/denominator;
            fracAlive=1-fracDead;
        }

        FracDeadVector.push_back(fracDead);

        if(i<20)
        {
            myfile<<i<<'\t'<<fracDead<<'\t'<<fracAlive<<'\t';
            myfile<<fracDead_PaperAnalysis<<'\t'<<fracAlive_PaperAnalysis<<'\t';
            myfile<<Dead[i]<<'\t'<<Alive[i]<<endl;
        }
        if(i<7 && ExpFracDead[i]!=0)
        {
            temp_SumSqFracDead+= pow((fracDead-(ExpFracDead[i])) ,2);
            n_AbsSumSq++;
            n++;


            //Checking AbsSumSq- gobbledygookykooky
            if(ExpFracDead[i]>1e-6)
            {
                MMH_Temp_AbsSumSq=MMH_Temp_AbsSumSq+((pow((fracDead-ExpFracDead[i]),2))/pow((ExpFracDead[i]) ,2));
                MMH_n_Abs++;
                MMH_n++;
            }
        }
    }

    if(AfterAT==0)
    {
        Results_SumSq[x].FracDeadCells=float(temp_SumSqFracDead/float(n));
        Results_SumSq[x].MMH_AbsoluteSumOfSq+=MMH_Temp_AbsSumSq;
        Results_SumSq[x].MMH_FracDeadCells=float(MMH_Temp_AbsSumSq/float(MMH_n));
        Results_SumSq[x].AbsoluteSumOfSq+=temp_SumSqFracDead;

        Observation.InfectedAlive.clear();
        Observation.InfectedDead.clear();

    }
    else if(AfterAT==1)
    {
        Results_SumSq[x].FracDeadCells_AfterAT=float(temp_SumSqFracDead/float(n));
        Results_SumSq[x].MMH_AbsoluteSumOfSq_AfterAT+=MMH_Temp_AbsSumSq;
        Results_SumSq[x].MMH_FracDeadCells_AfterAT=float(MMH_Temp_AbsSumSq/float(MMH_n));
        Results_SumSq[x].AbsoluteSumOfSq_AfterAT+=temp_SumSqFracDead;

        Observation.InfectedAliveAfterAT.clear();
        Observation.InfectedDeadAfterAT.clear();

    }
    else if(AfterAT==-1)
    {
        Observation.InfectedAlive0_240.clear();
        Observation.InfectedDead0_240.clear();
    }
    else if(AfterAT==-2)
    {
        Observation.InfectedAlive0_360.clear();
        Observation.InfectedDead0_360.clear();
    }
    else if(AfterAT==-3)
    {
        Observation.InfectedAlive0_480.clear();
        Observation.InfectedDead0_480.clear();
    }

    myfile.close();


}
















void CalcOutputVector(vector<float> &AvgOutputVector, vector<float> &StdDevOutputVector, vector<vector<float>> DataVector, sim_parameters & ParValue)
{

    for (unsigned int i=0;i<DataVector[0].size();i++)
    {
        vector<float>TempVector;

        for(unsigned int j=0;j<ParValue.SimPerParSet;j++)
        {
            float TempValue=DataVector[j][i];
            TempVector.push_back(TempValue);
        }
        float TempAvg=0.0, TempStdDev=0.0;
        TempAvg=CalcAvg(TempVector);
        AvgOutputVector.push_back(TempAvg);


        TempStdDev=CalcStdDev(TempVector);
        StdDevOutputVector.push_back(TempStdDev);

    }
}



/*
 * Based on the range of the parameter that is being altered among each simulation set and how many parameter sets' simulations are completed
 * the value of the variable parameter is computed.
 * */
void ParametersChanging(sim_parameters &ParValue, float TempRatio, float vmin, float vmax)
{
    int VP=ParValue.VariableParameter;

    float NewValue=vmin+(vmax-vmin)*TempRatio;

    ParValue.ChangePar(VP, NewValue);

    ParValue.ValChangingVariable=NewValue;

}

/*
 * A short .txt report is generated each time the program runs
 * Gives details of all the parameter values, whether multiple contacts are allowed or not etc
 * */
void WriteReport(string fn, sim_parameters & ParValue, DefaultParameterRange ParRange)
{

    std::string appAddress;         // Address of application (ending with "/")
    appAddress = currentDir();

    appAddress.append("/"+fn+"/");



    ofstream myfile(appAddress+"AAA_Report"+string(".txt"));
    //fn.open(appAddress.c_str());
    myfile<<"Report for "<<fn<<endl;

    myfile<<endl;
    myfile<<"If hypothesis is null or linearly increasing probability of death, then ";
    if(ParValue.DeathDecisionAtStart==0)
        myfile<<"decision to die is taken at end of contact."<<endl;
    else if(ParValue.DeathDecisionAtStart==1)
        myfile<<"decision to die is taken at start of contact."<<endl;


    if(ParValue.InteractDuringApoptosis==0)
        myfile<<"No interactions can take place when infected cell is in apoptotic state. (No Zombie Contacts)"<<endl;
    else if(ParValue.InteractDuringApoptosis==1)
        myfile<<"Interactions can even take place when infected cell is in apoptotic state but they will not contribute to an decisions being taken. (Zombie Contacts)"<<endl;
    myfile<<endl<<endl;

    myfile<<"Hypothesis being tested= "<<ParValue.Hypothesis<<endl<<endl;

    myfile<<"Hypothesis 1: Each cell has equal probability of death during each interaction."<<endl;
    myfile<<"Hypothesis 2: Exponential damage dependent on duration of interaction."<<endl;
    myfile<<"Hypothesis 3: Exponential damage and linear repair."<<endl;
    myfile<<"Hypothesis 4: CTLs become more potent with increasing number of contacts."<<endl;
    myfile<<"Hypothesis 5: Linear damage"<<endl;
    myfile<<"Hypothesis 6: Increasing probability of death with increasing contacts."<<endl;
    myfile<<"Hypothesis 7: Linear damage and repair."<<endl;
    myfile<<"Hypothesis 8: Heterogenous killing capacity (Gaussian) of CTLs."<<endl<<endl;

    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;


    myfile<<"The parameter(s) being varied is/are: ";//----------------------------------change

    if(ParValue.SimulationSets!=1)
    {
    int TempVar=ParValue.VariableParameter;
    myfile<<NameParameter(ParValue.VariableParameter)<<" in the following range ";
    myfile<<ParRange.MinValues[TempVar]<<" to "<<ParRange.MaxValues[TempVar]<<"."<<endl<<endl;
    }
    else
    {
        myfile<<"Only one simulation set. The parameter is not varied."<<endl<<endl;
    }



    myfile<<"Number of simulations per parameter set are "<<ParValue.SimPerParSet<<endl;
    myfile<<"Number of parameter sets are "<<ParValue.SimulationSets<<endl;
    myfile<<"The time bin for all observation sets is "<<ParValue.TimeBin<<endl;
    myfile<<"The time bin for dataset 4I is "<<ParValue.TimeBin_4I<<endl;

    myfile<<"NOTE: If a parameter is being varied and has been mentioned above and the same parameter's value is given below as a constant then ignore the constant value."<<endl<<endl;

    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;

    string DataDist[7]= {"Fixed", "Normal", "LogNormal", "Values taken from Data", "BiModal", "Exponential", "NBDIstribs"};



    myfile<<"Speed Values of CTLs:  "<<DataDist[int (ParValue.velocity_mag_Distribution)]<<endl;

    switch(int (ParValue.velocity_mag_Distribution)){
    case Fixed:{myfile<<"Speed value: "<<ParValue.velocity_mag_Mean<<endl;break;}

    case Normal:{myfile<<"Mean= "<<ParValue.velocity_mag_Mean<<endl<<"SD= "<<ParValue.velocity_mag_SD<<endl; break;}
    case LogNormal:{myfile<<"How to get these values??"<<endl;break;}

    case BiModal:{myfile<<"How to get these values??"<<endl; break;}
    case Exponential:{myfile<<"How to get these values??"<<endl; break;}

    }

    myfile<<"hdaouhfoahfaofaof"<<endl;
    myfile<<ParValue.DistributionSpeed.print()<<endl;
    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;



    myfile<<"Turning angles:  "<<DataDist[int (ParValue.TurningAngleDistribution)]<<endl;

    switch(int (ParValue.TurningAngleDistribution)){
    case Fixed:{myfile<<"Turning angle: "<<ParValue.TurningAngleMean<<endl;break;}

    case Normal:{myfile<<"Mean= "<<ParValue.TurningAngleMean<<endl<<"SD= "<<ParValue.TurningAngleSD<<endl; break;}
    case LogNormal:{myfile<<"How to get these values??"<<endl;break;}

    case BiModal:{myfile<<"How to get these values??"<<endl; break;}
    case Exponential:{myfile<<"How to get these values??"<<endl; break;}

    }
    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;



    myfile<<"Time taken for infected cell to die:  "<<DataDist[int (ParValue.TimetoDieDistribution)]<<endl;

    switch(int (ParValue.TimetoDieDistribution)){
    case Fixed:{myfile<<"Time taken to die: "<<ParValue.TimetoDieMean<<endl;break;}

    case Normal:{myfile<<"Mean= "<<ParValue.TimetoDieMean<<endl<<"SD= "<<ParValue.TimetoDieSD<<endl; break;}
    case LogNormal:{myfile<<"How to get these values??"<<endl;break;}

    case BiModal:{myfile<<"How to get these values??"<<endl; break;}
    case Exponential:{myfile<<"How to get these values??"<<endl; break;}

    }
    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;



    myfile<<"Duration of interactions:  "<<DataDist[int (ParValue.InteractTimeDistribution)]<<endl;

    switch(int (ParValue.InteractTimeDistribution)){
    case Fixed:{myfile<<"Duration: "<<ParValue.InteractTimeMean<<endl;break;}

    case Normal:{myfile<<"Mean= "<<ParValue.InteractTimeMean<<endl<<"SD= "<<ParValue.InteractTimeSD<<endl; break;}
    case LogNormal:{myfile<<"How to get these values??"<<endl;break;}

    case BiModal:{myfile<<"How to get these values??"<<endl; break;}
    case Exponential:{myfile<<"How to get these values??"<<endl; break;}

    }
    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;



    myfile<<"Persistent time:  "<<DataDist[int (ParValue.PersistentTimeDistribution)]<<endl;

    switch(int (ParValue.PersistentTimeDistribution)){
    case Fixed:{myfile<<"Persistent time: "<<ParValue.PersistentTimeMean<<endl;break;}

    case Normal:{myfile<<"Mean= "<<ParValue.PersistentTimeMean<<endl<<"SD= "<<ParValue.PersistentTimeSD<<endl; break;}
    case LogNormal:{myfile<<"How to get these values??"<<endl;break;}

    case BiModal:{myfile<<"How to get these values??"<<endl; break;}
    case Exponential:{myfile<<"How to get these values??"<<endl; break;}

    }

    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;

    myfile<<"NOTE: ";


    if(ParValue.MultipleCTLSingleInfected==0)
        myfile<<"One infected cell can interact with only one CTl at a vertain time point"<<endl;
    else if(ParValue.MultipleCTLSingleInfected==1)
        myfile<<"One infected cell can interact with multiple CTLs but not vice versa"<<endl;

    myfile<<endl;
    myfile<<"If hypothesis is null or linearly increasing probability of death, then ";
    if(ParValue.DeathDecisionAtStart==0)
        myfile<<"decision to die is taken at end of contact."<<endl;
    else if(ParValue.DeathDecisionAtStart==1)
        myfile<<"decision to die is taken at start of contact."<<endl;


    if(ParValue.InteractDuringApoptosis==0)
        myfile<<"No interactions can take place when infected cell is in apoptotic state. (No Zombie Contacts)"<<endl;
    else if(ParValue.InteractDuringApoptosis==1)
        myfile<<"Interactions can even take place when infected cell is in apoptotic state but they will not contribute to an decisions being taken. (Zombie Contacts)"<<endl;


    if(ParValue.LongContactsControlledDamage==0)
        myfile<<"Long contacts keep on imparting damage throughout the interaction."<<endl;
    else if(ParValue.LongContactsControlledDamage==1)
        myfile<<"After a certain point (threshold), the long contacts do not impart damage. Tha value of that threshold is: "<<ParValue.ThresholdContactDurDamage<<endl;


    myfile<<"X dimension= "<<ParValue.box_X<<endl;
    myfile<<"Y dimension= "<<ParValue.box_Y<<endl;
    myfile<<"Z dimension= "<<ParValue.box_Z<<endl;
    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;

    myfile<<"In order to closely mimic the real system, we don't analyze the whole 3 dimensional space that is defined above. ";
    myfile<<"Instead we define a range of values that define the volume that we wish to analyze."<<endl;
    myfile<<"The lower limit for x dimension that is used to define volume to be analyzed= "<<ParValue.X_LowerLim<<endl;
    myfile<<"Range along x dimension that will be analyzed= "<<ParValue.X_AnalyzeRange<<endl;
    myfile<<"The lower limit for y dimension that is used to define volume to be analyzed= "<<ParValue.Y_LowerLim<<endl;
    myfile<<"Range along y dimension that will be analyzed= "<<ParValue.Y_AnalyzeRange<<endl;
    myfile<<"The lower limit for z dimension that is used to define volume to be analyzed= "<<ParValue.Z_LowerLim<<endl;
    myfile<<"Range along z dimension that will be analyzed= "<<ParValue.Z_AnalyzeRange<<endl;
    for(int i=0;i<150;i++)
        myfile<<"-";
    myfile<<endl<<endl;


    myfile<<"The confinement index of infected cells= "<<ParValue.factorConfinementInfected<<endl;


    myfile<<"Total time= "<<ParValue.SimulationTotalTime<<endl;
    myfile<<"Analysis starts after "<<ParValue.AnalyzeAfter<<endl;
    myfile<<"Time Step= "<<ParValue.dT<<endl;
    myfile<<"Number of CTLs= "<<ParValue.num_CTL<<endl;
    myfile<<"Number of stromal cells= "<<ParValue.num_somatic<<endl;


   // myfile<<"Minimum time taken for stromal cell to die (non equal times of death)= "<<ParValue.TimetoDieMin<<endl;
   // myfile<<"Maximum time taken for stromal cell to die (non equal times of death)= "<<ParValue.TimetoDieSD<<endl;
   // myfile<<"Time taken for stromal cell to die= "<<ParValue.TimetoDieMean<<endl;
   // myfile<<"Magnitude of velocity= "<<ParValue.velocity_mag<<endl;
    myfile<<"Stromal cell radius= "<<ParValue.SomaticRadius<<endl;
    myfile<<"CTL cell radius= "<<ParValue.CTLRadius<<endl;




    myfile<<"Hyp1: Probability to die during interaction= "<<ParValue.ProbDie_Interact_Hyp1<<endl;
    myfile<<"Hyp2: Alpha in damage for hypothesis 2. Function for damage: (1-exp(-alpha*t))= "<<ParValue.AlphaDamage_Hyp2<<endl;
    myfile<<"Hyp3: Alpha in damage for hypothesis 3. Function for damage: (1-exp(-alpha*t))= "<<ParValue.AlphaDamage_Hyp3<<endl;
    myfile<<"Hyp3: Beta in repair for hypothesis 3. Function for repair: (1-exp(-beta*t))= "<<ParValue.BetaRepair_Hyp3<<endl;
    myfile<<"Hyp4: Alpha in damage for hypothesis 4. Function for damage: (1-exp(-alpha*Number of prev contacts*t))= "<<ParValue.AlphaDamage_Hyp4<<endl;
    myfile<<"Hyp5: Alpha in damage for hypothesis 5. Function for damage: (1-exp(-alpha*Number of prev contacts*t))= "<<ParValue.AlphaDamage_Hyp5<<endl;
    myfile<<"Hyp6: Gamma factor for hypothesis 6. Function to get probability of cell death: (gamma*No. of prior contacts)= "<<ParValue.GammaFactor_Hyp6<<endl;
    myfile<<"Hyp7: Alpha in damage for hypothesis 7. Function for damage: (1-exp(-alpha*t))= "<<ParValue.AlphaDamage_Hyp7<<endl;
    myfile<<"Hyp7: Beta in repair for hypothesis 7. Function for repair: (1-exp(-beta*t))= "<<ParValue.BetaRepair_Hyp7<<endl;
    myfile<<"Hyp8: Killing factor mean gamma*x= "<<ParValue.KillingFactorMean_Hyp8<<endl;
    myfile<<"Hyp8: Killing factor SD gamma*x= "<<ParValue.KillingFactorSD_Hyp8<<endl;

    myfile.close();
}



/*
 * Generates script in R whcih when runs fives all the plots
 * */
void GenerateRScript(string fn, sim_parameters & ParValue)
{
    std::string appAddress;         // Address of application (ending with "/")
    appAddress = currentDir();
    appAddress.append("/"+fn+"/");

    ofstream myfile(appAddress+"AAA_R_code"+string(".R"));
    //fn.open(appAddress.c_str());
    myfile<<"#R script for "<<fn<<endl<<endl;



    myfile<<"TotPlots<-"<<ParValue.SimulationSets<<endl;
    myfile<<"SimPerParSet<-"<<ParValue.SimPerParSet<<endl;
    myfile<<"address<-\""<<appAddress<<"\""<<endl;

    ifstream commonRcode(appAddress + "../../Cytokill/R_code.R");
    if(!commonRcode) {
        cerr << "ERR: GenerateRScript, the common R code (R_code.R) was not found" <<appAddress<< "/../../Cytokill/R_code.R" << endl;
    } else {
        stringstream buffer;
        buffer << commonRcode.rdbuf();
        myfile << buffer.str();
    }
}
