#ifndef i_distribution
#define i_distribution
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include <sstream>
#include <fstream>
using namespace std;

enum{Fixed, Normal, LogNormal, FromData, BiModal, Exponential, NBDIstribs};

string nameDistrib(int ID);

struct probaLawFromTable{
    int size;
    vector<double> Xs;
    vector<double> NormalizedDensities;     // to probabilities
    vector<double> cumulatedValuesHigh;      // cumulated probabilities
    double normalizingCoeff;                // just for remembering the normalisation

    bool uniformizeInsideClasses;
    vector<double> lowBoundsXs;
    vector<double> highBoundsXs;

    static double IndependentRandRealGen(); // If later someone wants to separate generators for seeding.
    virtual double getRandValue();
    int getRandClass();
    probaLawFromTable(vector<double> &_Xs, vector<double> &densities, bool _uniformizeInsideClasses);
    probaLawFromTable(string filename, bool _uniformizeInsideClasses);
    void init(vector<double> &_Xs, vector<double> &densities, bool _uniformizeInsideClasses);

    virtual string print();
    virtual ~probaLawFromTable(){}
};

struct probaLawStoringResults : public probaLawFromTable
{
    vector<int> nbOfTimesPerClass;
    int totalNbEvents;
    probaLawStoringResults(vector<double> & _Xs, vector<double> & densities, bool _uniformizeInsideClasses);
    double getRandValue();
    double getFrequency(int classIndex);
    int fitDoubleInAClass(double val);
    void clearRecord();
    string print();
    string TestprobaLawFromTable();
    virtual ~probaLawStoringResults(){}
};

// automatic
struct histogramFromDistrib{
    vector<double> lowBoundsXs;
    vector<double> highBoundsXs;
    vector<double> averageXs;
    vector<double> densities;
    double vmin;
    double vmax;
    double maxDens;
    string print(){
        stringstream ss;
        int S = lowBoundsXs.size();
        if((int) highBoundsXs.size() != S) cerr << "Wrong sizes in histogramFromDistrib " << endl;
        if((int) averageXs.size() != S) cerr << "Wrong sizes in histogramFromDistrib " << endl;
        if((int) densities.size() != S) cerr << "Wrong sizes in histogramFromDistrib " << endl;
        for(int i = 0; i < S; ++i){
            ss << averageXs[i] << "\t" << densities[i] << endl;
        }
        return ss.str();
    }
    histogramFromDistrib(vector<double> listValues, int nbBins){
        int S = listValues.size();
        if(S == 0) {cerr << "Histogram: empty data < 1\n"; return;}
        if(nbBins < 1) {cerr << "Histogram: nbBins < 1\n"; return;}
        lowBoundsXs.resize(nbBins, 0.);
        highBoundsXs.resize(nbBins, 0.);
        averageXs.resize(nbBins, 0.);
        densities.resize(nbBins, 0.);

        maxDens = 0;
        vmin = listValues[0];
        vmax = listValues[0];
        for(int i = 0; i < S; ++i){
            vmin = min(vmin, listValues[i]);
            vmax = max(vmax, listValues[i]);
        }
        double interval = (vmax - vmin) / ((double) nbBins - 1);
        if(fabs(interval) < 1e-9) {
            lowBoundsXs.clear();
            lowBoundsXs.resize(1, vmin);
            highBoundsXs.clear();
            highBoundsXs.resize(1, vmax);
            averageXs.clear();
            averageXs.resize(1, 0.5*vmin + 0.5 * vmax );
            densities.clear();
            densities.resize(1, 1.0);
            cerr << "WRN : Histogram: too few stddev\n"; return;}
        for(int i = 0; i < nbBins; ++i){
            lowBoundsXs[i] = vmin + interval * ((double) i);
            highBoundsXs[i] = vmin + interval * ((double) (i + 1));
            averageXs[i] = 0.5*lowBoundsXs[i] + 0.5*highBoundsXs[i];
        }
        for(int i = 0; i < S; ++i){
            double val = listValues[i];
            int bin = (val - vmin) / interval;
            if(bin < 0) bin = 0;
            if(bin >= nbBins) bin = nbBins-1;
            densities[bin] += 1.0 / (double) S;
        }
        for(int i = 0; i < (int) densities.size(); ++i){
            maxDens = max(maxDens, densities[i]);
        }
    }
};

struct Law{
    Law(){
        type = Fixed;
        mu1 = NAN;
        sigma1 = NAN;
        mu2 = NAN;
        sigma2 = NAN;
        weight = NAN;
        currentDistribFromFile = NULL;

        record = false;
        recordedValues.clear();
    }
    void operator=(const Law &toCopy){
        type = toCopy.type;
        mu1 = toCopy.mu1;
        mu2 = toCopy.mu2;
        sigma1 = toCopy.sigma1;
        sigma2 = toCopy.sigma2;
        weight = toCopy.weight;
        cerr << "Copy" << endl;
        if(toCopy.sourceDataFile.size() > 0) currentDistribFromFile = new probaLawFromTable(toCopy.sourceDataFile, true);
        sourceDataFile = toCopy.sourceDataFile;

        record = false;
        recordedValues.clear();
    }

    Law(const Law &toCopy){
        type = toCopy.type;
        mu1 = toCopy.mu1;
        mu2 = toCopy.mu2;
        sigma1 = toCopy.sigma1;
        sigma2 = toCopy.sigma2;
        weight = toCopy.weight;
        //cerr << "Copy2, type=" << nameDistrib(type) << endl;
        if(toCopy.sourceDataFile.size() > 0) currentDistribFromFile = new probaLawFromTable(toCopy.sourceDataFile, true);
        sourceDataFile = toCopy.sourceDataFile;

        record = false;
        recordedValues.clear();
    }

    int type;
    double mu1;
    double sigma1;
    double mu2;
    double sigma2;
    double weight;
    probaLawFromTable* currentDistribFromFile;

    bool record;
    vector<double> recordedValues;
    string sourceDataFile; // in case of 'fromData'

    void parse(ifstream& fileBeingRead){
        string typeDistrib;
        fileBeingRead >> typeDistrib;
        type = -1;
        for(int i = 0; i < NBDIstribs-1; ++i){
            if(!typeDistrib.compare(nameDistrib(i))){
                type = i;
            }
        }
        // if the type is not correct we suppose this is a file name
        if(type < 0){
            set(typeDistrib); // typeDistrib will become sourceDataFile
        }

        int nbArguments = 0;
        switch (type){
            case FromData: nbArguments = 0; break;
            case Fixed: case Exponential: nbArguments = 1; break;
            case Normal: case LogNormal: nbArguments = 2; break;
            case BiModal: nbArguments = 5; break;
            default: cerr << "ERR: Law:parse, problem with finding the correct distribution. Should never happen" << endl;
        }
        if(nbArguments >= 1){fileBeingRead >> mu1;}
        if(nbArguments >= 2){fileBeingRead >> sigma1;}
        if(nbArguments >= 3){fileBeingRead >> mu2;}
        if(nbArguments >= 4){fileBeingRead >> sigma2;}
        if(nbArguments >= 5){fileBeingRead >> weight;}
        resetRecords();
    }


    string print(){
        if(type == FromData) return sourceDataFile;
        stringstream res;
        res << nameDistrib(type) << endl;
        int nbArguments = 0;
        switch (type){
            case FromData: nbArguments = 0; break;
            case Fixed: case Exponential: nbArguments = 1; break;
            case Normal: case LogNormal: nbArguments = 2; break;
            case BiModal: nbArguments = 5; break;
            default: cerr << "ERR: Law:parse, problem with outputting the correct distribution. Should never happen" << endl;
        }
        if(nbArguments >= 1){res << mu1 << "/t";}
        if(nbArguments >= 2){res << sigma1 << "/t";}
        if(nbArguments >= 3){res << mu2 << "/t";}
        if(nbArguments >= 4){res << sigma2 << "/t";}
        if(nbArguments >= 5){res << weight << "/t";}
        return res.str();
    }


    // two ways to initiate a distribution: from a file or from a distribution
    void set(string fileName){
        if(currentDistribFromFile) delete currentDistribFromFile;
        //cout << "Reading distribution file " << fileName << endl;
        currentDistribFromFile = new probaLawFromTable(fileName, true);
        type = FromData;
        sourceDataFile = fileName;
        resetRecords();
    }
    void set(int _type, double _mu1, double _sigma1 = NAN, double _weight = NAN, double _mu2 = NAN, double _sigma2 = NAN){
        type = _type;
        mu1 = _mu1;
        mu2 = _mu2;
        sigma1 = _sigma1;
        sigma2 = _sigma2;
        weight  =_weight;
        resetRecords();
    }

    void setRecord(bool _record){
        record = _record;
        resetRecords();
    }
    void resetRecords(){
        recordedValues.clear();
    }

    // then, calling the distrib
    double getRandValue(){
        double res;
        switch(type){
        case Fixed:{res = mu1; break;}
        case Normal:{res = random::normal(mu1, sigma1); break;}
        case LogNormal:{res = random::logNormal(mu1, sigma1);break;}
        case FromData:{res = (currentDistribFromFile) ? currentDistribFromFile->getRandValue() : NAN; break;}
        case BiModal:{res = random::biModal(mu1, sigma1, weight, mu2, sigma2); break;}
        case Exponential:{res = random::exponential(mu1); break;}
        default:{res = NAN; std::cerr << "ERR: Type of distribution unknown";}
        }
        if(record) recordedValues.push_back(res);
        return res;
    }

    //std::pair<vector<double>, vector<double>> getHistogram();
    histogramFromDistrib getHistogram(){
       return histogramFromDistrib(recordedValues, recordedValues.size() / 50);
    }
};

#endif
