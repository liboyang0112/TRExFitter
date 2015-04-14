#include "TtHFitter/Common.h"

#include "TtHFitter/TthPlot.h"
#include "TtHFitter/FitResults.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/Region.h"

#ifndef __TtHFit__
#define __TtHFit__

class Region;
class Sample;
class Systematic;

class TtHFit {
public:
    
    enum FitType {
        ControlRegion = 1,
        ControlSignalRegion = 2
    };
    
    TtHFit(string name="MyMeasurement");
    ~TtHFit();
    
    void SetPOI(string name="SigXsecOverSM");
    void SetStatErrorConfig(bool useIt=true, float thres=0.05, string fStatErrCons="Gaussian");
    void SetLumiErr(float err);
    void SetLumi(const float lumi);
    void SetFitType(FitType type);
    
    Sample* NewSample(string name,int type=0);
    Systematic* NewSystematic(string name);
    Region* NewRegion(string name);
    
    // ntuple stuff
    void AddNtuplePath(string path);
    void SetMCweight(string weight);
    void SetSelection(string selection);
    void SetNtupleName(string name);
    
    // histogram stuff
    void AddHistoPath(string path);
    void SetHistoName(string name);
    
    // systematic handle
    void SmoothSystematics(string syst="all");
    
    // create new root file with all the histograms
    void WriteHistos(string fileName="",bool recreate=true);
    
    // read from ..
    void ReadNtuples();
    void ReadHistograms();
    void ReadHistos(string fileName="");
    void ReadAll(bool readNtuples=true,string fileName="");
    
    void DrawAndSaveAll(string opt="");
    
    TthPlot* DrawSummary(string opt="");
    
    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols,int nRows);
    void DrawSignalRegionsPlot(int nRows,int nCols, std::vector < Region* > &regions);
    void DrawPieChartPlot();
    
    // turn to RooStat::HistFactory
    void ToRooStat(bool createWorkspace=true, bool exportOnly=true);
    
    // fit etc...
    void Fit();
    void PlotFittedNP();
    void GetLimit();
    
    // get fit results from txt file
    void ReadFitResults(string fileName);
    
    void Print();
    
    string fName;
    string fResultsFolder;
    FitType fFitType;
    
    std::vector < Region* > fRegions;
    std::vector < Sample* > fSamples;
    std::vector < Systematic* > fSystematics;
    
    int fNRegions;
    int fNSamples;
    int fNSyst;
    string fPOI;
    bool fUseStatErr;
    float fStatErrThres;
    string fStatErrCons;
    
    float fLumi;
    float fLumiErr;
    
    float fThresholdSystPruning_Normalisation;
    float fThresholdSystPruning_Shape;
    
    vector<string> fNtuplePaths;
    string fMCweight;
    string fSelection;
    string fNtupleName;
    
    vector<string> fHistoPaths;
    string fHistoName;
    
    FitResults *fFitResults;
};

#endif
