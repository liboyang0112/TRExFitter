#include "TtHFitter/Common.h"

#include "TtHFitter/TthPlot.h"
#include "TtHFitter/FitResults.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/Region.h"
#include "TtHFitter/ConfigParser.h"

#ifndef __TtHFit__
#define __TtHFit__

class Region;
class Sample;
class Systematic;
class ConfigParser;

class TtHFit {
public:
    
    enum FitType {
        SPLUSB = 1,
        BONLY = 2
    };
    
    enum FitRegion {
        CRONLY = 1,
        CRSR = 2
    };
    
    enum InputType {
        HIST = 0,
        NTUP = 1
    };
    
    TtHFit(string name="MyMeasurement");
    ~TtHFit();
    
    void SetPOI(string name="SigXsecOverSM");
    void SetStatErrorConfig(bool useIt=true, float thres=0.05, string fStatErrCons="Gaussian");
    void SetLumiErr(float err);
    void SetLumi(const float lumi);
    void SetFitType(FitType type);
    void SetFitRegion(FitRegion region);
    
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
    
    void SmoothSystematics(string syst="all");
    
    // create new root file with all the histograms
    void WriteHistos(string fileName="");
    
    void DrawSystPlots();
    
    // config file
    void ReadConfigFile(string fileName,string options="");
    
    // read from ..
    void ReadNtuples();
    void ReadHistograms();
    void ReadHistos(string fileName="");
    
    void DrawAndSaveAll(string opt="");
    
    TthPlot* DrawSummary(string opt="");
    void BuildYieldTable(string opt="");
    
    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols=0,int nRows=0);
    void DrawSignalRegionsPlot(int nRows,int nCols, std::vector < Region* > &regions);
    void DrawPieChartPlot();
    
    // turn to RooStat::HistFactory
    void ToRooStat(bool createWorkspace=true, bool exportOnly=true);
    
    void DrawPruningPlot();
    
    // fit etc...
    void Fit();
    void PlotFittedNP();
    void PlotCorrelationMatrix();
    void GetLimit();
    void GetSignificance();
    
    // get fit results from txt file
    void ReadFitResults(string fileName);
    
    void Print();
    
    Region* GetRegion(string name);
    
    string fName;
    string fLabel;
    string fResultsFolder;
    
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
    
    std::vector<string> fNtuplePaths;
    string fMCweight;
    string fSelection;
    string fNtupleName;
    
    std::vector<string> fHistoPaths;
    string fHistoName;
    string fHistoFile;
    
    FitResults *fFitResults;
    
    int fIntCode_overall;
    int fIntCode_shape;
    
    int fInputType; // 0: histo, 1: ntup
    
    ConfigParser *fConfig;
    
    bool fSystControlPlots;
    
    std::vector<string> fRegionsToPlot;

    bool fHistoCheckCrash;
    
    string fLumiLabel;
    string fCmeLabel;
    
    string fSuffix;
    string fSaveSuf;
    string fLoadSuf;
    
    bool fUpdate;
    
    //
    // Fit caracteristics
    //
    FitType fFitType;
    FitRegion fFitRegion;
    std::vector< std::string > fFitRegionsToFit;
    std::map< std::string, double > fFitNPValues;
    double fFitPOIAsimov;
    bool fFitIsBlind;
};

#endif
