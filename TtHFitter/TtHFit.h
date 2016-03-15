#include "TtHFitter/Common.h"

#include "TGaxis.h"
#include "TPad.h"
#include "TPie.h"
#include "TF1.h"

#include "TtHFitter/TthPlot.h"
#include "TtHFitter/FitResults.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/NormFactor.h" // KERIM
#include "TtHFitter/Region.h"
#include "TtHFitter/ConfigParser.h"

#ifndef __TtHFit__
#define __TtHFit__

class Region;
class Sample;
class Systematic;
class NormFactor; // KERIM
class ConfigParser;
class RooDataSet;
class RooWorkspace;

class TtHFit {
public:
    
    enum FitType {
        SPLUSB = 1,
        BONLY = 2
    };
    
    enum FitRegion {
        CRONLY = 1,
        CRSR = 2,
        USERSPECIFIC = 3
    };
    
    enum InputType {
        HIST = 0,
        NTUP = 1
    };
    
    enum LimitType {
        ASYMPTOTIC = 0,
        TOYS = 1
    };
    
    TtHFit(string name="MyMeasurement");
    ~TtHFit();
    
    void SetPOI(string name="SigXsecOverSM");
    void SetStatErrorConfig(bool useIt=true, float thres=0.05, string fStatErrCons="Gaussian");
    void SetLumiErr(float err);
    void SetLumi(const float lumi);
    void SetFitType(FitType type);
    void SetLimitType( LimitType type );
    std::string CheckName( const std::string &name );
    void SetFitRegion(FitRegion region);
    
    Sample* NewSample(string name,int type=0);
    Systematic* NewSystematic(string name);
    Region* NewRegion(string name);
    
    // ntuple stuff
    void AddNtuplePath(string path);
    void SetMCweight(string weight);
    void SetSelection(string selection);
    void SetNtupleName(string name);
    void SetNtupleFile(string name);
    void ComputeBining(int regIter);
    
    // histogram stuff
    void AddHistoPath(string path);
    void SetHistoName(string name);
    
    void SmoothSystematics(string syst="all");
    
    // create new root file with all the histograms
    void WriteHistos(/*string fileName=""*/);
    
    void DrawSystPlots();
    void DrawSystPlotsSumSamples();

    // config file
    void ReadConfigFile(string fileName,string options="");
    
    // read from ..
    void ReadNtuples();
    void ReadHistograms();
    void ReadHistos(/*string fileName=""*/);
    
    void DrawAndSaveAll(string opt="");
   
    // separation plots
    void DrawAndSaveSeparationPlots();
    
    TthPlot* DrawSummary(string opt="");
    void BuildYieldTable(string opt="");
    
    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols=0,int nRows=0);
    void DrawSignalRegionsPlot(int nRows,int nCols, std::vector < Region* > &regions);
    void DrawPieChartPlot(const std::string &opt="", int nCols=0,int nRows=0);
    void DrawPieChartPlot(const std::string &opt, int nCols,int nRows, std::vector < Region* > &regions);
    
    // turn to RooStat::HistFactory
    void ToRooStat(bool createWorkspace=true, bool exportOnly=true);
    
    void DrawPruningPlot();
    
    // fit etc...
    void Fit();
    RooDataSet* DumpData( RooWorkspace *ws, std::map < std::string, int > &regionDataType, std::map < std::string, double > &npValues, const double poiValue);
    std::map < std::string, double > PerformFit( RooWorkspace *ws, std::vector < std::string > &regionsToFit, RooDataSet* inputData, FitType fitType=SPLUSB );
    RooWorkspace* PerformWorkspaceCombination( std::vector < std::string > &regionsToFit );

    void PlotFittedNP();
    void PlotCorrelationMatrix();
    void GetLimit();
    void GetSignificance();
    void GetLikelihoodScan( RooWorkspace *ws, string varName, RooDataSet* data);
    
    // get fit results from txt file
    void ReadFitResults(string fileName);
    
    void Print();
    
    Region* GetRegion(string name);
    Sample* GetSample(string name);
    
    void ProduceNPRanking(string NPnames="all");
    void PlotNPRanking();
    
    void PrintSystTables();
    
    // -------------------------
      
    string fName;
    string fLabel;
    string fResultsFolder;
    
    std::vector < Region* > fRegions;
    std::vector < Sample* > fSamples;
    std::vector < Systematic* > fSystematics;
    std::vector < NormFactor* > fNormFactors;
    std::vector < string > fSystematicNames;
    std::vector < string > fNormFactorNames;
    
    int fNRegions;
    int fNSamples;
    int fNSyst;
    int fNNorm;
    string fPOI;
    bool fUseStatErr;
    float fStatErrThres;
    string fStatErrCons;
    
    float fLumi;
    float fLumiScale;
    float fLumiErr;
    
    float fThresholdSystPruning_Normalisation;
    float fThresholdSystPruning_Shape;
    float fThresholdSystLarge;
    std::vector<string> fNtuplePaths;
    string fNtupleFile;
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
    bool fSystDataPlot_upFrame;
    bool fStatOnly;
    
    std::vector<string> fRegionsToPlot;
    std::vector<string> fSummaryPlotRegions;
    std::vector<string> fSummaryPlotLabels;
    std::vector<string> fSummaryPlotValidationRegions;
    std::vector<string> fSummaryPlotValidationLabels;

//     bool fHistoCheckCrash;
    
    string fLumiLabel;
    string fCmeLabel;
    
    string fSuffix;
    string fSaveSuf;
    string fLoadSuf;
    
    bool fUpdate;
    
    float fBlindingThreshold;
    
    int fRankingMaxNP;
    std::string fRankingOnly;
    
    std::string fImageFormat;
    
    //
    // Fit caracteristics
    //
    FitType fFitType;
    FitRegion fFitRegion;
    std::vector< std::string > fFitRegionsToFit;
    std::map< std::string, double > fFitNPValues;
    double fFitPOIAsimov;
    bool fFitIsBlind;
    vector<string> fVarNameLH;
    vector<string> fVarNameMinos;

    //
    // Limit parameters
    //
    LimitType fLimitType;
    bool fLimitIsBlind;
    double fLimitPOIAsimov;
};

#endif
