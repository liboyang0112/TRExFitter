#include "TtHFitter/Common.h"

#include "TtHFitter/TtHFit.h"

#ifndef __MultiFit__
#define __MultiFit__

class RooWorkspace;

class MultiFit {
public:

    MultiFit(string name="MyMultiFit");
    ~MultiFit();

    void ReadConfigFile(string configFile,string options);
    void AddFitFromConfig(string configFile,string options,string label,string loadSuf="",string wsFile="");
    RooWorkspace* CombineWS();
    void SaveCombinedWS();
    std::map < std::string, double > FitCombinedWS( int fitType=1, string inputData="" );
    void GetCombinedLimit(string inputData="obsData"); // or asimovData
    void GetCombinedSignificance(string inputData="obsData"); // or asimovData

    void ComparePOI(string POI);
    void CompareLimit();
    void ComparePulls(string caterogy="");
    void CompareNormFactors(string category="");
    void PlotCombinedCorrelationMatrix();
    void ProduceNPRanking(string NPnames="all");
    void PlotNPRankingManager();
    void PlotNPRanking(bool flagSysts=true, bool flagGammas=false);
    void PlotSummarySoverB();    
    void GetLikelihoodScan( RooWorkspace *ws, string varName, RooDataSet* data,bool recreate=true,bool compare=false);

    TH1F* Combine(std::vector<TH1F*>);
    TH1F* OrderBins(TH1F* h,std::vector<float> vec);
    TH1F* Rebin(TH1F* h,std::vector<float> vec, bool isData=true);

    std::vector< string > fFitNames;
    std::vector< TtHFit* > fFitList;
    std::vector< string > fFitLabels;
    std::vector< string > fFitSuffs;
    std::vector< string > fWsFiles;
    std::vector< string > fDirectory;
    std::vector< string > fInputName;

    std::vector< string > fNPCategories;

    bool fCombine;
    bool fCompare;
    bool fStatOnly;
    bool fIncludeStatOnly;

    bool fCompareLimits;
    bool fComparePOI;
    bool fComparePulls;
    bool fPlotCombCorrMatrix;

    string fName;
    string fDir;
    string fOutDir;
    string fLabel;
    bool fShowObserved;
    string fLimitTitle;
    string fPOITitle;
    string fRankingOnly;

    string fPOI;
    float fPOIMin;
    float fPOIMax;
    float fPOIVal;
    string fPOIPrecision;
    float fLimitMax;

    bool fUseRnd;
    float fRndRange;
    long int fRndSeed;

    string fLumiLabel;
    string fCmeLabel;

    ConfigParser *fConfig;

    string fSaveSuf;
    std::vector< bool > fFitShowObserved;

    string fDataName;
    int fFitType;
    bool fSignalInjection;

    bool fCombineChByCh;
    bool fFastFit;
    bool fFastFitForRanking;
    string fNuisParListFile;

    bool fPlotSoverB;
    string fSignalTitle;

    std::string fFitResultsFile;
    std::string fLimitsFile;
    std::vector<std::string> fLimitsFiles;
    std::string fBonlySuffix;

    bool fShowSystForPOI;
    bool fGetGoodnessOfFit;
    
    std::vector<std::string> fVarNameLH;
};

#endif
