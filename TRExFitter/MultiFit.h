#ifndef MULTIFIT_H
#define MULTIFIT_H

/// Framework includes
#include "TRExFitter/Common.h"

/// c++ includes
#include <map>

/// Forwards declaration
class ConfigParser;
class RooDataSet;
class RooWorkspace;
class TH1F;
class TRExFit;

class MultiFit {
public:

    MultiFit(std::string name="MyMultiFit");
    ~MultiFit();

    void AddFitFromConfig(std::string configFile,std::string options,std::string label,std::string loadSuf="",std::string wsFile="");
    RooWorkspace* CombineWS();
    void SaveCombinedWS();
    std::map < std::string, double > FitCombinedWS( int fitType=1, std::string inputData="", bool performFit=true );
    void GetCombinedLimit(std::string inputData="obsData"); // or asimovData
    void GetCombinedSignificance(std::string inputData="obsData"); // or asimovData

    void ComparePOI(std::string POI);
    void CompareLimit();
    void ComparePulls(std::string caterogy="");
    void CompareNormFactors(std::string category="");
    void PlotCombinedCorrelationMatrix();
    void ProduceNPRanking(std::string NPnames="all");
    void PlotNPRankingManager();
    void PlotNPRanking(bool flagSysts=true, bool flagGammas=false);
    void PlotSummarySoverB();
    void GetLikelihoodScan( RooWorkspace *ws, std::string varName, RooDataSet* data,bool recreate=true,bool compare=false);
    void BuildGroupedImpactTable();

    TH1F* Combine(std::vector<TH1F*>);
    TH1F* OrderBins(TH1F* h,std::vector<float> vec);
    TH1F* Rebin(TH1F* h,std::vector<float> vec, bool isData=true);
    
    std::vector< std::string > fFitNames;
    std::vector< TRExFit* > fFitList;
    std::vector< std::string > fFitLabels;
    std::vector< std::string > fFitSuffs;
    std::vector< std::string > fWsFiles;
    std::vector< std::string > fDirectory;
    std::vector< std::string > fInputName;

    std::vector< std::string > fNPCategories;

    bool fCombine;
    bool fCompare;
    bool fStatOnly;
    bool fIncludeStatOnly;

    bool fCompareLimits;
    bool fComparePOI;
    bool fComparePulls;
    bool fPlotCombCorrMatrix;

    std::string fName;
    std::string fDir;
    std::string fOutDir;
    std::string fLabel;
    bool fShowObserved;
    std::string fLimitTitle;
    std::string fPOITitle;
    std::string fRankingOnly;
    std::string fGroupedImpactCategory;

    std::string fPOI;
    float fPOIMin;
    float fPOIMax;
    float fPOIVal;
    std::string fPOIPrecision;
    float fLimitMax;

    bool fUseRnd;
    float fRndRange;
    long int fRndSeed;

    std::string fLumiLabel;
    std::string fCmeLabel;

    ConfigParser *fConfig;

    std::string fSaveSuf;
    std::vector< bool > fFitShowObserved;

    std::string fDataName;
    int fFitType;
    bool fSignalInjection;

    bool fCombineChByCh;
    bool fFastFit;
    bool fFastFitForRanking;
    std::string fNuisParListFile;

    bool fPlotSoverB;
    std::string fSignalTitle;

    std::string fFitResultsFile;
    std::string fLimitsFile;
    std::vector<std::string> fLimitsFiles;
    std::string fBonlySuffix;

    bool fShowSystForPOI;
    bool fGetGoodnessOfFit;

    std::vector<std::string> fVarNameLH;
    float fLHscanMin; 
    float fLHscanMax;
    int fLHscanSteps; 
    bool fDoGroupedSystImpactTable;
    
    bool fRunROOTMacros;

    std::string fPOIName;
    float fPOINominal;
};

#endif