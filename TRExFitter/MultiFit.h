#ifndef MULTIFIT_H
#define MULTIFIT_H

/// Framework includes
#include "TRExFitter/Common.h"

/// c++ includes
#include <map>
#include <memory>

/// Forwards declaration
class ConfigParser;
class RooDataSet;
class RooWorkspace;
class TH1D;
class TRExFit;

class MultiFit {
public:

    explicit MultiFit(std::string name="MyMultiFit");

    ~MultiFit() = default;

    void AddFitFromConfig(const std::string& configFile, const std::string& opt, const std::string& options,
                          const std::string& label, std::string loadSuf="",std::string wsFile="");
    RooWorkspace* CombineWS() const;
    void SaveCombinedWS() const;
    std::map < std::string, double > FitCombinedWS( int fitType, const std::string& inputData, bool doLHscanOnly ) const;
    void GetCombinedLimit(std::string inputData="obsData") const; // or asimovData
    void GetCombinedSignificance(std::string inputData="obsData") const; // or asimovData

    void ComparePOI(const std::string& POI) const;
    void CompareLimit();
    void ComparePulls(std::string category="") const;
    void CompareNormFactors(std::string category="") const;
    void PlotCombinedCorrelationMatrix() const;
    void ProduceNPRanking(std::string NPnames="all") const;
    void PlotNPRankingManager() const;
    void PlotNPRanking(bool flagSysts=true, bool flagGammas=false) const;
    void PlotSummarySoverB() const;
    void GetLikelihoodScan( RooWorkspace *ws, const std::string& varName, RooDataSet* data,bool recreate=true) const;
    void Get2DLikelihoodScan( RooWorkspace *ws, const std::vector<std::string>& varName, RooDataSet* data) const;
    void BuildGroupedImpactTable() const;

    TH1D* Combine(std::vector<TH1D*> hists) const;
    TH1D* OrderBins(TH1D* h, std::vector<double> vec) const;
    TH1D* Rebin(TH1D* h, const std::vector<double>& vec, bool isData=true) const;

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
    double fPOIMin;
    double fPOIMax;
    double fPOIVal;
    std::string fPOIPrecision;
    double fLimitMax;

    bool fUseRnd;
    double fRndRange;
    long int fRndSeed;

    std::string fLumiLabel;
    std::string fCmeLabel;

    std::unique_ptr<ConfigParser> fConfig;

    std::string fSaveSuf;
    std::vector< bool > fFitShowObserved;

    std::string fDataName;
    int fFitType;

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
    std::vector<std::vector<std::string> > fVarName2DLH;
    double fLHscanMin;
    double fLHscanMax;
    int fLHscanSteps;
    double fLHscanMinY;
    double fLHscanMaxY;
    int fLHscanStepsY;
    bool fParal2D;
    int fParal2Dstep;
    bool fDoGroupedSystImpactTable;

    std::string fPOIName;
    double fPOINominal;

    //
    // Limit parameters
    //
    bool fLimitIsBlind;
    double fLimitPOIAsimov;
    bool fSignalInjection;
    double fSignalInjectionValue;
    std::string fLimitParamName;
    double fLimitParamValue;
    std::string fLimitOutputPrefixName;
    double fLimitsConfidence;

    //
    // Significance parameters
    //
    bool fSignificanceIsBlind;
    double fSignificancePOIAsimov;
    std::string fSignificanceParamName;
    double fSignificanceParamValue;
    std::string fSignificanceOutputPrefixName;

    bool fShowTotalOnly;

    bool fuseGammasForCorr;

    double fPOIInitial;
    std::vector<std::string> fConfigPaths;
};

#endif
