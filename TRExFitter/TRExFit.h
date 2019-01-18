#ifndef TRExFIT_H
#define TRExFIT_H

/// Framework includes
#include "TRExFitter/HistoTools.h"

// RooFit
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"

/// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forwards class declaration
class ConfigParser;
class FitResults;
class NormFactor;
class RooDataSet;
class RooWorkspace;
class Region;
class Sample;
class SampleHist;
class ShapeFactor;
class Systematic;
class TRExPlot;
class TFile;

class TRExFit {
public:

    enum FitType {
        UNDEFINED = 0,
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

    enum TemplateInterpolationOption{
        LINEAR = 0,
        SMOOTHLINEAR = 1,
        SQUAREROOT = 2
    };

    enum PrunningType{
        SEPARATESAMPLE = 0,
        BACKGROUNDREFERENCE = 1,
        COMBINEDREFERENCE = 2
    };

    struct TemplateWeight{
        std::string function;
        std::string range;
        std::string name;
        float value;
    };

    enum BlindingType {
      SOVERB = 1,
      SOVERSPLUSB = 2,
      SOVERSQRTB = 3,
      SOVERSQRTSPLUSB = 4
    };

    TRExFit(std::string name="MyMeasurement");
    ~TRExFit();

    void SetPOI(std::string name="SigXsecOverSM");
    void SetStatErrorConfig(bool useIt=true, float thres=0.05, std::string cons="Poisson");
    void SetLumiErr(float err);
    void SetLumi(const float lumi);
    void SetFitType(FitType type);
    void SetLimitType( LimitType type );
    void SetFitRegion(FitRegion region);

    Sample* NewSample(const std::string& name,int type=0);
    Systematic* NewSystematic(const std::string& name);
    Region* NewRegion(const std::string& name);

    // ntuple stuff
    void AddNtuplePath(const std::string& path);
    void SetMCweight(const std::string& weight);
    void SetSelection(const std::string& selection);
    void SetNtupleName(const std::string& name);
    void SetNtupleFile(const std::string& name);
    void ComputeBinning(int regIter);
    void DefineVariable(int regIter);

    // histogram stuff
    void AddHistoPath(const std::string& path);
    void SetHistoName(const std::string& name);

    void SmoothSystematics(std::string syst="all");

    // create new root file with all the histograms
    void CreateRootFiles();
    void WriteHistos(bool reWriteOrig=true) const;

    void DrawSystPlots() const;
    void DrawSystPlotsSumSamples() const;

    // read from ..
    void ReadNtuples();
    void ReadHistograms();
    void ReadHistos(/*std::string fileName=""*/);
    void CloseInputFiles();
    void CorrectHistograms();

    void DrawAndSaveAll(std::string opt="");

    // separation plots
    void DrawAndSaveSeparationPlots() const;

    TRExPlot* DrawSummary(std::string opt="", TRExPlot* = 0) ;
    void DrawMergedPlot(std::string opt="",std::string group="") const;
    void BuildYieldTable(std::string opt="",std::string group="") const;

    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols=0,int nRows=0) const;
    void DrawSignalRegionsPlot(int nRows,int nCols, std::vector < Region* > &regions) const;
    void DrawPieChartPlot(const std::string &opt="", int nCols=0,int nRows=0) const;
    void DrawPieChartPlot(const std::string &opt, int nCols,int nRows, std::vector < Region* > &regions) const;

    void CreateCustomAsimov() const;

    // turn to RooStat::HistFactory
    void ToRooStat(bool createWorkspace=true, bool exportOnly=true);

    void SystPruning() const;
    void DrawPruningPlot() const;

    // fit etc...
    void Fit();
    RooDataSet* DumpData( RooWorkspace *ws, std::map < std::string, int > &regionDataType, std::map < std::string, double > &npValues, const double poiValue);
    std::map < std::string, double > PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType=SPLUSB, bool save=false, int debugLevel=1 );
    RooWorkspace* PerformWorkspaceCombination( std::vector < std::string > &regionsToFit ) const;

    void PlotFittedNP();
    void PlotCorrelationMatrix();
    void GetLimit();
    void GetSignificance();
    void GetLikelihoodScan( RooWorkspace *ws, std::string varName, RooDataSet* data) const;

    // get fit results from txt file
    void ReadFitResults(const std::string& fileName);

    void PrintConfigSummary() const;

    Region* GetRegion(const std::string& name) const;
    Sample* GetSample(const std::string& name) const;
    std::size_t GetSampleIndex(const std::string& name) const;

    void ProduceNPRanking(std::string NPnames="all");
    void PlotNPRanking(bool flagSysts=true, bool flagGammas=true) const;
    void PlotNPRankingManager() const;

    void PrintSystTables(std::string opt="") const;

    void MergeSystematics(); // this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter

    // for template fitting
    void AddTemplateWeight(const std::string& name, float);

    std::vector<TemplateWeight> GetTemplateWeightVec(const TemplateInterpolationOption& opt);

    std::string GetWeightFunction(std::vector<std::pair<float,std::string> > templatePair, unsigned int itemp, const TemplateInterpolationOption& opt) const;

    /**
     * Function that returns string that represents smoothed abs value function
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSmoothLinearInterpolation(unsigned int itemp) const;

    /**
     * Helper function to calualte numerical correction to the smoothed linear function
     * @param parameter in the argument of the hyperbolic tangent function
     * @param size of the x axis interval
     * @param central position of the function
     * @param left position of the function on x axis
     * @param right position of the function on x axis
     * @param parameter of iteration, set to 1 for the first iteration
     * @return correction
     */
    double GetCorrection(float k, float width, float x_mean, float x_left, float init = 1) const;

    /**
     * Helper function to approximate absolute value by sqrt(x^2+e)
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSquareRootLinearInterpolation(unsigned int itemp) const;

    /**
     * Helper function to apply correction to square root aproximation
     * @param will return value for a from -a*sqrt(x^2+epsilon) +b
     * @param will return value for b from -a*sqrt(x^2+epsilon) +b
     * @param central position of the function
     * @param left position of the function on x axis
     * @param epsilon = precision of the approximation
     */
    void GetSquareCorrection(double *a, double *b, float x_i, float x_left, float epsilon) const;

    /**
     * Helper function to smooth morphing templates according to a given functional form (bin-by-bin)
     * @param name of the morphing parameter
     * @param functional form
     * @param pointer to an array of parameter values
     */
    void SmoothMorphTemplates(const std::string& name,const std::string& formula="pol1",double *p=0x0) const;

    /**
     * Helper function that draws plots with morphijng templates
     * @param name of the morphing parameter
     */
    void DrawMorphingPlots(const std::string& name) const;

    bool MorphIsAlreadyPresent(const std::string& name, const float value) const;

    // for grouped impact evaluation
    void ProduceSystSubCategoryMap();

    void BuildGroupedImpactTable() const;

    /**
     * Helper function that runs toys experiments
     * @param A pointer to a workspace needed to run the fit
     */
    void RunToys(RooWorkspace* ws);
    
    /**
     * Helper function to compute the variable string to be used when reading ntuples, for a given region, sample combination
     * @param pointer to the Region
     * @param pointer to the Sample
     */
    std::string Variable(Region *reg,Sample *smp);
    
    /**
     * Helper function to compute the selection string to be used when reading ntuples, for a given region, sample combination
     * @param pointer to the Region
     * @param pointer to the Sample
     */
    std::string FullSelection(Region *reg,Sample *smp);
    
    /**
     * Helper function to compute the weight string to be used when reading ntuples, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = NULL)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::string FullWeight(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true);
    
    /**
     * Helper function to compute the full paths to be used when reading ntuples, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = NULL)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::vector<std::string> FullNtuplePaths(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true);
    
    /**
     * Helper function to compute the full paths to be used when reading histograms, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = NULL)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::vector<std::string> FullHistogramPaths(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true);

    /**
     *
     * A helper function to read singleHistogram from a file
     * @param vector of histo paths
     * @param pointer to systematic
     * @param index of channel
     * @param index of sample
     * @param flag if the systematic is up or down
     * @param flag if we are reading MC or Data
     * @return the read histogram
     */ 
    TH1D* ReadSingleHistogram(const std::vector<std::string>& fullPaths, Systematic* syst,
        int i_ch, int i_smp, bool isUp, bool isMC);

    /**
    * A helper function to get SampleHisto from a region that matches a name of the sample
    * @param Region
    * @@param name
    * @return SampleHist
    */
    SampleHist* GetSampleHistFromName(const Region* const reg, const std::string& name) const;

    /**
     * A helper function to Copy a smoothing from a reference histogram to other histograms bin by bin
     * @param SampleHist
     * @param nominal histogram
     * @param up variation of histogram
     * @param down variation
     * @param flag if we want up or down variation returned
     * @return Up/Down variation for the new histogram
     */
    TH1* CopySmoothedHisto(const SampleHist* const sh, const TH1* const nominal, const TH1* const up, const TH1* const down, const bool isUp) const;

    /**
     * A helper function to get an index of systemaati variation that matches some name
     * @param SampleHist
     * @param name of the systematics
     * @return index
     */ 
    int GetSystIndex(const SampleHist* const sh, const std::string& name) const;
    
    // -------------------------

    std::string fName;
    std::string fDir;
    std::string fLabel;
    std::string fResultsFolder;
    std::string fInputFolder;
    std::string fInputName;
    std::string fFitResultsFile;

    std::vector < TFile* > fFiles;

    std::vector < Region* > fRegions;
    std::vector < Sample* > fSamples;
    std::vector < Systematic* > fSystematics;
    std::vector < NormFactor* > fNormFactors;
    std::vector < ShapeFactor* > fShapeFactors;
    std::vector < std::string > fSystematicNames;
    std::vector < std::string > fNormFactorNames;
    std::vector < std::string > fShapeFactorNames;

    int fNRegions;
    int fNSamples;
    int fNSyst;
    int fNNorm;
    int fNShape;
    std::string fPOI;
    std::string fPOIunit;
    bool fUseStatErr;
    float fStatErrThres;
    std::string fStatErrCons;
    bool fUseGammaPulls;

    float fLumi;
    float fLumiScale;
    float fLumiErr;

    float fThresholdSystPruning_Normalisation;
    float fThresholdSystPruning_Shape;
    float fThresholdSystLarge;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleNames;
    std::string fMCweight;
    std::string fSelection;

    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoNames;

    FitResults *fFitResults;

    bool fWithPullTables;

    int fIntCode_overall;
    int fIntCode_shape;

    int fInputType; // 0: histo, 1: ntup

    ConfigParser *fConfig;

    bool fSystControlPlots;
    bool fSystDataPlot_upFrame;
    bool fStatOnly;
    bool fGammasInStatOnly;
    bool fStatOnlyFit;
    bool fFixNPforStatOnlyFit;

    std::vector<std::string> fRegionsToPlot;
    std::vector<std::string> fSummaryPlotRegions;
    std::vector<std::string> fSummaryPlotLabels;
    std::vector<std::string> fSummaryPlotValidationRegions;
    std::vector<std::string> fSummaryPlotValidationLabels;

    float fYmin;
    float fYmax;
    float fRatioYmin;
    float fRatioYmax;
    float fRatioYminPostFit;
    float fRatioYmaxPostFit;
    std::string fRatioYtitle;
    std::string fRatioType;

    std::string fLumiLabel;
    std::string fCmeLabel;

    std::string fSuffix;
    std::string fSaveSuffix;

    bool fUpdate;
    bool fKeepPruning;

    float fBlindingThreshold;
    BlindingType fBlindingType;

    int fRankingMaxNP;
    std::string fRankingOnly;
    std::string fRankingPlot;
    std::string fImageFormat;
    std::string fAtlasLabel;

    bool fDoSummaryPlot;
    bool fDoMergedPlot;
    bool fDoTables;
    bool fDoSignalRegionsPlot;
    bool fDoPieChartPlot;

    std::string fGroupedImpactCategory;

    std::string fSummaryPrefix;

    //
    // Fit caracteristics
    //
    FitType fFitType;
    FitRegion fFitRegion;
    std::vector< std::string > fFitRegionsToFit;
    std::map< std::string, double > fFitNPValues;
    std::map< std::string, double > fFitFixedNPs;
    double fFitPOIAsimov;
    bool fFitIsBlind;
    bool fUseRnd;
    float fRndRange;
    long int fRndSeed;
    std::vector<std::string> fVarNameLH;
    float fLHscanMin;
    float fLHscanMax;
    int fLHscanSteps;
    std::vector<std::string> fVarNameMinos;
    std::vector<std::string> fVarNameHide;
    std::string fWorkspaceFileName;
    bool fDoGroupedSystImpactTable;
    std::map<std::string, std::string> fSubCategoryImpactMap;

    //
    // Limit parameters
    //
    LimitType fLimitType;
    bool fLimitIsBlind;
    double fLimitPOIAsimov;
    bool fSignalInjection;
    float fSignalInjectionValue;
    std::string fLimitParamName;
    float fLimitParamValue;
    std::string fLimitOutputPrefixName;
    float fLimitsConfidence;

    //
    // Significance parameters
    //
    bool fSignificanceIsBlind;
    double fSignificancePOIAsimov;
    std::string fSignificanceParamName;
    float fSignificanceParamValue;
    std::string fSignificanceOutputPrefixName;

    bool fCleanTables;
    bool bles;
    bool fSystCategoryTables;

    std::vector< std::string > fRegionGroups;

    bool fKeepPrefitBlindedBins;
    TH1D* fBlindedBins;

    std::string fCustomAsimov;

    std::string fTableOptions;

    bool fGetGoodnessOfFit;
    int fGetChi2;

    bool fTtresSmoothing;

    HistoTools::SmoothOption fSmoothOption;

    bool fSuppressNegativeBinWarnings;

    std::vector<std::string> fCustomFunctions;

    std::vector<std::string> fMorphParams;
    std::vector<std::pair<float,std::string> > fTemplatePair;
    std::vector<TRExFit::TemplateWeight> fTemplateWeightVec;
    TemplateInterpolationOption fTemplateInterpolationOption;

    std::string fBootstrap;
    int fBootstrapIdx;

    std::vector<std::string> fDecorrSysts;
    std::string fDecorrSuff;

    bool fDoNonProfileFit;
    float fNonProfileFitSystThreshold;
    int fFitToys;
    float fToysHistoMin;
    float fToysHistoMax;
    int fToysHistoNbins;
    std::string fSmoothMorphingTemplates;
    int fPOIPrecision;

    std::string fRankingPOIName;
    bool fUseATLASRounding;
    bool fUseATLASRoundingTxt;
    bool fUseATLASRoundingTex;
    bool fuseGammasForCorr;
    bool fPropagateSystsForMorphing;
    PrunningType fPrunningType;

    std::vector<int> fPrePostFitCanvasSize;
    std::vector<int> fSummaryCanvasSize;
    std::vector<int> fMergeCanvasSize;
    std::vector<int> fPieChartCanvasSize;
    std::vector<int> fNPRankingCanvasSize;
    
    std::vector<std::string> fBlindedParameters;
    
    float fLabelX;
    float fLabelY;
    float fLegendX1;
    float fLegendX2;
    float fLegendY;
    
    float fLabelXSummary;
    float fLabelYSummary;
    float fLegendX1Summary;
    float fLegendX2Summary;
    float fLegendYSummary;
    
    float fLabelXMerge;
    float fLabelYMerge;
    float fLegendX1Merge;
    float fLegendX2Merge;
    float fLegendYMerge;
    
    int fLegendNColumns;
    int fLegendNColumnsSummary;
    int fLegendNColumnsMerge;
    
    bool fShowRatioPad;
    bool fShowRatioPadSummary;
    bool fShowRatioPadMerge;
};

#endif
