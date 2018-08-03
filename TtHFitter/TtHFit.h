#ifndef TTHFIT_H
#define TTHFIT_H

#include "TtHFitter/Common.h"

#include "TGaxis.h"
#include "TPad.h"
#include "TPie.h"
#include "TF1.h"
#include "TRandom3.h"

#include "TtHFitter/TthPlot.h"
#include "TtHFitter/FitResults.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/NormFactor.h"
#include "TtHFitter/ShapeFactor.h"
#include "TtHFitter/ConfigParser.h"
#include "TtHFitter/HistoTools.h"
#include "TtHFitter/SampleHist.h"

class RooDataSet;
class RooWorkspace;
class Region;

class TtHFit {
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

    struct TemplateWeight{
        std::string function;
        std::string range;
        std::string name;
        float value;
    };

    TtHFit(std::string name="MyMeasurement");
    ~TtHFit();

    void SetPOI(std::string name="SigXsecOverSM");
    void SetStatErrorConfig(bool useIt=true, float thres=0.05, std::string cons="Poisson");
    void SetLumiErr(float err);
    void SetLumi(const float lumi);
    void SetFitType(FitType type);
    void SetLimitType( LimitType type );
    void SetFitRegion(FitRegion region);

    Sample* NewSample(std::string name,int type=0);
    Systematic* NewSystematic(std::string name);
    Region* NewRegion(std::string name);

    // ntuple stuff
    void AddNtuplePath(std::string path);
    void SetMCweight(std::string weight);
    void SetSelection(std::string selection);
    void SetNtupleName(std::string name);
    void SetNtupleFile(std::string name);
    void ComputeBining(int regIter);
    void defineVariable(int regIter);

    // histogram stuff
    void AddHistoPath(std::string path);
    void SetHistoName(std::string name);

    void SmoothSystematics(std::string syst="all");

    // create new root file with all the histograms
    void CreateRootFiles();
//     void CloseRootFiles();
    void WriteHistos();

    void DrawSystPlots();
    void DrawSystPlotsSumSamples();

    // read from ..
    void ReadNtuples();
    void ReadHistograms();
    void ReadHistos(/*std::string fileName=""*/);
    void CloseInputFiles();
    void CorrectHistograms();

    void DrawAndSaveAll(std::string opt="");

    // separation plots
    void DrawAndSaveSeparationPlots();

    TthPlot* DrawSummary(std::string opt="", TthPlot* = 0);
    void DrawMergedPlot(std::string opt="",std::string group="");
    void BuildYieldTable(std::string opt="",std::string group="");

    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols=0,int nRows=0);
    void DrawSignalRegionsPlot(int nRows,int nCols, std::vector < Region* > &regions);
    void DrawPieChartPlot(const std::string &opt="", int nCols=0,int nRows=0);
    void DrawPieChartPlot(const std::string &opt, int nCols,int nRows, std::vector < Region* > &regions);

    void CreateCustomAsimov();

    // turn to RooStat::HistFactory
    void ToRooStat(bool createWorkspace=true, bool exportOnly=true);

    void DrawPruningPlot();

    // fit etc...
    void Fit();
    RooDataSet* DumpData( RooWorkspace *ws, std::map < std::string, int > &regionDataType, std::map < std::string, double > &npValues, const double poiValue);
    std::map < std::string, double > PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType=SPLUSB, bool save=false, int debugLevel=1 );
    RooWorkspace* PerformWorkspaceCombination( std::vector < std::string > &regionsToFit );

    void PlotFittedNP();
    void PlotCorrelationMatrix();
    void GetLimit();
    void GetSignificance();
    void GetLikelihoodScan( RooWorkspace *ws, std::string varName, RooDataSet* data);

    // get fit results from txt file
    void ReadFitResults(const std::string& fileName);

    void Print();

    Region* GetRegion(std::string name);
    Sample* GetSample(std::string name);

    void ProduceNPRanking(std::string NPnames="all");
    void PlotNPRanking(bool flagSysts=true, bool flagGammas=true);
    void PlotNPRankingManager();

    void PrintSystTables(std::string opt="");

    void MergeSystematics(); // this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter

    // for template fitting
    void AddTemplateWeight(const std::string& name, float);
    std::vector<TemplateWeight> GetTemplateWeightVec(const TemplateInterpolationOption& opt);
    std::string GetWeightFunction(unsigned int itemp, const TemplateInterpolationOption& opt) const;

    /*
     * Function that returns string that represents smoothed abs value function
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSmoothLinearInterpolation(unsigned int itemp) const;

    /*
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

    /*
     * Helper function to approximate absolute value by sqrt(x^2+e)
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSquareRootLinearInterpolation(unsigned int itemp) const;

    /*
     * Helper function to apply correction to square root aproximation
     * @param will return value for a from -a*sqrt(x^2+epsilon) +b
     * @param will return value for b from -a*sqrt(x^2+epsilon) +b
     * @param central position of the function
     * @param left position of the function on x axis
     * @param epsilon = precision of the approximation
     */
    void GetSquareCorrection(double *a, double *b, float x_i, float x_left, float epsilon) const;

    void SmoothMorphTemplates(std::string name);
    bool MorphIsAlreadyPresent(const std::string& name, const float value) const;

    // for grouped impact evaluation
    void ProduceSystSubCategoryMap();
    void BuildGroupedImpactTable();

    /*
     * Helper function to calculate nominal scale for morphed samples
     * @param A pointer to SampleHist for which we need to calculate the scale factor
     * @return A scale factor
     */
    float GetNominalMorphScale(const SampleHist* const sh) const;

    /*
     * Helper function that runs toys experiments
     * @param A pointer to a workspace needed to run the fit
     */
    void RunToys(RooWorkspace* ws);
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
    std::string fNtupleFile;
    std::string fMCweight;
    std::string fSelection;
    std::string fNtupleName;

    std::vector<std::string> fHistoPaths;
    std::string fHistoName;
    std::string fHistoFile;

    FitResults *fFitResults;

    bool fWithPullTables;

    int fIntCode_overall;
    int fIntCode_shape;

    int fInputType; // 0: histo, 1: ntup

    ConfigParser *fConfig;

    bool fSystControlPlots;
    bool fSystDataPlot_upFrame;
    bool fStatOnly;
    bool fStatOnlyFit;
    bool fFixNPforStatOnlyFit;

    bool fRunROOTMacros;

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

    std::string fLumiLabel;
    std::string fCmeLabel;

    std::string fSuffix;
    std::string fSaveSuffix;

    bool fUpdate;
    bool fKeepPruning;

    float fBlindingThreshold;

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

    //
    // Significance parameters
    //
    bool fSignificanceIsBlind;
    double fSignificancePOIAsimov;

    bool fCleanTables;
    bool fSystCategoryTables;

    std::vector< std::string > fRegionGroups;

    bool fKeepPrefitBlindedBins;
    TH1F* fBlindedBins;

    std::string fCustomAsimov;

    int fRandomPOISeed;

    std::string fTableOptions;

    bool fGetGoodnessOfFit;
    int fGetChi2;

    bool fTtresSmoothing;

    HistoTools::SmoothOption fSmoothOption;

    bool fSuppressNegativeBinWarnings;

    std::vector<std::string> fCustomFunctions;

//     bool fRunMorphing;
    std::vector<std::string> fMorphParams;
    std::vector<std::pair<float,std::string> > fTemplatePair;
    std::vector<TtHFit::TemplateWeight> fTemplateWeightVec;
    TemplateInterpolationOption fTemplateInterpolationOption;

    std::string fBootstrap;
    int fBootstrapIdx;

    std::vector<std::string> fDecorrSysts;
    std::string fDecorrSuff;

    bool fDoNonProfileFit;
    int fFitToys;
    bool fSmoothMorphingTemplates;
    int fPOIPrecision;

    std::string fRankingPOIName;
};

#endif
