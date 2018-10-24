#ifndef REGION_H
#define REGION_H

/// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/TRExFit.h"

/// c++ includes
#include <map>
#include <string>
#include <vector>

/// Forwards class declaration
class FitResults;
class Sample;
class Systematic;
class TH1;
class THStack;
class TGraphAsymmErrors;
class TRExPlot;
class TRExFit;
class SampleHist;
class ShapeFactor;
class CorrelationMatrix;
class TFile;

class Region {
public:

    enum RegionType {
        CONTROL = 1,
        VALIDATION = 2,
        SIGNAL = 3
    };

    enum DataType {
        REALDATA = 0,
        ASIMOVDATA = 1
    };

    Region(std::string name);
    ~Region();

    // -------
    // Methods
    // -------

    SampleHist* SetSampleHist(Sample *sample, std::string histoName, std::string fileName);
    SampleHist* SetSampleHist(Sample *sample, TH1* hist );
    SampleHist* GetSampleHist(const std::string &sampleName) const;

    void BuildPreFitErrorHist();
    TRExPlot* DrawPreFit(const std::vector<int>& canvasSize, std::string opt="");
    double GetMultFactors( FitResults *fitRes, 
                                std::ofstream& pullTex,
                                const int i /*sample*/, const int i_bin /*bin number*/,
                                const double binContent0,
                                const std::string &syst = "",
                                const bool isUp = true) const;

    void BuildPostFitErrorHist(FitResults *fitRes, const std::vector<std::string>& morph_names);
    TRExPlot* DrawPostFit(FitResults *fitRes,std::ofstream & pullTex, const std::vector<std::string>& morph_names, const std::vector<int>& canvasSize, std::string opt="");

    void SetBinning(int N, double *bins);
    void Rebin(int N);
    void SetRebinning(int N, double *bins);
    void SetRegionType(RegionType type);
    void SetRegionDataType( DataType type );
    void AddSample(Sample *sample);

    void AddSelection(const std::string& selection);
    void AddMCweight(const std::string& weight);
    void SetVariable(const std::string& variable,int nbin,float xmin,float xmax,std::string corrVar1="",std::string corrVar2="");
    void SetAlternativeVariable(const std::string& variable, const std::string& sample);
    bool UseAlternativeVariable(const std::string& sample);
    std::string GetAlternativeVariable(const std::string& sample) const;
    void SetAlternativeSelection(const std::string& selection, const std::string& sample);
    bool UseAlternativeSelection(const std::string& sample);
    std::string GetAlternativeSelection(const std::string& sample) const;

    void AddSystematic(Systematic *syst);

    // cosmetics
    void SetVariableTitle(const std::string& name);
    void SetLabel(const std::string& label,std::string shortLabel="");

    // log
    void Print() const;

    void PrintSystTable(FitResults* fitRes,std::string opt="") const;

    /**
      * Helper function to get postfit scales of "normalization" parameters used for morphing
      * @param pointer to FitResults class that stores the fit output
      * @param dummy parameter that will be filled
      * @param dummy parameter that will be filled
      */
    void PrepareMorphScales(FitResults *fitRes, std::vector<double> *morph_scale, std::vector<double> *morph_scale_nominal) const;
    // -------
    // Members
    // -------

    std::string fName;
    std::string fVariableTitle;
    std::string fYTitle;
    std::string fLabel; // something like "e/mu + 6 j, >=4 b b"
    std::string fShortLabel; // something like "6j,3b"
    std::string fTexLabel;
    std::string fFitName;
    RegionType fRegionType;
    DataType fRegionDataType;
    bool fHasData;
    SampleHist *fData;
    bool fHasSig;
    int fNSig;
    SampleHist *fSig[MAXsamples];
    int fNBkg;
    SampleHist *fBkg[MAXsamples];
    int fNSamples;
    std::vector < SampleHist* > fSampleHists;
    std::vector < Sample* > fSamples;
    float fYmaxScale;
    float fYmin;
    float fYmax;
    float fRatioYmin;
    float fRatioYmax;
    float fRatioYminPostFit;
    float fRatioYmaxPostFit;

    // to draw
    THStack *fStack;
    TH1* fTot;
    TGraphAsymmErrors *fErr;
    TH1* fTotUp[MAXsyst];
    TH1* fTotDown[MAXsyst];

    // post fit
    THStack *fStack_postFit;
    TH1* fTot_postFit;
    TGraphAsymmErrors *fErr_postFit;
    TH1* fTotUp_postFit[MAXsyst];
    TH1* fTotDown_postFit[MAXsyst];

    // ntuple stuff
    std::string fBinTransfo;
    double fTransfoDzBkg;
    double fTransfoDzSig;
    double fTransfoFzBkg;
    double fTransfoFzSig;
    double fTransfoJpar1;
    double fTransfoJpar2;
    double fTransfoJpar3;
    std::vector<std::string> fAutoBinBkgsInSig;
    std::string fVariable;
    std::map<std::string, std::string> fAlternativeVariables;
    std::map<std::string, std::string> fAlternativeSelections;
    std::string fCorrVar1;
    std::string fCorrVar2;
    int fNbins;
    float fXmin;
    float fXmax;
    std::string fSelection;
    std::string fMCweight;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtuplePathSuffs;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleFileSuffs;
    std::vector<std::string> fNtupleNames;
    std::vector<std::string> fNtupleNameSuffs;

    // histogram stuff
    double *fHistoBins;
    int fHistoNBinsRebin;
    double *fHistoBinsPost;
    int fHistoNBinsRebinPost;
    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoPathSuffs;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoFileSuffs;
    std::vector<std::string> fHistoNames;
    std::vector<std::string> fHistoNameSuffs;

    int fNSyst;
    std::vector < Systematic* > fSystematics;
    int fNNorm;
    std::vector < NormFactor* >  fNormFactors;
    int fNShape;
    std::vector < ShapeFactor* > fShapeFactors;

    // plot objects
    TRExPlot *fPlotPreFit;
    TRExPlot *fPlotPostFit;

    bool fUseStatErr;

    int fIntCode_overall;
    int fIntCode_shape;

    std::vector< std::string > fSystNames;
    std::vector< std::string > fNpNames;

    TRExFit::FitType fFitType;
    std::string fPOI;
    std::string fFitLabel;

    std::string fLumiLabel;
    std::string fCmeLabel;

    float fLumiScale;

    bool fLogScale;

    float fBinWidth;

    float fBlindingThreshold;

    bool fSkipSmoothing;

    std::string fATLASlabel;
    std::string fSuffix;

    std::string fGroup; // used to split yield tables

    TH1D* fBlindedBins;
    bool fKeepPrefitBlindedBins;
    int fGetChi2;

    std::vector<int> fDropBins;

    std::vector<std::string> fBinLabels;

    float fChi2val;
    int fNDF;
    float fChi2prob;

    bool fUseGammaPulls;

    std::vector<float> fXaxisRange;
};


// Functions

// for post-fit plots
double GetDeltaN(double alpha, double Iz, double Ip, double Imi, int intCode=4);
std::map < int , double > GetDeltaNForUncertainties(double alpha, double alpha_errUp, double alpha_errDown, double Iz, double Ip, double Imi, int intCode);

// To build the total error band
TGraphAsymmErrors* BuildTotError( TH1* h_nominal, std::vector< TH1* > h_up, std::vector< TH1* > h_down, std::vector< std::string > systNames, CorrelationMatrix *matrix=0x0 );

std::pair<double,int> GetChi2Test( TH1* h_data, TH1* h_nominal, std::vector< TH1* > h_up, std::vector< std::string > fSystNames, CorrelationMatrix *matrix=0x0 );

#endif
