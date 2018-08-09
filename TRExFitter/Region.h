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
    SampleHist* GetSampleHist(std::string &sampleName);

    void BuildPreFitErrorHist();
    TRExPlot* DrawPreFit(std::string opt="");
    double GetMultFactors( FitResults *fitRes, 
                                std::ofstream& pullTex,
                                const int i /*sample*/, const int i_bin /*bin number*/,
                                const double binContent0,
                                const std::string &syst = "",
                                const bool isUp = true);

    void BuildPostFitErrorHist(FitResults *fitRes, const std::vector<std::string>& morph_names);
    TRExPlot* DrawPostFit(FitResults *fitRes,std::ofstream & pullTex, const std::vector<std::string>& morph_names, std::string opt="");

    void SetBinning(int N, double *bins);
    void Rebin(int N);
    void SetRegionType(RegionType type);
    void SetRegionDataType( DataType type );
    void AddSample(Sample *sample);

    void AddSelection(std::string selection);
    void AddMCweight(std::string weight);
    void SetVariable(std::string variable,int nbin,float xmin,float xmax,std::string corrVar1="",std::string corrVar2="");
    void SetAlternativeVariable(std::string variable,std::string sample);
    bool UseAlternativeVariable(std::string sample);
    std::string GetAlternativeVariable(std::string sample);
    void SetAlternativeSelection(std::string selection,std::string sample);
    bool UseAlternativeSelection(std::string sample);
    std::string GetAlternativeSelection(std::string sample);

    void SetHistoName(std::string name); // name of the histogram to read (the same for each sample)
    void AddSystematic(Systematic *syst);

    // cosmetics
    void SetVariableTitle(std::string name);
    void SetLabel(std::string label,std::string shortLabel="");

    // log
    void Print();

    void PrintSystTable(FitResults* fitRes,std::string opt="");

    /**
      * Helper function to get postfit scales of "normalization" parameters used for morphing
      * @param pointer to FitResults class that stores the fit output
      * @param dummy parameter that will be filled
      * @param dummy parameter that will be filled
      */
    void PrepareMorphScales(FitResults *fitRes, std::vector<double> *morph_scale, std::vector<double> *morph_scale_nominal);
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
    float fXmin, fXmax;
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

    TH1F* fBlindedBins;
    bool fKeepPrefitBlindedBins;
    int fGetChi2;

    std::vector<int> fDropBins;

    std::vector<std::string> fBinLabels;

    float fChi2val;
    int fNDF;
    float fChi2prob;

    bool fUseGammaPulls;
};


// Functions

// for post-fit plots
float GetDeltaN(float alpha, float Iz, float Ip, float Imi, int intCode=4);
std::map < int , double > GetDeltaNForUncertainties(float alpha, float alpha_errUp, float alpha_errDown, float Iz, float Ip, float Imi, int intCode);


// To build the total error band
TGraphAsymmErrors* BuildTotError( TH1* h_nominal, std::vector< TH1* > h_up, std::vector< TH1* > h_down, std::vector< std::string > systNames, CorrelationMatrix *matrix=0x0 );

std::pair<double,int> GetChi2Test( TH1* h_data, TH1* h_nominal, std::vector< TH1* > h_up, std::vector< std::string > fSystNames, CorrelationMatrix *matrix=0x0 );

#endif