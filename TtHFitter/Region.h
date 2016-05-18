#ifndef __Region__
#define __Region__

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TBox.h"

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

#include "TtHFitter/Common.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/SampleHist.h"
#include "TtHFitter/FitResults.h"

#include "TtHFitter/TthPlot.h"
#include "TtHFitter/TtHFit.h"


class Sample;
class Systematic;
class SampleHist;

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
    
    Region(string name);
    ~Region();

    // -------
    // Methods
    // -------
    
    SampleHist* SetSampleHist(Sample *sample, string histoName, string fileName);
    SampleHist* SetSampleHist(Sample *sample, TH1* hist );
    SampleHist* GetSampleHist(string &sampleName);

    void BuildPreFitErrorHist();
    TthPlot* DrawPreFit(string opt="");
    void BuildPostFitErrorHist(FitResults *fitRes);
    TthPlot* DrawPostFit(FitResults *fitRes,string opt="");
      
    void SetBinning(int N, double *bins);
    void Rebin(int N);
    void SetRegionType(RegionType type);
    void SetRegionDataType( DataType type );
    void AddSample(Sample *sample);
    
    void AddSelection(string selection);
    void AddMCweight(string weight);
    void SetVariable(string variable,int nbin,float xmin,float xmax,string corrVar1="",string corrVar2="");

    void SetHistoName(string name); // name of the histogram to read (the same for each sample)
    void AddSystematic(Systematic *syst);

    // cosmetics
    void SetVariableTitle(string name);
    void SetLabel(string label,string shortLabel="");
    
    // log
    void Print();
    
    void PrintSystTable(string opt="");
    
    // -------
    // Members
    // -------
    
    string fName;
    string fVariableTitle;
    string fLabel; // something like "e/µ + 6 j, >=4 b b"
    string fShortLabel; // something like "6j,3b"
    string fFitName;
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
    string fBinTransfo;
    double fTransfoDzBkg;
    double fTransfoDzSig;
    double fTransfoFzBkg;
    double fTransfoFzSig;
    double fTransfoJpar1;
    double fTransfoJpar2;
    double fTransfoJpar3;
    string fVariable;
    string fCorrVar1;
    string fCorrVar2;
    int fNbins;
    float fXmin, fXmax;
    string fSelection;
    string fMCweight;
    std::vector<string> fNtuplePaths;
    std::vector<string> fNtuplePathSuffs;
    std::vector<string> fNtupleFiles;
    std::vector<string> fNtupleFileSuffs;
    std::vector<string> fNtupleNames;
    std::vector<string> fNtupleNameSuffs;

    // histogram stuff
//     string fHistoName;
    double *fHistoBins;
    int fHistoNBinsRebin;
    std::vector<string> fHistoPaths;
    std::vector<string> fHistoPathSuffs;
    std::vector<string> fHistoFiles;
    std::vector<string> fHistoFileSuffs;
    std::vector<string> fHistoNames;
    std::vector<string> fHistoNameSuffs;
    
    int fNSyst;
    std::vector < Systematic* > fSystematics;
    int fNNorm;
    std::vector < NormFactor* > fNormFactors;
    
    // plot objects
    TthPlot *fPlotPreFit;
    TthPlot *fPlotPostFit;
    
    bool fUseStatErr;
    
    int fIntCode_overall;
    int fIntCode_shape;
    
    std::vector< string > fSystNames;
    std::vector< string > fNpNames;
    
    TtHFit::FitType fFitType;
    string fPOI;
    string fFitLabel;
    
    string fLumiLabel;
    string fCmeLabel;
    
    float fLumiScale;
    
    bool fLogScale;
    
    float fBinWidth;
    
    float fBlindingThreshold;
    
    bool fSkipSmoothing;
};


// Functions

// for post-fit plots
float GetDeltaN(float alpha, float Iz, float Ip, float Imi, int intCode=4);

// To build the total error band
TGraphAsymmErrors* BuildTotError( TH1* h_nominal, std::vector< TH1* > h_up, std::vector< TH1* > h_down, std::vector< string > systNames, CorrelationMatrix *matrix=0x0 );

#endif
