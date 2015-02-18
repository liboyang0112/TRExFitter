#include <string>

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TMath.h"

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

#include "TtHFitter/Common.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/SampleHist.h"
#include "TtHFitter/FitResults.h"

#include "TtHFitter/TthPlot.h"

#ifndef __Region__
#define __Region__

class Sample;
class Systematic;
class SampleHist;

class Region {
public:
  Region(string name);
  ~Region();

  // -------
  // Methods
  // -------
  
  SampleHist* SetSampleHist(Sample *sample, string histoName, string fileName);
  SampleHist* SetSampleHist(Sample *sample, TH1* hist );
  SampleHist* GetSampleHist(string sampleName);  

  void BuildPreFitErrorHist();
  TthPlot* DrawPreFit(string opt="");
  void BuildPostFitErrorHist(FitResults *fitRes);
  TthPlot* DrawPostFit(FitResults *fitRes,string opt="");
  
  void AddSample(Sample *sample);
  void AddSelection(string selection);
  void AddMCweight(string weight);
  void SetVariable(string variable,int nbin,float xmin,float xmax);
  void SetAllSamples(bool readNtuples=true,string fileName="MyMeasurement_histos.root");
  
  void AddSystematic(Systematic *syst);

  // cosmetics
  void SetVariableTitle(string name);
  void SetLabel(string label,string shortLabel="");
  
  // log
  void Print();
  
  // -------
  // Members
  // -------
  
  string fName;
  string fVariableTitle;
  string fLabel; // something like "e/µ + 6 j, >=4 b b"
  string fShortLabel; // something like "6j,3b"
  string fFitName;
  bool fHasData;
  SampleHist *fData;
  bool fHasSig;
  SampleHist *fSig;
  int fNBkg;
  SampleHist *fBkg[MAXsamples];
  int fNSamples;
  SampleHist* fSampleHists[MAXsamples];
  Sample* fSamples[MAXsamples];
  
  // to draw
  THStack *fStack;
  TH1* fTot;
  TGraphAsymmErrors *fErr;

  // post fit
  THStack *fStack_postFit;
  TH1* fTot_postFit;
  TGraphAsymmErrors *fErr_postFit;
  
  // ntuple stuff
  string fVariable;
  int fNbins;
  float fXmin, fXmax;
  string fSelection;
  string fMCweight;
  vector<string> fNtuplePaths;
  vector<string> fNtuplePathSuffs;
  vector<string> fNtupleFiles;
  vector<string> fNtupleFileSuffs;
  vector<string> fNtupleNames;
  vector<string> fNtupleNameSuffs;

  // systematics & norm.factors
  int fNSyst;
  Systematic* fSystematics[MAXsyst];
  int fNNorm;
  NormFactor* fNormFactors[MAXnorm]; 
  
  // plot objects
  TthPlot *fPlotPreFit;
  TthPlot *fPlotPostFit;

  bool fUseStatErr;
};

// for post-fit plots
float GetDeltaN(float alpha, float Iz, float Ip, float Imi);

#endif
