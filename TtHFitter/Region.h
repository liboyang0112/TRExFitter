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
  
  SampleHist* SetDataHist(Sample *sample, string histoName, string fileName);
  SampleHist* SetDataHist(Sample *sample, TH1* hist );
  SampleHist* SetSigHist(Sample *sample, string histoName, string fileName);
  SampleHist* SetSigHist(Sample *sample, TH1* hist );
  SampleHist* AddBkgHist(Sample *sample, string histoName, string fileName);
  SampleHist* AddBkgHist(Sample *sample, TH1* hist);  
  void AddSystematic(Systematic *syst);
  SampleHist* GetSampleHist(string sampleName);  
  void BuildPreFitErrorHist();
  TCanvas* DrawPreFit();
  void BuildPostFitErrorHist(FitResults *fitRes);
  TCanvas* DrawPostFit(FitResults *fitRes);
  
  void AddSample(Sample *sample);
  void AddSelection(string selection);
  void SetVariable(string variable,int nbin,float xmin,float xmax);
  void SetAllSamples(bool readNtuples=true,string fileName="MyMeasurement_histos.root");
  
  void SetVariableTitle(string name);
  void SetLabel(string label);
  
  string fName;
  string fVariableTitle;
  string fLabel;
  string fFitName;
  bool fHasData;
  SampleHist *fData;
  bool fHasSig;
  SampleHist *fSig;
  int fNBkg;
  SampleHist *fBkg[MAXsamples];
  int fNSamples;
  Sample* fSamples[MAXsamples];
  
  // to draw
  THStack *fStack;
  TH1* fTot;
  TGraphAsymmErrors *fErr;
  
  // ntuple stuff
  string fVariable;
  int fNbins;
  float fXmin, fXmax;
  string fSelection;
  string fMCweight;
  vector<string> fNtuplePaths;
  vector<string> fNtuplePathSuffs;
  vector<string> fNtupleNames;
  vector<string> fNtupleNameSuffs;

  // systematics & norm.factors
  int fNSyst;
  Systematic* fSystematics[MAXsyst];
  int fNNorm;
  NormFactor* fNormFactors[MAXnorm];  
};

float GetDeltaN(float alpha, float Iz, float Ip, float Imi);

#endif
