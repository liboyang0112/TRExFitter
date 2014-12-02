#include <string>

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

#include "TtHFitter/Common.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/NormFactor.h"


#ifndef __Sample__
#define __Sample__

class Sample {
public:
  Sample(string name);
  ~Sample();
  
  void SetIsData(bool isData=true);
  void SetIsSignal(bool isSignal=true);
  void SetTitle(string title);
  void SetFillColor(int color);
  void SetLineColor(int color);
  void AddNormFactor(string name,float nominal,float min,float max);
  void AddNormFactor(NormFactor *factor);
  
  void AddNtupleName(string name);
  void AddNtuplePath(string path);
  void AddSystematic(Systematic *syst);
  void SetMCweight(string weight);
  
  string fName;
  string fFitName;
  string fTitle;
  bool fIsData;
  bool fIsSignal;
  int fFillColor;
  int fLineColor;
  
  // to read from ntuples
  string fSelection;
  string fMCweight;
  string fTreeName;
  vector<string> fNtuplePaths;
  vector<string> fNtuplePathSuffs;
  vector<string> fNtupleNames;
  vector<string> fNtupleNameSuffs;
  
  // systematics & norm.factors
  int fNSyst;
  Systematic* fSystematics[MAXsyst];
  int fNNorm;
  NormFactor* fNormFactors[MAXnorm];
  
//   // to read from histograms
//   vector<string> fHistoPaths;
//   vector<string> fHistoFiles;
//   vector<string> fHistoNames;
};

#endif

