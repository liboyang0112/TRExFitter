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

#ifndef __TtHFit__
#define __TtHFit__

class Region;
class Sample;
class Systematic;

class TtHFit {
public:
  TtHFit(string name="MyMeasurement");
  ~TtHFit();

  void SetPOI(string name="SigXsecOverSM");
  void SetLumiErr(float err);
  Sample* NewSample(string name);
  Systematic* NewSystematic(string name);
  Region* NewRegion(string name);
  
  // ntuple stuff
  void AddNtuplePath(string path);
  void SetMCweight(string weight);
  void SetSelection(string selection);
  void SetTreeName(string name);

  // create new root file with all the histograms
  void WriteHistos(string fileName="MyMeasurement_histos.root",bool recreate=true);
  
  // read from ..
  void ReadAll(bool readNtuples=true,string fileName="MyMeasurement_histos.root");
  
  void DrawAndSaveAll();

  // turn to RooStat::HistFactory
  void ToRooStat(bool createWorkspace=true, bool exportOnly=true);
  
  string fName;
  string fResultsFolder;
  
  Region *fRegions[MAXregions];
  Sample *fSamples[MAXsamples];
  Systematic *fSystematics[MAXsyst];
  int fNRegions;
  int fNSamples;
  int fNSyst;
  string fPOI;
  
  float fLumiErr;
  
  vector<string> fNtuplePaths;
  string fMCweight;
  string fSelection;
  string fTreeName;
  
  FitResults *fFitResults;
};

#endif
