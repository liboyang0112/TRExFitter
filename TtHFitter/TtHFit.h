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
  void SetStatErrorConfig(bool useIt=true, float thres=0.05);
  void SetLumiErr(float err);
  Sample* NewSample(string name,int type=0);
  Systematic* NewSystematic(string name);
  Region* NewRegion(string name);
  
  // ntuple stuff
  void AddNtuplePath(string path);
  void SetMCweight(string weight);
  void SetSelection(string selection);
  void SetNtupleName(string name);

  // histogram stuff
  void AddHistoPath(string path);
  void SetHistoName(string name);

  // systematic handle
  void SmoothSystematics(string syst="all");
  
  // create new root file with all the histograms
  void WriteHistos(string fileName="",bool recreate=true);
  
  // read from ..
  void ReadNtuples();
  void ReadHistograms();
  void ReadHistos(string fileName="");
  void ReadAll(bool readNtuples=true,string fileName="");
  
  void DrawAndSaveAll(string opt="");

  TthPlot* DrawSummary(string opt="");
  
  void DrawSystPlots(string syst="all");
  
  // regions examples:
  // ...
  void DrawSignalRegionsPlot(int nRows,int nCols,Region *regions[MAXregions]=0x0);
  void DrawPieChartPlot();

  // turn to RooStat::HistFactory
  void ToRooStat(bool createWorkspace=true, bool exportOnly=true);

  // fit etc...
  void Fit();
  void PlotFittedNP();
  void GetLimit();
  
  // get fit results from txt file
  void ReadFitResults(string fileName);

  void Print();
  
  
  string fName;
  string fResultsFolder;
  
  Region *fRegions[MAXregions];
  Sample *fSamples[MAXsamples];
  Systematic *fSystematics[MAXsyst];
  int fNRegions;
  int fNSamples;
  int fNSyst;
  string fPOI;
  bool fUseStatErr;
  float fStatErrThres;
  
  float fLumiErr;
  
  vector<string> fNtuplePaths;
  string fMCweight;
  string fSelection;
  string fNtupleName;
  
  vector<string> fHistoPaths;
  string fHistoName;

  FitResults *fFitResults;
};

#endif
