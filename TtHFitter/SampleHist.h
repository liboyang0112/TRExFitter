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


#ifndef __SampleHist__
#define __SampleHist__

class SampleHist {
public:
  SampleHist(Sample *sample,TH1 *hist, bool isData=false, bool isSig=false);
  SampleHist(Sample *sample, string histoName, string fileName, bool isData=false, bool isSig=false);
  ~SampleHist();
  
  TH1* GetHist();
  Sample* GetSample();
  SystematicHisto* AddOverallSyst(string name,float up,float down);
  SystematicHisto* AddHistoSyst(string name,TH1* h_up,TH1* h_down);
  SystematicHisto* AddHistoSyst(string name,string histoName_up, string fileName_up,string histoName_down, string fileName_down);
  SystematicHisto* GetSystematic(string systName);
  NormFactor* AddNormFactor(string name,float nominal, float min, float max);
  NormFactor* AddNormFactor(NormFactor *normFactor);
  NormFactor* GetNormFactor(string name);

  bool HasSyst(string name);
  
  void WriteToFile();
  void ReadFromFile();
  
  string fName;
  Sample *fSample;
  TH1 *fHist;
  string fFileName;
  string fHistoName;
  bool fIsData;
  bool fIsSig;

  int fNSyst;
  SystematicHisto* fSyst[MAXsyst];
  int fNNorm;
  NormFactor* fNormFactors[MAXnorm];
};

#endif

