#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

#include "TtHFitter/Common.h"

#include "TtHFitter/SystematicHist.h"
#include "TtHFitter/Sample.h"

#ifndef __SampleHist__
#define __SampleHist__

class SampleHist {
public:
  SampleHist();
  SampleHist(Sample *sample,TH1 *hist);
  SampleHist(Sample *sample, string histoName, string fileName);
  ~SampleHist();

  TH1* GetHist();
  Sample* GetSample();
  SystematicHist* AddOverallSyst(string name,float up,float down);
  SystematicHist* AddStatSyst(string name,int i_bin);
  SystematicHist* AddHistoSyst(string name,TH1* h_up,TH1* h_down);
  SystematicHist* AddHistoSyst(string name,string histoName_up, string fileName_up,string histoName_down, string fileName_down, int pruned=0);
  SystematicHist* GetSystematic(string systName);
  SystematicHist* GetSystFromNP(string NuisParName);
  NormFactor* AddNormFactor(string name,float nominal, float min, float max);
  NormFactor* AddNormFactor(NormFactor *normFactor);
  NormFactor* GetNormFactor(string name);
  ShapeFactor* AddShapeFactor(string name,float nominal, float min, float max);
  ShapeFactor* AddShapeFactor(ShapeFactor *shapeFactor);
  ShapeFactor* GetShapeFactor(string name);

  bool HasSyst(string name);
  bool HasNorm(string name);
  bool HasShapeFactor(string name);

  void WriteToFile(TFile *f=0x0);
  void ReadFromFile();

  void FixEmptyBins(const bool suppress);

  void Print();

  void Rebin(int ngroup = 2, const Double_t* xbins = 0);
  void DrawSystPlot( const string &syst="all", TH1* h_data=0x0, bool SumAndData=false, bool bothPanels=false );
  void SmoothSyst(string syst="all",bool force=false, bool TtresSmoothing = false);

  void Divide(  SampleHist* sh);
  void Multiply(SampleHist* sh);
  void Add(     SampleHist* sh,float scale=1.);

  void SampleHistAdd(SampleHist* h);
  void CloneSampleHist(SampleHist* h, std::set<std::string> names);

//   void SmoothNominal(); // new
//   void Rebin(); // new

  string fName;
  Sample *fSample;
  TH1 *fHist;
  TH1 *fHist_orig;  // new
  TH1 *fHist_regBin;  // new
  TH1 *fHist_postFit;
  string fFileName;
  string fHistoName;
  bool fIsData;
  bool fIsSig;

  int fNSyst;
  std::vector < SystematicHist* > fSyst;
  int fNNorm;
  std::vector < NormFactor* > fNormFactors;
  int fNShape;
  std::vector < ShapeFactor* > fShapeFactors;

  // other useful info
  string fFitName;
  string fRegionName;
  string fRegionLabel;
  string fVariableTitle;
  bool fSystSmoothed;
};

#endif

