// c++ stuff
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

// ROOT stuff
#include "TArrow.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TNamed.h"
#include "TObject.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

// RooStats stuff
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

// ATLAS stuff
#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "AtlasUtils.h"

using namespace std;

#ifndef __Common__
#define __Common__

struct SampleType{
  enum {
    Background, // 0
    Signal, // 1
    Data // 2
  };
};

struct SystType{
  enum {
    Overall, // 0
    Shape, // 1
    Histo // 2
  };
};

namespace TtHFitter{
  extern int DEBUGLEVEL;
  void SetDebugLevel(int level=0);
};

const int MAXregions = 10;
const int MAXsamples = 10;
const int MAXsyst = 150;
const int MAXnorm = 3;

TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight);
TH1* HistFromFile(string fullName);
TH1* HistFromFile(string fileName,string histoName);
void WriteHistToFile(TH1* h,string fileName,string option="UPDATE");
vector<string> CreatePathsList( vector<string> paths, vector<string> pathSufs, 
                                vector<string> files, vector<string> fileSufs, 
                                vector<string> names, vector<string> nameSufs);
vector<string> CombinePathSufs( vector<string> pathSufs, vector<string> newPathSufs );
vector<string> ToVec(string s);
// string RemovePrefix(string s,string prefix);
std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace);


#endif
