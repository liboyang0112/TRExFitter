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

namespace TtHFitter{
    extern int DEBUGLEVEL;
    void SetDebugLevel(int level=0);
    extern bool SHOWYIELDS; // flag to show or not yields in plots
    extern bool SYSTCONTROLPLOTS;
    extern float CORRELATIONTHRESHOLD;
};

const int MAXregions = 100;
const int MAXsamples = 100;
const int MAXsyst = 200;
const int MAXnorm = 3;

TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight);
TH1* HistFromFile(string fullName);
TH1* HistFromFile(string fileName,string histoName);
void WriteHistToFile(TH1* h,string fileName,string option="UPDATE");
std::vector<string> CreatePathsList( std::vector<string> paths, std::vector<string> pathSufs,
                                    std::vector<string> files, std::vector<string> fileSufs,
                                    std::vector<string> names, std::vector<string> nameSufs);
std::vector<string> CombinePathSufs( std::vector<string> pathSufs, std::vector<string> newPathSufs );
std::vector<string> ToVec(string s);
// string RemovePrefix(string s,string prefix);
string ReplaceString(string subject, const string& search,
                     const string& replace);

int FindInStringVector(std::vector<string> v, string s);

#endif
