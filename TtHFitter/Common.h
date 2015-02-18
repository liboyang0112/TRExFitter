#include <iostream>
#include <sstream>
#include <fstream>

#include "TChain.h"
#include "TH1F.h"

#ifndef __Common__
#define __Common__

using namespace std;

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
  int DEBUGLEVEL;
  void SetDebugLevel(int level=0);
};

const int MAXregions = 10;
const int MAXsamples = 10;
const int MAXsyst = 100;
const int MAXnorm = 3;

TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight);
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
void SmoothSystHistos(TH1* h_nominal,TH1* h_syst_up,TH1* h_syst_down); // h_syst_up/down will be overwritten (!!)
float Separation(TH1* h1,TH1* h2);

#endif
