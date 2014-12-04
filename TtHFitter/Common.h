#include <iostream>

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

#endif
