#include <iostream>

#include "TChain.h"
#include "TH1F.h"

#ifndef __Common__
#define __Common__

using namespace std;

const int MAXregions = 10;
const int MAXsamples = 10;
const int MAXsyst = 100;
const int MAXnorm = 3;

TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight);
TH1* HistFromFile(string fileName,string histoName);
void WriteHistToFile(TH1* h,string fileName,string option="UPDATE");

#endif
