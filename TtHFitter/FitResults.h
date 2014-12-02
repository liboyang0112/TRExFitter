#include "TtHFitter/NuisParameter.h"
#include "TtHFitter/CorrelationMatrix.h"

#ifndef __FitResults__
#define __FitResults__

class FitResults {
public:
  FitResults();
  ~FitResults();
  
  // ReadFromROOT(string rootFileName);
  // ReadFromTXT(string txtFileName);
  
  void AddNuisPar(NuisParameter *par);
  float GetNuisParValue(string p);
  float GetNuisParErrUp(string p);
  float GetNuisParErrDown(string p);
  
  vector<string> fNuisParNames;
  map<string,int> fNuisParIdx;
  map<string,bool> fNuisParIsThere;
  NuisParameter *fNuisPar[MAXsyst];
  CorrelationMatrix *fCorrMatrix;

};

#endif
