#include "TtHFitter/FitResults.h"

FitResults::FitResults(){
}
FitResults::~FitResults(){
}

void FitResults::AddNuisPar(NuisParameter *par){
  fNuisPar[(int)fNuisParNames.size()] = par;
  string p = par->fName;
  fNuisParIdx[p] = (int)fNuisParNames.size();
  fNuisParNames.push_back(p);
  fNuisParIsThere[p] = true;
}

float FitResults::GetNuisParValue(string p){
  int idx = -1;
  if(fNuisParIsThere[p]){
    idx = fNuisParIdx[p];
  }
  else{
    cout << "  WARNING: NP " << p << " not found... Returning 0." << endl;
    return 0.;
  }
  return fNuisPar[idx]->fFitValue;
}

float FitResults::GetNuisParErrUp(string p){
  int idx = -1;
  if(fNuisParIsThere[p]){
    idx = fNuisParIdx[p];
  }
  else{
    cout << "  WARNING: NP " << p << " not found... Returning error = 1." << endl;
    return 1.;
  }
  return fNuisPar[idx]->fPostFitUp;
}

float FitResults::GetNuisParErrDown(string p){
  int idx = -1;
  if(fNuisParIsThere[p]){
    idx = fNuisParIdx[p];
  }
  else{
    cout << "  WARNING: NP " << p << " not found... Returning error = 1." << endl;
    return 1.;
  }
  return fNuisPar[idx]->fPostFitDown;
}
