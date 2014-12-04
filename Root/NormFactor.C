#include "TtHFitter/NormFactor.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

NormFactor::NormFactor(){}
NormFactor::NormFactor(string name, float nominal, float min, float max){
  Set(name,nominal,min,max);
}
NormFactor::~NormFactor(){}
void NormFactor::Set(string name, float nominal, float min, float max){
  fName = name;
  fNominal = nominal;
  fMin = min;
  fMax = max;
}

void NormFactor::Print(){
  cout << "        NormFactor: " << fName << "\t" << fNominal << ", " << fMin << ", " << fMax << endl;
}