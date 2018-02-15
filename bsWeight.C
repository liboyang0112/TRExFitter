#include "TRandom3.h"

float bsWeight(int index=0,int eventNumber=0,int runNumber=0){
  gRandom->SetSeed(index+eventNumber*1e3+runNumber*1e6);
  return gRandom->Poisson(1);
};
