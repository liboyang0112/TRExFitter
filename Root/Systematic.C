#include "TtHFitter/Systematic.h"

// -------------------------------------------------------------------------------------------------
// Systematic

Systematic::Systematic(string name){
  fName = name;
  //
  fWeightUp = "";
  fWeightSufUp = "";  
  fNtuplePathsUp.clear();
  fNtuplePathSufUp = "";
  fNtupleFilesUp.clear();
  fNtupleFileSufUp = "";
  fTreeNameUp = "";
  fTreeNameSufUp = "";
  //
  fWeightDown = "";
  fWeightSufDown = "";  
  fNtuplePathsDown.clear();
  fNtuplePathSufDown = "";
  fNtupleFilesDown.clear();
  fNtupleFileSufDown = "";
  fTreeNameDown = "";
  fTreeNameSufDown = "";

}
Systematic::~Systematic(){}

bool Systematic::IsOverallOnly(){
  bool res = true;
  if(fWeightUp!="") res = false;
  if(fWeightSufUp!="") res = false;
  if(fNtuplePathsUp.size()!=0) res = false;
  if(fNtuplePathSufUp!="") res = false;
  if(fNtupleFilesUp.size()!=0) res = false;
  if(fNtupleFileSufUp!="") res = false;
  if(fTreeNameUp!="") res = false;
  if(fTreeNameSufUp!="") res = false;
  //
  if(fWeightDown!="") res = false;
  if(fWeightSufDown!="") res = false;
  if(fNtuplePathsDown.size()!=0) res = false;
  if(fNtuplePathSufDown!="") res = false;
  if(fNtupleFilesDown.size()!=0) res = false;
  if(fNtupleFileSufDown!="") res = false;
  if(fTreeNameDown!="") res = false;
  if(fTreeNameSufDown!="") res = false;
  return res;
}
