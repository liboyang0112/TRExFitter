#include "TtHFitter/Systematic.h"

// -------------------------------------------------------------------------------------------------
// Systematic

Systematic::Systematic(string name,int type,float up,float down){
  fName = name;
  fType = type;
  //
  fOverallUp = up;
  fOverallDown = down;
  //
  fWeightUp = "";
  fWeightSufUp = "";
  fNtuplePathsUp.clear();
  fNtuplePathSufUp = "";
  fNtupleFilesUp.clear();
  fNtupleFileSufUp = "";
  fNtupleNamesUp.clear();
  fNtupleNameSufUp = "";
  //
  fWeightDown = "";
  fWeightSufDown = "";  
  fNtuplePathsDown.clear();
  fNtuplePathSufDown = "";
  fNtupleFilesDown.clear();
  fNtupleFileSufDown = "";
  fNtupleNamesDown.clear();
  fNtupleNameSufDown = "";
  //
  fHistoPathsUp.clear();
  fHistoPathSufUp = "";
  fHistoFilesUp.clear();
  fHistoFileSufUp = "";
  fHistoNamesUp.clear();
  fHistoNameSufUp = "";
  //
  fHistoPathsDown.clear();
  fHistoPathSufDown = "";
  fHistoFilesDown.clear();
  fHistoFileSufDown = "";
  fHistoNamesDown.clear();
  fHistoNameSufDown = "";
}
Systematic::~Systematic(){}
