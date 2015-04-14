#include "TtHFitter/Systematic.h"

// -------------------------------------------------------------------------------------------------
// Systematic

Systematic::Systematic(string name,int type,float up,float down){
  fName = name;
  fType = type;
  fSmoothType = 0;
  fSymmetrisationType = 0;
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
Systematic::~Systematic(){
    fNtuplePathsUp.clear();
    fNtupleFilesUp.clear();
    fNtupleNamesUp.clear();
    fNtuplePathsDown.clear();
    fNtupleFilesDown.clear();
    fNtupleNamesDown.clear();
    fHistoPathsUp.clear();
    fHistoFilesUp.clear();
    fHistoNamesUp.clear(); 
    fHistoPathsDown.clear();
    fHistoFilesDown.clear();
    fHistoNamesDown.clear();
}
