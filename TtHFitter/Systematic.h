#include "TtHFitter/Common.h"

#ifndef __Systematic__
#define __Systematic__

class Systematic {
public:
  Systematic(string name,int type=0,float up=0,float down=0);
  ~Systematic();
  
  string fName;
  string fTitle;
  int fType;
  int fSmoothType;
  int fSymmetrisationType;
    
  float fOverallUp;
  float fOverallDown;
  
  // from ntuples - up
  string fWeightUp;
  string fWeightSufUp;  
  vector<string> fNtuplePathsUp;
  string fNtuplePathSufUp;
  vector<string> fNtupleFilesUp;
  string fNtupleFileSufUp;
  vector<string> fNtupleNamesUp;
  string fNtupleNameSufUp;

  // from ntuples - down
  string fWeightDown;
  string fWeightSufDown;  
  vector<string> fNtuplePathsDown;
  string fNtuplePathSufDown;
  vector<string> fNtupleFilesDown;
  string fNtupleFileSufDown;
  vector<string> fNtupleNamesDown;
  string fNtupleNameSufDown;

  // from histos - up
  vector<string> fHistoPathsUp;
  string fHistoPathSufUp;
  vector<string> fHistoFilesUp;
  string fHistoFileSufUp;
  vector<string> fHistoNamesUp;
  string fHistoNameSufUp;
  
  // from histos - down
  vector<string> fHistoPathsDown;
  string fHistoPathSufDown;
  vector<string> fHistoFilesDown;
  string fHistoFileSufDown;
  vector<string> fHistoNamesDown;
  string fHistoNameSufDown;
};

#endif
