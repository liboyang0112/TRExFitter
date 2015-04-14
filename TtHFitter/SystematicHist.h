#include "TtHFitter/Common.h"

#ifndef __SystematicHist__
#define __SystematicHist__

// #include "TtHFitter/Sample.h"

class SystematicHist {
public:
  SystematicHist(string name);
  ~SystematicHist();

  void WriteToFile();
  void ReadFromFile();
  bool IsShape();
  
  void Print();
  
  string fName;

  bool fIsOverall;
  bool fIsShape;
  int fSmoothType;
  int fSymmetrisationType;

  TH1* fHistUp;
  TH1* fHistShapeUp;
  float fNormUp;
  string fFileNameUp;
  string fHistoNameUp;
  string fFileNameShapeUp;
  string fHistoNameShapeUp;
  TH1* fHistUp_original;

  TH1* fHistDown;
  TH1* fHistShapeDown;
  float fNormDown;
  string fFileNameDown;
  string fHistoNameDown;
  string fFileNameShapeDown;
  string fHistoNameShapeDown;
  TH1* fHistDown_original;
}; 

#endif
