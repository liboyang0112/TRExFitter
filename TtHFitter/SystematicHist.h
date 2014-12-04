#include "TH1.h"

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
//   Sample *fSample;

  bool fIsOverall;
  bool fIsShape;

  TH1* fHistUp;
  TH1* fHistShapeUp;
  float fNormUp;
  string fFileNameUp;
  string fHistoNameUp;
  string fFileNameShapeUp;
  string fHistoNameShapeUp;

  TH1* fHistDown;
  TH1* fHistShapeDown;
  float fNormDown;
  string fFileNameDown;
  string fHistoNameDown;
  string fFileNameShapeDown;
  string fHistoNameShapeDown;
}; 

#endif
