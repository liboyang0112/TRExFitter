#include "TH1.h"

#ifndef __SystematicHisto__
#define __SystematicHisto__

// #include "TtHFitter/Sample.h"

class SystematicHisto {
public:
  SystematicHisto(string name);
  ~SystematicHisto();

  void WriteToFile();
  void ReadFromFile();
  bool IsShape();
  
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
