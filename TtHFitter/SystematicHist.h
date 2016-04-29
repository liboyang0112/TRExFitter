#include "TtHFitter/Common.h"
#include "TtHFitter/Systematic.h"

#ifndef __SystematicHist__
#define __SystematicHist__

class SystematicHist {
public:
  SystematicHist(string name);
  ~SystematicHist();

  void WriteToFile();
  void ReadFromFile();
  bool IsShape();
  
  void Print();

  void Divide(TH1* h);
  void Divide(SystematicHist *syh);
  void Multiply(TH1* h);
  void Multiply(SystematicHist *syh);
  
  string fName;
  Systematic *fSystematic;

  bool fIsOverall;
  bool fIsShape;
  int fSmoothType;
  int fSymmetrisationType;

  TH1* fHistUp;
  TH1* fHistUp_orig;
  TH1* fHistShapeUp;
  float fNormUp;
  string fFileNameUp;
  string fHistoNameUp;
  string fFileNameShapeUp;
  string fHistoNameShapeUp;
  TH1* fHistUp_original;
  TH1* fHistUp_postFit;

  TH1* fHistDown;
  TH1* fHistDown_orig;
  TH1* fHistShapeDown;
  float fNormDown;
  string fFileNameDown;
  string fHistoNameDown;
  string fFileNameShapeDown;
  string fHistoNameShapeDown;
  TH1* fHistDown_original;
  TH1* fHistDown_postFit;
}; 

#endif
