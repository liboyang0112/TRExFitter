#include "TtHFitter/SystematicHist.h"

// -------------------------------------------------------------------------------------------------
// SystematicHist

SystematicHist::SystematicHist(string name){
  fName = name;

  fIsOverall = false;
  fIsShape = false;

  fHistUp = 0x0;
  fHistShapeUp = 0x0;
  fNormUp = 0;
  fFileNameUp = "";
  fHistoNameUp = "";
  fFileNameShapeUp = "";
  fHistoNameShapeUp = "";

  fHistDown = 0x0;
  fHistShapeDown = 0x0;
  fNormDown = 0;
  fFileNameDown = "";
  fHistoNameDown = "";
  fFileNameShapeDown = "";
  fHistoNameShapeDown = "";
}
SystematicHist::~SystematicHist(){}

void SystematicHist::WriteToFile(){
  WriteHistToFile(fHistUp,fFileNameUp);
  WriteHistToFile(fHistDown,fFileNameDown);
  if(fIsShape){
    WriteHistToFile(fHistShapeUp,fFileNameShapeUp);
    WriteHistToFile(fHistShapeDown,fFileNameShapeDown);
  }
}

void SystematicHist::ReadFromFile(){
  fHistUp = HistFromFile(fFileNameUp,fHistoNameUp);
  fHistShapeUp = HistFromFile(fFileNameShapeUp,fHistoNameShapeUp);
  fHistDown = HistFromFile(fFileNameDown,fHistoNameDown);
  fHistShapeDown = HistFromFile(fFileNameShapeDown,fHistoNameShapeDown);
}

bool SystematicHist::IsShape(){
  if(fHistUp!=0x0 || fHistDown!=0x0) return true;
  return false;
}

void SystematicHist::Print(){
  cout << "        Systematic: " << fName;
  if(fHistShapeUp==0x0 && fHistShapeDown==0x0 && fHistUp==0x0 && fHistDown==0x0) cout << Form("\toverall (%.3f,%.3f)",fNormUp,fNormDown) << endl;
  else cout << endl;
}
