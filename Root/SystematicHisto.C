#include "TtHFitter/SystematicHisto.h"

// -------------------------------------------------------------------------------------------------
// SystematicHisto

SystematicHisto::SystematicHisto(string name){
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
SystematicHisto::~SystematicHisto(){}

void SystematicHisto::WriteToFile(){
  WriteHistToFile(fHistUp,fFileNameUp);
  WriteHistToFile(fHistDown,fFileNameDown);
  if(fIsShape){
    WriteHistToFile(fHistShapeUp,fFileNameShapeUp);
    WriteHistToFile(fHistShapeDown,fFileNameShapeDown);
  }
}

void SystematicHisto::ReadFromFile(){
  fHistUp = HistFromFile(fFileNameUp,fHistoNameUp);
  fHistShapeUp = HistFromFile(fFileNameShapeUp,fHistoNameShapeUp);
  fHistDown = HistFromFile(fFileNameDown,fHistoNameDown);
  fHistShapeDown = HistFromFile(fFileNameShapeDown,fHistoNameShapeDown);
}

bool SystematicHisto::IsShape(){
  if(fHistUp!=0x0 || fHistDown!=0x0) return true;
  return false;
}
