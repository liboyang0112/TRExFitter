#include "TtHFitter/SystematicHist.h"
#include "TtHFitter/HistoTools.h"

// -------------------------------------------------------------------------------------------------
// SystematicHist

//_____________________________________________________________________________
//
SystematicHist::SystematicHist(string name){
    fName = name;
    fSystematic = 0x0;

    fIsOverall = false;
    fIsShape = false;
    fSmoothType = 0;
    fSymmetrisationType = 0;

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

//_____________________________________________________________________________
//
SystematicHist::~SystematicHist(){
    if(fHistUp) delete fHistUp;
    if(fHistShapeUp) delete fHistShapeUp;
    if(fHistDown) delete fHistDown;
    if(fHistShapeDown) delete fHistShapeDown;
}

//_____________________________________________________________________________
//
void SystematicHist::WriteToFile(){
    WriteHistToFile(fHistUp,fFileNameUp);
    WriteHistToFile(fHistDown,fFileNameDown);
    if(fIsShape){
        WriteHistToFile(fHistShapeUp,fFileNameShapeUp);
        WriteHistToFile(fHistShapeDown,fFileNameShapeDown);
        WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeUp),fFileNameShapeUp);
        WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeDown),fFileNameShapeDown);
    }
}

//_____________________________________________________________________________
//
void SystematicHist::ReadFromFile(){
    fHistUp = HistFromFile(fFileNameUp,fHistoNameUp);
    fHistShapeUp = HistFromFile(fFileNameShapeUp,fHistoNameShapeUp);
    fHistDown = HistFromFile(fFileNameDown,fHistoNameDown);
    fHistShapeDown = HistFromFile(fFileNameShapeDown,fHistoNameShapeDown);
}

//_____________________________________________________________________________
//
bool SystematicHist::IsShape(){
    if(fHistUp!=0x0 || fHistDown!=0x0) return true;
    return false;
}

//_____________________________________________________________________________
//
void SystematicHist::Print(){
    cout << "        Systematic: " << fName;
    if(fHistShapeUp==0x0 && fHistShapeDown==0x0 && fHistUp==0x0 && fHistDown==0x0) cout << Form("\toverall (%.3f,%.3f)",fNormUp,fNormDown) << endl;
    else cout << endl;
}
