// Class include
#include "TRExFitter/SystematicHist.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"

// ROOT includes
#include "TFile.h"
#include "TH1F.h"

using namespace std;

// -------------------------------------------------------------------------------------------------
// SystematicHist

//_____________________________________________________________________________
//
SystematicHist::SystematicHist(string name){
    fName = name;
    fSystematic = nullptr;

    fIsOverall = false;
    fIsShape = false;
    fSmoothType = 0;
    fSymmetrisationType = 0;

    fShapePruned = false;
    fNormPruned  = false;
    fBadShape    = false;
    fBadNorm     = false;

    fHistUp = nullptr;
    fHistUp_orig = nullptr;
    fHistShapeUp = nullptr;
    fNormUp = 0;
    fFileNameUp = "";
    fHistoNameUp = "";
    fFileNameShapeUp = "";
    fHistoNameShapeUp = "";
    fHistUp_postFit = nullptr;

    fHistDown = nullptr;
    fHistDown_orig = nullptr;
    fHistShapeDown = nullptr;
    fNormDown = 0;
    fFileNameDown = "";
    fHistoNameDown = "";
    fFileNameShapeDown = "";
    fHistoNameShapeDown = "";
    fHistDown_postFit = nullptr;

    fScaleUp   = 1.;
    fScaleDown = 1.;
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
void SystematicHist::WriteToFile(TFile *f){
    if(f==nullptr){
        WriteHistToFile(fHistUp,fFileNameUp);
        WriteHistToFile(fHistDown,fFileNameDown);
        WriteHistToFile(fHistUp_orig,fFileNameUp);
        WriteHistToFile(fHistDown_orig,fFileNameDown);
        if(fIsShape){
            WriteHistToFile(fHistShapeUp,fFileNameShapeUp);
            WriteHistToFile(fHistShapeDown,fFileNameShapeDown);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeUp),fFileNameShapeUp);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeDown),fFileNameShapeDown);
        }
        if(fSystematic->fType==Systematic::SHAPE){
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistUp),fFileNameUp);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistDown),fFileNameDown);
        }
    }
    else{
        WriteHistToFile(fHistUp,f);
        WriteHistToFile(fHistDown,f);
        WriteHistToFile(fHistUp_orig,f);
        WriteHistToFile(fHistDown_orig,f);
        if(fIsShape){
            WriteHistToFile(fHistShapeUp,f);
            WriteHistToFile(fHistShapeDown,f);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeUp),f);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeDown),f);
        }
        if(fSystematic->fType==Systematic::SHAPE){
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistUp),f);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistDown),f);
        }
    }
}

//_____________________________________________________________________________
//
void SystematicHist::ReadFromFile(){
    fHistUp      = HistFromFile(fFileNameUp,fHistoNameUp);
    fHistUp_orig = HistFromFile(fFileNameUp,fHistoNameUp+"_orig");
    if(fHistUp_orig==nullptr) fHistUp_orig = fHistUp;
    fHistShapeUp = HistFromFile(fFileNameShapeUp,fHistoNameShapeUp);
    fHistDown      = HistFromFile(fFileNameDown,fHistoNameDown);
    fHistDown_orig = HistFromFile(fFileNameDown,fHistoNameDown+"_orig");
    if(fHistDown_orig==nullptr) fHistDown_orig = fHistDown;
    fHistShapeDown = HistFromFile(fFileNameShapeDown,fHistoNameShapeDown);
}

//_____________________________________________________________________________
//
bool SystematicHist::IsShape(){
    if(fHistUp!=nullptr || fHistDown!=nullptr) return true;
    return false;
}

//_____________________________________________________________________________
//
void SystematicHist::Print(){
    std::string temp = "        Systematic: " + fName;
    if(fHistShapeUp==nullptr && fHistShapeDown==nullptr && fHistUp==nullptr && fHistDown==nullptr) temp + Form("\toverall (%.3f,%.3f)",fNormUp,fNormDown);
    WriteInfoStatus("SystematicHist::Print", temp);
}

//_____________________________________________________________________________
//
void SystematicHist::Divide(TH1 *h){
    fHistUp->Divide(h);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Divide(h);
    fHistDown->Divide(h);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Divide(h);
}

//_____________________________________________________________________________
//
void SystematicHist::Divide(SystematicHist *syh){
    fHistUp->Divide(       syh->fHistUp);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Divide(  syh->fHistShapeUp);
    fHistDown->Divide(     syh->fHistDown);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Divide(syh->fHistShapeDown);
}

//_____________________________________________________________________________
//
void SystematicHist::Multiply(TH1 *h){
    fHistUp->Multiply(h);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Multiply(h);
    fHistDown->Multiply(h);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Multiply(h);
}

//_____________________________________________________________________________
//
void SystematicHist::Multiply(SystematicHist *syh){
    fHistUp->Multiply(       syh->fHistUp);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Multiply(  syh->fHistShapeUp);
    fHistDown->Multiply(     syh->fHistDown);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Multiply(syh->fHistShapeDown);
}

//_____________________________________________________________________________
//
void SystematicHist::Add(TH1 *h,float scale){
    fHistUp->Add(h,scale);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Add(h,scale);
    fHistDown->Add(h,scale);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Add(h,scale);
}

//_____________________________________________________________________________
//
void SystematicHist::Add(SystematicHist *syh,float scale){
    fHistUp->Add(       syh->fHistUp,scale);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Add(  syh->fHistShapeUp,scale);
    fHistDown->Add(     syh->fHistDown,scale);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Add(syh->fHistShapeDown,scale);
}
