// Class include
#include "TRExFitter/SystematicHist.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"

// -------------------------------------------------------------------------------------------------
// SystematicHist

//_____________________________________________________________________________
//
SystematicHist::SystematicHist(const std::string& name) :
    fName(name),
    fSystematic(nullptr),
    fIsOverall(false),
    fIsShape(false),
    fSmoothType(0),
    fSymmetrisationType(HistoTools::NOSYMMETRIZATION),
    fShapePruned(false),
    fNormPruned(false),
    fBadShape(false),
    fBadNorm(false),
    fHistUp(nullptr),
    fHistUp_orig(nullptr),
    fHistUp_preSmooth(nullptr),
    fHistShapeUp(nullptr),
    fNormUp(0),
    fFileNameUp(""),
    fHistoNameUp(""),
    fFileNameShapeUp(""),
    fHistoNameShapeUp(""),
    fHistUp_postFit(nullptr),
    fHistDown(nullptr),
    fHistDown_orig(nullptr),
    fHistDown_preSmooth(nullptr),
    fHistShapeDown(nullptr),
    fNormDown(0),
    fFileNameDown(""),
    fHistoNameDown(""),
    fFileNameShapeDown(""),
    fHistoNameShapeDown(""),
    fHistDown_postFit(nullptr),
    fScaleUp(1.),
    fScaleDown(1.) {
}

//_____________________________________________________________________________
//
SystematicHist::~SystematicHist(){
    delete fHistUp;
    delete fHistShapeUp;
    delete fHistDown;
    delete fHistShapeDown;
}

//_____________________________________________________________________________
//
void SystematicHist::WriteToFile(TFile *f,bool reWriteOrig) const{
    if(f==nullptr){
        WriteHistToFile(fHistUp,fFileNameUp);
        WriteHistToFile(fHistDown,fFileNameDown);
        if(reWriteOrig) WriteHistToFile(fHistUp_orig,fFileNameUp);
        if(reWriteOrig) WriteHistToFile(fHistDown_orig,fFileNameDown);
        if(fIsShape){
            WriteHistToFile(fHistShapeUp,fFileNameShapeUp);
            WriteHistToFile(fHistShapeDown,fFileNameShapeDown);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeUp).get(),fFileNameShapeUp);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeDown).get(),fFileNameShapeDown);
        }
        if(fSystematic->fType==Systematic::SHAPE){
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistUp).get(),fFileNameUp);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistDown).get(),fFileNameDown);
        }
    }
    else{
        WriteHistToFile(fHistUp,f);
        WriteHistToFile(fHistDown,f);
        if(reWriteOrig) WriteHistToFile(fHistUp_orig,f);
        if(reWriteOrig) WriteHistToFile(fHistDown_orig,f);
        if(fIsShape){
            WriteHistToFile(fHistShapeUp,f);
            WriteHistToFile(fHistShapeDown,f);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeUp).get(),f);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistShapeDown).get(),f);
        }
        if(fSystematic->fType==Systematic::SHAPE){
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistUp).get(),f);
            WriteHistToFile(HistoTools::TranformHistogramBinning(fHistDown).get(),f);
        }
    }
}

//_____________________________________________________________________________
//
void SystematicHist::ReadFromFile(){
    fHistUp      = HistFromFile(fFileNameUp,fHistoNameUp).release();
    fHistUp_orig = HistFromFile(fFileNameUp,fHistoNameUp+"_orig").release();
    if(fHistUp_orig==nullptr) fHistUp_orig = fHistUp;
    fHistShapeUp = HistFromFile(fFileNameShapeUp,fHistoNameShapeUp).release();
    fHistDown      = HistFromFile(fFileNameDown,fHistoNameDown).release();
    fHistDown_orig = HistFromFile(fFileNameDown,fHistoNameDown+"_orig").release();
    if(fHistDown_orig==nullptr) fHistDown_orig = fHistDown;
    fHistShapeDown = HistFromFile(fFileNameShapeDown,fHistoNameShapeDown).release();
}

//_____________________________________________________________________________
//
bool SystematicHist::IsShape() const{
    if(fHistUp!=nullptr || fHistDown!=nullptr) return true;
    return false;
}

//_____________________________________________________________________________
//
void SystematicHist::Print() const{
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
void SystematicHist::Add(TH1 *h,double scale){
    fHistUp->Add(h,scale);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Add(h,scale);
    fHistDown->Add(h,scale);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Add(h,scale);
}

//_____________________________________________________________________________
//
void SystematicHist::Add(SystematicHist *syh,double scale){
    fHistUp->Add(       syh->fHistUp,scale);
    if(fHistShapeUp!=nullptr)   {
        if (syh->fHistShapeUp != nullptr) fHistShapeUp->Add(  syh->fHistShapeUp,scale);
        else if (syh->fHistUp != nullptr){
          // the other syst is overall, while this is not: get by hand its dummy shape syst
          std::unique_ptr<TH1> htemp (static_cast<TH1*>(syh->fHistUp->Clone("hDummyShapeUp")));
          htemp->Scale(1.0/(1.0+syh->fNormUp));
          fHistShapeUp->Add(  htemp.get(),scale);
          htemp.reset(nullptr);
        }
    }
    fHistDown->Add(     syh->fHistDown,scale);
    if(fHistShapeDown!=nullptr)   {
        if (syh->fHistShapeDown != nullptr) fHistShapeDown->Add(  syh->fHistShapeDown,scale);
        else if (syh->fHistDown != nullptr){// the other syst is overall, while this is not
          std::unique_ptr<TH1> htemp (static_cast<TH1*>(syh->fHistDown->Clone("hDummyShapeDown")));
          htemp->Scale(1.0/(1.0+syh->fNormDown));
          fHistShapeDown->Add(  htemp.get(),scale);
          htemp.reset(nullptr);
        }
    }
}
