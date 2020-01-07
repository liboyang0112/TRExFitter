// Class include
#include "TRExFitter/SampleHist.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/HistoTools.h"

// ROOT includes
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

// c++ includes
#include <iostream>
#include <memory>

using namespace std;

// -------------------------------------------------------------------------------------------------
// SampleHist

//_____________________________________________________________________________
//
SampleHist::SampleHist() :
    fName(""),
    fSample(nullptr),
    fHist(nullptr),
    fHist_orig(nullptr),
    fHist_regBin(nullptr),
    fHist_preSmooth(nullptr),
    fHist_postFit(nullptr),
    fFileName(""),
    fHistoName(""),
    fIsData(false),
    fIsSig(false),
    fNSyst(0),
    fNNorm(0),
    fNShape(0),
    fFitName(""),
    fRegionName("Region"),
    fRegionLabel("Region"),
    fVariableTitle("Variable"),
    fSystSmoothed(false) {
}

//_____________________________________________________________________________
//
SampleHist::SampleHist(Sample *sample,TH1 *hist) : 
    fName(sample->fName),
    fSample(sample),
    fHist(nullptr),
    fHist_orig(nullptr),
    fHist_regBin(nullptr),
    fHist_preSmooth(nullptr),
    fHist_postFit(nullptr),
    fFileName(""),
    fHistoName(""),
    fIsData(false),
    fIsSig(false),
    fNSyst(0),
    fNNorm(0),
    fNShape(0),
    fFitName(""),
    fRegionName("Region"),
    fRegionLabel("Region"),
    fVariableTitle("Variable"),
    fSystSmoothed(false) {
    
    fHist = unique_ptr<TH1>(static_cast<TH1*>(hist->Clone(Form("h_%s",fName.c_str()))));
    fHist->SetFillColor(fSample->fFillColor);
    fHist->SetLineColor(fSample->fLineColor);
    fHist->SetLineWidth(1);

    fHist_orig = std::unique_ptr<TH1>(static_cast<TH1*>(fHist->Clone(Form("%s_orig",fHist->GetName()))));
    
    fIsMorph = fSample->fIsMorph;
}

//_____________________________________________________________________________
//
SampleHist::SampleHist(Sample *sample, const std::string& histoName, const std::string& fileName) :
    fName(sample->fName),
    fSample(sample),
    fHist(nullptr),
    fHist_orig(nullptr),
    fHist_regBin(nullptr),
    fHist_preSmooth(nullptr),
    fHist_postFit(nullptr),
    fFileName(fileName),
    fHistoName(histoName),
    fIsData(false),
    fIsSig(false),
    fNSyst(0),
    fNNorm(0),
    fNShape(0),
    fFitName(""),
    fRegionName("Region"),
    fRegionLabel("Region"),
    fVariableTitle("Variable"),
    fSystSmoothed(false) {

    fHist = Common::HistFromFile(fileName,histoName);

    if (fHist == nullptr) {
        WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
        exit(EXIT_FAILURE);
    }
    fHist->SetFillColor(fSample->fFillColor);
    fHist->SetLineColor(fSample->fLineColor);
    fHist->SetLineWidth(1);

    fHist_orig = Common::HistFromFile(fileName,histoName+"_orig");
    if(fHist_orig==nullptr){
        fHist_orig = std::unique_ptr<TH1>(static_cast<TH1*>(fHist->Clone(Form("%s_orig",fHist->GetName()))));
    }

    fIsMorph = fSample->fIsMorph;
}

//_____________________________________________________________________________
//
SampleHist::~SampleHist(){
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddOverallSyst(const std::string& name,const std::string& storedName,double up,double down){
    SystematicHist *syh = GetSystematic(name);
    // ... and if not create a new one
    if(!syh){
        fSyst.emplace_back(new SystematicHist(name));
        fNSyst ++;
        syh = fSyst.back().get();
    }
    //
    syh->fHistUp.reset( static_cast<TH1*>(fHist->Clone(Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistDown.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistUp  ->Scale(1.+up);
    syh->fHistDown->Scale(1.+down);
    syh->fHistUp_orig.reset(static_cast<TH1*>(fHist_orig->Clone(Form("%s_%s_%s_Up_orig",  fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistDown_orig.reset(static_cast<TH1*>(fHist_orig->Clone(Form("%s_%s_%s_Down_orig",fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistUp_orig  ->Scale(1.+up);
    syh->fHistDown_orig->Scale(1.+down);
    syh->fIsOverall = true;
    syh->fIsShape   = false;
    syh->fNormUp   = up;
    syh->fNormDown = down;
    return syh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddStatSyst(const std::string& name,const std::string& storedName, int i_bin) {
    int bin = i_bin+1; // counting of bins in Root starts with 1, in TRExFitter with 0
    SystematicHist *syh = GetSystematic(name);
    // ... and if not create a new one
    if(!syh){
        fSyst.emplace_back(new SystematicHist(name));
        fNSyst ++;
        syh = fSyst.back().get();
    }
    const double binContent = fHist->GetBinContent(bin);
    const double binError = binContent > 1e-4 ? fHist->GetBinError(bin) : 1e-7;
    syh->fHistUp.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistDown.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistShapeUp.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_%s_Shape_Up",  fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistShapeDown.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_%s_Shape_Down",fRegionName.c_str(),fSample->fName.c_str(),storedName.c_str()))));
    syh->fHistShapeUp  ->SetBinContent(bin, binContent + binError);
    syh->fHistShapeDown->SetBinContent(bin, binContent - binError);
    syh->fHistUp  ->SetBinContent(bin, binContent + binError);
    syh->fHistDown->SetBinContent(bin, binContent - binError);
    syh->fHistUp_orig.reset(static_cast<TH1*>(syh->fHistUp  ->Clone(Form("%s_orig",syh->fHistUp  ->GetName()))));
    syh->fHistDown_orig.reset(static_cast<TH1*>(syh->fHistDown->Clone(Form("%s_orig",syh->fHistDown->GetName()))));
    syh->fHistShapeUp  ->Scale(Common::EffIntegral(fHist.get()) / Common::EffIntegral(syh->fHistShapeUp.get()));
    syh->fHistShapeDown->Scale(Common::EffIntegral(fHist.get()) / Common::EffIntegral(syh->fHistShapeDown.get()));
    syh->fIsOverall = true;
    syh->fIsShape   = true;
    syh->fNormUp   = ( Common::EffIntegral(syh->fHistUp.get())   -  Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
    syh->fNormDown = ( Common::EffIntegral(syh->fHistDown.get()) - Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
    return syh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddHistoSyst(const std::string& name,const std::string& storedName,TH1* h_up,TH1* h_down){
    
    // before doing anything else, check if the sampleHist can be created
    if(h_up  ==nullptr) return nullptr;
    if(h_down==nullptr) return nullptr;

    SystematicHist *syh = GetSystematic(name);
    // ... and if not create a new one
    if(!syh){
        fSyst.emplace_back(new SystematicHist(name));
        fNSyst ++;
        syh = fSyst.back().get();
    }
    //
    syh->fHistUp.reset(static_cast<TH1*>(h_up  ->Clone(Form("%s_%s_Up",  fHist->GetName(),storedName.c_str()))));
    syh->fHistDown.reset(static_cast<TH1*>(h_down->Clone(Form("%s_%s_Down",fHist->GetName(),storedName.c_str()))));
    syh->fHistUp_orig.reset(static_cast<TH1*>(h_up  ->Clone(Form("%s_%s_Up_orig",  fHist->GetName(),storedName.c_str()))));
    syh->fHistDown_orig.reset(static_cast<TH1*>(h_down->Clone(Form("%s_%s_Down_orig",fHist->GetName(),storedName.c_str()))));
    syh->fHistUp_preSmooth.reset(static_cast<TH1*>(h_up  ->Clone(Form("%s_%s_Up_preSmooth",  fHist->GetName(),storedName.c_str()))));
    syh->fHistDown_preSmooth.reset(static_cast<TH1*>(h_down->Clone(Form("%s_%s_Down_preSmooth",fHist->GetName(),storedName.c_str()))));
    syh->fHistShapeUp.reset(static_cast<TH1*>(h_up  ->Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),storedName.c_str()))));
    syh->fHistShapeDown.reset(static_cast<TH1*>(h_down->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),storedName.c_str()))));
    if(Common::EffIntegral(syh->fHistShapeUp.get()) > 0. ){
        syh->fHistShapeUp  ->Scale(Common::EffIntegral(fHist.get()) / Common::EffIntegral(syh->fHistShapeUp.get()));
    } else {
        syh->fHistShapeUp.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),storedName.c_str()))));
    }
    if(Common::EffIntegral(syh->fHistShapeDown.get()) > 0. ){
        syh->fHistShapeDown->Scale(Common::EffIntegral(fHist.get()) / Common::EffIntegral(syh->fHistShapeDown.get()));
    } else {
        syh->fHistShapeDown.reset(static_cast<TH1*>(fHist->Clone(Form("%s_%s_Shape_Down",  fHist->GetName(),storedName.c_str()))));
    }

    syh->fIsOverall = true;
    syh->fIsShape   = true;
    syh->fNormUp   = ( Common::EffIntegral(syh->fHistUp.get())   -  Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
    syh->fNormDown = ( Common::EffIntegral(syh->fHistDown.get()) - Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
    if(syh->fNormUp == 0 && syh->fNormDown == 0) syh->fIsOverall = false;
    return syh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddHistoSyst(const std::string& name,
                                         const std::string& storedName,
                                         const std::string& histoName_up,
                                         const std::string& fileName_up,
                                         const std::string& histoName_down,
                                         const std:: string& fileName_down,
                                         int pruned/*1: norm only, 2: shape only*/){

    // before doing anything else, check if the sampleHist can be created
    std::unique_ptr<TH1> hUp   = Common::HistFromFile(fileName_up,  histoName_up);
    std::unique_ptr<TH1> hDown = Common::HistFromFile(fileName_down,histoName_down);
    if(hUp  ==nullptr) return nullptr;
    if(hDown==nullptr) return nullptr;

    SystematicHist *sh = GetSystematic(name);
    // ... and if not create a new one
    if(!sh){
        fSyst.emplace_back(new SystematicHist(name));
        fNSyst ++;
        sh = fSyst.back().get();
    }
    
    const bool normOnly  = (pruned==1);
    const bool shapeOnly = (pruned==2);
    
    sh->fFileNameUp   = fileName_up;
    sh->fFileNameDown = fileName_down;
    sh->fHistoNameUp   = histoName_up;
    sh->fHistoNameDown = histoName_down;
    sh->fHistUp   = Common::HistFromFile(sh->fFileNameUp,  sh->fHistoNameUp);
    sh->fHistDown = Common::HistFromFile(sh->fFileNameDown,sh->fHistoNameDown);
    sh->fHistUp_orig   = Common::HistFromFile(sh->fFileNameUp,  sh->fHistoNameUp  +"_orig");
    sh->fHistDown_orig = Common::HistFromFile(sh->fFileNameDown,sh->fHistoNameDown+"_orig");
    if(sh->fHistUp   == nullptr) return nullptr;
    if(sh->fHistDown == nullptr) return nullptr;
    if(sh->fHistUp_orig  ==nullptr) sh->fHistUp_orig.reset(static_cast<TH1D*>(sh->fHistUp->Clone(Form("%s_orig",sh->fHistUp->GetName()))));
    if(sh->fHistDown_orig==nullptr) sh->fHistDown_orig.reset(static_cast<TH1D*>(sh->fHistDown->Clone(Form("%s_orig",sh->fHistDown->GetName()))));
    //
    if(normOnly){
        sh->fIsShape   = false;
    }
    else{
        sh->fHistShapeUp.reset(static_cast<TH1*>(sh->fHistUp  ->Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),storedName.c_str()))));
        sh->fHistShapeDown.reset(static_cast<TH1*>(sh->fHistDown->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),storedName.c_str()))));
        if(Common::EffIntegral(sh->fHistShapeUp.get()) > 0. ){
            sh->fHistShapeUp  -> Scale(Common::EffIntegral(fHist.get()) / Common::EffIntegral(sh->fHistShapeUp.get()));
        } else {
            sh->fHistShapeUp.reset(static_cast<TH1*>(fHist -> Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),storedName.c_str()))));
        }

        if(Common::EffIntegral(sh->fHistShapeDown.get()) > 0. ){
            sh->fHistShapeDown->Scale(Common::EffIntegral(fHist.get()) / Common::EffIntegral(sh->fHistShapeDown.get()));
        } else {
            sh->fHistShapeDown.reset(static_cast<TH1*>(fHist -> Clone(Form("%s_%s_Shape_Down",  fHist->GetName(),storedName.c_str()))));
        }
        sh->fIsShape   = true;
    }
    //
    if(shapeOnly){
        sh->fIsOverall = false;
        sh->fNormUp   = 0;
        sh->fNormDown = 0;
    }
    else{
        sh->fNormUp   = ( Common::EffIntegral(sh->fHistUp.get())   -  Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
        sh->fNormDown = ( Common::EffIntegral(sh->fHistDown.get()) - Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
        sh->fIsOverall = true;
    }
    if(sh->fNormUp == 0 && sh->fNormDown == 0) sh->fIsOverall = false;
    //
    return sh;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::AddNormFactor(NormFactor *normFactor){
    NormFactor *norm = GetNormFactor(normFactor->fName);
    if(!norm){
        fNormFactors.emplace_back(std::move(normFactor));
        fNNorm ++;
        norm = fNormFactors.back().get();
    }
    else{
        norm = normFactor;
    }
    return norm;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::AddNormFactor(const std::string& name,double nominal, double min, double max){
    NormFactor *norm = GetNormFactor(name);
    if(!norm){
        fNormFactors.emplace_back(new NormFactor(name,nominal,min,max));
        fNNorm ++;
        norm = fNormFactors.back().get();
    }
    return norm;
}

//_____________________________________________________________________________
//
ShapeFactor* SampleHist::AddShapeFactor(ShapeFactor *shapeFactor){
    ShapeFactor *shape = GetShapeFactor(shapeFactor->fName);
    if(!shape){
        fShapeFactors.emplace_back(std::move(shapeFactor));
        fNShape ++;
        shape = fShapeFactors.back().get();
    } else {
        shape = shapeFactor;
    }
    return shape;
}

//_____________________________________________________________________________
//
ShapeFactor* SampleHist::AddShapeFactor(const std::string& name,double nominal, double min, double max){
    ShapeFactor *shape = GetShapeFactor(name);
    if(!shape){
        fShapeFactors.emplace_back(new ShapeFactor(name,nominal,min,max));
        fNShape ++;
        shape = fShapeFactors.back().get();
    }
    return shape;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::GetSystematic(const std::string& systName) const{
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(systName == fSyst[i_syst]->fName) return fSyst[i_syst].get();
    }
    return nullptr;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::GetSystFromNP(const std::string& NuisParName) const{
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(NuisParName == fSyst[i_syst]->fSystematic->fNuisanceParameter) return fSyst[i_syst].get();
    }
    return nullptr;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::GetNormFactor(const std::string& name) const{
    for(int i_syst=0;i_syst<fNNorm;i_syst++){
        if(name == fNormFactors[i_syst]->fName) return fNormFactors[i_syst].get();
    }
    return nullptr;
}

//_____________________________________________________________________________
//
ShapeFactor* SampleHist::GetShapeFactor(const std::string& name) const{
    for(int i_syst=0;i_syst<fNShape;i_syst++){
        if(name == fShapeFactors[i_syst]->fName) return fShapeFactors[i_syst].get();
    }
    return nullptr;
}

//_____________________________________________________________________________
//
bool SampleHist::HasSyst(const std::string& name) const{
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSyst[i_syst]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasNorm(const std::string& name) const{
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        if(fNormFactors[i_norm]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasShapeFactor(const std::string& name) const{
    for(int i_shape=0;i_shape<fNShape;i_shape++){
        if(fShapeFactors[i_shape]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
void SampleHist::WriteToFile(TFile *f,bool reWriteOrig){
    if(f==nullptr){
        if(fHist_orig!=nullptr && reWriteOrig)   Common::WriteHistToFile(fHist_orig.get(), fFileName);
        if(fHist!=nullptr)        Common::WriteHistToFile(fHist.get(), fFileName);
    }
    else{
        if(fHist_orig!=nullptr && reWriteOrig)   Common::WriteHistToFile(fHist_orig.get(), f);
        if(fHist!=nullptr)        Common::WriteHistToFile(fHist.get(), f);
    }
    // create the regular binning histogram
    fHist_regBin = std::unique_ptr<TH1>(HistoTools::TranformHistogramBinning(fHist.get()));
    if(fHist_regBin!=nullptr) Common::WriteHistToFile(fHist_regBin.get(),f);
    //
    // save separate gammas as histograms
    if(fSample->fSeparateGammas){
        TH1 *htempUp   = static_cast<TH1*>(fHist->Clone());
        TH1 *htempDown = static_cast<TH1*>(fHist->Clone());
        for(int i_bin=1;i_bin<=fHist->GetNbinsX();++i_bin) {
            htempUp  ->AddBinContent(i_bin, 1.*fHist->GetBinError(i_bin));
            htempDown->AddBinContent(i_bin,-1.*fHist->GetBinError(i_bin));
        }
        TH1 *htempUp_orig   = static_cast<TH1*>(fHist_orig->Clone());
        TH1 *htempDown_orig = static_cast<TH1*>(fHist_orig->Clone());
        for(int i_bin=1;i_bin<=fHist_orig->GetNbinsX();++i_bin) {
            htempUp_orig  ->AddBinContent(i_bin, 1.*fHist_orig->GetBinError(i_bin));
            htempDown_orig->AddBinContent(i_bin,-1.*fHist_orig->GetBinError(i_bin));
        }
        std::string systName = "stat_"+fSample->fName;
        Systematic *gamma = nullptr;
        if(GetSystematic(systName)) gamma = GetSystematic(systName)->fSystematic;  //GetSystematic(systName);
        if(gamma==nullptr) gamma = new Systematic(systName,Systematic::SHAPE);
        WriteDebugStatus("SampleHist::WriteToFile", "adding separate gammas as SHAPE systematic " + systName);
        gamma->fRegions.clear();
        gamma->fRegions.push_back(fRegionName);
        SystematicHist *syh = AddHistoSyst(systName,systName,htempUp,htempDown);
        if (!syh) {
            WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
            exit(EXIT_FAILURE);
        }
        syh->fHistUp_orig.reset(htempUp_orig);
        syh->fHistDown_orig.reset(htempDown_orig);
        gamma->fNuisanceParameter = gamma->fName;
        TRExFitter::NPMAP[gamma->fName] = gamma->fNuisanceParameter;
        syh->fSystematic = gamma;
        // setting histo name and file
        syh->fHistoNameUp = fRegionName+"_"+fSample->fName+"_stat_"+fSample->fName+"_Up";
        syh->fFileNameUp = fFileName;
    }
    //
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        // make sure they all have the correct name!
        fSyst[i_syst]->fHistUp  ->SetName( Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fSystematic->fStoredName.c_str()) );
        fSyst[i_syst]->fHistDown->SetName( Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fSystematic->fStoredName.c_str()) );
        fSyst[i_syst]->fHistUp_orig  ->SetName( Form("%s_%s_%s_Up_orig",  fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fSystematic->fStoredName.c_str()) );
        fSyst[i_syst]->fHistDown_orig->SetName( Form("%s_%s_%s_Down_orig",fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fSystematic->fStoredName.c_str()) );
        if(f==nullptr) fSyst[i_syst]->WriteToFile(nullptr,reWriteOrig);
        else           fSyst[i_syst]->WriteToFile(f,      reWriteOrig);
        // for shape hist, save also the syst(up)-nominal (to feed HistFactory)
        if(fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
            TH1* hVar = HistoTools::TranformHistogramBinning(
              static_cast<TH1*>(fSyst[i_syst]->fHistUp->Clone(Form("%s_%s_%s_Up_Var",  fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fSystematic->fStoredName.c_str())))
            );
            hVar->Add(fHist_regBin.get(),-1);
            hVar->Divide(fHist_regBin.get());
            // no negative bins here!
            for(int i_bin=1;i_bin<=hVar->GetNbinsX();i_bin++){
                if(hVar->GetBinContent(i_bin)<0) hVar->SetBinContent(i_bin,-1.*hVar->GetBinContent(i_bin));
            }
            if(f==nullptr) Common::WriteHistToFile(hVar,fFileName);
            else           Common::WriteHistToFile(hVar,f);
        }
    }
}

//_____________________________________________________________________________
//
void SampleHist::ReadFromFile(){
    fHist      = Common::HistFromFile(fFileName,fHistoName);
    fHist_orig = Common::HistFromFile(fFileName,fHistoName+"_orig");
}

//_____________________________________________________________________________
//
void SampleHist::NegativeTotalYieldWarning(TH1* hist, double yield) const{
    std::string temp = hist->GetName();
    WriteWarningStatus("SampleHist::NegativeTotalYieldWarning", "The total yield in " + temp + " is negative: " + std::to_string(yield));
    WriteWarningStatus("SampleHist::NegativeTotalYieldWarning", "    --> unable to preserve normalization while fixing bins with negative yields!");
}

//_____________________________________________________________________________
//
void SampleHist::FixEmptyBins(const bool suppress){
    //
    // store yields (nominal and systs)
    double initialYield = fHist->Integral();
    if (initialYield<0) { NegativeTotalYieldWarning(fHist.get(), initialYield); } // warning if total yield is negative, in which case normalization cannot be preserved
    vector<double> yieldUp;
    vector<double> yieldDown;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        SystematicHist* syh = fSyst[i_syst].get();
        if(syh==nullptr) continue;
        if(syh->fHistUp  ==nullptr) continue;
        if(syh->fHistDown==nullptr) continue;
        double tmpYieldUp   = syh->fHistUp->Integral();
        double tmpYieldDown = syh->fHistDown->Integral();
        yieldUp.push_back(tmpYieldUp);
        yieldDown.push_back(tmpYieldDown);
        // warnings if total yield in systematic variations is negative, in which case normalization cannot be preserved
        if (tmpYieldUp  <0){
            NegativeTotalYieldWarning(syh->fHistUp.get(),   tmpYieldUp);
        }
        if (tmpYieldDown<0){
            NegativeTotalYieldWarning(syh->fHistDown.get(), tmpYieldDown);
        }
    }
    //
    // store minimum stat unc for non-zero bins
    double minStat = -1;
    for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
        double content = fHist->GetBinContent(i_bin);
        double error   = fHist->GetBinError(  i_bin);
        if(content>0 && error>0){
            if(minStat<0 || error<minStat) minStat = error;
        }
    }
    //
    // loop o bins looking for negatives or zeros
    for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
        double content = fHist->GetBinContent(i_bin);
        double error   = fHist->GetBinError(  i_bin);
        if(content<=0){
            std::string temp = fHist->GetName();
            if (!suppress){
                WriteWarningStatus("SampleHist::FixEmptyBins", "Checking your nominal histogram " +temp + ", the bin " + std::to_string(i_bin) +
                " has a null/negative bin content (content = " + std::to_string(content) + ") ! You should have a look at this !");
                WriteWarningStatus("SampleHist::FixEmptyBins", "    --> For now setting this bin to 1e-06  +/- 1e-06!!! ");
            }
            // set nominal to 10^-6
            fHist->SetBinContent(i_bin,1e-6);
            if(error>0) {
                // if error defined, use it
                // this should in general make sense, even if the yield is negative and gets scaled to 1e-6
                fHist -> SetBinError(i_bin, error);
            } else {
                // try to guess stat. uncertainty if GUESSMCSTATERROR is enabled
                if(TRExFitter::GUESSMCSTATERROR){
                    if(minStat>0){
                        // if there was at least one bin with meaningful error, use the smallest
                        fHist -> SetBinError(i_bin, minStat);
                    } else {
                        // if not, give up and assign a meaningless error ;)
                        fHist -> SetBinError(i_bin, 1e-06);
                    }
                } else {
                    // no guessing, assign a 100% uncertainty
                    fHist -> SetBinError(i_bin, 1e-06);
                }
            }

            // loop on systematics and set them accordingly
            // uncertainties are not changed!
            for(int i_syst=0;i_syst<fNSyst;i_syst++){
                SystematicHist* syh = fSyst[i_syst].get();
                if(syh->fHistUp  ->GetBinContent(i_bin)<=0) syh->fHistUp  ->SetBinContent(i_bin,1e-06);
                if(syh->fHistDown->GetBinContent(i_bin)<=0) syh->fHistDown->SetBinContent(i_bin,1e-06);
            }
        }
    }
    // at this stage the integral should be strictly positive
    if(fHist->Integral()<0){
        WriteErrorStatus("SampleHist::FixEmptyBins", "Protection against negative yields failed, this should not happen");
        exit(EXIT_FAILURE);
    }
    // correct the overall normalization again, so it agrees with the initial status
    if(fHist->Integral()!=initialYield){
        if (initialYield>0) {
            // keep the original overall normalisation if the initial yield was positive
            double tmpScalingFactor = initialYield/fHist->Integral();
            for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
                fHist->SetBinContent(i_bin, fHist->GetBinContent(i_bin)*tmpScalingFactor);
            }
        } else if (TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL){
            // if the initial yield was negative, scale such that the total integral is 1e-06
            double tmpScalingFactor = 1e-6/fHist->Integral();
            for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
                fHist->SetBinContent(i_bin, fHist->GetBinContent(i_bin)*tmpScalingFactor);
            }
        }
    }
    // TODO apply the same logic also for the systematics histograms!
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        SystematicHist* syh = fSyst[i_syst].get();
        if(syh->fHistUp  ->Integral()!=yieldUp[i_syst]  ) syh->fHistUp  ->Scale(yieldUp[i_syst]  /syh->fHistUp  ->Integral());
        if(syh->fHistDown->Integral()!=yieldDown[i_syst]) syh->fHistDown->Scale(yieldDown[i_syst]/syh->fHistDown->Integral());
    }
}

//_____________________________________________________________________________
//
void SampleHist::Print() const{
    std::string temp = fHist->GetName();
    WriteDebugStatus("SampleHist::Print", "      Sample: " + fName + "\t" + temp);
    if(fNSyst>0){
        temp = "        Systematics:   ";
        for(int i_syst=0;i_syst<fNSyst;i_syst++){
            temp+= " " + fSyst[i_syst]->fName;
        }
        WriteDebugStatus("SampleHist::Print", temp);
    }
    if(fNNorm>0){
        temp = "        NormFactor(s): ";
        for(int i_norm=0;i_norm<fNNorm;i_norm++){
            temp+= " " + fNormFactors[i_norm]->fName;
        }
        WriteDebugStatus("SampleHist::Print", temp);
    }
    if(fNShape>0){
        temp = "        ShapeFactor(s): ";
        for(int i_shape=0;i_shape<fNShape;i_shape++){
            temp+= " " + fShapeFactors[i_shape]->fName;
        }
        WriteDebugStatus("SampleHist::Print", temp);
    }
}

//_____________________________________________________________________________
//
void SampleHist::Rebin(int ngroup, const Double_t* xbins){
    fHist->Rebin(ngroup,"",xbins);
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSyst[i_syst]->fHistUp!=nullptr) fSyst[i_syst]->fHistUp->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistDown!=nullptr) fSyst[i_syst]->fHistDown->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistShapeUp!=nullptr) fSyst[i_syst]->fHistShapeUp->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistShapeDown!=nullptr) fSyst[i_syst]->fHistShapeDown->Rebin(ngroup,"",xbins);
    }
}

//_____________________________________________________________________________
// this draws the control plots (for each systematic) with the syst variations for this region & all sample
void SampleHist::DrawSystPlot( const string &syst, TH1* const h_data, bool SumAndData, bool bothPanels ) const{
    if (SumAndData && h_data == nullptr){
        WriteWarningStatus("SampleHist::DrawSystPlot", "Data histogram passed is nullptr and you want to plot syst effect on data, returning.");
        WriteWarningStatus("SampleHist::DrawSystPlot", "Maybe you do not have data sample defined?");
        return;
    }

    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;
        
        TCanvas c("c","c",800,600);
        TPad pad0("pad0","pad0",0,0.30,1,1,0,0,0);
        pad0.SetTickx(true);
        pad0.SetTicky(true);
        pad0.SetTopMargin(0.05);
        pad0.SetBottomMargin(0.115);
        pad0.SetLeftMargin(0.14);
        pad0.SetRightMargin(0.04);
        pad0.SetFrameBorderMode(0);
        //
        TPad pad1("pad1","pad1",0,0,1,0.38,0,0,0);
        pad1.SetTickx(true);
        pad1.SetTicky(true);
        pad1.SetTopMargin(0.0);
        pad1.SetBottomMargin(0.27);
        pad1.SetLeftMargin(0.14);
        pad1.SetRightMargin(0.04);
        pad1.SetFrameBorderMode(0);
        //
        pad0.Draw();
        pad1.Draw();
        pad0.cd();

        std::unique_ptr<TH1> nominal(static_cast<TH1*>(fHist->Clone("nominal")));
        std::unique_ptr<TH1> nominal_orig(static_cast<TH1*>(fHist_preSmooth->Clone("nominal_orig")));
        std::unique_ptr<TH1> syst_up(static_cast<TH1*>(fSyst[i_syst]->fHistUp->Clone()));
        std::unique_ptr<TH1> syst_up_orig(static_cast<TH1*>(fSyst[i_syst]->fHistUp_preSmooth->Clone()));
        std::unique_ptr<TH1> syst_down(static_cast<TH1*>(fSyst[i_syst]->fHistDown->Clone()));
        std::unique_ptr<TH1> syst_down_orig(static_cast<TH1*>(fSyst[i_syst]->fHistDown_preSmooth->Clone()));
        std::unique_ptr<TH1> data(nullptr);
        if (SumAndData) data = std::unique_ptr<TH1>(static_cast<TH1*>(h_data->Clone("nominal")));
        std::unique_ptr<TH1> tmp(static_cast<TH1*>(nominal->Clone()));

        // drop shape or norm (for cases where this is not yet done in the stored histogrmas, i.e. in case of pruning or decorrelation)
        if (fSyst[i_syst] != nullptr && fSyst[i_syst]->fSystematic != nullptr) {
            if(fSyst[i_syst]->fSystematic->fIsNormOnly){
                Common::DropShape(syst_up.get(),syst_down.get(),nominal.get());
            }
            if(fSyst[i_syst]->fSystematic->fIsShapeOnly){
                Common::DropNorm(syst_up.get(),syst_down.get(),nominal.get());
            }
        }
        
        // Cosmetics
        nominal->SetLineColor(kBlack);
        nominal->SetLineWidth(2);
        nominal->SetFillColor(0);
        nominal_orig->SetLineColor(kBlack);
        nominal_orig->SetLineStyle(2);
        nominal_orig->SetLineWidth(2);
        nominal_orig->SetFillColor(0);
        nominal->SetMinimum(0);
        syst_up->SetLineColor(kRed);
        syst_up->SetLineWidth(2);
        syst_up->SetLineStyle(1);
        syst_up->SetFillStyle(0);
        syst_down->SetLineColor(kBlue);
        syst_down->SetLineWidth(2);
        syst_down->SetLineStyle(1);
        syst_down->SetFillStyle(0);
        syst_up_orig->SetLineColor(kRed);
        syst_up_orig->SetLineWidth(2);
        syst_up_orig->SetLineStyle(2);
        syst_up_orig->SetFillStyle(0);
        syst_down_orig->SetLineColor(kBlue);
        syst_down_orig->SetLineWidth(2);
        syst_down_orig->SetLineStyle(2);
        syst_down_orig->SetFillStyle(0);
        tmp->Scale(0);
        tmp->SetFillColor(0);
        if (SumAndData) data->SetMarkerColor(kBlack);

        // make copies for ratio
        std::unique_ptr<TH1> nominal_ratio(static_cast<TH1*>(nominal->Clone()));      
        std::unique_ptr<TH1> nominal_orig_ratio(static_cast<TH1*>(nominal_orig->Clone()));      
        std::unique_ptr<TH1> syst_up_ratio(static_cast<TH1*>(syst_up->Clone()));      
        std::unique_ptr<TH1> syst_up_orig_ratio(static_cast<TH1*>(syst_up_orig->Clone()));      
        std::unique_ptr<TH1> syst_down_ratio(static_cast<TH1*>(syst_down->Clone()));      
        std::unique_ptr<TH1> syst_down_orig_ratio(static_cast<TH1*>(syst_down_orig->Clone()));      
        std::unique_ptr<TH1> data_ratio(nullptr);
        if (SumAndData) data_ratio = std::unique_ptr<TH1>(static_cast<TH1*>(data->Clone()));      
        std::unique_ptr<TH1> tmp_ratio(static_cast<TH1*>(tmp->Clone()));      
        
        DrawSystPlotUpper(&pad0,
                          nominal.get(),
                          nominal_orig.get(),
                          syst_up.get(),
                          syst_up_orig.get(),
                          syst_down.get(),
                          syst_down_orig.get(),
                          data.get(),
                          tmp.get(),
                          SumAndData,
                          bothPanels);

        // Draw laels
        TLatex tex{};
        tex.SetNDC();
        if(SumAndData) {
            if(fSyst[i_syst]->fSystematic) tex.DrawLatex(0.17,0.79,Form("%s",fSyst[i_syst]->fSystematic->fTitle.c_str()));
            else                           tex.DrawLatex(0.17,0.79,Form("%s",fSyst[i_syst]->fName.c_str()));
        } else{
            if(fSyst[i_syst]->fSystematic) tex.DrawLatex(0.17,0.79,Form("%s, %s",fSyst[i_syst]->fSystematic->fTitle.c_str(),fSample->fTitle.c_str()));
            else                           tex.DrawLatex(0.17,0.79,Form("%s, %s",fSyst[i_syst]->fName.c_str(),fSample->fTitle.c_str()));
        }
        tex.DrawLatex(0.17,0.72,fRegionLabel.c_str());

        std::unique_ptr<TLegend> leg(nullptr);
        if(SumAndData) leg = std::make_unique<TLegend>(0.7,0.71,0.9,0.9);
        else           leg = std::make_unique<TLegend>(0.7,0.71,0.9,0.85);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetMargin(0.2);

        const float yield_nominal = Common::CorrectIntegral(nominal.get());
        const float yield_up      = Common::CorrectIntegral(syst_up.get());
        const float yield_down    = Common::CorrectIntegral(syst_down.get());
        float yield_data(0);
        if (SumAndData) yield_data = Common::CorrectIntegral(data.get());
        const float acc_up = (yield_up-yield_nominal)/yield_nominal;
        const float acc_down = (yield_down-yield_nominal)/yield_nominal;
        std::string sign_up =  "+";
        if(acc_up<0) sign_up = "-";
        std::string sign_down =  "+";
        if(acc_down<0) sign_down = "-";
        leg->AddEntry(syst_up.get(),  Form("+ 1 #sigma (%s%.1f %%)",sign_up.c_str(),  std::fabs(acc_up  *100)),"l");
        leg->AddEntry(syst_down.get(),Form(" - 1 #sigma (%s%.1f %%)",sign_down.c_str(),std::fabs(acc_down*100)),"l");
        leg->Draw("same");

        //Legend to define the line style
        TLegend leg2(0.65,0.64,0.9,0.7);
        leg2.SetFillStyle(0);
        leg2.SetBorderSize(0);
        leg2.SetNColumns(2);
        leg2.SetTextSize(gStyle->GetTextSize());
        leg2.SetTextFont(gStyle->GetTextFont());
        std::unique_ptr<TH1D> syst_up_black(static_cast<TH1D*>(syst_up->Clone()));
        std::unique_ptr<TH1D> syst_up_origin_black(static_cast<TH1D*>(syst_up_orig->Clone()));
        syst_up_black->SetLineColor(kBlack);
        syst_up_origin_black->SetLineColor(kBlack);
        leg2.AddEntry(syst_up_origin_black.get(),"Original","l");
        leg2.AddEntry(syst_up_black.get(),"Modified","l");
        leg2.Draw("same");

        std::unique_ptr<TLegend> leg3(nullptr);
        if(SumAndData){
            float acc_data = 0.;
            if (std::fabs(yield_nominal) > 1e-6) acc_data = (yield_data-yield_nominal)/yield_nominal;
            else acc_data = 99999999;
            std::string sign_data =  "+";
            if(acc_data<0) sign_data = "-";
            leg3 = std::make_unique<TLegend>(0.7,0.43,0.9,0.62);
            leg3->SetFillStyle(0);
            leg3->SetBorderSize(0);
            leg3->SetTextSize(gStyle->GetTextSize());
            leg3->SetTextFont(gStyle->GetTextFont());
            leg3->SetMargin(0.2);
            leg3->AddEntry(data.get(),"Data","p");
            leg3->AddEntry(nominal.get(),"Total prediction","l");
            leg3->Draw("same");
        }
        
 
        DrawSystPlotRatio(&pad1,
                          nominal_ratio.get(),
                          nominal_orig_ratio.get(),
                          syst_up_ratio.get(),
                          syst_up_orig_ratio.get(),
                          syst_down_ratio.get(),
                          syst_down_orig_ratio.get(),
                          data_ratio.get(),
                          tmp_ratio.get(),
                          SumAndData);

        const float xmin = nominal->GetBinLowEdge(1);
        const float xmax = nominal->GetBinLowEdge(nominal->GetNbinsX()+1);
        TLine one(xmin,0.,xmax,0.);
        one.SetLineColor(kBlack);
        one.SetLineWidth(2);
        one.Draw("same HIST");
        
        TLine line(0.01,1,0.1,1);
        line.SetLineColor(kWhite);
        line.SetLineWidth(20);
        line.DrawLineNDC(0.07,1,0.135,1);
        
        /// Make folders
        gSystem->mkdir(fFitName.c_str());
        gSystem->mkdir((fFitName+"/Systematics").c_str());
        gSystem->mkdir((fFitName+"/Systematics/"+fSyst[i_syst]->fName).c_str());

        for(std::size_t i_format=0; i_format < TRExFitter::IMAGEFORMAT.size(); ++i_format){
            if(SumAndData) {
                c.SaveAs(Form("%s/Systematics/%s/%s_%s.%s",fFitName.c_str(),fSyst[i_syst]->fName.c_str(), fName.c_str(), fSyst[i_syst]->fName.c_str(), TRExFitter::IMAGEFORMAT[i_format].c_str()));
            } else { 
                c.SaveAs(Form("%s/Systematics/%s/%s_%s.%s",fFitName.c_str(),fSyst[i_syst]->fName.c_str(),fHist->GetName(), fSyst[i_syst]->fName.c_str(), TRExFitter::IMAGEFORMAT[i_format].c_str()));
            }
        }

    }
}

//_____________________________________________________________________________
//
void SampleHist::SmoothSyst(const HistoTools::SmoothOption &smoothOpt, string syst, bool force){
    if(fSystSmoothed && !force) return;
    TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
    TH1* h_syst_up;
    TH1* h_syst_down;

    for(int i_syst=0;i_syst<fNSyst;i_syst++){

        if(syst!="all" && fSyst[i_syst]->fName != syst) continue;

        if(fSyst[i_syst]->fHistUp  ==nullptr) continue;
        if(fSyst[i_syst]->fHistDown==nullptr) continue;
        if(fSyst[i_syst]->fSystematic==nullptr) continue;

        h_syst_up = static_cast<TH1*>(fSyst[i_syst]->fHistUp->Clone());
        h_syst_down = static_cast<TH1*>(fSyst[i_syst]->fHistDown->Clone());

        if(fSyst[i_syst]->fSmoothType + fSyst[i_syst]->fSymmetrisationType<=0){
            HistoTools::Scale(fSyst[i_syst]->fHistUp.get(),   fHist.get(),fSyst[i_syst]->fScaleUp);
            HistoTools::Scale(fSyst[i_syst]->fHistDown.get(), fHist.get(),fSyst[i_syst]->fScaleDown);
            continue;
        }

        //
        // Pre-smoothing
        //
        // (do smoothing and symmetrization before pre-smoothing in case of two-sided systematics)
        if(fSyst[i_syst]->fSymmetrisationType==HistoTools::SYMMETRIZETWOSIDED){
            if(fSyst[i_syst]->fIsShape){
                if(fSyst[i_syst]->fSystematic->fSampleSmoothing){
                    HistoTools::ManageHistograms(   fSyst[i_syst]->fSmoothType,
                                                    fSyst[i_syst]->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal,//nominal histogram
                                                    fSyst[i_syst]->fHistUp.get(),
                                                    fSyst[i_syst]->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    fSyst[i_syst]->fScaleUp,
                                                    fSyst[i_syst]->fScaleDown, // scale factors
                                                    fSyst[i_syst]->fSystematic->fSampleSmoothOption // overwrite smoothing option
                                                );
                }
                else{
                    HistoTools::ManageHistograms(   fSyst[i_syst]->fSmoothType,
                                                    fSyst[i_syst]->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal,//nominal histogram
                                                    fSyst[i_syst]->fHistUp.get(),
                                                    fSyst[i_syst]->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    fSyst[i_syst]->fScaleUp,
                                                    fSyst[i_syst]->fScaleDown, // scale factors
                                                    smoothOpt
                                                );
                }
            }
            //
            // need to ad these lines to make sure overall only systematics get scaled as well
            else{
                HistoTools::Scale(fSyst[i_syst]->fHistUp.get(),  fHist.get(),fSyst[i_syst]->fScaleUp);
                HistoTools::Scale(fSyst[i_syst]->fHistDown.get(),fHist.get(),fSyst[i_syst]->fScaleDown);
            }
        }

        if(fSyst[i_syst]->fSystematic->fPreSmoothing){
            TH1* h_tmp_up   = h_syst_up!=nullptr   ? (TH1*)h_syst_up  ->Clone() : nullptr;
            TH1* h_tmp_down = h_syst_down!=nullptr ? (TH1*)h_syst_down->Clone() : nullptr;
            if(h_tmp_up!=nullptr || h_tmp_down!=nullptr){
                TH1* h_tmp_nominal = (TH1*)h_nominal  ->Clone();
                for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++){
                    h_tmp_nominal->GetBinError(i_bin,0);
                }
                if(h_tmp_up!=nullptr){
                    double tmp_nom_up = h_tmp_up->Integral();
                    h_tmp_up->Add(h_tmp_nominal,-1);
                    h_tmp_up->Divide(h_tmp_nominal);
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_up->AddBinContent(i_bin, 100.);
                    h_tmp_up->Smooth();
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_up->AddBinContent(i_bin,-100.);
                    h_tmp_up->Multiply(h_tmp_nominal);
                    h_tmp_up->Add(h_tmp_nominal, 1);
                    h_tmp_up->Scale(tmp_nom_up/h_tmp_up->Integral());
                    h_syst_up = (TH1*)h_tmp_up->Clone();
                }
                if(h_tmp_down!=nullptr){
                    double tmp_nom_down = h_tmp_down->Integral();
                    h_tmp_down->Add(h_tmp_nominal,-1);
                    h_tmp_up->Divide(h_tmp_nominal);
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_down->AddBinContent(i_bin, 100.);
                    h_tmp_down->Smooth();
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_down->AddBinContent(i_bin,-100.);
                    h_tmp_up->Multiply(h_tmp_nominal);
                    h_tmp_down->Add(h_tmp_nominal, 1);
                    h_tmp_down->Scale(tmp_nom_down/h_tmp_down->Integral());
                    h_syst_down = (TH1*)h_tmp_down->Clone();
                }
                delete h_tmp_nominal;
            }
            delete h_tmp_up;
            delete h_tmp_down;
        }

        //
        // Call the function for smoothing and symmetrisation
        //
        if(fSyst[i_syst]->fSymmetrisationType!=HistoTools::SYMMETRIZETWOSIDED){
            if(fSyst[i_syst]->fIsShape){
                if(fSyst[i_syst]->fSystematic->fSampleSmoothing){
                    HistoTools::ManageHistograms(   fSyst[i_syst]->fSmoothType,
                                                    fSyst[i_syst]->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal,//nominal histogram
                                                    fSyst[i_syst]->fHistUp.get(),
                                                    fSyst[i_syst]->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    fSyst[i_syst]->fScaleUp,
                                                    fSyst[i_syst]->fScaleDown, // scale factors
                                                    fSyst[i_syst]->fSystematic->fSampleSmoothOption // overwrite smoothing option
                                                );
                }
                else{
                    HistoTools::ManageHistograms(   fSyst[i_syst]->fSmoothType,
                                                    fSyst[i_syst]->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal,//nominal histogram
                                                    fSyst[i_syst]->fHistUp.get(),
                                                    fSyst[i_syst]->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    fSyst[i_syst]->fScaleUp,
                                                    fSyst[i_syst]->fScaleDown, // scale factors
                                                    smoothOpt
                                                );
                }
            }
            //
            // need to ad these lines to make sure overall only systematics get scaled as well
            else{
                HistoTools::Scale(fSyst[i_syst]->fHistUp.get(),  fHist.get(),fSyst[i_syst]->fScaleUp);
                HistoTools::Scale(fSyst[i_syst]->fHistDown.get(),fHist.get(),fSyst[i_syst]->fScaleDown);
            }
        }

        //
        // keep the variation below 100% in each bin, if the option Smooth is set for the sample
        //
        if(fSample->fSmooth){
            if(h_syst_up!=nullptr){
                for(int iBin = 1; iBin <= h_syst_up  ->GetNbinsX(); ++iBin ){
                    double relDiff = (h_syst_up->GetBinContent(iBin) - h_nominal->GetBinContent(iBin))/ h_nominal->GetBinContent(iBin);
                    if(relDiff>=1. ) h_syst_up->SetBinContent(iBin, (1.+0.99)*h_nominal->GetBinContent(iBin) );
                    if(relDiff<=-1.) h_syst_up->SetBinContent(iBin, (1.-0.99)*h_nominal->GetBinContent(iBin) );
                }
            }
            if(h_syst_down!=nullptr){
                for(int iBin = 1; iBin <= h_syst_down  ->GetNbinsX(); ++iBin ){
                    double relDiff = (h_syst_down->GetBinContent(iBin) - h_nominal->GetBinContent(iBin))/ h_nominal->GetBinContent(iBin);
                    if(relDiff>=1. ) h_syst_down->SetBinContent(iBin, (1.+0.99)*h_nominal->GetBinContent(iBin) );
                    if(relDiff<=-1.) h_syst_down->SetBinContent(iBin, (1.-0.99)*h_nominal->GetBinContent(iBin) );
                }
            }
        }

        //
        // Save stuff
        //
        fSyst[i_syst]->fHistUp.reset(static_cast<TH1*>(h_syst_up->Clone(fSyst[i_syst]->fHistUp->GetName())));
        fSyst[i_syst]->fHistDown.reset(static_cast<TH1*>(h_syst_down->Clone(fSyst[i_syst]->fHistDown->GetName())));

        //
        // Perform a check of the output histograms (check for 0 bins and other pathologic behaviour)
        //
        HistoTools::CheckHistograms( h_nominal /*nominal*/, fSyst[i_syst].get() /*systematic*/, fSample -> fType != Sample::SIGNAL, TRExFitter::HISTOCHECKCRASH /*cause crash if problem*/);

        //
        // Normalisation component first
        //
        if(h_nominal->Integral()!=0){
            fSyst[i_syst]->fNormUp   = fSyst[i_syst]->fHistUp  ->Integral()/h_nominal->Integral() - 1.;
            fSyst[i_syst]->fNormDown = fSyst[i_syst]->fHistDown->Integral()/h_nominal->Integral() - 1.;
        } else {
            WriteErrorStatus("SampleHist::SmoothSyst", "A nominal histogram with 0 integral has been found. Please check ! ");
            WriteErrorStatus("SampleHist::SmoothSyst", "            -> Sample: " + fName);
        }

        if(fSyst[i_syst]->fIsShape){
            // update shape hists as well
            fSyst[i_syst]->fHistShapeUp.reset(static_cast<TH1*>(h_syst_up  ->Clone(fSyst[i_syst]->fHistShapeUp->GetName())));
            fSyst[i_syst]->fHistShapeDown.reset(static_cast<TH1*>(h_syst_down->Clone(fSyst[i_syst]->fHistShapeDown->GetName())));
            if(fSyst[i_syst]->fHistShapeUp  ->Integral()>0){
                fSyst[i_syst]->fHistShapeUp  ->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeUp  ->Integral());
            } else {
                fSyst[i_syst]->fHistShapeUp.reset(static_cast<TH1*>(fHist ->Clone(fSyst[i_syst]->fHistShapeUp->GetName())));
            }

            if(fSyst[i_syst]->fHistShapeDown->Integral() > 0.){
                fSyst[i_syst]->fHistShapeDown->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeDown->Integral());
            } else {
                fSyst[i_syst]->fHistShapeDown.reset(static_cast<TH1*>(fHist ->Clone(fSyst[i_syst]->fHistShapeDown->GetName())));
            }
        }
    }
    fSystSmoothed = true;
}

//_____________________________________________________________________________
//
void SampleHist::CloneSampleHist(SampleHist* h, const std::set<std::string>& names, double scale){
    fName = h->fName;
    fHist           = std::unique_ptr<TH1>(static_cast<TH1*>(h->fHist->Clone()));
    fHist_preSmooth = std::unique_ptr<TH1>(static_cast<TH1*>(h->fHist_preSmooth->Clone()));
    fHist_orig      = std::unique_ptr<TH1>(static_cast<TH1*>(h->fHist_orig->Clone()));
    fHist->Scale(scale);
    fHist_preSmooth->Scale(scale);
    fHist_orig->Scale(scale);
    fFileName = h->fFileName;
    fHistoName = h->fHistoName;
    fIsData = h->fIsData;
    fIsSig = h->fIsSig;
    fNSyst = h->fNSyst;
    for(const auto& systname : names){
        bool notFound=true;
        for(int i_syst=0; i_syst<h->fNSyst; i_syst++){
            SystematicHist* syst_tmp = new SystematicHist("tmp");
            if(systname!=h->fSyst[i_syst]->fName) continue;
            TH1* tmp = static_cast<TH1*>(h->fSyst[i_syst]->fHistUp->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistUp.reset(tmp);

            tmp = static_cast<TH1*>(h->fSyst[i_syst]->fHistUp_preSmooth->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistUp_preSmooth.reset(tmp);

            tmp = static_cast<TH1*>(h->fSyst[i_syst]->fHistUp_orig->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistUp_orig.reset(tmp);

            tmp = static_cast<TH1*>(h->fSyst[i_syst]->fHistDown->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistDown.reset(tmp);

            tmp = static_cast<TH1*>(h->fSyst[i_syst]->fHistDown_preSmooth->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistDown_preSmooth.reset(tmp);

            tmp = static_cast<TH1*>(h->fSyst[i_syst]->fHistDown_orig->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistDown_orig.reset(tmp);

            syst_tmp->fName = h->fSyst[i_syst]->fName;
            fSyst.emplace_back(std::move(syst_tmp));
            notFound=false;
        }
        if(notFound){
            SystematicHist* syst_tmp = new SystematicHist("tmp");
            ++fNSyst;
            syst_tmp->fHistUp.reset(static_cast<TH1*>(h->fHist->Clone()));
            syst_tmp->fHistUp_orig.reset(static_cast<TH1*>(h->fHist_orig->Clone()));
            syst_tmp->fHistDown.reset(static_cast<TH1*>(h->fHist->Clone()));
            syst_tmp->fHistDown_orig.reset(static_cast<TH1*>(h->fHist_orig->Clone()));
            syst_tmp->fName = systname;
            fSyst.emplace_back(std::move(syst_tmp));
        }
    }

    fFitName = h->fFitName;
    fRegionName = h->fRegionName;
    fRegionLabel = h->fRegionLabel;
    fVariableTitle = h->fVariableTitle;
    fSystSmoothed = h->fSystSmoothed;
}

//_____________________________________________________________________________
//
void SampleHist::SampleHistAdd(SampleHist* h, double scale){
    fHist          ->Add(h->fHist.get(),          scale);
    fHist_preSmooth->Add(h->fHist_preSmooth.get(),scale);
    fHist_orig     ->Add(h->fHist_orig.get(),     scale);
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        bool wasIn = false;
        for(int j_syst=0;j_syst<h->fNSyst;j_syst++){
            if(fSyst[i_syst]->fName==h->fSyst[j_syst]->fName){
                fSyst[i_syst]->fHistUp  ->Add(h->fSyst[j_syst]->fHistUp.get(),  scale);
                fSyst[i_syst]->fHistDown->Add(h->fSyst[j_syst]->fHistDown.get(),scale);
                if(fSyst[i_syst]->fHistUp_preSmooth!=nullptr)   fSyst[i_syst]->fHistUp_preSmooth  ->Add(h->fSyst[j_syst]->fHistUp_preSmooth.get(),  scale);
                else                                            fSyst[i_syst]->fHistUp_preSmooth.reset(static_cast<TH1*>(fHist_preSmooth->Clone()));
                if(fSyst[i_syst]->fHistDown_preSmooth!=nullptr) fSyst[i_syst]->fHistDown_preSmooth->Add(h->fSyst[j_syst]->fHistDown_preSmooth.get(),scale);
                else                                            fSyst[i_syst]->fHistDown_preSmooth.reset(static_cast<TH1*>(fHist_preSmooth->Clone()));
                fSyst[i_syst]->fHistUp_orig  ->Add(h->fSyst[j_syst]->fHistUp_orig.get(),  scale);
                fSyst[i_syst]->fHistDown_orig->Add(h->fSyst[j_syst]->fHistDown_orig.get(),scale);
                wasIn = true;
            }
        }
        if(wasIn) continue;
        fSyst[i_syst]->fHistUp  ->Add(h->fHist.get(),scale);
        fSyst[i_syst]->fHistDown->Add(h->fHist.get(),scale);
        if(fSyst[i_syst]->fHistUp_preSmooth!=nullptr)   fSyst[i_syst]->fHistUp_preSmooth  ->Add(h->fHist_preSmooth.get(),scale);
        else                                            fSyst[i_syst]->fHistUp_preSmooth.reset(static_cast<TH1*>(fHist_preSmooth->Clone()));
        if(fSyst[i_syst]->fHistDown_preSmooth!=nullptr) fSyst[i_syst]->fHistDown_preSmooth->Add(h->fHist_preSmooth.get(),scale);
        else                                            fSyst[i_syst]->fHistDown_preSmooth.reset(static_cast<TH1*>(fHist_preSmooth->Clone()));
        fSyst[i_syst]->fHistUp_orig  ->Add(h->fHist_orig.get(),scale );
        fSyst[i_syst]->fHistDown_orig->Add(h->fHist_orig.get(),scale);
    }
}

//_____________________________________________________________________________
//
void SampleHist::Divide(SampleHist *sh){
    if (sh->fHist!=nullptr) fHist->Divide( sh->fHist.get() );
    else  {
       if (TRExFitter::HISTOCHECKCRASH) {
            WriteErrorStatus("SampleHist::Divide", "Sample "+sh->fName+ " not found when trying to divide "+fSample->fName+" by it");
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("SampleHist::Divide", "Sample "+sh->fName+ " not found when trying to divide "+fSample->fName+" by it");
        }
    }

    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        const std::string systName = fSyst[i_syst]->fName;
        const std::string NuisParName = fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = sh->GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Divide", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Divide", "Using its nominal. ");
            fSyst[i_syst]->Divide( sh->fHist.get() );
        }
        else{
            WriteDebugStatus("SampleHist::Divide", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Divide", "Properly computing with that. ");
            fSyst[i_syst]->Divide( syh );
        }
    }
    // loop on all the systematics in the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        const std::string systName = sh->fSyst[i_syst]->fName;
        const std::string NuisParName = sh->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Divide", "Adding syst "+ NuisParName + " (through syst "+ systName + ") to sample "+ fName);
            std::unique_ptr<TH1> hUp(static_cast<TH1*>(fHist->Clone("h_tmp_up")));
            std::unique_ptr<TH1> hDown(static_cast<TH1*>(fHist->Clone("h_tmp_down")));
            hUp  ->Divide(  sh->fHist.get() );
            hUp  ->Multiply(sh->fSyst[i_syst]->fHistUp.get());
            hUp  ->Scale(-1);
            hUp  ->Add(fHist.get(),2);
            //
            hDown->Divide(  sh->fHist.get() );
            hDown->Multiply(sh->fSyst[i_syst]->fHistDown.get());
            hDown->Scale(-1);
            hDown->Add(fHist.get(),2);
            //
            syh = AddHistoSyst(NuisParName,NuisParName,hUp.get(),hDown.get());
            if (syh == nullptr) {
                WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
                exit(EXIT_FAILURE);
            }
            syh->fHistUp_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistUp_orig  ->GetName())));
            syh->fHistDown_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistDown_orig->GetName())));
            Systematic* tmpsyst = new Systematic(*(sh->fSyst[i_syst]->fSystematic));
            // want to inherit the triggering systematic, to follow one (and -only one-) convention:
            tmpsyst->fName = NuisParName;
            tmpsyst->fStoredName = NuisParName;
            if (tmpsyst->fType == Systematic::OVERALL ) {
                  tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                  tmpsyst->fIsNormOnly = false;
            }
            syh->fSystematic = tmpsyst;
            fSample->AddSystematic(syh->fSystematic);
        }
    }
}

//_____________________________________________________________________________
//
void SampleHist::Multiply(SampleHist *sh){
    std::unique_ptr<TH1> hOrig( static_cast<TH1*>(fHist->Clone("h_tmp_orig")));
    if (sh->fHist!=nullptr) fHist->Multiply( sh->fHist.get() );
    else  {
       if (TRExFitter::HISTOCHECKCRASH) {
            WriteErrorStatus("SampleHist::Multiply", "Sample "+sh->fName+ " not found when trying to multiply it to "+fSample->fName);
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("SampleHist::Multiply", "Sample "+sh->fName+ " not found when trying to multiply it to "+fSample->fName);
        }
    }

    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        const std::string systName = fSyst[i_syst]->fName;
        const std::string NuisParName = fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = sh->GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Multiply", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Multiply", "Using its nominal. ");
            fSyst[i_syst]->Multiply( sh->fHist.get() );
        }
        else{
            WriteDebugStatus("SampleHist::Multiply", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Multiply", "Properly computing with that. ");
            fSyst[i_syst]->Multiply( syh );
        }
    }
    // loop on all the systematics in the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        const std::string systName = sh->fSyst[i_syst]->fName;
        const std::string NuisParName = sh->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Multiply", "Adding syst "+ NuisParName + " (through syst "+ systName + ") to sample "+ fName);
            std::unique_ptr<TH1> hUp(static_cast<TH1*>(hOrig->Clone("h_tmp_up")));
            std::unique_ptr<TH1> hDown(static_cast<TH1*>(hOrig->Clone("h_tmp_down")));
            hUp  ->Multiply(sh->fSyst[i_syst]->fHistUp.get());
            hDown->Multiply(sh->fSyst[i_syst]->fHistDown.get());
            syh = AddHistoSyst(NuisParName,NuisParName,hUp.get(),hDown.get());
            if (syh == nullptr) {
                WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
                exit(EXIT_FAILURE);
            }
            syh->fHistUp_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistUp_orig  ->GetName())));
            syh->fHistDown_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistDown_orig->GetName())));
            Systematic* tmpsyst = new Systematic(*(sh->fSyst[i_syst]->fSystematic));
            // want to inherit the triggering systematic, to follow one (and -only one-) convention:
            tmpsyst->fName = NuisParName;
            tmpsyst->fStoredName = NuisParName;
            if (tmpsyst->fType == Systematic::OVERALL ) {
                  tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                  tmpsyst->fIsNormOnly = false;
            }
            syh->fSystematic = tmpsyst;
            fSample->AddSystematic(syh->fSystematic);
        }
    }
}

//_____________________________________________________________________________
//
void SampleHist::Add(SampleHist *sh,double scale){
    std::unique_ptr<TH1> hOrig(static_cast<TH1*>(fHist->Clone("h_tmp_orig")));
    if (sh->fHist != nullptr) fHist->Add( sh->fHist.get(), scale );
    else  {
       if (TRExFitter::HISTOCHECKCRASH) {
            WriteErrorStatus("SampleHist::Add", "Sample "+sh->fName+ " not found when trying to add it to "+fSample->fName);
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("SampleHist::Add", "Sample "+sh->fName+ " not found when trying to add it to "+fSample->fName);
        }
    }

    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        const std::string systName = fSyst[i_syst]->fName;
        const std::string NuisParName = fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = sh->GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Add", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Add", "Using its nominal. ");
            fSyst[i_syst]->Add( sh->fHist.get(), scale );
        }
        else{
            WriteDebugStatus("SampleHist::Add", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Add", "Properly computing with that. ");
            fSyst[i_syst]->Add( syh, scale );
        }
    }
    // loop on all the systematics of the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        const std::string systName = sh->fSyst[i_syst]->fName;
        const std::string NuisParName = sh->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Add", "Adding syst "+ NuisParName + " (through syst "+ systName + ") to sample "+ fName);
            std::unique_ptr<TH1> hUp(static_cast<TH1*>(hOrig->Clone("h_tmp_up")));
            std::unique_ptr<TH1> hDown(static_cast<TH1*>(hOrig->Clone("h_tmp_down")));
            if (sh->fSyst[i_syst]->fHistUp == nullptr) {
               if (TRExFitter::HISTOCHECKCRASH) {
                    WriteErrorStatus("SampleHist::Add", "Systematic "+sh->fSyst[i_syst]->fName+ " up var. not found when trying to adding it to "+fSample->fName);
                    exit(EXIT_FAILURE);
               } else {
                    WriteWarningStatus("SampleHist::Add", "Systematic "+sh->fSyst[i_syst]->fName+ " up var. not found when trying to adding it to "+fSample->fName);
               }
            }
            else  hUp  ->Add( sh->fSyst[i_syst]->fHistUp.get(), scale);
            if (sh->fSyst[i_syst]->fHistDown == nullptr) {
               if (TRExFitter::HISTOCHECKCRASH) {
                    WriteErrorStatus("SampleHist::Add", "Systematic "+sh->fSyst[i_syst]->fName+ " down var. not found when trying to adding it to "+fSample->fName);
                    exit(EXIT_FAILURE);
               } else {
                    WriteWarningStatus("SampleHist::Add", "Systematic "+sh->fSyst[i_syst]->fName+ " down var. not found when trying to adding it to "+fSample->fName);
               }
            }
            else hDown->Add( sh->fSyst[i_syst]->fHistDown.get(),scale );
            syh = AddHistoSyst(NuisParName,NuisParName,hUp.get(),hDown.get());
            if (syh == nullptr) {
                WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
                exit(EXIT_FAILURE);
            }
            syh->fHistUp_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistUp_orig  ->GetName())));
            syh->fHistDown_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistDown_orig->GetName())));
            Systematic* tmpsyst = new Systematic(*(sh->fSyst[i_syst]->fSystematic));
            // want to inherit the triggering systematic, to follow one (and -only one-) convention:
            tmpsyst->fName = NuisParName;
            tmpsyst->fStoredName = NuisParName;
            if (tmpsyst->fType == Systematic::OVERALL ) {
                  tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                  tmpsyst->fIsNormOnly = false;
            }
            syh->fSystematic = tmpsyst;
            fSample->AddSystematic(syh->fSystematic);
        }
    }
}

//_____________________________________________________________________________
//
void SampleHist::Scale(double scale){
    fHist->Scale( scale );
    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        fSyst[i_syst]->fHistUp->Scale( scale );
        if(fSyst[i_syst]->fHistShapeUp!=nullptr)   fSyst[i_syst]->fHistShapeUp->Scale( scale );
        fSyst[i_syst]->fHistDown->Scale( scale );
        if(fSyst[i_syst]->fHistShapeDown!=nullptr) fSyst[i_syst]->fHistShapeDown->Scale( scale );
    }
}

//_____________________________________________________________________________
//
void SampleHist::SystPruning(PruningUtil *pu,TH1* hTot){
    for(auto& syh : fSyst){
        if(!syh) continue;
        if(!syh->fHistUp) continue;
        if(!syh->fHistDown) continue;
        const int pruningResult = pu->CheckSystPruning(syh->fHistUp.get(), syh->fHistDown.get(), fHist.get(), hTot);
        syh->fBadShape = (pruningResult==-3 || pruningResult==-4);
        syh->fBadNorm = (pruningResult==-2 || pruningResult==-4);
        syh->fShapePruned = (pruningResult==1 || pruningResult==3 || syh->fBadShape);
        syh->fNormPruned = (pruningResult==2 || pruningResult==3 || syh->fBadNorm);
    }
}

//_____________________________________________________________________________
//
void SampleHist::DrawSystPlotUpper(TPad* pad0,
                                   TH1* nominal,
                                   TH1* nominal_orig,
                                   TH1* syst_up,
                                   TH1* syst_up_orig,
                                   TH1* syst_down,
                                   TH1* syst_down_orig,
                                   TH1* data,
                                   TH1* tmp,
                                   bool SumAndData,
                                   bool bothPanels) const {
    pad0->cd();
    nominal->SetLineStyle(1);
    double max = std::max({std::fabs(nominal->GetMaximum()),
                          std::fabs(nominal_orig->GetMaximum()),
                          std::fabs(syst_up->GetMaximum()),
                          std::fabs(syst_up_orig->GetMaximum()),
                          std::fabs(syst_down->GetMaximum()),
                          std::fabs(syst_down_orig->GetMaximum()),
                          std::fabs(syst_up->GetMinimum()),
                          std::fabs(syst_up_orig->GetMinimum()),
                          std::fabs(syst_down->GetMinimum()),
                          std::fabs(syst_down_orig->GetMinimum())});

    if (SumAndData) max = std::max({max, std::fabs(data->GetMaximum()), std::fabs(data->GetMinimum())});

    tmp->GetYaxis()->SetTitle("Number of events");
    tmp->SetMinimum(1e-5);
    tmp->SetMaximum(max * 2.0);
    tmp->GetXaxis()->SetTitle(fVariableTitle.c_str());
    tmp->Draw("HIST");

    if(TRExFitter::SYSTERRORBARS){
        syst_down_orig->SetMarkerSize(0);
        syst_up_orig->SetMarkerSize(0);
        syst_down_orig->DrawCopy("same E");
        syst_up_orig->DrawCopy("same E");
    } else {
        syst_down_orig->DrawCopy("same HIST");
        syst_up_orig->DrawCopy("same HIST");
    }
    syst_down->DrawCopy("same HIST");
    syst_up->DrawCopy("same HIST");
    nominal->DrawCopy("same HIST");
    nominal->SetFillStyle(3005);
    nominal->SetFillColor(kBlue);
    nominal->SetMarkerSize(0);
    nominal->DrawCopy("e2same");
    nominal_orig->DrawCopy("same HIST");
    if(bothPanels && SumAndData) data->Draw("EX0same");
}

//_____________________________________________________________________________
//
void SampleHist::DrawSystPlotRatio(TPad* pad1,
                                   TH1* nominal,
                                   TH1* nominal_orig,
                                   TH1* syst_up,
                                   TH1* syst_up_orig,
                                   TH1* syst_down,
                                   TH1* syst_down_orig,
                                   TH1* data,
                                   TH1* tmp,
                                   bool SumAndData) const {

    pad1->cd();

    nominal->SetLineStyle(2);

    syst_up->Add(nominal, -1);
    syst_down->Add(nominal, -1);
    if (SumAndData) data->Add(nominal, -1);
    syst_up->Divide(nominal);
    syst_down->Divide(nominal);
    if (SumAndData) data->Divide(nominal);
    for(int i_bin=1; i_bin <= nominal->GetNbinsX(); ++i_bin){
        if(nominal->GetBinContent(i_bin)<1e-5){
            syst_up  ->SetBinContent(i_bin,0.);
            syst_down->SetBinContent(i_bin,0.);
            if (SumAndData) data->SetBinContent(i_bin,0.);
        }
    }
    syst_up->Scale(100);
    syst_down->Scale(100);
    if(SumAndData) data->Scale(100);
    
    syst_up_orig->Add(nominal_orig, -1);
    syst_down_orig->Add(nominal_orig, -1);
    syst_up_orig->Divide(nominal_orig);
    syst_down_orig->Divide(nominal_orig);
    for(int i_bin=1; i_bin <= nominal_orig->GetNbinsX(); ++i_bin){
        if(nominal_orig->GetBinContent(i_bin)<1e-5){
            syst_up_orig  ->SetBinContent(i_bin,0.);
            syst_down_orig->SetBinContent(i_bin,0.);
        }
    }
    syst_up_orig->Scale(100);
    syst_down_orig->Scale(100);

    tmp->GetYaxis()->SetTitle("#frac{Syst.-Nom.}{Nom.} [%]");
    tmp->GetYaxis()->SetTitleOffset(1.6);
    tmp->GetXaxis()->SetTitleOffset(3.);

    double max = std::max({std::fabs(syst_up->GetMaximum()),
                          std::fabs(syst_up_orig->GetMaximum()),
                          std::fabs(syst_down->GetMaximum()),
                          std::fabs(syst_down_orig->GetMaximum()),
                          std::fabs(syst_up->GetMinimum()),
                          std::fabs(syst_up_orig->GetMinimum()),
                          std::fabs(syst_down->GetMinimum()),
                          std::fabs(syst_down_orig->GetMinimum())});

    if (SumAndData) max = std::max({max, std::fabs(data->GetMaximum()), std::fabs(data->GetMinimum())});
    
    if(TRExFitter::OPTION["SystPlotRatioRange"]!=0){
        tmp->SetMinimum(-TRExFitter::OPTION["SystPlotRatioRange"]);
        tmp->SetMaximum( TRExFitter::OPTION["SystPlotRatioRange"]);
    }
    else{
        tmp->SetMinimum(-max*1.5);
        tmp->SetMaximum( max*1.5);
    }

    tmp->GetXaxis()->SetTitle(fVariableTitle.c_str());
    tmp->Draw("HIST");
    
    if(TRExFitter::SYSTERRORBARS){
        syst_down_orig->SetMarkerSize(0);
        syst_up_orig->SetMarkerSize(0);
        syst_down_orig->DrawCopy("same E");
        syst_up_orig->DrawCopy("same E");
    }
    else{
        syst_down_orig->DrawCopy("same HIST");
        syst_up_orig->DrawCopy("same HIST");
    }
    syst_down->DrawCopy("same HIST");
    syst_up->DrawCopy("same HIST");
    nominal->SetFillStyle(3005);
    nominal->SetFillColor(kBlue);
    nominal->SetMarkerSize(0);
    for (int i=1; i <= nominal->GetNbinsX(); ++i) {
        nominal->SetBinError(i, nominal->GetBinError(i)*100. / nominal->GetBinContent(i));
        nominal->SetBinContent(i,0);
    }
    nominal->DrawCopy("e2same");
    if(SumAndData) data->Draw("EX0same");
}
