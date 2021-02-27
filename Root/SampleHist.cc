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
#include "TFile.h"
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
    fHist->SetDirectory(nullptr);
    fHist_orig->SetDirectory(nullptr);
    
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
    fHist->SetDirectory(nullptr);
    fHist_orig->SetDirectory(nullptr);

    fIsMorph = fSample->fIsMorph;
}

//_____________________________________________________________________________
//
SampleHist::~SampleHist(){
}

//_____________________________________________________________________________
//
std::shared_ptr<SystematicHist> SampleHist::AddOverallSyst(const std::string& name,const std::string& storedName,double up,double down){
    std::shared_ptr<SystematicHist> syh = GetSystematic(name);
    // ... and if not create a new one
    if(!syh){
        fSyst.emplace_back(new SystematicHist(name));
        syh = fSyst.back();
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
    syh->fHistUp->SetDirectory(nullptr);
    syh->fHistDown->SetDirectory(nullptr);
    syh->fHistUp_orig->SetDirectory(nullptr);
    syh->fHistDown_orig->SetDirectory(nullptr);
    syh->fIsOverall = true;
    syh->fIsShape   = false;
    syh->fNormUp   = up;
    syh->fNormDown = down;
    return syh;
}

//_____________________________________________________________________________
//
std::shared_ptr<SystematicHist> SampleHist::AddStatSyst(const std::string& name,const std::string& storedName, int i_bin) {
    int bin = i_bin+1; // counting of bins in Root starts with 1, in TRExFitter with 0
    std::shared_ptr<SystematicHist> syh = GetSystematic(name);
    // ... and if not create a new one
    if(!syh){
        fSyst.emplace_back(new SystematicHist(name));
        syh = fSyst.back();
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
    syh->fHistUp->SetDirectory(nullptr);
    syh->fHistDown->SetDirectory(nullptr);
    syh->fHistShapeUp->SetDirectory(nullptr);
    syh->fHistShapeDown->SetDirectory(nullptr);
    syh->fHistUp_orig->SetDirectory(nullptr);
    syh->fHistDown_orig->SetDirectory(nullptr);
    return syh;
}

//_____________________________________________________________________________
//
std::shared_ptr<SystematicHist> SampleHist::AddHistoSyst(const std::string& name,const std::string& storedName,TH1* h_up,TH1* h_down){
    
    // before doing anything else, check if the sampleHist can be created
    if(h_up  ==nullptr) return nullptr;
    if(h_down==nullptr) return nullptr;

    std::shared_ptr<SystematicHist> syh = GetSystematic(name);
    // ... and if not create a new one
    if(!syh){
        fSyst.emplace_back(new SystematicHist(name));
        syh = fSyst.back();
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

    syh->fHistUp->SetDirectory(nullptr);
    syh->fHistDown->SetDirectory(nullptr);
    syh->fHistShapeUp->SetDirectory(nullptr);
    syh->fHistShapeDown->SetDirectory(nullptr);
    syh->fHistUp_orig->SetDirectory(nullptr);
    syh->fHistDown_orig->SetDirectory(nullptr);
    syh->fHistUp_preSmooth->SetDirectory(nullptr);
    syh->fHistDown_preSmooth->SetDirectory(nullptr);
    syh->fIsOverall = true;
    syh->fIsShape   = true;
    syh->fNormUp   = ( Common::EffIntegral(syh->fHistUp.get())   -  Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
    syh->fNormDown = ( Common::EffIntegral(syh->fHistDown.get()) - Common::EffIntegral(fHist.get()) ) / Common::EffIntegral(fHist.get());
    if(syh->fNormUp == 0 && syh->fNormDown == 0) syh->fIsOverall = false;
    return syh;
}

//_____________________________________________________________________________
//
std::shared_ptr<SystematicHist> SampleHist::AddHistoSyst(const std::string& name,
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

    std::shared_ptr<SystematicHist> sh = GetSystematic(name);
    // ... and if not create a new one
    if(!sh){
        fSyst.emplace_back(new SystematicHist(name));
        sh = fSyst.back();
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
    sh->fHistUp->SetDirectory(nullptr);
    sh->fHistDown->SetDirectory(nullptr);
    sh->fHistShapeUp->SetDirectory(nullptr);
    sh->fHistShapeDown->SetDirectory(nullptr);
    sh->fHistUp_orig->SetDirectory(nullptr);
    sh->fHistDown_orig->SetDirectory(nullptr);
    //
    return sh;
}

//_____________________________________________________________________________
//
std::shared_ptr<NormFactor> SampleHist::AddNormFactor(std::shared_ptr<NormFactor> normFactor) {
    std::shared_ptr<NormFactor> norm = GetNormFactor(normFactor->fName);
    if(!norm){
        fSample->fNormFactors.emplace_back(normFactor);
        norm = fSample->fNormFactors.back();
    }
    else{
        norm = normFactor;
    }
    return norm;
}

//_____________________________________________________________________________
//
std::shared_ptr<NormFactor> SampleHist::AddNormFactor(const std::string& name,double nominal, double min, double max){
    std::shared_ptr<NormFactor> norm = GetNormFactor(name);
    if(!norm){
        fSample->fNormFactors.emplace_back(new NormFactor(name,nominal,min,max));
        norm = fSample->fNormFactors.back();
    }
    return norm;
}

//_____________________________________________________________________________
//
std::shared_ptr<ShapeFactor> SampleHist::AddShapeFactor(std::shared_ptr<ShapeFactor> shapeFactor){
    std::shared_ptr<ShapeFactor> shape = GetShapeFactor(shapeFactor->fName);
    if(!shape){
        fSample->fShapeFactors.emplace_back(shapeFactor);
        shape = fSample->fShapeFactors.back();
    } else {
        shape = shapeFactor;
    }
    return shape;
}

//_____________________________________________________________________________
//
std::shared_ptr<ShapeFactor> SampleHist::AddShapeFactor(const std::string& name,double nominal, double min, double max){
    std::shared_ptr<ShapeFactor> shape = GetShapeFactor(name);
    if(!shape){
        fSample->fShapeFactors.emplace_back(new ShapeFactor(name,nominal,min,max));
        shape = fSample->fShapeFactors.back();
    }
    return shape;
}

//_____________________________________________________________________________
//
std::shared_ptr<SystematicHist> SampleHist::GetSystematic(const std::string& systName) const{
    for(const auto& isyst : fSyst) {
        if(systName == isyst->fName) return isyst;
    }
    return nullptr;
}

//_____________________________________________________________________________
//
std::shared_ptr<SystematicHist> SampleHist::GetSystFromNP(const std::string& NuisParName) const{
    for(const auto& isyst : fSyst) {
        if(NuisParName == isyst->fSystematic->fNuisanceParameter) return isyst;
    }
    return nullptr;
}

//_____________________________________________________________________________
//
std::shared_ptr<NormFactor> SampleHist::GetNormFactor(const std::string& name) const{
    if (!fSample) {
        WriteErrorStatus("SampleHist::GetNormFactor", "Nullptr for fSample!");
        exit(EXIT_FAILURE);
    }

    for(const auto& inorm : fSample->fNormFactors) {
        if(name == inorm->fName) return inorm;
    }
    return nullptr;
}

//_____________________________________________________________________________
//
std::shared_ptr<ShapeFactor> SampleHist::GetShapeFactor(const std::string& name) const{
    if (!fSample) {
        WriteErrorStatus("SampleHist::GetShapeFactor", "Nullptr for fSample!");
        exit(EXIT_FAILURE);
    }
    for(const auto& ishape : fSample->fShapeFactors) {
        if(name == ishape->fName) return ishape;
    }
    return nullptr;
}

//_____________________________________________________________________________
//
bool SampleHist::HasSyst(const std::string& name) const{
    for(const auto& isyst : fSyst) {
        if(isyst->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasNorm(const std::string& name) const{
    for(const auto& inorm : fSample->fNormFactors) {
        if(inorm->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasShapeFactor(const std::string& name) const{
    for(const auto& ishape : fSample->fShapeFactors) {
        if(ishape->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
void SampleHist::WriteToFile(const std::vector<int>& blindedBins, std::shared_ptr<TFile> f,bool reWriteOrig){
    if(f==nullptr){
        if(fHist_orig!=nullptr && reWriteOrig)   Common::WriteHistToFile(fHist_orig.get(), fFileName);
        if(fHist!=nullptr)        Common::WriteHistToFile(fHist.get(), fFileName);
    }
    else{
        if(fHist_orig!=nullptr && reWriteOrig)   Common::WriteHistToFile(fHist_orig.get(), f);
        if(fHist!=nullptr)        Common::WriteHistToFile(fHist.get(), f);
    }
    // create the regular binning histogram
    fHist_regBin = HistoTools::TranformHistogramBinning(fHist.get(), blindedBins);
    if(fHist_regBin!=nullptr) Common::WriteHistToFile(fHist_regBin.get(),f);
    //
    // save separate gammas as histograms
    if(fSample->fSeparateGammas){
        std::unique_ptr<TH1> htempUp(static_cast<TH1*>(fHist->Clone()));
        std::unique_ptr<TH1> htempDown(static_cast<TH1*>(fHist->Clone()));
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
        std::shared_ptr<Systematic> gamma = nullptr;
        if(GetSystematic(systName)) gamma = GetSystematic(systName)->fSystematic;  //GetSystematic(systName);
        if(gamma==nullptr) gamma = std::make_shared<Systematic>(systName,Systematic::SHAPE);
        WriteDebugStatus("SampleHist::WriteToFile", "adding separate gammas as SHAPE systematic " + systName);
        gamma->fRegions.clear();
        gamma->fRegions.push_back(fRegionName);
        std::shared_ptr<SystematicHist> syh = AddHistoSyst(systName,systName,htempUp.get(),htempDown.get());
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
    for(auto& isyst : fSyst) {
        // make sure they all have the correct name!
        isyst->fHistUp  ->SetName( Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),isyst->fSystematic->fStoredName.c_str()) );
        isyst->fHistDown->SetName( Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),isyst->fSystematic->fStoredName.c_str()) );
        isyst->fHistUp_orig  ->SetName( Form("%s_%s_%s_Up_orig",  fRegionName.c_str(),fSample->fName.c_str(),isyst->fSystematic->fStoredName.c_str()) );
        isyst->fHistDown_orig->SetName( Form("%s_%s_%s_Down_orig",fRegionName.c_str(),fSample->fName.c_str(),isyst->fSystematic->fStoredName.c_str()) );
        if(f==nullptr) isyst->WriteToFile(blindedBins, nullptr,reWriteOrig);
        else           isyst->WriteToFile(blindedBins, f,      reWriteOrig);
        // for shape hist, save also the syst(up)-nominal (to feed HistFactory)
        if(isyst->fSystematic->fType==Systematic::SHAPE){
            std::unique_ptr<TH1> tmp(static_cast<TH1*>(isyst->fHistUp->Clone(Form("%s_%s_%s_Up_Var",  fRegionName.c_str(),fSample->fName.c_str(),isyst->fSystematic->fStoredName.c_str()))));
            std::unique_ptr<TH1> hVar = HistoTools::TranformHistogramBinning(tmp.get(), blindedBins);
            hVar->Add(fHist_regBin.get(),-1);
            hVar->Divide(fHist_regBin.get());
            // no negative bins here!
            for(int i_bin=1;i_bin<=hVar->GetNbinsX();i_bin++){
                if(hVar->GetBinContent(i_bin)<0) hVar->SetBinContent(i_bin,-1.*hVar->GetBinContent(i_bin));
            }
            if(f==nullptr) Common::WriteHistToFile(hVar.get(),fFileName);
            else           Common::WriteHistToFile(hVar.get(),f);
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
    for(const auto& isyst : fSyst) {
        if(isyst==nullptr) continue;
        if(isyst->fHistUp  ==nullptr) continue;
        if(isyst->fHistDown==nullptr) continue;
        const double tmpYieldUp   = isyst->fHistUp->Integral();
        const double tmpYieldDown = isyst->fHistDown->Integral();
        yieldUp.push_back(tmpYieldUp);
        yieldDown.push_back(tmpYieldDown);
        // warnings if total yield in systematic variations is negative, in which case normalization cannot be preserved
        if (tmpYieldUp  <0){
            NegativeTotalYieldWarning(isyst->fHistUp.get(),   tmpYieldUp);
        }
        if (tmpYieldDown<0){
            NegativeTotalYieldWarning(isyst->fHistDown.get(), tmpYieldDown);
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
            for(auto& isyst : fSyst) {
                if(isyst->fHistUp  ->GetBinContent(i_bin)<=0) isyst->fHistUp  ->SetBinContent(i_bin,1e-06);
                if(isyst->fHistDown->GetBinContent(i_bin)<=0) isyst->fHistDown->SetBinContent(i_bin,1e-06);
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
    for(std::size_t i_syst = 0; i_syst < fSyst.size(); ++i_syst) {
        std::shared_ptr<SystematicHist> syh = fSyst[i_syst];
        if(syh->fHistUp  ->Integral()!=yieldUp[i_syst]  ) syh->fHistUp  ->Scale(yieldUp[i_syst]  /syh->fHistUp  ->Integral());
        if(syh->fHistDown->Integral()!=yieldDown[i_syst]) syh->fHistDown->Scale(yieldDown[i_syst]/syh->fHistDown->Integral());
    }
}

//_____________________________________________________________________________
//
void SampleHist::Print() const{
    std::string temp = fHist->GetName();
    WriteDebugStatus("SampleHist::Print", "      Sample: " + fName + "\t" + temp);
    if(fSyst.size() >0){
        temp = "        Systematics:   ";
        for(const auto& isyst : fSyst) {
            temp+= " " + isyst->fName;
        }
        WriteDebugStatus("SampleHist::Print", temp);
    }
    if(fSample->fNormFactors.size() > 0) {
        temp = "        NormFactor(s): ";
        for(const auto& inorm : fSample->fNormFactors) {
            temp+= " " + inorm->fName;
        }
        WriteDebugStatus("SampleHist::Print", temp);
    }
    if(fSample->fShapeFactors.size()> 0) {
        temp = "        ShapeFactor(s): ";
        for(const auto& ishape : fSample-> fShapeFactors) {
            temp+= " " + ishape->fName;
        }
        WriteDebugStatus("SampleHist::Print", temp);
    }
}

//_____________________________________________________________________________
//
void SampleHist::Rebin(int ngroup, const Double_t* xbins){
    fHist->Rebin(ngroup,"",xbins);
    for(auto& isyst : fSyst) {
        if(isyst->fHistUp!=nullptr)        isyst->fHistUp->Rebin(ngroup,"",xbins);
        if(isyst->fHistDown!=nullptr)      isyst->fHistDown->Rebin(ngroup,"",xbins);
        if(isyst->fHistShapeUp!=nullptr)   isyst->fHistShapeUp->Rebin(ngroup,"",xbins);
        if(isyst->fHistShapeDown!=nullptr) isyst->fHistShapeDown->Rebin(ngroup,"",xbins);
    }
}

//_____________________________________________________________________________
// this draws the control plots (for each systematic) with the syst variations for this region & all sample
void SampleHist::DrawSystPlot( const string &syst, TH1* const h_data, bool SumAndData, bool bothPanels ) const{
    for(auto& isyst : fSyst) {
        if(syst!="all" && isyst->fName.find(syst)==string::npos) continue;
        
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
        std::unique_ptr<TH1> syst_up(static_cast<TH1*>(isyst->fHistUp->Clone()));
        std::unique_ptr<TH1> syst_up_orig(static_cast<TH1*>(isyst->fHistUp_preSmooth->Clone()));
        std::unique_ptr<TH1> syst_down(static_cast<TH1*>(isyst->fHistDown->Clone()));
        std::unique_ptr<TH1> syst_down_orig(static_cast<TH1*>(isyst->fHistDown_preSmooth->Clone()));
        std::unique_ptr<TH1> data(nullptr);
        if (SumAndData && h_data) data = std::unique_ptr<TH1>(static_cast<TH1*>(h_data->Clone("nominal")));
        std::unique_ptr<TH1> tmp(static_cast<TH1*>(nominal->Clone()));

        // drop shape or norm (for cases where this is not yet done in the stored histogrmas, i.e. in case of pruning or decorrelation)
        if (isyst != nullptr && isyst->fSystematic != nullptr) {
            if(isyst->fSystematic->fIsNormOnly){
                Common::DropShape(syst_up.get(),syst_down.get(),nominal.get());
            }
            if(isyst->fSystematic->fIsShapeOnly){
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
        if (SumAndData&&h_data) data->SetMarkerColor(kBlack);

        // make copies for ratio
        std::unique_ptr<TH1> nominal_ratio(static_cast<TH1*>(nominal->Clone()));      
        std::unique_ptr<TH1> nominal_orig_ratio(static_cast<TH1*>(nominal_orig->Clone()));      
        std::unique_ptr<TH1> syst_up_ratio(static_cast<TH1*>(syst_up->Clone()));      
        std::unique_ptr<TH1> syst_up_orig_ratio(static_cast<TH1*>(syst_up_orig->Clone()));      
        std::unique_ptr<TH1> syst_down_ratio(static_cast<TH1*>(syst_down->Clone()));      
        std::unique_ptr<TH1> syst_down_orig_ratio(static_cast<TH1*>(syst_down_orig->Clone()));      
        std::unique_ptr<TH1> data_ratio(nullptr);
        if (SumAndData&&h_data) data_ratio = std::unique_ptr<TH1>(static_cast<TH1*>(data->Clone()));      
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
                          SumAndData&&h_data,
                          bothPanels);

        // Draw laels
        TLatex tex{};
        tex.SetNDC();
        if(SumAndData) {
            if(isyst->fSystematic) tex.DrawLatex(0.17,0.79,Form("%s",isyst->fSystematic->fTitle.c_str()));
            else                   tex.DrawLatex(0.17,0.79,Form("%s",isyst->fName.c_str()));
        } else{
            if(isyst->fSystematic) tex.DrawLatex(0.17,0.79,Form("%s, %s",isyst->fSystematic->fTitle.c_str(),fSample->fTitle.c_str()));
            else                   tex.DrawLatex(0.17,0.79,Form("%s, %s",isyst->fName.c_str(),fSample->fTitle.c_str()));
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
        if (SumAndData&&h_data) yield_data = Common::CorrectIntegral(data.get());
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
        if(SumAndData&&h_data){
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
                          SumAndData&&h_data);

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
        gSystem->mkdir((fFitName+"/Systematics/"+isyst->fName).c_str());

        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            if(SumAndData) {
                c.SaveAs(Form("%s/Systematics/%s/%s_%s.%s",fFitName.c_str(),isyst->fName.c_str(), fName.c_str(), isyst->fName.c_str(), format.c_str()));
            } else { 
                c.SaveAs(Form("%s/Systematics/%s/%s_%s.%s",fFitName.c_str(),isyst->fName.c_str(),fHist->GetName(), isyst->fName.c_str(), format.c_str()));
            }
        }

    }
}

//_____________________________________________________________________________
//
void SampleHist::SmoothSyst(const HistoTools::SmoothOption &smoothOpt, const bool useAlternativeShapeHistFactory, string syst, bool force){
    if(fSystSmoothed && !force) return;
    std::unique_ptr<TH1> h_nominal(static_cast<TH1*>(fHist->Clone("h_nominal")));
    if (!h_nominal->GetSumw2()) h_nominal->Sumw2();

    for(auto& isyst : fSyst) {

        if(syst!="all" && isyst->fName != syst) continue;

        if(isyst->fHistUp  ==nullptr) continue;
        if(isyst->fHistDown==nullptr) continue;
        if(isyst->fSystematic==nullptr) continue;

        TH1* h_syst_up = nullptr;
        TH1* h_syst_down = nullptr;

        if(isyst->fSmoothType + isyst->fSymmetrisationType<=0){
            HistoTools::Scale(isyst->fHistUp.get(),   fHist.get(),isyst->fScaleUp);
            HistoTools::Scale(isyst->fHistDown.get(), fHist.get(),isyst->fScaleDown);
            continue;
        }

        TH1* shape_up   = nullptr;
        TH1* shape_down = nullptr;
        std::unique_ptr<TH1> shapeOriginUp  (static_cast<TH1*>(isyst->fHistUp->Clone()));
        std::unique_ptr<TH1> shapeOriginDown(static_cast<TH1*>(isyst->fHistDown->Clone()));
        const double nominal = h_nominal->Integral();
        if (std::fabs(shapeOriginUp->Integral()) > 1e-6) {
            shapeOriginUp  ->Scale(nominal/shapeOriginUp->Integral());
        }
        if (std::fabs(shapeOriginDown->Integral()) > 1e-6) {
            shapeOriginDown->Scale(nominal/shapeOriginDown->Integral());
        }

        //
        // Pre-smoothing
        //
        // (do smoothing and symmetrization before pre-smoothing in case of two-sided systematics)
        if(isyst->fSymmetrisationType==HistoTools::SYMMETRIZETWOSIDED){
            if(isyst->fIsShape){
                if(isyst->fSystematic->fSampleSmoothing){
                    HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                    isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal.get(),//nominal histogram
                                                    isyst->fHistUp.get(),
                                                    isyst->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    isyst->fScaleUp,
                                                    isyst->fScaleDown, // scale factors
                                                    isyst->fSystematic->fSampleSmoothOption // overwrite smoothing option
                                                );

                    if (useAlternativeShapeHistFactory) {
                        HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                        isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                        h_nominal.get(),//nominal histogram
                                                        shapeOriginUp.get(),
                                                        shapeOriginDown.get(),//original histograms
                                                        shape_up, shape_down, //modified histograms
                                                        1.0,
                                                        1.0, // scale factors
                                                        isyst->fSystematic->fSampleSmoothOption // overwrite smoothing option
                                                    );
                    }
                }
                else{
                    HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                    isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal.get(),//nominal histogram
                                                    isyst->fHistUp.get(),
                                                    isyst->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    isyst->fScaleUp,
                                                    isyst->fScaleDown, // scale factors
                                                    smoothOpt
                                                );
                    if (useAlternativeShapeHistFactory) {
                        HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                        isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                        h_nominal.get(),//nominal histogram
                                                        shapeOriginUp.get(),
                                                        shapeOriginDown.get(),//original histograms
                                                        shape_up, shape_down, //modified histograms
                                                        1.0,
                                                        1.0, // scale factors
                                                        smoothOpt
                                                    );
                    }
                }
            }
            //
            // need to ad these lines to make sure overall only systematics get scaled as well
            else{
                HistoTools::Scale(isyst->fHistUp.get(),  fHist.get(),isyst->fScaleUp);
                HistoTools::Scale(isyst->fHistDown.get(),fHist.get(),isyst->fScaleDown);
            }
        }

        if(isyst->fSystematic->fPreSmoothing){
            std::unique_ptr<TH1> h_tmp_up(nullptr);
            std::unique_ptr<TH1> h_tmp_down(nullptr);
            if (h_syst_up) {
                h_tmp_up.reset(static_cast<TH1*>(h_syst_up->Clone()));
                h_tmp_down.reset(static_cast<TH1*>(h_syst_down->Clone()));
            }
            if(h_tmp_up!=nullptr || h_tmp_down!=nullptr){
                std::unique_ptr<TH1> h_tmp_nominal(static_cast<TH1*>(h_nominal->Clone()));
                for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++){
                    h_tmp_nominal->GetBinError(i_bin,0);
                }
                if(h_tmp_up!=nullptr){
                    const double tmp_nom_up = h_tmp_up->Integral();
                    h_tmp_up->Add(h_tmp_nominal.get(),-1);
                    h_tmp_up->Divide(h_tmp_nominal.get());
                    for(int i_bin = 1; i_bin <= h_tmp_nominal->GetNbinsX(); ++i_bin) {
                        h_tmp_up->AddBinContent(i_bin, 100.);
                    }
                    h_tmp_up->Smooth();
                    for(int i_bin = 1; i_bin <= h_tmp_nominal->GetNbinsX(); ++i_bin) {
                        h_tmp_up->AddBinContent(i_bin,-100.);
                    }
                    h_tmp_up->Multiply(h_tmp_nominal.get());
                    h_tmp_up->Add(h_tmp_nominal.get(), 1);
                    h_tmp_up->Scale(tmp_nom_up/h_tmp_up->Integral());
                    delete h_syst_up;
                    h_syst_up = static_cast<TH1*>(h_tmp_up->Clone());
                }
                if(h_tmp_down!=nullptr){
                    const double tmp_nom_down = h_tmp_down->Integral();
                    h_tmp_down->Add(h_tmp_nominal.get(),-1);
                    h_tmp_up->Divide(h_tmp_nominal.get());
                    for(int i_bin = 1; i_bin <= h_tmp_nominal->GetNbinsX(); ++i_bin) {
                        h_tmp_down->AddBinContent(i_bin, 100.);
                    }
                    h_tmp_down->Smooth();
                    for(int i_bin = 1; i_bin <= h_tmp_nominal->GetNbinsX(); ++i_bin) {
                        h_tmp_down->AddBinContent(i_bin,-100.);
                    }
                    h_tmp_up->Multiply(h_tmp_nominal.get());
                    h_tmp_down->Add(h_tmp_nominal.get(), 1);
                    h_tmp_down->Scale(tmp_nom_down/h_tmp_down->Integral());
                    delete h_syst_down;
                    h_syst_down = static_cast<TH1*>(h_tmp_down->Clone());
                }
            }
        }

        //
        // Call the function for smoothing and symmetrisation
        //
        if(isyst->fSymmetrisationType!=HistoTools::SYMMETRIZETWOSIDED){
            if(isyst->fIsShape){
                if(isyst->fSystematic->fSampleSmoothing){
                    HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                    isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal.get(),//nominal histogram
                                                    isyst->fHistUp.get(),
                                                    isyst->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    isyst->fScaleUp,
                                                    isyst->fScaleDown, // scale factors
                                                    isyst->fSystematic->fSampleSmoothOption // overwrite smoothing option
                                                );
                    if (useAlternativeShapeHistFactory) {
                        HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                        isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                        h_nominal.get(),//nominal histogram
                                                        shapeOriginUp.get(),
                                                        shapeOriginDown.get(),//original histograms
                                                        shape_up, shape_down, //modified histograms
                                                        1.0,
                                                        1.0, // scale factors
                                                        isyst->fSystematic->fSampleSmoothOption // overwrite smoothing option
                                                    );
                    }
                }
                else{
                    HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                    isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                    h_nominal.get(),//nominal histogram
                                                    isyst->fHistUp.get(),
                                                    isyst->fHistDown.get(),//original histograms
                                                    h_syst_up, h_syst_down, //modified histograms
                                                    isyst->fScaleUp,
                                                    isyst->fScaleDown, // scale factors
                                                    smoothOpt
                                                );
                    if (useAlternativeShapeHistFactory) {
                        HistoTools::ManageHistograms(   isyst->fSmoothType,
                                                        isyst->fSymmetrisationType,//parameters of the histogram massaging
                                                        h_nominal.get(),//nominal histogram
                                                        shapeOriginUp.get(),
                                                        shapeOriginDown.get(),//original histograms
                                                        shape_up, shape_down, //modified histograms
                                                        1.0,
                                                        1.0, // scale factors
                                                        smoothOpt
                                                    );
                    }
                }
            }
            //
            // need to ad these lines to make sure overall only systematics get scaled as well
            else{
                HistoTools::Scale(isyst->fHistUp.get(),  fHist.get(),isyst->fScaleUp);
                HistoTools::Scale(isyst->fHistDown.get(),fHist.get(),isyst->fScaleDown);
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
        isyst->fHistUp.reset(static_cast<TH1*>(h_syst_up->Clone(isyst->fHistUp->GetName())));
        isyst->fHistDown.reset(static_cast<TH1*>(h_syst_down->Clone(isyst->fHistDown->GetName())));

        //
        // Perform a check of the output histograms (check for 0 bins and other pathologic behaviour)
        //
        HistoTools::CheckHistograms( h_nominal.get() /*nominal*/, isyst.get() /*systematic*/, fSample -> fType != Sample::SIGNAL, TRExFitter::HISTOCHECKCRASH /*cause crash if problem*/);

        //
        // Normalisation component first
        //
        if(h_nominal->Integral()!=0){
            isyst->fNormUp   = isyst->fHistUp  ->Integral()/h_nominal->Integral() - 1.;
            isyst->fNormDown = isyst->fHistDown->Integral()/h_nominal->Integral() - 1.;
        } else {
            WriteErrorStatus("SampleHist::SmoothSyst", "A nominal histogram with 0 integral has been found. Please check ! ");
            WriteErrorStatus("SampleHist::SmoothSyst", "            -> Sample: " + fName);
        }

        if(isyst->fIsShape) {
            if (useAlternativeShapeHistFactory) {
                isyst->fHistShapeUp.reset(static_cast<TH1*>(shape_up  ->Clone(  isyst->fHistShapeUp->GetName())));
                isyst->fHistShapeDown.reset(static_cast<TH1*>(shape_down->Clone(isyst->fHistShapeDown->GetName())));
            } else {
                // update shape hists as well
                isyst->fHistShapeUp.reset(static_cast<TH1*>(h_syst_up    ->Clone(isyst->fHistShapeUp->GetName())));
                isyst->fHistShapeDown.reset(static_cast<TH1*>(h_syst_down->Clone(isyst->fHistShapeDown->GetName())));
                if(isyst->fHistShapeUp  ->Integral()>0){
                    isyst->fHistShapeUp  ->Scale(fHist->Integral() / isyst->fHistShapeUp  ->Integral());
                } else {
                    isyst->fHistShapeUp.reset(static_cast<TH1*>(fHist ->Clone(isyst->fHistShapeUp->GetName())));
                }

                if(isyst->fHistShapeDown->Integral() > 0.){
                    isyst->fHistShapeDown->Scale(fHist->Integral() / isyst->fHistShapeDown->Integral());
                } else {
                    isyst->fHistShapeDown.reset(static_cast<TH1*>(fHist ->Clone(isyst->fHistShapeDown->GetName())));
                }
            }
        }
        delete h_syst_up;
        delete h_syst_down;

        delete shape_up;
        delete shape_down;
    }
    fSystSmoothed = true;
}

//_____________________________________________________________________________
//
void SampleHist::CloneSampleHist(SampleHist* h, const std::set<std::string>& names, double scale){
    fName = h->fName;
    if (h->fHist)           fHist           .reset(static_cast<TH1*>(h->fHist->Clone()));
    if (h->fHist_preSmooth) fHist_preSmooth .reset(static_cast<TH1*>(h->fHist_preSmooth->Clone()));
    if (h->fHist_orig)      fHist_orig      .reset(static_cast<TH1*>(h->fHist_orig->Clone()));
    if (fHist) fHist->Scale(scale);
    if (fHist_preSmooth) fHist_preSmooth->Scale(scale);
    if (fHist_orig) fHist_orig->Scale(scale);
    fFileName = h->fFileName;
    fHistoName = h->fHistoName;
    fIsData = h->fIsData;
    fIsSig = h->fIsSig;
    for(const auto& systname : names){
        bool notFound=true;
        for(const auto& isyst : h->fSyst) {
            std::shared_ptr<SystematicHist> syst_tmp = std::make_shared<SystematicHist>("tmp");
            if(systname != isyst->fName) continue;
            TH1* tmp = static_cast<TH1*>(isyst->fHistUp->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistUp.reset(tmp);

            tmp = static_cast<TH1*>(isyst->fHistUp_preSmooth->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistUp_preSmooth.reset(tmp);

            tmp = static_cast<TH1*>(isyst->fHistUp_orig->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistUp_orig.reset(tmp);

            tmp = static_cast<TH1*>(isyst->fHistDown->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistDown.reset(tmp);

            tmp = static_cast<TH1*>(isyst->fHistDown_preSmooth->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistDown_preSmooth.reset(tmp);

            tmp = static_cast<TH1*>(isyst->fHistDown_orig->Clone());
            tmp->Scale(scale);
            syst_tmp->fHistDown_orig.reset(tmp);

            syst_tmp->fName = isyst->fName;
            fSyst.emplace_back(std::move(syst_tmp));
            notFound=false;
        }
        if(notFound){
            std::shared_ptr<SystematicHist> syst_tmp = std::make_shared<SystematicHist>("tmp");
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
    for(std::size_t i_syst = 0; i_syst < fSyst.size(); ++i_syst){
        bool wasIn = false;
        for(std::size_t j_syst = 0; j_syst < h->fSyst.size(); ++j_syst) {
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
void SampleHist::SampleHistAddNominal(SampleHist* h, double scale) {
    if (h->fHist)           fHist          ->Add(h->fHist.get(),          scale);
    if (h->fHist_preSmooth) fHist_preSmooth->Add(h->fHist_preSmooth.get(),scale);
    if (fHist_orig)         fHist_orig     ->Add(h->fHist_orig.get(),     scale);
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
    for(auto& isyst : fSyst) {
        if(!fSample->fUseSystematics) break;
        const std::string systName = isyst->fName;
        const std::string NuisParName = isyst->fSystematic->fNuisanceParameter;
        std::shared_ptr<SystematicHist> syh = sh->GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Divide", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Divide", "Using its nominal. ");
            isyst->Divide( sh->fHist.get() );
        }
        else{
            WriteDebugStatus("SampleHist::Divide", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Divide", "Properly computing with that. ");
            isyst->Divide( syh.get() );
        }
    }
    // loop on all the systematics in the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(const auto& isyst : sh->fSyst) {
        if(!fSample->fUseSystematics) break;
        const std::string systName = isyst->fName;
        const std::string NuisParName = isyst->fSystematic->fNuisanceParameter;
        std::shared_ptr<SystematicHist> syh = GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Divide", "Adding syst "+ NuisParName + " (through syst "+ systName + ") to sample "+ fName);
            std::unique_ptr<TH1> hUp(static_cast<TH1*>(fHist->Clone("h_tmp_up")));
            std::unique_ptr<TH1> hDown(static_cast<TH1*>(fHist->Clone("h_tmp_down")));
            hUp  ->Divide(  sh->fHist.get() );
            hUp  ->Multiply(isyst->fHistUp.get());
            hUp  ->Scale(-1);
            hUp  ->Add(fHist.get(),2);
            //
            hDown->Divide(  sh->fHist.get() );
            hDown->Multiply(isyst->fHistDown.get());
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
            std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*(isyst->fSystematic));
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
    for(auto& isyst : fSyst) {
        if(!fSample->fUseSystematics) break;
        const std::string systName = isyst->fName;
        const std::string NuisParName = isyst->fSystematic->fNuisanceParameter;
        std::shared_ptr<SystematicHist> syh = sh->GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Multiply", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Multiply", "Using its nominal. ");
            isyst->Multiply( sh->fHist.get() );
        }
        else{
            WriteDebugStatus("SampleHist::Multiply", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Multiply", "Properly computing with that. ");
            isyst->Multiply( syh.get() );
        }
    }
    // loop on all the systematics in the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(const auto& isyst : sh->fSyst) {
        if(!fSample->fUseSystematics) break;
        const std::string systName = isyst->fName;
        const std::string NuisParName = isyst->fSystematic->fNuisanceParameter;
        std::shared_ptr<SystematicHist> syh = GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Multiply", "Adding syst "+ NuisParName + " (through syst "+ systName + ") to sample "+ fName);
            std::unique_ptr<TH1> hUp(static_cast<TH1*>(hOrig->Clone("h_tmp_up")));
            std::unique_ptr<TH1> hDown(static_cast<TH1*>(hOrig->Clone("h_tmp_down")));
            hUp  ->Multiply(isyst->fHistUp.get());
            hDown->Multiply(isyst->fHistDown.get());
            syh = AddHistoSyst(NuisParName,NuisParName,hUp.get(),hDown.get());
            if (syh == nullptr) {
                WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
                exit(EXIT_FAILURE);
            }
            syh->fHistUp_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistUp_orig  ->GetName())));
            syh->fHistDown_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistDown_orig->GetName())));
            std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*(isyst->fSystematic));
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
    for(auto& isyst : fSyst){
        if(!fSample->fUseSystematics) break;
        const std::string systName = isyst->fName;
        const std::string NuisParName = isyst->fSystematic->fNuisanceParameter;
        std::shared_ptr<SystematicHist> syh = sh->GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Add", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Add", "Using its nominal. ");
            isyst->Add( sh->fHist.get(), scale );
        }
        else{
            WriteDebugStatus("SampleHist::Add", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Add", "Properly computing with that. ");
            isyst->Add( syh.get(), scale );
        }
    }
    // loop on all the systematics of the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(const auto& isyst : sh->fSyst) {
        if(!fSample->fUseSystematics) break;
        const std::string systName = isyst->fName;
        const std::string NuisParName = isyst->fSystematic->fNuisanceParameter;
        std::shared_ptr<SystematicHist> syh = GetSystFromNP( NuisParName );
        if(syh==nullptr){
            WriteDebugStatus("SampleHist::Add", "Adding syst "+ NuisParName + " (through syst "+ systName + ") to sample "+ fName);
            std::unique_ptr<TH1> hUp(static_cast<TH1*>(hOrig->Clone("h_tmp_up")));
            std::unique_ptr<TH1> hDown(static_cast<TH1*>(hOrig->Clone("h_tmp_down")));
            if (isyst->fHistUp == nullptr) {
               if (TRExFitter::HISTOCHECKCRASH) {
                    WriteErrorStatus("SampleHist::Add", "Systematic "+isyst->fName+ " up var. not found when trying to adding it to "+fSample->fName);
                    exit(EXIT_FAILURE);
               } else {
                    WriteWarningStatus("SampleHist::Add", "Systematic "+isyst->fName+ " up var. not found when trying to adding it to "+fSample->fName);
               }
            }
            else  hUp  ->Add(isyst->fHistUp.get(), scale);
            if (isyst->fHistDown == nullptr) {
               if (TRExFitter::HISTOCHECKCRASH) {
                    WriteErrorStatus("SampleHist::Add", "Systematic "+isyst->fName+ " down var. not found when trying to adding it to "+fSample->fName);
                    exit(EXIT_FAILURE);
               } else {
                    WriteWarningStatus("SampleHist::Add", "Systematic "+isyst->fName+ " down var. not found when trying to adding it to "+fSample->fName);
               }
            }
            else hDown->Add(isyst->fHistDown.get(),scale);
            syh = AddHistoSyst(NuisParName,NuisParName,hUp.get(),hDown.get());
            if (syh == nullptr) {
                WriteErrorStatus("TRExFit::SampleHist", "Histo pointer is nullptr, cannot continue running the code");
                exit(EXIT_FAILURE);
            }
            syh->fHistUp_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistUp_orig  ->GetName())));
            syh->fHistDown_orig.reset(static_cast<TH1*>(fHist_orig->Clone(syh->fHistDown_orig->GetName())));
            std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*(isyst->fSystematic));
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
    for(auto& isyst : fSyst) {
        if(!fSample->fUseSystematics) break;
        isyst->fHistUp->Scale( scale );
        if(isyst->fHistShapeUp!=nullptr) isyst->fHistShapeUp->Scale( scale );
        isyst->fHistDown->Scale( scale );
        if(isyst->fHistShapeDown!=nullptr) isyst->fHistShapeDown->Scale( scale );
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
