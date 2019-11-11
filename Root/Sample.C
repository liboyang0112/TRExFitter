// Class include
#include "TRExFitter/Sample.h"

// Framework includes
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"


// -------------------------------------------------------------------------------------------------
// Sample

//__________________________________________________________________________________
//
Sample::Sample(const std::string& name,int type) :
    fName(name),
    fType(type),
    fFitName(name),
    fTitle(name),
    fTexTitle(""),
    fGroup(""),
    fFillColor(kWhite),
    fLineColor(kBlack),
    fNormalizedByTheory(true),
    fIgnoreSelection(""),
    fIgnoreWeight(""),
    fUseMCStat(true),
    fUseSystematics(true),
    fDivideBy(""),
    fMultiplyBy(""),
    fNormToSample(""),
    fSmooth(false),
    fBuildPullTable(0),
    fSelection("1"),
    fMCweight("1"),
    fNSyst(0),
    fNNorm(0),
    fNShape(0),
    fSeparateGammas(false),
    fMCstatScale(1.),
    fCorrelateGammasWithSample(""),
    fSystFromSample("") {
}

//__________________________________________________________________________________
//
Sample::~Sample(){
    for(auto isf : fShapeFactors) {
        delete isf;
    }
}

//__________________________________________________________________________________
// cosmetics
void Sample::SetTitle(const std::string& title){
    fTitle = title;
}

//__________________________________________________________________________________
//
void Sample::SetFillColor(int color){
    fFillColor = color;
}

//__________________________________________________________________________________
//
void Sample::SetLineColor(int color){
    fLineColor = color;
}

//__________________________________________________________________________________
//
void Sample::NormalizedByTheory(const bool norm ){
    fNormalizedByTheory = norm;
}

//__________________________________________________________________________________
// read from ntuples
void Sample::SetMCweight(const std::string& weight){
    fMCweight = weight;
}

//__________________________________________________________________________________
//
void Sample::SetSelection(const std::string& selection){
    fSelection = selection;
}

//__________________________________________________________________________________
//
void Sample::AddNtuplePath(const std::string& path){
    fNtuplePaths.push_back(path);
}

//__________________________________________________________________________________
//
void Sample::AddNtupleFile(const std::string& file){
    fNtupleFiles.push_back(file);
}

//__________________________________________________________________________________
//
void Sample::AddNtupleName(const std::string& name){
    fNtupleNames.push_back(name);
}

//__________________________________________________________________________________
// read from histograms
void Sample::AddHistoPath(const std::string& path){
    fHistoPaths.push_back(path);
}

//__________________________________________________________________________________
//
void Sample::AddHistoFile(const std::string& file){
    fHistoFiles.push_back(file);
}

//__________________________________________________________________________________
//
void Sample::AddHistoName(const std::string& name){
    fHistoNames.push_back(name);
}

//__________________________________________________________________________________
// norm factors and systs
void Sample::AddNormFactor(NormFactor* normFactor){
    fNormFactors.emplace_back(normFactor);
    fNNorm ++;
}

//__________________________________________________________________________________
//
void Sample::AddShapeFactor(ShapeFactor* shapeFactor){
    fShapeFactors.push_back(shapeFactor);
    fNShape ++;
}

//__________________________________________________________________________________
//
void Sample::AddSystematic(Systematic* syst){
    fSystematics.emplace_back(syst);
    fNSyst++;
}

//__________________________________________________________________________________
//
bool Sample::HasSystematic(const std::string& name) const{
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSystematics[i_syst]->fName==name) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
bool Sample::HasNuisanceParameter(const std::string& name) const{
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSystematics[i_syst]->fNuisanceParameter==name) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
bool Sample::HasNormFactor(const std::string& name) const{
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        if(fNormFactors[i_norm]->fName==name) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
NormFactor* Sample::AddNormFactor(const std::string& name,double nominal,double min,double max,bool isConst){
    fNormFactors.emplace_back(new NormFactor(name,nominal,min,max,isConst));
    fNNorm ++;
    return fNormFactors[fNNorm-1].get();
}

//__________________________________________________________________________________
//
ShapeFactor* Sample::AddShapeFactor(const std::string& name,double nominal,double min,double max,bool isConst){
    fShapeFactors.push_back(new ShapeFactor(name,nominal,min,max,isConst));
    fNShape ++;
    return fShapeFactors[fNShape-1];
}

//__________________________________________________________________________________
//
Systematic* Sample::AddSystematic(const std::string& name,int type,double up,double down){
    fSystematics.emplace_back(new Systematic(name,type,up,down));
    fNSyst++;
    return fSystematics[fNSyst-1].get();
}
