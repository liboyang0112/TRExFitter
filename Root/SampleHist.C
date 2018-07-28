#include "TtHFitter/HistoTools.h"
#include "TtHFitter/SampleHist.h"
#include "TtHFitter/StatusLogbook.h"

using namespace std;

// -------------------------------------------------------------------------------------------------
// SampleHist

//_____________________________________________________________________________
//
SampleHist::SampleHist(){
    fName = "";
    //
    fSample = 0x0;
    //
    fHist = 0x0;
    fHist_orig = 0x0;
    //
    fHist_postFit = 0x0;
    //
    fHistoName = "";
    fFileName = "";
    fFitName = "";
    fNSyst = 0;
    fNNorm = 0;
    fNShape = 0;
    fRegionName = "Region";
    fRegionLabel = "Region";
    fVariableTitle = "Variable";
    fSystSmoothed = false;
    fIsMorph.clear();
    //
    fSyst.clear();
}

//_____________________________________________________________________________
//
SampleHist::SampleHist(Sample *sample,TH1 *hist){
    fSample = sample;
    fName = fSample->fName;
    fIsMorph = fSample->fIsMorph;
    //
    //Nominal histogram configuration
    fHist = (TH1*)hist->Clone(Form("h_%s",fName.c_str()));
    fHist->SetFillColor(fSample->fFillColor);
    fHist->SetLineColor(fSample->fLineColor);
    fHist_orig = (TH1*)fHist->Clone(Form("%s_orig",fHist->GetName()));
    //
    fHist_postFit = 0x0;
    //
    fHistoName = "";
    fFileName = "";
    fFitName = "";
    fNSyst = 0;
    fNNorm = 0;
    fNShape = 0;
    fRegionName = "Region";
    fRegionLabel = "Region";
    fVariableTitle = "Variable";
    fSystSmoothed = false;
}

//_____________________________________________________________________________
//
SampleHist::SampleHist(Sample *sample, string histoName, string fileName){
    fSample = sample;
    fHist = HistFromFile(fileName,histoName);
    fHist_orig = HistFromFile(fileName,histoName+"_orig");
    if(fHist_orig==0x0) fHist_orig = (TH1*)fHist->Clone(Form("%s_orig",fHist->GetName()));
//     if(fHist_orig==0x0) fHist_orig = (TH1*)fHist->Clone( (TString)fHist->GetName() + "_orig" );
    fHist->SetFillColor(fSample->fFillColor);
    fHist->SetLineColor(fSample->fLineColor);
    fHist->SetLineWidth(1);
    fName = fSample->fName;
    fIsMorph = fSample->fIsMorph;
    fHistoName = histoName;
    fFileName = fileName;
    fHist_postFit = 0x0;
    fNSyst = 0;
    fNNorm = 0;
    fNShape = 0;
    fRegionName = "Region";
    fVariableTitle = "Variable";
    fSystSmoothed = false;
}

//_____________________________________________________________________________
//
SampleHist::~SampleHist(){
    delete fHist;
    delete fHist_postFit;
    for(unsigned int i = 0; i<fSyst.size(); ++i){
        if(fSyst[i]){
            delete fSyst[i];
        }
    }
    fSyst.clear();
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddOverallSyst(string name,float up,float down){
    SystematicHist *syh;
    // try if it's already there...
    syh = GetSystematic(name);
    // ... and if not create a new one
    if(syh==0x0){
        syh = new SystematicHist(name);
        fSyst.push_back(syh);
        fNSyst ++;
    }
    //
//     syh->fHistUp   = (TH1*)fHist->Clone(Form("%s_%s_Up",fHist->GetName(),name.c_str()));
//     syh->fHistDown = (TH1*)fHist->Clone(Form("%s_%s_Down",fHist->GetName(),name.c_str()));
    syh->fHistUp   = (TH1*)fHist->Clone(Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),name.c_str()));
    syh->fHistDown = (TH1*)fHist->Clone(Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),name.c_str()));
    syh->fHistUp  ->Scale(1.+up);
    syh->fHistDown->Scale(1.+down);
    syh->fHistUp_orig   = (TH1*)syh->fHistUp  ->Clone(Form("%s_orig",syh->fHistUp  ->GetName()));
    syh->fHistDown_orig = (TH1*)syh->fHistDown->Clone(Form("%s_orig",syh->fHistDown->GetName()));
    syh->fIsOverall = true;
    syh->fIsShape   = false;
    syh->fNormUp   = up;
    syh->fNormDown = down;
    return syh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddStatSyst(string name, int i_bin) {
    int bin = i_bin+1; // counting of bins in Root starts with 1, in TRExFitter with 0
    SystematicHist *syh;
    // try if it's already there...
    syh = GetSystematic(name);
    // ... and if not create a new one
    if(syh==0x0){
        syh = new SystematicHist(name);
        fSyst.push_back(syh);
        fNSyst ++;
    }
    double binContent = fHist->GetBinContent(bin);
    double binError = binContent > 1e-4 ? fHist->GetBinError(bin) : 1e-7;
    syh->fHistUp   = (TH1*)fHist->Clone(Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),name.c_str()));
    syh->fHistDown = (TH1*)fHist->Clone(Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),name.c_str()));
    syh->fHistShapeUp   = (TH1*)fHist->Clone(Form("%s_%s_%s_Shape_Up",  fRegionName.c_str(),fSample->fName.c_str(),name.c_str()));
    syh->fHistShapeDown = (TH1*)fHist->Clone(Form("%s_%s_%s_Shape_Down",fRegionName.c_str(),fSample->fName.c_str(),name.c_str()));
    syh->fHistShapeUp  ->SetBinContent(bin, binContent + binError);
    syh->fHistShapeDown->SetBinContent(bin, binContent - binError);
    syh->fHistUp  ->SetBinContent(bin, binContent + binError);
    syh->fHistDown->SetBinContent(bin, binContent - binError);
    syh->fHistUp_orig   = (TH1*)syh->fHistUp  ->Clone(Form("%s_orig",syh->fHistUp  ->GetName()));
    syh->fHistDown_orig = (TH1*)syh->fHistDown->Clone(Form("%s_orig",syh->fHistDown->GetName()));
    syh->fHistShapeUp  ->Scale(fHist->Integral() / syh->fHistShapeUp  ->Integral());
    syh->fHistShapeDown->Scale(fHist->Integral() / syh->fHistShapeDown->Integral());
    syh->fIsOverall = true;
    syh->fIsShape   = true;
    syh->fNormUp   = ( syh->fHistUp->Integral()   - fHist->Integral() ) / fHist->Integral();
    syh->fNormDown = ( syh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
    return syh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddHistoSyst(string name,TH1* h_up,TH1* h_down){

    // before doing anything else, check if the sampleHist can be created
    if(h_up  ==0x0) return 0x0;
    if(h_down==0x0) return 0x0;

    SystematicHist *syh;
    // try if it's already there...
    syh = GetSystematic(name);
    // ... and if not create a new one
    if(syh==0x0){
        syh = new SystematicHist(name);
        fSyst.push_back(syh);
        fNSyst ++;
    }
    //
    syh->fHistUp   = (TH1*)h_up  ->Clone(Form("%s_%s_Up",  fHist->GetName(),name.c_str()));
    syh->fHistDown = (TH1*)h_down->Clone(Form("%s_%s_Down",fHist->GetName(),name.c_str()));
    syh->fHistUp_orig   = (TH1*)h_up  ->Clone(Form("%s_%s_Up_orig",  fHist->GetName(),name.c_str()));
    syh->fHistDown_orig = (TH1*)h_down->Clone(Form("%s_%s_Down_orig",fHist->GetName(),name.c_str()));
    syh->fHistShapeUp   = (TH1*)h_up  ->Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),name.c_str()));
    syh->fHistShapeDown = (TH1*)h_down->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),name.c_str()));
    if(syh->fHistShapeUp  ->Integral() > 0. ){
        syh->fHistShapeUp  ->Scale(fHist->Integral() / syh->fHistShapeUp  ->Integral());
    } else {
        syh->fHistShapeUp  = (TH1*)fHist->Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),name.c_str()));
    }
    if(syh->fHistShapeDown  ->Integral() > 0. ){
        syh->fHistShapeDown->Scale(fHist->Integral() / syh->fHistShapeDown->Integral());
    } else {
        syh->fHistShapeDown  = (TH1*)fHist->Clone(Form("%s_%s_Shape_Down",  fHist->GetName(),name.c_str()));
    }

    syh->fIsOverall = true;
    syh->fIsShape   = true;
    syh->fNormUp   = ( syh->fHistUp->Integral()   - fHist->Integral() ) / fHist->Integral();
    syh->fNormDown = ( syh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
    if(syh->fNormUp == 0 && syh->fNormDown == 0) syh->fIsOverall = false;
    return syh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddHistoSyst(string name,string histoName_up, string fileName_up,string histoName_down, string fileName_down, int pruned/*1: norm only, 2: shape only*/){

    // before doing anything else, check if the sampleHist can be created
    TH1* hUp   = HistFromFile(fileName_up,  histoName_up);
    TH1* hDown = HistFromFile(fileName_down,histoName_down);
    if(hUp  ==0x0) return 0x0;
    if(hDown==0x0) return 0x0;

    SystematicHist *sh;
    // try if it's already there...
    sh = GetSystematic(name);
    // ... and if not create a new one
    if(sh==0x0){
        sh = new SystematicHist(name);
        fSyst.push_back(sh);
        fNSyst ++;
    }
    //
    bool normOnly  = (pruned==1);
    bool shapeOnly = (pruned==2);
    //
    sh->fFileNameUp   = fileName_up;
    sh->fFileNameDown = fileName_down;
    sh->fHistoNameUp   = histoName_up;
    sh->fHistoNameDown = histoName_down;
    sh->fHistUp   = HistFromFile(sh->fFileNameUp,  sh->fHistoNameUp);
    sh->fHistDown = HistFromFile(sh->fFileNameDown,sh->fHistoNameDown);
    sh->fHistUp_orig   = HistFromFile(sh->fFileNameUp,  sh->fHistoNameUp  +"_orig");
    sh->fHistDown_orig = HistFromFile(sh->fFileNameDown,sh->fHistoNameDown+"_orig");
    if(sh->fHistUp   == 0x0) return 0x0;
    if(sh->fHistDown == 0x0) return 0x0;
    if(sh->fHistUp_orig  ==0x0) sh->fHistUp_orig   = (TH1F*)sh->fHistUp->Clone(  Form("%s_orig",sh->fHistUp->GetName()  ));
    if(sh->fHistDown_orig==0x0) sh->fHistDown_orig = (TH1F*)sh->fHistDown->Clone(Form("%s_orig",sh->fHistDown->GetName()));
    //
    if(normOnly){
        sh->fIsShape   = false;
    }
    else{
        sh->fHistShapeUp   = (TH1*)sh->fHistUp  ->Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),name.c_str()));
        sh->fHistShapeDown = (TH1*)sh->fHistDown->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),name.c_str()));
        if(sh->fHistShapeUp  ->Integral() > 0. ){
            sh->fHistShapeUp  -> Scale(fHist->Integral() / sh->fHistShapeUp  ->Integral());
        } else {
            sh->fHistShapeUp = (TH1*)fHist -> Clone(Form("%s_%s_Shape_Up",  fHist->GetName(),name.c_str()));
        }

        if(sh->fHistShapeDown  ->Integral() > 0. ){
            sh->fHistShapeDown->Scale(fHist->Integral() / sh->fHistShapeDown->Integral());
        } else {
            sh->fHistShapeDown = (TH1*)fHist -> Clone(Form("%s_%s_Shape_Down",  fHist->GetName(),name.c_str()));
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
        sh->fNormUp   = ( sh->fHistUp->Integral()   - fHist->Integral() ) / fHist->Integral();
        sh->fNormDown = ( sh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
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
    if(norm==0x0){
        fNormFactors.push_back(normFactor);
        fNNorm ++;
    }
    else{
        norm = normFactor;
    }
    return norm;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::AddNormFactor(string name,float nominal, float min, float max){
    NormFactor *norm = GetNormFactor(name);
    if(norm==0x0){
        fNormFactors.push_back(new NormFactor(name,nominal,min,max));
        fNNorm ++;
    }
    else{
        norm = new NormFactor(name,nominal,min,max);;
    }
    return norm;
}

//_____________________________________________________________________________
//
ShapeFactor* SampleHist::AddShapeFactor(ShapeFactor *shapeFactor){
    ShapeFactor *shape = GetShapeFactor(shapeFactor->fName);
    if(shape==0x0){
        fShapeFactors.push_back(shapeFactor);
        fNShape ++;
    }
    else{
        shape = shapeFactor;
    }
    return shape;
}

//_____________________________________________________________________________
//
ShapeFactor* SampleHist::AddShapeFactor(string name,float nominal, float min, float max){
    ShapeFactor *shape = GetShapeFactor(name);
    if(shape==0x0){
        fShapeFactors.push_back(new ShapeFactor(name,nominal,min,max));
        fNShape ++;
    }
    else{
        shape = new ShapeFactor(name,nominal,min,max);;
    }
    return shape;
}

//_____________________________________________________________________________


//
SystematicHist* SampleHist::GetSystematic(string systName){
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(systName == fSyst[i_syst]->fName) return fSyst[i_syst];
    }
    return 0x0;
}



//
SystematicHist* SampleHist::GetSystFromNP(string NuisParName){
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(NuisParName == fSyst[i_syst]->fSystematic->fNuisanceParameter) return fSyst[i_syst];
    }
    return 0x0;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::GetNormFactor(string name){
    for(int i_syst=0;i_syst<fNNorm;i_syst++){
        if(name == fNormFactors[i_syst]->fName) return fNormFactors[i_syst];
    }
    return 0x0;
}

//_____________________________________________________________________________
//
ShapeFactor* SampleHist::GetShapeFactor(string name){
    for(int i_syst=0;i_syst<fNShape;i_syst++){
        if(name == fShapeFactors[i_syst]->fName) return fShapeFactors[i_syst];
    }
    return 0x0;
}

//_____________________________________________________________________________
//
bool SampleHist::HasSyst(string name){
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSyst[i_syst]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasNorm(string name){
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        if(fNormFactors[i_norm]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasShapeFactor(string name){
    for(int i_shape=0;i_shape<fNShape;i_shape++){
        if(fShapeFactors[i_shape]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
void SampleHist::WriteToFile(TFile *f){
    if(f==0x0){
        if(fHist_orig!=0x0)   WriteHistToFile(fHist_orig,  fFileName);
        if(fHist!=0x0)        WriteHistToFile(fHist,       fFileName);
    }
    else{
        if(fHist_orig!=0x0)   WriteHistToFile(fHist_orig,  f);
        if(fHist!=0x0)        WriteHistToFile(fHist,       f);
    }
    // create the regular binning histogram
    fHist_regBin = HistoTools::TranformHistogramBinning(fHist);
    if(fHist_regBin!=0x0) WriteHistToFile(fHist_regBin,f);
    //
    // save separate gammas as histograms
    if(fSample->fSeparateGammas){
        TH1 *htempUp   = (TH1*)fHist->Clone();
        TH1 *htempDown = (TH1*)fHist->Clone();
        for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
            htempUp  ->AddBinContent(i_bin, 1.*fHist->GetBinError(i_bin));
            htempDown->AddBinContent(i_bin,-1.*fHist->GetBinError(i_bin));
        }
        std::string systName = "stat_"+fSample->fName;
        Systematic *gamma = 0x0;
        if(GetSystematic(systName)) gamma = GetSystematic(systName)->fSystematic;  //GetSystematic(systName);
        if(gamma==0x0) gamma = new Systematic(systName,Systematic::SHAPE);
        WriteDebugStatus("SampleHist::WriteToFile", "adding separate gammas as SHAPE systematic " + systName);
        gamma->fRegions.clear();
        gamma->fRegions.push_back(fRegionName);
        SystematicHist *syh = AddHistoSyst(systName,htempUp,htempDown);
        gamma->fNuisanceParameter = gamma->fName;
        TtHFitter::NPMAP[gamma->fName] = gamma->fNuisanceParameter;
        syh->fSystematic = gamma;
    }
    //
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        // make sure they all have the correct name!
        fSyst[i_syst]->fHistUp  ->SetName( Form("%s_%s_%s_Up",  fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fName.c_str()) );
        fSyst[i_syst]->fHistDown->SetName( Form("%s_%s_%s_Down",fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fName.c_str()) );
        fSyst[i_syst]->fHistUp_orig  ->SetName( Form("%s_%s_%s_Up_orig",  fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fName.c_str()) );
        fSyst[i_syst]->fHistDown_orig->SetName( Form("%s_%s_%s_Down_orig",fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fName.c_str()) );
        if(f==0x0) fSyst[i_syst]->WriteToFile();
        else       fSyst[i_syst]->WriteToFile(f);
        // for shape hist, save also the syst(up)-nominal (to feed HistFactory)
        if(fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
            TH1* hVar = HistoTools::TranformHistogramBinning(
              (TH1*)fSyst[i_syst]->fHistUp->Clone(Form("%s_%s_%s_Up_Var",  fRegionName.c_str(),fSample->fName.c_str(),fSyst[i_syst]->fName.c_str()))
            );
            hVar->Add(fHist_regBin,-1);
            hVar->Divide(fHist_regBin);
            // no negative bins here!
            for(int i_bin=1;i_bin<=hVar->GetNbinsX();i_bin++){
                if(hVar->GetBinContent(i_bin)<0) hVar->SetBinContent(i_bin,-1.*hVar->GetBinContent(i_bin));
            }
            if(f==0x0) WriteHistToFile(hVar,fFileName);
            else       WriteHistToFile(hVar,f);
        }
    }
}

//_____________________________________________________________________________
//
void SampleHist::ReadFromFile(){
    fHist      = HistFromFile(fFileName,fHistoName);
    fHist_orig = HistFromFile(fFileName,fHistoName+"_orig");
}

//_____________________________________________________________________________
//
void SampleHist::FixEmptyBins(const bool suppress){
    //
    // store yields (nominal and systs)
    float yield = fHist->Integral();
    vector<float> yieldUp;
    vector<float> yieldDown;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        SystematicHist* syh = fSyst[i_syst];
        if(syh==0x0) continue;
        if(syh->fHistUp  ==0x0) continue;
        if(syh->fHistDown==0x0) continue;
        yieldUp.push_back(syh->fHistUp->Integral());
        yieldDown.push_back(syh->fHistDown->Integral());
    }
    //
    // store minimum stat unc for non-zero bins
    float minStat = -1;
    for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
        float content = fHist->GetBinContent(i_bin);
        float error   = fHist->GetBinError(  i_bin);
        if(content>0 && error>0){
            if(minStat<0 || error<minStat) minStat = error;
        }
    }
    //
    // loop o bins looking for negatives or zeros
    for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
        float content = fHist->GetBinContent(i_bin);
        float error   = fHist->GetBinError(  i_bin);
        if(content<=0){
            std::string temp = fHist->GetName();
            if (!suppress){
                WriteWarningStatus("SampleHist::FixEmptyBins", "Checking your nominal histogram " +temp + ", the bin " + std::to_string(i_bin) +
                " has a null/negative bin content (content = " + std::to_string(content) + ") ! You should have a look at this !");
                WriteWarningStatus("SampleHist::FixEmptyBins", "    --> For now setting this bin to 1e-06  +/- 1e-06!!! ");
            }
            // set nominal to 10^-6
            fHist->SetBinContent(i_bin,1e-6);
            if(!TtHFitter::GUESSMCSTATERROR){
                fHist -> SetBinError(i_bin, 1e-06);
            } else {
                // if error defined, use it
                if(error>0)        fHist -> SetBinError(i_bin, error);
                // if not, if there was at least one bin with meaningful error, use the smallest
                else if(minStat>0) fHist -> SetBinError(i_bin, minStat);
                // if not, give up and assign a meaningless error ;)
                else               fHist -> SetBinError(i_bin, 1e-06);
            }
            // loop on systematics and set them accordingly
            for(int i_syst=0;i_syst<fNSyst;i_syst++){
                SystematicHist* syh = fSyst[i_syst];
                if(syh->fHistUp  ->GetBinContent(i_bin)<=0) syh->fHistUp  ->SetBinContent(i_bin,1e-06);
                if(syh->fHistDown->GetBinContent(i_bin)<=0) syh->fHistDown->SetBinContent(i_bin,1e-06);
            }
        }
    }
    // set to 0 if negative overall
    if(fHist->Integral()<0){
        for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
            fHist->SetBinContent(i_bin,1e-6);
            fHist->SetBinError(i_bin,1e-6);
        }
    }
    // keep the original overall Normalisation
    else if(fHist->Integral()!=yield){
        fHist->Scale(yield/fHist->Integral());
    }
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        SystematicHist* syh = fSyst[i_syst];
        if(syh->fHistUp  ->Integral()!=yieldUp[i_syst]  ) syh->fHistUp  ->Scale(yieldUp[i_syst]  /syh->fHistUp  ->Integral());
        if(syh->fHistDown->Integral()!=yieldDown[i_syst]) syh->fHistDown->Scale(yieldDown[i_syst]/syh->fHistDown->Integral());
    }
}

//_____________________________________________________________________________
//
void SampleHist::Print(){
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
        if(fSyst[i_syst]->fHistUp!=0x0) fSyst[i_syst]->fHistUp->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistDown!=0x0) fSyst[i_syst]->fHistDown->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistShapeUp!=0x0) fSyst[i_syst]->fHistShapeUp->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistShapeDown!=0x0) fSyst[i_syst]->fHistShapeDown->Rebin(ngroup,"",xbins);
    }
}

//_____________________________________________________________________________
// this draws the control plots (for each systematic) with the syst variations for this region & all sample
void SampleHist::DrawSystPlot( const string &syst, TH1* h_data, bool SumAndData, bool bothPanels ){

    //
    // Draw the distributions for nominal, syst (before and after smoothing)
    //
    float yield_syst_up = 0;
	float yield_syst_down = 0;
	float yield_nominal = 0;
	float yield_data = 0;
    TCanvas *c = new TCanvas("c","c",800,600);
    //
    TPad* pad0 = new TPad("pad0","pad0",0,0.30,1,1,0,0,0);
    pad0->SetTickx(true);
    pad0->SetTicky(true);
    pad0->SetTopMargin(0.05);
    pad0->SetBottomMargin(0.115);
    pad0->SetLeftMargin(0.14);
    pad0->SetRightMargin(0.04);
    pad0->SetFrameBorderMode(0);
    //
    TPad* pad1 = new TPad("pad1","pad1",0,0,1,0.38,0,0,0);
    pad1->SetTickx(true);
    pad1->SetTicky(true);
    pad1->SetTopMargin(0.0);
    pad1->SetBottomMargin(0.27);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.04);
    pad1->SetFrameBorderMode(0);
    //
    pad0->Draw();
    pad1->Draw();
    pad0->cd();

    TH1* h_nominal = nullptr;
    TH1* h_nominal_orig = nullptr;
    TH1* h_dataCopy = nullptr;
    TH1* h_syst_up = nullptr;
    TH1* h_syst_down = nullptr;
    TH1* h_syst_up_orig = nullptr;
    TH1* h_syst_down_orig = nullptr;

    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;
        std::vector < bool > drawRatio;
        drawRatio.push_back(false);
        drawRatio.push_back(true);

        for ( const bool ratioON : drawRatio ){

            if(ratioON) pad1->cd();
            else pad0 -> cd();

            if(SumAndData) h_dataCopy = (TH1*)h_data->Clone();
            h_nominal = (TH1*)fHist->Clone("h_nominal");
            h_nominal->SetLineColor(kBlack);
            h_nominal->SetLineWidth(2);
            h_nominal_orig = (TH1*)fHist_orig->Clone("h_nominal_orig");
            h_nominal_orig->SetLineColor(kBlack);
            h_nominal_orig->SetLineStyle(2);
            h_nominal_orig->SetLineWidth(2);
            if(SumAndData) h_dataCopy->SetMarkerColor(kBlack);
            if(ratioON)h_nominal->SetLineStyle(2);
            else h_nominal->SetLineStyle(1);

            h_nominal->SetFillStyle(0);
            h_nominal_orig->SetFillStyle(0);

            TH1* h_1 = (TH1*)h_nominal->Clone();
            h_nominal->SetMinimum(0);
            h_nominal->SetMaximum(h_nominal->GetMaximum());
            h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
            h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
            h_syst_up_orig = (TH1*)fSyst[i_syst]->fHistUp_orig->Clone();
            h_syst_down_orig = (TH1*)fSyst[i_syst]->fHistDown_orig->Clone();
//             if(FindInStringVector( fSyst[i_syst]->fSystematic->fDropNormIn, "all" )>=0 ||
//                FindInStringVector( fSyst[i_syst]->fSystematic->fDropNormIn, fRegionName )>=0){ // FIXME
//                h_syst_up->Scale(h_nominal->Integral()/h_syst_up->Integral());
//                h_syst_down->Scale(h_nominal->Integral()/h_syst_down->Integral());
//                h_syst_up_orig->Scale(h_nominal->Integral()/h_syst_up_orig->Integral());
//                h_syst_down_orig->Scale(h_nominal->Integral()/h_syst_down_orig->Integral());
//             }
            h_syst_up->SetLineColor(kRed);
            h_syst_up->SetLineWidth(2);
            h_syst_up->SetLineStyle(1);
            h_syst_up->SetFillStyle(0);
            h_syst_down->SetLineColor(kBlue);
            h_syst_down->SetLineWidth(2);
            h_syst_down->SetLineStyle(1);
            h_syst_down->SetFillStyle(0);
            h_syst_up_orig->SetLineColor(kRed);
            h_syst_up_orig->SetLineWidth(2);
            h_syst_up_orig->SetLineStyle(2);
            h_syst_up_orig->SetFillStyle(0);
            h_syst_down_orig->SetLineColor(kBlue);
            h_syst_down_orig->SetLineWidth(2);
            h_syst_down_orig->SetLineStyle(2);
            h_syst_down_orig->SetFillStyle(0);

//             yield_nominal = h_nominal->Integral();
//             yield_syst_up = h_syst_up->Integral();
//             yield_syst_down = h_syst_down->Integral();
//             if(SumAndData) yield_data = h_dataCopy->Integral();
            yield_nominal = CorrectIntegral(h_nominal);
            yield_syst_up = CorrectIntegral(h_syst_up);
            yield_syst_down = CorrectIntegral(h_syst_down);
            if(SumAndData) yield_data = CorrectIntegral(h_dataCopy);

            // draw Relative difference
            h_1->Scale(0);

            if(ratioON){
                h_syst_up->Add(h_nominal,-1);
                h_syst_down->Add(h_nominal,-1);
                if(SumAndData) h_dataCopy->Add(h_nominal,-1);
                h_syst_up->Divide(h_nominal);
                h_syst_down->Divide(h_nominal);
                if(SumAndData) h_dataCopy->Divide(h_nominal);
                // fix empty bins
                for(int i_bin=1;i_bin<=h_nominal->GetNbinsX();i_bin++){
                    if(h_nominal->GetBinContent(i_bin)<1e-5){
                        h_syst_up  ->SetBinContent(i_bin,0.);
                        h_syst_down->SetBinContent(i_bin,0.);
                        if(SumAndData) h_dataCopy->SetBinContent(i_bin,0.);
                    }
                }
                h_syst_up->Scale(100);
                h_syst_down->Scale(100);
                if(SumAndData) h_dataCopy->Scale(100);
                h_syst_up_orig->Add(h_nominal_orig,-1);
                h_syst_down_orig->Add(h_nominal_orig,-1);
                h_syst_up_orig->Divide(h_nominal_orig);
                h_syst_down_orig->Divide(h_nominal_orig);
                // fix empty bins
                for(int i_bin=1;i_bin<=h_nominal_orig->GetNbinsX();i_bin++){
                    if(h_nominal_orig->GetBinContent(i_bin)<1e-5){
                        h_syst_up_orig  ->SetBinContent(i_bin,0.);
                        h_syst_down_orig->SetBinContent(i_bin,0.);
                    }
                }
                h_syst_up_orig->Scale(100);
                h_syst_down_orig->Scale(100);
            }

            double ymax = 0;
            if(!ratioON) ymax = TMath::Max( ymax,TMath::Abs(h_nominal->GetMaximum()));
            if((ratioON || bothPanels) && SumAndData ) ymax = TMath::Max( ymax,TMath::Abs(h_dataCopy->GetMaximum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_up->GetMaximum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_down->GetMaximum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_up->GetMinimum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_down->GetMinimum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_up_orig->GetMaximum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_down_orig->GetMaximum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_up_orig->GetMinimum()));
            ymax = TMath::Max( ymax,TMath::Abs(h_syst_down_orig->GetMinimum()));
            if(!ratioON) {
                h_1->GetYaxis()->SetTitle("Number of events");
                h_1->SetMinimum(1e-05);
                h_1->SetMaximum( ymax*2.0 );
            }
            else {
                h_1->GetYaxis()->SetTitle("#frac{Syst.-Nom.}{Nom.} [%]");
                h_1->GetYaxis()->SetTitleOffset(1.6);
                h_1->GetXaxis()->SetTitleOffset(3.);
                if(TtHFitter::OPTION["SystPlotRatioRange"]!=0){
                    h_1->SetMinimum(-TtHFitter::OPTION["SystPlotRatioRange"]);
                    h_1->SetMaximum( TtHFitter::OPTION["SystPlotRatioRange"]);
                }
                else{
                    h_1->SetMinimum(-ymax*1.5);
                    h_1->SetMaximum( ymax*1.5);
                }
            }
            h_1->GetXaxis()->SetTitle(fVariableTitle.c_str());

            h_1->Draw("HIST");
            if(TtHFitter::SYSTERRORBARS){
                h_syst_down_orig->SetMarkerSize(0);
                h_syst_up_orig->SetMarkerSize(0);
                h_syst_down_orig->DrawCopy("same E");
                h_syst_up_orig->DrawCopy("same E");
            }
            else{
                h_syst_down_orig->DrawCopy("same HIST");
                h_syst_up_orig->DrawCopy("same HIST");
            }
            h_syst_down->DrawCopy("same HIST");
            h_syst_up->DrawCopy("same HIST");
            if(!ratioON){
                h_nominal->DrawCopy("same HIST");
                h_nominal->SetFillStyle(3005);
                h_nominal->SetFillColor(kBlue);
                h_nominal->SetMarkerSize(0);
                h_nominal->DrawCopy("e2same");
                h_nominal_orig->DrawCopy("same HIST");
            }
            else {
	        double xmin=h_nominal->GetBinLowEdge(1);
	        double xmax=h_nominal->GetBinLowEdge(h_nominal->GetNbinsX()+1);
	        TLine* one= new TLine(xmin,0.,xmax,0.);
	        one->SetLineColor(kBlack);
	        //one->SetLineStyle(7);
	        one->SetLineWidth(2);
	        one->Draw("same HIST");
                h_nominal -> SetFillStyle(3005);
                h_nominal -> SetFillColor(kBlue);
                h_nominal -> SetMarkerSize(0);
                for ( int i=1; i<=h_nominal->GetNbinsX(); ++i ) {
                    h_nominal -> SetBinError( i, h_nominal -> GetBinError( i ) * 100. / h_nominal -> GetBinContent( i ) );
                    h_nominal -> SetBinContent( i, 0 );
                }
                h_nominal -> DrawCopy("e2same");
            }
            if((ratioON || bothPanels) && SumAndData ) h_dataCopy->Draw("EX0same");

            if(!ratioON){
                // Creates a legend for the plot
                TLatex *tex = new TLatex();
                tex->SetNDC();
                if(SumAndData) {
                    if(fSyst[i_syst]->fSystematic!=0x0) tex->DrawLatex(0.17,0.79,Form("%s",fSyst[i_syst]->fSystematic->fTitle.c_str()));
                    else                                tex->DrawLatex(0.17,0.79,Form("%s",fSyst[i_syst]->fName.c_str()));
                    //if(fSyst[i_syst]->fSystematic!=0x0) tex->DrawLatex(0.17,0.79,Form("%s, %s",fSyst[i_syst]->fSystematic->fTitle.c_str(),"all samples"));
                    //else                                tex->DrawLatex(0.17,0.79,Form("%s, %s",fSyst[i_syst]->fName.c_str(),"all samples"));
                }
                else{
                    if(fSyst[i_syst]->fSystematic!=0x0) tex->DrawLatex(0.17,0.79,Form("%s, %s",fSyst[i_syst]->fSystematic->fTitle.c_str(),fSample->fTitle.c_str()));
                    else                                tex->DrawLatex(0.17,0.79,Form("%s, %s",fSyst[i_syst]->fName.c_str(),fSample->fTitle.c_str()));
                }
                tex->DrawLatex(0.17,0.72,fRegionLabel.c_str());

                //Legend of the histograms
		TLegend *leg;
		if(SumAndData) leg = new TLegend(0.7,0.71,0.9,0.9);
                else leg = new TLegend(0.7,0.71,0.9,0.85);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                leg->SetTextSize(gStyle->GetTextSize());
                leg->SetTextFont(gStyle->GetTextFont());
                leg->SetMargin(0.2);

                float acc_up = (yield_syst_up-yield_nominal)/yield_nominal;
                string sign_up =  "+";
                if(acc_up<0) sign_up = "-";
                float acc_down = (yield_syst_down-yield_nominal)/yield_nominal;
                string sign_down =  "+";
                if(acc_down<0) sign_down = "-";
                leg->AddEntry(h_syst_up,  Form("+ 1 #sigma (%s%.1f %%)",sign_up.c_str(),  TMath::Abs(acc_up  *100)),"l");
                leg->AddEntry(h_syst_down,Form(" - 1 #sigma (%s%.1f %%)",sign_down.c_str(),TMath::Abs(acc_down*100)),"l");
                leg->Draw();

                //Legend to define the line style
//                 TLegend *leg2 = new TLegend(0.605,0.69,0.9,0.74);
                TLegend *leg2 = new TLegend(0.65,0.64,0.9,0.7);
                leg2->SetFillStyle(0);
                leg2->SetBorderSize(0);
                leg2->SetNColumns(2);
                leg2->SetTextSize(gStyle->GetTextSize());
                leg2->SetTextFont(gStyle->GetTextFont());
                TH1F* h_syst_up_black = (TH1F*)h_syst_up -> Clone();
                h_syst_up_black -> SetLineColor(kBlack);
                TH1F* h_syst_up_origin_black = (TH1F*)h_syst_up_orig -> Clone();
                h_syst_up_origin_black -> SetLineColor(kBlack);
                leg2->AddEntry(h_syst_up_origin_black,"Original","l");
                leg2->AddEntry(h_syst_up_black,"Modified","l");
                leg2 -> Draw();
                if(SumAndData){
                    float acc_data = 0;
					if (yield_nominal != 0) acc_data = (yield_data-yield_nominal)/yield_nominal;
					else acc_data = 99999999;
                    string sign_data =  "+";
                    if(acc_data<0) sign_data = "-";
                    TLegend *leg3 = new TLegend(0.7,0.43,0.9,0.62);
                    leg3->SetFillStyle(0);
                    leg3->SetBorderSize(0);
                    leg3->SetTextSize(gStyle->GetTextSize());
                    leg3->SetTextFont(gStyle->GetTextFont());
                    leg3->SetMargin(0.2);
		    		leg3->AddEntry(h_dataCopy,"Data","p");
		    		leg3->AddEntry(h_nominal,"Total prediction","l");
                    //leg3->AddEntry(h_dataCopy,Form("data (%s%.1f %%)",sign_data.c_str(),TMath::Abs(acc_data*100)),"l");
                    leg3 -> Draw();
                }
            } else {
                TLine line(0.01,1,0.1,1);
                line.SetLineColor(kWhite);
                line.SetLineWidth(20);
                line.DrawLineNDC(0.07,1,0.135,1);
            }
        }

        gSystem->mkdir(fFitName.c_str());
        gSystem->mkdir((fFitName+"/Systematics").c_str());
        gSystem->mkdir((fFitName+"/Systematics/"+fSyst[i_syst]->fName).c_str());

        for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++){
            if(SumAndData) c->SaveAs(Form("%s/Systematics/%s/%s.%s",fFitName.c_str(),fSyst[i_syst]->fName.c_str(), fName.c_str(), TtHFitter::IMAGEFORMAT[i_format].c_str()));
            else c->SaveAs(Form("%s/Systematics/%s/%s.%s",fFitName.c_str(),fSyst[i_syst]->fName.c_str(),fHist->GetName(), TtHFitter::IMAGEFORMAT[i_format].c_str()));
        }
    }
    delete c;
}

//_____________________________________________________________________________
//
void SampleHist::SmoothSyst(const HistoTools::SmoothOption &smoothOpt, string syst,bool force, bool TtresSmoothing){
    if(fSystSmoothed && !force) return;
    TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
    TH1* h_syst_up;
    TH1* h_syst_down;

    for(int i_syst=0;i_syst<fNSyst;i_syst++){

        if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;

        if(fSyst[i_syst]->fHistUp  ==0x0) continue;
        if(fSyst[i_syst]->fHistDown==0x0) continue;
        if(fSyst[i_syst]->fSystematic==0x0) continue;

        h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
        h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
        // h_syst_up = (TH1*)fSyst[i_syst]->fHistUp_orig->Clone();
        // h_syst_down = (TH1*)fSyst[i_syst]->fHistDown_orig->Clone();

        if(fSyst[i_syst]->fSystematic->fPreSmoothing){
            TH1* h_tmp_up   = h_syst_up!=0x0   ? (TH1*)h_syst_up  ->Clone() : 0x0;
            TH1* h_tmp_down = h_syst_down!=0x0 ? (TH1*)h_syst_down->Clone() : 0x0;
            if(h_tmp_up!=0x0 || h_tmp_down!=0x0){
                TH1* h_tmp_nominal = (TH1*)h_nominal  ->Clone();
                for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++){
                    h_tmp_nominal->GetBinError(i_bin,0);
                }
                if(h_tmp_up!=0x0){
                    float tmp_nom_up = h_tmp_up->Integral();
                    h_tmp_up->Add(h_tmp_nominal,-1);
                    h_tmp_up->Divide(h_tmp_nominal);
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_up->AddBinContent(i_bin, 100.);
//                     SmoothHistogram(h_tmp_up,-1,3);
                    h_tmp_up->Smooth();
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_up->AddBinContent(i_bin,-100.);
                    h_tmp_up->Multiply(h_tmp_nominal);
                    h_tmp_up->Add(h_tmp_nominal, 1);
                    h_tmp_up->Scale(tmp_nom_up/h_tmp_up->Integral());
                    h_syst_up = (TH1*)h_tmp_up->Clone();
                }
                if(h_tmp_down!=0x0){
                    float tmp_nom_down = h_tmp_down->Integral();
                    h_tmp_down->Add(h_tmp_nominal,-1);
                    h_tmp_up->Divide(h_tmp_nominal);
                    for(int i_bin=1;i_bin<=h_tmp_nominal->GetNbinsX();i_bin++) h_tmp_down->AddBinContent(i_bin, 100.);
                    h_tmp_down->Smooth();
//                     SmoothHistogram(h_tmp_down,-1,3);
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

        if(fSyst[i_syst]->fSmoothType + fSyst[i_syst]->fSymmetrisationType<=0){
            HistoTools::Scale(fSyst[i_syst]->fHistUp,  fHist,fSyst[i_syst]->fScaleUp);
            HistoTools::Scale(fSyst[i_syst]->fHistDown,fHist,fSyst[i_syst]->fScaleDown);
            continue;
        }

        //
        // Call the function for smoothing and symmetrisation
        //
        if(fSyst[i_syst]->fIsShape){
            HistoTools::ManageHistograms(   fSyst[i_syst]->fSmoothType + fSyst[i_syst]->fSymmetrisationType,//parameters of the histogram massaging
                                            h_nominal,//nominal histogram
                                            fSyst[i_syst]->fHistUp, fSyst[i_syst]->fHistDown,//original histograms
                                            h_syst_up, h_syst_down, //modified histograms
                                            fSyst[i_syst]->fScaleUp,fSyst[i_syst]->fScaleDown, // scale factors
                                            smoothOpt,
                                            TtresSmoothing // alternative smoothing
                                         );
        }
        //
        // need to ad these lines to make sure overall only systematics get scaled as well
        else{
            HistoTools::Scale(fSyst[i_syst]->fHistUp,  fHist,fSyst[i_syst]->fScaleUp);
            HistoTools::Scale(fSyst[i_syst]->fHistDown,fHist,fSyst[i_syst]->fScaleDown);
        }

        //
        // keep the variation below 100% in each bin, if the option Smooth is set for the sample
        //
        if(fSample->fSmooth){
            if(h_syst_up!=0x0){
                for(int iBin = 1; iBin <= h_syst_up  ->GetNbinsX(); ++iBin ){
                    float relDiff = (h_syst_up->GetBinContent(iBin) - h_nominal->GetBinContent(iBin))/ h_nominal->GetBinContent(iBin);
                    if(relDiff>=1. ) h_syst_up->SetBinContent(iBin, (1.+0.99)*h_nominal->GetBinContent(iBin) );
                    if(relDiff<=-1.) h_syst_up->SetBinContent(iBin, (1.-0.99)*h_nominal->GetBinContent(iBin) );
                }
            }
            if(h_syst_down!=0x0){
                for(int iBin = 1; iBin <= h_syst_down  ->GetNbinsX(); ++iBin ){
                    float relDiff = (h_syst_down->GetBinContent(iBin) - h_nominal->GetBinContent(iBin))/ h_nominal->GetBinContent(iBin);
                    if(relDiff>=1. ) h_syst_down->SetBinContent(iBin, (1.+0.99)*h_nominal->GetBinContent(iBin) );
                    if(relDiff<=-1.) h_syst_down->SetBinContent(iBin, (1.-0.99)*h_nominal->GetBinContent(iBin) );
                }
            }
        }

        //
        // Save stuff
        //
//         fSyst[i_syst]->fHistUp_orig = (TH1*)fSyst[i_syst]->fHistUp->Clone();
        fSyst[i_syst]->fHistUp = (TH1*)h_syst_up->Clone(fSyst[i_syst]->fHistUp->GetName());
//         fSyst[i_syst]->fHistDown_orig = (TH1*)fSyst[i_syst]->fHistDown->Clone();
        fSyst[i_syst]->fHistDown = (TH1*)h_syst_down->Clone(fSyst[i_syst]->fHistUp->GetName());

        //
        // Perform a check of the output histograms (check for 0 bins and other pathologic behaviours)
        //
//         HistoTools::CheckHistograms( h_nominal /*nominal*/, fSyst[i_syst] /*systematic*/, fSample -> fType != Sample::SIGNAL, true /*cause crash if problem*/);
        HistoTools::CheckHistograms( h_nominal /*nominal*/, fSyst[i_syst] /*systematic*/, fSample -> fType != Sample::SIGNAL, TtHFitter::HISTOCHECKCRASH /*cause crash if problem*/);

        //
        // Normalisation component first
        //
        if(h_nominal->Integral()!=0){
            fSyst[i_syst]->fNormUp   = fSyst[i_syst]->fHistUp  ->Integral()/h_nominal->Integral() - 1.;
            fSyst[i_syst]->fNormDown = fSyst[i_syst]->fHistDown->Integral()/h_nominal->Integral() - 1.;
        } else {
            WriteErrorStatus("SampleHist::SmoothSyst", "A nominal histogram with 0 intergral has been found. Please check ! ");
            WriteErrorStatus("SampleHist::SmoothSyst", "            -> Sample: " + fName);
        }

        if(fSyst[i_syst]->fIsShape){
            // update shape hists as well
            fSyst[i_syst]->fHistShapeUp   = (TH1*)h_syst_up  ->Clone(fSyst[i_syst]->fHistShapeUp  ->GetName());
            fSyst[i_syst]->fHistShapeDown = (TH1*)h_syst_down->Clone(fSyst[i_syst]->fHistShapeDown->GetName());
            if(fSyst[i_syst]->fHistShapeUp  ->Integral()>0){
                fSyst[i_syst]->fHistShapeUp  ->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeUp  ->Integral());
            } else {
                fSyst[i_syst]->fHistShapeUp  = (TH1*)fHist ->Clone(fSyst[i_syst]->fHistShapeUp  ->GetName());
            }

            if(fSyst[i_syst]->fHistShapeDown->Integral() > 0.){
                fSyst[i_syst]->fHistShapeDown->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeDown->Integral());
            } else {
                fSyst[i_syst]->fHistShapeDown  = (TH1*)fHist ->Clone(fSyst[i_syst]->fHistShapeDown  ->GetName());
            }
        }
    }
    fSystSmoothed = true;
}

//_____________________________________________________________________________
//
void SampleHist::CloneSampleHist(SampleHist* h, std::set<std::string> names, float scale){
    fName = h->fName;
    fHist = (TH1*)h->fHist->Clone();
    fHist_orig = (TH1*)h->fHist_orig->Clone();
    if (scale >= 0 && scale <=1){
        fHist->Scale(scale*fHist->Integral());
        fHist_orig->Scale(scale*fHist_orig->Integral());
    } else {
        WriteDebugStatus("SampleHist::CloneSampleHist", "Not scaling the histograms because the scale was set to: " +std::to_string(scale));
    }
    fFileName = h->fFileName;
    fHistoName = h->fHistoName;
    fIsData = h->fIsData;
    fIsSig = h->fIsSig;
    fNSyst = h->fNSyst;
    for(auto systname : names ){
        bool notFound=true;
        for(int i_syst=0; i_syst<h->fNSyst; i_syst++){
            SystematicHist* syst_tmp = new SystematicHist("tmp");
            if(systname!=h->fSyst[i_syst]->fName) continue;
            syst_tmp->fHistUp = (TH1*)h->fSyst[i_syst]->fHistUp->Clone();
            syst_tmp->fHistUp_orig = (TH1*)h->fSyst[i_syst]->fHistUp_orig->Clone();
            syst_tmp->fHistDown = (TH1*)h->fSyst[i_syst]->fHistDown->Clone();
            syst_tmp->fHistDown_orig = (TH1*)h->fSyst[i_syst]->fHistDown_orig->Clone();
            syst_tmp->fName = h->fSyst[i_syst]->fName;
            fSyst.push_back(syst_tmp);
            notFound=false;
        }
        if(notFound){
            SystematicHist* syst_tmp = new SystematicHist("tmp");
            ++fNSyst;
            syst_tmp->fHistUp = (TH1*)h->fHist->Clone();
            syst_tmp->fHistUp_orig = (TH1*)h->fHist->Clone();
            syst_tmp->fHistDown = (TH1*)h->fHist->Clone();
            syst_tmp->fHistDown_orig = (TH1*)h->fHist->Clone();
            syst_tmp->fName = systname;
            fSyst.push_back(syst_tmp);
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
void SampleHist::SampleHistAdd(SampleHist* h){
    fHist->Add(h->fHist);
    fHist_orig->Add(h->fHist_orig);
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        bool wasIn=false;
        for(int j_syst=0;j_syst<h->fNSyst;j_syst++){
            if(fSyst[i_syst]->fName==h->fSyst[j_syst]->fName){
                fSyst[i_syst]->fHistUp->Add(h->fSyst[j_syst]->fHistUp);
                fSyst[i_syst]->fHistDown->Add(h->fSyst[j_syst]->fHistDown);
                fSyst[i_syst]->fHistUp_orig->Add(h->fSyst[j_syst]->fHistUp_orig);
                fSyst[i_syst]->fHistDown_orig->Add(h->fSyst[j_syst]->fHistDown_orig);
                wasIn=true;
            }
        }
        if(wasIn) continue;
        fSyst[i_syst]->fHistUp->Add(h->fHist);
        fSyst[i_syst]->fHistDown->Add(h->fHist);
        fSyst[i_syst]->fHistUp_orig->Add(h->fHist);
        fSyst[i_syst]->fHistDown_orig->Add(h->fHist);
    }
}

//_____________________________________________________________________________
//
void SampleHist::Divide(SampleHist *sh){
    TH1* hOrig = (TH1*)fHist->Clone("h_tmp_orig");
    fHist->Divide( sh->fHist );
    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        string systName = fSyst[i_syst]->fName;
//        SystematicHist *syh = sh->GetSystematic( systName );
        string NuisParName = fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = sh->GetSystFromNP( NuisParName );
        if(syh==0x0){
            fSyst[i_syst]->Divide( sh->fHist );
            WriteDebugStatus("SampleHist::Divide", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Divide", "Using its nominal. ");
        }
        else{
            fSyst[i_syst]->Divide( syh );
            WriteDebugStatus("SampleHist::Divide", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Divide", "Properly computing with that. ");
        }
    }
    // loop on all the systematics in the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
     for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        string systName = sh->fSyst[i_syst]->fName;
//        SystematicHist *syh = GetSystematic( systName );
        string NuisParName = sh->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = GetSystFromNP( NuisParName );
        if(syh==0x0){
            WriteDebugStatus("SampleHist::Divide", "Adding syst "+ systName + " to sample "+ fName);
            TH1* hUp   = (TH1*)fHist->Clone("h_tmp_up"  );
            TH1* hDown = (TH1*)fHist->Clone("h_tmp_down");
            hUp  ->Divide(   sh->fHist );
            hUp  ->Multiply( sh->fSyst[i_syst]->fHistUp  );
            hUp  ->Scale(-1);
            hUp  ->Add(fHist,2);
            //
            hDown->Divide(   sh->fHist );
            hDown->Multiply( sh->fSyst[i_syst]->fHistDown);
            hDown->Scale(-1);
            hDown->Add(fHist,2);
            //
            syh = AddHistoSyst(systName,hUp,hDown);
            syh->fHistUp_orig   = (TH1*)fHist_orig->Clone(syh->fHistUp_orig  ->GetName());
            syh->fHistDown_orig = (TH1*)fHist_orig->Clone(syh->fHistDown_orig->GetName());
            syh->fSystematic = sh->fSyst[i_syst]->fSystematic;
            fSample->AddSystematic(sh->fSyst[i_syst]->fSystematic);
            delete hUp;
            delete hDown;
        }
    }
    delete hOrig;
}

//_____________________________________________________________________________
//
void SampleHist::Multiply(SampleHist *sh){
    TH1* hOrig = (TH1*)fHist->Clone("h_tmp_orig");
    fHist->Multiply( sh->fHist );
    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        string systName = fSyst[i_syst]->fName;
//        SystematicHist *syh = sh->GetSystematic( systName );
        string NuisParName = fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = sh->GetSystFromNP( NuisParName );
        if(syh==0x0){
            fSyst[i_syst]->Multiply( sh->fHist );
            WriteDebugStatus("SampleHist::Multiply", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Multiply", "Using its nominal. ");
        }
        else{
            fSyst[i_syst]->Multiply( syh );
            WriteDebugStatus("SampleHist::Multiply", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Multiply", "Properly computing with that. ");
        }
    }
    // loop on all the systematics in the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        string systName = sh->fSyst[i_syst]->fName;
//        SystematicHist *syh = GetSystematic( systName );
        string NuisParName = sh->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = GetSystFromNP( NuisParName );
        if(syh==0x0){
            WriteDebugStatus("SampleHist::Multiply", "Adding syst "+ systName + " to sample "+ fName);
            TH1* hUp   = (TH1*)hOrig->Clone("h_tmp_up"  );
            TH1* hDown = (TH1*)hOrig->Clone("h_tmp_down");
            hUp  ->Multiply( sh->fSyst[i_syst]->fHistUp   );
            hDown->Multiply( sh->fSyst[i_syst]->fHistDown );
            syh = AddHistoSyst(systName,hUp,hDown);
            syh->fHistUp_orig   = (TH1*)fHist_orig->Clone(syh->fHistUp_orig  ->GetName());
            syh->fHistDown_orig = (TH1*)fHist_orig->Clone(syh->fHistDown_orig->GetName());
            syh->fSystematic = sh->fSyst[i_syst]->fSystematic;
            fSample->AddSystematic(sh->fSyst[i_syst]->fSystematic);
            delete hUp;
            delete hDown;
        }
    }
    delete hOrig;
}

//_____________________________________________________________________________
//
void SampleHist::Add(SampleHist *sh,float scale){
    TH1* hOrig = (TH1*)fHist->Clone("h_tmp_orig");
    fHist->Add( sh->fHist, scale );
    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        string systName = fSyst[i_syst]->fName;
//        SystematicHist *syh = sh->GetSystematic( systName );
        string NuisParName = fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = sh->GetSystFromNP( NuisParName );
        if(syh==0x0){
            fSyst[i_syst]->Add( sh->fHist, scale );
            WriteDebugStatus("SampleHist::Add", "Syst. "+ systName +"(" + NuisParName +")"+ " not present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Add", "Using its nominal. ");
        }
        else{
            fSyst[i_syst]->Add( syh, scale );
            WriteDebugStatus("SampleHist::Add", "Syst. "+ systName +"(" + NuisParName +")"+ " present in  "+ sh->fName);
            WriteDebugStatus("SampleHist::Add", "Properly computing with that. ");
        }
    }
    // loop on all the systematics the the other SampleHist, and see if some of them are NOT in this
    // if so, add a new SystematicHist
    for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        string systName = sh->fSyst[i_syst]->fName;
//        SystematicHist *syh = GetSystematic( systName );
        string NuisParName = sh->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        SystematicHist *syh = GetSystFromNP( NuisParName );
        if(syh==0x0){
            WriteDebugStatus("SampleHist::Add", "Adding syst "+ systName + " to sample "+ fName);
            TH1* hUp   = (TH1*)hOrig->Clone("h_tmp_up"  );
            TH1* hDown = (TH1*)hOrig->Clone("h_tmp_down");
            hUp  ->Add( sh->fSyst[i_syst]->fHistUp  ,scale );
            hDown->Add( sh->fSyst[i_syst]->fHistDown,scale );
            syh = AddHistoSyst(systName,hUp,hDown);
            syh->fHistUp_orig   = (TH1*)fHist_orig->Clone(syh->fHistUp_orig  ->GetName());
            syh->fHistDown_orig = (TH1*)fHist_orig->Clone(syh->fHistDown_orig->GetName());
            syh->fSystematic = sh->fSyst[i_syst]->fSystematic;
            fSample->AddSystematic(sh->fSyst[i_syst]->fSystematic);
            delete hUp;
            delete hDown;
        }
    }
    delete hOrig;
}

//_____________________________________________________________________________
//
void SampleHist::Scale(float scale){
    fHist->Scale( scale );
    // loop on all the systematics in this SampleHist
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(!fSample->fUseSystematics) break;
        fSyst[i_syst]->fHistUp->Scale( scale );
        fSyst[i_syst]->fHistUp_orig->Scale( scale );
        fSyst[i_syst]->fHistDown->Scale( scale );
        fSyst[i_syst]->fHistDown_orig->Scale( scale );
    }
}
