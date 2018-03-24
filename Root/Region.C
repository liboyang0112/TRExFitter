#include <iomanip>
#include "TtHFitter/Region.h"
#include "TtHFitter/StatusLogbook.h"

#include "TMatrixD.h"

using namespace std;

// -------------------------------------------------------------------------------------------------
// class Region

//__________________________________________________________________________________
//
Region::Region(string name){
    fName = name;
    fLabel = name;
    fTexLabel = "";
    fShortLabel = name;
    fRegionType = CONTROL;
    fNBkg = 0;
    fNSig = 0;
    fNSyst = 0;
    fNNorm = 0;
    fNShape = 0;
    fHasData = false;
    fHasSig = false;
    fNSamples = 0;
    fUseStatErr = false;
    fHistoBins = 0;
    fHistoNBinsRebin = -1;
    fYTitle = "";
    fYmaxScale = 0;
    fYmax = 0;
    fRatioYmax = 2.;
    fRatioYmin = 0;
    fRatioYmaxPostFit = 1.5;
    fRatioYminPostFit = 0.5;

    string cName = "c_"+fName;
    int canvasWidth = 600;
    int canvasHeight = 700;
    if(TtHFitter::OPTION["CanvasWidth"]!=0)  canvasWidth  = TtHFitter::OPTION["CanvasWidth"];
    if(TtHFitter::OPTION["CanvasHeight"]!=0) canvasHeight = TtHFitter::OPTION["CanvasHeight"];
    if(canvasWidth!=600 || canvasHeight!=700) fPlotPreFit = new TthPlot(cName,canvasWidth,canvasHeight);
    else                                      fPlotPreFit = new TthPlot(cName);
    fPlotPreFit->fShowYields = TtHFitter::SHOWYIELDS;
    cName = "c_"+fName+"_postFit";
    if(canvasWidth!=600 || canvasHeight!=700) fPlotPostFit = new TthPlot(cName,canvasWidth,canvasHeight);
    else                                      fPlotPostFit = new TthPlot(cName);
    fPlotPostFit->fShowYields = TtHFitter::SHOWYIELDS;
    
    fSampleHists.clear();
    fSamples.clear();
    fSystematics.clear();
    fNormFactors.clear();
    fShapeFactors.clear();

    fAlternativeVariables.clear();
    
    fIntCode_overall = 4;
    fIntCode_shape = 0;

    fCorrVar1 = "";
    fCorrVar2 = "";
    
    fFitName = "";
    fFitType = TtHFit::SPLUSB;
    fPOI = "";
    fFitLabel = "";
    
    fSelection = "1";
    fMCweight = "1";
    
    fLogScale = false;
    fSkipSmoothing = false;
    
    fBinWidth = 0;
    
    fLumiScale = 1.;
    
    fBlindingThreshold = -1;
    
    fRegionDataType = REALDATA;

    fBinTransfo = "";
    fTransfoDzBkg = 0.;
    fTransfoDzSig = 0.;
    fTransfoFzBkg = 0.;
    fTransfoFzSig = 0.;
    fTransfoJpar1 = 0.;
    fTransfoJpar2 = 0.;
    fTransfoJpar3 = 0.;
    fAutoBinBkgsInSig.clear();
    
    fATLASlabel = "Internal";
    
    fSuffix = "";
    fGroup = "";
    
    fBlindedBins = 0x0;
    fKeepPrefitBlindedBins = false;
    fGetChi2 = 0;
    
    fDropBins.clear();

    fBinLabels.clear();
    
    fChi2val = -1;
    fNDF = -1;
    fChi2prob = -1;
    
    fUseGammaPulls = false;
    
    fTot = 0x0;
}

//__________________________________________________________________________________
//
Region::~Region(){
    if(fPlotPreFit) delete fPlotPreFit;
    if(fPlotPostFit) delete fPlotPostFit;
    
    for(unsigned int i = 0; i < fSampleHists.size(); ++i){
        if(fSampleHists[i]){
            delete fSampleHists[i];
        }
    }
    fSampleHists.clear();
}

//__________________________________________________________________________________
//
SampleHist* Region::SetSampleHist(Sample *sample, string histoName, string fileName){
    fSampleHists.push_back(new SampleHist( sample, histoName, fileName ));
    if(sample->fType==Sample::DATA){
        fHasData = true;
        fData = fSampleHists[fNSamples];
    }
    else if(sample->fType==Sample::SIGNAL){
        fHasSig = true;
        fSig[fNSig] = fSampleHists[fNSamples];
        fNSig ++;
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg[fNBkg] = fSampleHists[fNSamples];
        fNBkg ++;
    }
    else if(sample->fType==Sample::GHOST){
        WriteDebugStatus("Region::SetSampleHist", "Adding GHOST sample.");
    }
    else{
        WriteErrorStatus("Region::SetSampleHist", "SampleType not supported.");
    }
    fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists[fNSamples]->fRegionName = fName;
    fSampleHists[fNSamples]->fRegionLabel = fLabel;
    fSampleHists[fNSamples]->fFitName = fFitName;
    fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
    fNSamples++;
    return fSampleHists[fNSamples-1];
}

//__________________________________________________________________________________
//
SampleHist* Region::SetSampleHist(Sample *sample, TH1* hist ){
    fSampleHists.push_back(new SampleHist( sample, hist ));
    if(sample->fType==Sample::DATA){
        fHasData = true;
        fData = fSampleHists[fNSamples];
    }
    else if(sample->fType==Sample::SIGNAL){
        fHasSig = true;
        fSig[fNSig] = fSampleHists[fNSamples];
        fNSig ++;
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg[fNBkg] = fSampleHists[fNSamples];
        fNBkg ++;
    }
    else if(sample->fType==Sample::GHOST){
        WriteDebugStatus("Region::SetSampleHist", "Adding GHOST sample.");
    }
    else{
        WriteErrorStatus("Region::SetSampleHist", "SampleType not supported.");
    }
    fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists[fNSamples]->fRegionName = fName;
    fSampleHists[fNSamples]->fRegionLabel = fLabel;
    fSampleHists[fNSamples]->fFitName = fFitName;
    fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
    fNSamples++;
    return fSampleHists[fNSamples-1];
}

//__________________________________________________________________________________
//
void Region::AddSample(Sample* sample){
    fSamples.push_back(sample);
    fNSamples++;
}

//__________________________________________________________________________________
//
void Region::AddSystematic(Systematic *syst){
    fSystematics.push_back(syst);
    fNSyst++;
}

//__________________________________________________________________________________
//
void Region::SetBinning(int N, double *bins){
    fNbins = fHistoNBinsRebin = N;
    fHistoBins = new double[N+1];
    for(int i=0; i<=N; ++i) fHistoBins[i] = bins[i];
}

//__________________________________________________________________________________
//
void Region::Rebin(int N){
    fNbins = fHistoNBinsRebin = N;
}

//__________________________________________________________________________________
//
void Region::SetRegionType( RegionType type ){
    fRegionType = type;
}

//__________________________________________________________________________________
//
void Region::SetRegionDataType( DataType type ){
    fRegionDataType = type;
}

//__________________________________________________________________________________
//
SampleHist* Region::GetSampleHist(string &sampleName){
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSampleHists[i_smp]->fName == sampleName) return fSampleHists[i_smp];
    }
    return 0x0;
}

//__________________________________________________________________________________
//
void Region::BuildPreFitErrorHist(){
    WriteInfoStatus("Region::BuildPreFitErrorHist", "Building pre-fit plot for region " + fName + " ...");
    //
    float yieldNominal(0.), yieldUp(0.), yieldDown(0.);
    float diffUp(0.), diffDown(0.);
    //
    fSystNames.clear();
    fNpNames.clear();
    std::map<string,bool> systIsThere;
    systIsThere.clear();
    TH1* hUp = 0x0;
    TH1* hDown = 0x0;
    string systName = "";
    string systNuisPar = "";
    SystematicHist *sh = 0x0;
    
    //
    // Collect all the systematics on all the samples
    // (in parallel also colect all the nuisance parameters)
    //
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType == Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType == Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType == Sample::SIGNAL && !TtHFitter::SHOWSTACKSIG) continue;
        
        // NOTE: probably we don't need NormFactors here, as they don't introduce prefit uncertainty...
//         //
//         // Norm factors
//         //
//         for(int i_norm=0;i_norm<fSampleHists[i]->fNNorm;i_norm++){
//             systName = fSampleHists[i]->fNormFactors[i_norm]->fName;
//             // skip POI if B-only fit
//             if(fFitType==TtHFit::BONLY && systName==fPOI) continue;
//             if(!systIsThere[systName]){
//                 fSystNames.push_back(systName);
//                 systIsThere[systName] = true;
//             }
//             systNuisPar = systName;
// //             if(fSampleHists[i]->fSyst[i_syst]->fSystematic!=0x0) systNuisPar = fSampleHists[i]->fSyst[i_syst]->fSystematic->fNuisanceParameter;
//             if(FindInStringVector(fNpNames,systNuisPar)<0){
//                 fNpNames.push_back(systNuisPar);
//             }
//         }
        
        //
        // Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
            systName = fSampleHists[i]->fSyst[i_syst]->fName;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
            if(fSampleHists[i]->fSyst[i_syst]->fSystematic!=0x0)
                systNuisPar = fSampleHists[i]->fSyst[i_syst]->fSystematic->fNuisanceParameter;
            if(FindInStringVector(fNpNames,systNuisPar)<0){
                fNpNames.push_back(systNuisPar);
            }
        }
    }
    
    //
    // Build pre-fit error hists (for each sample and each systematic)
    //

    // - loop on samples
    //
    for(int i=0;i<fNSamples;i++){
        
        // skip data
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !TtHFitter::SHOWSTACKSIG) continue;
        
        WriteDebugStatus("Region::BuildPreFitErrorHist", "  Sample: " + fSampleHists[i]->fName);
        
        // - loop on systematics
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            WriteDebugStatus("Region::BuildPreFitErrorHist", "    Systematic: " + fSystNames[i_syst]);
            systName = fSystNames[i_syst];
            
            // get SystematicHist
            sh = fSampleHists[i]->GetSystematic(systName);
            
            // hack: add a systematic hist if not there... FIXME
            if(sh==0x0){
                fSampleHists[i]->AddHistoSyst(systName,fSampleHists[i]->fHist,fSampleHists[i]->fHist);
                sh = fSampleHists[i]->GetSystematic(systName);
                // initialize the up and down variation histograms
                // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
                sh->fHistUp   = (TH1*)fSampleHists[i]->fHist->Clone(Form("%s_%s_Up",  fSampleHists[i]->fHist->GetName(),systName.c_str()));
                sh->fHistDown = (TH1*)fSampleHists[i]->fHist->Clone(Form("%s_%s_Down",fSampleHists[i]->fHist->GetName(),systName.c_str()));
            }
            //
            // - loop on bins
            for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
                diffUp = 0.;
                diffDown = 0.;
                hUp = 0x0;
                hDown = 0x0;
                yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
                // if it's a systematic (NB: skip Norm-Factors!!)
                if(fSampleHists[i]->HasSyst(fSystNames[i_syst])){
                    hUp   = sh->fHistUp;
                    hDown = sh->fHistDown;
                    if(hUp!=0x0)    yieldUp     = hUp  ->GetBinContent(i_bin);
                    else            yieldUp     = yieldNominal;
                    if(hDown!=0x0)  yieldDown   = hDown->GetBinContent(i_bin);
                    else            yieldDown   = yieldNominal;
                    diffUp   += yieldUp   - yieldNominal;
                    diffDown += yieldDown - yieldNominal;
                }
                WriteDebugStatus("Region::BuildPreFitErrorHist", "        Bin " + std::to_string(i_bin) + ":  " + " \t +" + std::to_string(100*diffUp/yieldNominal)
                 + "%\t " + std::to_string(100*diffDown/yieldNominal) + "%");
            }
        }
    }
    
    // at this point all the sample-by-sample pre-fit variation histograms should be filled
    //
    // Now build the total prediction variations, for each systematic
    // - loop on systematics
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        systName = fSystNames[i_syst];
        
        // initialize the tot variation hists
        fTotUp[i_syst]   = (TH1*)fTot->Clone(Form("h_%s_tot_%s_Up",  fName.c_str(), systName.c_str()));
        fTotDown[i_syst] = (TH1*)fTot->Clone(Form("h_%s_tot_%s_Down",fName.c_str(), systName.c_str()));
        // - loop on bins
        for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
            diffUp = 0.;
            diffDown = 0.;
            // - loop on samples
            for(int i=0;i<fNSamples;i++){
                //
                // scale according to NormFactors
                float scale = 1.;
                for(unsigned int i_nf=0;i_nf<fSampleHists[i]->fSample->fNormFactors.size();i_nf++){
                    NormFactor *nf = fSampleHists[i]->fSample->fNormFactors[i_nf];
                    // if this norm factor is a morphing one
                    if(nf->fName.find("morph_")!=string::npos){
                        std::string formula = TtHFitter::SYSTMAP[nf->fName];
                        std::string name = TtHFitter::NPMAP[nf->fName];
                        formula = ReplaceString(formula,name,"x");
                        TF1* f_morph = new TF1("f_morph",formula.c_str(),nf->fMin,nf->fMax);
                        scale *= f_morph->Eval(nf->fNominal);
                    }
                    else {
                        scale *= fSampleHists[i]->fSample->fNormFactors[i_nf]->fNominal;
                    }
                }
                //
                // skip data
                if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
                if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
                if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !TtHFitter::SHOWSTACKSIG) continue;
                // get SystematicHist
                sh = fSampleHists[i]->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);
                diffUp   += (sh->fHistUp  ->GetBinContent(i_bin) - yieldNominal)*scale;
                diffDown += (sh->fHistDown->GetBinContent(i_bin) - yieldNominal)*scale;
            }
            // add the proper bin content to the variation hists
            fTotUp[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }
    
    WriteDebugStatus("Region::BuildPreFitErrorHist", "----");
    
    //
    // build the vectors of variations (sum histograms for systematics with the same NP)
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        // look for all the already stored systematics, to find if one had the same NP
        bool found = false;
        for(int j_syst=0;j_syst<i_syst;j_syst++){
            if(TtHFitter::NPMAP[fSystNames[i_syst]]==TtHFitter::NPMAP[fSystNames[j_syst]]){
                found = true;
//                 h_up[   FindInStringVector(fNpNames,TtHFitter::NPMAP[fSystNames[i_syst]]) ]->Add(fTotUp[  i_syst]);
//                 h_down[ FindInStringVector(fNpNames,TtHFitter::NPMAP[fSystNames[i_syst]]) ]->Add(fTotDown[i_syst]);
                int whichsyst = FindInStringVector(fNpNames,TtHFitter::NPMAP[fSystNames[i_syst]]);
                TH1*h_diff_up   = (TH1*) h_up[whichsyst]  ->Clone(Form("%s_%s","clone_",h_up[whichsyst]  ->GetName()));
                TH1*h_diff_down = (TH1*) h_down[whichsyst]->Clone(Form("%s_%s","clone_",h_down[whichsyst]->GetName()));
                h_diff_up  ->Add(fTotUp[  i_syst],fTot,1,-1);
                h_diff_down->Add(fTotDown[i_syst],fTot,1,-1);
                h_up[   FindInStringVector(fNpNames,TtHFitter::NPMAP[fSystNames[i_syst]])]->Add(h_diff_up);
                h_down[ FindInStringVector(fNpNames,TtHFitter::NPMAP[fSystNames[i_syst]])]->Add(h_diff_down);
                break;
            }
        }
        if(found) continue;
        //
        // if shape-only systematics to show, normalize each syst variation
        if(TtHFitter::OPTION["ShapeOnlySystBand"]>0){
            fTotUp[i_syst]  ->Scale(fTot->Integral()/fTotUp[i_syst]  ->Integral());
            fTotDown[i_syst]->Scale(fTot->Integral()/fTotDown[i_syst]->Integral());
        }
        h_up.  push_back( fTotUp[i_syst]   );
        h_down.push_back( fTotDown[i_syst] );
    }
    fErr = BuildTotError( fTot, h_up, h_down, fNpNames );
    fErr->SetName("g_totErr");
    // at this point fTot and fErr should be ready
    
    //
    // Goodness of pre-fit
    if(fGetChi2){
        // remove blinded bins
        TH1* h_data = (TH1*)fData->fHist->Clone();
        if(fBlindedBins!=0x0){
            for(int i_bin=1;i_bin<=h_data->GetNbinsX();i_bin++){
                if(fBlindedBins->GetBinContent(i_bin)>0) h_data->SetBinContent(i_bin,-1);
                if(find(fDropBins.begin(),fDropBins.end(),i_bin-1)!=fDropBins.end()) h_data->SetBinContent(i_bin,-1);
            }
        }
        if(fGetChi2==1) fNpNames.clear();
        std::pair<double,int> res = GetChi2Test( h_data, fTot, h_up, h_down, fNpNames );
        fChi2val = res.first;
        fNDF = res.second;
        fChi2prob = ROOT::Math::chisquared_cdf_c( res.first, res.second);
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- PRE-FIT AGREEMENT EVALUATION -----------------------");
        if(fGetChi2==1)
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- -------- STAT-ONLY --------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "--- REGION " + fName + ":");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  chi2        = " + std::to_string(fChi2val));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  ndof        = " + std::to_string(fNDF));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  probability = " + std::to_string(fChi2prob));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
    }
}

//__________________________________________________________________________________
//
TthPlot* Region::DrawPreFit(string opt){
    
    TthPlot *p = fPlotPreFit;
    if(fYmaxScale==0) p->SetYmaxScale(1.8);
    else              p->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) p->fYmax = fYmax;
    if(fYmin!=0) p->fYmin = fYmin;
    p->fRatioYmax = fRatioYmax;
    p->fRatioYmin = fRatioYmin;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") p->SetYaxis(fYTitle);
    //
    // For 4-top-style plots
    if(TtHFitter::OPTION["FourTopStyle"]>0){
        if(fRegionType==CONTROL) p->AddLabel("#font[62]{Control Region}");
        else if(fRegionType==SIGNAL) p->AddLabel("#font[62]{Signal Region}");
        else if(fRegionType==VALIDATION) p->AddLabel("#font[62]{Validation Region}");
        else p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        p->AddLabel("#font[52]{Pre-fit}");
        p->fLegendNColumns = TtHFitter::OPTION["LegendNColumns"];
    }
    //
    // old-style plots
    else{
        p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        if(TtHFitter::OPTION["NoPrePostFit"]==0) p->AddLabel("Pre-Fit");
        if(TtHFitter::OPTION["LegendNColumns"]!=0) p->fLegendNColumns = TtHFitter::OPTION["LegendNColumns"];
    }
    //
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    if(fBlindingThreshold>=0) p->SetBinBlinding(true,fBlindingThreshold);

    if(fBinLabels.size() && ((int)fBinLabels.size()==fNbins)) {
      for(int i_bin=0; i_bin<fNbins; i_bin++) {
        p->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
      }
    }
    
    //
    // build h_tot
    //
    fTot = 0x0;
    string title;
    TH1* h = 0x0;
    if(fHasData && opt.find("blind")==string::npos) p->SetData(fData->fHist,fData->fSample->fTitle);
    for(int i=0;i<fNSig;i++){
        title = fSig[i]->fSample->fTitle;
        if(fSig[i]->fSample->fGroup != "") title = fSig[i]->fSample->fGroup;
        h = (TH1*)fSig[i]->fHist->Clone();
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if(!fSig[i]->fSample->fUseMCStat && !fSig[i]->fSample->fSeparateGammas){
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++) h->SetBinError(i_bin,0.);
        }
        // else still check if the value is reasonable
        else{
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++){
                if(h->GetBinError(i_bin)>10*h->GetBinContent(i_bin)) h->SetBinError(i_bin,h->GetBinContent(i_bin));
            }
        }
        // scale it according to NormFactors
        for(unsigned int i_nf=0;i_nf<fSig[i]->fSample->fNormFactors.size();i_nf++){
            NormFactor *nf = fSig[i]->fSample->fNormFactors[i_nf];
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos){
                std::string formula = TtHFitter::SYSTMAP[nf->fName];
                std::string name = TtHFitter::NPMAP[nf->fName];
                formula = ReplaceString(formula,name,"x");
                TF1* f_morph = new TF1("f_morph",formula.c_str(),nf->fMin,nf->fMax);
                float scale = f_morph->Eval(nf->fNominal);
                h->Scale(scale);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fSig[i]->fSample->fName + " by " + std::to_string(scale));
            }
            else{
                h->Scale(nf->fNominal);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fSig[i]->fSample->fName + " by " + std::to_string(fSig[i]->fSample->fNormFactors[i_nf]->fNominal));
            }
        }
        if(TtHFitter::SHOWSTACKSIG)   p->AddSignal(    h,title);
        if(TtHFitter::SHOWNORMSIG){
            if( (TtHFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
             || !TtHFitter::OPTION["NormSigSRonly"] )
                p->AddNormSignal(h,title);
        }
        else{
            if(TtHFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(h,title);
        }
        if(TtHFitter::SHOWOVERLAYSIG) p->AddOverSignal(h,title);
        if(TtHFitter::SHOWSTACKSIG){
            if(fTot==0x0) fTot = (TH1*)h->Clone("h_tot");
            else          fTot->Add(h);
        }
    }
    for(int i=0;i<fNBkg;i++){
        title = fBkg[i]->fSample->fTitle;
        if(fBkg[i]->fSample->fGroup != "") title = fBkg[i]->fSample->fGroup;
        h = (TH1*)fBkg[i]->fHist->Clone();
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if(!fBkg[i]->fSample->fUseMCStat && !fBkg[i]->fSample->fSeparateGammas){
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++) h->SetBinError(i_bin,0.);
        }
        // else still check if the value is reasonable
        else{
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++){
                if(h->GetBinError(i_bin)>10*h->GetBinContent(i_bin)) h->SetBinError(i_bin,h->GetBinContent(i_bin));
            }
        }
        // scale it according to NormFactors
        for(unsigned int i_nf=0;i_nf<fBkg[i]->fSample->fNormFactors.size();i_nf++){
            NormFactor *nf = fBkg[i]->fSample->fNormFactors[i_nf];
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos){
                std::string formula = TtHFitter::SYSTMAP[nf->fName];
                std::string name = TtHFitter::NPMAP[nf->fName];
                formula = ReplaceString(formula,name,"x");
                TF1* f_morph = new TF1("f_morph",formula.c_str(),nf->fMin,nf->fMax);
                float scale = f_morph->Eval(nf->fNominal);
                h->Scale(scale);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fBkg[i]->fSample->fName + " by " + std::to_string(scale));
            }
            else{
                h->Scale(nf->fNominal);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fBkg[i]->fSample->fName + " by " + std::to_string(fBkg[i]->fSample->fNormFactors[i_nf]->fNominal));
            }
        }
        p->AddBackground(h,title);
        if(fTot==0x0) fTot = (TH1*)h->Clone("h_tot");
        else          fTot->Add(h);
    }
    
    //
    // set error to 0 if no MCstat
    //
    if(!fUseStatErr){
        for(int i_bin=1;i_bin<=fTot->GetNbinsX();i_bin++){
            fTot->SetBinError(i_bin,0);
        }
    }
    else{
        for(int i_bin=0;i_bin<fTot->GetNbinsX()+2;i_bin++){
            if(fTot->GetBinError(i_bin)>10*fTot->GetBinContent(i_bin)) fTot->SetBinError(i_bin,fTot->GetBinContent(i_bin));
        }
    }
    
    p->SetTotBkg((TH1*)fTot);
    p->BlindData();
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);
    fBlindedBins = p->h_blinding;
    
    //
    // Computes the uncertainty bands arround the h_tot histogram
    //
    BuildPreFitErrorHist();
    
    //
    // Print chi2 info
    //
    if(fGetChi2 && TtHFitter::SHOWCHI2)  p->SetChi2KS(fChi2prob,-1,fChi2val,fNDF);
    
    //
    // Sets the last ingredients in the TthPlot object
    //
    p->SetTotBkgAsym(fErr);
    p->fATLASlabel = fATLASlabel;
    if(fLogScale) opt += " log";
    p->Draw(opt);
    return p;
}

//__________________________________________________________________________________
//
double Region::GetMultFactors( FitResults *fitRes, std::ofstream& pullTex,
                                const int i /*sample*/, const int i_bin /*bin number*/, 
                                const double binContent0,
                                const std::string &var_syst_name,
                                const bool isUp ){

    double multNorm(1.), multShape(0.);
    float systValue = 0;
    for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
        std::string systName = fSampleHists[i]->fSyst[i_syst]->fName;
        TString systNameNew(systName); // used in pull tables
        if(fSampleHists[i]->fSyst[i_syst]->fSystematic!=0x0)
            systName = fSampleHists[i]->fSyst[i_syst]->fSystematic->fNuisanceParameter;
        if(fSampleHists[i]->fSyst[i_syst]->fSystematic!=0x0){
            if(fSampleHists[i]->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
                continue;
            }
            else {
                if(var_syst_name!="" && systName==var_syst_name){
                    systValue = fitRes->GetNuisParValue(systName) + (isUp ? fitRes->GetNuisParErrUp(systName) : fitRes->GetNuisParErrDown(systName));
                } else {
                    systValue = fitRes->GetNuisParValue(systName);
                }
            }
        } else{
            systValue = fitRes->GetNuisParValue(systName);
        }

        //
        // Normalisation component: use the exponential interpolation and the multiplicative combination
        //
        if(fSampleHists[i]->fSyst[i_syst]->fIsOverall){
            float binContentUp   = (fSampleHists[i]->fSyst[i_syst]->fNormUp+1) * binContent0;
            float binContentDown = (fSampleHists[i]->fSyst[i_syst]->fNormDown+1) * binContent0;
            multNorm *= (GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall));
            if (fSampleHists[i]->fSample->fBuildPullTable>0){
                if ((((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall)) > 1.01) || ((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall)) < 0.99)) && (i_bin==1)) {
                    WriteDebugStatus("Region::DrawPostFit", "Syst " + systName +" in bin " + std::to_string(i_bin) + " has norm effect " +
                                  std::to_string( ((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall)) -1 )*100));
                    pullTex  << setprecision(2) << "norm "<< systNameNew.ReplaceAll("_","-") << "&"<< ((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall)) -1 )*100  << " \\% \\\\\n"<< endl;
                }
            }
        }

        
        //
        // Shape component: use the linear interpolation and the additive combination
        //
        if(fSampleHists[i]->fSyst[i_syst]->fIsShape){
            float binContentUp   = fSampleHists[i]->fSyst[i_syst]->fHistShapeUp->GetBinContent(i_bin);
            float binContentDown = fSampleHists[i]->fSyst[i_syst]->fHistShapeDown->GetBinContent(i_bin);
            multShape += (GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape) -1 );
            if (fSampleHists[i]->fSample->fBuildPullTable==2){
                if (((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape) -1 ) > 0.03) || ((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape) -1 ) < - 0.03)) {
                    WriteDebugStatus("Region::DrawPostFit", "Syst " + systName +" in bin " + std::to_string(i_bin) + " has shape effect " +
                                  std::to_string( ((GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape)) -1 )*100));
                    pullTex << setprecision(2) <<"shape "<<  systNameNew.ReplaceAll("_","-") <<" bin " << i_bin << "&"<<(GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape) -1 )*100  << " \\% \\\\\n"<< endl;
                }
            }
        }
    }
    return multNorm*(multShape+1.);
}

//__________________________________________________________________________________
//
void Region::BuildPostFitErrorHist(FitResults *fitRes){
    
    WriteInfoStatus("Region::BuildPostFitErrorHist", "Building post-fit plot for region " + fName + " ...");
    float diffUp(0.), diffDown(0.);
    
    //
    // 0) Collect all the systematics on all the samples
    //
    fSystNames.clear();
    std::map<string,bool> systIsThere;
    systIsThere.clear();
    float systValue(0.);
    float systErrUp(0.);
    float systErrDown(0.);
    //TH1* hUp = 0x0;
    //TH1* hDown = 0x0;
    string systName = "";
    string systNuisPar = "";
    SystematicHist *sh = 0x0;
    string systNameSF = "";
    int iBinSF = 0;
    
    for(int i_sample=0;i_sample<fNSamples;i_sample++){
        if(fSampleHists[i_sample]->fSample->fType == Sample::DATA) continue;
        if(fSampleHists[i_sample]->fSample->fType == Sample::GHOST) continue;
        if(fSampleHists[i_sample]->fSample->fType == Sample::SIGNAL && !TtHFitter::SHOWSTACKSIG) continue;
        
        //
        // Norm factors
        //
        for(int i_norm=0;i_norm<fSampleHists[i_sample]->fNNorm;i_norm++){
            NormFactor *nf = fSampleHists[i_sample]->fNormFactors[i_norm];
            systName = nf->fName;
//             // if this norm factor is a morphing one => save the nuis.par
//             if(nf->fName.find("morph_")!=string::npos){
//                 systName = TtHFitter::NPMAP[nf->fName];
//             }
            // skip POI if B-only fit FIXME
            if(fFitType==TtHFit::BONLY && systName==fPOI) continue;
            if(nf->fConst) continue;            
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }

        //
        // Shape factors
        //
        // extract number of bins
        TH1* hSFTmp;
        hSFTmp = (TH1*)fSampleHists[i_sample]->fHist->Clone();
        // loop over shape factors
        for(int i_shape=0;i_shape<fSampleHists[i_sample]->fNShape;i_shape++){
            systName = fSampleHists[i_sample]->fShapeFactors[i_shape]->fName;
//             systName = fSampleHists[i_sample]->fShapeFactors[i_shape]->fNuisanceParameter;
            // add syst name for each bin
            for(int i_bin = 0; i_bin < hSFTmp->GetNbinsX(); i_bin++){     
                systNameSF = systName + "_bin_" + std::to_string(i_bin);
                // the shape factor naming used i_bin - 1 for the first bin
                // add it as one syst per bin
                if(!systIsThere[systNameSF]){
                    fSystNames.push_back(systNameSF);
                    systIsThere[systNameSF] = true;
                }
            }
        }
        hSFTmp->~TH1();
        
        //
        // Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i_sample]->fNSyst;i_syst++){
            systName = fSampleHists[i_sample]->fSyst[i_syst]->fName;
//             if(fSampleHists[i_sample]->fSyst[i_syst]->fSystematic!=0x0)
//                 systName = fSampleHists[i_sample]->fSyst[i_syst]->fSystematic->fNuisanceParameter;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }
        
        //
        // Gammas
        //
        if(fUseGammaPulls && (fSampleHists[i_sample]->fSample->fUseMCStat || fSampleHists[i_sample]->fSample->fSeparateGammas)){
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                std::string gammaName = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                if(fSampleHists[i_sample]->fSample->fSeparateGammas)
                            gammaName = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i_sample]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                if(!systIsThere[gammaName]){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                }
            }
        }
    }
    
    //
    // 1) Build post-fit error hists (for each sample and each systematic):
    //
    
    // - loop on samples
    for(int i=0;i<fNSamples;i++){
        
        // skip data
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        WriteVerboseStatus("Region::BuildPostFitErrorHist", "  Sample: " + fSampleHists[i]->fName);
        
        // - loop on systematics
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            
            WriteVerboseStatus("Region::BuildPostFitErrorHist", "    Systematic: " + fSystNames[i_syst]);
            
            //
            // Get fit result
            //
            systName    = fSystNames[i_syst];
            if(TtHFitter::NPMAP[systName]=="") TtHFitter::NPMAP[systName] = systName;
            systValue   = fitRes->GetNuisParValue(TtHFitter::NPMAP[systName]);
            systErrUp   = fitRes->GetNuisParErrUp(TtHFitter::NPMAP[systName]);
            systErrDown = fitRes->GetNuisParErrDown(TtHFitter::NPMAP[systName]);
            
            WriteVerboseStatus("Region::BuildPostFitErrorHist", "      alpha = " + std::to_string(systValue) + " +" + std::to_string(systErrUp) + " " + std::to_string(systErrDown));
            
            // this to include (prefit) error from SHAPE syst
            if(fSampleHists[i]->GetSystematic(systName)){
                if(fSampleHists[i]->GetSystematic(systName)->fSystematic!=0x0){
                    if(fSampleHists[i]->GetSystematic(systName)->fSystematic->fType==Systematic::SHAPE){
                        systValue   = 0;
                        systErrUp   = 1;
                        systErrDown = -1;
                    }
                }
            }
            
            //
            // Get SystematicHist
            //
            sh = fSampleHists[i]->GetSystematic(systName);
            
            // hack: add a systematic hist if not there
            if(sh==0x0){
                fSampleHists[i]->AddHistoSyst(systName,fSampleHists[i]->fHist,fSampleHists[i]->fHist);
                sh = fSampleHists[i]->GetSystematic(systName);
            }
            
            //
            // initialize the up and down variation histograms
            // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
            //
            sh->fHistUp_postFit   = (TH1*)fSampleHists[i]->fHist_postFit->Clone(Form("%s_%s_Up_postFit",  fSampleHists[i]->fHist->GetName(),systName.c_str()));
            sh->fHistDown_postFit = (TH1*)fSampleHists[i]->fHist_postFit->Clone(Form("%s_%s_Down_postFit",fSampleHists[i]->fHist->GetName(),systName.c_str()));

            // - loop on bins
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                diffUp = 0.;
                diffDown = 0.;
                //hUp = 0x0;
                //hDown = 0x0;
                double yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
                double yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);  // store nominal yield for this bin, but do it post fit
                double yieldNominal_postFit_nfOnly = yieldNominal;
                for(auto nf : fSampleHists[i]->fNormFactors){
                    yieldNominal_postFit_nfOnly *= fitRes->GetNuisParValue(TtHFitter::NPMAP[nf->fName]);
                }
                
                size_t posTmp = systName.find("_bin_");
                std::string gammaName      = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                std::string gammaNameShape = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                //
                // if it's a gamma
                if(gammaName==fSystNames[i_syst] && fSampleHists[i]->fSample->fUseMCStat && !fSampleHists[i]->fSample->fSeparateGammas){
                    diffUp   += yieldNominal_postFit*systErrUp;
                    diffDown += yieldNominal_postFit*systErrDown;
                }
                //
                // if it's a specific-sample gamma
                else if(gammaName==fSystNames[i_syst] && fSampleHists[i]->fSample->fSeparateGammas){
                    diffUp   += yieldNominal_postFit*systErrUp;
                    diffDown += yieldNominal_postFit*systErrDown;
                }
                //
                // if it's a norm factor
                else if(fSampleHists[i]->HasNorm(fSystNames[i_syst])){
                    // if this norm factor is a morphing one
                    if(fSystNames[i_syst].find("morph_")!=string::npos){
                        std::string formula = TtHFitter::SYSTMAP[fSystNames[i_syst]];
                        std::string name = TtHFitter::NPMAP[fSystNames[i_syst]];
                        formula = ReplaceString(formula,name,"x");
                        TF1* f_morph = new TF1("f_morph",formula.c_str(),fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fMin,fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fMax);
                        // float scaleUp   = f_morph->Eval(systValue+systErrUp);
                        // float scaleDown = f_morph->Eval(systValue-systErrDown);
                        // diffUp   += yieldNominal*(scaleUp  -1);
                        // diffDown += yieldNominal*(scaleDown-1);
                        float scaleNom  = f_morph->Eval(systValue);
                        float scaleUp   = f_morph->Eval(systValue+systErrUp);
                        float scaleDown = f_morph->Eval(systValue+systErrDown);
                        diffUp   += yieldNominal_postFit*(scaleUp/scaleNom  -1);
                        diffDown += yieldNominal_postFit*(scaleDown/scaleNom-1);
                    }
                    else{
                        diffUp   += yieldNominal_postFit*systErrUp/systValue;
                        diffDown += yieldNominal_postFit*systErrDown/systValue;
                    }
                }
                //
                // ShapeFactor have to get NP per bin
                else if(posTmp != std::string::npos){
                    // get the shape factor name without bin index
                    systNameSF = systName.substr(0, posTmp);
                    // get the shape factor bin as integer
                    iBinSF = std::atoi(systName.substr(posTmp + 5).c_str()) + 1;
                    // FIXME could still be a problem with pruning?
                    if(fSampleHists[i]->HasShapeFactor(systNameSF) && iBinSF == i_bin){
                        // diffUp   += yieldNominal*systErrUp;
                        // diffDown += yieldNominal*systErrDown;
                        diffUp   += yieldNominal_postFit*systErrUp;
                        diffDown += yieldNominal_postFit*systErrDown;
                    }
                }
                //
                // Systematics treatment 
                //
                else if(fSampleHists[i]->HasSyst(fSystNames[i_syst])){
                    std::ofstream dummy;
                    double multNom = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal );
                    double multUp = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal, fSystNames[i_syst], true);
                    double multDown = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal, fSystNames[i_syst], false);
                    diffUp   += (multUp/multNom - 1.)*yieldNominal_postFit;
                    diffDown += (multDown/multNom - 1.)*yieldNominal_postFit;
                }
                
                WriteVerboseStatus("Region::BuildPostFitErrorHist", "        Bin " + std::to_string(i_bin) + ":   " + "\t +" + std::to_string(100*diffUp/yieldNominal)
                    + "%\t " + std::to_string(100*diffDown/yieldNominal) + "%");

                //
                // Add the proper bin content to the variation hists (coming from post-fit total histogram)
                //
                sh->fHistUp_postFit  ->AddBinContent( i_bin, diffUp   );
                sh->fHistDown_postFit->AddBinContent( i_bin, diffDown );
            }
        }
    }
    
    // at this point all the sample-by-sample post-fit variation histograms should be filled
    //
    // Now build the total prediction variations, for each systematic
    //
    
    // - loop on systematics
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        systName = fSystNames[i_syst];
        //
        // Initialize the tot variation hists
        //
        fTotUp_postFit[i_syst]   = (TH1*)fTot_postFit->Clone(Form("h_tot_%s_Up_postFit",  systName.c_str()));
        fTotUp_postFit[i_syst] -> Scale(0); //initialising the content to 0
        fTotDown_postFit[i_syst] = (TH1*)fTot_postFit->Clone(Form("h_tot_%s_Down_postFit",systName.c_str()));
        fTotDown_postFit[i_syst] -> Scale(0); //initialising the content to 0
        
        // - loop on bins
        for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
            diffUp = 0.;
            diffDown = 0.;
            // - loop on samples
            for(int i=0;i<fNSamples;i++){
                // skip data
                if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
                if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
                if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !TtHFitter::SHOWSTACKSIG) continue;
                // skip signal if Bkg only
                if(fFitType==TtHFit::BONLY && fSampleHists[i]->fSample->fType==Sample::SIGNAL) continue;
                // get SystematicHist
                sh = fSampleHists[i]->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                // yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);
                if(!sh) continue;
                diffUp   += sh->fHistUp_postFit  ->GetBinContent(i_bin);
                diffDown += sh->fHistDown_postFit->GetBinContent(i_bin);
                // if(fSampleHists[i]->fSample->fIsMorph){
                //     if (std::fabs(yieldNominal_postFit) > 1e-6){
                //         if(sh->fHistUp_postFit)   diffUp   += sh->fHistUp_postFit  ->GetBinContent(i_bin) - yieldNominal_postFit;
                //         if(sh->fHistDown_postFit) diffDown += sh->fHistDown_postFit->GetBinContent(i_bin) - yieldNominal_postFit;
                //     }
                // } else {
                //     if(sh->fHistUp_postFit)   diffUp   += sh->fHistUp_postFit  ->GetBinContent(i_bin) - yieldNominal_postFit;
                //     if(sh->fHistDown_postFit) diffDown += sh->fHistDown_postFit->GetBinContent(i_bin) - yieldNominal_postFit;
                // }
            }
            // add the proper bin content to the variation hists
            fTotUp_postFit[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown_postFit[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }
    
    // at this point all the total expectation post-fit variation histograms should be filled
    //
    // Build the vectors of variations (necessary to call BuildTotError)
    //
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    std::vector<std::string> systNuisPars;
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        h_up.  push_back( fTotUp_postFit[i_syst]   );
        h_down.push_back( fTotDown_postFit[i_syst] );
        systNuisPars.push_back(TtHFitter::NPMAP[fSystNames[i_syst]]);
    }
    fErr_postFit = BuildTotError( fTot_postFit, h_up, h_down, systNuisPars, fitRes->fCorrMatrix );
    fErr_postFit->SetName("g_totErr_postFit");
    // at this point fTot and fErr _postFit should be ready
    
    //
    // Goodness of post-fit
    if(fGetChi2){
        // remove blinded bins
        TH1* h_data = (TH1*)fData->fHist->Clone();
        if(fBlindedBins!=0x0){
            for(int i_bin=1;i_bin<=h_data->GetNbinsX();i_bin++){
                if(fBlindedBins->GetBinContent(i_bin)>0) h_data->SetBinContent(i_bin,-1);
                if(find(fDropBins.begin(),fDropBins.end(),i_bin-1)!=fDropBins.end()) h_data->SetBinContent(i_bin,-1);
            }
        }
        if(fGetChi2==1) fSystNames.clear();
        std::pair<double,int> res = GetChi2Test( h_data, fTot_postFit, h_up, h_down, fSystNames, fitRes->fCorrMatrix );
        fChi2val = res.first;
        fNDF = res.second;
        fChi2prob = ROOT::Math::chisquared_cdf_c( res.first, res.second);
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- POST-FIT AGREEMENT EVALUATION -----------------------");
        if(fGetChi2==1)
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- -------- STAT-ONLY --------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "--- REGION " + fName + ":");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  chi2        = " + std::to_string(fChi2val));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  ndof        = " + std::to_string(fNDF));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  probability = " + std::to_string(fChi2prob));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
    }
    
}

//__________________________________________________________________________________
//
TthPlot* Region::DrawPostFit(FitResults *fitRes,ofstream& pullTex,string opt){

    if(TtHFitter::PREFITONPOSTFIT){
        fPlotPostFit->h_tot_bkg_prefit = (TH1*)fPlotPreFit->GetTotBkg()->Clone("h_tot_bkg_prefit");
    }
  
    TthPlot *p = fPlotPostFit;
    if(fYmaxScale==0) p->SetYmaxScale(1.8);
    else              p->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) p->fYmax = fYmax;
    if(fYmin!=0) p->fYmin = fYmin;
    p->fRatioYmax = fRatioYmaxPostFit;
    p->fRatioYmin = fRatioYminPostFit;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") p->SetYaxis(fYTitle);
    //
    // For 4-top-style plots
    if(TtHFitter::OPTION["FourTopStyle"]>0){
        if(fRegionType==CONTROL) p->AddLabel("#font[62]{Control Region}");
        else if(fRegionType==SIGNAL) p->AddLabel("#font[62]{Signal Region}");
        else if(fRegionType==VALIDATION) p->AddLabel("#font[62]{Validation Region}");
        else p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
//         p->AddLabel("#font[52]{B-only fit}");
        p->AddLabel("#font[52]{Post-fit}");
        p->fLegendNColumns = TtHFitter::OPTION["LegendNColumns"];
    }
    //
    // old-style plots
    else{
        p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        if(TtHFitter::OPTION["NoPrePostFit"]==0) p->AddLabel("Post-Fit");
        if(TtHFitter::OPTION["LegendNColumns"]!=0) p->fLegendNColumns = TtHFitter::OPTION["LegendNColumns"];
    }
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);

    if(fBinLabels.size() && ((int)fBinLabels.size()==fNbins)) {
        for(int i_bin=0; i_bin<fNbins; i_bin++) {
            p->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
        }
    }
    
    //
    // 0) Create a new hist for each sample
    //
    TH1* hSmpNew[MAXsamples];
    for(int i=0;i<fNSamples;i++){
        hSmpNew[i] = (TH1*)fSampleHists[i]->fHist->Clone();
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if((!fSampleHists[i]->fSample->fUseMCStat && !fSampleHists[i]->fSample->fSeparateGammas) || fUseGammaPulls){
            for(int i_bin=0;i_bin<hSmpNew[i]->GetNbinsX()+2;i_bin++) hSmpNew[i]->SetBinError(i_bin,0.);
        }
    }
    
    //
    // 1) Propagates the post-fit NP values to the central value (pulls)
    //
    string systName;
    TH1* hNew;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if (fSampleHists[i]->fSample->fBuildPullTable>0){
            WriteDebugStatus("Region::DrawPostFit", "Propagating post-fit to Sample " + fSampleHists[i]->fSample->fTitle);
            TString sampleTex= fSampleHists[i]->fSample->fTexTitle;
            pullTex << "\\hline\n" << endl;
            pullTex << "{\\color{blue}{$\\rightarrow \\,$ "<< sampleTex << "}} & \\\\\n"<< endl;
        }
        hNew = (TH1*)hSmpNew[i]->Clone();
        for(int i_bin=1;i_bin<=hNew->GetNbinsX();i_bin++){
            double binContent0 = hSmpNew[i]->GetBinContent(i_bin);
            double mult_factor = GetMultFactors(fitRes, pullTex, i, i_bin, binContent0);
            
            //
            // Final computation
            //
            double binContentNew = binContent0*mult_factor;
            
            // 
            // gammas
            if(fUseGammaPulls && (fSampleHists[i]->fSample->fUseMCStat || fSampleHists[i]->fSample->fSeparateGammas)){
                // find the gamma for this bin of this distribution in the fit results
                std::string gammaName = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                if(fSampleHists[i]->fSample->fSeparateGammas)
                            gammaName = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                WriteDebugStatus("Region::DrawPostFit", "Looking for gamma " + gammaName);
                float gammaValue = fitRes->GetNuisParValue(gammaName);
                WriteDebugStatus("Region::DrawPostFit", "  -->  pull = " + std::to_string(gammaValue));
                // linear effect
                if(gammaValue>0) binContentNew *= gammaValue;
            }
            
            //
            // Setting to new values
            //
            hNew->SetBinContent(i_bin,binContentNew);
        }
        hSmpNew[i] = (TH1*)hNew->Clone();
        fSampleHists[i]->fHist_postFit = hSmpNew[i];
        hNew->~TH1();
    }
    delete hNew;
    
    //
    // 2) Scale all samples by norm factors
    //    Done after the propagation of the NP (avoids nans due to "0" value of some NormFactors)
    //    Seems consistent with application in Roostats
    //
    string nfName;
    float nfValue;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        for(int i_norm=0;i_norm<fSampleHists[i]->fNNorm;i_norm++){
            NormFactor *nf = fSampleHists[i]->fNormFactors[i_norm];
            nfName = nf->fName;
//             nfName = nf->fNuisanceParameter; // not implemeted yet...
            if(nf->fConst) nfValue = nf->fNominal;
            else           nfValue = fitRes->GetNuisParValue(TtHFitter::NPMAP[nfName]);
//             if(nfValue==0) nfValue = 0.0001;  // FIXME
//             if(nfName=="SigXsecOverSM" && nfValue==0) nfValue = 0.0001;   // FIXME
//             if(nfName=="SigXsecOverSM" && nfValue==0) nfValue = 1;   // FIXME
//             if(nfName=="SigXsecOverSM") nfValue = 1;   // FIXME
            //
            // if this norm factor is a morphing one
            if(nfName.find("morph_")!=string::npos){
                std::string formula = TtHFitter::SYSTMAP[nfName];
                std::string name = TtHFitter::NPMAP[nfName];
                formula = ReplaceString(formula,name,"x");
                TF1* f_morph = new TF1("f_morph",formula.c_str(),nf->fMin,nf->fMax);
                float scale = f_morph->Eval(nfValue);
                hSmpNew[i]->Scale(scale);
            }
            else{
                hSmpNew[i]->Scale(nfValue);
            }
        }
    }
    
    
    //
    // 2)b) Scale all samples by shape factors factors
    //

    string sfName;
    string sfNameBin;
    float sfValue;
    double binContentSFNew = 0;
    int iBinSF = 0;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        for(int i_shape=0;i_shape<fSampleHists[i]->fNShape;i_shape++){
            ShapeFactor *sf = fSampleHists[i]->fShapeFactors[i_shape];
            sfName = sf->fName;
            // loop over bins
            // there should be a NP per bin in the fit file
            // already checked by GetNuisParValue()
            // the shape factor naming used i_bin - 1 for the first bin
            for(int i_bin = 1; i_bin <= hSmpNew[i]->GetNbinsX(); i_bin++){
                iBinSF = i_bin - 1;
                sfNameBin = sfName + "_bin_" + std::to_string(iBinSF);
                if(sf->fConst) sfValue = sf->fNominal; 
                else           sfValue = fitRes->GetNuisParValue(sfNameBin);
                // scale bin content by shape factor
                binContentSFNew = hSmpNew[i]->GetBinContent(i_bin) * sfValue;
                // Setting to new value
                hSmpNew[i]->SetBinContent(i_bin, binContentSFNew);
            }
        }
    }
    
    //
    // 3) Add the new Sig and Bkg to plot
    //
    string title;
    TH1* hBkgNew[MAXsamples];
    TH1* hSigNew[MAXsamples];
    for(int i=0, i_bkg=0, i_sig=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::BACKGROUND){
            hBkgNew[i_bkg] = hSmpNew[i];
            i_bkg++;
        }
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL){
            hSigNew[i_sig] = hSmpNew[i];
            i_sig++;
        }
    }
    if(fHasData && opt.find("blind")==string::npos) p->SetData(fData->fHist,fData->fSample->fTitle);
    for(int i=0;i<fNSig;i++){
        title = fSig[i]->fSample->fTitle;
        if(fSig[i]->fSample->fGroup != "") title = fSig[i]->fSample->fGroup;
        if(TtHFitter::SHOWSTACKSIG)    p->AddSignal(    hSigNew[i],title);
//         if(TtHFitter::SHOWNORMSIG)     p->AddNormSignal(hSigNew[i],title);
//         if(TtHFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(hSigNew[i],title);
        if(TtHFitter::SHOWNORMSIG){
            if( (TtHFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
             || !TtHFitter::OPTION["NormSigSRonly"] )
                p->AddNormSignal(hSigNew[i],title);
        }
        else{
            if(TtHFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(hSigNew[i],title);
        }
        if(TtHFitter::SHOWOVERLAYSIG)  p->AddOverSignal(hSigNew[i],title);
    }
    for(int i=0;i<fNBkg;i++){
        title = fBkg[i]->fSample->fTitle;
        if(fBkg[i]->fSample->fGroup != "") title = fBkg[i]->fSample->fGroup;
        p->AddBackground(hBkgNew[i],title);
    }
    
    //
    // 4) Build post-fit error band
    //    Build total histogram (fTot_postFit)
    //
    int j = 0;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !TtHFitter::SHOWSTACKSIG) continue;
        if(j==0) fTot_postFit = (TH1*)hSmpNew[i]->Clone("h_tot_postFit");
        else fTot_postFit->Add(hSmpNew[i]);
        j++;
    }
    
    //
    // MCstat set to 0 if disabled
    //
    if(!fUseStatErr || fUseGammaPulls){
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
            fTot_postFit->SetBinError(i_bin,0);
        }
    }
    
    p->SetTotBkg(fTot_postFit);
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);
    
    //
    // blinding bins
    //
    if(fBlindingThreshold>=0){
        p->SetBinBlinding(true,fBlindingThreshold);
        if(fKeepPrefitBlindedBins && fBlindedBins!=0x0) p->SetBinBlinding(true,fBlindedBins);
    }
    p->BlindData();
    
    //
    // Build error band
    //
    BuildPostFitErrorHist(fitRes);
    
    //
    // Print chi2 info
    //
    if(fGetChi2 && TtHFitter::SHOWCHI2)  p->SetChi2KS(fChi2prob,-1,fChi2val,fNDF);
    
    //
    // 5) Finishes configuration of TthPlot objects
    //
    p->SetTotBkgAsym(fErr_postFit);
    p->fATLASlabel = fATLASlabel;
    if(fLogScale) opt += " log";
    
    p->Draw(opt);
    
    //
    // Print bin content and errors
    //
    WriteDebugStatus("Region::DrawPostFit", "--------------------");
    WriteDebugStatus("Region::DrawPostFit", "Final bin contents");
    WriteDebugStatus("Region::DrawPostFit", "--------------------");
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
        WriteDebugStatus("Region::DrawPostFit", std::to_string(i_bin) + ":\t" + std::to_string(fTot_postFit->GetBinContent(i_bin)) + " +" +
        std::to_string( fErr_postFit->GetErrorYhigh(i_bin-1)) + " -" + std::to_string(fErr_postFit->GetErrorYlow(i_bin-1)));
    }
    
    //
    // Save in a root file...
    //
    gSystem->mkdir((fFitName+"/Histograms").c_str());
    WriteInfoStatus("Region::DrawPostFit", "Writing file " + fFitName+"/Histograms/"+fName+fSuffix+"_postFit.root");
    TFile *f = new TFile((fFitName+"/Histograms/"+fName+fSuffix+"_postFit.root").c_str(),"RECREATE");
    f->cd();
    fErr_postFit->Write("",TObject::kOverwrite);
    fTot_postFit->Write("",TObject::kOverwrite);
    if( fPlotPostFit->h_tot_bkg_prefit)
      fPlotPostFit->h_tot_bkg_prefit->Write("",TObject::kOverwrite);
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        if(fTotUp_postFit[i_syst])   fTotUp_postFit[i_syst]  ->Write("",TObject::kOverwrite);
        if(fTotDown_postFit[i_syst]) fTotDown_postFit[i_syst]->Write("",TObject::kOverwrite);
    }
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType == Sample::DATA){
            fSampleHists[i]->fHist->Write(Form("h_%s",fSampleHists[i]->fName.c_str()),TObject::kOverwrite);
            continue;
        }
        if(fSampleHists[i]->fHist_postFit){
            fSampleHists[i]->fHist_postFit->Write(Form("h_%s_postFit",fSampleHists[i]->fName.c_str()),TObject::kOverwrite);
            for(unsigned int i_syst=0;i_syst<fSampleHists[i]->fSyst.size();i_syst++){
                if(fSampleHists[i]->fSyst[i_syst]) 
                    if(fSampleHists[i]->fSyst[i_syst]->fHistUp_postFit)
                        fSampleHists[i]->fSyst[i_syst]->fHistUp_postFit  ->Write(
                          Form("h_%s_%s_Up_postFit",fSampleHists[i]->fName.c_str(),fSampleHists[i]->fSyst[i_syst]->fName.c_str()), TObject::kOverwrite);
                if(fSampleHists[i]->fSyst[i_syst])
                    if(fSampleHists[i]->fSyst[i_syst]->fHistDown_postFit)
                        fSampleHists[i]->fSyst[i_syst]->fHistDown_postFit->Write(
                          Form("h_%s_%s_Down_postFit",fSampleHists[i]->fName.c_str(),fSampleHists[i]->fSyst[i_syst]->fName.c_str()),TObject::kOverwrite);
            }
        }
    }
    f->Close();
    delete f;
    //
    return p;
}

//__________________________________________________________________________________
//
void Region::AddSelection(string selection){
    if(selection=="") return;
    if(fSelection=="1" || fSelection=="") fSelection = selection;
    else fSelection += " && "+selection;
}

//__________________________________________________________________________________
//
void Region::AddMCweight(string weight){
    if(weight=="") return;
    if(fMCweight=="1" || fMCweight=="") fMCweight = weight;
    else fMCweight += " * "+weight;
}

//__________________________________________________________________________________
//
void Region::SetVariable(string variable,int nbin,float xmin,float xmax,string corrVar1,string corrVar2){
    fVariable = variable;
    fCorrVar1 = corrVar1;
    fCorrVar2 = corrVar2;
    fNbins = nbin;
    fXmin = xmin;
    fXmax = xmax;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeVariable(string variable,string sample){
    fAlternativeVariables[sample] = variable;
}


//__________________________________________________________________________________
//
bool Region::UseAlternativeVariable(string sample){
    if (fAlternativeVariables.find(sample)==fAlternativeVariables.end())
        return false;
    else
        return true;
}

//__________________________________________________________________________________
//
void Region::SetHistoName(string name){
//     fHistoName = name;
    fHistoNames.clear();
    fHistoNames.push_back(name);
}

//__________________________________________________________________________________
//
void Region::SetVariableTitle(string name){
    fVariableTitle = name;
}

//__________________________________________________________________________________
//
void Region::SetLabel(string label,string shortLabel){
    fLabel = label;
    if(shortLabel=="") fShortLabel = label;
    else fShortLabel = shortLabel;
}

//__________________________________________________________________________________
//
void Region::Print(){
    WriteInfoStatus("Region::Print", "    Region: " + fName);
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        fSampleHists[i_smp]->Print();
    }
}

//__________________________________________________________________________________
//
void Region::PrintSystTable(FitResults *fitRes, string opt){
    bool isPostFit  = false; if(opt.find("post")!=string::npos)     isPostFit  = true;
    bool doClean    = false; if(opt.find("clean")!=string::npos)    doClean    = true;
    bool doCategory = false; if(opt.find("category")!=string::npos) doCategory = true;
    //
    ofstream out;
    ofstream texout;
    ofstream out_cat;
    ofstream texout_cat;
    gSystem->mkdir(fFitName.c_str());
    gSystem->mkdir((fFitName+"/Tables").c_str());
    if(isPostFit){
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.tex").c_str());
        }
    }
    else{
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.tex").c_str());
        }
    }
    Sample *s = 0x0;
    SampleHist *sh = 0x0;
    SystematicHist *syh = 0x0;
    
    
    out << " | ";
    texout << "\\documentclass[10pt]{article}" << endl;
    texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
    texout << "\\begin{document}" << endl;
    texout << "\\begin{table}[htbp]" << endl;
    texout << "\\begin{center}" << endl;
    texout << "\\begin{tabular}{|c" ;

    if(doCategory){
        out_cat << " | ";
        texout_cat << "\\documentclass[10pt]{article}" << endl;
        texout_cat << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
        texout_cat << "\\begin{document}" << endl;
        texout_cat << "\\begin{table}[htbp]" << endl;
        texout_cat << "\\begin{center}" << endl;
        texout_cat << "\\begin{tabular}{|c" ;
    }

    
    float Ncol = 2.;
    for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
        sh = fSampleHists[i_smp];
        s = sh->fSample;
        if(s->fType==Sample::DATA) continue;
        if(s->fType==Sample::GHOST) continue;
        texout << "|c";
        if(doCategory) texout_cat << "|c";
        Ncol+=1;
    }
    texout << "|}" << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        texout_cat << "|}" << endl;
        texout_cat << "\\hline " << endl;
    }
    //
    float i_col = 1;
    for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
        sh = fSampleHists[i_smp];
        s = sh->fSample;
        if(s->fType==Sample::DATA) continue;
        if(s->fType==Sample::GHOST) continue;
        std::string title = s->fTitle;
        if(s->fTexTitle!="") title = s->fTexTitle;
        out << "      | " << s->fTitle;
        texout << "      & " << title;
        if(doCategory){
            out_cat << "      | " << s->fTitle;
            texout_cat << "      & " << title;
        }
        i_col+=1;
    }
    out << " |" << endl;
    texout << " \\\\ " << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        out_cat << " |" << endl;
        texout_cat << " \\\\ " << endl;
        texout_cat << "\\hline " << endl;
    }

    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        if(TtHFitter::SYSTMAP[fSystNames[i_syst]]!="") out << " | " << TtHFitter::SYSTMAP[fSystNames[i_syst]];
        else                                           out << " | " << fSystNames[i_syst];
        if(TtHFitter::SYSTTEX[fSystNames[i_syst]]!="")      texout << "  " << TtHFitter::SYSTTEX[fSystNames[i_syst]];
        else if(TtHFitter::SYSTMAP[fSystNames[i_syst]]!=""){
            string fixedTitle = TtHFitter::SYSTMAP[fSystNames[i_syst]];
            fixedTitle = ReplaceString(fixedTitle,"#geq","$\\geq$");
            texout << "  " << fixedTitle;
        }
        else                                                texout << " " << fSystNames[i_syst];
        for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
            sh = fSampleHists[i_smp];
            s = sh->fSample;
            if(s->fType==Sample::DATA) continue;
            if(s->fType==Sample::GHOST) continue;
            syh = sh->GetSystematic(fSystNames[i_syst]);
            if(syh==0x0){
                out << " |    nan   ";
                texout << " &    nan   ";
            }
            else{
              float normUp;
              float normDown;
              if(isPostFit){
                normUp   = (syh->fHistUp_postFit->Integral()   - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
                normDown = (syh->fHistDown_postFit->Integral() - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
              }
              else{
                normUp = syh->fNormUp;
                normDown = syh->fNormDown;
              }

              out << " | " << normUp;
              texout << setprecision(3) << " & " << normUp;
              out << " / " << normDown;
              texout << setprecision(3) << " / " << normDown;
              //                 pt[i_smp]->AddText(Form("%.2f / %.2f",syh->fNormUp,syh->fNormDown) );
            }
        }

        out << " |" << endl;
        texout << " \\\\ " << endl;
    }


    if(doCategory){
        //--- Systematic tables per category:
        //--- Get systematic names per category and sample
        std::set<std::string> category_names;
        std::map<std::string, std::map<std::string, std::vector<std::string> > > category_syst_names;
        for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
            sh = fSampleHists[i_smp];
            s = sh->fSample;
            for(int i_samplesyst=0; i_samplesyst<(int)s->fSystematics.size();i_samplesyst++){
                if(s->fType==Sample::DATA) continue;
                if(s->fType==Sample::GHOST) continue;
                std::string category = s->fSystematics[i_samplesyst]->fCategory;
                if (category!=""){
                    category_names.insert(category);
                    std::vector<std::string> systePerSample;
                    // Check for systematics existing not in all regions...
                    if (std::find(fSystNames.begin(), fSystNames.end(), s->fSystematics.at(i_samplesyst)->fName)!=fSystNames.end()){
                        category_syst_names[s->fName][category].push_back(s->fSystematics.at(i_samplesyst)->fName);
                    }
                }
            }
        }
        //--- Loop over categories and samples, to get systematic effects using BuildTotError function
        for (auto category : category_names){
            out_cat << " | " << category; 
            texout_cat << " " << category;
            for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
                sh = fSampleHists[i_smp];
                s = sh->fSample;
                if(s->fType==Sample::DATA) continue;
                if(s->fType==Sample::GHOST) continue;

                std::vector<TH1*> category_histo_up;
                std::vector<TH1*> category_histo_down;
                std::vector<std::string> sample_syste = category_syst_names[s->fName][category];
                for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
                    if(!sh->HasSyst(fSystNames[i_syst]))
                      continue;

                    syh = sh->GetSystematic(fSystNames[i_syst]);
                    if (std::find(sample_syste.begin(), sample_syste.end(), fSystNames.at(i_syst)) == sample_syste.end())
                      continue;

                    if(isPostFit){
                        TH1 *h_up = (TH1*) syh->fHistUp_postFit;
                        h_up->SetDirectory(0);
                        TH1 *h_down = (TH1*) syh->fHistDown_postFit;
                        h_down->SetDirectory(0);
                        category_histo_up.push_back(h_up);
                        category_histo_down.push_back(h_down);
                    }
                    else{
                        TH1 *h_up = (TH1*) syh->fHistUp;
                        h_up->SetDirectory(0);
                        TH1 *h_down = (TH1*) syh->fHistDown;
                        h_down->SetDirectory(0);
                        category_histo_up.push_back(h_up);
                        category_histo_down.push_back(h_down);
                    }
                }

                TGraphAsymmErrors *g_err;
                double err = 0.;
                if (isPostFit){
                    g_err = BuildTotError(sh->fHist_postFit, category_histo_up, category_histo_down, category_syst_names[s->fName][category], fitRes->fCorrMatrix);
                    if (category_histo_up.size()>0 && sh->fHist_postFit->Integral()>0.){
                        for (int ibin=1; ibin<sh->fHist_postFit->GetNbinsX()+1; ibin++){
                            if (pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist_postFit->GetBinError(ibin), 2)>=0.){ // dummy check
                                //Need to substract the statistical unc.
                                err += (sqrt(pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist_postFit->GetBinError(ibin), 2)))/sh->fHist_postFit->Integral();
                            }
                        }
                    }
                }
                else{
                    g_err = BuildTotError(sh->fHist, category_histo_up, category_histo_down, category_syst_names[s->fName][category]);
                    if (category_histo_up.size()>0 && sh->fHist->Integral()>0.){
                        for (int ibin=1; ibin<sh->fHist->GetNbinsX()+1; ibin++){
                            if (pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist->GetBinError(ibin), 2)>=0.){
                                //Need to substract the statistical unc.
                                err += (sqrt(pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist->GetBinError(ibin), 2)))/sh->fHist->Integral();
                            }
                        }
                    }
                }
                out_cat << setprecision(3) << " | " << err;
                texout_cat << setprecision(3) << " & " << err;
            }
            out_cat << " |" << endl;
            texout_cat << " \\\\ " << endl;
        }
    }

    texout << "\\hline " << endl;
    texout << "\\end{tabular} " << endl;
    texout << "\\caption{Relative effect of each systematic on the yields.} " << endl;
    texout << "\\end{center} " << endl;
    texout << "\\end{table} " << endl;
    texout << "\\end{document}" << endl;

    if(doCategory){
        texout_cat << "\\hline " << endl;
        texout_cat << "\\end{tabular} " << endl;
        texout_cat << "\\caption{Realtive effect of each group of systematics on the yields.} " << endl;
        texout_cat << "\\end{center} " << endl;
        texout_cat << "\\end{table} " << endl;
        texout_cat << "\\end{document}" << endl;
    }

    if(doClean){
        std::string shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
        if(isPostFit) shellcommand += "_postFit";
        shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
        shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst_clean.tex";
        gSystem->Exec(shellcommand.c_str());
        if(doCategory){
            shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
            shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += "_clean.tex";
            gSystem->Exec(shellcommand.c_str());
        }
    }
}

// //__________________________________________________________________________________
// //
// void Region::PrintSystTable(string opt){
//     bool isPostFit = false;
//     if(opt.find("post")!=string::npos){
//         isPostFit = true;
//     }
//     //
//     ofstream out;
//     ofstream texout;
//     gSystem->mkdir(fFitName.c_str());
//     gSystem->mkdir((fFitName+"/Tables").c_str());
//     if(isPostFit){
//         out.open((fFitName+"/Tables/"+fName+"_syst_postFit.txt").c_str());
//         texout.open((fFitName+"/Tables/"+fName+"_syst_postFit.tex").c_str());
//     }
//     else{
//         out.open((fFitName+"/Tables/"+fName+"_syst.txt").c_str());
//         texout.open((fFitName+"/Tables/"+fName+"_syst.tex").c_str());
//     }
//     Sample *s = 0x0;
//     SampleHist *sh = 0x0;
//     SystematicHist *syh = 0x0;
//     out << " | ";
//     
//     texout << "\\documentclass[10pt]{article}" << endl;
//     texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
//     texout << "\\begin{document}" << endl;
//     //texout << "\\small" << endl;
//     
//     texout << "\\begin{table}[htbp]" << endl;
//     texout << "\\begin{center}" << endl;
//     //texout << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << endl;
//     texout << "\\begin{tabular}{|c" ;
//     for (int i =0; i<fSampleHists.size(); i++)
//       {
//     texout << "|c";      
//       }
//     texout << "|}" << endl;
// 
//     texout << "\\hline " << endl;
//     
//     float Ncol = 2.;
//     for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
//         sh = fSampleHists[i_smp];
//         s = sh->fSample;
//         if(s->fType==Sample::DATA) continue;
//         if(s->fType==Sample::GHOST) continue;
//         Ncol+=1;
//     }
//     // plot a table with ROOT ;)
//     int width = (Ncol)*120;
//     int height = fSystNames.size()*25;
//     //
// //     TCanvas *c = new TCanvas("c","c",width,height);
// //     TPaveText *pt0;
// //     TPaveText *pt[MAXsamples];
// //     pt0 = new TPaveText(0,0,2./Ncol,1);
// //     pt0->SetTextSize(gStyle->GetTextSize());
// //     pt0->SetFillStyle(0);
//     //
//     float i_col = 1;
//     for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
// //         pt[i_smp] = new TPaveText((i_col+1.)/Ncol,0,(i_col+2.)/Ncol,1); // ,"blNDC"
// //         pt[i_smp]->SetTextSize(gStyle->GetTextSize());
// //         pt[i_smp]->SetFillStyle(0);
//         sh = fSampleHists[i_smp];
//         s = sh->fSample;
//         if(s->fType==Sample::DATA) continue;
//         if(s->fType==Sample::GHOST) continue;
//         out << "      | " << s->fTitle;
//         texout << "      & " << s->fTitle;
// //         pt[i_smp]->AddText(s->fTitle.c_str());
//         i_col+=1;
//     }
//     out << " |" << endl;
//     texout << " \\\\ " << endl;
//     texout << "\\hline " << endl;
// 
// //     pt0->AddText(" ");
//     //
// //     TPave *b[100];
// //     int i_gray = 0;
//     for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
// //         if(i_syst%2==0) {
// //             b[i_gray] = new TPave(0,(1.*fSystNames.size()-i_syst)/(fSystNames.size()+1.),1,(1.*fSystNames.size()-i_syst-1)/(fSystNames.size()+1.),0);
// //             b[i_gray]->SetFillColor(kGray);
// //             b[i_gray]->Draw("NB");
// //             i_gray++;
// //         }
//         //
//         if(TtHFitter::SYSTMAP[fSystNames[i_syst]]!="") {out << " | " << TtHFitter::SYSTMAP[fSystNames[i_syst]]; texout << "  " << TtHFitter::SYSTMAP[fSystNames[i_syst]];}
//         else                                           {out << " | " << fSystNames[i_syst]; texout << " " << fSystNames[i_syst];}
// //         if(TtHFitter::SYSTMAP[fSystNames[i_syst]]!="") pt0->AddText(TtHFitter::SYSTMAP[fSystNames[i_syst]].c_str());
// //         else                                           pt0->AddText(fSystNames[i_syst].c_str());
// //         if(i_syst==0) pt0->AddLine(0,(1.*fSystNames.size()-i_syst)/(fSystNames.size()+1.),1,(1.*fSystNames.size()-i_syst)/(fSystNames.size()+1.));
//         for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
// //             if(i_syst==0) pt[i_smp]->AddLine(0,(1.*fSystNames.size()-i_syst)/(fSystNames.size()+1.),1,(1.*fSystNames.size()-i_syst)/(fSystNames.size()+1.));
//             sh = fSampleHists[i_smp];
//             s = sh->fSample;
//             if(s->fType==Sample::DATA) continue;
//             if(s->fType==Sample::GHOST) continue;
//             syh = sh->GetSystematic(fSystNames[i_syst]);
//             if(syh==0x0){
//                 out << " |    nan   ";
//                 texout << " &    nan   ";
// //                 pt[i_smp]->AddText(" - ");
//             }
// //             sh = GetSampleHist(fSamples[i_smp]->fName);
//             else{
//                 float normUp;
//                 float normDown;
//                 if(isPostFit){
//                     normUp   = (syh->fHistUp_postFit->Integral()   - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
//                     normDown = (syh->fHistDown_postFit->Integral() - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
//                 }
//                 else{
//                     normUp = syh->fNormUp;
//                     normDown = syh->fNormDown;
//                 }
//                 out << " | " << normUp;
//                 texout << setprecision(3) << " & " << normUp;
//                 out << " / " << normDown;
//                 texout << setprecision(3) << " / " << normDown;
// //                 pt[i_smp]->AddText(Form("%.2f / %.2f",syh->fNormUp,syh->fNormDown) );
//             }
// //             if(i_syst==(int)fSystNames.size()-1) pt[i_smp]->Draw("NB");
//         }
//         out << " |" << endl;
//         texout << " \\\\ " << endl;
// //         if(i_syst==(int)fSystNames.size()-1) pt0->Draw("NB");
//     }
// 
//     texout << "\\hline " << endl;
// 
//     texout << "\\end{tabular} " << endl;
//     texout << "\\caption{Yields of the analysis} " << endl;
//     texout << "\\end{center} " << endl;
//     texout << "\\end{table} " << endl;
//     texout << "\\end{document}" << endl;
// 
//     //TString shellcommand = "cat l3_syst.tex|sed -e \"s/\\#/ /g\" > truc.tex";
//     //texout.open((fFitName+"/Tables/"+fName+"_syst.tex").c_str());
// 
//     TString shellcommand = "cat "+fFitName+"/Tables/"+fName+"_syst.tex|sed -e \"s/\\#/ /g\" > "+fFitName+"/Tables/"+fName+"_syst_clean.tex";
//     //TString shellcommand = "cat "+fFitName+"/Tables/"+fName+"_syst.tex|sed -e \"s/t\\#bar\\{t\\}/t\\$\\bar\\{t\\}\\$/g\" > "+fFitName+"/Tables/"+fName+"_syst_truc.tex";
//     gSystem->Exec(shellcommand);
//     
//     //cat l3_syst.tex|sed -e "s/\#/ /g" > truc.tex
//     
// //     c->SaveAs((fFitName+"/Tables/"+fName+"_syst.pdf").c_str());
// }

// --------------- Functions --------------- //

float GetDeltaN(float alpha, float Iz, float Ip, float Imi, int intCode){
    // protection against negative values
    if(Ip<=0)  Ip  = 0.00000001*Iz;
    if(Imi<=0) Imi = 0.00000001*Iz;
    
    float deltaN = 0.;
    if(alpha>0)      deltaN = Ip;
    else if(alpha<0) deltaN = Imi;
    else             return 1.; 

    if(intCode==4){
        
        //////////////////////////////////////////////////////////////
        // NORMALISATION
        // =============
        // Exponential-polynomial
        // Equation solved with Mathematica
        //////////////////////////////////////////////////////////////
        
        if(TMath::Abs(alpha)>1){
            deltaN /= Iz; // divide h_tmp by the nominal
            deltaN = pow( deltaN, TMath::Abs(alpha) );  // d -> d^(|a|)
        } else {
            // polinomial: equations solved with Mathematica
            float a1 = -(15*Imi - 15*Ip - 7*Imi*TMath::Log(Imi/Iz) + Imi*pow(TMath::Log(Imi/Iz),2) + 7*Ip*TMath::Log(Ip/Iz) - Ip*pow(TMath::Log(Ip/Iz),2))/(16.*Iz);
            float a2 = -3 + (3*Imi)/(2.*Iz) + (3*Ip)/(2.*Iz) - (9*Imi*TMath::Log(Imi/Iz))/(16.*Iz) + (Imi*pow(TMath::Log(Imi/Iz),2))/(16.*Iz) -
            (9*Ip*TMath::Log(Ip/Iz))/(16.*Iz) + (Ip*pow(TMath::Log(Ip/Iz),2))/(16.*Iz);
            float a3 = (5*Imi)/(8.*Iz) - (5*Ip)/(8.*Iz) - (5*Imi*TMath::Log(Imi/Iz))/(8.*Iz) + (Imi*pow(TMath::Log(Imi/Iz),2))/(8.*Iz) + (5*Ip*TMath::Log(Ip/Iz))/(8.*Iz) -
            (Ip*pow(TMath::Log(Ip/Iz),2))/(8.*Iz);
            float a4 = 3 - (3*Imi)/(2.*Iz) - (3*Ip)/(2.*Iz) + (7*Imi*TMath::Log(Imi/Iz))/(8.*Iz) -
            (Imi*pow(TMath::Log(Imi/Iz),2))/(8.*Iz) + (7*Ip*TMath::Log(Ip/Iz))/(8.*Iz) - (Ip*pow(TMath::Log(Ip/Iz),2))/(8.*Iz);
            float a5 = (-3*Imi)/(16.*Iz) + (3*Ip)/(16.*Iz) + (3*Imi*TMath::Log(Imi/Iz))/(16.*Iz) - (Imi*pow(TMath::Log(Imi/Iz),2))/(16.*Iz) -
            (3*Ip*TMath::Log(Ip/Iz))/(16.*Iz) + (Ip*pow(TMath::Log(Ip/Iz),2))/(16.*Iz);
            float a6 = -1 + Imi/(2.*Iz) + Ip/(2.*Iz) - (5*Imi*TMath::Log(Imi/Iz))/(16.*Iz) + (Imi*pow(TMath::Log(Imi/Iz),2))/(16.*Iz) - (5*Ip*TMath::Log(Ip/Iz))/(16.*Iz) +
            (Ip*pow(TMath::Log(Ip/Iz),2))/(16.*Iz);
            float a = alpha;
            deltaN = 1 + a1*a + a2*a*a + a3*a*a*a + a4*a*a*a*a + a5*a*a*a*a*a + a6*a*a*a*a*a*a;
        }
        
    }
    else if (intCode==0) {
        
        //////////////////////////////////////////////////////////////
        // SHAPES
        // ======
        // Linear-polynomial
        // Extracted and adapted from RooFit
        //////////////////////////////////////////////////////////////
        
        if(TMath::Abs(alpha)>1){
            deltaN = 1. + TMath::Abs(alpha)*(deltaN - Iz)/Iz;
        } else {
            double eps_plus = Ip - Iz;
            double eps_minus = Iz - Imi;
            double S = 0.5 * (eps_plus + eps_minus);
            double A = 0.0625 * (eps_plus - eps_minus);
            double val = Iz + alpha * (S + alpha * A * ( 15 + alpha * alpha * (-10 + alpha * alpha * 3  ) ) );
            if (val < 0) val = 0.;
            if(Iz == Iz && Iz>0) deltaN = val/Iz;
            else deltaN = 0.;
        }
    }
    
    if(deltaN!=deltaN) deltaN = 1;  // to avoid nan
    if(deltaN<=0) deltaN = 0; //protection against negative values (can happen for linear extrapolation)

    return deltaN;
}

//___________________________________________________________
//
std::map < int , double > GetDeltaNForUncertainties(float alpha, float alpha_errUp, float alpha_errDown, float Iz, float Ip, float Imi, int intCode){
    double nominal = GetDeltaN(alpha, Iz, Ip, Imi, intCode);
    double up = GetDeltaN(alpha+alpha_errUp, Iz, Ip, Imi, intCode);
    double down = GetDeltaN(alpha+alpha_errDown, Iz, Ip, Imi, intCode);
    return {{-1,down},{0,nominal},{1,up}};
}


//--------------- ~ ---------------

// function to get pre/post-fit agreement
std::pair<double,int> GetChi2Test( TH1* h_data, TH1* h_nominal, std::vector< TH1* > h_up, std::vector< TH1* > h_down, std::vector< string > fSystNames, CorrelationMatrix *matrix ){
    unsigned int nbins = h_nominal->GetNbinsX();
    unsigned int nsyst = fSystNames.size();
    int ndf = 0;
    for(unsigned int i=0;i<nbins;i++){
        double ydata_i = h_data->GetBinContent(i+1);
        if(ydata_i<0) continue; // skip dropped / blinded bins
        ndf ++;
    }
    TMatrixD C(ndf,ndf);
    //
    int ibin = 0;
    int jbin = 0;
    for(unsigned int i=0;i<nbins;i++){
        double ydata_i = h_data->GetBinContent(i+1);
        if(ydata_i<0) continue; // skip dropped / blinded bins
        double ynom_i = h_nominal->GetBinContent(i+1);
        jbin = 0;
        for(unsigned int j=0;j<nbins;j++){
            double sum = 0.;
            double ydata_j = h_data->GetBinContent(j+1);
            if(ydata_j<0) continue; // skip dropped / blinded bins
            double ynom_j = h_nominal->GetBinContent(j+1);
            for(unsigned int n=0;n<nsyst;n++){
                double ysyst_i_n = h_up[n]->GetBinContent(i+1);
                for(unsigned int m=0;m<nsyst;m++){
                    double ysyst_j_m = h_up[m]->GetBinContent(j+1);
                    // more than Bill's suggestion: add correlation between systematics!!
                    double corr = 0.;
                    if(n==m) corr = 1.;
                    else if(matrix!=0x0) corr = matrix->GetCorrelation(fSystNames[n],fSystNames[m]);
                    else continue;
                    sum += (ysyst_i_n-ynom_i) * corr * (ysyst_j_m-ynom_j);
                }
            }
            if(i==j && ynom_i>0) sum += ynom_i;  // add stat uncertainty to diagonal
            if(i==j && ynom_i>0) sum += pow(h_nominal->GetBinError(i+1),2); // add MC stat as well
            C[ibin][jbin] = sum;
            jbin++;
        }
        ibin++;
    }
    //
    if(TtHFitter::DEBUGLEVEL) C.Print();
    //
    // Invert the matrix
    C.Invert();
    if(TtHFitter::DEBUGLEVEL) C.Print();
    //
    double chi2 = 0.;
    ibin = 0;
    jbin = 0;
    for(unsigned int i=0;i<nbins;i++){
        double ydata_i = h_data->GetBinContent(i+1);
        if(ydata_i<0) continue; // skip dropped / blinded bins
        double ynom_i = h_nominal->GetBinContent(i+1);
        jbin = 0;
        for(unsigned int j=0;j<nbins;j++){
            double ydata_j = h_data->GetBinContent(j+1);
            if(ydata_j<0) continue; // skip dropped / blinded bins
            double ynom_j = h_nominal->GetBinContent(j+1);
            chi2 += (ydata_i - ynom_i)*C[ibin][jbin]*(ydata_j-ynom_j);
            jbin++;
        }
        ibin++;
    }
    //
    return std::make_pair(chi2,ndf);
}

//--------------- ~ ---------------

// function to bouild total error band from:
// - a nominal histo (tot exepcted)
// - syst variation histos (eventually already scaled by post-fit pulls)
// - correlation matrix
// Note: if matrix = 0x0 => no correlation, i.e. matrix = 1 (used for pre-fit, or to neglect correlation)
TGraphAsymmErrors* BuildTotError( TH1* h_nominal, std::vector< TH1* > h_up, std::vector< TH1* > h_down, std::vector< string > fSystNames, CorrelationMatrix *matrix ){
    TGraphAsymmErrors *g_totErr = new TGraphAsymmErrors( h_nominal );
    float finalErrPlus(0.);
    float finalErrMinus(0.);
    float corr(0.);
    float errUp_i(0.),   errUp_j(0.);
    float errDown_i(0.), errDown_j(0.);
    //
    // - loop on bins
    for(int i_bin=1;i_bin<h_nominal->GetNbinsX()+1;i_bin++){
        finalErrPlus = 0;
        finalErrMinus = 0;
        corr = 0;
        // yieldNominal = h_nominal->GetBinContent(i_bin);
        // - loop on the syst, two by two, to include the correlations
//         std::vector<string> systNames_unique;
        for(unsigned int i_syst=0;i_syst<fSystNames.size();i_syst++){
//             if (std::find(systNames_unique.begin(), systNames_unique.end(), fSystNames[i_syst]) == systNames_unique.end())
//                 systNames_unique.push_back(fSystNames[i_syst]);
//             else continue;
            for(unsigned int j_syst=0;j_syst<fSystNames.size();j_syst++){
//                 if (fSystNames[i_syst]==fSystNames[j_syst] && i_syst!=j_syst) continue;
                if(matrix!=0x0){
                    corr = matrix->GetCorrelation(fSystNames[i_syst],fSystNames[j_syst]);
                }
                else{
                    if(fSystNames[i_syst]==fSystNames[j_syst]) corr = 1.;
                    else               corr = 0.;
                }
                errUp_i   = h_up[i_syst]  ->GetBinContent(i_bin);// - yieldNominal;
                errDown_i = h_down[i_syst]->GetBinContent(i_bin);// - yieldNominal;
                errUp_j   = h_up[j_syst]  ->GetBinContent(i_bin);// - yieldNominal;
                errDown_j = h_down[j_syst]->GetBinContent(i_bin);// - yieldNominal;
                
                //
                // Symmetrize (seems to be done in Roostats ??)
                //
                double err_i = (errUp_i - errDown_i)/2.;
                double err_j = (errUp_j - errDown_j)/2.;
                
                //
                // Compute the + and - variations
                //
                finalErrPlus  += err_i * err_j * corr;
                finalErrMinus += err_i * err_j * corr;
            }
        }
        // add stat uncertainty, which should have been stored as orignal bin errors in the h_nominal (if fUseStatErr is true)
        finalErrPlus  += pow( h_nominal->GetBinError(i_bin), 2 );
        finalErrMinus += pow( h_nominal->GetBinError(i_bin), 2 );
        
        g_totErr->SetPointEYhigh(i_bin-1,sqrt(TMath::Abs(finalErrPlus )));
        g_totErr->SetPointEYlow( i_bin-1,sqrt(TMath::Abs(finalErrMinus)));
    }
    
    return g_totErr;
}

//--------------- ~ ---------------
 
