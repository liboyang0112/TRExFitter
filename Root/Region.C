#include "TtHFitter/Region.h"

// -------------------------------------------------------------------------------------------------
// class Region

//__________________________________________________________________________________
//
Region::Region(string name){
    fName = name;
    fLabel = name;
    fShortLabel = name;
    fRegionType = CONTROL;
    fNBkg = 0;
    fNSyst = 0;
    fHasData = false;
    fHasSig = false;
    fNSamples = 0;
    fUseStatErr = false;
    fHistoBins = 0;
    fHistoNBinsRebin = -1;

    string cName = "c_"+fName;
    fPlotPreFit = new TthPlot(cName);
    fPlotPreFit->fShowYields = TtHFitter::SHOWYIELDS;
    cName = "c_"+fName+"_postFit";
    fPlotPostFit = new TthPlot(cName);
    fPlotPostFit->fShowYields = TtHFitter::SHOWYIELDS;
    
    fSampleHists.clear();
    fSamples.clear();
    fSystematics.clear();
    fNormFactors.clear();
    
    fIntCode_overall = 4;
    fIntCode_shape = 0;
    
    fFitName = "";
    fFitType = TtHFit::SPLUSB;
    fPOI = "";
    fFitLabel = "";
    
    fSelection = "1";
    fMCweight = "1";
    
    fLogScale = false;
    
    fBinWidth = 0;
    
    fLumiScale = 1.;
    
    fSig = 0x0;
//     fBkg = 0;
    
    fBlindingThreshold = -1;
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
        fSig = fSampleHists[fNSamples];
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg[fNBkg] = fSampleHists[fNSamples];
        fNBkg ++;
    }
    else{
        cout << "ERROR: SampleType not supported." << endl;
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
        fSig = fSampleHists[fNSamples];
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg[fNBkg] = fSampleHists[fNSamples];
        fNBkg ++;
    }
    else{
        cout << "ERROR: SampleType not supported." << endl;
    }
    fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists[fNSamples]->fRegionName = fName;
    fSampleHists[fNSamples]->fRegionLabel = fLabel;
    fSampleHists[fNSamples]->fFitName = fFitName;
    fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
    fNSamples++;
    return fSampleHists[fNSamples-1];
}

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
    fHistoNBinsRebin = N;
    fHistoBins = new double[N+1];
    for(unsigned int i=0; i<=N; ++i) fHistoBins[i] = bins[i];    
}

//__________________________________________________________________________________
//
void Region::Rebin(int N){
    fHistoNBinsRebin = N;
}

//__________________________________________________________________________________
//
void Region::SetRegionType( RegionType type ){
    fRegionType = type;
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
    if(TtHFitter::DEBUGLEVEL>0){
        cout << "-----------------------------------------------" << endl;
        cout << "->     Pre-Fit Plot for Region " << fName << endl;
        cout << "-----------------------------------------------" << endl;
    }
    else {
        cout << "Building pre-fit plot for region " << fName << " ..." << endl;
    }
    //
    float yieldNominal(0.), yieldUp(0.), yieldDown(0.);
    float diffUp(0.), diffDown(0.);
    //
    fSystNames.clear();
    std::map<string,bool> systIsThere;
    systIsThere.clear();
    TH1* hUp = 0x0;
    TH1* hDown = 0x0;
    string systName = "";
    SystematicHist *sh = 0x0;
    
    //
    // Collect all the systematics on all the samples
    //
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType == Sample::DATA) continue;
        
        //
        // Norm factors
        //
        for(int i_norm=0;i_norm<fSampleHists[i]->fNNorm;i_norm++){
            systName = fSampleHists[i]->fNormFactors[i_norm]->fName;
            // skip POI if B-only fit
            if(fFitType==TtHFit::BONLY && systName==fPOI) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }
        
        //
        // Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
            systName = fSampleHists[i]->fSyst[i_syst]->fName;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
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
        
        if(TtHFitter::DEBUGLEVEL>0) cout << "  Sample: " << fSampleHists[i]->fName << endl;
        
        // - loop on systematics
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            if(TtHFitter::DEBUGLEVEL>0) cout << "    Systematic: " << fSystNames[i_syst] << endl;
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
            
            // - loop on bins
            for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
                if(TtHFitter::DEBUGLEVEL>0) cout << "        Bin " << i_bin << ":  ";// << endl;
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
                if(TtHFitter::DEBUGLEVEL>0) cout << "\t +" << 100*diffUp/yieldNominal << "%\t " << 100*diffDown/yieldNominal << "%" << endl;
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
                // skip data
                if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
                // get SystematicHist
                sh = fSampleHists[i]->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);
                diffUp   += sh->fHistUp  ->GetBinContent(i_bin) - yieldNominal;
                diffDown += sh->fHistDown->GetBinContent(i_bin) - yieldNominal;
            }
            // add the proper bin content to the variation hists
            fTotUp[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }
    
    //
    // build the vectors of variations
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        h_up.  push_back( fTotUp[i_syst]   );
        h_down.push_back( fTotDown[i_syst] );
    }
    fErr = BuildTotError( fTot, h_up, h_down, fSystNames );
    fErr->SetName("g_totErr");
    
    // at this point fTot and fErr should be ready
}

//__________________________________________________________________________________
//
TthPlot* Region::DrawPreFit(string opt){
    
    TthPlot *p = fPlotPreFit;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    p->AddLabel(fFitLabel);
    p->AddLabel(fLabel);
    p->AddLabel("Pre-Fit");
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    if(fBlindingThreshold>=0) p->SetBinBlinding(true,fBlindingThreshold);
    
    //
    // build h_tot
    //
    fTot = 0x0;
    string title;
    if(fHasData && opt.find("blind")==string::npos) p->SetData(fData->fHist,fData->fSample->fTitle);
    if(fHasSig){
        title = fSig->fSample->fTitle;
        if(fSig->fSample->fGroup != "") title = fSig->fSample->fGroup;
        p->AddSignal(fSig->fHist,title);
        if(TtHFitter::SHOWNORMSIG) p->AddNormSignal(fSig->fHist,title+" (norm)");
        if(fTot==0x0) fTot = (TH1*)fSig->fHist->Clone("h_tot");
        else          fTot->Add(fSig->fHist);
    }
    for(int i=0;i<fNBkg;i++){
        title = fBkg[i]->fSample->fTitle;
        if(fBkg[i]->fSample->fGroup != "") title = fBkg[i]->fSample->fGroup;
        p->AddBackground(fBkg[i]->fHist,title);
        if(fTot==0x0) fTot = (TH1*)fBkg[i]->fHist->Clone("h_tot");
        else          fTot->Add(fBkg[i]->fHist);
    }
    
    //
    // set error to 0 if no MCstat
    //
    if(!fUseStatErr){
        for(int i_bin=1;i_bin<=fTot->GetNbinsX();i_bin++){
            fTot->SetBinError(i_bin,0);
        }
    }
    
    //
    // Computes the uncertainty bands arround the h_tot histogram
    //
    BuildPreFitErrorHist();
    
    //
    // Sets the last ingredients in the TthPlot object
    //
    p->SetTotBkg((TH1*)fTot);
    p->SetTotBkgAsym(fErr);
    p->fATLASlabel = "Internal";
    if(fLogScale) opt += " log";
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);
    p->Draw(opt);
    return p;
}

//__________________________________________________________________________________
//
void Region::BuildPostFitErrorHist(FitResults *fitRes){
    
    if(TtHFitter::DEBUGLEVEL>0){
        cout << "-----------------------------------------------" << endl;
        cout << "->     Post-Fit Plot for Region " << fName << endl;
        cout << "-----------------------------------------------" << endl;
    } else {
        cout << "Building post-fit plot for region " << fName << " ..." << endl;
    }
    
    float yieldNominal(0.), yieldUp(0.), yieldDown(0.);
    float yieldNominal_postFit = 0.;
    float diffUp(0.), diffDown(0.);
    float deltaN = 0.;
    float totYield = 0.;
    
    //
    // 0) Collect all the systematics on all the samples
    //
    fSystNames.clear();
    std::map<string,bool> systIsThere;
    systIsThere.clear();
    float systValue(0.);
    float systErrUp(0.);
    float systErrDown(0.);
    TH1* hUp = 0x0;
    TH1* hDown = 0x0;
    string systName = "";
    SystematicHist *sh = 0x0;
    
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType == Sample::DATA) continue;
        
        //
        // Norm factors
        //
        for(int i_norm=0;i_norm<fSampleHists[i]->fNNorm;i_norm++){
            systName = fSampleHists[i]->fNormFactors[i_norm]->fName;
            // skip POI if B-only fit FIXME
            if(fFitType==TtHFit::BONLY && systName==fPOI) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }
        
        //
        // Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
            systName = fSampleHists[i]->fSyst[i_syst]->fName;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
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
        if(TtHFitter::DEBUGLEVEL>0) cout << "  Sample: " << fSampleHists[i]->fName << endl;
        
        // - loop on systematics
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            
            if(TtHFitter::DEBUGLEVEL>0) cout << "    Systematic: " << fSystNames[i_syst] << endl;
            
            //
            // Get fit result
            //
            systName = fSystNames[i_syst];
            systValue   = fitRes->GetNuisParValue(systName);
            systErrUp   = fitRes->GetNuisParErrUp(systName);
            systErrDown = fitRes->GetNuisParErrDown(systName);
            
            if(TtHFitter::DEBUGLEVEL>0) cout << "      alpha = " << systValue << " +" << systErrUp << " " << systErrDown << endl;
            
            //
            // Get SystematicHist
            //
            sh = fSampleHists[i]->GetSystematic(systName);
            
            // hack: add a systematic hist if not there... FIXME
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
                
                if(TtHFitter::DEBUGLEVEL>0) cout << "        Bin " << i_bin << ":  ";// << endl;
                diffUp = 0.;
                diffDown = 0.;
                hUp = 0x0;
                hDown = 0x0;
                yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
                double yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);  // store nominal yield for this bin, but do it post fit
                
                // if it's a norm-factor ==> FIXME
                if(fSampleHists[i]->HasNorm(fSystNames[i_syst])){
                    diffUp   += yieldNominal*systErrUp;
                    diffDown += yieldNominal*systErrDown;
                }
                
                //
                // Systematics treatment
                //
                else if(fSampleHists[i]->HasSyst(fSystNames[i_syst])){
                    hUp   = sh->fHistUp;
                    hDown = sh->fHistDown;
                    if(hUp!=0x0)    yieldUp     = hUp ->GetBinContent(i_bin);
                    else            yieldUp     = yieldNominal;
                    if(hDown!=0x0)  yieldDown   = hDown -> GetBinContent(i_bin);
                    else            yieldDown   = yieldNominal;
                    
                    //
                    // Compute updated relative syst. variations (linear scaling)
                    //
                    double scaleUp = (yieldUp/yieldNominal-1.)*(systErrUp-systValue);
                    double scaleDown = (yieldDown/yieldNominal-1.)*(systValue-systErrDown);

                    diffUp += scaleUp*yieldNominal_postFit;
                    diffDown += scaleDown*yieldNominal_postFit;
                }
                
                if(TtHFitter::DEBUGLEVEL>0) cout << "\t +" << 100*diffUp/yieldNominal << "%\t " << 100*diffDown/yieldNominal << "%" << endl;

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
        fTotDown_postFit[i_syst] = (TH1*)fTot_postFit->Clone(Form("h_tot_%s_Down_postFit",systName.c_str()));
        
        // - loop on bins
        for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
            diffUp = 0.;
            diffDown = 0.;
            // - loop on samples
            for(int i=0;i<fNSamples;i++){
                // skip data
                if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
                // get SystematicHist
                sh = fSampleHists[i]->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);
                diffUp   += sh->fHistUp_postFit  ->GetBinContent(i_bin) - yieldNominal_postFit;
                diffDown += sh->fHistDown_postFit->GetBinContent(i_bin) - yieldNominal_postFit;
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
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        h_up.  push_back( fTotUp_postFit[i_syst]   );
        h_down.push_back( fTotDown_postFit[i_syst] );
    }
    fErr_postFit = BuildTotError( fTot_postFit, h_up, h_down, fSystNames, fitRes->fCorrMatrix );
    fErr_postFit->SetName("g_totErr_postFit");
    // at this point fTot and fErr _postFit should be ready
}

//__________________________________________________________________________________
//
TthPlot* Region::DrawPostFit(FitResults *fitRes,string opt){
    
    TthPlot *p = fPlotPostFit;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    p->AddLabel(fFitLabel);
    p->AddLabel(fLabel);
    p->AddLabel("Post-Fit");
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    if(fBlindingThreshold>=0) p->SetBinBlinding(true,fBlindingThreshold);
    
    //
    // 0) Create a new hist for each sample
    //
    TH1* hSmpNew[MAXsamples];
    for(int i=0;i<fNSamples;i++){
        hSmpNew[i] = (TH1*)fSampleHists[i]->fHist->Clone();
    }
    
    //
    // 1) Propagates the post-fit NP values to the central value (pulls)
    //
    string systName;
    float systValue;
    float systErrUp;
    float systErrDown;
    TH1* hSyst;
    TH1* hNew;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        hNew = (TH1*)hSmpNew[i]->Clone();
        for(int i_bin=1;i_bin<=hNew->GetNbinsX();i_bin++){
            double binContent0 = hSmpNew[i]->GetBinContent(i_bin);
            double multShape = 0.;
            double multNorm = 1.;
            for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
                systName = fSampleHists[i]->fSyst[i_syst]->fName;
                systValue = fitRes->GetNuisParValue(systName);
                
                //
                // Normalisation component: use the exponential interpolation and the multiplicative combination
                //
                if(fSampleHists[i]->fSyst[i_syst]->fIsOverall){
                    float binContentUp   = (fSampleHists[i]->fSyst[i_syst]->fNormUp+1) * binContent0;
                    float binContentDown = (fSampleHists[i]->fSyst[i_syst]->fNormDown+1) * binContent0;
                    multNorm *= (GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall));
                }
                
                //
                // Shape component: use the linear interpolation and the additive combination
                //
                if(fSampleHists[i]->fSyst[i_syst]->fIsShape){
                    float binContentUp   = fSampleHists[i]->fSyst[i_syst]->fHistShapeUp->GetBinContent(i_bin);
                    float binContentDown = fSampleHists[i]->fSyst[i_syst]->fHistShapeDown->GetBinContent(i_bin);
                    multShape += (GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape) -1 );
                }
            }
            
            //
            // Final computation
            //
            double binContentNew = binContent0*multNorm*(multShape+1.);
            
            //
            // Setting to new values
            //
            hNew->SetBinContent(i_bin,binContentNew);
        }
        hSmpNew[i] = (TH1*)hNew->Clone();
        fSampleHists[i]->fHist_postFit = hSmpNew[i];
        hNew->~TH1();
    }
    
    //
    // 2) Scale all samples by norm factors
    //    Done after the propagation of the NP (avoids nans due to "0" value of some NormFactors)
    //    Seems consistent with application in Roostats
    //
    string nfName;
    float nfValue;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        for(int i_norm=0;i_norm<fSampleHists[i]->fNNorm;i_norm++){
            nfName = fSampleHists[i]->fNormFactors[i_norm]->fName;
            nfValue = fitRes->GetNuisParValue(nfName);
            hSmpNew[i]->Scale(nfValue);
        }
    }
    
    //
    // 3) Add the new Sig and Bkg to plot
    //
    string title;
    TH1* hBkgNew[MAXsamples];
    TH1* hSigNew;
    for(int i=0, i_bkg=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::BACKGROUND){
            hBkgNew[i_bkg] = hSmpNew[i];
            i_bkg++;
        }
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL)
            hSigNew = hSmpNew[i];
    }
    if(fHasData) p->SetData(fData->fHist,fData->fSample->fTitle);
    if(fHasSig){
        title = fSig->fSample->fTitle;
        if(fSig->fSample->fGroup != "") title = fSig->fSample->fGroup;
        p->AddSignal(hSigNew,title);
        if(TtHFitter::SHOWNORMSIG){
            p->AddNormSignal(hSigNew,title+" (norm)");
        }
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
        if(j==0) fTot_postFit = (TH1*)hSmpNew[i]->Clone("h_tot_postFit");
        else fTot_postFit->Add(hSmpNew[i]);
        j++;
    }
    
    //
    // MCstat set to 0 if disabled
    //
    if(!fUseStatErr){
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
            fTot_postFit->SetBinError(i_bin,0);
        }
    }
    
    //
    // Build error band
    //
    BuildPostFitErrorHist(fitRes);
    
    //
    // 5) Finishes configuration of TthPlot objects
    //
    p->SetTotBkg(fTot_postFit);
    p->SetTotBkgAsym(fErr_postFit);
    p->fATLASlabel = "Internal"; //FIXME
    if(fLogScale) opt += " log";
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);
    p->Draw(opt);
    
    //
    // Print bin content and errors
    //
    if(TtHFitter::DEBUGLEVEL>0){
        cout << "--------------------" << endl;
        cout << "Final bin contents" << endl;
        cout << "--------------------" << endl;
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
            cout << i_bin << ":\t";
            cout << fTot_postFit->GetBinContent(i_bin);
            cout << " +";
            cout << fErr_postFit->GetErrorYhigh(i_bin-1);
            cout << " -";
            cout << fErr_postFit->GetErrorYlow(i_bin-1);
            cout << endl;
        }
    }
    
    //
    // Save in a root file...
    //
    cout << "Writing file " << fFitName+"/Histograms/"+fName+"_postFit.root" << endl;
    TFile *f = new TFile((fFitName+"/Histograms/"+fName+"_postFit.root").c_str(),"RECREATE");
    fErr_postFit->Write("",TObject::kOverwrite);
    fTot_postFit->Write("",TObject::kOverwrite);
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
                if(fSampleHists[i]->fSyst[i_syst]) fSampleHists[i]->fSyst[i_syst]->fHistUp_postFit  ->Write(
                    Form("h_%s_%s_Up_postFit",fSampleHists[i]->fName.c_str(),fSampleHists[i]->fSyst[i_syst]->fName.c_str()), TObject::kOverwrite);
                if(fSampleHists[i]->fSyst[i_syst]) fSampleHists[i]->fSyst[i_syst]->fHistDown_postFit->Write(
                    Form("h_%s_%s_Down_postFit",fSampleHists[i]->fName.c_str(),fSampleHists[i]->fSyst[i_syst]->fName.c_str()),TObject::kOverwrite);
            }
        }
    }
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
void Region::SetVariable(string variable,int nbin,float xmin,float xmax){
    fVariable = variable;
    fNbins = nbin;
    fXmin = xmin;
    fXmax = xmax;
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
    cout << "    Region: " << fName << endl;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        fSampleHists[i_smp]->Print();
    }
}


// --------------- Functions --------------- //

float GetDeltaN(float alpha, float Iz, float Ip, float Imi, int intCode){
    // protection against negative values
    //   if(Ip<0)  Ip  = 0.00001*Iz;
    //   if(Imi<0) Imi = 0.00001*Iz;
    
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
    float yieldNominal(0.);
    float errUp_i(0.),   errUp_j(0.);
    float errDown_i(0.), errDown_j(0.);
    //
    // - loop on bins
    for(int i_bin=1;i_bin<h_nominal->GetNbinsX()+1;i_bin++){
        finalErrPlus = 0;
        finalErrMinus = 0;
        corr = 0;
        yieldNominal = h_nominal->GetBinContent(i_bin);
        double err =h_nominal->GetBinError(i_bin);
        
        // - loop on the syst, two by two, to include the correlations
        for(unsigned int i_syst=0;i_syst<fSystNames.size();i_syst++){
            for(unsigned int j_syst=0;j_syst<fSystNames.size();j_syst++){
                if(matrix!=0x0){
                    corr = matrix->GetCorrelation(fSystNames[i_syst],fSystNames[j_syst]);
                } else {
                    if(i_syst==j_syst) corr = 1.;
                    else               corr = 0.;
                }
                errUp_i   = h_up[i_syst]  ->GetBinContent(i_bin) - yieldNominal;
                errDown_i = h_down[i_syst]->GetBinContent(i_bin) - yieldNominal;
                errUp_j   = h_up[j_syst]  ->GetBinContent(i_bin) - yieldNominal;
                errDown_j = h_down[j_syst]->GetBinContent(i_bin) - yieldNominal;
                
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

                // Michele's note: at this point errPlus should contain the positive variation, errMinus the negative
//                 finalErrPlus  += corr*errPlus[i_syst]*errPlus[j_syst];
//                 finalErrMinus += corr*errMinus[i_syst]*errMinus[j_syst];
//                 // TEST:
 //                finalErrPlus  += corr * errUp_i   * errUp_j;
 //                finalErrMinus += corr * errDown_i * errDown_j;
                
//
//                
//                // TEST 1:
//                
//                if(corr>=0) finalErrPlus  +=  corr * errUp_i   * errUp_j;
//                else        finalErrPlus  += -corr * errUp_i   * errDown_j;
//                if(corr>=0) finalErrMinus +=  corr * errDown_i * errDown_j;
//                else        finalErrMinus += -corr * errDown_i * errUp_j;
//
//                 if(i_bin==1) cout << sqrt(TMath::Abs(finalErrPlus)) << endl;
                // TEST 2:
//                 if(errUp[i_syst]>=0 && errUp[j_syst]>=0) finalErrPlus  += corr*errUp[i_syst]*errUp[j_syst];
//                 else if(errUp[i_syst]>=0 && errDown[j_syst]>=0) finalErrPlus  += -corr*errUp[i_syst]*errDown[j_syst];
//                 else if(errDown[i_syst]>=0 && errUp[j_syst]>=0) finalErrPlus  += -corr*errDown[i_syst]*errUp[j_syst];
//                 else if(errDown[i_syst]>=0 && errDown[j_syst]>=0) finalErrPlus  += corr*errDown[i_syst]*errDown[j_syst];
//                 if(errDown[i_syst]<=0 && errDown[j_syst]<=0) finalErrMinus += corr*errDown[i_syst]*errDown[j_syst];
//                 else if(errDown[i_syst]<=0 && errUp[j_syst]<=0) finalErrMinus += -corr*errDown[i_syst]*errUp[j_syst];
//                 else if(errUp[i_syst]<=0 && errDown[j_syst]<=0) finalErrMinus += -corr*errUp[i_syst]*errDown[j_syst];
//                 else if(errUp[i_syst]<=0 && errUp[j_syst]<=0) finalErrMinus += corr*errUp[i_syst]*errUp[j_syst];
            }
        }
        // add stat uncertainty, which should have been stored as orignal bin errors in the h_nominal (if fUseStatErr is true)
        finalErrPlus  += pow( h_nominal->GetBinError(i_bin), 2 );
        finalErrMinus += pow( h_nominal->GetBinError(i_bin), 2 );
        //
        g_totErr->SetPointEYhigh(i_bin-1,sqrt(TMath::Abs(finalErrPlus )));
        g_totErr->SetPointEYlow( i_bin-1,sqrt(TMath::Abs(finalErrMinus)));
    }
    //
    return g_totErr;
}

//--------------- ~ ---------------
