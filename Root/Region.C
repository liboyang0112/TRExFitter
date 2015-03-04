#include "TtHFitter/Region.h"

// -------------------------------------------------------------------------------------------------
// class Region

Region::Region(string name){
  fName = name;
  fLabel = name;
  fShortLabel = name;
  fNBkg = 0;
  fNSyst = 0;
  fHasData = false;
  fHasSig = false;
  fNSamples = 0;
  fUseStatErr = false;
  //
  string cName = "c_"+fName;
  fPlotPreFit = new TthPlot(cName);
  cName = "c_"+fName+"_postFit";
  fPlotPostFit = new TthPlot(cName);
}
Region::~Region(){}

SampleHist* Region::SetSampleHist(Sample *sample, string histoName, string fileName){
  fSampleHists[fNSamples] = new SampleHist( sample, histoName, fileName );
  if(sample->fType==SampleType::Data){
    fHasData = true;
    fData = fSampleHists[fNSamples];
  }
  else if(sample->fType==SampleType::Signal){
    fHasSig = true;
    fSig = fSampleHists[fNSamples];
  }
  else if(sample->fType==SampleType::Background){
    fBkg[fNBkg] = fSampleHists[fNSamples];
    fNBkg ++;
  }
  else{
    cout << "ERROR: SampleType not supported." << endl;
  }
  fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
  fSampleHists[fNSamples]->fRegionName = fName;
  fSampleHists[fNSamples]->fFitName = fFitName;
  fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
  fNSamples++;
  return fSampleHists[fNSamples-1];
}

SampleHist* Region::SetSampleHist(Sample *sample, TH1* hist ){
  fSampleHists[fNSamples] = new SampleHist( sample, hist );
  if(sample->fType==SampleType::Data){
    fHasData = true;
    fData = fSampleHists[fNSamples];
  }
  else if(sample->fType==SampleType::Signal){
    fHasSig = true;
    fSig = fSampleHists[fNSamples];
  }
  else if(sample->fType==SampleType::Background){
    fBkg[fNBkg] = fSampleHists[fNSamples];
    fNBkg ++;
  }
  else{
    cout << "ERROR: SampleType not supported." << endl;
  }
  fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
  fSampleHists[fNSamples]->fRegionName = fName;
  fSampleHists[fNSamples]->fFitName = fFitName;
  fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
  fNSamples++;
  return fSampleHists[fNSamples-1];
}

// SampleHist* Region::SetDataHist(Sample *sample, string histoName, string fileName){
//   fData = new SampleHist( sample, histoName, fileName, true, false );
//   fHasData = true;
//   fData->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
//   return fData;
// }
// SampleHist* Region::SetDataHist(Sample *sample, TH1* hist ){
//   fData = new SampleHist( sample, hist, true, false );
//   fHasData = true;
//   fData->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
//   return fData;
// }
// 
// SampleHist* Region::SetSigHist(Sample *sample, string histoName, string fileName){
//   fSig = new SampleHist( sample, histoName, fileName, false, true );
//   fHasSig = true;
//   fSig->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
//   return fSig;
// }
// SampleHist* Region::SetSigHist(Sample *sample, TH1* hist ){
//   fSig = new SampleHist( sample, hist, false, true );
//   fHasSig = true;
//   fSig->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
//   return fSig;
// }
// 
// SampleHist* Region::AddBkgHist(Sample *sample, string histoName, string fileName){
//   fBkg[fNBkg] = new SampleHist( sample, histoName, fileName, false, false );
//   fBkg[fNBkg]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
//   fNBkg ++;
//   return fBkg[fNBkg-1];
// }
// SampleHist* Region::AddBkgHist(Sample *sample, TH1* hist){
//   fBkg[fNBkg] = new SampleHist( sample, hist, false, false );
//   fBkg[fNBkg]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
//   fNBkg ++;
//   return fBkg[fNBkg-1];
// }

void Region::AddSample(Sample* sample){
  fSamples[fNSamples] = sample;
  fNSamples++;
}

void Region::AddSystematic(Systematic *syst){
  fSystematics[fNSyst] = syst;
  fNSyst++;
}

SampleHist* Region::GetSampleHist(string sampleName){
  for(int i_smp=0;i_smp<fNSamples;i_smp++){
//     if(fSampleHists[i_smp]->fSample->fName == sampleName) return fSampleHists[i_smp];
    if(fSampleHists[i_smp]->fName == sampleName) return fSampleHists[i_smp];
  }
  return 0x0;
}

// (internal)
void Region::BuildPreFitErrorHist(){
  if(TtHFitter::DEBUGLEVEL>0){
    cout << "-----------------------------------------------" << endl;
    cout << "->     Pre-Fit Plot for Region " << fName << endl;
    cout << "-----------------------------------------------" << endl;
  }
  else
    cout << "Building pre-fit plot for region " << fName << " ..." << endl;
  // build hTot
  for(int i=0;i<fNBkg;i++){
    if(i==0) fTot = (TH1*)fBkg[i]->fHist->Clone("hTot");
    else fTot->Add(fBkg[i]->fHist);
  }
  if(fHasSig){
    fTot->Add(fSig->fHist);
  }
  //
  // build error band
  float yieldNominal, yieldUp, yieldDown;
  float diffUp, diffDown;
  float diffPlus, diffMinus;
  float errPlus = 0.;
  float errMinus = 0.;
  float totYield = 0.;
  fErr = new TGraphAsymmErrors(fTot);
  //
  // collect all the systematics on all the samples
  vector<string> systNames;
  systNames.clear();
  map<string,bool> systIsThere;
  systIsThere.clear();
  string systName = "";
  if(TtHFitter::DEBUGLEVEL>0){
    cout << "Building syst list..." << endl;
  }
  for(int i=0;i<fNBkg;i++){
    for(int i_syst=0;i_syst<fBkg[i]->fNSyst;i_syst++){
      systName = fBkg[i]->fSyst[i_syst]->fName;
      if(!systIsThere[systName]){
        if(TtHFitter::DEBUGLEVEL>0){
          cout << " - " << systName << endl;
        }
        systNames.push_back(systName);
        systIsThere[systName] = true;
      }
    }
  }
  if(fHasSig){
    for(int i_syst=0;i_syst<fSig->fNSyst;i_syst++){
      systName = fSig->fSyst[i_syst]->fName;
      if(!systIsThere[systName]){
        if(TtHFitter::DEBUGLEVEL>0){
          cout << " - " << systName << endl;
        }
        systNames.push_back(systName);
        systIsThere[systName] = true;
      }
    }
  }
  //
  for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
    if(TtHFitter::DEBUGLEVEL>0){
      cout << "Bin " << i_bin << ":" << endl;
    }
    errPlus = 0.;
    errMinus = 0.;
    totYield = 0.;
    // syst = -1 is the stat unc.
    for(int i_syst=-1;i_syst<(int)systNames.size();i_syst++){
      if(i_syst<0){
        // add stat unc.
        if(TtHFitter::DEBUGLEVEL>0){
          cout << "  Adding stat uncertainty";
        }
        for(int i=0;i<fNBkg;i++){
          errPlus  += fBkg[i]->fHist->GetBinError(i_bin)*fBkg[i]->fHist->GetBinError(i_bin);
          errMinus += fBkg[i]->fHist->GetBinError(i_bin)*fBkg[i]->fHist->GetBinError(i_bin);
          totYield += fBkg[i]->fHist->GetBinContent(i_bin);
        }
        if(TtHFitter::DEBUGLEVEL>0){
          cout << "\t +" << 100*sqrt(errPlus)/totYield << "%";
          cout << "\t -" << 100*sqrt(errMinus)/totYield << "%";
          cout << endl;
        }
        continue;
      }
      if(TtHFitter::DEBUGLEVEL>0){
        cout << "  Adding syst " << systNames[i_syst];
      }
      diffUp = 0.;
      diffDown = 0.;
      TH1* hUp = 0x0;
      TH1* hDown = 0x0;
      for(int i=0;i<fNBkg;i++){
        yieldNominal = fBkg[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
        hUp = 0x0;
        hDown = 0x0;
        if(fBkg[i]->HasSyst(systNames[i_syst])){
          hUp   = fBkg[i]->GetSystematic(systNames[i_syst])->fHistUp;
          hDown = fBkg[i]->GetSystematic(systNames[i_syst])->fHistDown;
        }
        if(hUp!=0x0)    yieldUp     = hUp->GetBinContent(i_bin);
        else            yieldUp     = yieldNominal;
        if(hDown!=0x0)  yieldDown   = hDown->GetBinContent(i_bin);
        else            yieldDown   = yieldNominal;
        diffUp   += yieldUp   - yieldNominal;
        diffDown += yieldDown - yieldNominal;
      }
      if(fHasSig){
        yieldNominal = fSig->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
        hUp = 0x0;
        hDown = 0x0;
        if(fSig->HasSyst(systNames[i_syst])){
          hUp   = fSig->GetSystematic(systNames[i_syst])->fHistUp;
          hDown = fSig->GetSystematic(systNames[i_syst])->fHistDown;
        }
        if(hUp!=0x0)    yieldUp     = hUp->GetBinContent(i_bin);
        else            yieldUp     = yieldNominal;
        if(hDown!=0x0)  yieldDown   = hDown->GetBinContent(i_bin);
        else            yieldDown   = yieldNominal;
        diffUp   += yieldUp   - yieldNominal;
        diffDown += yieldDown - yieldNominal;
      }
      //
      if(diffUp>=0 && diffDown<=0)     errPlus += diffUp*diffUp;
      else if(diffDown>0 && diffUp<=0) errPlus += diffDown*diffDown;
      else if(diffUp>0 && diffDown>0)  errPlus += TMath::Max(diffUp,diffDown)*TMath::Max(diffUp,diffDown);
      if(diffUp<0 && diffDown>=0)      errMinus += diffUp*diffUp;
      else if(diffDown<0 && diffUp>=0) errMinus += diffDown*diffDown;
      else if(diffUp<0 && diffDown<0)  errMinus += TMath::Min(diffUp,diffDown)*TMath::Min(diffUp,diffDown);
      if(TtHFitter::DEBUGLEVEL>0){
        cout << "\t +" << 100*diffUp/totYield << "%";
        cout << "\t " << 100*diffDown/totYield << "%";
        cout << endl;
      }
    }
    errPlus = sqrt(errPlus);
    errMinus = sqrt(errMinus);

    fErr->SetPointEYhigh(i_bin-1,errPlus);
    fErr->SetPointEYlow(i_bin-1,errMinus);
  }
  // at this point fTot and fErr should be ready
}

TthPlot* Region::DrawPreFit(string opt){
//   string cName = "c_"+fName;
//   TthPlot *p = new TthPlot(cName);
  TthPlot *p = fPlotPreFit;
  p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
  p->SetChannel(fLabel);
  //
  if(fHasData && opt.find("blind")==string::npos) p->SetData(fData->fHist,fData->fSample->fTitle);
  if(fHasSig){
    p->AddSignal(fSig->fHist,fSig->fSample->fTitle);
    p->AddNormSignal(fSig->fHist,fSig->fSample->fTitle+"(norm)");
  }
  for(int i=0;i<fNBkg;i++)
    p->AddBackground(fBkg[i]->fHist,fBkg[i]->fSample->fTitle);
  //
  BuildPreFitErrorHist();
  //
  p->SetTotBkg((TH1*)fTot);
  p->SetTotBkgAsym(fErr);
  p->fATLASlabel = "Internal";
  p->Draw(opt);
  //
//   return p->GetCanvas();
//   fPlotPreFit = p;
  return p;
}


void Region::BuildPostFitErrorHist(FitResults *fitRes){
  if(TtHFitter::DEBUGLEVEL>0){
    cout << "-----------------------------------------------" << endl;
    cout << "->     Post-Fit Plot for Region " << fName << endl;
    cout << "-----------------------------------------------" << endl;
  }
  else
    cout << "Building post-fit plot for region " << fName << " ..." << endl;
  //
  float yieldNominal, yieldUp, yieldDown;
  float diffUp, diffDown;
  float diffPlus, diffMinus;
  float errPlus[MAXsyst];
  float errMinus[MAXsyst];
  float deltaN;
  fErr_postFit = new TGraphAsymmErrors(fTot_postFit);
  //
  // collect all the systematics on all the samples
  vector<string> systNames;
  systNames.clear();
  map<string,bool> systIsThere;
  systIsThere.clear();
  float systValue;
  float systErrUp;
  float systErrDown;
  TH1* hSyst;
  TH1* hNew;
  string systName = "";
//   cout << "Building syst list..." << endl;
  // backgrounds
  for(int i=0;i<fNBkg;i++){
    // norm factors
    for(int i_norm=0;i_norm<fBkg[i]->fNNorm;i_norm++){
      systName = fBkg[i]->fNormFactors[i_norm]->fName;
      if(!systIsThere[systName]){
//         cout << " " << systName << endl;
        systNames.push_back(systName);
        systIsThere[systName] = true;
      }
    }
    // syst
    for(int i_syst=0;i_syst<fBkg[i]->fNSyst;i_syst++){
      systName = fBkg[i]->fSyst[i_syst]->fName;
      if(!systIsThere[systName]){
        cout << " " << systName << endl;
        systNames.push_back(systName);
        systIsThere[systName] = true;
      }
    }
  }
  // signal
  if(fHasSig){
    // norm factors
    for(int i_norm=0;i_norm<fSig->fNNorm;i_norm++){
      systName = fSig->fNormFactors[i_norm]->fName;
      if(!systIsThere[systName]){
//         cout << " " << systName << endl;
        systNames.push_back(systName);
        systIsThere[systName] = true;
      }
    }
    // syst
    for(int i_syst=0;i_syst<fSig->fNSyst;i_syst++){
      systName = fSig->fSyst[i_syst]->fName;
      if(!systIsThere[systName]){
//         cout << " " << systName << endl;
        systNames.push_back(systName);
        systIsThere[systName] = true;
      }
    }
  }
  //
  for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
    if(TtHFitter::DEBUGLEVEL>0) cout << "Bin " << i_bin << ":" << endl;
    for(int i_syst=0;i_syst<(int)systNames.size();i_syst++){
      if(TtHFitter::DEBUGLEVEL>0) cout << "  Adding syst " << systNames[i_syst] << endl;
      systName = systNames[i_syst];
      systValue = fitRes->GetNuisParValue(systName);
      systErrUp = fitRes->GetNuisParErrUp(systName);
      systErrDown = fitRes->GetNuisParErrDown(systName);
      if(TtHFitter::DEBUGLEVEL>0) cout << "    alpha = " << systValue << " +" << systErrUp << " " << systErrDown << endl;
      diffUp = 0.;
      diffDown = 0.;
      TH1* hUp = 0x0;
      TH1* hDown = 0x0;
      for(int i=0;i<fNBkg;i++){
        if(TtHFitter::DEBUGLEVEL>0) cout << "    Sample " << fBkg[i]->fName << endl;
        yieldNominal = fBkg[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
//         yieldNominal = fBkg[i]->fHist_postFit->GetBinContent(i_bin);  // store nominal yield for this bin, but do it post fit!
        hUp = 0x0;
        hDown = 0x0;
        // norm
        if(fBkg[i]->HasNorm(systNames[i_syst])){
//           diffUp += yieldNominal*(systValue+systErrUp)-yieldNominal*systValue;
          diffUp += yieldNominal*systErrUp;
          diffDown += yieldNominal*systErrDown;
          if(TtHFitter::DEBUGLEVEL>0) cout << "\t +" << 100*diffUp/yieldNominal << "%\t " << 100*diffDown/yieldNominal << "%" << endl;
        }
        // syst
        if(fBkg[i]->HasSyst(systNames[i_syst])){
          hUp   = fBkg[i]->GetSystematic(systNames[i_syst])->fHistUp;
          hDown = fBkg[i]->GetSystematic(systNames[i_syst])->fHistDown;
          if(hUp!=0x0)    yieldUp     = hUp->GetBinContent(i_bin);
          else            yieldUp     = yieldNominal;
          if(hDown!=0x0)  yieldDown   = hDown->GetBinContent(i_bin);
          else            yieldDown   = yieldNominal;
          deltaN = GetDeltaN( systValue, yieldNominal,yieldUp,yieldDown);
          diffUp += yieldNominal*( GetDeltaN( systValue+systErrUp, yieldNominal,yieldUp,yieldDown) - deltaN );
          diffDown += yieldNominal*( GetDeltaN( systValue+systErrDown, yieldNominal,yieldUp,yieldDown) - deltaN );
          if(TtHFitter::DEBUGLEVEL>0) cout << "\t +" << 100*diffUp/yieldNominal << "%\t " << 100*diffDown/yieldNominal << "%" << endl;
        }
      }
      if(fHasSig){
        if(TtHFitter::DEBUGLEVEL>0) cout << "    Sample " << fSig->fName << endl;
        yieldNominal = fSig->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
//         yieldNominal = fSig->fHist_postFit->GetBinContent(i_bin);  // store nominal yield for this bin, but do it post-fit (wrong I think...)
        hUp = 0x0;
        hDown = 0x0;
        // norm
        if(fSig->HasNorm(systNames[i_syst])){
          diffUp += yieldNominal*systErrUp;
          diffDown += yieldNominal*systErrDown;
          if(TtHFitter::DEBUGLEVEL>0) cout << "\t +" << 100*diffUp/yieldNominal << "%\t " << 100*diffDown/yieldNominal << "%" << endl;
        }
        // syst
        if(fSig->HasSyst(systNames[i_syst])){
          hUp   = fSig->GetSystematic(systNames[i_syst])->fHistUp;
          hDown = fSig->GetSystematic(systNames[i_syst])->fHistDown;
          if(hUp!=0x0)    yieldUp     = hUp->GetBinContent(i_bin);
          else            yieldUp     = yieldNominal;
          if(hDown!=0x0)  yieldDown   = hDown->GetBinContent(i_bin);
          else            yieldDown   = yieldNominal;
          deltaN = GetDeltaN( systValue, yieldNominal,yieldUp,yieldDown);
          diffUp   += yieldNominal*( GetDeltaN( systValue+systErrUp, yieldNominal,yieldUp,yieldDown) - deltaN );
          diffDown += yieldNominal*( GetDeltaN( systValue+systErrDown, yieldNominal,yieldUp,yieldDown) - deltaN );
          if(TtHFitter::DEBUGLEVEL>0) cout << "\t +" << 100*diffUp/yieldNominal << "%\t " << 100*diffDown/yieldNominal << "%" << endl;
        }
      }
      // store errors up and down
      if(diffUp>=0 && diffDown<=0)     errPlus[i_syst] = diffUp;
      else if(diffDown>0 && diffUp<=0) errPlus[i_syst] = diffDown;
      else if(diffUp>0 && diffDown>0)  errPlus[i_syst] = TMath::Max(diffUp,diffDown);
      if(diffUp<0 && diffDown>=0)      errMinus[i_syst] = diffUp;
      else if(diffDown<0 && diffUp>=0) errMinus[i_syst] = diffDown;
      else if(diffUp<0 && diffDown<0)  errMinus[i_syst] = TMath::Min(diffUp,diffDown);
    }
    // Loop again on all the syst, two by two, to include the correlations
    float finalErrPlus = 0;
    float finalErrMinus = 0;
    float corr;
    for(int i_syst=0;i_syst<(int)systNames.size();i_syst++){
      for(int j_syst=0;j_syst<(int)systNames.size();j_syst++){
        corr = fitRes->fCorrMatrix->GetCorrelation(systNames[i_syst],systNames[j_syst]);
        finalErrPlus  += corr*errPlus[i_syst]*errPlus[j_syst];
        finalErrMinus += corr*errMinus[i_syst]*errMinus[j_syst];
      }
    }
    // add stat unc
    if(fUseStatErr){
      cout << "  Adding stat uncertainty" << endl;
      for(int i=0;i<fNBkg;i++){
        finalErrPlus  += pow(fBkg[i]->fHist->GetBinError(i_bin),2);
        finalErrMinus += pow(fBkg[i]->fHist->GetBinError(i_bin),2);
      }
      if(fHasSig){
        finalErrPlus  += pow(fSig->fHist->GetBinError(i_bin),2);
        finalErrMinus += pow(fSig->fHist->GetBinError(i_bin),2);
      }
    }
    //
    fErr_postFit->SetPointEYhigh(i_bin-1,sqrt(finalErrPlus));
    fErr_postFit->SetPointEYlow(i_bin-1,sqrt(finalErrMinus));
  }
  // at this point fTot and fErr _postFit should be ready
}

TthPlot* Region::DrawPostFit(FitResults *fitRes,string opt){
  TthPlot *p = fPlotPostFit;
  p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
  p->SetChannel(fLabel);
  //
  // 0) Create a new hist for each sample
  TH1* hBkgNew[MAXsamples];
  TH1* hSigNew;
  for(int i=0;i<fNBkg;i++){
    hBkgNew[i] = (TH1*)fBkg[i]->fHist->Clone();
  }
  if(fHasSig){
    hSigNew = (TH1*)fSig->fHist->Clone();
  }
  // 1) Scale all samples by norm factors (FIXME: before of after the syst? here's before...)
  string nfName;
  float nfValue;
  for(int i=0;i<fNBkg;i++){
    for(int i_norm=0;i_norm<fBkg[i]->fNNorm;i_norm++){
      nfName = fBkg[i]->fNormFactors[i_norm]->fName;
      nfValue = fitRes->GetNuisParValue(nfName);
      hBkgNew[i]->Scale(nfValue);
    }
  }
  if(fHasSig){
    for(int i_norm=0;i_norm<fSig->fNNorm;i_norm++){
      nfName = fSig->fNormFactors[i_norm]->fName;
      nfValue = fitRes->GetNuisParValue(nfName);
      hSigNew->Scale( nfValue );
    }
  }
  // 2) Scale all samples by the syst
  string systName;
  float systValue;
  float systErrUp;
  float systErrDown;
  TH1* hSyst;
  TH1* hNew;
  float binContent0;
  float binContentNew;
  float binContentUp;
  float binContentDown;
  for(int i=0;i<fNBkg;i++){
    hNew = (TH1*)hBkgNew[i]->Clone();
    for(int i_bin=1;i_bin<=hNew->GetNbinsX();i_bin++){
      binContent0 = hBkgNew[i]->GetBinContent(i_bin);
      binContentNew = binContent0;
      for(int i_syst=0;i_syst<fBkg[i]->fNSyst;i_syst++){
        systName = fBkg[i]->fSyst[i_syst]->fName;
        systValue = fitRes->GetNuisParValue(systName);
        binContentUp = fBkg[i]->fSyst[i_syst]->fHistUp->GetBinContent(i_bin);
        binContentDown = fBkg[i]->fSyst[i_syst]->fHistDown->GetBinContent(i_bin);
        binContentNew += (GetDeltaN(systValue,binContent0,binContentUp,binContentDown) - 1.)*binContent0;
//         binContentNew += (GetDeltaN(systValue,binContentNew,binContentUp,binContentDown) - 1.)*binContentNew;
//         binContentNew += (GetDeltaN(systValue,binContent0,binContentUp,binContentDown) - 1.)*binContentNew;
//         binContentNew += (GetDeltaN(systValue,binContentNew,binContentUp,binContentDown) - 1.)*binContent0;
//         if(systValue>0) binContentNew += systValue*(binContentUp-binContentNew); // linear - TEST
//         if(systValue<0) binContentNew -= systValue*(binContentDown-binContentNew); // linear - TEST
      }
      hNew->SetBinContent(i_bin,binContentNew);
    }
    hBkgNew[i] = (TH1*)hNew->Clone();
    fBkg[i]->fHist_postFit = hBkgNew[i];
    hNew->~TH1();
  }
  if(fHasSig){
    hNew = (TH1*)hSigNew->Clone();
    for(int i_bin=1;i_bin<=hNew->GetNbinsX();i_bin++){
      binContent0 = hSigNew->GetBinContent(i_bin);
      binContentNew = binContent0;
      for(int i_syst=0;i_syst<fSig->fNSyst;i_syst++){
        systName = fSig->fSyst[i_syst]->fName;
        systValue = fitRes->GetNuisParValue(systName);
        binContentUp = fSig->fSyst[i_syst]->fHistUp->GetBinContent(i_bin);
        binContentDown = fSig->fSyst[i_syst]->fHistDown->GetBinContent(i_bin);
        binContentNew += (GetDeltaN(systValue,binContent0,binContentUp,binContentDown) - 1.)*binContent0;
      }
      hNew->SetBinContent(i_bin,binContentNew);
    }
    hSigNew = (TH1*)hNew->Clone();
    fSig->fHist_postFit = hSigNew;
    hNew->~TH1();
  }
  
  // 3) Add the new Sig and Bkg to plot
  if(fHasData) p->SetData(fData->fHist,fData->fSample->fTitle);
  if(fHasSig)
    p->AddSignal(hSigNew,fSig->fSample->fTitle);
  for(int i=0;i<fNBkg;i++)
    p->AddBackground(hBkgNew[i],fBkg[i]->fSample->fTitle);
  //
  // 4) Build post-fit error band
  // build hTot
  for(int i=0;i<fNBkg;i++){
    if(i==0) fTot_postFit = (TH1*)hBkgNew[i]->Clone("hTot");
    else fTot_postFit->Add(hBkgNew[i]);
  }
  if(fHasSig){
    fTot_postFit->Add(hSigNew);
  }
  // Build error band
  BuildPostFitErrorHist(fitRes);
  //
  p->SetTotBkg(fTot_postFit);
  p->SetTotBkgAsym(fErr_postFit);
  p->Draw();
  //
  // print bin content and errors
  if(TtHFitter::DEBUGLEVEL>0){
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
//   return p->GetCanvas();
//   fPlotPostFit = p;
  return p;
}

void Region::AddSelection(string selection){
  if(fSelection=="") fSelection = selection;
  else fSelection += " && "+selection;
}

void Region::AddMCweight(string weight){
  if(fMCweight=="") fMCweight = weight;
  else fMCweight += " * "+weight;
}

void Region::SetVariable(string variable,int nbin,float xmin,float xmax){
  fVariable = variable;
  fNbins = nbin;
  fXmin = xmin;
  fXmax = xmax;
}

void Region::SetHistoName(string name){
  fHistoName = name;
}

void Region::SetVariableTitle(string name){
  fVariableTitle = name;
}

void Region::SetLabel(string label,string shortLabel){
  fLabel = label;
  if(shortLabel=="") fShortLabel = label;
  else fShortLabel = shortLabel;
}


void Region::SetAllSamples(bool readNtuples,string fileName){
//   TH1F* h;
//   TH1F* h_up;
//   TH1F* h_down;
//   TH1F* htmp;
//   string fullPath;
//   TFile *f;
//   if(!readNtuples)  f = new TFile(fileName.c_str());
//   for(int i_smp=0;i_smp<fNSamples;i_smp++){
//     SampleHist* sh;
//     h = 0x0;
//     if(readNtuples){
//       cout << "Building  " << Form("h_%s_%s",fSamples[i_smp]->fName.c_str(),fName.c_str()) << endl;
//       for(int i_path=0;i_path<(int)fSamples[i_smp]->fNtuplePaths.size();i_path++){
//         for(int i_name=0;i_name<(int)fSamples[i_smp]->fNtupleFiles.size();i_name++){
//           fullPath = fSamples[i_smp]->fNtuplePaths[i_path] + fSamples[i_smp]->fNtupleFiles[i_name] + "/" + fSamples[i_smp]->fNtupleName;
//           htmp = HistFromNtuple( fullPath, fVariable, fNbins, fXmin, fXmax, 
//                                  fSelection, fSamples[i_smp]->fIsData ? "1" : fSamples[i_smp]->fMCweight );
//           if(h==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fSamples[i_smp]->fName.c_str(),fName.c_str()));
//           else h->Add(htmp);
//           htmp->~TH1F();
//         }
//       }
//     }
//     else{
//       cout << "Getting " << Form("h_%s_%s",fSamples[i_smp]->fName.c_str(),fName.c_str()) << endl;
//       h = (TH1F*)f->Get( Form( "%s_%s", fName.c_str(), fSamples[i_smp]->fName.c_str() ) );
//     }
//     if(fSamples[i_smp]->fIsData)        sh = SetDataHist(fSamples[i_smp],(TH1*)h->Clone());
//     else if(fSamples[i_smp]->fIsSignal) sh = SetSigHist(fSamples[i_smp],(TH1*)h->Clone());
//     else                                sh = AddBkgHist(fSamples[i_smp],(TH1*)h->Clone());
//     h->~TH1F();
//     //
//     // add systematics
//     for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
//       // overall systs
//       if(fSamples[i_smp]->fSystematics[i_syst]->IsOverallOnly()){
//         cout << "Adding OverallSyst: " << fSamples[i_smp]->fSystematics[i_syst]->fName;
//         cout << " " << fSamples[i_smp]->fSystematics[i_syst]->fOverallUp;
//         cout << " " << fSamples[i_smp]->fSystematics[i_syst]->fOverallDown;
//         cout << endl;
//         sh->AddOverallSyst( fSamples[i_smp]->fSystematics[i_syst]->fName, 
//                             fSamples[i_smp]->fSystematics[i_syst]->fOverallUp, 
//                             fSamples[i_smp]->fSystematics[i_syst]->fOverallDown );
//       }
//       // histo syst
//       else{
//         cout << "Adding HistoSyst: " << fSamples[i_smp]->fSystematics[i_syst]->fName << endl;
//         if(readNtuples){
//           // up:
//           // - tree name
//           string treeNameUp = fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameUp;
//           if(treeNameUp=="") treeNameUp = fSamples[i_smp]->fNtupleName;
//           treeNameUp += fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameSufUp;
//           // - ntuple weight
//           string weightUp = fSamples[i_smp]->fSystematics[i_syst]->fWeightUp;
//           if(weightUp=="") weightUp = fSamples[i_smp]->fMCweight;
//           else             weightUp = fSamples[i_smp]->fMCweight + " * " + weightUp;
//           // - ntuple names
//           vector<string> ntupleNamesUp;
//           for(int i_name=0;i_name<(int)fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesUp.size();i_name++){
//             ntupleNamesUp.push_back(fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesUp[i_name]);
//           }
//           if(ntupleNamesUp.size()==0){
//             for(int i_name=0;i_name<(int)fSamples[i_smp]->fNtupleFiles.size();i_name++){
//               ntupleNamesUp.push_back(fSamples[i_smp]->fNtupleFiles[i_name]);
//             }
//           }
//           for(int i_name=0;i_name<(int)ntupleNamesUp.size();i_name++){
//             ntupleNamesUp[i_name] += fSamples[i_smp]->fSystematics[i_syst]->fNtupleFileSufUp;
//           }
//           // - ntuple paths
//           vector<string> ntuplePathsUp;
//           for(int i_path=0;i_path<(int)fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsUp.size();i_path++){
//             ntuplePathsUp.push_back(fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsUp[i_path]);
//           }
//           if(ntuplePathsUp.size()==0){
//             for(int i_path=0;i_path<(int)fSamples[i_smp]->fNtuplePaths.size();i_path++){
//               ntuplePathsUp.push_back(fSamples[i_smp]->fNtuplePaths[i_path]);
//             }
//           }
//           for(int i_path=0;i_path<(int)ntuplePathsUp.size();i_path++){
//             ntuplePathsUp[i_path] += fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathSufUp;
//           }
//           // down:
//           // - tree name
//           string treeNameDown = fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameDown;
//           if(treeNameDown=="") treeNameDown = fSamples[i_smp]->fNtupleName;
//           treeNameDown += fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameSufDown;
//           // - ntuple weight
//           // - ntuple weight
//           string weightDown = fSamples[i_smp]->fSystematics[i_syst]->fWeightDown;
//           if(weightDown=="") weightDown = fSamples[i_smp]->fMCweight;
//           else             weightDown = fSamples[i_smp]->fMCweight + " * " + weightDown;
// //           string weightDown = fSamples[i_smp]->fSystematics[i_syst]->fWeightDown;
// //           if(weightDown=="") weightDown = fSamples[i_smp]->fMCweight;
// //           weightDown += fSamples[i_smp]->fSystematics[i_syst]->fWeightSufDown;
//           // - ntuple names
//           vector<string> ntupleNamesDown;
//           for(int i_name=0;i_name<(int)fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesDown.size();i_name++){
//             ntupleNamesDown.push_back(fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesDown[i_name]);
//           }
//           if(ntupleNamesDown.size()==0){
//             for(int i_name=0;i_name<(int)fSamples[i_smp]->fNtupleFiles.size();i_name++){
//               ntupleNamesDown.push_back(fSamples[i_smp]->fNtupleFiles[i_name]);
//             }
//           }
//           for(int i_name=0;i_name<(int)ntupleNamesDown.size();i_name++){
//             ntupleNamesDown[i_name] += fSamples[i_smp]->fSystematics[i_syst]->fNtupleFileSufDown;
//           }
//           // - ntuple paths
//           vector<string> ntuplePathsDown;
//           for(int i_path=0;i_path<(int)fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsDown.size();i_path++){
//             ntuplePathsDown.push_back(fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsDown[i_path]);
//           }
//           if(ntuplePathsDown.size()==0){
//             for(int i_path=0;i_path<(int)fSamples[i_smp]->fNtuplePaths.size();i_path++){
//               ntuplePathsDown.push_back(fSamples[i_smp]->fNtuplePaths[i_path]);
//             }
//           }
//           for(int i_path=0;i_path<(int)ntuplePathsDown.size();i_path++){
//             ntuplePathsDown[i_path] += fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathSufDown;
//           }
//           //
//           cout << "\t" << ntuplePathsUp[0] << "\t" << ntupleNamesUp[0] << "\t" << treeNameUp << "\t" << weightUp << endl;
//           cout << "\t" << ntuplePathsDown[0] << "\t" << ntupleNamesDown[0] << "\t" << treeNameDown << "\t" << weightDown << endl;
//           //
//           // building the syst samples
//           h_up = 0x0;
//           h_down = 0x0;
//           // up
//           for(int i_path=0;i_path<(int)ntuplePathsUp.size();i_path++){
//             for(int i_name=0;i_name<(int)ntupleNamesUp.size();i_name++){
//               fullPath = ntuplePathsUp[i_path] + ntupleNamesUp[i_name] + "/" + treeNameUp;
//               htmp = HistFromNtuple( fullPath, fVariable, fNbins, fXmin, fXmax, 
//                                     fSelection, weightUp );
//               if(htmp==0x0){
//                 cout << "ERROR: cannot extract histo... " << endl;
//                 continue;
//               }
//               if(h_up==0) h_up = (TH1F*)htmp->Clone(Form("h_%s_%s_%s_Up",
//                                                 fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
//               else h_up->Add(htmp);
//               htmp->~TH1F();
//             }
//           }
//           // down
//           for(int i_path=0;i_path<(int)ntuplePathsDown.size();i_path++){
//             for(int i_name=0;i_name<(int)ntupleNamesDown.size();i_name++){
//               fullPath = ntuplePathsDown[i_path] + ntupleNamesDown[i_name] + "/" + treeNameDown;
//               htmp = HistFromNtuple( fullPath, fVariable, fNbins, fXmin, fXmax, 
//                                       fSelection, weightDown );
//               if(htmp==0x0){
//                 cout << "ERROR: cannot extract histo... " << endl;
//                 continue;
//               }
//               if(h_down==0) h_down = (TH1F*)htmp->Clone(Form("h_%s_%s_%s_Down",
//                                                   fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
//               else h_down->Add(htmp);
//               htmp->~TH1F();
//             }
//           }
//         }
//         else{
//           cout << "Getting " << Form("h_%s_%s_%s_Up",fName.c_str(), fSamples[i_smp]->fName.c_str(), fSamples[i_smp]->fSystematics[i_syst]->fName.c_str() ) << endl;
//           h_up = (TH1F*)f->Get( Form( "%s_%s_%s_Up", fName.c_str(), fSamples[i_smp]->fName.c_str(), fSamples[i_smp]->fSystematics[i_syst]->fName.c_str() ) );
//           cout << "Getting " << Form("h_%s_%s_%s_Down",fName.c_str(), fSamples[i_smp]->fName.c_str(), fSamples[i_smp]->fSystematics[i_syst]->fName.c_str() ) << endl;
//           h_down = (TH1F*)f->Get( Form( "%s_%s_%s_Down", fName.c_str(), fSamples[i_smp]->fName.c_str(), fSamples[i_smp]->fSystematics[i_syst]->fName.c_str() ) );
//         }
//         sh->AddHistoSyst( fSamples[i_smp]->fSystematics[i_syst]->fName, h_up, h_down );
//         h_up->~TH1F();
//         h_down->~TH1F();
//       }
//     }
//     //
//     // add norm factors
//     for(int i_norm=0;i_norm<fSamples[i_smp]->fNNorm;i_norm++){
//       cout << "Adding NormFactor: " << fSamples[i_smp]->fNormFactors[i_norm]->fName;
//       cout << " " << fSamples[i_smp]->fNormFactors[i_norm]->fNominal;
//       cout << "[" << fSamples[i_smp]->fNormFactors[i_norm]->fMin;
//       cout << "-" << fSamples[i_smp]->fNormFactors[i_norm]->fMax;
//       cout << "]" << endl;
//       sh->AddNormFactor( fSamples[i_smp]->fNormFactors[i_norm] );
//     }
//   }
}

void Region::Print(){
  cout << "    Region: " << fName << endl;
  for(int i_smp=0;i_smp<fNSamples;i_smp++){
    fSampleHists[i_smp]->Print();
  }
}



/////////////
float GetDeltaN(float alpha, float Iz, float Ip, float Imi){
  // protection against negative values
  if(Ip<0)  Ip  = 0.00001*Iz;
  if(Imi<0) Imi = 0.00001*Iz;
//   cout << "Running GetDeltaN." << endl;
//   cout << "  alpha = " << alpha << endl;
//   cout << "  I0 = " << Iz << endl;
//   cout << "  Ip = " << Ip << endl;
//   cout << "  Im = " << Imi << endl;
  float deltaN;
  if(alpha>0)      deltaN = Ip;
  else if(alpha<0) deltaN = Imi;
  else             return 0.;
  if(TMath::Abs(alpha)>=1){
    // exponential
    deltaN /= Iz; // divde h_tmp by the nominal
    deltaN = pow( deltaN, TMath::Abs(alpha) );  // d -> d^(|a|)
  }
  else{
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
    float a = alpha; //systValue[systName[i_sys]];
    deltaN = 1 + a1*a + a2*a*a + a3*a*a*a + a4*a*a*a*a + a5*a*a*a*a*a + a6*a*a*a*a*a*a;
  }
  if(deltaN!=deltaN) deltaN = 1;  // to avoid nan
//   cout << "  Returning = " << deltaN*Iz << endl;
  return deltaN;
}
/////////////
