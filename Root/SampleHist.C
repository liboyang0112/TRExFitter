#include "TtHFitter/SampleHist.h"

// -------------------------------------------------------------------------------------------------
// SampleHist

SampleHist::SampleHist(Sample *sample,TH1 *hist){
  fSample = sample;
  fHist = (TH1*)hist->Clone(Form("h_%s",fName.c_str()));
  fHist->SetFillColor(fSample->fFillColor);
  fHist->SetLineColor(fSample->fLineColor);
  fName = fSample->fName;
  fHistoName = "";
  fFileName = "";
  fNSyst = 0;
  fNNorm = 0;
  // add overall systematics and normFactors from sample
  for(int i_syst=0;i_syst<sample->fNSyst;i_syst++){
    if(sample->fSystematics[i_syst]->fType == SystType::Overall)
      AddOverallSyst(sample->fSystematics[i_syst]->fName,sample->fSystematics[i_syst]->fOverallUp,sample->fSystematics[i_syst]->fOverallDown);
  }
  for(int i_norm=0;i_norm<sample->fNNorm;i_norm++){
    AddNormFactor(sample->fNormFactors[i_norm]);
  }
}
SampleHist::SampleHist(Sample *sample, string histoName, string fileName){
  fSample = sample;
  fHist = HistFromFile(fileName,histoName);
  fHist->SetFillColor(fSample->fFillColor);
  fHist->SetLineColor(fSample->fLineColor);
  fName = fSample->fName;
  fHistoName = histoName;
  fFileName = fileName;
  fNSyst = 0;
  fNNorm = 0;
  // add overall systematics and normFactors from sample
  for(int i_syst=0;i_syst<sample->fNSyst;i_syst++){
    if(sample->fSystematics[i_syst]->fType == SystType::Overall)
      AddOverallSyst(sample->fSystematics[i_syst]->fName,sample->fSystematics[i_syst]->fOverallUp,sample->fSystematics[i_syst]->fOverallDown);
  }
  for(int i_norm=0;i_norm<sample->fNNorm;i_norm++){
    AddNormFactor(sample->fNormFactors[i_norm]);
  }
}
SampleHist::~SampleHist(){}

SystematicHist* SampleHist::AddOverallSyst(string name,float up,float down){
  SystematicHist *sh;
  // try if it's already there...
  sh = GetSystematic(name);
  // ... and if not create a new one
  if(sh==0x0){
    sh = new SystematicHist(name);
    fSyst[fNSyst] = sh;
    fNSyst ++;
  }
  //
  sh->fHistUp   = (TH1*)fHist->Clone(Form("%s_%s_Up",fHist->GetName(),name.c_str()));
  sh->fHistDown = (TH1*)fHist->Clone(Form("%s_%s_Down",fHist->GetName(),name.c_str()));
  sh->fHistUp->Scale(1.+up);
  sh->fHistDown->Scale(1.+down);
  sh->fIsOverall = true;
  sh->fIsShape   = false;
  sh->fNormUp   = up;
  sh->fNormDown = down;
  return sh;
}

SystematicHist* SampleHist::AddHistoSyst(string name,TH1* h_up,TH1* h_down){
  SystematicHist *sh;
  // try if it's already there...
  sh = GetSystematic(name);
  // ... and if not create a new one
  if(sh==0x0){
    sh = new SystematicHist(name);
    fSyst[fNSyst] = sh;
    fNSyst ++;
  }
  //
  sh->fHistUp   = (TH1*)h_up->Clone(Form("%s_%s_Up",fHist->GetName(),name.c_str()));
  sh->fHistDown = (TH1*)h_down->Clone(Form("%s_%s_Down",fHist->GetName(),name.c_str()));
  sh->fHistShapeUp = (TH1*)h_up->Clone(Form("%s_%s_Shape_Up",fHist->GetName(),name.c_str()));
  sh->fHistShapeDown = (TH1*)h_down->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),name.c_str()));
  sh->fHistShapeUp->Scale(fHist->Integral() / sh->fHistShapeUp->Integral());
  sh->fHistShapeDown->Scale(fHist->Integral() / sh->fHistShapeDown->Integral());
  sh->fIsOverall = true;
  sh->fIsShape   = true;
  sh->fNormUp   = ( sh->fHistUp->Integral() - fHist->Integral() ) / fHist->Integral();
  sh->fNormDown = ( sh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
  if(sh->fNormUp == 0 && sh->fNormDown == 0) sh->fIsOverall = false;
  return sh;
}

SystematicHist* SampleHist::AddHistoSyst(string name,string histoName_up, string fileName_up,string histoName_down, string fileName_down){
  SystematicHist *sh;
  // try if it's already there...
  sh = GetSystematic(name);
  // ... and if not create a new one
  if(sh==0x0){
    sh = new SystematicHist(name);
    fSyst[fNSyst] = sh;
    fNSyst ++;
  }
  //
  sh->fFileNameUp = fileName_up;
  sh->fFileNameDown = fileName_down;
  sh->fHistoNameUp = histoName_up;
  sh->fHistoNameDown = histoName_down;
  sh->fHistUp = HistFromFile(sh->fFileNameUp,sh->fHistoNameUp);
  sh->fHistDown = HistFromFile(sh->fFileNameDown,sh->fHistoNameDown);
  sh->fHistShapeUp = (TH1*)sh->fHistUp->Clone(Form("%s_%s_Shape_Up",fHist->GetName(),name.c_str()));
  sh->fHistShapeDown = (TH1*)sh->fHistDown->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),name.c_str()));
  sh->fHistShapeUp->Scale(fHist->Integral() / sh->fHistShapeUp->Integral());
  sh->fHistShapeDown->Scale(fHist->Integral() / sh->fHistShapeDown->Integral());
  sh->fIsOverall = true;
  sh->fIsShape   = true;
  sh->fNormUp   = ( sh->fHistUp->Integral() - fHist->Integral() ) / fHist->Integral();
  sh->fNormDown = ( sh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
  if(sh->fNormUp == 0 && sh->fNormDown == 0) sh->fIsOverall = false;
  return sh;
}

NormFactor* SampleHist::AddNormFactor(NormFactor *normFactor){
  NormFactor *norm = GetNormFactor(normFactor->fName);
  if(norm==0x0){
    fNormFactors[fNNorm] = normFactor;
    fNNorm ++;
  }
  else{
    norm = normFactor;
  }
}

NormFactor* SampleHist::AddNormFactor(string name,float nominal, float min, float max){
  NormFactor *norm = GetNormFactor(name);
  if(norm==0x0){
    fNormFactors[fNNorm] = new NormFactor(name,nominal,min,max);
    fNNorm ++;
  }
  else{
    norm = new NormFactor(name,nominal,min,max);;
  }
}

SystematicHist* SampleHist::GetSystematic(string systName){
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    if(systName == fSyst[i_syst]->fName) return fSyst[i_syst];
  }
  return 0x0;
}

NormFactor* SampleHist::GetNormFactor(string name){
  for(int i_syst=0;i_syst<fNNorm;i_syst++){
    if(name == fNormFactors[i_syst]->fName) return fNormFactors[i_syst];
  }
  return 0x0;
}

bool SampleHist::HasSyst(string name){
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    if(fSyst[i_syst]->fName == name) return true;
  }
  return false;
}

bool SampleHist::HasNorm(string name){
  for(int i_norm=0;i_norm<fNNorm;i_norm++){
    if(fNormFactors[i_norm]->fName == name) return true;
  }
  return false;
}

void SampleHist::WriteToFile(){
  WriteHistToFile(fHist,fFileName);
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    fSyst[i_syst]->WriteToFile();
  }
}

void SampleHist::ReadFromFile(){
  fHist = HistFromFile(fFileName,fHistoName);
}

void SampleHist::FixEmptyBins(){
  for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
    if(fHist->GetBinContent(i_bin)<=0) fHist->SetBinContent(i_bin,1e-3);
  }
}

void SampleHist::Print(){
  cout << "      Sample: " << fName << "\t" << fHist->GetName() << endl;
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    fSyst[i_syst]->Print();
  }
  for(int i_norm=0;i_norm<fNNorm;i_norm++){
    fNormFactors[i_norm]->Print();
  }
}

void SampleHist::Rebin(int ngroup, const Double_t* xbins){
  fHist->Rebin(ngroup,"",xbins);
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    if(fSyst[i_syst]->fHistUp!=0x0) fSyst[i_syst]->fHistUp->Rebin(ngroup,"",xbins);
    if(fSyst[i_syst]->fHistDown!=0x0) fSyst[i_syst]->fHistDown->Rebin(ngroup,"",xbins);
    if(fSyst[i_syst]->fHistShapeUp!=0x0) fSyst[i_syst]->fHistShapeUp->Rebin(ngroup,"",xbins);
    if(fSyst[i_syst]->fHistShapeDown!=0x0) fSyst[i_syst]->fHistShapeDown->Rebin(ngroup,"",xbins);
  }
}

void SampleHist::Smooth(int ntimes){
  fHist->Smooth(ntimes);
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    if(fSyst[i_syst]->fHistUp!=0x0) fSyst[i_syst]->fHistUp->Smooth(ntimes);
    if(fSyst[i_syst]->fHistDown!=0x0) fSyst[i_syst]->fHistDown->Smooth(ntimes);
    if(fSyst[i_syst]->fHistShapeUp!=0x0) fSyst[i_syst]->fHistShapeUp->Smooth(ntimes);
    if(fSyst[i_syst]->fHistShapeDown!=0x0) fSyst[i_syst]->fHistShapeDown->Smooth(ntimes);
  }
}

// this draws the control plots (for each systematic) with the syst variations for this region & sample
void SampleHist::DrawSystPlot(){
  TCanvas *c = new TCanvas("c","c",600,600);
  c->Divide(1,2);
  TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
    h_nominal->SetLineColor(kBlack);
    h_nominal->SetLineWidth(2);
    h_nominal->SetLineStyle(2);
    h_nominal->SetFillStyle(0);
  TH1* h_1 = (TH1*)h_nominal->Clone("h_1");
    h_1->Divide(h_1);
  TH1* h_syst_up;
  TH1* h_syst_down;
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    c->cd(1);
    gPad->SetBottomMargin(0);
    h_nominal->Draw("HIST");
    h_nominal->SetMinimum(0);
    h_nominal->SetMaximum(h_nominal->GetMaximum()*2);
    h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
    h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
      h_syst_up->SetLineColor(kRed);
      h_syst_up->SetLineWidth(2);
      h_syst_up->SetLineStyle(1);
      h_syst_up->SetFillStyle(0);
      h_syst_down->SetLineColor(kBlue);
      h_syst_down->SetLineWidth(2);
      h_syst_down->SetLineStyle(1);
      h_syst_down->SetFillStyle(0);
    h_syst_down->DrawClone("same HIST");
    h_syst_up->DrawClone("same HIST");
    c->cd(2);
//     gPad = new TPad("c_2", "c_2", .1, .1, .9, .3);
//     gPad->Draw();
    gPad->SetTopMargin(0);
    h_syst_up->Divide(h_nominal);
    h_syst_down->Divide(h_nominal);
    h_1->Draw("HIST");
    h_syst_down->Draw("same HIST");
    h_syst_up->Draw("same HIST");
    h_1->SetMinimum(0);
    h_1->SetMaximum(2);
    c->SaveAs(Form("systPlots/%s_%s.png",fHist->GetName(),fSyst[i_syst]->fName.c_str()));
  }
}