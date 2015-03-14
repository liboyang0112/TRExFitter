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
  fFitName = "";
  fNSyst = 0;
  fNNorm = 0;
  fRegionName = "Region";
  fVariableTitle = "Variable";
  fSystSmoothed = false;
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
  fHist->SetLineWidth(1);
  fName = fSample->fName;
  fHistoName = histoName;
  fFileName = fileName;
  fNSyst = 0;
  fNNorm = 0;
  fRegionName = "Region";
  fVariableTitle = "Variable";
  fSystSmoothed = false;
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
void SampleHist::DrawSystPlot(string syst){
  // smooth here FIXME
  SmoothSyst(syst);
  //
  float yield_syst_up, yield_syst_down, yield_nominal;
  TCanvas *c = new TCanvas("c","c",800,600);
  TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
  h_nominal->SetLineColor(kBlack);
  h_nominal->SetLineWidth(2);
  h_nominal->SetLineStyle(2);
  h_nominal->SetFillStyle(0);
  TH1* h_1 = (TH1*)h_nominal->Clone("h_1");
  TH1* h_syst_up;
  TH1* h_syst_down;
  TH1* h_syst_up_orig;
  TH1* h_syst_down_orig;
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;
      h_nominal->SetMinimum(0);
      h_nominal->SetMaximum(h_nominal->GetMaximum()*2);
      h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
      h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
      h_syst_up_orig = (TH1*)fSyst[i_syst]->fHistUp_original->Clone();
      h_syst_down_orig = (TH1*)fSyst[i_syst]->fHistDown_original->Clone();
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
    yield_nominal = h_nominal->Integral();
    yield_syst_up = h_syst_up->Integral();
    yield_syst_down = h_syst_down->Integral();
    // draw Relative difference
    h_1->Scale(0);
    h_syst_up->Add(h_nominal,-1);
    h_syst_down->Add(h_nominal,-1);
    h_syst_up->Divide(h_nominal);
    h_syst_down->Divide(h_nominal);
    h_syst_up->Scale(100);
    h_syst_down->Scale(100);
    h_syst_up_orig->Add(h_nominal,-1);
    h_syst_down_orig->Add(h_nominal,-1);
    h_syst_up_orig->Divide(h_nominal);
    h_syst_down_orig->Divide(h_nominal);
    h_syst_up_orig->Scale(100);
    h_syst_down_orig->Scale(100);
    h_1->Draw("HIST");
    h_syst_down_orig->Draw("same HIST");
    h_syst_up_orig->Draw("same HIST");
    h_syst_down->Draw("same HIST");
    h_syst_up->Draw("same HIST");
    double ymax = 0;
    ymax = TMath::Max( ymax,TMath::Abs(h_syst_up->GetMaximum()));
    ymax = TMath::Max( ymax,TMath::Abs(h_syst_down->GetMaximum()));
    ymax = TMath::Max( ymax,TMath::Abs(h_syst_up->GetMinimum()));
    ymax = TMath::Max( ymax,TMath::Abs(h_syst_down->GetMinimum()));
    h_1->SetMinimum(-ymax*1.8);
    h_1->SetMaximum( ymax*1.8);
    h_1->GetYaxis()->SetTitle("Relative difference [%]");
    h_1->GetXaxis()->SetTitle(fVariableTitle.c_str());
    //
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->DrawLatex(0.2,0.89,Form("%s, %s",fSyst[i_syst]->fName.c_str(),fSample->fName.c_str()));
    tex->DrawLatex(0.2,0.84,fRegionName.c_str());
    TLegend *leg = new TLegend(0.6,0.8,0.9,0.94);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    float acc_up = (yield_syst_up-yield_nominal)/yield_nominal;
    string sign_up =  "+";
    if(acc_up<0) sign_up = "-";
    float acc_down = (yield_syst_down-yield_nominal)/yield_nominal;
    string sign_down =  "+";
    if(acc_down<0) sign_down = "-";
    leg->AddEntry(h_syst_up,  Form("+1#sigma (%s%.1f%)",sign_up.c_str(),  TMath::Abs(acc_up  *100)),"l");
    leg->AddEntry(h_syst_down,Form("-1#sigma (%s%.1f%)",sign_down.c_str(),TMath::Abs(acc_down*100)),"l");
    leg->Draw();
    //
//     c->SaveAs(Form("systPlots/%s_%s.png",fHist->GetName(),fSyst[i_syst]->fName.c_str()));
    cout << fFitName << endl;
    gSystem->mkdir(fFitName.c_str());
    gSystem->mkdir((fFitName+"/systPlots/").c_str());
//     string saveName = fFitName+"/systPlots/"+fHist->GetName()+"_"+fSyst[i_syst]->fName+".png";
    const char* saveName = Form("%s/systPlots/%s_%s.png",fFitName.c_str(),fHist->GetName(),fSyst[i_syst]->fName.c_str());
    c->SaveAs(saveName);
  }
  delete c;
}

void SampleHist::SmoothSyst(string syst,bool force){
  if(fSystSmoothed && !force) return;
  TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
  TH1* h_syst_up;
  TH1* h_syst_down;
  //
  for(int i_syst=0;i_syst<fNSyst;i_syst++){
    cout << fSyst[i_syst]->fName << endl;
    if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;
    h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
    h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
    //
    SmoothSystHistos(h_nominal,h_syst_up,h_syst_down); // see Root/Commmon.C
    //
    // save stuff
    fSyst[i_syst]->fHistUp_original = (TH1*)fSyst[i_syst]->fHistUp->Clone();
    fSyst[i_syst]->fHistUp = h_syst_up;
    fSyst[i_syst]->fHistDown_original = (TH1*)fSyst[i_syst]->fHistDown->Clone();
    fSyst[i_syst]->fHistDown = h_syst_down;
    if(fSyst[i_syst]->fIsShape){
      // update shpae hists as well
      fSyst[i_syst]->fHistShapeUp = (TH1*)h_syst_up->Clone(fSyst[i_syst]->fHistShapeUp->GetName());
      fSyst[i_syst]->fHistShapeDown = (TH1*)h_syst_down->Clone(fSyst[i_syst]->fHistShapeDown->GetName());
      fSyst[i_syst]->fHistShapeUp->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeUp->Integral());
      fSyst[i_syst]->fHistShapeDown->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeDown->Integral());
    }
  }
  fSystSmoothed = true;
}
