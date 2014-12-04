#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TDirectory.h"

#include "Root/Common.C"
#include "Root/NuisParameter.C"
#include "Root/CorrelationMatrix.C"
#include "Root/FitResults.C"
#include "Root/Systematic.C"
#include "Root/SystematicHisto.C"
#include "Root/NormFactor.C"
#include "Root/Sample.C"
#include "Root/Region.C"
#include "Root/SampleHist.C"
#include "Root/TtHFit.C"

#include "Root/TthPlot.C"


void FitExample(){
    
//   // create dummy data
//   TDirectory *dir = gDirectory;
//   TFile *f = new TFile("inputs.root","RECREATE");
//   TH1F *h_data = new TH1F("h_data","h_data",10,0,1200);
//     h_data->FillRandom("pol0",100);
//     h_data->FillRandom("pol1",10);
// //     h_data->FillRandom("pol1",110);
//   TH1F *h_bkg = new TH1F("h_bkg","h_bkg",10,0,1200);
//     h_bkg->FillRandom("pol0",100000);
//     h_bkg->Scale(100./h_bkg->Integral());
//   TH1F *h_bkg_jesUp = new TH1F("h_bkg_jesUp","h_bkg_jesUp",10,0,1200);
//     h_bkg_jesUp->FillRandom("pol0",100000);
//     h_bkg_jesUp->FillRandom("pol1",10000);
//     h_bkg_jesUp->Scale(100./100000);
//   TH1F *h_bkg_jesDown = new TH1F("h_bkg_jesDown","h_bkg_jesDown",10,0,1200);
//     h_bkg_jesDown->FillRandom("pol1",10000);
//     h_bkg_jesDown->Scale(-1);
//     h_bkg_jesDown->FillRandom("pol0",100000);
//     h_bkg_jesDown->Scale(100./100000);
//   TH1F *h_sig = new TH1F("h_sig","h_sig",10,0,1200);
//     h_sig->FillRandom("pol1",1000);
//     h_sig->Scale(10./h_sig->Integral());
//   h_data->Write("",TObject::kOverwrite);
//   h_bkg->Write("",TObject::kOverwrite);
//   h_bkg_jesUp->Write("",TObject::kOverwrite);
//   h_bkg_jesDown->Write("",TObject::kOverwrite);
//   h_sig->Write("",TObject::kOverwrite);
//   dir->cd();
//   f->Close();
//   f->~TFile();
//   delete f;
//   return;
  
  // create the fit object
  TtHFit *myFit = new TtHFit();
    myFit->SetPOI("SigXsecOverSM");
    
  // create the samples
  Sample *data = myFit->NewSample("Data");
    data->SetTitle("Data 2012");
    data->SetIsData();
  Sample *bkg = myFit->NewSample("Background");
    bkg->SetFillColor(kWhite);
    bkg->SetLineColor(kBlack);
  Sample *sig = myFit->NewSample("Signal");
    sig->SetFillColor(kRed);
    sig->SetLineColor(kRed);
    sig->AddNormFactor("SigXsecOverSM",1,0,5);

  // create the systematics
  Systematic *BkgXsec = myFit->NewSystematic("BkgXsec");
  Systematic *JES = myFit->NewSystematic("JES");
    
  // create fit regions
  Region *SR_6j4b = myFit->NewRegion("SR_6j4b");
    SampleHist *hS_data = SR_6j4b->SetDataHist(data,"h_data","inputs.root");
    SampleHist *hS_bkg = SR_6j4b->AddBkgHist(bkg,"h_bkg","inputs.root");
    SampleHist *hS_sig = SR_6j4b->SetSigHist(sig,"h_sig","inputs.root");

  // add some systematics to the regions
  hS_bkg->AddOverallSyst("BkgXsec",0.10,-0.05);
  hS_bkg->AddHistoSyst("JES","h_bkg_jesUp","inputs.root","h_bkg_jesDown","inputs.root");
  
  TCanvas *c = SR_6j4b->DrawPreFit();
  
  c->SaveAs("Test.png");
  
  myFit->WriteHistos();  
  myFit->ReadAll();  
}
