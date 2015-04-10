#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "Root/Common.C"
#include "Root/NuisParameter.C"
#include "Root/CorrelationMatrix.C"
#include "Root/FitResults.C"
#include "Root/Systematic.C"
#include "Root/SystematicHist.C"
#include "Root/NormFactor.C"
#include "Root/Sample.C"
#include "Root/Region.C"
#include "Root/SampleHist.C"
#include "Root/TtHFit.C"

#include "Root/TthPlot.C"

#include "AtlasStyle.C"
#include "AtlasLabels.C"
#include "AtlasUtils.C"


// Simplest example:
//  - a single region is used
//  - input (dummy) histograms are created on the flight and set manually

void FitExample(){

  TtHFitter::SetDebugLevel(1);
  
  TH1F *h_data = new TH1F("h_data","h_data",10,0,1200);
    h_data->Sumw2();
    h_data->FillRandom("pol0",1000);
    h_data->FillRandom("pol1",100);
//     h_data->FillRandom("pol1",110);
  TH1F *h_bkg = new TH1F("h_bkg","h_bkg",10,0,1200);
    h_bkg->Sumw2();
    h_bkg->FillRandom("pol0",1000000);
    h_bkg->Scale(1000./h_bkg->Integral());
  TH1F *h_sig = new TH1F("h_sig","h_sig",10,0,1200);
    h_sig->Sumw2();
    h_sig->FillRandom("pol1",10000);
    h_sig->Scale(100./h_sig->Integral());
  
  // a histo syst...
  TH1F *h_bkg_jesUp = (TH1F*)h_bkg->Clone("h_bkg_jesUp");
    h_bkg_jesUp->SetBinContent(1,h_bkg_jesUp->GetBinContent(1)*1.2);
  TH1F *h_bkg_jesDo = (TH1F*)h_bkg->Clone("h_bkg_jesDo");
    h_bkg_jesDo->SetBinContent(1,h_bkg_jesDo->GetBinContent(1)*0.8);

    
  TtHFit *myFit = new TtHFit("FitExample0");

    // samples used in the following are declared here
    Sample *Signal = myFit->NewSample("Signal",SampleType::Signal);
      Signal->SetFillColor(kRed);
    Sample *Background = myFit->NewSample("Background",SampleType::Background);
    Sample *Data = myFit->NewSample("Data",SampleType::Data);

//     Region *CR1 = myFit->NewRegion("CR1");
//       SampleHist *CR1_Background = CR1->SetSampleHist(Background,h_CR1_Background);

    Region *SR1 = myFit->NewRegion("SR1");
      SampleHist *SR1_Data = SR1->SetSampleHist(Data,h_data);
      SampleHist *SR1_Background = SR1->SetSampleHist(Background,h_bkg);
        SR1_Background->AddOverallSyst("BkgXsec",0.10,-0.10);
        SR1_Background->AddHistoSyst("JES",h_bkg_jesUp,h_bkg_jesDo);
//         SR1_Background->AddSystematic("BkgXsec",SystsType::Overall);
      SampleHist *SR1_Signal = SR1->SetSampleHist(Signal,h_sig);
        SR1_Signal->AddNormFactor("mu",1,0,5);

  // print on the screen what's inside this TtHFit: regions, samples, systematics...
//   myFit->Print();
  
  // draw all pre-fit
//   myFit->DrawAndSaveAll();

//   SR1_Background->DrawSystPlot();
  
  // saves all in a root file for later usage
  myFit->WriteHistos();
  
  // export to RooStat:
  //  - creates xml
  //  - creates a workspace
  //  - make a quick fit
  myFit->SetPOI("mu");
//   myFit->SetStatErrorConfig(true,0.05,"Poisson");
  myFit->SetStatErrorConfig(true,0.00,"Gaussian");
  myFit->ToRooStat(true,false);
  
  
    
}