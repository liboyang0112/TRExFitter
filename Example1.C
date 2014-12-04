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


// Simplest example

void Example1(){

  TtHFitter::SetDebugLevel(0);
  
  TH1F *h_data = new TH1F("h_data","h_data",10,0,1200);
    h_data->Sumw2();
    h_data->FillRandom("pol0",100);
    h_data->FillRandom("pol1",10);
//     h_data->FillRandom("pol1",110);
  TH1F *h_bkg = new TH1F("h_bkg","h_bkg",10,0,1200);
    h_bkg->Sumw2();
    h_bkg->FillRandom("pol0",100000);
    h_bkg->Scale(100./h_bkg->Integral());
  TH1F *h_sig = new TH1F("h_sig","h_sig",10,0,1200);
    h_sig->Sumw2();
    h_sig->FillRandom("pol1",1000);
    h_sig->Scale(10./h_sig->Integral());
  
  
  TtHFit *myFit = new TtHFit();

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
//         SR1_Background->AddSystematic("BkgXsec",SystsType::Overall);
      SampleHist *SR1_Signal = SR1->SetSampleHist(Signal,h_sig);
        SR1_Signal->AddNormFactor("mu",1,0,10);

  myFit->Print();
  
  myFit->DrawAndSaveAll();
    
}