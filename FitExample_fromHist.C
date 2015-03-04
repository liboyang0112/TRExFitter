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

// opt:
// "t" -> read ntuples
// "h" -> read histograms
// "w" -> create workspace
// "f" -> fit
// "d" -> draw pre-fit
// "p" -> draw post-fit

void FitExample_fromHist(string opt="h",bool update=false){
  TtHFitter::SetDebugLevel(1);    

  // interpret opt
  bool readHistograms  = opt.find("h")!=string::npos;
  bool createWorkspace = opt.find("w")!=string::npos;
  bool doFit           = opt.find("f")!=string::npos;
  bool doLimit         = opt.find("l")!=string::npos;
  bool drawPreFit      = opt.find("d")!=string::npos;
  bool drawPostFit     = opt.find("p")!=string::npos;
  bool systSmoothing   = opt.find("s")!=string::npos;
  
  // create the fit object
  TtHFit *myFit = new TtHFit("FitExample1");
    // histogram stuff
    myFit->AddHistoPath("ExampleInputs");
    
  // create the samples
  Sample *data = myFit->NewSample("Data",SampleType::Data);
    data->SetTitle("Data 2012");
    data->AddHistoFile("data");
    
  Sample *bkg1 = myFit->NewSample("Bkg1",SampleType::Background);
    bkg1->SetTitle("Backgr.1");
    bkg1->SetFillColor(kYellow);
    bkg1->SetLineColor(kBlack);
    bkg1->AddHistoFile("bkg1");

    Systematic *JES_bkg1 = bkg1->AddSystematic("JES",SystType::Histo);
      JES_bkg1->fHistoNameSufUp = "_jesUp";
      JES_bkg1->fHistoNameSufDown = "_jesDown";

  Sample *bkg2 = myFit->NewSample("Bkg2",SampleType::Background);
    bkg2->SetTitle("Backg.2");
    bkg2->SetFillColor(kBlue-9);
    bkg2->SetLineColor(kBlack);
    bkg2->AddHistoFile("bkg2");
    bkg2->AddSystematic("BkgXsec",SystType::Overall,0.10,-0.10);
      
  Sample *sig = myFit->NewSample("Signal",SampleType::Signal);
    sig->SetTitle("Signal");
    sig->SetFillColor(kRed);
    sig->SetLineColor(kRed);
    sig->AddNormFactor("SigXsecOverSM",1,0,100);
    sig->AddHistoFile("sig");
    
  // signal region
  Region *SR_1 = myFit->NewRegion("SR_1");
    SR_1->SetHistoName("HTj");
    SR_1->SetVariableTitle("H_{T} [GeV]");
    SR_1->SetLabel("Signal Region 1", "SR 1");
    
  //
  // do actual things
  //
  if(readHistograms){
    myFit->ReadHistograms();
    myFit->Print();
    myFit->WriteHistos("FitExample1_histos.root",!update);
  }
  else{
    myFit->ReadHistos("FitExample1_histos.root");
  }
  
  if(systSmoothing){
    myFit->SmoothSystematics("all");
    myFit->WriteHistos("FitExample1_histos.root",!update);
  }
  
  if(drawPreFit){
    myFit->DrawAndSaveAll();
    myFit->DrawSummary();
    myFit->DrawSystPlots();
    myFit->DrawSignalRegionsPlot(2,2);
  }

  if(createWorkspace){
    myFit->SetPOI("SigXsecOverSM");
    myFit->SetLumiErr(0.);
    myFit->SetStatErrorConfig(false,0.05);
    myFit->ToRooStat(true,true);
  }

  // use the external tool FitCrossCheckForLimits fir fitting
  if(doFit){
//     myFit->Fit(); // with FitCrossCheckForLimits
    myFit->PlotFittedNP();
  }
  
  if(doLimit){
    myFit->GetLimit();
  }
  
  if(drawPostFit){
    myFit->DrawAndSaveAll("post");
  }
  
}
 
