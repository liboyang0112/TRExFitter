#include "TtHFitter/Common.h"

#include "TtHFitter/NuisParameter.h"
#include "TtHFitter/CorrelationMatrix.h"
#include "TtHFitter/FitResults.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/SystematicHist.h"
#include "TtHFitter/NormFactor.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/Region.h"
#include "TtHFitter/SampleHist.h"
#include "TtHFitter/TtHFit.h"
#include "TtHFitter/TthPlot.h"

// -------------------------------------------------------

void FitExample_fromHist(string opt="h",bool update=false){
  TtHFitter::SetDebugLevel(1); 
  
//   gROOT->ProcessLine(".L AtlasStyle.C");
//   gROOT->ProcessLine(".L AtlasLabels.C");
//   gROOT->ProcessLine(".L AtlasUtils.C");

  SetAtlasStyle();
  
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

// -------------------------------------------------------
// main function

int main(int argc, char **argv){
  string opt="h";
  bool update=false;
  
  if(argc>1) opt    = argv[1];
  if(argc>2) update = atoi(argv[2])!=0;

  // call the function
  FitExample_fromHist(opt,update);
  
  return 0;
}
