#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

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

void Fit4t(bool readNtuples=true,bool update=false){
    

  // create the fit object
  TtHFit *myFit = new TtHFit();
    myFit->SetPOI("SigXsecOverSM");
    // ntuple stuff
    myFit->AddNtuplePath("/gpfs/atlas/public/TopMiniSLs/TRC_23/NTUP_ELE_V1_L2/nomin/");
//     myFit->AddNtuplePath("/gpfs/atlas/public/TopMiniSLs/TRC_23/NTUP_MUON_V1_L2/nomin/");
    myFit->SetMCweight("FinalWeight * weight_BTag");
    myFit->SetSelection("(lep_trigMatch&3)>0");
    myFit->SetTreeName("mini");
// //     myFit->SetNominalPathSuf("");

  // create norm the factor for the POI
  NormFactor *SigXsecOverSM = new NormFactor("SigXsecOverSM",1,0,100);
    
  // create systematics
  Systematic *TtXsec = new Systematic("tt_xsec");
    TtXsec->fOverallUp = 0.10;
    TtXsec->fOverallDown = -0.10;
//   Systematic *TtMC = new Systematic("tt_mc");
//     TtMC->fNtupleFilesUp.push_back("ttbar_MadGraph.root");
//     TtMC->fNtupleFilesDown.push_back("ttbar_hdampHeraPdf.root");
  Systematic *TtNjet = new Systematic("tt_njet");
    TtNjet->fWeightUp = "(1+(jet_n - 4)*0.05)";
    TtNjet->fWeightDown = "(1-(jet_n - 4)*0.05)";
  
  // create the samples
  Sample *data = myFit->NewSample("Data");
    data->SetTitle("Data 2012");
    data->SetIsData();
    data->AddNtupleName("data.root");
    
  Sample *ttbar = myFit->NewSample("ttbar");
    ttbar->SetTitle("t#bar{t}+jets");
    ttbar->SetFillColor(kWhite);
    ttbar->SetLineColor(kBlack);
// //     ttbar->AddNtupleName("ttbar_hdampHeraPdf.root");
    ttbar->AddNtupleName("ttbar_sherpa.root");
    ttbar->AddSystematic(TtXsec);
//     ttbar->AddSystematic(TtMC);
    ttbar->AddSystematic(TtNjet);
//     ttbar->AddOverallSyst("tt_xsec",0.10,-0.10);
//     // 24% uncertainty every jet bin above 6 (kind of Berend scaling)
//     ttbar->AddWeightSyst("tt_jet_ratio","(1+(jet_n - 6)*0.24)","(1-(jet_n - 6)*0.24)");
// //     ttbar->AddNtupleSyst("tt_mc","ttbar_MadGraph.root","ttbar_hdampHeraPdf.root");
//     
  Sample *sig = myFit->NewSample("Signal");
    sig->SetIsSignal();
    sig->SetFillColor(kRed);
    sig->SetLineColor(kRed);
    sig->AddNormFactor(SigXsecOverSM);
    sig->AddNtupleName("tttt_SM.root");
//     
//     
  // create fit regions
  Region *CR_ge4j2b = myFit->NewRegion("CR_ge4j2b");
    CR_ge4j2b->AddSelection("jet_n>=4 && nJet_tagged==2");
    CR_ge4j2b->SetVariable("jet_n",7,4,11);
    CR_ge4j2b->SetLabel("e/#mu + #geq4 jets, 2 b-tags");
  Region *CR_ge4j3b = myFit->NewRegion("CR_ge4j3b");
    CR_ge4j3b->AddSelection("jet_n>=4 && nJet_tagged==3");
    CR_ge4j3b->SetVariable("jet_n",7,4,11);
    CR_ge4j3b->SetLabel("e/#mu + #geq4 jets, 3 b-tags");
  Region *CR_ge4jge4b = myFit->NewRegion("CR_ge4jge4b");
    CR_ge4jge4b->AddSelection("jet_n>=4 && nJet_tagged>=4");
    CR_ge4jge4b->SetVariable("jet_n",7,4,11);
    CR_ge4jge4b->SetLabel("e/#mu + #geq4 jets, #geq4 b-tags");
// 
//     
// //   Region *SR_7j3b = myFit->NewRegion("SR_7j3b");
// //     SR_7j3b->AddSelection("jet_n==7 && nJet_tagged==3");
// //     SR_7j3b->SetVariable("(HTj+lep_pt+met_et)/1e3",8,400,1200);
// // 
// //   Region *SR_8j3b = myFit->NewRegion("SR_8j3b");
// //     SR_8j3b->AddSelection("jet_n==8 && nJet_tagged==3");
// //     SR_8j3b->SetVariable("(HTj+lep_pt+met_et)/1e3",8,400,1200);
// // 
//   Region *SR_9j3b = myFit->NewRegion("SR_9j3b");
//     SR_9j3b->AddSelection("jet_n==9 && nJet_tagged==3");
//     SR_9j3b->SetVariable("(HTj+lep_pt+met_et)/1e3",8,400,1200);
// 
//   Region *SR_ge10j3b = myFit->NewRegion("SR_ge10j3b");
//     SR_ge10j3b->AddSelection("jet_n>=10 && nJet_tagged==3");
//     SR_ge10j3b->SetVariable("(HTj+lep_pt+met_et)/1e3",8,400,1200);

// 
// //   Region *SR_7jge4b = myFit->NewRegion("SR_7jge4b");
// //     SR_7jge4b->AddSelection("jet_n==7 && nJet_tagged>=4");
// //     SR_7jge4b->SetVariable("(HTj+lep_pt+met_et)/1e3",8,400,1200);
// // 
// //   Region *SR_8jge4b = myFit->NewRegion("SR_8jge4b");
// //     SR_8jge4b->AddSelection("jet_n==8 && nJet_tagged>=4");
// //     SR_8jge4b->SetVariable("(HTj+lep_pt+met_et)/1e3",8,400,1200);

//   Region *SR_9jge4b = myFit->NewRegion("SR_9jge4b");
//     SR_9jge4b->AddSelection("jet_n==9 && nJet_tagged>=4");
//     SR_9jge4b->SetVariable("(HTj+lep_pt+met_et)/1e3",4,400,1200);

//   Region *SR_ge10jge4b = myFit->NewRegion("SR_ge10jge4b");
//     SR_ge10jge4b->AddSelection("jet_n>=10 && nJet_tagged>=4");
//     SR_ge10jge4b->SetVariable("(HTj+lep_pt+met_et)/1e3",4,400,1200);
// 
//     
  if(readNtuples){
    myFit->ReadAll(true);
    myFit->WriteHistos("MyMeasurement_histos.root",!update);
  }
  else{
    myFit->ReadAll(false);
  }
  
  CR_ge4j2b->SetVariableTitle("Number of jets");
  CR_ge4j3b->SetVariableTitle("Number of jets");
  CR_ge4jge4b->SetVariableTitle("Number of jets");
  
  myFit->DrawAndSaveAll();
  
//   if(!readNtuples)  myFit->ToRooStat(false);
//   if(!readNtuples)  myFit->ToRooStat(true);
//   if(!readNtuples)  myFit->ToRooStat(true,false);
  
  
// return;  
  
  // Fit Results
  FitResults* fitRes = new FitResults();
  NuisParameter* np0 = new NuisParameter("SigXsecOverSM");
    np0->fFitValue = 100;
    np0->fPostFitUp = 200;
    np0->fPostFitDown = -100;
  NuisParameter* np1 = new NuisParameter("tt_xsec");
    np1->fFitValue = 1.53;// 0.41;
    np1->fPostFitUp = 0.043;//0.83;
    np1->fPostFitDown = -0.043;//-0.83;
  NuisParameter* np2 = new NuisParameter("tt_njet");
    np2->fFitValue = -1.6;//-0.0018;
    np2->fPostFitUp = 0.06;//0.03;
    np2->fPostFitDown = -0.06;//-0.03;
  CorrelationMatrix* matrix = new CorrelationMatrix();
    matrix->SetCorrelation("tt_xsec","tt_xsec",1);
    matrix->SetCorrelation("tt_xsec","tt_njet",-0.005);
    matrix->SetCorrelation("tt_njet","tt_njet",1);
    matrix->SetCorrelation("tt_njet","tt_xsec",-0.005);
  fitRes->AddNuisPar(np0);
  fitRes->AddNuisPar(np1);
  fitRes->AddNuisPar(np2);
  fitRes->fCorrMatrix = matrix;
  
  myFit->fFitResults = fitRes;
  
  CR_ge4j2b->DrawPostFit(fitRes) -> SaveAs("CR_ge4j2b_PostFit.png");
  CR_ge4j3b->DrawPostFit(fitRes) -> SaveAs("CR_ge4j3b_PostFit.png");
  CR_ge4jge4b->DrawPostFit(fitRes) -> SaveAs("CR_ge4jge4b_PostFit.png");
  
}
 
