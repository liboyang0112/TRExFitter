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
// "w" -> create workspace
// "f" -> fit
// "d" -> draw pre-fit
// "p" -> draw post-fit

void Fit4t(string opt="t",bool update=false){
  bool useTRF = true;
  
  TtHFitter::SetDebugLevel(1);    

  // interpret opt
  bool readNtuples     = opt.find("t")!=string::npos;
  bool createWorkspace = opt.find("w")!=string::npos;
  bool doFit           = opt.find("f")!=string::npos;
  bool doLimit         = opt.find("l")!=string::npos;
  bool drawPreFit      = opt.find("d")!=string::npos;
  bool drawPostFit     = opt.find("p")!=string::npos;
  bool systSmoothing   = opt.find("s")!=string::npos;
  
  // create the fit object
  TtHFit *myFit = new TtHFit();
    // ntuple stuff
    myFit->AddNtuplePath("/gpfs/atlas/public/TopMiniSLs/TRC_23/NTUP_ELE_V1_L2/nomin/");
    myFit->AddNtuplePath("/gpfs/atlas/public/TopMiniSLs/TRC_23/NTUP_MUON_V1_L2/nomin/");
    if(!useTRF) myFit->SetMCweight("FinalWeight * weight_BTag");
    else myFit->SetMCweight("FinalWeight");
    myFit->SetSelection("(lep_trigMatch&3)>0");
    myFit->SetNtupleName("mini");
    
  // create the samples
  Sample *data = myFit->NewSample("Data",SampleType::Data);
    data->SetTitle("Data 2012");
    data->AddNtupleFile("data");
    
  Sample *ttH = myFit->NewSample("ttH",SampleType::Background);
    ttH->SetTitle("t#bar{t}+H");
    ttH->SetFillColor(kGray);
    ttH->SetLineColor(kBlack);
    ttH->AddNtupleFile("ttH");
  Sample *ttV = myFit->NewSample("ttV",SampleType::Background);
    ttV->SetTitle("t#bar{t}+W/Z");
    ttV->SetFillColor(kGray+2);
    ttV->SetLineColor(kBlack);
    ttV->AddNtupleFile("ttbar_V");

  Sample *ttbar = myFit->NewSample("ttbar",SampleType::Background);
    ttbar->SetTitle("t#bar{t}+jets");
    ttbar->SetFillColor(kWhite);
    ttbar->SetLineColor(kBlack);
    ttbar->AddNtupleFile("ttbar_hdampHeraPdf");
//     ttbar->AddNtupleFile("ttbar_sherpa");
//     ttbar->AddNtupleFile("ttbar_MadGraph");

    // norm uncertainties
    ttbar->AddSystematic("ttXsec",SystType::Overall,0.10,-0.10);
//     ttbar->AddSystematic("ttXsec",SystType::Overall,0.0001,-0.0001);
    ttH->AddSystematic("tthXsec",SystType::Overall,0.50,-0.50);
    ttV->AddSystematic("ttvXsec",SystType::Overall,0.50,-0.50);

    // 24% uncertainty every jet bin above 6 (kind of Berend scaling)
    Systematic *ttJetRatio = ttbar->AddSystematic("ttJetRatio",SystType::Histo);
      ttJetRatio->fWeightSufUp = "(1+(jet_n - 4)*0.05)";
      ttJetRatio->fWeightSufDown = "(1-(jet_n - 4)*0.05)";

    Systematic *ttMCsherpa = ttbar->AddSystematic("ttMCsherpa",SystType::Histo);
      ttMCsherpa->fNtupleFilesUp.push_back("ttbar_sherpa");
    Systematic *ttMCmg = ttbar->AddSystematic("ttMCmg",SystType::Histo);
      ttMCmg->fNtupleFilesUp.push_back("ttbar_MadGraph");
      
  Sample *sig = myFit->NewSample("Signal",SampleType::Signal);
    sig->SetTitle("SM t#bar{t}t#bar{t}");
    sig->SetFillColor(kRed);
    sig->SetLineColor(kRed);
    sig->AddNormFactor("SigXsecOverSM",1,0,100);
    sig->AddNtupleFile("tttt_SM");
    
    
  // create fit regions
  Region *SR_ge10jge4b = myFit->NewRegion("SR_ge10jge4b");
    if(!useTRF) SR_ge10jge4b->AddSelection("jet_n>=10 && nJet_tagged>=4");
    else{
      SR_ge10jge4b->AddSelection("jet_n>=10 && (channelNumber!=runNumber || nJet_tagged>=4)");
      SR_ge10jge4b->AddMCweight("weight_TRF_ge4b");
    }
    SR_ge10jge4b->SetVariable("(HTj+lep_pt+met_et)/1e3",7,400,1800);
    SR_ge10jge4b->SetVariableTitle("H_{T} [GeV]");
    SR_ge10jge4b->SetLabel("e/#mu + #geq10 j, #geq4 b", "#geq10 j, #geq4 b");

  Region *SR_ge10j3b = myFit->NewRegion("SR_ge10j3b");
    if(!useTRF) SR_ge10j3b->AddSelection("jet_n>=10 && nJet_tagged==3");
    else{
      SR_ge10j3b->AddSelection("jet_n>=10 && (channelNumber!=runNumber || nJet_tagged==3)");
      SR_ge10j3b->AddMCweight("weight_TRF_3b");
    }
    SR_ge10j3b->SetVariable("(HTj+lep_pt+met_et)/1e3",7,400,1800);
    SR_ge10j3b->SetVariableTitle("H_{T} [GeV]");
    SR_ge10j3b->SetLabel("e/#mu + #geq10 j, 3 b", "#geq10 j, 3 b");
    
  Region *SR_9jge4b = myFit->NewRegion("SR_9jge4b");
    if(!useTRF) SR_9jge4b->AddSelection("jet_n==9 && nJet_tagged>=4");
    else{
      SR_9jge4b->AddSelection("jet_n==9 && (channelNumber!=runNumber || nJet_tagged>=4)");
      SR_9jge4b->AddMCweight("weight_TRF_ge4b");
    }
    SR_9jge4b->SetVariable("(HTj+lep_pt+met_et)/1e3",7,400,1800);
    SR_9jge4b->SetVariableTitle("H_{T} [GeV]");
    SR_9jge4b->SetLabel("e/#mu + 9 j, #geq4 b", "9 j, #geq3 b");

  // control region
  Region *CR_ge4j2b = myFit->NewRegion("CR_ge4j2b");
    if(!useTRF) CR_ge4j2b->AddSelection("jet_n>=4 && nJet_tagged==2");
    else{
      CR_ge4j2b->AddSelection("jet_n>=4 && (channelNumber!=runNumber || nJet_tagged==2)");
      CR_ge4j2b->AddMCweight("weight_TRF_2b");
    }
    CR_ge4j2b->SetVariable("jet_n",7,4,11);
    CR_ge4j2b->SetVariableTitle("Number of jets");
    CR_ge4j2b->SetLabel("e/#mu + #geq4 j, 2 b", "#geq4 j, 2 b");
    
    
  //
  // do actual things
  //
  if(readNtuples){
    myFit->ReadNtuples();
    myFit->Print();
    myFit->WriteHistos("MyMeasurement_histos.root",!update);
  }
  else{
    myFit->ReadHistos("MyMeasurement_histos.root");
  }
  
  if(systSmoothing){
    myFit->SmoothSystematics("all");
    myFit->WriteHistos("MyMeasurement_histos.root",!update);
  }
  
  if(drawPreFit){
    myFit->DrawAndSaveAll();
//     myFit->DrawAndSaveAll("blind");
    myFit->DrawSystPlots();
//     Region *regions[4];
//     regions[0] = CR_ge4j2b;
//     regions[1] = SR_9jge4b;
//     regions[2] = SR_ge10j3b;
//     regions[3] = SR_ge10jge4b;
//     myFit->DrawSignalRegionsPlot(2,2,regions);
  }

  if(createWorkspace){
    myFit->SetPOI("SigXsecOverSM");
//     myFit->SetLumiErr(0.037);
    myFit->SetLumiErr(0.);
//     myFit->SetStatErrorConfig(true,0.05);
    myFit->SetStatErrorConfig(false,0.05);
    myFit->ToRooStat(true,true);
//     myFit->ToRooStat(true,false);
  }

  // use the external tool FitCrossCheckForLimits fir fitting
  if(doFit){
    // PlotHistosBeforeFit=0,
    // PlotMorphingControlPlots=1, 
    // PlotHistosAfterFitEachSubChannel=2, 
    // PlotHistosAfterFitGlobal=3, 
    // PlotsNuisanceParametersVSmu=4, 
    // PlotsStatisticalTest=5
    int algo = 3;
//     int algo = 0;
    string workspace = "results/MyMeasurement_combined_MyMeasurement_model.root";
    string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 0,\"%s\",\"./xcheckResults/\",\"combined\",\"ModelConfig\",\"obsData\")'",algo,workspace.c_str());
    gSystem->Exec(cmd.c_str());
    //
    // plot the NP fit plot
    gSystem->Exec("python plotNP.py \"xcheckResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt\"");
  }
  
  if(doLimit){
    string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/MyMeasurement_combined_MyMeasurement_model.root\",\"combined\",\"ModelConfig\",\"obsData\")'";
    gSystem->Exec(cmd.c_str());
  }
  
  if(drawPostFit){
    myFit->ReadFitResults("xcheckResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt");
//     SR_ge10jge4b->DrawPostFit(myFit->fFitResults)->SaveAs("Test.png");
    myFit->DrawAndSaveAll("post");
  }
  
}
 
