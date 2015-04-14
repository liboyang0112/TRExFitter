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
#include "TtHFitter/ConfigParser.h"

#include <string>

// -------------------------------------------------------

void FitExample_fromHist(string opt="h",bool update=false){
  TtHFitter::SetDebugLevel(1); 

  // options
  bool useMCstat = false;
  
  SetAtlasStyle();
  
  // interpret opt
  bool readHistograms  = opt.find("h")!=string::npos;
  bool createWorkspace = opt.find("w")!=string::npos;
  bool doFit           = opt.find("f")!=string::npos;
  bool doLimit         = opt.find("l")!=string::npos;
  bool drawPreFit      = opt.find("d")!=string::npos;
  bool drawPostFit     = opt.find("p")!=string::npos;
    
  // Read the config file
  ConfigParser *myConfig = new ConfigParser();
  myConfig->ReadFile("util/myFit.config");
  ConfigSet *cs; // to store stuff later
  
  // set the stuff accordingly...
  int type;
  
  // set fit
  cs = myConfig->GetConfigSet("Fit");
  TtHFit *myFit = new TtHFit( cs->GetValue() );
  myFit->AddHistoPath( cs->Get("HistoPath") );
  if(cs->Get("LumiScale")!="") myFit -> SetLumi( atof(cs->Get("LumiScale").c_str()) );
  if(cs->Get("FitType")!=""){
      if(cs->Get("FitType")=="ControlSignalRegion") myFit -> SetFitType(TtHFit::ControlSignalRegion);
      else if(cs->Get("FitType")=="ControlRegion")  myFit -> SetFitType(TtHFit::ControlRegion);
      else {
          std::cerr << "Unknown FitType argument : " << cs->Get("FitType") << std::endl;
      }
  }
  if(cs->Get("SystPruningShape")!="") myFit->fThresholdSystPruning_Shape = atof(cs->Get("SystPruningShape").c_str());
  if(cs->Get("SystPruningNorm")!="") myFit->fThresholdSystPruning_Normalisation = atof(cs->Get("SystPruningNorm").c_str());

  // set regions
  int nReg = 0;
  Region *reg;
  while(true){
    cs = myConfig->GetConfigSet("Region",nReg);
    if(cs==0x0) break;
    reg = myFit->NewRegion(cs->GetValue());
    reg->SetHistoName(cs->Get("HistoName"));
    reg->SetVariableTitle(cs->Get("VariableTitle"));
    reg->SetLabel(cs->Get("Label"),cs->Get("ShortLabel"));
    //Potential rebinning
    if(cs->Get("Rebin")!="") reg -> Rebin(atoi(cs->Get("Rebin").c_str()));
    if(cs->Get("Binning")!=""){
        std::vector < string > vec_bins = Vectorize(cs->Get("Binning"), ',');
        const int nBounds = vec_bins.size();
        double bins[nBounds];
        for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
            bins[iBound] = atof(vec_bins[iBound].c_str());
        }
        reg -> SetBinning(nBounds-1,bins);
    }
    nReg++;
  }
  
  // set samples 
  int nSmp = 0;
  Sample *smp;
  while(true){
    cs = myConfig->GetConfigSet("Sample",nSmp);
    if(cs==0x0) break;
    type = SampleType::Background;
    if(cs->Get("Type")=="signal") type = SampleType::Signal;
    if(cs->Get("Type")=="data") type = SampleType::Data;
    smp = myFit->NewSample(cs->GetValue(),type);
    smp->SetTitle(cs->Get("Title"));
    smp->AddHistoFile(cs->Get("HistoFile"));
    if(cs->Get("FillColor")!="")
      smp->SetFillColor(atoi(cs->Get("FillColor").c_str()));
    if(cs->Get("LineColor")!="")
      smp->SetLineColor(atoi(cs->Get("LineColor").c_str()));
    if(cs->Get("NormFactor")!="")
      smp->AddNormFactor(
        Vectorize(cs->Get("NormFactor"),',')[0],
        atof(Vectorize(cs->Get("NormFactor"),',')[1].c_str()),
        atof(Vectorize(cs->Get("NormFactor"),',')[2].c_str()),
        atof(Vectorize(cs->Get("NormFactor"),',')[3].c_str())
      );
    // ...
    nSmp++;
  }
  
  
  // set systs
  int nSys = 0;
  Systematic *sys;
  Sample *sam;
  while(true){
    cs = myConfig->GetConfigSet("Systematic",nSys);
    if(cs==0x0) break;
    string samples_str = cs->Get("Samples");
    if(samples_str=="") samples_str = "all";
    vector<string> samples = Vectorize(samples_str,',');
    type = SystType::Histo;
    if(cs->Get("Type")=="overall" || cs->Get("Type")=="Overall")
      type = SystType::Overall;
    for(int i_smp=0;i_smp<myFit->fNSamples;i_smp++){
      sam = myFit->fSamples[i_smp];
      if(sam->fType == SampleType::Data) continue;
      if(samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() ){
//         cout << "Adding syst " << cs->GetValue() << " for sample "<< sam->fName << endl;
        sys = sam->AddSystematic(cs->Get("Title"),type);
        if(type==SystType::Histo){
          sys->fHistoNameSufUp = cs->Get("HistoNameSufUp");
          sys->fHistoNameSufDown = cs->Get("HistoNameSufDown");
          // ...
        }
        else if(type==SystType::Overall){
          sys->fOverallUp = atof( cs->Get("NormUp").c_str() );
          sys->fOverallDown = atof( cs->Get("NormDown").c_str() );
        }
      }
    }
    // ...
    nSys++;
  }
  
  
//     Systematic *JES_bkg1 = bkg1->AddSystematic("JES",SystType::Histo);
//       JES_bkg1->fHistoNameSufUp = "_jesUp";
//       JES_bkg1->fHistoNameSufDown = "_jesDown";
    
  //
  // do actual things
  //

  myFit->SetStatErrorConfig(useMCstat,0.05);

  if(readHistograms){
    myFit->ReadHistograms();
    myFit->Print();
    myFit->WriteHistos("",!update);
  }
  else{
    myFit->ReadHistos();
  }
    
  if(drawPreFit){
    myFit->DrawAndSaveAll();
    myFit->DrawSummary();
    myFit->DrawSignalRegionsPlot(2,2);
  }

  if(createWorkspace){
    myFit->SetPOI("SigXsecOverSM");
    myFit->SetLumiErr(0.);
    myFit->ToRooStat(true,true);
  }

  // use the external tool FitCrossCheckForLimits fir fitting
  if(doFit){
    myFit->Fit(); // with FitCrossCheckForLimits
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
