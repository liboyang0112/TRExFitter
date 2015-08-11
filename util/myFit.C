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
#include "TtHFitter/HistoTools.h"

#include <string>

// -------------------------------------------------------
// -------------------------------------------------------

void FitExample(string opt="h",string configFile="util/myFit.config",bool update=false){
    SetAtlasStyle();
    
    // interpret opt
    bool readHistograms  = opt.find("h")!=string::npos;
    bool readNtuples     = opt.find("n")!=string::npos;
    bool createWorkspace = opt.find("w")!=string::npos;
    bool doFit           = opt.find("f")!=string::npos;
    bool doLimit         = opt.find("l")!=string::npos;
    bool doSignificance  = opt.find("s")!=string::npos;
    bool drawPreFit      = opt.find("d")!=string::npos;
    bool drawPostFit     = opt.find("p")!=string::npos;
    
    TtHFit *myFit = new TtHFit();
    myFit->ReadConfigFile(configFile);
    
    // check compatibility between run option and config file
    if(readHistograms && myFit->fInputType!=TtHFit::HIST){
        std::cerr << "ERROR: Option \"h\" asked but no HISTO InputType speficied in the configuration file. Aborting." << std::endl;
        return;
    }
    if(readNtuples && myFit->fInputType!=TtHFit::NTUP){
        std::cerr << "ERROR: Option \"n\" asked but no NTUP InputType speficied in the configuration file. Aborting." << std::endl;
        return;
    }
      
    // -------------------------------------------------------

    if(readHistograms){
        myFit->ReadHistograms();
        myFit->Print();
        myFit->SmoothSystematics("all");
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        myFit->WriteHistos("",!update);
    }
    else if(readNtuples){
        myFit->ReadNtuples();
        myFit->Print();
        myFit->SmoothSystematics("all");
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        myFit->WriteHistos("",!update);
    }
    else{
        myFit->ReadHistos();
    }
      
    if(drawPreFit){
//         if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        myFit->DrawAndSaveAll();
        myFit->DrawSummary("log");
        myFit->BuildYieldTable();
        int nCols = 2;
        int nRows = 2;
        if(myFit->fNRegions>4){
            nCols = (int)sqrt(myFit->fNRegions);
            if(sqrt(myFit->fNRegions)>nCols) nCols++;
            nRows = (int)sqrt(myFit->fNRegions);
            if(nCols*nRows < myFit->fNRegions) nRows++;
        }
        myFit->DrawSignalRegionsPlot(nCols,nRows);
    }

    if(createWorkspace){
        myFit->DrawPruningPlot();
        myFit->SetLumiErr(0.);
        myFit->ToRooStat(true,true);
    }

    // use the external tool FitCrossCheckForLimits fir fitting
    if(doFit){
        myFit->Fit(); // with FitCrossCheckForLimits
//         myFit->PlotFittedNP();
//         myFit->PlotCorrelationMatrix();
    }
    
    if(doLimit){
        myFit->GetLimit();
    }
    
    if(doSignificance){
        myFit->GetSignificance();
    }
    
    if(drawPostFit){
        myFit->DrawAndSaveAll("post");
        myFit->PlotFittedNP();
        myFit->PlotCorrelationMatrix();
        myFit->DrawSummary("log post");
        myFit->BuildYieldTable("post");
    }
}

// -------------------------------------------------------
// -------------------------------------------------------
// main function
// -------------------------------------------------------
// -------------------------------------------------------

int main(int argc, char **argv){
  string opt="h";
  string config="util/myFit.config";
  bool update=false;
  
  if(argc>1) opt    = argv[1];
  if(argc>2) config = argv[2];
  if(argc>3) update = atoi(argv[3])!=0;

  // call the function
  FitExample(opt,config,update);
  
  return 0;
}
