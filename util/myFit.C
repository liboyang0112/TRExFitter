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
#include "TtHFitter/MultiFit.h"

#include <string>

// -------------------------------------------------------
// -------------------------------------------------------

void FitExample(string opt="h",string configFile="util/myl2tau.config",string options=""){
    SetAtlasStyle();
    
    // multi-fit
    bool isMultiFit      = opt.find("m")!=string::npos;    
    if(isMultiFit){
        MultiFit *myMultiFit = new MultiFit();
        myMultiFit->ReadConfigFile(configFile,options);
        myMultiFit->ComparePulls();
        myMultiFit->ComparePOI("SigXsecOverSM");
        myMultiFit->CompareLimit();
        return;
    }
    
    // interpret opt
    bool readHistograms  = opt.find("h")!=string::npos;
    bool readNtuples     = opt.find("n")!=string::npos;
    bool createWorkspace = opt.find("w")!=string::npos;
    bool doFit           = opt.find("f")!=string::npos;
    bool doRanking       = opt.find("r")!=string::npos;
    bool doLimit         = opt.find("l")!=string::npos;
    bool doSignificance  = opt.find("s")!=string::npos;
    bool drawPreFit      = opt.find("d")!=string::npos;
    bool drawPostFit     = opt.find("p")!=string::npos;
    bool drawSeparation  = opt.find("a")!=string::npos;
    
    TtHFit *myFit = new TtHFit();
    myFit->ReadConfigFile(configFile,options);
    
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
        myFit->WriteHistos();
    }
    else if(readNtuples){
        myFit->ReadNtuples();
        myFit->Print();
        myFit->SmoothSystematics("all");
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        myFit->WriteHistos();
    }
    else{
        if(drawPreFit || drawPostFit || createWorkspace || drawSeparation) myFit->ReadHistos();
    }

    if(createWorkspace){
        myFit->DrawPruningPlot();
        myFit->SetLumiErr(0.);
        myFit->ToRooStat(true,true);
    }

    // use the external tool FitCrossCheckForLimits fir fitting
    if(doFit){
        myFit->Fit();
        myFit->PlotFittedNP();
        myFit->PlotCorrelationMatrix();
    }
    if(doRanking){
        if(myFit->fRankingOnly!="plot")  myFit->ProduceNPRanking( myFit->fRankingOnly );
        if(myFit->fRankingOnly=="all" || myFit->fRankingOnly=="plot")  myFit->PlotNPRanking();
    }
    
    if(doLimit){
        myFit->GetLimit();
    }
    
    if(doSignificance){
        myFit->GetSignificance();
    }
    
    if(drawPreFit){
        myFit->DrawAndSaveAll();
        myFit->DrawSummary("log");
        myFit->DrawSummary("logvalid");
        myFit->BuildYieldTable();
        myFit->PrintSystTables();
        int nCols = 2;
        int nRows = 2;
        if(myFit->fNRegions>4){
            nCols = (int)sqrt(myFit->fNRegions);
            if(sqrt(myFit->fNRegions)>nCols) nCols++;
            nRows = (int)sqrt(myFit->fNRegions);
            if(nCols*nRows < myFit->fNRegions) nRows++;
        }
        myFit->DrawSignalRegionsPlot(nCols,nRows);
        myFit->DrawPieChartPlot("pre",nCols,nRows);
    }
    
    if(drawPostFit){
        myFit->DrawAndSaveAll("post");
        myFit->DrawSummary("log post");
        myFit->DrawSummary("log post valid");
        myFit->BuildYieldTable("post");
        int nCols = 2;
        int nRows = 2;
        if(myFit->fNRegions>4){
            nCols = (int)sqrt(myFit->fNRegions);
            if(sqrt(myFit->fNRegions)>nCols) nCols++;
            nRows = (int)sqrt(myFit->fNRegions);
            if(nCols*nRows < myFit->fNRegions) nRows++;
        }
        myFit->DrawPieChartPlot("post",nCols,nRows);
    }

    if(drawSeparation){
        myFit->DrawAndSaveSeparationPlots();
    //    myFit->ListOfBestSeparationVariables(); // for the future list of best separation variables
    //    myFit->ListOfBestDataMCVariables();     // for the future list of best data-mc agreement variables based on KS test
    }
  
}

// -------------------------------------------------------
// -------------------------------------------------------
// main function
// -------------------------------------------------------
// -------------------------------------------------------

int main(int argc, char **argv){
  string opt="h";
  string config="util/myl2tau,config";
  string options="";
  
  if(argc>1) opt     = argv[1];
  if(argc>2) config  = argv[2];
  if(argc>3) options = argv[3];

  // call the function
  FitExample(opt,config,options);
  
  return 0;
}
