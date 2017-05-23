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

void FitExample(string opt="h",string configFile="util/myFit.config",string options=""){
    SetAtlasStyle();
    
    // interpret opt
    bool readHistograms  = opt.find("h")!=string::npos;
    bool readNtuples     = opt.find("n")!=string::npos;
    bool rebinAndSmooth  = opt.find("b")!=string::npos; // new: separated from the input creation
    bool createWorkspace = opt.find("w")!=string::npos;
    bool doFit           = opt.find("f")!=string::npos;
    bool doRanking       = opt.find("r")!=string::npos;
    bool doLimit         = opt.find("l")!=string::npos;
    bool doSignificance  = opt.find("s")!=string::npos;
    bool drawPreFit      = opt.find("d")!=string::npos;
    bool drawPostFit     = opt.find("p")!=string::npos;
    bool drawSeparation  = opt.find("a")!=string::npos;
    
    // multi-fit
    bool isMultiFit      = opt.find("m")!=string::npos;
    if(isMultiFit){
        MultiFit *myMultiFit = new MultiFit();
        myMultiFit->ReadConfigFile(configFile,options);
        //
        if(myMultiFit->fCombine){
            if(createWorkspace){
                myMultiFit->SaveCombinedWS();
            }
            if(doFit){
                myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName );
            }
            if(doLimit){
                myMultiFit->GetCombinedLimit( myMultiFit->fDataName );
            }
            if(doSignificance){
                myMultiFit->GetCombinedSignificance( myMultiFit->fDataName );
            }
            if(doRanking){
                if(myMultiFit->fRankingOnly!="plot")  myMultiFit->ProduceNPRanking( myMultiFit->fRankingOnly );
                if(myMultiFit->fRankingOnly=="all" || myMultiFit->fRankingOnly=="plot")  myMultiFit->PlotNPRanking();
            }
        }
        //
        if(myMultiFit->fCompare){
            if(myMultiFit->fComparePulls){
                for(unsigned int i_cat=0;i_cat<myMultiFit->fNPCategories.size();i_cat++){
                    myMultiFit->ComparePulls(myMultiFit->fNPCategories[i_cat]);
                }
                myMultiFit->CompareNormFactors("");
            }
            if(myMultiFit->fPlotCombCorrMatrix) myMultiFit->PlotCombinedCorrelationMatrix();
            if(myMultiFit->fComparePOI)    myMultiFit->ComparePOI(myMultiFit->fPOI);
            if(myMultiFit->fCompareLimits) myMultiFit->CompareLimit();
        }
        //
        if(myMultiFit->fPlotSoverB)    myMultiFit->PlotSummarySoverB();
        //
        return;
    }
    
    // proceed if not multi-fit
    
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
        myFit->CreateRootFiles();
        myFit->ReadHistograms();
        myFit->Print();
        myFit->CorrectHistograms(); // apply rebinning, smoothing etc...
        myFit->WriteHistos();
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TtHFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();
    }
    else if(readNtuples){
        myFit->CreateRootFiles();
        myFit->ReadNtuples();
        myFit->Print();
        myFit->CorrectHistograms(); // apply rebinning, smoothing etc...
        myFit->WriteHistos();
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TtHFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();
    }
    else{
        if(drawPreFit || drawPostFit || createWorkspace || drawSeparation || rebinAndSmooth) myFit->ReadHistos();
    }
    
    // new
    if(rebinAndSmooth){
        bool udpate = myFit->fUpdate;
        myFit->fUpdate = true;
        myFit->CreateRootFiles();  // ?
        myFit->fUpdate = udpate;
        myFit->CorrectHistograms(); // apply rebinning, smoothing etc...
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TtHFitter::SYSTDATAPLOT) myFit->DrawSystPlotsSumSamples();
        myFit->WriteHistos();
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
        if(myFit->fRankingOnly=="all" || myFit->fRankingOnly=="plot")  myFit->PlotNPRankingManager();
    }
    
    if(doLimit){
        myFit->GetLimit();
    }
    
    if(doSignificance){
        myFit->GetSignificance();
    }
    
    if(drawPreFit){
        if(TtHFitter::OPTION["PrefitRatioMax"]==2){
            myFit->DrawAndSaveAll("prefit");
            myFit->DrawSummary("log prefit");
            myFit->DrawSummary("log valid prefit");
        }
        else{
            myFit->DrawAndSaveAll();
            myFit->DrawSummary("log");
            myFit->DrawSummary("logvalid");
        }
        myFit->BuildYieldTable();
        for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
            myFit->BuildYieldTable("",myFit->fRegionGroups[i_gr]);
        }
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
        for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
            myFit->BuildYieldTable("post",myFit->fRegionGroups[i_gr]);
        }
        myFit->PrintSystTables("post");
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
  string config="util/myFit,config";
  string options="";
  
  if(argc>1) opt     = argv[1];
  if(argc>2) config  = argv[2];
  if(argc>3) options = argv[3];

  // call the function
  FitExample(opt,config,options);
  
  return 0;
}
