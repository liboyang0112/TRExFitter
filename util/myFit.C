#include "TtHFitter/Common.h"

#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/ConfigReader.h"
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

void FitExample(std::string opt="h",std::string configFile="util/myFit.config",std::string options=""){
    SetAtlasStyle();
    
    RooStats::UseNLLOffset(true);
    
    // interpret opt
    bool readHistograms  = opt.find("h")!=std::string::npos;
    bool readNtuples     = opt.find("n")!=std::string::npos;
    bool rebinAndSmooth  = opt.find("b")!=std::string::npos; // new: separated from the input creation
    bool createWorkspace = opt.find("w")!=std::string::npos;
    bool doFit           = opt.find("f")!=std::string::npos;
    bool doRanking       = opt.find("r")!=std::string::npos;
    bool doLimit         = opt.find("l")!=std::string::npos;
    bool doSignificance  = opt.find("s")!=std::string::npos;
    bool drawPreFit      = opt.find("d")!=std::string::npos;
    bool drawPostFit     = opt.find("p")!=std::string::npos;
    bool drawSeparation  = opt.find("a")!=std::string::npos;
    
    if(!readNtuples && !rebinAndSmooth){
        TH1::AddDirectory(kFALSE); // FIXME: it would be nice to have a solution which works always
    }
    
    // multi-fit
    bool isMultiFit      = opt.find("m")!=std::string::npos;
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
                if(myMultiFit->fRankingOnly=="all" || myMultiFit->fRankingOnly=="plot")  myMultiFit->PlotNPRankingManager();
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
            if(myMultiFit->fVarNameLH.size()>0) myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName, false );
        }
        //
        if(myMultiFit->fPlotSoverB)    myMultiFit->PlotSummarySoverB();
        //
        return;
    }
    
    // proceed if not multi-fit
    
    TtHFit *myFit = new TtHFit();

    // initialize config reader 
    ConfigReader reader(myFit);

    // read the actual config
    int sc = reader.ReadFullConfig(configFile,options);
    if(sc!=0){
        WriteErrorStatus("myFit::FitExample", "Failed to read the config file.");
        return;
    }
    
    if (TtHFitter::DEBUGLEVEL < 2){
        gErrorIgnoreLevel = kError;
        RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }

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
        myFit->MergeSystematics();
        myFit->CreateCustomAsimov();
        myFit->WriteHistos();
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TtHFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();
    }
    else if(readNtuples){
        myFit->CreateRootFiles();
        myFit->ReadNtuples();
        myFit->Print();
        myFit->CorrectHistograms(); // apply rebinning, smoothing etc...
        myFit->MergeSystematics();
        myFit->CreateCustomAsimov();
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
        myFit->MergeSystematics();
        myFit->CreateCustomAsimov();
        myFit->WriteHistos();
        if(TtHFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TtHFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();
    }

    if(createWorkspace){
        myFit->DrawPruningPlot();
        myFit->SetLumiErr(0.);
        myFit->ToRooStat(true,true);
    }

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

    TthPlot* prefit_plot = 0;
    TthPlot* prefit_plot_valid = 0;
    if( drawPostFit and TtHFitter::PREFITONPOSTFIT ) {
      drawPreFit = true;
      //DrawSummary depends on DrawAndSaveAll to have been executed first. so might as well do the whole drawprefit block to avoid duplication
    }
    
    if(drawPreFit){
        if(TtHFitter::OPTION["PrefitRatioMax"]==2){
            myFit->DrawAndSaveAll("prefit");
            if(myFit->fDoMergedPlot){
                if(myFit->fRegionGroups.size()==0)
                    myFit->DrawMergedPlot("prefit");
                for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
                    myFit->DrawMergedPlot("prefit",myFit->fRegionGroups[i_gr]);
                }
            }
            if(myFit->fDoSummaryPlot){
                prefit_plot       = myFit->DrawSummary("log prefit");
                prefit_plot_valid = myFit->DrawSummary("log valid prefit");
            }
        }
        else{
            myFit->DrawAndSaveAll();
            if(myFit->fDoMergedPlot){
                if(myFit->fRegionGroups.size()==0)
                    myFit->DrawMergedPlot("");
                for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
                    myFit->DrawMergedPlot("",myFit->fRegionGroups[i_gr]);
                }
            }
            if(myFit->fDoSummaryPlot){
                prefit_plot       = myFit->DrawSummary("log");
                prefit_plot_valid = myFit->DrawSummary("log valid");
            }
        }
        if(myFit->fDoTables){
            myFit->BuildYieldTable();
            for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
                myFit->BuildYieldTable("",myFit->fRegionGroups[i_gr]);
            }
            myFit->PrintSystTables();
        }
        int nCols = 2;
        int nRows = 2;
        if(myFit->fNRegions>4){
            nCols = (int)sqrt(myFit->fNRegions);
            if(sqrt(myFit->fNRegions)>nCols) nCols++;
            nRows = (int)sqrt(myFit->fNRegions);
            if(nCols*nRows < myFit->fNRegions) nRows++;
        }
        if(myFit->fDoSignalRegionsPlot) myFit->DrawSignalRegionsPlot(nCols,nRows);
        if(myFit->fDoPieChartPlot)      myFit->DrawPieChartPlot("pre",nCols,nRows);
    }
    
    if(drawPostFit){
        myFit->DrawAndSaveAll("post");
        if(myFit->fDoMergedPlot){
            if(myFit->fRegionGroups.size()==0)
                myFit->DrawMergedPlot("post");
            for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
                myFit->DrawMergedPlot("post",myFit->fRegionGroups[i_gr]);
            }
        }
        if(myFit->fDoSummaryPlot){
            myFit->DrawSummary("log post",      prefit_plot);
            myFit->DrawSummary("log post valid",prefit_plot_valid);
        }
        if(myFit->fDoTables){
            myFit->BuildYieldTable("post");
            for(unsigned int i_gr=0;i_gr<myFit->fRegionGroups.size();i_gr++){
                myFit->BuildYieldTable("post",myFit->fRegionGroups[i_gr]);
            }
            myFit->PrintSystTables("post");
        }
        int nCols = 2;
        int nRows = 2;
        if(myFit->fNRegions>4){
            nCols = (int)sqrt(myFit->fNRegions);
            if(sqrt(myFit->fNRegions)>nCols) nCols++;
            nRows = (int)sqrt(myFit->fNRegions);
            if(nCols*nRows < myFit->fNRegions) nRows++;
        }
        if(myFit->fDoPieChartPlot)      myFit->DrawPieChartPlot("post",nCols,nRows);
    }

    if(drawSeparation){
        myFit->DrawAndSaveSeparationPlots();
    //    myFit->ListOfBestSeparationVariables(); // for the future list of best separation variables
    //    myFit->ListOfBestDataMCVariables();     // for the future list of best data-mc agreement variables based on KS test
    }
    
    if(drawPreFit || drawPostFit || createWorkspace || drawSeparation || rebinAndSmooth) myFit->CloseInputFiles();
  
}

// -------------------------------------------------------
// -------------------------------------------------------
// main function
// -------------------------------------------------------
// -------------------------------------------------------

int main(int argc, char **argv){
  std::string opt="h";
  std::string config="util/myFit,config";
  std::string options="";
  
  if(argc>1) opt     = argv[1];
  if(argc>2) config  = argv[2];
  if(argc>3) options = argv[3];

  // call the function
  FitExample(opt,config,options);
  
  return 0;
}
