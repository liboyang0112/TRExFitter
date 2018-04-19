#include "TtHFitter/Common.h"

#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/ConfigReader.h"
#include "TtHFitter/ConfigReaderMulti.h"
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

// trick to suppress the RooFit banner
// int doBanner(){ return 0; }

void FitExample(std::string opt="h",std::string configFile="util/myFit.config",std::string options=""){
    
    // pre-read the config just to extract the debug level
    std::string debugStr = ReadValueFromConfig(configFile,"DebugLevel");
    if(debugStr=="") WriteWarningStatus("", "Not able to pre-read the DebugLevel => keep the defaul (1).");
    else if(debugStr!="1") TtHFitter::DEBUGLEVEL = atoi(debugStr.c_str());

    // pre-read the logo option
    std::string logoStr = ReadValueFromConfig(configFile,"Logo");
    if(logoStr=="TRUE"){
        std::ifstream logoFile("$TREXFITTER_HOME/logo.txt");
        std::string str;
        std::string logo = "";
        while(getline(logoFile,str)){
            if(!logoFile.good()) break;
            logo+=str;
            logo+="\n";
        }
        std::cout << logo << std::endl;
    }
    
    if(TtHFitter::DEBUGLEVEL<=0)      gErrorIgnoreLevel = kError;
    else if(TtHFitter::DEBUGLEVEL<=1) gErrorIgnoreLevel = kWarning;
    
    // now can set ATLAS style
    if(TtHFitter::DEBUGLEVEL<=0) std::cout.setstate(std::ios_base::failbit);
    SetAtlasStyle();
    if(TtHFitter::DEBUGLEVEL<=0) std::cout.clear();
    
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
    bool groupedImpact   = opt.find("i")!=std::string::npos;
    
    if(!readNtuples && !rebinAndSmooth){
        TH1::AddDirectory(kFALSE); // FIXME: it would be nice to have a solution which works always
    }
    
    // multi-fit
    bool isMultiFit      = opt.find("m")!=std::string::npos;
    if(isMultiFit){
        MultiFit *myMultiFit = new MultiFit();
        ConfigReaderMulti confReaderMulti(myMultiFit);
        int sc = confReaderMulti.ReadFullConfig(configFile,options) ;
    
        if (sc != 0){
            WriteErrorStatus("myFit::FitExample", "Failed to read the config file for multifit.");
            exit(EXIT_FAILURE);
        }

        if(myMultiFit->fCombine){
            if(createWorkspace){
                std::cout << "Combining workspaces..." << std::endl;
                myMultiFit->SaveCombinedWS();
            }
            if(doFit){
                std::cout << "Fitting combining workspace..." << std::endl;
                myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName );
            }
            if(doLimit){
                std::cout << "Getting combined limit..." << std::endl;
                myMultiFit->GetCombinedLimit( myMultiFit->fDataName );
            }
            if(doSignificance){
                std::cout << "Getting combined significance..." << std::endl;
                myMultiFit->GetCombinedSignificance( myMultiFit->fDataName );
            }
            if(doRanking){
                std::cout << "Getting combined ranking..." << std::endl;
                if(myMultiFit->fRankingOnly!="plot")  myMultiFit->ProduceNPRanking( myMultiFit->fRankingOnly );
                if(myMultiFit->fRankingOnly=="all" || myMultiFit->fRankingOnly=="plot")  myMultiFit->PlotNPRankingManager();
            }
            if(groupedImpact){
                std::cout << "Getting combined grouped systematic impact..." << std::endl;
                myMultiFit->fDoGroupedSystImpactTable = true;
                if(myMultiFit->fGroupedImpactCategory!="combine") myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName ); // this calls TtHFit::PerformFit(), which then does the calculation if fDoGroupedSystImpactTable==true
                else                                              myMultiFit->BuildGroupedImpactTable(); // combine the results into one table with option "combine"
            }
        }
        //
        if(myMultiFit->fCompare){
            std::cout << "Comparing fits..." << std::endl;
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
        if(myMultiFit->fPlotSoverB){
            std::cout << "Comparing fits..." << std::endl;
            myMultiFit->PlotSummarySoverB();
        }
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
        exit(EXIT_FAILURE);
    }
    
    WriteInfoStatus("myFit::FitExample", "Successfully read config file.");
    
    if (TtHFitter::DEBUGLEVEL < 2){
        gErrorIgnoreLevel = kError;
        RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }

    // check compatibility between run option and config file
    if(readHistograms && myFit->fInputType!=TtHFit::HIST){
        WriteErrorStatus("myFit::FitExample", "Option \"h\" asked but no HISTO InputType specified in the configuration file. Aborting.");
        return;
    }
    if(readNtuples && myFit->fInputType!=TtHFit::NTUP){
        WriteErrorStatus("myFit::FitExample", "Option \"n\" asked but no NTUP InputType specified in the configuration file. Aborting.");
        return;
    }
      
    // -------------------------------------------------------

    if(readHistograms){
        std::cout << "Reading histograms..." << std::endl;
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
        std::cout << "Reading ntuples..." << std::endl;
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
        std::cout << "Rebinning and smoothing..." << std::endl;
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
        std::cout << "Creating workspace..." << std::endl;
        myFit->DrawPruningPlot();
        myFit->SetLumiErr(0.);
        myFit->ToRooStat(true,true);
    }

    if(doFit){
        std::cout << "Fitting..." << std::endl;
        myFit->Fit();
        myFit->PlotFittedNP();
        myFit->PlotCorrelationMatrix();
    }
    if(doRanking){
        std::cout << "Doing ranking..." << std::endl;
        if(myFit->fRankingOnly!="plot")  myFit->ProduceNPRanking( myFit->fRankingOnly );
        if(myFit->fRankingOnly=="all" || myFit->fRankingOnly=="plot")  myFit->PlotNPRankingManager();
    }
    
    if(doLimit){
        std::cout << "Extracting limit..." << std::endl;
        myFit->GetLimit();
    }
    
    if(doSignificance){
        std::cout << "Extracting significance..." << std::endl;
        myFit->GetSignificance();
    }

    if(groupedImpact){
        std::cout << "Doing grouped systematics impact table..." << std::endl;
        myFit->fDoGroupedSystImpactTable = true;
        if(myFit->fGroupedImpactCategory!="combine") myFit->Fit(); // this calls TtHFit::PerformFit(), which then does the calculation if fDoGroupedSystImpactTable==true
        else                                         myFit->BuildGroupedImpactTable(); // combine the results into one table with option "combine"
    }

    TthPlot* prefit_plot = 0;
    TthPlot* prefit_plot_valid = 0;
    if( drawPostFit and TtHFitter::PREFITONPOSTFIT ) {
      drawPreFit = true;
      //DrawSummary depends on DrawAndSaveAll to have been executed first. so might as well do the whole drawprefit block to avoid duplication
    }
    
    if(drawPreFit){
        std::cout << "Drawing pre-fit plots..." << std::endl;
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
        std::cout << "Drawing post-fit plots..." << std::endl;
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
        std::cout << "Drawing separation plots..." << std::endl;
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
  std::string version = "3.23";
  std::cout << "\033[1mTRExFitter v" << version << " -- Developed by Michele Pinamonti, Loic Valery, Alexander Held, Tomas Dado\033[0m" << std::endl;
  std::cout << "                    No rights reserved, feel free to use and modify it ;)" << std::endl;
  
  std::string opt="h";
  std::string config="util/myFit.config";
  std::string options="";
  
  if(argc>1) opt     = argv[1];
  if(argc>2) config  = argv[2];
  if(argc>3) options = argv[3];

  // call the function
  FitExample(opt,config,options);
  
  return 0;
}
