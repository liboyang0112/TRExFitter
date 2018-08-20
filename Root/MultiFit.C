// Class include
#include "TRExFitter/MultiFit.h"

// Framework onncludes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/TRExFit.h"
#include "TRExFitter/RunSig.h"
#include "TRExFitter/RunAsymptoticsCLs.h"
#include "TRExFitter/RunAsymptoticsCLs_inject.h"

// Roofit includes
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"

// HistFactory includes
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"

// ATLAS stuff
#include "AtlasUtils/AtlasStyle.h"
#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// -------------------------------------------------------------------------------------------------
// class MultiFit

//__________________________________________________________________________________
//
MultiFit::MultiFit(string name){
    fFitList.clear();
    fName = name;
    fDir = "";
    fOutDir = "";
    fLabel = "";
    fShowObserved = false;
    fLimitTitle = "95% CL limit on #sigma/#sigma_{SM}(t#bar{t}H) at m_{H} = 125 GeV";
    fPOITitle = "best fit #mu = #sigma^{t#bar{t}H}/#sigma^{t#bar{t}H}_{SM} for m_{H} = 125 GeV";
    fConfig = new ConfigParser();
    fSaveSuf = "";
    fFitShowObserved.clear();
    fPOI = "";
    fPOIMin = 0;
    fPOIMax = 10;
    fPOIVal = 1;
    fPOIPrecision = "1";
    fLimitMax = 0;
    //
    fCombine       = false;
    fCompare       = false;
    fCompareLimits = true;
    fComparePOI    = true;
    fComparePulls  = true;
    fPlotCombCorrMatrix  = false;
    fStatOnly      = false;
    fIncludeStatOnly = false;
    //
    fDataName      = "obsData";
    fFitType       = 1; // 1: S+B, 2: B-only
    fSignalInjection = false;
    //
    fCombineChByCh = false;
    //
    fNPCategories.clear();
    fNPCategories.push_back("");
    //
    fRndRange = 0.1;
    fRndSeed = -999;
    fUseRnd = false;
    //
    fRankingOnly = "all";
    fGroupedImpactCategory = "all";
    fFastFit = false;
    fFastFitForRanking = true;
    fNuisParListFile = "";
    //
    fPlotSoverB = false;
    fSignalTitle = "t#bar{t}H";
    //
    fFitResultsFile = "";
    fLimitsFile = "";
    fLimitsFiles.clear();
    fBonlySuffix = "";
    fShowSystForPOI = false;
    fVarNameLH.clear();
    fLHscanMin = 999999;
    fLHscanMax = -999999;
    fLHscanSteps = 30;
    //
    fDoGroupedSystImpactTable = false;
    //
    fRunROOTMacros = true; // FIXME: had to temporary set it to true by default, otherwise it crashes...
    fPOIName = "#mu";
    fPOINominal = 1;
    
    fShowTotalOnly = false;
}

//__________________________________________________________________________________
//
MultiFit::~MultiFit(){
    fFitList.clear();
}

//__________________________________________________________________________________
//
void MultiFit::AddFitFromConfig(string configFile,string options,string label,string loadSuf,string wsFile){
    // keep debug level
    int debug = TRExFitter::DEBUGLEVEL;

    fFitList.push_back(new TRExFit());

    // initialize config reader 
    ConfigReader reader(fFitList[fFitList.size()-1]);

    if (reader.ReadFullConfig(configFile,options) != 0){
        WriteErrorStatus("MultiFit::AddFitFromConfig", "Failed to read the config file.");
        exit(EXIT_FAILURE);
}

    fFitLabels.push_back(label);
    fFitSuffs.push_back(loadSuf);
    fWsFiles.push_back(wsFile);
    
    TRExFitter::DEBUGLEVEL = debug;
}

//__________________________________________________________________________________
//
RooWorkspace* MultiFit::CombineWS(){
    WriteInfoStatus("MultiFit::CombineWS", "....................................");
    WriteInfoStatus("MultiFit::CombineWS", "Combining workspaces...");

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);

    std::vector < RooWorkspace* > vec_ws;
    std::vector < std::string > vec_chName;
    RooStats::HistFactory::Measurement *measurement = nullptr;

    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        std::string fitName = fFitList[i_fit]->fInputName;
        std::string fitDir = fFitList[i_fit]->fName;
        WriteDebugStatus("MultiFit::CombineWS", "Adding Fit: " + fitName + ", " + fFitLabels[i_fit] + ", " + fFitSuffs[i_fit] + fitDir);

        RooStats::HistFactory::Measurement *meas;
        std::string fileName = fitDir + "/RooStats/" + fitName + "_combined_" + fitName + fFitSuffs[i_fit] + "_model.root";
        if(fWsFiles[i_fit]!="") fileName = fWsFiles[i_fit];
        WriteDebugStatus("MultiFit::CombineWS", "Opening file " + fileName );
        TFile *rootFile = new TFile(fileName.c_str(),"read");
        RooWorkspace* m_ws = (RooWorkspace*) rootFile->Get("combined");
        WriteDebugStatus("MultiFit::CombineWS", "Getting " + fitName+fFitSuffs[i_fit] );
        meas = (RooStats::HistFactory::Measurement*) rootFile -> Get( (fitName+fFitSuffs[i_fit]).c_str());
        //
        // import measurement if not there yet
        if(!measurement){
            measurement = meas;
        }

        if(!fCombineChByCh){
            //
            // Combine combined workspaces directly
            std::vector<RooStats::HistFactory::Channel> chVec = meas->GetChannels();
            for(unsigned int i_ch=0;i_ch<chVec.size();i_ch++){
                vec_ws.push_back(m_ws);
                vec_chName.push_back(chVec[i_ch].GetName());
            }
        }

        //
        // Alternative way: combine the individual workspaces for the different channels
        // Loop on all the regions in each fit
        if(fCombineChByCh){
            for(unsigned int i_reg=0;i_reg<fFitList[i_fit]->fRegions.size();i_reg++){
                Region *reg = fFitList[i_fit]->fRegions[i_reg];
                if(reg->fRegionType==Region::VALIDATION) continue;
                std::string fileName_tmp = fitDir + "/RooStats/" + fitName + "_" + reg->fName + "_" + fitName + fFitSuffs[i_fit] + "_model.root";
                WriteDebugStatus("MultiFit::CombineWS", "  Opening file " + fileName_tmp );
                TFile *rootFile_tmp = new TFile(fileName_tmp.c_str(),"read");
                RooWorkspace* m_ws_tmp = (RooWorkspace*) rootFile_tmp->Get(reg->fName.c_str());
                WriteDebugStatus("MultiFit::CombineWS", "  Getting " + reg->fName );
                vec_ws.push_back(m_ws_tmp);
                vec_chName.push_back(reg->fName);
            }
        }
        //
    }

    //
    // Create the HistoToWorkspaceFactoryFast object to perform safely the combination
    //
    if(!measurement){
        WriteErrorStatus("MultiFit::CombineWS", "The measurement object has not been retrieved ! Please check.");
        return 0;
    }
    RooStats::HistFactory::HistoToWorkspaceFactoryFast factory(*measurement);

    // Creating the combined model
    RooWorkspace* ws = factory.MakeCombinedModel( vec_chName, vec_ws );

    WriteInfoStatus("MultiFit::CombineWS", "....................................");

    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws, *measurement );

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    return ws;
}

//__________________________________________________________________________________
//
void MultiFit::SaveCombinedWS(){
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    //
    // Creating the rootfile
    //
    TFile *f = new TFile( (fOutDir+"/ws_combined"+fSaveSuf+".root").c_str() , "recreate" );
    //
    // Creating the workspace
    //
    RooWorkspace *ws = CombineWS();
    //
    // Save the workspace
    //
    f->cd();
    ws->Write("combWS");
    f->Close();
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
}

//__________________________________________________________________________________
//
std::map < std::string, double > MultiFit::FitCombinedWS(int fitType, string inputData, bool performFit){
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    TFile *f = new TFile((fOutDir+"/ws_combined"+fSaveSuf+".root").c_str() );
    RooWorkspace *ws = (RooWorkspace*)f->Get("combWS");

    std::map < std::string, double > result;

    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////

    //
    // Fit configuration (1: SPLUSB or 2: BONLY)
    //
    FittingTool *fitTool = new FittingTool();
    fitTool->SetDebug(TRExFitter::DEBUGLEVEL);
    if(fitType==2){
        fitTool -> ValPOI(0.);
        fitTool -> ConstPOI(true);
    } else if(fitType==1){
        fitTool -> ValPOI(fPOIVal);
        fitTool -> ConstPOI(false);
    }
    if(fUseRnd) fitTool -> SetRandomNP(fRndRange, fUseRnd, fRndSeed);

    //
    // Fit starting from custom point
    if(fFitResultsFile!=""){
        fFitList[0]->ReadFitResults(fFitResultsFile);
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(unsigned int i_np=0;i_np<fFitList[0]->fFitResults->fNuisPar.size();i_np++){
            npNames.push_back(  fFitList[0]->fFitResults->fNuisPar[i_np]->fName );
            npValues.push_back( fFitList[0]->fFitResults->fNuisPar[i_np]->fFitValue );
        }
        fitTool -> SetNPs( npNames,npValues );
    }

    std::vector<std::string> vVarNameMinos; vVarNameMinos.clear();
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        for(unsigned int i_minos=0;i_minos<fFitList[i_fit]->fVarNameMinos.size();i_minos++){
            if(FindInStringVector(vVarNameMinos,fFitList[i_fit]->fVarNameMinos[i_minos])<0){
                vVarNameMinos.push_back( fFitList[i_fit]->fVarNameMinos[i_minos] );
            }
        }
    }

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    if(vVarNameMinos.size()>0){
        WriteDebugStatus("MultiFit::FitCombinedWS", "Setting the variables to use MINOS with:");
        for(unsigned int i_minos=0;i_minos<vVarNameMinos.size();i_minos++){
            WriteDebugStatus("MultiFit::FitCombinedWS",  "  " + vVarNameMinos[i_minos]);
        }
        fitTool -> UseMinos(vVarNameMinos);
    }

    //
    // Gets needed objects for the fit
    //
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    //
    // Creates the data object
    //
    RooDataSet* data = 0;
    if(inputData=="asimovData"){
        RooArgSet empty;// = RooArgSet();
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        data = (RooDataSet*)RooStats::AsymptoticCalculator::MakeAsimovData( (*mc), RooArgSet(ws->allVars()), (RooArgSet&)empty);
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    }
    else if(inputData!=""){
        data = (RooDataSet*)ws->data( inputData.c_str() );
    } else {
        WriteWarningStatus("MultiFit::FitCombinedWS", "You didn't specify inputData => will try with observed data !");
        data = (RooDataSet*)ws->data("obsData");
        if(!data){
            WriteWarningStatus("MultiFit::FitCombinedWS", "Observed data not present => will use with asimov data !");
            data = (RooDataSet*)ws->data("asimovData");
        }
    }
    
    if (!data){
        WriteErrorStatus("MultiFit::FitCombinedWS", "Data returns null ptr, probably wrong name in DataName?");
        exit(EXIT_FAILURE);
    }

    // Performs the fit
    gSystem -> mkdir((fOutDir+"/Fits/").c_str(),true);

    if(performFit){
        //
        // Get initial ikelihood value from Asimov
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        double nll0 = 0.;
        if(fGetGoodnessOfFit) nll0 = fitTool -> FitPDF( mc, simPdf, (RooDataSet*)ws->data("asimovData"), false, true );
        
        //
        // Get number of degrees of freedom
        // - number of bins
        int ndof = data->numEntries();
        // - minus number of free & non-constant parameters
        std::vector<std::string> nfList;
        for(auto fit : fFitList){
            for(auto nf : fit->fNormFactors){
                if(nf->fConst) continue;
                if(FindInStringVector(nfList,nf->fName)>=0) continue;
                if(fFitType==2 && fPOI==nf->fName) continue;
                nfList.push_back(nf->fName);
            }
        }
        ndof -= nfList.size();
        
        fitTool -> MinimType("Minuit2");

        // Full fit
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        double nll = fitTool -> FitPDF( mc, simPdf, data, fFastFit );
        fitTool -> ExportFitResultInTextFile(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
        result = fitTool -> ExportFitResultInMap();

        //
        // Goodness of fit
        if(fGetGoodnessOfFit){
            double deltaNLL = nll-nll0;
            double prob = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndof);
            WriteInfoStatus("MultiFit::FitCombinedWS", "----------------------- -------------------------- -----------------------");
            WriteInfoStatus("MultiFit::FitCombinedWS", "----------------------- GOODNESS OF FIT EVALUATION -----------------------");
            WriteInfoStatus("MultiFit::FitCombinedWS", "  NLL0        = " + std::to_string(nll0));
            WriteInfoStatus("MultiFit::FitCombinedWS", "  NLL         = " + std::to_string(nll));
            WriteInfoStatus("MultiFit::FitCombinedWS", "  ndof        = " + std::to_string(ndof));
            WriteInfoStatus("MultiFit::FitCombinedWS", "  dNLL        = " + std::to_string(deltaNLL));
            WriteInfoStatus("MultiFit::FitCombinedWS", "  2dNLL/nof   = " + std::to_string(2.*deltaNLL/ndof));
            WriteInfoStatus("MultiFit::FitCombinedWS", "  probability = " + std::to_string(prob));
            WriteInfoStatus("MultiFit::FitCombinedWS", "----------------------- -------------------------- -----------------------");
            WriteInfoStatus("MultiFit::FitCombinedWS", "----------------------- -------------------------- -----------------------");
        }

        //
        // grouped systematics impact
        if(fDoGroupedSystImpactTable){
            std::string outNameGroupedImpact = fOutDir+"/Fits/GroupedImpact"+fSaveSuf;
            if(fGroupedImpactCategory!="all") outNameGroupedImpact += "_"+fGroupedImpactCategory;
            outNameGroupedImpact += ".txt";
            // need to create a merged list of fSubCategoryImpactMap from the include Fits
            std::map<std::string, std::string> mergedMap;
            mergedMap.clear();
            for(auto fit : fFitList){
                fit->ProduceSystSubCategoryMap();
                for(auto m : fit->fSubCategoryImpactMap){
                    if(mergedMap[m.first]=="") mergedMap[m.first] = m.second;
                    else if(mergedMap[m.first]!=m.second){
                        WriteWarningStatus("MultiFit::FitCombinedWS","Systematics assigned to different SubCategory in the different included Fits. Keeping first Fit convention.");
                    }
                }
            }
            fitTool -> SetSystMap( mergedMap );
            fitTool -> GetGroupedImpact( mc, simPdf, data, ws, fGroupedImpactCategory, outNameGroupedImpact);
        }
    }
    
    //
    // Calls the  function to create LH scan with respect to a parameter
    //
    if(fVarNameLH.size()>0){
        if (fVarNameLH[0]=="all"){
            for(map<string,string>::iterator it=TRExFitter::SYSTMAP.begin(); it!=TRExFitter::SYSTMAP.end(); ++it){
                GetLikelihoodScan( ws, it->first, data, true, fCompare);
            }
        }
        else{
            for(unsigned int i=0; i<fVarNameLH.size(); ++i){
                GetLikelihoodScan( ws, fVarNameLH[i], data, true, fCompare);
            }
        }
    }
    
    if(performFit){
        // Stat-only fit:
        // - read fit resutls
        // - fix all NP to fitted ones before fitting
        if(fIncludeStatOnly){
            WriteInfoStatus("MultiFit::FitCombinedWS", "Fitting stat-only: reading fit results from full fit from file:");
            WriteInfoStatus("MultiFit::FitCombinedWS", "  " + (fOutDir+"/Fits/"+fName+fSaveSuf+".txt"));
            fFitList[0]->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
            std::vector<std::string> npNames;
            std::vector<double> npValues;
            for(unsigned int i_np=0;i_np<fFitList[0]->fFitResults->fNuisPar.size();i_np++){
                bool isNF = false;
                for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
                    if(!fFitList[i_fit]->fFixNPforStatOnlyFit &&
                        FindInStringVector(fFitList[i_fit]->fNormFactorNames,fFitList[0]->fFitResults->fNuisPar[i_np]->fName)>=0){
                        isNF = true;
                        break;
                    }
                }
                if(isNF) continue;
                npNames.push_back(  fFitList[0]->fFitResults->fNuisPar[i_np]->fName );
                npValues.push_back( fFitList[0]->fFitResults->fNuisPar[i_np]->fFitValue );
            }
            fitTool -> FixNPs(npNames,npValues);
            fitTool -> FitPDF( mc, simPdf, data );
            fitTool -> ExportFitResultInTextFile(fOutDir+"/Fits/"+fName+fSaveSuf+"_statOnly.txt");
        }
    }
    
    return result;
}
//__________________________________________________________________________________
//
void MultiFit::GetCombinedLimit(string inputData){
    WriteInfoStatus("MultiFit::GetCombinedLimit", "Runing runAsymptoticsCLs macro...");
    
    string wsFileName = fOutDir+"/ws_combined"+fSaveSuf+".root";
    string cmd;
    if(fSignalInjection){
        if (!fRunROOTMacros){
            LimitsCLs_inject::RunAsymptoticsCLs_inject(wsFileName.c_str(), "combWS", "ModelConfig", inputData.c_str(), "asimovData_0", (fOutDir+"/Limits/").c_str(),(fName+fSaveSuf).c_str(),0.95);
        }
        cmd = "root -l -b -q 'runAsymptoticsCLs_inject.C+(\""+wsFileName+"\",\"combWS\",\"ModelConfig\",\""+inputData+"\",\"asimovData_0\",\""+fOutDir+"/Limits/\",\""+fName+fSaveSuf+"\",0.95)'";
    }
    else{
        if (!fRunROOTMacros){
            LimitsCLs::RunAsymptoticsCLs(wsFileName.c_str(), "combWS", "ModelConfig", inputData.c_str(), "asimovData_0", (fOutDir+"/Limits/").c_str(),(fName+fSaveSuf).c_str(),0.95);
        }
        cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+wsFileName+"\",\"combWS\",\"ModelConfig\",\""+inputData+"\",\"asimovData_0\",\""+fOutDir+"/Limits/\",\""+fName+fSaveSuf+"\",0.95)'";
    }
    //
    // Finally computing the limit
    //
    if (fRunROOTMacros) gSystem->Exec(cmd.c_str());
}
//__________________________________________________________________________________
//
void MultiFit::GetCombinedSignificance(string inputData){
    WriteInfoStatus("MultiFit::GetCombinedSignificance", "Runing runSig macro...");
  
    string wsFileName = fOutDir+"/ws_combined"+fSaveSuf+".root";
    string cmd;
    cmd = "root -l -b -q 'runSig.C(\""+wsFileName+"\",\"combWS\",\"ModelConfig\",\""+inputData+"\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+fSaveSuf+"\",\""+fOutDir+"/Significance\")'";

    //
    // Finally computing the significance
    //
    if (!fRunROOTMacros){
        RunSig(wsFileName.c_str(), "combWS", "ModelConfig", inputData.c_str(), "asimovData_1", "conditionalGlobs_1", "nominalGlobs", (fName+fSaveSuf).c_str(), (fOutDir+"/Significance").c_str());
    }
    else gSystem->Exec(cmd.c_str());
}
//__________________________________________________________________________________
//
void MultiFit::ComparePOI(string POI){
    float xmin = 0;
    float xmax = 2;

    xmax = fPOIMax + (fPOIMax-fPOIMin);
    xmin = fPOIMin;

    string process = fLabel;

    // Fit titles
    vector<string> names;
    vector<string> dirs;
    vector<string> suffs;
    vector<string> titles;
    vector<string> pois;
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        WriteInfoStatus("MultiFit::ComparePOI", "Adding Fit: " + fFitList[i_fit]->fInputName + ", " + fFitLabels[i_fit] + ", " + fFitSuffs[i_fit]);
        names.push_back( fFitList[i_fit]->fInputName );
        dirs.push_back( fFitList[i_fit]->fName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
        if(fFitList[i_fit]->fPOI!="") pois.push_back( fFitList[i_fit]->fPOI );
        else                          pois.push_back( POI );
    }
    if(fCombine){
        WriteInfoStatus("MultiFit::ComparePOI", "Adding Combined Fit");
        names.push_back( fName );
        dirs.push_back( fOutDir );
        titles.push_back( "Combined" );
        suffs.push_back( "" );
        pois.push_back( POI );
    }

    int N = names.size();

    float ymin = -0.5;
    float ymax = N+1-0.5;

    TCanvas *c = new TCanvas("c","c",700,500);
    gStyle->SetEndErrorSize(6.);

    TGraph *g_central    = new TGraph(N);
    TGraphAsymmErrors *g_stat = new TGraphAsymmErrors(N);
    TGraphAsymmErrors *g_tot  = new TGraphAsymmErrors(N);

    int Ndiv = N+1;

    NuisParameter *par;
    bool found = false;

    bool isComb = false;

    // get values
    TRExFit *fit = nullptr;
    for(int i=0;i<N;i++){
        if(fCombine && i==N-1) isComb = true;
        else                   isComb = false;
        //
        if(!isComb) fit = fFitList[i];
        if(!isComb){
            if(fit->fFitResultsFile=="") fit->ReadFitResults(dirs[i]+"/Fits/"+names[i]+suffs[i]+".txt");
            else                         fit->ReadFitResults(fit->fFitResultsFile);
        }
        else{
            if(fFitResultsFile=="") fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
            else                    fit->ReadFitResults(fFitResultsFile);
        }
        found = false;
        for(unsigned int j = 0; j<fit->fFitResults->fNuisPar.size(); ++j){
            par = fit->fFitResults->fNuisPar[j];
            if( pois[i] == par->fName ){
                g_central->SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_stat   ->SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_tot    ->SetPoint(N-i-1,par->fFitValue,N-i-1);
                //
                // temporary put the full uncertainty
                g_stat->SetPointEXhigh(N-i-1,par->fPostFitUp);
                g_stat->SetPointEXlow(N-i-1,-par->fPostFitDown);
                g_stat->SetPointEYhigh(N-i-1,0);
                g_stat->SetPointEYlow(N-i-1,0);
                //
                g_tot->SetPointEXhigh(N-i-1,par->fPostFitUp);
                g_tot->SetPointEXlow(N-i-1,-par->fPostFitDown);
                g_tot->SetPointEYhigh(N-i-1,0);
                g_tot->SetPointEYlow(N-i-1,0);
                //
                found = true;
                break;
            }
        }
        if(!found){
            g_central->SetPoint(N-i-1,-10,N-i-1);
            g_stat->SetPoint(N-i-1,-10,N-i-1);
            g_tot->SetPoint(N-i-1,-10,N-i-1);
            g_stat->SetPointError(N-i-1,0,0,0,0);
            g_tot->SetPointError(N-i-1,0,0,0,0);
        }
    }
    // stat error
    if (!fShowTotalOnly){
        for(int i=0;i<N;i++){
            if(fCombine && i==N-1) isComb = true;
            else                   isComb = false;
            //
            if(!isComb) fit = fFitList[i];
            if(!isComb){
                if(fit->fFitResultsFile=="") fit->ReadFitResults(dirs[i]+"/Fits/"+names[i]+suffs[i]+"_statOnly.txt");
                else                         fit->ReadFitResults(ReplaceString(fit->fFitResultsFile,".txt","_statOnly.txt"));
            }
            else{
                if(fFitResultsFile=="")      fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+"_statOnly.txt");
                else                         fit->ReadFitResults(ReplaceString(fFitResultsFile,".txt","_statOnly.txt"));
            }
            found = false;
            for(unsigned int j = 0; j<fit->fFitResults->fNuisPar.size(); ++j){
                par = fit->fFitResults->fNuisPar[j];
                if( pois[i] == par->fName ){
                    g_stat->SetPointEXhigh(N-i-1,par->fPostFitUp);
                    g_stat->SetPointEXlow(N-i-1,-par->fPostFitDown);
                    g_stat->SetPointEYhigh(N-i-1,0);
                    g_stat->SetPointEYlow(N-i-1,0);
                    found = true;
                    break;
                }
            }
        }
    }

    g_tot->SetLineWidth(3);
    g_stat->SetLineWidth(3);
    g_stat->SetMarkerStyle(kOpenCircle);
    g_tot->SetLineWidth(3);
    g_tot->SetMarkerStyle(kOpenCircle);
    if(TRExFitter::OPTION["FourTopStyle"]){
        g_tot->SetLineColor(kAzure);
        g_tot->SetFillColor(kAzure);
        g_stat->SetLineColor(kCyan);
        g_stat->SetFillColor(kCyan);
        g_stat->SetMarkerStyle(kFullCircle);
        g_stat->SetMarkerColor(kAzure);
        g_stat->SetMarkerSize(1.5);
        //
        for(int i=0;i<N;i++){
            g_tot->SetPointEYlow(i,0.02);
            g_tot->SetPointEYhigh(i,0.02);
            g_stat->SetPointEYlow(i,0.04);
            g_stat->SetPointEYhigh(i,0.04);
        }
    }
    else{
        g_stat->SetLineColor(kGreen-8);
        g_stat->SetMarkerSize(0);
        g_tot->SetLineColor(kBlack);
    }
    if(TRExFitter::OPTION["FourTopStyle"]){
        g_central->SetMarkerColor(kWhite);
        g_central->SetMarkerStyle(kFullCircle);
        g_central->SetMarkerSize(1.);
    }
    else{
        g_central->SetMarkerStyle(kFullCircle);
        g_central->SetMarkerColor(kRed);
        g_central->SetMarkerSize(1.5);
    }
    g_tot->SetMarkerSize(0);

    TH1D* h_dummy = new TH1D("h_dummy","h_dummy",1,xmin,xmax);
    h_dummy->Draw();
    h_dummy->SetMinimum(ymin);
    h_dummy->SetMaximum(ymax);
    h_dummy->SetLineColor(kWhite);
    h_dummy->GetYaxis()->Set(N+1,ymin,ymax);
    h_dummy->GetYaxis()->SetNdivisions(Ndiv);

    TLatex *tex = new TLatex();

    for(int i=0;i<N;i++){
        h_dummy->GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
        if(fShowSystForPOI){
            tex->SetTextSize(gStyle->GetTextSize()*1.2);
            tex->DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form(("#font[62]{%." + fPOIPrecision + "f}").c_str(),g_central->GetX()[N-i-1]));
            tex->DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{^{#plus%." + fPOIPrecision + "f}}").c_str(),g_tot->GetErrorXhigh(N-i-1)));
            tex->DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{_{#minus%." + fPOIPrecision + "f}}").c_str(),g_tot->GetErrorXlow(N-i-1)));
            tex->DrawLatex(xmin+0.69*(xmax-xmin),N-i-1,"(");
            if (!fShowTotalOnly){
	            tex->DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + fPOIPrecision + "f}}").c_str(),g_stat->GetErrorXhigh(N-i-1)));
	            tex->DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + fPOIPrecision + "f}}").c_str(),g_stat->GetErrorXlow(N-i-1)));
                tex->DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + fPOIPrecision + "f}}").c_str(),
                    sqrt( pow(g_tot->GetErrorXhigh(N-i-1),2) - pow(g_stat->GetErrorXhigh(N-i-1),2) ) ) );
                tex->DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + fPOIPrecision + "f}}").c_str(),
                    sqrt( pow(g_tot->GetErrorXlow(N-i-1) ,2) - pow(g_stat->GetErrorXlow(N-i-1) ,2) ) ) );
                tex->DrawLatex(xmin+0.94*(xmax-xmin),N-i-1,")");
            }
        }
        else{
            tex->DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form((fPOIName+" = %." + fPOIPrecision + "f").c_str(),g_central->GetX()[N-i-1]));
            tex->DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("^{#plus%." + fPOIPrecision + "f}").c_str(),g_tot->GetErrorXhigh(N-i-1)));
            tex->DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("_{#minus%." + fPOIPrecision + "f}").c_str(),g_tot->GetErrorXlow(N-i-1)));
            if (!fShowTotalOnly){
                tex->DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("^{#plus%." + fPOIPrecision + "f}").c_str(),g_stat->GetErrorXhigh(N-i-1)));
                tex->DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("_{#minus%." + fPOIPrecision + "f}").c_str(),g_stat->GetErrorXlow(N-i-1)));
            }
        }
    }

    TLine *l_SM = new TLine(fPOINominal,-0.5,fPOINominal,N-0.5);
    l_SM->SetLineWidth(2);
    l_SM->SetLineColor(kGray);
    l_SM->Draw("same");

    if(fCombine){
        TLine *l_h = new TLine(xmin,0.5,xmax,0.5);
        l_h->SetLineWidth(2);
        l_h->SetLineColor(kBlack);
        l_h->SetLineStyle(kDashed);
        l_h->Draw("same");
    }

    if(TRExFitter::OPTION["FourTopStyle"]){
        g_tot->Draw("E2 same");
        if (!fShowTotalOnly){
            g_stat->Draw("PE2 same");
        }
        g_central->Draw("P same");
    }
    else{
        g_tot->Draw("E same");
        if (!fShowTotalOnly){
            g_stat->Draw("E same");
        }
        g_central->Draw("P same");
    }

    gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
    gPad->SetBottomMargin( 1.15*gPad->GetBottomMargin() );
    gPad->SetTopMargin( 1.8*gPad->GetTopMargin() );
    h_dummy->GetXaxis()->SetTitle(fPOITitle.c_str());
    h_dummy->GetYaxis()->SetTickSize(0);

    c->RedrawAxis();

    ATLASLabel(0.32,0.93,fFitList[0]->fAtlasLabel.c_str(),kBlack);
    myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(process!="") myText(0.94,0.85,kBlack,Form("#kern[-1]{%s}",process.c_str()));

    TLegend *leg;
    leg = new TLegend(0.35,0.775,0.7,0.9);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(g_tot,"tot.","l");
    if (!fShowTotalOnly){
        leg->AddEntry(g_stat,"stat.","l");
    }
    leg->Draw();

    if(fShowSystForPOI){
        tex->DrawLatex(xmin+0.6*(xmax-xmin),N-0.4,"#font[62]{tot}");
        tex->DrawLatex(xmin+0.69*(xmax-xmin),N-0.4,"(");
        if (!fShowTotalOnly){
            tex->DrawLatex(xmin+0.72*(xmax-xmin),N-0.4,"stat");
        }
        tex->DrawLatex(xmin+0.83*(xmax-xmin),N-0.4,"syst");
        tex->DrawLatex(xmin+0.94*(xmax-xmin),N-0.4,")");
    }
    else{
        tex->DrawLatex(xmin+(0.7-0.02)*(xmax-xmin),N-0.4,"( tot )");
        if (!fShowTotalOnly){
            tex->DrawLatex(xmin+(0.85-0.02)*(xmax-xmin),N-0.4,"( stat )");
        }
    }

    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
        c->SaveAs( (fOutDir+"/POI"+fSaveSuf+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::CompareLimit(){
    float xmax = 2;
    string process = fLabel;
    gStyle->SetEndErrorSize(0.);

    // ---

    // Fit titles
    vector<string> dirs;   dirs.clear();
    vector<string> names;  names.clear();
    vector<string> suffs;  suffs.clear();
    vector<string> titles; titles.clear();
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        WriteInfoStatus("MultiFit::CompareLimit", "Adding Fit: " + fFitList[i_fit]->fInputName + ", " + fFitLabels[i_fit] + ", " + fFitSuffs[i_fit]);
        dirs.push_back( fFitList[i_fit]->fName );
        names.push_back( fFitList[i_fit]->fInputName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
    }
    if(fCombine){
        WriteInfoStatus("MultiFit::CompareLimit", "Adding combined limit");
        dirs.push_back( fOutDir );
        names.push_back( fName );
        titles.push_back( "Combined" );
        suffs.push_back( fSaveSuf );
        if(fShowObserved) fFitShowObserved.push_back(true);
    }

    // ---

    bool showObs = fShowObserved;

    unsigned int N = names.size();

    float ymin = -0.5;
    float ymax = N-0.5;

    TCanvas *c = new TCanvas("c","c",700,500);

    TGraphErrors *g_obs = new TGraphErrors(N);
    TGraphErrors *g_exp = new TGraphErrors(N);
    TGraphErrors *g_inj = new TGraphErrors(N);
    TGraphAsymmErrors *g_1s = new TGraphAsymmErrors(N);
    TGraphAsymmErrors *g_2s = new TGraphAsymmErrors(N);

    int Ndiv = N+1;

    TFile *f;
    TH1* h;
    TH1* h_old = nullptr;

    // get values
    for(unsigned int i=0;i<N;i++){
        if(i>fLimitsFiles.size()-1){
            if(fLimitsFile!="") fLimitsFiles.push_back(fLimitsFile);
            else                fLimitsFiles.push_back("");
        }
        if(fLimitsFiles[i]==""){
            if(fSignalInjection){
                f = new TFile(Form("%s/Limits/%s_injection.root",dirs[i].c_str(),(names[i]+suffs[i]).c_str()) );
                WriteInfoStatus("MultiFit::CompareLimit", "Reading file " + dirs[i] + "/Limits/" + (names[i]+suffs[i]) + "_injection.root");
            }
            else{
                f = new TFile(Form("%s/Limits/%s.root",dirs[i].c_str(),(names[i]+suffs[i]).c_str()) );
                WriteInfoStatus("MultiFit::CompareLimit", "Reading file " + dirs[i] + "/Limits/" + (names[i]+suffs[i]) + ".root");
            }
        }
        else{
            f = new TFile(fLimitsFiles[i].c_str());
            WriteInfoStatus("MultiFit::CompareLimit", "Reading file " + fLimitsFiles[i]);
        }
        h = (TH1*)f->Get("limit");
        if(fSignalInjection) h_old = (TH1*)f->Get("limit_old");

        WriteDebugStatus("MultiFit::CompareLimit", "bin 1 content: " + std::to_string(h->GetBinContent(1)));
        if(fFitShowObserved[i]) g_obs->SetPoint(N-i-1,h->GetBinContent(1),N-i-1);
        else g_obs->SetPoint(N-i-1,-1,N-i-1);
        g_exp->SetPoint(N-i-1,h->GetBinContent(2),N-i-1);
        if(fSignalInjection) g_inj->SetPoint(N-i-1,h_old->GetBinContent(7),N-i-1);
        g_1s->SetPoint(N-i-1,h->GetBinContent(2),N-i-1);
        g_2s->SetPoint(N-i-1,h->GetBinContent(2),N-i-1);
        g_obs->SetPointError(N-i-1,0,0.5);
        g_exp->SetPointError(N-i-1,0,0.5);
        g_inj->SetPointError(N-i-1,0,0.5);
        g_1s->SetPointError(N-i-1,h->GetBinContent(2)-h->GetBinContent(5),h->GetBinContent(4)-h->GetBinContent(2),0.5,0.5);
        g_2s->SetPointError(N-i-1,h->GetBinContent(2)-h->GetBinContent(6),h->GetBinContent(3)-h->GetBinContent(2),0.5,0.5);

        if(h->GetBinContent(1)>xmax) xmax = h->GetBinContent(1);
        if(h->GetBinContent(2)>xmax) xmax = h->GetBinContent(2);
        if(h->GetBinContent(3)>xmax) xmax = h->GetBinContent(3);
        if(h->GetBinContent(4)>xmax) xmax = h->GetBinContent(4);
        if(h->GetBinContent(5)>xmax) xmax = h->GetBinContent(5);
        if(h->GetBinContent(6)>xmax) xmax = h->GetBinContent(6);
    }

    g_obs->SetLineWidth(3);
    g_exp->SetLineWidth(3);
    g_exp->SetLineStyle(2);
    g_inj->SetLineWidth(3);
    g_inj->SetLineStyle(2);
    g_inj->SetLineColor(kRed);
    g_1s->SetFillColor(kGreen);
    g_1s->SetLineWidth(3);
    g_1s->SetLineStyle(2);
    g_2s->SetFillColor(kYellow);
    g_2s->SetLineWidth(3);
    g_2s->SetLineStyle(2);

    g_2s->SetMarkerSize(0);
    g_1s->SetMarkerSize(0);
    g_exp->SetMarkerSize(0);
    g_obs->SetMarkerSize(0);
    g_inj->SetMarkerSize(0);

    if(fLimitMax!=0) xmax = fLimitMax;

    TH1D* h_dummy = new TH1D("h_dummy","h_dummy",1,0,xmax);
    h_dummy->Draw();
    h_dummy->SetMinimum(ymin);
    h_dummy->SetMaximum(ymax);
    h_dummy->SetLineColor(kWhite);
    h_dummy->GetYaxis()->Set(N,ymin,ymax);
    h_dummy->GetYaxis()->SetNdivisions(Ndiv);
    for(unsigned int i=0;i<N;i++){
        h_dummy->GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
    }

    g_2s->Draw("E2 same");
    g_1s->Draw("E2 same");
    g_exp->Draw("E same");
    if(showObs) g_obs->Draw("E same");
    if(fSignalInjection) g_inj->Draw("E same");

    TLine *l_SM = new TLine(fPOINominal,-0.5,fPOINominal,N-0.5);
    l_SM->SetLineWidth(2);
    l_SM->SetLineColor(kGray);
    l_SM->Draw("same");

    c->RedrawAxis();

    gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
    gPad->SetBottomMargin( 1.15*gPad->GetBottomMargin() );
    gPad->SetTopMargin( 1.8*gPad->GetTopMargin() );
    h_dummy->GetXaxis()->SetTitle(fLimitTitle.c_str());

    ATLASLabel(0.32,0.93,fFitList[0]->fAtlasLabel.c_str(),kBlack);
    myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(process!="") myText(0.94,0.85,kBlack,Form("#kern[-1]{%s}",process.c_str()));

    TLegend *leg;
    if(showObs) leg = new TLegend(0.65,0.2,0.95,0.40);
    else        leg = new TLegend(0.65,0.2,0.95,0.35);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(g_1s,"Expected #pm 1#sigma","lf");
    leg->AddEntry(g_2s,"Expected #pm 2#sigma","lf");
    if(showObs) leg->AddEntry(g_obs,"Observed","l");
    if(fSignalInjection) leg->AddEntry(g_inj,("Expected ("+fPOIName+"=1)").c_str(),"l");
    leg->Draw();

    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
        c->SaveAs( (fOutDir+"/Limits" + fSaveSuf +  + "."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::ComparePulls(string category){
    float ydist = 0.2;

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    vector<float>  yshift;
    int color[] = {kBlack,kRed,kBlue,kViolet};

    int style[] = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown};

    unsigned int N = fFitList.size();
    if(fCombine) N++;

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            WriteInfoStatus("MultiFit::ComparePulls", "Adding Combined Fit");
            dirs.push_back( fOutDir );
            names.push_back( fName );
            titles.push_back( "Combined" );
            suffs.push_back( "" );
        }
        else{
            dirs.push_back( fFitList[i_fit]->fName );
            names.push_back( fFitList[i_fit]->fInputName );
            titles.push_back( fFitLabels[i_fit] );
            suffs.push_back( fFitSuffs[i_fit] );
        }
        yshift.push_back( 0. - ydist*N/2. + ydist*i_fit );
    }

    float xmin = -2.9;
    float xmax = 2.9;
    float max = 0;
    string npToExclude[] = {"gamma_","stat_"};
    bool brazilian = true;

    // create a list of Systematics
    std::vector< string > Names;  Names.clear();
    std::vector< string > Titles; Titles.clear();
    std::vector< string > Categories; Categories.clear();
    string systName;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        for(int i_syst=0;i_syst<fFitList[i_fit]->fNSyst;i_syst++){
            systName = fFitList[i_fit]->fSystematics[i_syst]->fNuisanceParameter;
            if(systName == "") systName = fFitList[i_fit]->fSystematics[i_syst]->fName;
            if(FindInStringVector(Names,systName)<0){
                Names.push_back(systName);
                Titles.push_back(fFitList[i_fit]->fSystematics[i_syst]->fTitle);
                if(FindInStringVector(fFitList[i_fit]->fDecorrSysts,fFitList[i_fit]->fSystematics[i_syst]->fName)>=0){
                    Titles[Titles.size()-1] += fFitList[i_fit]->fDecorrSuff;
                }
                Categories.push_back(fFitList[i_fit]->fSystematics[i_syst]->fCategory);
            }
        }
    }
    unsigned int Nsyst = Names.size();

    // read fit resutls
    NuisParameter *par;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        if(fFitList[i_fit]->fFitResultsFile!="") fFitList[i_fit]->ReadFitResults(fFitList[i_fit]->fFitResultsFile);
        else                                     fFitList[i_fit]->ReadFitResults(dirs[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".txt");
    }

    // exclude unused systematics
    std::vector<string> NamesNew; NamesNew.clear();
    std::vector<string> TitlesNew; TitlesNew.clear();
    std::vector<string> CategoriesNew; CategoriesNew.clear();
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        FitResults *fitRes;
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            fitRes = fFitList[i_fit]->fFitResults;
            for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
                par = fitRes->fNuisPar[j];
                systName = par->fName;
                if(systName==Names[i_syst]){
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
        if(found){
            if(category=="" || category==Categories[i_syst]){
                NamesNew.push_back(Names[i_syst]);
                TitlesNew.push_back(Titles[i_syst]);
                CategoriesNew.push_back(Categories[i_syst]);
            }
        }
    }
    //
    Nsyst = NamesNew.size();
    Names.clear();
    Titles.clear();
    Categories.clear();
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        Names.push_back(NamesNew[i_syst]);
        Titles.push_back(TitlesNew[i_syst]);
        CategoriesNew.push_back(Categories[i_syst]);
    }
    if(fNuisParListFile!=""){
        //
        // reorder NPs
        Names.clear();
        Titles.clear();
        Categories.clear();
        ifstream in;
        in.open(fNuisParListFile.c_str());
        while(true){
            in >> systName;
            if(!in.good()) break;
            WriteDebugStatus("MultiFit::ComparePulls", "Looking for " + systName);
            std::string temp_string = "";
            for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
                if(NamesNew[i_syst]==systName){
                    temp_string+=  "found";
                    temp_string += ", title = " + TitlesNew[i_syst];
                    Names.push_back(NamesNew[i_syst]);
                    Titles.push_back(TitlesNew[i_syst]);
                    Categories.push_back(CategoriesNew[i_syst]);
                    break;
                }
            }
            WriteDebugStatus("MultiFit::ComparePulls", temp_string);
        }
        in.close();
    }
    Nsyst = Names.size();

    // fill stuff
    std::vector< TGraphAsymmErrors* > g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,float> centralMap;
        std::map<string,float> errUpMap;
        std::map<string,float> errDownMap;
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            if(fFitResultsFile!="") fitRes->ReadFromTXT(fFitResultsFile);
            else                    fitRes->ReadFromTXT(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults;
        }
        for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
            par = fitRes->fNuisPar[j];
            systName = par->fName;
            centralMap[systName] = par->fFitValue;
            errUpMap[systName]   = par->fPostFitUp;
            errDownMap[systName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.push_back( new TGraphAsymmErrors(Nsyst) );
        for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
            systName = Names[i_syst];
            if(centralMap[systName]!=0 || (errUpMap[systName]!=0 || errDownMap[systName]!=0)){
                g[i_fit]->SetPoint(i_syst,centralMap[systName],(Nsyst-i_syst-1)+0.5+yshift[i_fit]);
                g[i_fit]->SetPointEXhigh(i_syst,  errUpMap[systName]);
                g[i_fit]->SetPointEXlow( i_syst, -errDownMap[systName]);
            }
            else{
                g[i_fit]->SetPoint(i_syst,-10,-10);
                g[i_fit]->SetPointEXhigh(i_syst, 0);
                g[i_fit]->SetPointEXlow( i_syst, 0);
            }
        }
    }

    max = Nsyst;

    int lineHeight = 20;
    int offsetUp = 50;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas *c = new TCanvas("c","c",800,newHeight);
    c->SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D *h_dummy = new TH1D("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy->SetMaximum(max);
    h_dummy->SetLineWidth(0);
    h_dummy->SetFillStyle(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->SetFillColor(kWhite);
    h_dummy->SetMinimum(0.);
    h_dummy->GetYaxis()->SetLabelSize(0);
    h_dummy->Draw();
    h_dummy->GetYaxis()->SetNdivisions(0);

    TLine l0;
    TBox b1, b2;
    if(brazilian){
        l0 = TLine(0,0,0,max);
        l0.SetLineStyle(7);
        l0.SetLineColor(kBlack);
        b1 = TBox(-1,0,1,max);
        b2 = TBox(-2,0,2,max);
        b1.SetFillColor(kGreen);
        b2.SetFillColor(kYellow);
        b2.Draw("same");
        b1.Draw("same");
        l0.Draw("same");
    }

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        g[i_fit]->SetLineColor(color[i_fit]);
        g[i_fit]->SetMarkerColor(color[i_fit]);
        g[i_fit]->SetMarkerStyle(style[i_fit]);
        g[i_fit]->Draw("P same");
    }

    TLatex *systs = new TLatex();
    systs->SetTextSize( systs->GetTextSize()*0.8 );
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        systs->DrawLatex(3.,(Nsyst-i_syst-1)+0.25,Titles[i_syst].c_str());
    }
    h_dummy->GetXaxis()->SetLabelSize( h_dummy->GetXaxis()->GetLabelSize()*0.9 );

    TLegend *leg;
    leg = new TLegend(0.01,1.-0.03*(30./max),0.75,0.99);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        leg->AddEntry(g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg->Draw();

    gPad->RedrawAxis();

    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
        if(category=="") c->SaveAs((fOutDir+"/NuisPar_comp"+fSaveSuf+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
        else             c->SaveAs((fOutDir+"/NuisPar_comp"+fSaveSuf+"_"+category+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
    }
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::CompareNormFactors(string category){
    float ydist = 0.2;

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    vector<float>  yshift;

    int color[] = {kBlack,kRed,kBlue,kViolet};
    int style[] = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown};

    unsigned int N = fFitList.size();
    if(fCombine) N++;

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            WriteInfoStatus("MultiFit::CompareNormFactors", "Adding Combined Fit");
            dirs.push_back( fOutDir );
            names.push_back( fName );
            titles.push_back( "Combined" );
            suffs.push_back( "" );
        }
        else{
            dirs.push_back( fFitList[i_fit]->fName );
            names.push_back( fFitList[i_fit]->fInputName );
            titles.push_back( fFitLabels[i_fit] );
            suffs.push_back( fFitSuffs[i_fit] );
        }
        yshift.push_back( 0. - ydist*N/2. + ydist*(N-i_fit-1) );
    }

    float xmin = -1;
    float xmax = 10;
    float max = 0;

    // create a list of Norm Factors
    std::vector< string > Names;  Names.clear();
    std::vector< string > Titles; Titles.clear();
    std::vector< string > Categories; Categories.clear();
    string normName;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        for(int i_norm=0;i_norm<fFitList[i_fit]->fNNorm;i_norm++){
            normName = fFitList[i_fit]->fNormFactors[i_norm]->fName;
            if(normName==fPOI) continue;
            if(FindInStringVector(Names,normName)<0){
                Names.push_back(normName);
                Titles.push_back(fFitList[i_fit]->fNormFactors[i_norm]->fTitle);
                Categories.push_back(fFitList[i_fit]->fNormFactors[i_norm]->fCategory);
            }
        }
    }
    unsigned int Nnorm = Names.size();

    // read fit resutls
    NuisParameter *par;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        fFitList[i_fit]->ReadFitResults(dirs[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".txt");
    }

    // exclude norm factors
    std::vector<string> NamesNew; NamesNew.clear();
    std::vector<string> TitlesNew; TitlesNew.clear();
    std::vector<string> CategoriesNew; CategoriesNew.clear();
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        FitResults *fitRes;
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            fitRes = fFitList[i_fit]->fFitResults;
            for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
                par = fitRes->fNuisPar[j];
                normName = par->fName;
                if(normName==Names[i_norm]){
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
        if(found){
            if(category=="" || category==Categories[i_norm]){
                NamesNew.push_back(Names[i_norm]);
                TitlesNew.push_back(Titles[i_norm]);
                CategoriesNew.push_back(Categories[i_norm]);
            }
        }
    }
    //
    Nnorm = NamesNew.size();
    Names.clear();
    Titles.clear();
    Categories.clear();
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        Names.push_back(NamesNew[i_norm]);
        Titles.push_back(TitlesNew[i_norm]);
        CategoriesNew.push_back(Categories[i_norm]);
    }
    if(fNuisParListFile!=""){
        //
        // reorder NPs
        Names.clear();
        Titles.clear();
        Categories.clear();
        ifstream in;
        in.open(fNuisParListFile.c_str());
        while(true){
            in >> normName;
            if(!in.good()) break;
            WriteDebugStatus("MultiFit::CompareNormFactors", "Looking for " + normName + "... ");
            std::string temp_string = "";
            for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
                if(NamesNew[i_norm]==normName){
                    temp_string+= "found";
                    temp_string+= ", title = " + TitlesNew[i_norm];
                    Names.push_back(NamesNew[i_norm]);
                    Titles.push_back(TitlesNew[i_norm]);
                    Categories.push_back(CategoriesNew[i_norm]);
                    break;
                }
            }
            WriteDebugStatus("MultiFit::CompareNormFactors", temp_string);
        }
        in.close();
    }
    Nnorm = Names.size();

    // fill stuff
    std::vector< TGraphAsymmErrors* > g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,float> centralMap;
        std::map<string,float> errUpMap;
        std::map<string,float> errDownMap;
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            if(fFitResultsFile!="") fitRes->ReadFromTXT(fFitResultsFile);
            else                    fitRes->ReadFromTXT(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults;
        }
        for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
            par = fitRes->fNuisPar[j];
            normName = par->fName;
            centralMap[normName] = par->fFitValue;
            errUpMap[normName]   = par->fPostFitUp;
            errDownMap[normName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.push_back( new TGraphAsymmErrors(Nnorm) );
        for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
            normName = Names[i_norm];
            if(centralMap[normName]!=0 || (errUpMap[normName]!=0 || errDownMap[normName]!=0)){
                g[i_fit]->SetPoint(i_norm,centralMap[normName],(Nnorm-i_norm-1)+0.5+yshift[i_fit]);
                g[i_fit]->SetPointEXhigh(i_norm,  errUpMap[normName]);
                g[i_fit]->SetPointEXlow( i_norm, -errDownMap[normName]);
            }
            else{
                g[i_fit]->SetPoint(i_norm,-10,-10);
                g[i_fit]->SetPointEXhigh(i_norm, 0);
                g[i_fit]->SetPointEXlow( i_norm, 0);
            }
        }
    }

    max = Nnorm;

    int lineHeight = 50;
    int offsetUp = 50;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas *c = new TCanvas("c","c",800,newHeight);
    c->SetTicks(1,0);
    gPad->SetLeftMargin(0.2/(8./6.));
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D *h_dummy = new TH1D("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy->SetMaximum(max);
    h_dummy->SetLineWidth(0);
    h_dummy->SetFillStyle(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->SetFillColor(kWhite);
    h_dummy->SetMinimum(0.);
    h_dummy->GetYaxis()->SetLabelSize(0);
    h_dummy->Draw();
    h_dummy->GetYaxis()->SetNdivisions(0);

    TLine l1;
    l1 = TLine(1,0,1,max);
    l1.SetLineColor(kGray);
    l1.Draw("same");

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        g[i_fit]->SetLineColor(color[i_fit]);
        g[i_fit]->SetMarkerColor(color[i_fit]);
        g[i_fit]->SetMarkerStyle(style[i_fit]);
        g[i_fit]->Draw("P same");
    }

    TLatex *norms = new TLatex();
    norms->SetTextSize( norms->GetTextSize()*0.8 );
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        norms->DrawLatex(xmin-0.15*(xmax-xmin),(Nnorm-i_norm-1)+0.25,Titles[i_norm].c_str());
    }
    TLatex *values = new TLatex();
    values->SetTextSize( values->GetTextSize()*0.8 );
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            values->DrawLatex((xmin+(xmax-xmin)*(0.45+i_fit*0.55/N)),(Nnorm-i_norm-1)+0.25,
                              Form("%.2f^{+%.2f}_{-%.2f}",g[i_fit]->GetX()[i_norm],g[i_fit]->GetErrorXhigh(i_norm),g[i_fit]->GetErrorXlow(i_norm)));
        }
    }
    h_dummy->GetXaxis()->SetLabelSize( h_dummy->GetXaxis()->GetLabelSize()*0.9 );

    TLegend *leg;
    leg = new TLegend(0.45,1.-0.02*(30./max),0.98,0.99);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        leg->AddEntry(g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg->Draw();

    gPad->RedrawAxis();

    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
        if(category=="") c->SaveAs((fOutDir+"/NormFactors_comp"+fSaveSuf+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
        else             c->SaveAs((fOutDir+"/NormFactors_comp"+fSaveSuf+"_"+category+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
    }
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::PlotCombinedCorrelationMatrix(){
    TRExFit *fit = fFitList[0];
    if(fit->fStatOnly){
        WriteInfoStatus("MultiFit::PlotCombinedCorrelationMatrix", "Stat only fit => No Correlation Matrix generated.");
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TRExFitter::CORRELATIONTHRESHOLD)
    fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
    if(fit->fFitResults){
        for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
            fit->fFitResults->DrawCorrelationMatrix(fOutDir+"/CorrMatrix_comb"+fSaveSuf+"."+TRExFitter::IMAGEFORMAT[i_format],TRExFitter::CORRELATIONTHRESHOLD);
    }
}

//____________________________________________________________________________________
//
void MultiFit::ProduceNPRanking( string NPnames/*="all"*/ ){
    WriteInfoStatus("MultiFit::ProduceNPRanking", "....................................");
    WriteInfoStatus("MultiFit::ProduceNPRanking", "Producing Ranking...");

    if(fFitType==2){
        WriteErrorStatus("MultiFit::ProduceNPRanking", "For ranking plots, the SPLUSB FitType is needed.");
        abort();
    }

    string inputData = fDataName;
    unsigned int N = fFitList.size();

    // create a list of Systematics
    std::vector< Systematic* > vSystematics;  vSystematics.clear();
    std::vector< std::string > Names;  Names.clear();
    string systName;
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        for(int i_syst=0;i_syst<fFitList[i_fit]->fNSyst;i_syst++){
            systName = fFitList[i_fit]->fSystematics[i_syst]->fNuisanceParameter;
            if(FindInStringVector(Names,systName)<0){
                Names.push_back(systName);
                vSystematics.push_back(fFitList[i_fit]->fSystematics[i_syst]);
            }
        }
    }
    unsigned int Nsyst = Names.size();

    // create a list of norm factors
    std::vector< NormFactor* > vNormFactors;  vNormFactors.clear();
    std::vector< std::string > nfNames;  nfNames.clear();
    string normName;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        for(int i_norm=0;i_norm<fFitList[i_fit]->fNNorm;i_norm++){
            normName = fFitList[i_fit]->fNormFactors[i_norm]->fName;
            if(FindInStringVector(nfNames,systName)<0){
                nfNames.push_back(normName);
                vNormFactors.push_back(fFitList[i_fit]->fNormFactors[i_norm]);
            }
        }
    }
    unsigned int Nnorm = nfNames.size();

    //
    // List of systematics to check
    //
    std::vector< string > nuisPars;
    std::vector< bool > isNF;
    std::vector<string> systNames_unique;
    for(int i_syst=0;i_syst< (int) Nsyst;i_syst++){
        if(NPnames=="all" || NPnames==vSystematics[i_syst]->fNuisanceParameter ||
            ( atoi(NPnames.c_str())==i_syst && (atoi(NPnames.c_str())>0 || strcmp( NPnames.c_str(), "0") ==0 ) )
            ){
            if(vSystematics[i_syst]->fType == Systematic::SHAPE) continue;
            if (std::find(systNames_unique.begin(), systNames_unique.end(),
                vSystematics[i_syst]->fNuisanceParameter) == systNames_unique.end()){
                systNames_unique.push_back(vSystematics[i_syst]->fNuisanceParameter);
            }
            else {
                continue;
            }
            nuisPars.push_back( vSystematics[i_syst]->fNuisanceParameter );
            isNF.push_back( false );
        }
    }
    for(int i_norm=0;i_norm<(int)Nnorm;i_norm++){
        if(fPOI==vNormFactors[i_norm]->fName) continue;
        if(NPnames=="all" || NPnames==vNormFactors[i_norm]->fName ||
            ( ((atoi(NPnames.c_str())-(int)Nnorm) == i_norm) && (atoi(NPnames.c_str())>0 || strcmp( NPnames.c_str(), "0")==0) )
            ){
            nuisPars.push_back( vNormFactors[i_norm]->fName );
            isNF.push_back( true );
        }
    }

    //
    // Text files containing information necessary for drawing of ranking plot
    //
    string outName = fOutDir+"/Fits/NPRanking";
    if(NPnames!="all") outName += "_"+NPnames;
    outName += ".txt";
    ofstream outName_file(outName.c_str());
    //
    double central;
    double up;
    double down;
    double muhat;
    std::map< string,double > muVarUp;
    std::map< string,double > muVarDown;
    std::map< string,double > muVarNomUp;
    std::map< string,double > muVarNomDown;

    //
    // Get the combined model
    //
    TFile *f = new TFile((fOutDir+"/ws_combined"+fSaveSuf+".root").c_str() );
    RooWorkspace *ws = (RooWorkspace*)f->Get("combWS");

    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());

    //
    // Creates the data object
    //
    RooDataSet* data = 0;
    if(inputData=="asimovData"){
        RooArgSet empty;// = RooArgSet();
        data = (RooDataSet*)RooStats::AsymptoticCalculator::MakeAsimovData( (*mc), RooArgSet(ws->allVars()), (RooArgSet&)empty);
    }
    else if(inputData!=""){
        data = (RooDataSet*)ws->data( inputData.c_str() );
    } else {
        WriteWarningStatus("MultiFit::ProduceNPRanking", "You didn't specify inputData => will try with observed data !");
        data = (RooDataSet*)ws->data("obsData");
        if(!data){
            WriteWarningStatus("MultiFit::ProduceNPRanking", "Observed data not present => will use with asimov data !");
            data = (RooDataSet*)ws->data("asimovData");
        }
    }

    // Loop on NPs to find gammas and add to the list to be ranked
    if(NPnames=="all" || NPnames.find("gamma")!=string::npos || (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0)){
        RooRealVar* var = NULL;
        RooArgSet* nuis = (RooArgSet*) mc->GetNuisanceParameters();
        if(nuis){
            TIterator* it2 = nuis->createIterator();
            int i_gamma = 0;
            while( (var = (RooRealVar*) it2->Next()) ){
                string np = var->GetName();
                if(np.find("gamma")!=string::npos){
                    // add the nuisance parameter to the list nuisPars if it's there in the ws
                    // remove "gamma"...
                    if(np==NPnames || (((atoi(NPnames.c_str())-(int)Nsyst-(int)Nnorm)==i_gamma) && (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0)) || NPnames=="all"){
                        nuisPars.push_back(ReplaceString(np,"gamma_",""));
                        isNF.push_back( true );
                        if(NPnames!="all") break;
                    }
                    i_gamma++;
                }
            }
        }
    }

    //
    // Create snapshot to keep inital values
    //
    ws -> saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));

    //
    // Initialize the FittingTool object
    //
    FittingTool *fitTool = new FittingTool();
    fitTool -> SetDebug(TRExFitter::DEBUGLEVEL);
    fitTool -> ValPOI(1.);
    fitTool -> ConstPOI(false);

    TRExFit *fit = fFitList[fFitList.size()-1];
    fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
    muhat = fit->fFitResults -> GetNuisParValue( fPOI );

    for(unsigned int i=0;i<nuisPars.size();i++){

        //Getting the postfit values of the nuisance parameter
        central = fit->fFitResults -> GetNuisParValue(   nuisPars[i] );
        up      = fit->fFitResults -> GetNuisParErrUp(   nuisPars[i] );
        down    = fit->fFitResults -> GetNuisParErrDown( nuisPars[i] );
        // for gammas
        if( (NPnames=="all" && nuisPars[i].find("_bin_")!=string::npos) ){
            nuisPars[i] = "gamma_" + nuisPars[i];
        }
        outName_file <<  nuisPars[i] << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";

        //Set the NP to its post-fit *up* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool -> ResetFixedNP();
        fitTool -> FixNP( nuisPars[i], central + TMath::Abs(up  ) );
        fitTool -> FitPDF( mc, simPdf, data, fFastFitForRanking );
        muVarUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
        //
        //Set the NP to its post-fit *down* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool -> ResetFixedNP();
        fitTool -> FixNP( nuisPars[i], central - TMath::Abs(down) );
        fitTool -> FitPDF( mc, simPdf, data, fFastFitForRanking );
        muVarDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        outName_file << muVarUp[nuisPars[i]]-muhat << "   " <<  muVarDown[nuisPars[i]]-muhat<< "  ";

        if(isNF[i]){
            muVarNomUp[   nuisPars[i] ] = muhat;
            muVarNomDown[ nuisPars[i] ] = muhat;
        }
        else{
            //Set the NP to its pre-fit *up* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            float prefitUp   = 1.;
            float prefitDown = 1.;
            fitTool -> ResetFixedNP();
            fitTool -> FixNP( nuisPars[i], central + prefitUp );
            fitTool -> FitPDF( mc, simPdf, data, fFastFitForRanking );
            muVarNomUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
            //
            //Set the NP to its pre-fit *down* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool -> ResetFixedNP();
            fitTool -> FixNP( nuisPars[i], central - prefitDown );
            fitTool -> FitPDF( mc, simPdf, data, fFastFitForRanking );
            //
            muVarNomDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        }
        outName_file << muVarNomUp[nuisPars[i]]-muhat << "   " <<  muVarNomDown[nuisPars[i]]-muhat<< " "<<endl;

    }
    outName_file.close();
    ws->loadSnapshot("tmp_snapshot");

    f->Close();
    delete f;
    delete fitTool;
}

//____________________________________________________________________________________
//
void MultiFit::PlotNPRankingManager(){
  if(fFitList[0]->fRankingPlot=="Merge"  || fFitList[0]->fRankingPlot=="all") PlotNPRanking(true,true);
  if(fFitList[0]->fRankingPlot=="Systs"  || fFitList[0]->fRankingPlot=="all") PlotNPRanking(true,false);
  if(fFitList[0]->fRankingPlot=="Gammas" || fFitList[0]->fRankingPlot=="all") PlotNPRanking(false,true);
}

//____________________________________________________________________________________
//
void MultiFit::PlotNPRanking(bool flagSysts, bool flagGammas){
    WriteInfoStatus("MultiFit::PlotNPRanking", "....................................");
    WriteInfoStatus("MultiFit::PlotNPRanking", "Plotting Ranking...");
    //
    string fileToRead = fOutDir+"/Fits/NPRanking.txt";
    //
    // trick to merge the ranking outputs produced in parallel:
    string cmd = " if [[ `ls "+fOutDir+"/Fits/NPRanking_*` != \"\" ]] ; then";
    cmd       += " if [[ `ls "+fOutDir+"/Fits/NPRanking.txt` == \"\" ]] ; then";
    cmd       += " cat "+fName+"/Fits/NPRanking_* > "+fileToRead+" ; ";
    cmd       += " fi ;";
    cmd       += " fi ;";
    gSystem->Exec(cmd.c_str());
    //
    unsigned int maxNP = fFitList[0]->fRankingMaxNP;
    //
    string paramname;
    double nuiphat;
    double nuiperrhi;
    double nuiperrlo;
    double PoiUp;
    double PoiDown;
    double PoiNomUp;
    double PoiNomDown;
    std::vector<string> parname;
    std::vector<double> nuhat;
    std::vector<double> nuerrhi;
    std::vector<double> nuerrlo;
    std::vector<double> poiup;
    std::vector<double> poidown;
    std::vector<double> poinomup;
    std::vector<double> poinomdown;
    std::vector<double> number;

    ifstream fin( fileToRead.c_str() );
    fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    if (paramname=="Luminosity"){
        WriteErrorStatus("MultiFit::PlotNPRanking", "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"");
        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    }
    while (!fin.eof()){
        if(paramname.find("stat")!=string::npos && !flagGammas){
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            if (paramname=="Luminosity"){
                WriteErrorStatus("MultiFit::PlotNPRanking", "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"");
                fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            }
            continue;
        }
        if(paramname.find("stat")==string::npos && !flagSysts){
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            if (paramname=="Luminosity"){
                WriteErrorStatus("MultiFit::PlotNPRanking", "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"");
                fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            }
            continue;
        }
        parname.push_back(paramname);
        nuhat.push_back(nuiphat);
        nuerrhi.push_back(nuiperrhi);
        nuerrlo.push_back(nuiperrlo);
        poiup.push_back(PoiUp);
        poidown.push_back(PoiDown);
        poinomup.push_back(PoiNomUp);
        poinomdown.push_back(PoiNomDown);
        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        if (paramname=="Luminosity"){
            WriteErrorStatus("MultiFit::PlotNPRanking", "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"");
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        }
    }

    unsigned int SIZE = parname.size();
    number.push_back(0.5);
    for (unsigned int i=1;i<SIZE;i++){
        number.push_back(i+0.5);
        double sumi = 0.0;
        int index=-1;
        sumi += TMath::Max( TMath::Abs(poiup[i]),TMath::Abs(poidown[i]) );
        for (unsigned int j=1;j<=i;j++){
            double sumii = 0.0;
            sumii += TMath::Max( TMath::Abs(poiup[i-j]),TMath::Abs(poidown[i-j]) );
            if (sumi<sumii){
                if (index==-1){
                    swap(poiup[i],poiup[i-j]);
                    swap(poidown[i],poidown[i-j]);
                    swap(poinomup[i],poinomup[i-j]);
                    swap(poinomdown[i],poinomdown[i-j]);
                    swap(nuhat[i],nuhat[i-j]);
                    swap(nuerrhi[i],nuerrhi[i-j]);
                    swap(nuerrlo[i],nuerrlo[i-j]);
                    swap(parname[i],parname[i-j]);
                    index=i-j;
                }
                else{
                    swap(poiup[index],poiup[i-j]);
                    swap(poidown[index],poidown[i-j]);
                    swap(poinomup[index],poinomup[i-j]);
                    swap(poinomdown[index],poinomdown[i-j]);
                    swap(nuhat[index],nuhat[i-j]);
                    swap(nuerrhi[index],nuerrhi[i-j]);
                    swap(nuerrlo[index],nuerrlo[i-j]);
                    swap(parname[index],parname[i-j]);
                    index=i-j;
                }
            }
            else{
                break;
            }
        }
    }
    number.push_back(parname.size()-0.5);

    double poimax = 0;
    for (unsigned int i=0;i<SIZE;i++) {
        poimax = TMath::Max(poimax,TMath::Max( TMath::Abs(poiup[i]),TMath::Abs(poidown[i]) ));
        poimax = TMath::Max(poimax,TMath::Max( TMath::Abs(poinomup[i]),TMath::Abs(poinomdown[i]) ));
        nuerrlo[i] = TMath::Abs(nuerrlo[i]);
    }
    poimax *= 1.2;

    for (unsigned int i=0;i<SIZE;i++) {
        poiup[i]     *= (2./poimax);
        poidown[i]   *= (2./poimax);
        poinomup[i]  *= (2./poimax);
        poinomdown[i]*= (2./poimax);
    }

    // Resttrict to the first N
    if(SIZE>maxNP) SIZE = maxNP;

    // Graphical part - rewritten taking DrawPulls in TRExFitter
    float lineHeight  =  30;
    float offsetUp    =  60; // external
    float offsetDown  =  60;
    float offsetUp1   = 100; // internal
    float offsetDown1 =  15;
    int offset = offsetUp + offsetDown + offsetUp1 + offsetDown1;
    int newHeight = offset + SIZE*lineHeight;

    float xmin = -2;
    float xmax =  2;
    float max  =  0;

    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    TGraphAsymmErrors *g1 = new TGraphAsymmErrors();
    TGraphAsymmErrors *g2 = new TGraphAsymmErrors();
    TGraphAsymmErrors *g1a = new TGraphAsymmErrors();
    TGraphAsymmErrors *g2a = new TGraphAsymmErrors();

    int idx = 0;
    std::vector< string > Names;
    Names.clear();
    string parTitle;

    for(unsigned int i = parname.size()-SIZE; i<parname.size(); ++i){
        g->SetPoint(idx, nuhat[i],  idx+0.5);
        g->SetPointEXhigh(      idx, nuerrhi[i]);
        g->SetPointEXlow(       idx, nuerrlo[i]);

        g1->SetPoint(      idx, 0.,idx+0.5);
        g1->SetPointEXhigh(idx, poiup[i]);
        g1->SetPointEXlow( idx, 0.);
        g1->SetPointEYhigh(idx, 0.4);
        g1->SetPointEYlow( idx, 0.4);

        g2->SetPoint(      idx, 0.,idx+0.5);
        g2->SetPointEXhigh(idx, poidown[i]);
        g2->SetPointEXlow( idx, 0.);
        g2->SetPointEYhigh(idx, 0.4);
        g2->SetPointEYlow( idx, 0.4);

        g1a->SetPoint(      idx, 0.,idx+0.5);
        g1a->SetPointEXhigh(idx, poinomup[i]);
        g1a->SetPointEXlow( idx, 0.);
        g1a->SetPointEYhigh(idx, 0.4);
        g1a->SetPointEYlow( idx, 0.4);

        g2a->SetPoint(      idx, 0.,idx+0.5);
        g2a->SetPointEXhigh(idx, poinomdown[i]);
        g2a->SetPointEXlow( idx, 0.);
        g2a->SetPointEYhigh(idx, 0.4);
        g2a->SetPointEYlow( idx, 0.4);
        if(parname[i].find("gamma")!=string::npos || parname[i].find("stat_")!=string::npos){
            // get name of the region
            std::vector<std::string> tmpVec = Vectorize(parname[i],'_');
            int nWords = tmpVec.size();
            std::string regName = tmpVec[2];
            for(int i_word=3;i_word<nWords-2;i_word++){
                regName += tmpVec[i_word];
            }
            // find the short label of this region
            std::string regTitle = regName;
            for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
                for( int i_ch = 0; i_ch < fFitList[i_fit]->fNRegions; i_ch++ ){
                    if(fFitList[i_fit]->fRegions[i_ch]->fName==regName){
                        regTitle = fFitList[i_fit]->fRegions[i_ch]->fShortLabel;
                        break;
                    }
                }
            }
            // build the title of the nuis par
            parTitle = "#gamma (" + regTitle + " bin " + tmpVec[nWords-1] + ")";
        }
        else parTitle = TRExFitter::SYSTMAP[ parname[i] ];

        if(parTitle==""){
            parTitle = parname[i];
        }

        Names.push_back(parTitle);

        idx ++;
        if(idx > max)  max = idx;
    }

    TCanvas *c = new TCanvas("c","c",600,newHeight);
    c->SetTicks(0,0);
    gPad->SetLeftMargin(0.4);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D *h_dummy = new TH1D("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy->SetMaximum( SIZE + offsetUp1/lineHeight   );
    h_dummy->SetMinimum(      - offsetDown1/lineHeight );
    h_dummy->SetLineWidth(0);
    h_dummy->SetFillStyle(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->SetFillColor(kWhite);
    h_dummy->GetYaxis()->SetLabelSize(0);
    h_dummy->Draw();
    h_dummy->GetYaxis()->SetNdivisions(0);
    for(int i_bin=0;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
        h_dummy->SetBinContent(i_bin,-10);
    }

    g1->SetFillColor(kAzure-4);
    g2->SetFillColor(kCyan);
    g1->SetLineColor(g1->GetFillColor());
    g2->SetLineColor(g2->GetFillColor());

    g1a->SetFillColor(kWhite);
    g2a->SetFillColor(kWhite);
    g1a->SetLineColor(kAzure-4);
    g2a->SetLineColor(kCyan);
    g1a->SetFillStyle(0);
    g2a->SetFillStyle(0);
    g1a->SetLineWidth(1);
    g2a->SetLineWidth(1);

    g->SetLineWidth(2);

    g1a->Draw("5 same");
    g2a->Draw("5 same");
    g1->Draw("2 same");
    g2->Draw("2 same");
    g->Draw("p same");

    TLatex *systs = new TLatex();
    systs->SetTextAlign(32);
    systs->SetTextSize( systs->GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs->DrawLatex(xmin-0.1,i+0.5,Names[i].c_str());
    }
    h_dummy->GetXaxis()->SetLabelSize( h_dummy->GetXaxis()->GetLabelSize()*0.9 );
    h_dummy->GetXaxis()->CenterTitle();
    h_dummy->GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy->GetXaxis()->SetTitleOffset(1.2);

    TGaxis *axis_up = new TGaxis( -2, SIZE + (offsetUp1)/lineHeight, 2, SIZE + (offsetUp1)/lineHeight, -poimax,poimax, 510, "-" );
    axis_up->SetLabelOffset( 0.01 );
    axis_up->SetLabelSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetLabelFont(   gStyle->GetTextFont() );
    axis_up->Draw();
    axis_up->CenterTitle();
    axis_up->SetTitle(("#Delta"+fPOIName).c_str());
    if(SIZE==20) axis_up->SetTitleOffset(1.5);
    if(SIZE==10) axis_up->SetTitleOffset(1.25);
    axis_up->SetTitleSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetTitleFont(   gStyle->GetTextFont() );

    TPad *pad1 = new TPad("p1","Pad High",0,(newHeight-offsetUp-offsetUp1)/newHeight,0.4,1);
    pad1->Draw();

    pad1->cd();
    TLegend *leg1 = new TLegend(0.02,0.7,1,1.0,("Pre-fit impact on "+fPOIName+":").c_str());
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetMargin(0.25);
    leg1->SetNColumns(2);
    leg1->SetTextFont(gStyle->GetTextFont());
    leg1->SetTextSize(gStyle->GetTextSize());
    leg1->AddEntry(g1a,"#theta = #hat{#theta}+#Delta#theta","f");
    leg1->AddEntry(g2a,"#theta = #hat{#theta}-#Delta#theta","f");
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.02,0.32,1,0.62,("Post-fit impact on "+fPOIName+":").c_str());
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.25);
    leg2->SetNColumns(2);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetTextSize(gStyle->GetTextSize());
    leg2->AddEntry(g1,"#theta = #hat{#theta}+#Delta#hat{#theta}","f");
    leg2->AddEntry(g2,"#theta = #hat{#theta}-#Delta#hat{#theta}","f");
    leg2->Draw();

    TLegend *leg0 = new TLegend(0.02,0.1,1,0.25);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetMargin(0.2);
    leg0->SetTextFont(gStyle->GetTextFont());
    leg0->SetTextSize(gStyle->GetTextSize());
    leg0->AddEntry(g,"Nuis. Param. Pull","lp");
    leg0->Draw();

    c->cd();

    TLine l0;
    TLine l1;
    TLine l2;
    l0 = TLine(0,- offsetDown1/lineHeight,0,SIZE+0.5);// + offsetUp1/lineHeight);
    l0.SetLineStyle(kDashed);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    l1 = TLine(-1,- offsetDown1/lineHeight,-1,SIZE+0.5);// + offsetUp1/lineHeight);
    l1.SetLineStyle(kDashed);
    l1.SetLineColor(kBlack);
    l1.Draw("same");
    l2 = TLine(1,- offsetDown1/lineHeight,1,SIZE+0.5);// + offsetUp1/lineHeight);
    l2.SetLineStyle(kDashed);
    l2.SetLineColor(kBlack);
    l2.Draw("same");

    ATLASLabelNew(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.6*offsetUp1)/newHeight), (char*)fFitList[0]->fAtlasLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(       0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.3*offsetUp1)/newHeight), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    gPad->RedrawAxis();

    if(flagGammas && flagSysts){
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fOutDir+"/Ranking."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else if(flagGammas){
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fOutDir+"/RankingGammas."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else if(flagSysts){
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fOutDir+"/RankingSysts."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else{
      WriteWarningStatus("MultiFit::PlotNPRanking", "Your ranking plot felt in unknown category :s");
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fOutDir+"/RankingUnknown."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    //
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::GetLikelihoodScan( RooWorkspace *ws, string varName, RooDataSet* data,bool recreate,bool compare){
    WriteInfoStatus("MultiFit::GetLikelihoodScan", "Running likelihood scan for the parameter = " + varName);
    TString LHDir("LHoodPlots/");
    
    // shut-up RooFit!
    if(TRExFitter::DEBUGLEVEL<=1){
        if(TRExFitter::DEBUGLEVEL<=0) gErrorIgnoreLevel = kError;
        else if(TRExFitter::DEBUGLEVEL<=1) gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL) ;
        RooMsgService::instance().getStream(1).removeTopic(Generation) ;
        RooMsgService::instance().getStream(1).removeTopic(Plotting) ;
        RooMsgService::instance().getStream(1).removeTopic(LinkStateMgmt) ;
        RooMsgService::instance().getStream(1).removeTopic(Eval) ;  
        RooMsgService::instance().getStream(1).removeTopic(Caching) ;  
        RooMsgService::instance().getStream(1).removeTopic(Optimization) ;  
        RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;  
        RooMsgService::instance().getStream(1).removeTopic(InputArguments) ;  
        RooMsgService::instance().getStream(1).removeTopic(Tracing) ;  
        RooMsgService::instance().getStream(1).removeTopic(Contents) ;  
        RooMsgService::instance().getStream(1).removeTopic(DataHandling) ;  
        RooMsgService::instance().setStreamStatus(1,false);
    }

    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());

    bool isPoI = false;
    RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
    TString firstPOIname = (TString)firstPOI->GetName();
    if (firstPOIname.Contains(varName.c_str())) isPoI = true;

    RooRealVar* var = NULL;
    TString vname = "";
    std::string vname_s = "";
    bool foundSyst = false;
    Double_t minVal = -3;
    Double_t maxVal =  3;
    for(auto fit : fFitList){
        for(auto nf : fit->fNormFactors){
            if(nf->fName == varName){
                minVal = nf->fMin;
                maxVal = nf->fMax;
            }
        }
    }

    if (fLHscanMin < 99999) { // is actually set
        minVal = fLHscanMin;
    }
    
    if (fLHscanMax > -99999) { // is actually set
        maxVal = fLHscanMax;
    }

    if (isPoI){
        TIterator* it = mc->GetParametersOfInterest()->createIterator();
        while( (var = (RooRealVar*) it->Next()) ){
            vname=var->GetName();
            vname_s=var->GetName();
            if (vname.Contains(varName.c_str())) {
                WriteInfoStatus("MultiFit::GetLikelihoodScan", "GetLikelihoodScan for POI = " + vname_s);
                foundSyst=true;
                break;
            }
        }
    }
    else {
        TIterator* it = mc->GetNuisanceParameters()->createIterator();
        while( (var = (RooRealVar*) it->Next()) ){
            vname=var->GetName();
            vname_s=var->GetName();
            if (vname.Contains(varName.c_str())) {
                WriteInfoStatus("MultiFit::GetLikelihoodScan", "GetLikelihoodScan for NP = " + vname_s);
                foundSyst=true;
                break;
            }
        }
    }

    if(!foundSyst){
        WriteWarningStatus("MultiFit::GetLikelihoodScan", "systematic " + varName + " not found (most probably due to Pruning), skip LHscan !");
        return;
    }
    WriteInfoStatus("MultiFit::GetLikelihoodScan", "GetLikelihoodScan for parameter = " + vname_s);

    TCanvas* can = new TCanvas("NLLscan");
    can->SetTopMargin(0.1);
    RooCurve* curve;
    RooPlot* frameLH = var->frame(Title("-log(L) vs "+vname),Bins(fLHscanSteps),Range(minVal, maxVal));
    
    if(recreate){
        RooAbsReal* nll = simPdf->createNLL(*data, Constrain(*mc->GetNuisanceParameters()), Offset(1), NumCPU(TRExFitter::NCPU, RooFit::Hybrid));
        TString tag("");
        RooAbsReal* pll = nll->createProfile(*var);
        pll->plotOn(frameLH,RooFit::Precision(-1),LineColor(kRed), NumCPU(TRExFitter::NCPU));
        curve = frameLH->getCurve();
    }
    else{
        TFile *f = new TFile(fName+"/"+LHDir+"NLLscan_"+varName+"_curve.root");
        curve = (RooCurve*)f->Get("LHscan");
        curve->SetLineColor(kRed);
        curve->SetLineWidth(3);
    }
    curve->Draw();
    
    // take the LH curves also for other fits
    std::vector<RooCurve*> curve_fit;
    std::vector<RooCurve*> curve_fit_statOnly; // to implement
    TLegend *leg = new TLegend(0.5,0.85-0.06*(fFitList.size()+1),0.75,0.85);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    if(compare){
        for(auto fit : fFitList){
            TFile *f = nullptr;
            if(fit->fFitResultsFile!=""){
                std::vector<std::string> v = Vectorize(fit->fFitResultsFile,'/');
                f = new TFile(v[0]+"/"+LHDir+"NLLscan_"+varName+"_curve.root");
            }
            else{
                f = new TFile(fit->fName+"/"+LHDir+"NLLscan_"+varName+"_curve.root");
            }
            if(f!=nullptr) {
                curve_fit.push_back((RooCurve*)f->Get("LHscan"));
            }
            else {
                curve_fit.push_back(nullptr);
            }
        }
        //
        int idx = 0;
        for(auto crv : curve_fit){
            if (!crv) continue;
            if(idx==0) crv->SetLineColor(kBlue);
            if(idx==1) crv->SetLineColor(kGreen);
            frameLH->addPlotable(crv,"same");
            leg->AddEntry(crv,fFitList[idx]->fLabel.c_str(),"l");
            idx++;
        }
        leg->AddEntry(curve,"Combined","l");
    }

    frameLH->GetXaxis()->SetRangeUser(minVal,maxVal);

    // y axis
    frameLH->GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") frameLH->GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") frameLH->GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());

    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(vname);

    can->SetTitle(cname);
    can->SetName(cname);
    can->cd();
    frameLH->Draw();

    TLatex *tex = new TLatex();
    tex->SetTextColor(kGray+2);
    
    TLine *l1s = new TLine(minVal,0.5,maxVal,0.5);
    l1s->SetLineStyle(kDashed);
    l1s->SetLineColor(kGray);
    l1s->SetLineWidth(2);
    if(frameLH->GetMaximum()>2){
        l1s->Draw();
        tex->DrawLatex(maxVal,0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }
    
    if(isPoI){
        if(frameLH->GetMaximum()>2){
            TLine *l2s = new TLine(minVal,2,maxVal,2);
            l2s->SetLineStyle(kDashed);
            l2s->SetLineColor(kGray);
            l2s->SetLineWidth(2);
            l2s->Draw();
            tex->DrawLatex(maxVal,2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
        }
        //
        if(frameLH->GetMaximum()>4.5){
            TLine *l3s = new TLine(minVal,4.5,maxVal,4.5);
            l3s->SetLineStyle(kDashed);
            l3s->SetLineColor(kGray);
            l3s->SetLineWidth(2);
            l3s->Draw();
            tex->DrawLatex(maxVal,4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
        }
        //
        TLine *lv0 = new TLine(0,frameLH->GetMinimum(),0,frameLH->GetMaximum());
        lv0->Draw();
        //
        TLine *lh0 = new TLine(minVal,0,maxVal,0);
        lh0->Draw();
    }

    system(TString("mkdir -vp ")+fName+"/"+LHDir);

    if(compare){
        leg->Draw();
        ATLASLabel(0.15,0.93,fFitList[0]->fAtlasLabel.c_str(),kBlack);
        myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
        if(fLabel!="") myText(0.2,0.85,kBlack,Form("#kern[-1]{%s}",fLabel.c_str()));
    }
    
    frameLH->SetMinimum(0);
    can->RedrawAxis();
    curve->Draw("same");
    
    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        can->SaveAs( fName+"/"+LHDir+"NLLscan_"+varName+"."+TRExFitter::IMAGEFORMAT[i_format] );
    
    if(recreate){
        // write it to a ROOT file as well
        TFile *f = new TFile(fName+"/"+LHDir+"NLLscan_"+varName+"_curve.root","UPDATE");
        f->cd();
        curve->Write("LHscan",TObject::kOverwrite);
        f->Close();
        delete f;
    }
}

//____________________________________________________________________________________
//
void MultiFit::PlotSummarySoverB(){
    WriteInfoStatus("MultiFit::PlotSummarySoverB", "....................................");
    WriteInfoStatus("MultiFit::PlotSummarySoverB", "Producing S/B plot...");

    bool includeBonly = false;
    if(fBonlySuffix!="") includeBonly = true;

    fFitList[0]->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
    float muFit = fFitList[0]->fFitResults->GetNuisParValue(fPOI);
    if (HistFromFile( fOutDir+"/Limits/"+fName+fSaveSuf+".root/limit" ) == nullptr) {
        WriteWarningStatus("TRExFit::PlotSummarySoverB", "Histo pointer is nullptr, skipping plotting.");
        return;
    }
    float muLimit = HistFromFile( fOutDir+"/Limits/"+fName+fSaveSuf+".root/limit" )->GetBinContent(1);

    std::vector<string> fileNames; fileNames.clear();
    std::vector<string> fileNamesBonly; fileNamesBonly.clear();
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        for(unsigned int i_reg=0;i_reg<fFitList[i_fit]->fRegions.size();i_reg++){
            if(fFitList[i_fit]->fRegions[i_reg]->fRegionType==Region::VALIDATION) continue;
            fileNames.push_back(fFitList[i_fit]->fName+"/Histograms/"+fFitList[i_fit]->fRegions[i_reg]->fName+"_postFit.root");
            if(includeBonly)
                fileNamesBonly.push_back(fFitList[i_fit]->fName+"/Histograms/"+fFitList[i_fit]->fRegions[i_reg]->fName+fBonlySuffix+"_postFit.root");
        }
    }
    int Nhist = (int)fileNames.size();

    //
    // create a list of all the samples
    std::vector<string> sigList; sigList.clear();
    std::vector<string> bkgList; bkgList.clear();
    std::vector<string> dataList; dataList.clear();
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        for(unsigned int i_smp=0;i_smp<fFitList[i_fit]->fSamples.size();i_smp++){
            if(fFitList[i_fit]->fSamples[i_smp]->fType==Sample::SIGNAL && FindInStringVector(sigList,fFitList[i_fit]->fSamples[i_smp]->fName)<0)
                sigList.push_back(fFitList[i_fit]->fSamples[i_smp]->fName);
            if(fFitList[i_fit]->fSamples[i_smp]->fType==Sample::BACKGROUND && FindInStringVector(bkgList,fFitList[i_fit]->fSamples[i_smp]->fName)<0)
                bkgList.push_back(fFitList[i_fit]->fSamples[i_smp]->fName);
            if(fFitList[i_fit]->fSamples[i_smp]->fType==Sample::DATA && FindInStringVector(dataList,fFitList[i_fit]->fSamples[i_smp]->fName)<0)
                dataList.push_back(fFitList[i_fit]->fSamples[i_smp]->fName);
        }
    }
    //
    // create a list of all the systematics
    std::vector<string> systList; systList.clear();
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        // actual systematics
        for(unsigned int i_syst=0;i_syst<fFitList[i_fit]->fSystematics.size();i_syst++){
            if(FindInStringVector(systList,fFitList[i_fit]->fSystematics[i_syst]->fName)<0)
                systList.push_back(fFitList[i_fit]->fSystematics[i_syst]->fName);
        }
        // norm factors
        for(unsigned int i_norm=0;i_norm<fFitList[i_fit]->fNormFactors.size();i_norm++){
            if(FindInStringVector(systList,fFitList[i_fit]->fNormFactors[i_norm]->fName)<0)
                systList.push_back(fFitList[i_fit]->fNormFactors[i_norm]->fName);
        }
    }

    unsigned int Nsyst = systList.size();

    std::vector<TFile*> file;
    std::vector<TFile*> fileBonly;
    std::vector<TH1D* > h_sig;
    std::vector<TH1D* > h_bkg;
    std::vector<TH1D* > h_bkgBonly;
    std::vector<TH1D* > h_tot_bkg_prefit;
    std::vector<TH1D* > h_data;
    std::vector<std::vector<TH1D*> > h_syst_up  (Nsyst,std::vector<TH1D*>(Nhist));
    std::vector<std::vector<TH1D*> > h_syst_down(Nsyst,std::vector<TH1D*>(Nhist));

    // get histos
    for(int i_hist=0;i_hist<Nhist;i_hist++){
        TH1D* h_tmp = nullptr;
        TH1D* h_tmpBonly = nullptr;
        WriteDebugStatus("MultiFit::PlotSummarySoverB",  "Opening file " + fileNames[i_hist]);
        file.push_back(new TFile(fileNames[i_hist].c_str()));
        if(includeBonly){
            WriteDebugStatus("MultiFit::PlotSummarySoverB", "Opening file " + fileNamesBonly[i_hist]);
            fileBonly.push_back(new TFile(fileNamesBonly[i_hist].c_str()));
        }
        //
        // initialize null histogram pointers
        h_sig.push_back(nullptr);
        h_bkg.push_back(nullptr);
        h_bkgBonly.push_back(nullptr);
        h_tot_bkg_prefit.push_back(nullptr);
        h_data.push_back(nullptr);
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            h_syst_up  [i_syst][i_hist] = nullptr;
            h_syst_down[i_syst][i_hist] = nullptr;
        }
        //
        for(unsigned int i_sig=0;i_sig<sigList.size();i_sig++){
            WriteDebugStatus("MultiFit::PlotSummarySoverB", "  Getting histogram h_"+sigList[i_sig]+"_postFit");
            h_tmp = (TH1D*)file[i_hist]->Get( ("h_"+sigList[i_sig]+"_postFit").c_str() );
            if(h_tmp!=nullptr){
                WriteDebugStatus("MultiFit::PlotSummarySoverB", " ... FOUND");
                if(h_sig[i_hist]==nullptr) h_sig[i_hist] = h_tmp;
                else                   h_sig[i_hist]->Add(h_tmp);
            }
        }
        for(unsigned int i_bkg=0;i_bkg<bkgList.size();i_bkg++){
            WriteDebugStatus("MultiFit::PlotSummarySoverB", "  Getting histogram h_"+bkgList[i_bkg]+"_postFit");
            h_tmp = (TH1D*)file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_postFit").c_str() );
            if(includeBonly) h_tmpBonly = (TH1D*)fileBonly[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_postFit").c_str() );
            if(h_tmp!=nullptr){
                WriteDebugStatus("MultiFit::PlotSummarySoverB", " ... FOUND");
                if(h_bkg[i_hist]==nullptr) h_bkg[i_hist] = h_tmp;
                else                   h_bkg[i_hist]->Add(h_tmp);
            }
            if(h_tmpBonly!=nullptr){
                WriteDebugStatus("MultiFit::PlotSummarySoverB", " ... B-only FOUND");
                if(h_bkgBonly[i_hist]==nullptr) h_bkgBonly[i_hist] = h_tmpBonly;
                else                        h_bkgBonly[i_hist]->Add(h_tmpBonly);
            }
            // syst variations
            for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
                // up
                h_tmp = (TH1D*)file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_"+systList[i_syst]+"_Up_postFit").c_str() );
                if(h_tmp!=nullptr){
                    if(h_syst_up[i_syst][i_hist]==nullptr)   h_syst_up[i_syst][i_hist] = h_tmp;
                    else                                 h_syst_up[i_syst][i_hist]->Add(h_tmp);
                }
                // down
                if(h_tmp!=nullptr){
                    h_tmp = (TH1D*)file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_"+systList[i_syst]+"_Down_postFit").c_str() );
                    if(h_syst_down[i_syst][i_hist]==nullptr) h_syst_down[i_syst][i_hist] = h_tmp;
                    else                                 h_syst_down[i_syst][i_hist]->Add(h_tmp);
                }
            }
        }
        for(unsigned int i_data=0;i_data<dataList.size();i_data++){
            h_tmp = (TH1D*)file[i_hist]->Get( ("h_"+dataList[i_data]).c_str() );
            if(h_tmp!=nullptr){
                if(h_data[i_hist]==nullptr) h_data[i_hist] = h_tmp;
                else                    h_data[i_hist]->Add(h_tmp);
            }
        }

        if(TRExFitter::PREFITONPOSTFIT)
          h_tot_bkg_prefit[i_hist] = (TH1D*)file[i_hist]->Get("h_tot_bkg_prefit");

        //
        // Fix eventually empty histograms
        if(h_sig[i_hist] ==nullptr){
            h_sig[i_hist]  = (TH1D*)h_bkg[i_hist]->Clone(Form("h_sig[%d]", i_hist));
            h_sig[i_hist]->Scale(0.);
        }
        if(h_data[i_hist]==nullptr){
            h_data[i_hist] = (TH1D*)h_bkg[i_hist]->Clone(Form("h_data[%d]",i_hist));
            h_data[i_hist]->Scale(0.);
        }
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            // up
            if(h_syst_up[i_syst][i_hist]==nullptr){
                h_syst_up[i_syst][i_hist] = (TH1D*)h_bkg[i_hist]->Clone(Form("h_syst_up[%d][%d]", i_syst,i_hist));
                h_syst_up[i_syst][i_hist]->Scale(0.);
            }
            else{
                h_syst_up[i_syst][i_hist]->Add(h_bkg[i_hist],-1);
            }
            // down
            if(h_syst_down[i_syst][i_hist]==nullptr){
                h_syst_down[i_syst][i_hist] = (TH1D*)h_bkg[i_hist]->Clone(Form("h_syst_down[%d][%d]", i_syst,i_hist));
                h_syst_down[i_syst][i_hist]->Scale(0.);
            }
            else{
                h_syst_down[i_syst][i_hist]->Add(h_bkg[i_hist],-1);
            }
        }
    }

    // create combined histogram
    TH1D* h_bkg_comb  = Combine(h_bkg);
    TH1D* h_bkgBonly_comb = nullptr; if(includeBonly) h_bkgBonly_comb = Combine(h_bkgBonly);
    TH1D* h_sig_comb  = Combine(h_sig);
    TH1D* h_data_comb = Combine(h_data);
    TH1D* h_tot_bkg_prefit_comb = nullptr;
    if(TRExFitter::PREFITONPOSTFIT) h_tot_bkg_prefit_comb = Combine(h_tot_bkg_prefit);


    std::vector<TH1D*> h_syst_up_comb  (Nsyst);
    std::vector<TH1D*> h_syst_down_comb(Nsyst);
    for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
        h_syst_up_comb  [i_syst] = Combine(h_syst_up  [i_syst]);
        h_syst_down_comb[i_syst] = Combine(h_syst_down[i_syst]);
    }

    std::vector<float> SoverSqrtB;
    float sig, bkg;

    for(int i_bin=1;i_bin<=h_bkg_comb->GetNbinsX();i_bin++){
        sig = h_sig_comb->GetBinContent(i_bin);
        bkg = h_bkg_comb->GetBinContent(i_bin);
        SoverSqrtB.push_back(sig/bkg);
  //      SoverSqrtB.push_back(sig/sqrt(bkg));
    }

    TH1D* h_bkg_ord  = Rebin(h_bkg_comb,SoverSqrtB,false);
    TH1D* h_bkgBonly_ord = nullptr; if(includeBonly) h_bkgBonly_ord = Rebin(h_bkgBonly_comb,SoverSqrtB,false);
    TH1D* h_sig_ord  = Rebin(h_sig_comb,SoverSqrtB,false);
    TH1D* h_data_ord = Rebin(h_data_comb,SoverSqrtB);
    TH1D* h_tot_bkg_prefit_ord = nullptr;
    if(TRExFitter::PREFITONPOSTFIT) h_tot_bkg_prefit_ord = Rebin(h_tot_bkg_prefit_comb,SoverSqrtB,false);

    std::vector<TH1D*> h_syst_up_ord  (Nsyst);
    std::vector<TH1D*> h_syst_down_ord(Nsyst);
    for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
        h_syst_up_ord  [i_syst] = Rebin((TH1D*)(h_syst_up_comb  [i_syst]),SoverSqrtB,false);
        h_syst_down_ord[i_syst] = Rebin((TH1D*)(h_syst_down_comb[i_syst]),SoverSqrtB,false);
    }

    float errUp, errDown, err, err_tot;
    float corr;
    for(int i_bin=0;i_bin<h_bkg_ord->GetNbinsX()+2;i_bin++){
        err_tot = h_bkg_ord->GetBinError(i_bin); // this should be the stat unc
        errUp   = 0;
        errDown = 0;
        err = 0;
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            for(unsigned int j_syst=0;j_syst<systList.size();j_syst++){
                corr     = fFitList[0]->fFitResults->fCorrMatrix->GetCorrelation( systList[i_syst],systList[j_syst] );
                errUp   += corr * h_syst_up_ord  [i_syst]->GetBinContent(i_bin) * h_syst_up_ord  [j_syst]->GetBinContent(i_bin);
                errDown += corr * h_syst_down_ord[i_syst]->GetBinContent(i_bin) * h_syst_down_ord[j_syst]->GetBinContent(i_bin);
            }
        }
        errUp   = sqrt(errUp);
        errDown = sqrt(errDown);
        err = TMath::Abs(errUp+errDown)/2.;
        err_tot = sqrt(err_tot*err_tot + err*err);
        h_bkg_ord->SetBinError(i_bin,err_tot);
    }

    TCanvas *c = new TCanvas("c","c",600,600);

    TPad* pad0;
    TPad* pad1;
    pad0 = new TPad("pad0","pad0",0,0.28,1,1,0,0,0);
    pad0->SetTicks(1,1);
    pad0->SetTopMargin(0.05);
    pad0->SetBottomMargin(0);
    pad0->SetLeftMargin(0.14);
    pad0->SetRightMargin(0.05);
    pad0->SetFrameBorderMode(0);
    //
    pad1 = new TPad("pad1","pad1",0,0,1,0.28,0,0,0);
    pad1->SetTicks(1,1);
    pad1->SetTopMargin(0.0);
    pad1->SetBottomMargin(0.37);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->SetFrameBorderMode(0);

    pad1->Draw();
    pad0->Draw();
    pad0->cd();

    h_sig_ord->SetLineColor(kRed);
    h_sig_ord->SetFillColor(kRed);

    TH1D* h_sig_ord_lim = (TH1D*)h_sig_ord->Clone("h_sig_ord_lim");
    h_sig_ord_lim->Scale(muLimit/muFit);
    h_sig_ord_lim->SetFillColor(kOrange);
    h_sig_ord_lim->SetLineColor(kOrange);
    TH1D* h_sig_ord_lim_diff = (TH1D*)h_sig_ord_lim->Clone("h_sig_ord_lim_diff");
    h_sig_ord_lim_diff->Add(h_sig_ord,-1);

    THStack *h_s = new THStack();
    h_s->Add(h_bkg_ord);
    h_s->Add(h_sig_ord);
    h_s->Add(h_sig_ord_lim_diff);
    h_data_ord->Draw("EX0");

    h_s->Draw("HISTsame ][");
    TH1D* h_err = (TH1D*)h_bkg_ord->Clone("h_err");
    h_err->SetMarkerSize(0);
    h_err->SetFillColor(kBlack);
    h_err->SetFillStyle(3454);
    h_err->SetLineWidth(0);
    h_err->SetLineColor(kWhite);
    h_err->Draw("E2same");
    h_data_ord->Draw("EX0same");
    h_data_ord->SetMaximum(20*h_data_ord->GetMaximum());
    h_data_ord->SetMinimum(50);
    h_data_ord->SetLineWidth(2);
    h_data_ord->GetXaxis()->SetTitle("log_{10}(S/B)");
    h_data_ord->GetYaxis()->SetTitle("Events / 0.2");
    h_data_ord->GetYaxis()->SetTitleOffset(2.);
    h_data_ord->GetXaxis()->SetLabelSize(0);
    h_data_ord->GetXaxis()->SetTitleSize(0);

    if(includeBonly){
        h_bkgBonly_ord->SetLineColor(kBlack);
        h_bkgBonly_ord->SetLineStyle(kDashed);
        h_bkgBonly_ord->Draw("HISTsame ][");
    }

    if(TRExFitter::PREFITONPOSTFIT) {
      h_tot_bkg_prefit_ord->SetLineColor(kBlue);
      h_tot_bkg_prefit_ord->SetLineStyle(kDashed);
      h_tot_bkg_prefit_ord->Draw("HISTsame ][");
    }

    TLegend *leg;
    if(includeBonly) leg = new TLegend(0.6,0.50,0.90,0.92);
    else             leg = new TLegend(0.6,0.57,0.90,0.92);
    leg->SetFillStyle(0);
    leg->SetMargin(0.2);
    leg->SetBorderSize(0);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->AddEntry(h_data_ord,"Data","lep");
    leg->AddEntry(h_sig_ord_lim,Form(("%s ("+fPOIName+"_{95%% excl.}=%.1f)").c_str(),fSignalTitle.c_str(),muLimit),"f");
    leg->AddEntry(h_sig_ord,    Form(("%s ("+fPOIName+"_{fit}=%.1f)"       ).c_str(),fSignalTitle.c_str(),muFit),"f");
    leg->AddEntry(h_bkg_ord,"Background","f");
    leg->AddEntry(h_err,"Bkgd. Unc.","f");
    if(includeBonly) leg->AddEntry(h_bkgBonly_ord,("Bkgd. ("+fPOIName+"=0)").c_str(),"l");
    if(TRExFitter::PREFITONPOSTFIT) leg->AddEntry(h_tot_bkg_prefit_comb,"Pre-Fit Bkgd.","l");
    leg->Draw();

    ATLASLabelNew(0.17,0.87, (char*)fFitList[0]->fAtlasLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(0.17,0.80,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()) );
    if(fLabel!="") myText(0.17,0.18,kBlack,Form("%s Combined",fLabel.c_str()) );
    else           myText(0.17,0.18,kBlack,"Combined");
    std::string channels = "";
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        if(i_fit!=0){
            if(i_fit==fFitList.size()-1) channels += " and ";
            else                         channels += ", ";
        }
        channels += fFitList[i_fit]->fLabel;
    }
    myText(0.17,0.13,kBlack,channels.c_str());
    myText(0.17,0.05,kBlack,"Post-Fit");

    pad0->RedrawAxis();
    pad0->SetLogy();

    pad1->cd();
    pad1->GetFrame()->SetY1(2);
    TH1D *h_ratio   = (TH1D*)h_data_ord->Clone("h_ratio");
    TH1D *h_den     = (TH1D*)h_bkg_ord ->Clone("h_den");
    for(int i_bin=0;i_bin<h_den->GetNbinsX()+2;i_bin++){
        h_den->SetBinError(i_bin,0);
    }

    TH1D *h_ratioBonly = nullptr;
    if(includeBonly){
        h_ratioBonly = (TH1D*)h_bkgBonly_ord->Clone("h_ratioBonly");
        h_ratioBonly->Divide(h_den);
        h_ratioBonly->SetLineStyle(kDashed);
        h_ratioBonly->SetLineColor(kBlack);
    }

    TH1D* h_stackSig = (TH1D*)h_sig_ord ->Clone("h_sig_ratio");
    h_stackSig->Add(h_bkg_ord);
    h_stackSig->Divide(h_den);
    h_stackSig->SetFillColor(0);
    h_stackSig->SetFillStyle(0);
    h_stackSig->SetLineColor(kRed);

    TH1D* h_stackSigLim = (TH1D*)h_sig_ord_lim ->Clone("h_sig_lim_ratio");
    h_stackSigLim->Add(h_bkg_ord);
    h_stackSigLim->Divide(h_den);
    h_stackSigLim->SetFillColor(0);
    h_stackSigLim->SetFillStyle(0);
    h_stackSigLim->SetLineStyle(kDashed);
    h_stackSigLim->SetLineColor(kOrange+1);

    TH1D *h_ratio2  = (TH1D*)h_err->Clone("h_ratio2");
    h_ratio2->SetMarkerSize(0);
    h_ratio->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetTitle("Data / Bkgd.");
    h_ratio->GetYaxis()->SetTitleSize(20);
    h_ratio->GetYaxis()->SetTitleOffset(2.);
    h_ratio->GetYaxis()->SetLabelSize(20); // 0.04
    h_ratio ->Divide(h_den);
    h_ratio2->Divide(h_den);
    h_ratio->SetMarkerSize(1.2);
    h_ratio->SetLineWidth(2);
    gStyle->SetEndErrorSize(0.); // 4.
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->GetYaxis()->SetNdivisions(406);
    h_ratio->SetMinimum(0.6);
    h_ratio->SetMaximum(1.75);
    h_ratio->GetXaxis()->SetTitle(h_data_ord->GetXaxis()->GetTitle());
    h_ratio->GetXaxis()->SetTitleSize(20);
    h_ratio->GetXaxis()->SetTitleOffset(4.);
    h_ratio->GetXaxis()->SetLabelSize(20);
    TLine* hline = new TLine(h_ratio->GetXaxis()->GetXmin(),1,h_ratio->GetXaxis()->GetXmax(),1);
    h_ratio->Draw("E1");
    h_ratio2->Draw("same E2");
    hline->Draw();
    h_stackSig->Draw("same HIST ][");
    h_stackSigLim->Draw("same HIST ][");
    h_ratio->Draw("same E1");

    if(includeBonly) h_ratioBonly->Draw("same HIST ][");

    TLegend *leg2 = new TLegend(0.17,0.64,0.75,0.98);
    leg2->SetFillStyle(0);
    leg2->SetMargin(0.1);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(gStyle->GetTextSize());
    leg2->AddEntry(h_stackSigLim,Form(("%s ("+fPOIName+"_{95%% excl.}=%.1f) + Bkgd.").c_str(),fSignalTitle.c_str(),muLimit),"l");
    leg2->AddEntry(h_stackSig,   Form(("%s ("+fPOIName+"_{fit}=%.1f) + Bkgd."       ).c_str(),fSignalTitle.c_str(),muFit)  ,"l");
    leg2->Draw();

    if(includeBonly){
        TLegend *leg3 = new TLegend(0.17,0.4,0.75,0.5);
        leg3->SetFillStyle(0);
        leg3->SetMargin(0.1);
        leg3->SetBorderSize(0);
        leg3->SetTextSize(gStyle->GetTextSize());
        leg3->AddEntry(h_ratioBonly,"Bkgd. (from Bkgd-only fit)","l");
        leg3->Draw();
    }


    pad0->RedrawAxis();
    pad1->RedrawAxis();

    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fOutDir+"/SoverB_postFit."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );

    delete c;
}

//____________________________________________________________________________________
//
TH1D* MultiFit::Combine(vector<TH1D*> h){
    int Nbins = 0;
    int Nhist = h.size();
    for(int i_hist=0;i_hist<Nhist;i_hist++){
        if(h[i_hist]==nullptr) WriteWarningStatus("MultiFit::Combine", "empty histgram " + std::to_string(i_hist));
        else Nbins += h[i_hist]->GetNbinsX();
    }
    TH1D* h_new = new TH1D(Form("%s_comb",h[0]->GetName()),Form("%s_comb",h[0]->GetTitle()),Nbins,0,Nbins);
    int bin = 0;
    for(int i_hist=0;i_hist<Nhist;i_hist++){
        for(int i_bin=1;i_bin<=h[i_hist]->GetNbinsX();i_bin++){
            bin++;
            h_new->SetBinContent(bin,h[i_hist]->GetBinContent(i_bin));
            h_new->SetBinError(bin,h[i_hist]->GetBinError(i_bin));
        }
    }
    return h_new;
}

//____________________________________________________________________________________
// order bins of h acording to a[] (increasing order)
TH1D* MultiFit::OrderBins(TH1D* h,vector<float> vec){
    map<float,int> binIndex;
    int Nbins = h->GetNbinsX();
    for(int i_bin=1;i_bin<=Nbins;i_bin++){
        binIndex[vec[i_bin-1]] = i_bin;
    }
    sort(vec.begin(),vec.end());
    TH1D *h_new = (TH1D*)h->Clone();
    for(int i_bin=1;i_bin<=Nbins;i_bin++){
        h_new->SetBinContent(i_bin,h->GetBinContent(binIndex[vec[i_bin-1]]));
    }
    return h_new;
}

//____________________________________________________________________________________
// merge bins in bins of SoverSqrtB
TH1D* MultiFit::Rebin(TH1D* h,vector<float> vec, bool isData){
    TH1D* h_new = new TH1D(Form("%s_rebin",h->GetName()),Form("%s_rebin",h->GetTitle()),17,-3.8,-0.5);
    h_new->Sumw2();
    // new way
    for(int j_bin=1;j_bin<=h->GetNbinsX();j_bin++){
        float value=TMath::Log10(vec[j_bin-1]);
        if ( value<h_new->GetXaxis()->GetXmin() ) value=0.9999*h_new->GetXaxis()->GetXmin();
        if ( value>h_new->GetXaxis()->GetXmax() ) {
            double tmpvalue=1.0001*h_new->GetXaxis()->GetXmax();
            WriteDebugStatus("MultiFit::Rebin", "turning: " + std::to_string(value) + " in: " + std::to_string(tmpvalue));
            value=tmpvalue;
        }
        int i_bin=h_new->FindBin(value);
        h_new->SetBinContent(i_bin,h_new->GetBinContent(i_bin)+h->GetBinContent(j_bin));
        if (!isData) {
            h_new->SetBinError(i_bin,sqrt( pow(h_new->GetBinError(i_bin),2) + pow(h->GetBinError(j_bin),2) ) );
        }
    }
    if (isData) {
        for(int j_bin=1;j_bin<=h_new->GetNbinsX();j_bin++){
            h_new->SetBinError(j_bin, sqrt(h_new->GetBinContent(j_bin) ) );
        }
    }
    h_new->SetMinimum(1);
    return h_new;
}

//____________________________________________________________________________________
// combine individual results from grouped impact evaluation into one table
void MultiFit::BuildGroupedImpactTable(){
    WriteInfoStatus("TRExFit::BuildGroupedImpactTable", "merging grouped impact evaluations");
    std::string targetName = fOutDir+"/Fits/GroupedImpact"+fSaveSuf+".txt";

    if(std::ifstream(targetName).good()){
        WriteWarningStatus("TRExFit::BuildGroupedImpactTable","file " + targetName + " already exists, will not overwrite");
    }
    else{
        std::string cmd = " if [[ `ls "+fOutDir+"/Fits/GroupedImpact"+fSaveSuf+"_*` != \"\" ]] ; then";
        cmd            += " cat "+fOutDir+"/Fits/GroupedImpact_* > "+targetName+" ; ";
        cmd            += " fi ;";
        gSystem->Exec(cmd.c_str());
    }
}
