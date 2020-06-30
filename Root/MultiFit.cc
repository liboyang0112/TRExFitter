// Class include
#include "TRExFitter/MultiFit.h"

// Framework inncludes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/LikelihoodScanManager.h"
#include "TRExFitter/LimitEvent.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/RankingManager.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/TRExFit.h"
#include "TRExFitter/YamlConverter.h"

// CommonStatTools include
#include "CommonStatTools/runSig.h"
#include "CommonStatTools/runAsymptoticsCLs.h"

// Roofit includes
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooMinimizer.h"
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
#include "TGraph2D.h"
#include "TH2F.h"
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
MultiFit::MultiFit(const string& name) :
    fCombine(false),
    fCompare(false),
    fStatOnly(false),
    fIncludeStatOnly(false),
    fCompareLimits(true),
    fComparePOI(true),
    fComparePulls(true),
    fPlotCombCorrMatrix(false),
    fName(name),
    fDir(""),
    fOutDir(""),
    fLabel(""),
    fShowObserved(false),
    fLimitTitle("95% CL limit on XXX"),
    fPOITitle("best fit XXX"),
    fRankingOnly("all"),
    fGroupedImpactCategory("all"),
    fPOI(""),
    fPOIMin(0),
    fPOIMax(10),
    fPOIPrecision("1"),
    fLimitMax(0),
    fUseRnd(false),
    fRndRange(0.1),
    fRndSeed(-999),
    fLumiLabel(""),
    fCmeLabel(""),
    fCombiLabel("Combined"),
    fConfig(std::unique_ptr<ConfigParser>(new ConfigParser())),
    fSaveSuf(""),
    fDataName("obsData"),
    fFitType(1), // 1: S+B, 2: B-only
    fCombineChByCh(true),
    fFastFit(false),
    fFastFitForRanking(true),
    fNuisParListFile(""),
    fPlotSoverB(false),
    fSignalTitle("signal"),
    fFitResultsFile(""),
    fLimitsFile(""),
    fBonlySuffix(""),
    fShowSystForPOI(false),
    fGetGoodnessOfFit(false),
    fLHscanMin(999999),
    fLHscanMax(-999999),
    fLHscanSteps(30),
    fLHscanMinY(999999),
    fLHscanMaxY(-999999),
    fLHscanStepsY(30),
    fParal2D(false),
    fParal2Dstep(-1),
    fDoGroupedSystImpactTable(false),
    fPOIName("#mu"),
    fPOINominal(1),
    fPOIAsimov(1),
    fLimitIsBlind(false),
    fSignalInjection(false),
    fSignalInjectionValue(0),
    fLimitParamName("parameter"),
    fLimitParamValue(0),
    fLimitOutputPrefixName("myLimit"),
    fLimitsConfidence(0.95),
    fSignificanceIsBlind(false),
    fSignificancePOIAsimov(0),
    fSignificanceParamName("parameter"),
    fSignificanceParamValue(0),
    fSignificanceOutputPrefixName("mySignificance"),
    fShowTotalOnly(false),
    fuseGammasForCorr(false),
    fPOIInitial(1.),
    fHEPDataFormat(false),
    fFitStrategy(-1),
    fCPU(1),
    fBinnedLikelihood(false)
{
    fNPCategories.emplace_back("");
}

//__________________________________________________________________________________
//
void MultiFit::AddFitFromConfig(const std::string& configFile,
                                const std::string& opt,
                                const std::string& options,
                                const std::string& label,
                                const std::string& loadSuf,
                                const std::string& wsFile,
                                const bool useInFit,
                                const bool useInComparison) {

    // check if the config is not already processed (but it might be intended, if comparing different fits from same config)
    if (std::find(fConfigPaths.begin(), fConfigPaths.end(), configFile) != fConfigPaths.end()){
        WriteWarningStatus("MultiFit::AddFitFromConfig", "Config " + configFile + " is added twice. Make sure you know what you are doing."); // changed from error to warning, since in some cases one might want to include the same job twice (e.g. comparing fit results with different suffix)
    }

    fConfigPaths.emplace_back(configFile);

    // keep debug level
    const int debug = TRExFitter::DEBUGLEVEL;

    fFitList.emplace_back(new TRExFit());

    fFitList.back()->fUseInFit = useInFit;

    fFitList.back()->fUseInComparison = useInComparison;

    // initialize config reader
    ConfigReader reader(fFitList.back());

    if (reader.ReadFullConfig(configFile,opt,options) != 0){
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
RooWorkspace* MultiFit::CombineWS() const{
    WriteInfoStatus("MultiFit::CombineWS", "....................................");
    WriteInfoStatus("MultiFit::CombineWS", "Combining workspaces...");

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);

    gSystem->mkdir((fOutDir+"/RooStats").c_str());

    std::vector < RooWorkspace* > vec_ws;
    std::vector < std::string > vec_chName;
    RooStats::HistFactory::Measurement *measurement = nullptr;

    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        if (!fFitList.at(i_fit)->fUseInFit) continue;
        const std::string& fitName = fFitList[i_fit]->fInputName;
        const std::string& fitDir = fFitList[i_fit]->fName;
        WriteDebugStatus("MultiFit::CombineWS", "Adding Fit: " + fitName + ", " + fFitLabels[i_fit] + ", " + fFitSuffs[i_fit] + fitDir);

        std::string fileName = fitDir + "/RooStats/" + fitName + "_combined_" + fitName + fFitSuffs[i_fit] + "_model.root";
        if(fWsFiles[i_fit]!="") fileName = fWsFiles[i_fit];
        WriteDebugStatus("MultiFit::CombineWS", "Opening file " + fileName );
        TFile* rootFile = TFile::Open(fileName.c_str(),"read");
        if (!rootFile) {
            WriteErrorStatus("MultiFit::CombineWS", "Cannot open file: " + fileName);
            exit(EXIT_FAILURE);
        }
        RooWorkspace* m_ws = dynamic_cast<RooWorkspace*>(rootFile->Get("combined"));
        WriteDebugStatus("MultiFit::CombineWS", "Getting " + fitName+fFitSuffs[i_fit] );
        RooStats::HistFactory::Measurement* meas = dynamic_cast<RooStats::HistFactory::Measurement*>(rootFile -> Get( (fitName+fFitSuffs[i_fit]).c_str()));
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
                const std::string fileName_tmp = fitDir + "/RooStats/" + fitName + "_" + reg->fName + "_" + fitName + fFitSuffs[i_fit] + "_model.root";
                WriteDebugStatus("MultiFit::CombineWS", "  Opening file " + fileName_tmp );
                TFile* rootFile_tmp = TFile::Open(fileName_tmp.c_str(),"read");
                if (!rootFile_tmp) {
                    WriteErrorStatus("MultiFit::CombineWS", "Cannot open file: " + fileName_tmp);
                    exit(EXIT_FAILURE);
                }

                RooWorkspace* m_ws_tmp = dynamic_cast<RooWorkspace*>(rootFile_tmp->Get(reg->fName.c_str()));
                WriteDebugStatus("MultiFit::CombineWS", "  Getting " + reg->fName );
                vec_ws.push_back(m_ws_tmp);
                vec_chName.push_back(reg->fName);
            }
        }
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
    
    if (fBinnedLikelihood) {
        FitUtils::SetBinnedLikelihoodOptimisation(ws);
    }

    WriteInfoStatus("MultiFit::CombineWS", "....................................");

    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws, *measurement );

    measurement->PrintXML((fOutDir + "/RooStats/").c_str());

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    return ws;
}

//__________________________________________________________________________________
//
void MultiFit::SaveCombinedWS() const{
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    //
    // Creating the rootfile
    //
    std::unique_ptr<TFile> f(TFile::Open( (fOutDir+"/ws_combined"+fSaveSuf+".root").c_str(), "recreate"));
    if (!f) {
        WriteErrorStatus("MultiFit::SaveCombinedWS", "Cannot open file!");
        exit(EXIT_FAILURE);
    }
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
std::map < std::string, double > MultiFit::FitCombinedWS(int fitType, const std::string& inputData, bool doLHscanOnly) const {
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    std::unique_ptr<TFile> f(TFile::Open((fOutDir+"/ws_combined"+fSaveSuf+".root").c_str()));
    if (!f) {
        WriteErrorStatus("MultiFit::FitCombinedWS", "Cannot open file!");
        exit(EXIT_FAILURE);
    }
    RooWorkspace *ws = dynamic_cast<RooWorkspace*>(f->Get("combWS"));

    std::map < std::string, double > result;

    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////

    //
    // Fit configuration (1: SPLUSB or 2: BONLY)
    //
    FittingTool fitTool{};
    fitTool.SetUseHesse(true);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetDebug(TRExFitter::DEBUGLEVEL);
    fitTool.SetNCPU(fCPU);
    if(fitType==2){
        fitTool.ValPOI(0.);
        fitTool.ConstPOI(true);
    } else if(fitType==1){
        fitTool.ValPOI(fPOIInitial);
        fitTool.ConstPOI(false);
    }
    if(fUseRnd) fitTool.SetRandomNP(fRndRange, fUseRnd, fRndSeed);

    //
    // Fit starting from custom point
    if(fFitResultsFile!=""){
        fFitList[0]->ReadFitResults(fFitResultsFile);
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitList[0]->fFitResults->fNuisPar) {
            npNames.push_back( inp->fName);
            npValues.push_back(inp->fFitValue);
        }
        fitTool.SetNPs(npNames,npValues);
    }
    // Fix NPs that are specified in the individual configs
    for (const auto& ifit : fFitList){
        if (!ifit->fUseInFit) continue;
        if(ifit->fFitFixedNPs.size()>0){
            std::vector<std::string> npNames;
            std::vector<double> npValues;
            for(const auto& nuisParToFix : ifit->fFitFixedNPs){
                npNames.push_back( nuisParToFix.first );
                npValues.push_back( nuisParToFix.second );
            }
            fitTool.FixNPs(npNames,npValues);
        }
    }

    std::vector<std::string> vVarNameMinos;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& iminos : ifit->fVarNameMinos) {
            if(Common::FindInStringVector(vVarNameMinos,iminos) < 0){
                vVarNameMinos.push_back(iminos);
            }
        }
    }

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    if(vVarNameMinos.size()>0){
        WriteDebugStatus("MultiFit::FitCombinedWS", "Setting the variables to use MINOS with:");
        for(const auto& iminos : vVarNameMinos) {
            WriteDebugStatus("MultiFit::FitCombinedWS",  "  " + iminos);
        }
        fitTool.UseMinos(vVarNameMinos);
    }

    //
    // Gets needed objects for the fit
    //
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>((mc->GetPdf()));
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    if(inputData=="asimovData"){
        RooArgSet empty;// = RooArgSet();
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        RooRealVar * poi = static_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
        if (!poi){
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            WriteErrorStatus("TRExFit::DumpData", "Cannot find POI in workspace, exiting...");
            exit(EXIT_FAILURE);
        }
        poi->setVal(fPOIAsimov);
        data = static_cast<RooDataSet*>(RooStats::AsymptoticCalculator::MakeAsimovData( (*mc), RooArgSet(ws->allVars()), static_cast<RooArgSet&>(empty)));
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        
    }
    else if(inputData!=""){
        data = static_cast<RooDataSet*>(ws->data(inputData.c_str()));
    } else {
        WriteWarningStatus("MultiFit::FitCombinedWS", "You didn't specify inputData => will try with observed data !");
        data = static_cast<RooDataSet*>(ws->data("obsData"));
        if(!data){
            WriteWarningStatus("MultiFit::FitCombinedWS", "Observed data not present => will use with asimov data !");
            data = static_cast<RooDataSet*>(ws->data("asimovData"));
        }
    }

    if (!data){
        WriteErrorStatus("MultiFit::FitCombinedWS", "Data returns null ptr, probably wrong name in DataName?");
        exit(EXIT_FAILURE);
    }

    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, GetFitNormFactors());

    // Performs the fit
    gSystem -> mkdir((fOutDir+"/Fits/").c_str(),true);
    
    // save snapshot before fit
    ws->saveSnapshot("snapshot_BeforeFit_POI", *(mc->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_BeforeFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_BeforeFit_GO" , *(mc->GetGlobalObservables())    );
    
    double nll0 = 0.;
    if(fGetGoodnessOfFit) nll0 = fitTool.FitPDF( mc, simPdf, static_cast<RooDataSet*>(ws->data("asimovData")), false, true );
    
    // save snapshot after fit
    ws->saveSnapshot("snapshot_AfterFit_POI", *(mc->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_AfterFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(mc->GetGlobalObservables())    );
    
    //
    // Get number of degrees of freedom
    // - number of bins
    int ndof = data->numEntries();
    // - minus number of free & non-constant parameters
    std::vector<std::string> nfList;
    for(const auto& fit : fFitList){
        if (!fit->fUseInFit) continue;
        for(const auto& nf : fit->fNormFactors){
            if(nf->fConst) continue;
            if(Common::FindInStringVector(nfList,nf->fName)>=0) continue;
            if(fFitType==2 && fPOI==nf->fName) continue;
            nfList.push_back(nf->fName);
        }
    }
    ndof -= nfList.size();

    double nll = 0;
    if (!doLHscanOnly){
        nll = fitTool.FitPDF( mc, simPdf, data, fFastFit );
        std::vector<std::string> s_vec;
        fitTool.ExportFitResultInTextFile(fOutDir+"/Fits/"+fName+fSaveSuf+".txt", s_vec);
        result = fitTool.ExportFitResultInMap();
    }


    // Get initial ikelihood value from Asimov
    if(fGetGoodnessOfFit) {
        ws->loadSnapshot("snapshot_BeforeFit_POI");
        ws->loadSnapshot("snapshot_BeforeFit_GO");
        ws->loadSnapshot("snapshot_BeforeFit_NP");
        
        nll0 = fitTool.FitPDF( mc, simPdf, data, false, false, true);
    }

    //
    // Goodness of fit
    if(fGetGoodnessOfFit && !doLHscanOnly){
        const double deltaNLL = nll-nll0;
        const double prob = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndof);
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

    if (fGetGoodnessOfFit && !fDoGroupedSystImpactTable) {
        ws->loadSnapshot("snapshot_AfterFit_POI");
        ws->loadSnapshot("snapshot_AfterFit_GO");
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }

    //
    // grouped systematics impact
    if(fDoGroupedSystImpactTable && !doLHscanOnly){
        std::string outNameGroupedImpact = fOutDir+"/Fits/GroupedImpact"+fSaveSuf;
        if(fGroupedImpactCategory!="all") outNameGroupedImpact += "_"+fGroupedImpactCategory;
        outNameGroupedImpact += ".txt";
        // need to create a merged list of fSubCategoryImpactMap from the include Fits
        std::map<std::string, std::string> mergedMap;
        for(auto fit : fFitList){
            if (!fit->fUseInFit) continue;
            fit->ProduceSystSubCategoryMap();
            for(auto m : fit->fSubCategoryImpactMap) {
                if(mergedMap[m.first] == "") {
                    mergedMap[m.first] = m.second;
                } else if(mergedMap[m.first] != m.second){
                    WriteWarningStatus("MultiFit::FitCombinedWS","Systematics assigned to different SubCategory in the different included Fits. Keeping first Fit convention.");
                }
            }
        }
        fitTool.SetSystMap( mergedMap );
        fitTool.GetGroupedImpact( mc, simPdf, data, ws, fGroupedImpactCategory, outNameGroupedImpact);
    }

    //
    // Calls the  function to create LH scan with respect to a parameter
    //
    if(fVarNameLH.size()>0 && !doLHscanOnly && !fParal2D){ // Skip 1Dscan when paralelizing 2D
        if (fVarNameLH[0]=="all"){
            for(map<string,string>::iterator it=TRExFitter::SYSTMAP.begin(); it!=TRExFitter::SYSTMAP.end(); ++it){
                GetLikelihoodScan( ws, it->first, data, true);
            }
        } else{
            for(const auto& iLH : fVarNameLH){
                GetLikelihoodScan( ws, iLH, data, true);
            }
        }
    }
    if (doLHscanOnly && !fParal2D){ // Skip 1Dscan when paralelizing 2D
        if (fVarNameLH.size() == 0){
            WriteErrorStatus("MultiFit::MultiFit","Did not provide any LH scan parameter and running LH scan only. This is not correct.");
            exit(EXIT_FAILURE);
        }
        if (fVarNameLH[0]=="all"){
            WriteWarningStatus("MultiFit::MultiFit","You are running LHscan only option but running it for all parameters. Will not parallelize!.");
            for(map<string,string>::iterator it=TRExFitter::SYSTMAP.begin(); it!=TRExFitter::SYSTMAP.end(); ++it){
                GetLikelihoodScan( ws, it->first, data, true);
            }
        } else {
            GetLikelihoodScan( ws, fVarNameLH[0], data, true);
        }
    }
    // run 2D likelihood scan
    if(fVarName2DLH.size()>0){
        for (const auto & ipair : fVarName2DLH) {
            Get2DLikelihoodScan( ws, ipair, data);
        }
    }

    // Stat-only fit:
    // - read fit resutls
    // - fix all NP to fitted ones before fitting
    if(fIncludeStatOnly && !doLHscanOnly){
        WriteInfoStatus("MultiFit::FitCombinedWS", "Fitting stat-only: reading fit results from full fit from file:");
        WriteInfoStatus("MultiFit::FitCombinedWS", "  " + (fOutDir+"/Fits/"+fName+fSaveSuf+".txt"));
        fFitList[0]->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitList[0]->fFitResults->fNuisPar) {
            bool isNF = false;
            for(const auto& ifit : fFitList) {
                if (!ifit->fUseInFit) continue;
                if(!ifit->fFixNPforStatOnlyFit &&
                    Common::FindInStringVector(ifit->fNormFactorNames,inp->fName)>=0){
                    isNF = true;
                    break;
                }
            }
            if(isNF) continue;
            npNames.push_back( inp->fName);
            npValues.push_back(inp->fFitValue);
        }
        fitTool.ResetFixedNP();
        fitTool.FixNPs(npNames,npValues);
        fitTool.FitPDF( mc, simPdf, data );
        std::vector<std::string> s_vecTemp;
        fitTool.ExportFitResultInTextFile(fOutDir+"/Fits/"+fName+fSaveSuf+"_statOnly.txt", s_vecTemp);
    }

    return result;
}
//__________________________________________________________________________________
//
void MultiFit::GetCombinedLimit(string inputData) const{
    WriteInfoStatus("MultiFit::GetCombinedLimit", "Running runAsymptoticsCLs macro...");

    string wsFileName = fOutDir+"/ws_combined"+fSaveSuf+".root";
    int sigDebug = 3 - TRExFitter::DEBUGLEVEL;
    if (sigDebug < 0) sigDebug = 0;
    runAsymptoticsCLs(wsFileName.c_str(), "combWS", "ModelConfig", inputData.c_str(), fLimitParamName.c_str(), fLimitParamValue, fLimitOutputPrefixName.c_str(), (fOutDir+"/Limits/").c_str(), fLimitIsBlind, fLimitsConfidence, "asimovData_0", fSignalInjection, fSignalInjectionValue, sigDebug);
}
//__________________________________________________________________________________
//
void MultiFit::GetCombinedSignificance(string inputData) const{
    WriteInfoStatus("MultiFit::GetCombinedSignificance", "Running runSig macro...");

    string wsFileName = fOutDir+"/ws_combined"+fSaveSuf+".root";

    //
    // Finally computing the significance
    //
    int sigDebug = 3 - TRExFitter::DEBUGLEVEL;
    if (sigDebug < 0) sigDebug = 0;
    runSig(wsFileName.c_str(), "combWS", "ModelConfig", inputData.c_str(), fSignificanceParamName.c_str(), fSignificanceParamValue, fSignificanceOutputPrefixName.c_str(), (fOutDir+"/Significance").c_str(), fSignificanceIsBlind, "asimovData_1", "conditionalGlobs_1", "nominalGlobs", false, fSignificancePOIAsimov, sigDebug);
}
//__________________________________________________________________________________
//
void MultiFit::ComparePOI(const string& POI) const {
    double xmax = fPOIMax + (fPOIMax-fPOIMin);
    double xmin = fPOIMin;

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
        if(fFitList[i_fit]->fPOIs.size()>0) pois.push_back( fFitList[i_fit]->fPOIs[0] );
        else                                pois.push_back( POI );
    }
    if(fCombine){
        WriteInfoStatus("MultiFit::ComparePOI", "Adding Combined Fit");
        names.push_back( fName );
        dirs.push_back( fOutDir );
        titles.push_back( fCombiLabel );
        suffs.push_back( "" );
        pois.push_back( POI );
    }

    int N = names.size();

    double ymin = -0.5;
    double ymax = N+1-0.5;

    TCanvas c("c","c",700,500);
    gStyle->SetEndErrorSize(6.);

    TGraph g_central(N);
    TGraphAsymmErrors g_stat(N);
    TGraphAsymmErrors g_tot(N);

    int Ndiv = N+1;

    NuisParameter *par;

    // get values
    TRExFit *fit = nullptr;
    for(int i=0;i<N;i++){
        if(!(fCombine && i==N-1)){
            if (!fFitList.at(i)->fUseInComparison) continue;
        }
        const bool isComb = (fCombine && i==N-1) ? true : false;
        //
        if(!isComb){
            fit = fFitList[i];
            if(fit->fFitResultsFile=="") fit->ReadFitResults(dirs[i]+"/Fits/"+names[i]+suffs[i]+".txt");
            else                         fit->ReadFitResults(fit->fFitResultsFile);
        }
        else{
            if(fFitResultsFile=="") fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
            else                    fit->ReadFitResults(fFitResultsFile);
        }
        bool found(false);
        for(unsigned int j = 0; j<fit->fFitResults->fNuisPar.size(); ++j){
            par = fit->fFitResults->fNuisPar[j].get();
            if( pois[i] == par->fName ){
                g_central.SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_stat   .SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_tot    .SetPoint(N-i-1,par->fFitValue,N-i-1);
                //
                // temporary put the full uncertainty
                g_stat.SetPointEXhigh(N-i-1,par->fPostFitUp);
                g_stat.SetPointEXlow(N-i-1,-par->fPostFitDown);
                g_stat.SetPointEYhigh(N-i-1,0);
                g_stat.SetPointEYlow(N-i-1,0);
                //
                g_tot.SetPointEXhigh(N-i-1,par->fPostFitUp);
                g_tot.SetPointEXlow(N-i-1,-par->fPostFitDown);
                g_tot.SetPointEYhigh(N-i-1,0);
                g_tot.SetPointEYlow(N-i-1,0);
                //
                found = true;
                break;
            }
        }
        if(!found){
            g_central.SetPoint(N-i-1,-10,N-i-1);
            g_stat   .SetPoint(N-i-1,-10,N-i-1);
            g_tot    .SetPoint(N-i-1,-10,N-i-1);
            g_stat   .SetPointError(N-i-1,0,0,0,0);
            g_tot    .SetPointError(N-i-1,0,0,0,0);
        }
    }
    // stat error
    if (!fShowTotalOnly){
        for(int i=0;i<N;i++){
            const bool isComb = (fCombine && i==N-1) ? true : false;
            //
            if(!isComb){
                fit = fFitList[i];
                if(fit->fFitResultsFile=="") fit->ReadFitResults(dirs[i]+"/Fits/"+names[i]+suffs[i]+"_statOnly.txt");
                else                         fit->ReadFitResults(Common::ReplaceString(fit->fFitResultsFile,".txt","_statOnly.txt"));
            }
            else{
                if(fFitResultsFile=="")      fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+"_statOnly.txt");
                else                         fit->ReadFitResults(Common::ReplaceString(fFitResultsFile,".txt","_statOnly.txt"));
            }
            for(unsigned int j = 0; j<fit->fFitResults->fNuisPar.size(); ++j){
                par = fit->fFitResults->fNuisPar[j].get();
                if( pois[i] == par->fName ){
                    g_stat.SetPointEXhigh(N-i-1,par->fPostFitUp);
                    g_stat.SetPointEXlow(N-i-1,-par->fPostFitDown);
                    g_stat.SetPointEYhigh(N-i-1,0);
                    g_stat.SetPointEYlow(N-i-1,0);
                    break;
                }
            }
        }
    }

    g_tot .SetLineWidth(3);
    g_stat.SetLineWidth(3);
    g_stat.SetMarkerStyle(kOpenCircle);
    g_tot .SetLineWidth(3);
    g_tot .SetMarkerStyle(kOpenCircle);
    if(TRExFitter::OPTION["FourTopStyle"]){
        g_tot .SetLineColor(kAzure);
        g_tot .SetFillColor(kAzure);
        g_stat.SetLineColor(kCyan);
        g_stat.SetFillColor(kCyan);
        g_stat.SetMarkerStyle(kFullCircle);
        g_stat.SetMarkerColor(kAzure);
        g_stat.SetMarkerSize(1.5);
        //
        for(int i=0;i<N;i++){
            g_tot .SetPointEYlow(i,0.02);
            g_tot .SetPointEYhigh(i,0.02);
            g_stat.SetPointEYlow(i,0.04);
            g_stat.SetPointEYhigh(i,0.04);
        }
    }
    else{
        g_stat.SetLineColor(kGreen-8);
        g_stat.SetMarkerSize(0);
        g_tot .SetLineColor(kBlack);
    }
    if(TRExFitter::OPTION["FourTopStyle"]){
        g_central.SetMarkerColor(kWhite);
        g_central.SetMarkerStyle(kFullCircle);
        g_central.SetMarkerSize(1.);
    }
    else{
        g_central.SetMarkerStyle(kFullCircle);
        g_central.SetMarkerColor(kRed);
        g_central.SetMarkerSize(1.5);
    }
    g_tot.SetMarkerSize(0);

    TH1D h_dummy("h_dummy","h_dummy",1,xmin,xmax);
    h_dummy.Draw();
    h_dummy.SetMinimum(ymin);
    h_dummy.SetMaximum(ymax);
    h_dummy.SetLineColor(kWhite);
    h_dummy.GetYaxis()->Set(N+1,ymin,ymax);
    h_dummy.GetYaxis()->SetNdivisions(Ndiv);

    TLatex tex{};

    for(int i=0;i<N;i++){
        if(!(fCombine && i==N-1)){
            if (!fFitList.at(i)->fUseInComparison) continue;
        }
        h_dummy.GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
        if(fShowSystForPOI){
            tex.SetTextSize(gStyle->GetTextSize()*1.2);
            tex.DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form(("#font[62]{%." + fPOIPrecision + "f}").c_str(),g_central.GetX()[N-i-1]));
            tex.DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{^{#plus%." + fPOIPrecision + "f}}").c_str(),g_tot.GetErrorXhigh(N-i-1)));
            tex.DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{_{#minus%." + fPOIPrecision + "f}}").c_str(),g_tot.GetErrorXlow(N-i-1)));
            tex.DrawLatex(xmin+0.69*(xmax-xmin),N-i-1,"(");
            if (!fShowTotalOnly){
                tex.DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + fPOIPrecision + "f}}").c_str(),g_stat.GetErrorXhigh(N-i-1)));
                tex.DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + fPOIPrecision + "f}}").c_str(),g_stat.GetErrorXlow(N-i-1)));
                tex.DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + fPOIPrecision + "f}}").c_str(),
                    std::sqrt( g_tot.GetErrorXhigh(N-i-1) * g_tot.GetErrorXhigh(N-i-1) - g_stat.GetErrorXhigh(N-i-1)* g_stat.GetErrorXhigh(N-i-1) ) ) );
                tex.DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + fPOIPrecision + "f}}").c_str(),
                    std::sqrt( g_tot.GetErrorXlow(N-i-1)*g_tot.GetErrorXlow(N-i-1) - g_stat.GetErrorXlow(N-i-1)*g_stat.GetErrorXlow(N-i-1) ) ) );
                tex.DrawLatex(xmin+0.94*(xmax-xmin),N-i-1,")");
            }
        }
        else{
            tex.DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form((fPOIName+" = %." + fPOIPrecision + "f").c_str(),g_central.GetX()[N-i-1]));
            tex.DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("^{#plus%." + fPOIPrecision + "f}").c_str(),g_tot.GetErrorXhigh(N-i-1)));
            tex.DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("_{#minus%." + fPOIPrecision + "f}").c_str(),g_tot.GetErrorXlow(N-i-1)));
            if (!fShowTotalOnly){
                tex.DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("^{#plus%." + fPOIPrecision + "f}").c_str(),g_stat.GetErrorXhigh(N-i-1)));
                tex.DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("_{#minus%." + fPOIPrecision + "f}").c_str(),g_stat.GetErrorXlow(N-i-1)));
            }
        }
    }

    TLine l_SM(fPOINominal,-0.5,fPOINominal,N-0.5);
    l_SM.SetLineWidth(2);
    l_SM.SetLineColor(kGray);
    l_SM.Draw("same");

    TLine l_h(xmin,0.5,xmax,0.5);
    if(fCombine){
        l_h.SetLineWidth(2);
        l_h.SetLineColor(kBlack);
        l_h.SetLineStyle(kDashed);
        l_h.Draw("same");
    }

    if(TRExFitter::OPTION["FourTopStyle"]){
        g_tot.Draw("E2 same");
        if (!fShowTotalOnly){
            g_stat.Draw("PE2 same");
        }
        g_central.Draw("P same");
    }
    else{
        g_tot.Draw("E same");
        if (!fShowTotalOnly){
            g_stat.Draw("E same");
        }
        g_central.Draw("P same");
    }

    gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
    gPad->SetBottomMargin( 1.15*gPad->GetBottomMargin() );
    gPad->SetTopMargin( 1.8*gPad->GetTopMargin() );
    h_dummy.GetXaxis()->SetTitle(fPOITitle.c_str());
    h_dummy.GetYaxis()->SetTickSize(0);

    c.RedrawAxis();

    if (fFitList[0]->fAtlasLabel != "none") ATLASLabel(0.32,0.93,fFitList[0]->fAtlasLabel.c_str(),kBlack);
    myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(process!="") myText(0.94,0.85,kBlack,Form("#kern[-1]{%s}",process.c_str()));

    TLegend leg(0.35,0.775,0.7,0.9);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(&g_tot,"tot.","l");
    if (!fShowTotalOnly){
        leg.AddEntry(&g_stat,"stat.","l");
    }
    leg.Draw();

    if(fShowSystForPOI){
        tex.DrawLatex(xmin+0.6*(xmax-xmin),N-0.4,"#font[62]{tot}");
        tex.DrawLatex(xmin+0.69*(xmax-xmin),N-0.4,"(");
        if (!fShowTotalOnly){
            tex.DrawLatex(xmin+0.72*(xmax-xmin),N-0.4,"stat");
        }
        tex.DrawLatex(xmin+0.83*(xmax-xmin),N-0.4,"syst");
        tex.DrawLatex(xmin+0.94*(xmax-xmin),N-0.4,")");
    }
    else{
        tex.DrawLatex(xmin+(0.7-0.02)*(xmax-xmin),N-0.4,"( tot )");
        if (!fShowTotalOnly){
            tex.DrawLatex(xmin+(0.85-0.02)*(xmax-xmin),N-0.4,"( stat )");
        }
    }

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fOutDir+"/POI"+fSaveSuf+"."+format).c_str() );
    }
}

//__________________________________________________________________________________
//
void MultiFit::CompareLimit(){
    float xmax = 2.;
    const std::string process = fLabel;
    gStyle->SetEndErrorSize(0.);

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    for(unsigned int i_fit= 0; i_fit < fFitList.size(); ++i_fit) {
        WriteInfoStatus("MultiFit::CompareLimit", "Adding Fit: " + fFitList[i_fit]->fInputName + ", " + fFitLabels[i_fit] + ", " + fFitSuffs[i_fit]);
        dirs.emplace_back(fFitList[i_fit]->fName);
        names.emplace_back(fFitList[i_fit]->fInputName);
        titles.emplace_back(fFitLabels[i_fit]);
        suffs.emplace_back(fFitSuffs[i_fit]);
    }
    if(fCombine){
        WriteInfoStatus("MultiFit::CompareLimit", "Adding combined limit");
        dirs.emplace_back(fOutDir);
        names.emplace_back(fName);
        titles.emplace_back(fCombiLabel);
        suffs.emplace_back(fSaveSuf);
    }

    // ---

    const std::size_t N = names.size();

    const double ymin = -0.5;
    const double ymax = N-0.5;

    TCanvas c("c","c",700,500);

    TGraphErrors g_obs(N);
    TGraphErrors g_exp(N);
    TGraphErrors g_inj(N);
    TGraphAsymmErrors g_1s(N);
    TGraphAsymmErrors g_2s(N);

    const std::size_t Ndiv = N+1;

    std::unique_ptr<TFile> f = nullptr;

    // get values
    for(unsigned int i = 0; i < N; ++i) {
        if(i > fLimitsFiles.size()-1){
            if(fLimitsFile != "") fLimitsFiles.emplace_back(fLimitsFile);
            else                  fLimitsFiles.emplace_back("");
        }
        if(fLimitsFiles[i] == ""){
            const std::string name = fFitList.at(0)->fLimitOutputPrefixName + "_CL" + std::to_string(static_cast<int>(100*fFitList.at(0)->fLimitsConfidence));
            WriteInfoStatus("MultiFit::CompareLimit", "Reading file " + dirs[i] + "/Limits/asymptotics/" + name + ".root");
            f.reset(TFile::Open((dirs[i]+"/Limits/asymptotics/"+name+".root").c_str()));
        }
        else{
            WriteInfoStatus("MultiFit::CompareLimit", "Reading file " + fLimitsFiles[i]);
            f.reset(TFile::Open(fLimitsFiles[i].c_str()));
        }
        if (!f) {
            WriteWarningStatus("MultiFit::CompareLimit", "Cannot open the file!");
            continue;
        }
        TTree *tree = static_cast<TTree*>(f->Get("stats"));
        if (!tree) {
            WriteWarningStatus("MultiFit::CompareLimit", "Cannot read tree from the limit file");
            continue;
        }

        if (tree->GetEntries() == 0) {
            WriteWarningStatus("MultiFit::CompareLimit", "The tree has 0 entries");
            continue;
        }
        
        LimitEvent event(tree);

        // read the first event = read all info
        event.GetEntry(0);

        if(fShowObserved) g_obs.SetPoint(N-i-1,event.obs_upperlimit,N-i-1);
        else g_obs.SetPoint(N-i-1,-1,N-i-1);
        g_exp.SetPoint(N-i-1,event.exp_upperlimit,N-i-1);
        if(fSignalInjection) g_inj.SetPoint(N-i-1,event.inj_upperlimit,N-i-1);
        g_1s.SetPoint(N-i-1,event.exp_upperlimit,N-i-1);
        g_2s.SetPoint(N-i-1,event.exp_upperlimit,N-i-1);
        g_obs.SetPointError(N-i-1,0,0.5);
        g_exp.SetPointError(N-i-1,0,0.5);
        g_inj.SetPointError(N-i-1,0,0.5);
        g_1s.SetPointError(N-i-1,event.exp_upperlimit-event.exp_upperlimit_minus1,event.exp_upperlimit_plus1-event.exp_upperlimit,0.5,0.5);
        g_2s.SetPointError(N-i-1,event.exp_upperlimit-event.exp_upperlimit_minus2,event.exp_upperlimit_plus2-event.exp_upperlimit,0.5,0.5);

        const std::vector<float> tmp = {event.obs_upperlimit, event.exp_upperlimit,
                                        event.exp_upperlimit_plus1, event.exp_upperlimit_minus1,
                                        event.exp_upperlimit_plus2, event.exp_upperlimit_minus2, xmax};

        xmax = *std::max_element(tmp.begin(), tmp.end());

        delete tree;
        f->Close();
    }

    xmax*= 1.1;

    g_obs.SetLineWidth(3);
    g_exp.SetLineWidth(3);
    g_exp.SetLineStyle(2);
    g_inj.SetLineWidth(3);
    g_inj.SetLineStyle(2);
    g_inj.SetLineColor(kRed);
    g_1s.SetFillColor(kGreen);
    g_1s.SetLineWidth(3);
    g_1s.SetLineStyle(2);
    g_2s.SetFillColor(kYellow);
    g_2s.SetLineWidth(3);
    g_2s.SetLineStyle(2);

    g_2s.SetMarkerSize(0);
    g_1s.SetMarkerSize(0);
    g_exp.SetMarkerSize(0);
    g_obs.SetMarkerSize(0);
    g_inj.SetMarkerSize(0);

    if(fLimitMax!=0) xmax = fLimitMax;

    TH1D h_dummy("h_dummy","h_dummy",1,0,xmax);
    h_dummy.Draw();
    h_dummy.SetMinimum(ymin);
    h_dummy.SetMaximum(ymax);
    h_dummy.SetLineColor(kWhite);
    h_dummy.GetYaxis()->Set(N,ymin,ymax);
    h_dummy.GetYaxis()->SetNdivisions(Ndiv);
    for(unsigned int i = 0; i < N; ++i){
        h_dummy.GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
    }

    g_2s.Draw("E2 same");
    g_1s.Draw("E2 same");
    g_exp.Draw("E same");
    if(fShowObserved) g_obs.Draw("E same");
    if(fSignalInjection) g_inj.Draw("E same");

    TLine l_SM(fPOINominal,-0.5,fPOINominal,N-0.5);
    l_SM.SetLineWidth(2);
    l_SM.SetLineColor(kGray);
    l_SM.Draw("same");

    c.RedrawAxis();

    gPad->SetLeftMargin(2*gPad->GetLeftMargin());
    gPad->SetBottomMargin(1.15*gPad->GetBottomMargin());
    gPad->SetTopMargin(1.8*gPad->GetTopMargin());
    h_dummy.GetXaxis()->SetTitle(fLimitTitle.c_str());

    if (fFitList[0]->fAtlasLabel != "none") ATLASLabel(0.32,0.93,fFitList[0]->fAtlasLabel.c_str(),kBlack);
    myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(process!="") myText(0.94,0.85,kBlack,Form("#kern[-1]{%s}",process.c_str()));

    std::unique_ptr<TLegend> leg = nullptr;
    if(fShowObserved) leg = std::make_unique<TLegend>(0.65,0.2,0.95,0.40);
    else              leg = std::make_unique<TLegend>(0.65,0.2,0.95,0.35);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(&g_1s,"Expected #pm 1#sigma","lf");
    leg->AddEntry(&g_2s,"Expected #pm 2#sigma","lf");
    if(fShowObserved)    leg->AddEntry(&g_obs,"Observed","l");
    if(fSignalInjection) leg->AddEntry(&g_inj,("Expected ("+fPOIName+"=1)").c_str(),"l");
    leg->Draw();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs( (fOutDir+"/Limits" + fSaveSuf +  + "."+format).c_str() );
    }
}

//__________________________________________________________________________________
//
void MultiFit::ComparePulls(string category) const{
    double ydist = 0.2;

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    vector<double>  yshift;
    static const std::vector<int> color = {kBlack,kRed,kBlue,kViolet,kOrange};

    static const std::vector<int> style = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown,kOpenDiamond};

    unsigned int N = fFitList.size();
    if(fCombine) N++;

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            WriteInfoStatus("MultiFit::ComparePulls", "Adding Combined Fit");
            dirs.push_back( fOutDir );
            names.push_back( fName );
            titles.push_back( fCombiLabel );
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

    double xmin = -2.9;
    double xmax = 2.9;
    double max = 0.;

    // create a list of Systematics
    std::vector< string > Names;
    std::vector< string > Titles;
    std::vector< string > Categories;
    std::string systName;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        for(const auto& isyst : fFitList[i_fit]->fSystematics) {
            systName = isyst->fNuisanceParameter;
            if(systName == "") systName = isyst->fName;
            if(Common::FindInStringVector(Names,systName)<0){
                Names.push_back(systName);
                Titles.push_back(isyst->fTitle);
                if(Common::FindInStringVector(fFitList[i_fit]->fDecorrSysts,isyst->fName)>=0){
                    Titles[Titles.size()-1] += fFitList[i_fit]->fDecorrSuff;
                }
                Categories.push_back(isyst->fCategory);
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
    std::vector<string> NamesNew;
    std::vector<string> TitlesNew;
    std::vector<string> CategoriesNew;
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        FitResults *fitRes;
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            fitRes = fFitList[i_fit]->fFitResults;
            for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
                par = fitRes->fNuisPar[j].get();
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
    std::vector< TGraphAsymmErrors> g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,double> centralMap;
        std::map<string,double> errUpMap;
        std::map<string,double> errDownMap;
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            std::vector<std::string> s;
            if(fFitResultsFile!="") fitRes->ReadFromTXT(fFitResultsFile, s);
            else                    fitRes->ReadFromTXT(fOutDir+"/Fits/"+fName+fSaveSuf+".txt", s);
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults;
        }
        for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
            par = fitRes->fNuisPar[j].get();
            systName = par->fName;
            centralMap[systName] = par->fFitValue;
            errUpMap[systName]   = par->fPostFitUp;
            errDownMap[systName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.emplace_back(Nsyst);
        for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
            systName = Names[i_syst];
            if(centralMap[systName]!=0 || (errUpMap[systName]!=0 || errDownMap[systName]!=0)){
                g[i_fit].SetPoint(i_syst,centralMap[systName],(Nsyst-i_syst-1)+0.5+yshift[i_fit]);
                g[i_fit].SetPointEXhigh(i_syst,  errUpMap[systName]);
                g[i_fit].SetPointEXlow( i_syst, -errDownMap[systName]);
            }
            else{
                g[i_fit].SetPoint(i_syst,-10,-10);
                g[i_fit].SetPointEXhigh(i_syst, 0);
                g[i_fit].SetPointEXlow( i_syst, 0);
            }
        }
    }

    max = Nsyst;

    int lineHeight = 20;
    int offsetUp = 50;
    int offsetDown = 60;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0(0,0,0,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    TBox b1(-1,0,1,max);
    TBox b2(-2,0,2,max);
    b1.SetFillColor(kGreen);
    b2.SetFillColor(kYellow);
    b2.Draw("same");
    b1.Draw("same");
    l0.Draw("same");

    int i_style = 0;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if (!(fCombine && i_fit==N-1)){
            if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        g[i_fit].SetLineColor(color[i_style]);
        g[i_fit].SetMarkerColor(color[i_style]);
        g[i_fit].SetMarkerStyle(style[i_style]);
        g[i_fit].Draw("P same");
        i_style ++;
    }

    TLatex systs{};
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        systs.DrawLatex(3.,(Nsyst-i_syst-1)+0.25,Titles[i_syst].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    h_dummy.GetXaxis()->CenterTitle();
    h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy.GetXaxis()->SetTitleOffset(1.2);

    TLegend leg(0.01,1.-0.03*(30./max),0.75,0.99);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if (!(fCombine && i_fit==N-1)){
            if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        leg.AddEntry(&g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg.Draw();

    gPad->RedrawAxis();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        if(category=="") c.SaveAs((fOutDir+"/NuisPar_comp"+fSaveSuf+"."+format).c_str());
        else             c.SaveAs((fOutDir+"/NuisPar_comp"+fSaveSuf+"_"+category+"."+format).c_str());
    }
}

//__________________________________________________________________________________
//
void MultiFit::CompareNormFactors(string category) const{
    double ydist = 0.2;

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    vector<double>  yshift;

    static const std::vector<int> color = {kBlack,kRed,kBlue,kViolet,kOrange};
    static const std::vector<int> style = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown,kOpenDiamond};

    unsigned int N = fFitList.size();
    if(fCombine) N++;

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            WriteInfoStatus("MultiFit::CompareNormFactors", "Adding Combined Fit");
            dirs.push_back( fOutDir );
            names.push_back( fName );
            titles.push_back(fCombiLabel);
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

    double xmin = -1.;
    double xmax = 10.;
    double max = 0.;

    // create a list of Norm Factors
    std::vector< string > Names;
    std::vector< string > Titles;
    std::vector< string > Categories;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        for(const auto& inorm : fFitList[i_fit]->fNormFactors) {
            const std::string normName = inorm->fName;
            if(normName==fPOI) continue;
            if(Common::FindInStringVector(Names,normName)<0){
                Names.push_back(normName);
                Titles.push_back(inorm->fTitle);
                Categories.push_back(inorm->fCategory);
            }
        }
    }
    unsigned int Nnorm = Names.size();

    // read fit resutls
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        fFitList[i_fit]->ReadFitResults(dirs[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".txt");
    }

    // exclude norm factors
    std::vector<string> NamesNew;
    std::vector<string> TitlesNew;
    std::vector<string> CategoriesNew;
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            const FitResults* fitRes = fFitList[i_fit]->fFitResults;
            for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
                const NuisParameter* par = fitRes->fNuisPar[j].get();
                const std::string& normName = par->fName;
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
        std::string normName;
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
    std::vector< TGraphAsymmErrors > g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,double> centralMap;
        std::map<string,double> errUpMap;
        std::map<string,double> errDownMap;
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            std::vector<std::string> s;
            if(fFitResultsFile!="") fitRes->ReadFromTXT(fFitResultsFile, s);
            else                    fitRes->ReadFromTXT(fOutDir+"/Fits/"+fName+fSaveSuf+".txt", s);
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults;
        }
        for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
            NuisParameter* par = fitRes->fNuisPar[j].get();
            const std::string& normName = par->fName;
            centralMap[normName] = par->fFitValue;
            errUpMap[normName]   = par->fPostFitUp;
            errDownMap[normName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.emplace_back(Nnorm);
        for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
            const std::string& normName = Names[i_norm];
            if(centralMap[normName]!=0 || (errUpMap[normName]!=0 || errDownMap[normName]!=0)){
                g[i_fit].SetPoint(i_norm,centralMap[normName],(Nnorm-i_norm-1)+0.5+yshift[i_fit]);
                g[i_fit].SetPointEXhigh(i_norm,  errUpMap[normName]);
                g[i_fit].SetPointEXlow( i_norm, -errDownMap[normName]);
            }
            else{
                g[i_fit].SetPoint(i_norm,-10,-10);
                g[i_fit].SetPointEXhigh(i_norm, 0);
                g[i_fit].SetPointEXlow( i_norm, 0);
            }
        }
    }

    max = Nnorm;

    int lineHeight = 50;
    int offsetUp = 50;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.2/(8./6.));
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l1(1,0,1,max);
    l1.SetLineColor(kGray);
    l1.Draw("same");

    int i_style = 0;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(!(fCombine && i_fit==N-1)){
            if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        g[i_fit].SetLineColor(color[i_style]);
        g[i_fit].SetMarkerColor(color[i_style]);
        g[i_fit].SetMarkerStyle(style[i_style]);
        g[i_fit].Draw("P same");
        i_style++;
    }

    TLatex norms{};
    norms.SetTextSize( norms.GetTextSize()*0.8 );
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        norms.DrawLatex(xmin-0.15*(xmax-xmin),(Nnorm-i_norm-1)+0.25,Titles[i_norm].c_str());
    }
    TLatex values{};
    values.SetTextSize( values.GetTextSize()*0.8 );
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(!(fCombine && i_fit==N-1)){
                if (!fFitList.at(i_fit)->fUseInComparison) continue;
            }
            values.DrawLatex((xmin+(xmax-xmin)*(0.45+i_fit*0.55/N)),(Nnorm-i_norm-1)+0.25,
                             Form("%.2f^{+%.2f}_{-%.2f}",g[i_fit].GetX()[i_norm],g[i_fit].GetErrorXhigh(i_norm),g[i_fit].GetErrorXlow(i_norm)));
        }
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );

    TLegend leg(0.45,1.-0.02*(30./max),0.98,0.99);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(!(fCombine && i_fit==N-1)){
         if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        leg.AddEntry(&g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg.Draw();

    gPad->RedrawAxis();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        if(category=="") c.SaveAs((fOutDir+"/NormFactors_comp"+fSaveSuf+"."+format).c_str());
        else             c.SaveAs((fOutDir+"/NormFactors_comp"+fSaveSuf+"_"+category+"."+format).c_str());
    }
}

//__________________________________________________________________________________
//
void MultiFit::PlotCombinedCorrelationMatrix() const{
    TRExFit *fit = fFitList[0];
    if(fit->fStatOnly){
        WriteInfoStatus("MultiFit::PlotCombinedCorrelationMatrix", "Stat only fit => No Correlation Matrix generated.");
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TRExFitter::CORRELATIONTHRESHOLD)
    fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
    if(fit->fFitResults){
        fit->fFitResults->fOutFolder = fOutDir;
        std::vector<std::string> formats;
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            formats.emplace_back(fOutDir+"/CorrMatrix_comb"+fSaveSuf+"."+format);
        }
        fit->fFitResults->DrawCorrelationMatrix(formats,fuseGammasForCorr, fHEPDataFormat, TRExFitter::CORRELATIONTHRESHOLD);
    }
}

//____________________________________________________________________________________
//
void MultiFit::ProduceNPRanking(const  std::string& NPnames) const {
    WriteInfoStatus("MultiFit::ProduceNPRanking", "....................................");
    WriteInfoStatus("MultiFit::ProduceNPRanking", "Producing Ranking...");

    if(fFitType==2){
        WriteErrorStatus("MultiFit::ProduceNPRanking", "For ranking plots, the SPLUSB FitType is needed.");
        exit(EXIT_FAILURE);
    }

    const std::string& inputData = fDataName;

    // create a list of norm factors
    std::vector< std::shared_ptr<Systematic> > vSystematics = GetFitSystematics();
    std::vector< std::shared_ptr<NormFactor> > vNormFactors = GetFitNormFactors();
    
    //
    // List of systematics to check
    //
    RankingManager manager{};
    
    std::vector<string> systNames_unique;
    for(const auto& isyst : vSystematics) {
        if((NPnames=="all") || (NPnames == isyst->fNuisanceParameter)) {
            if(isyst->fType == Systematic::SHAPE) continue;
            if (std::find(systNames_unique.begin(), systNames_unique.end(),
                isyst->fNuisanceParameter) == systNames_unique.end()){
                systNames_unique.emplace_back(isyst->fNuisanceParameter);
            }
            else {
                continue;
            }
            manager.AddNuisPar(isyst->fNuisanceParameter, false);
        }
    }
    for(const auto& inorm: vNormFactors) {
        if(inorm->fName == fPOI) continue;
        if(NPnames=="all" || NPnames == inorm->fName){
            manager.AddNuisPar(inorm->fName, true);
        }
    }

    //
    // Text files containing information necessary for drawing of ranking plot
    //
    std::string outName = fOutDir+"/Fits/NPRanking";
    if(NPnames!="all") outName += "_"+NPnames;
    outName += ".txt";
    
    manager.SetOutputPath(outName);

    //
    // Get the combined model
    //
    std::unique_ptr<TFile> f(TFile::Open((fOutDir+"/ws_combined"+fSaveSuf+".root").c_str()));
    if (!f) {
        WriteErrorStatus("MultiFit::ProduceNPRanking", "Cannot open file!");
        exit(EXIT_FAILURE);
    }
    RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f->Get("combWS"));
    if (!ws) {
        WriteErrorStatus("MultiFit::ProduceNPRanking", "Cannot read workspace!");
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    if(inputData=="asimovData"){
        RooArgSet empty;// = RooArgSet();
        data = static_cast<RooDataSet*>(RooStats::AsymptoticCalculator::MakeAsimovData( (*mc), RooArgSet(ws->allVars()), (RooArgSet&)empty));
    }
    else if(inputData!=""){
        data = static_cast<RooDataSet*>(ws->data(inputData.c_str()));
    } else {
        WriteWarningStatus("MultiFit::ProduceNPRanking", "You didn't specify inputData => will try with observed data !");
        data = static_cast<RooDataSet*>(ws->data("obsData"));
        if(!data){
            WriteWarningStatus("MultiFit::ProduceNPRanking", "Observed data not present => will use with asimov data !");
            data = static_cast<RooDataSet*>(ws->data("asimovData"));
        }
    }

    // Loop on NPs to find gammas and add to the list to be ranked
    if(NPnames=="all" || NPnames.find("gamma")!=string::npos){
        RooRealVar* var = nullptr;
        const RooArgSet* nuis = static_cast<const RooArgSet*>(mc->GetNuisanceParameters());
        if(nuis){
            std::unique_ptr<TIterator> it2(nuis->createIterator());
            while( (var = static_cast<RooRealVar*>(it2->Next()))) {
                const std::string& np = var->GetName();
                if (np.find("saturated_model") != std::string::npos) continue;
                if(np.find("gamma")!=string::npos){
                    // add the nuisance parameter to the list nuisPars if it's there in the ws
                    // remove "gamma"...
                    if(np==NPnames || NPnames=="all"){
                        manager.AddNuisPar(Common::ReplaceString(np,"gamma_",""), true);
                        if(NPnames!="all") break;
                    }
                }
            }
        }
    }

    //
    // Create snapshot to keep inital values
    //
    ws -> saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));

    // find the last config that is used in a fit
    int validConfig(999999);
    for (int i = fFitList.size() - 1; i >= 0; --i) {
        if (fFitList.at(i)->fUseInFit) {
            validConfig = i;
            break;
        }
    }

    if (validConfig == 999999) {
        WriteErrorStatus("ProduceNPRanking", "Did not find a single config that should be used for a fit");
        exit(EXIT_FAILURE);
    }

    TRExFit *fit = fFitList[validConfig];
    fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
    
    manager.SetInjectGlobalObservables(fit->fInjectGlobalObservables);
    manager.SetNPValues(fit->fFitNPValues);
    manager.SetFixedNPs(GetFixedNPs());
    manager.SetFitStrategy(fFitStrategy);
    manager.SetPOIName(fPOI);
    manager.SetNCPU(fCPU);
    manager.SetRng(fit->fRndRange, fit->fUseRnd, fit->fRndSeed);
    manager.SetStatOnly(fStatOnly);

    manager.RunRanking(fit->fFitResults, ws, data, GetFitNormFactors());
    
    f->Close();
}

//____________________________________________________________________________________
//
void MultiFit::PlotNPRankingManager() const{
  if(fFitList[0]->fRankingPlot=="Merge"  || fFitList[0]->fRankingPlot=="all") PlotNPRanking(true,true);
  if(fFitList[0]->fRankingPlot=="Systs"  || fFitList[0]->fRankingPlot=="all") PlotNPRanking(true,false);
  if(fFitList[0]->fRankingPlot=="Gammas" || fFitList[0]->fRankingPlot=="all") PlotNPRanking(false,true);
}

//____________________________________________________________________________________
//
void MultiFit::PlotNPRanking(bool flagSysts, bool flagGammas) const {
    WriteInfoStatus("MultiFit::PlotNPRanking", "....................................");
    WriteInfoStatus("MultiFit::PlotNPRanking", "Plotting Ranking...");
    //
    const std::string fileToRead = fOutDir+"/Fits/NPRanking"+fSaveSuf+".txt";
    std::ifstream in(fileToRead.c_str());
    if (!in.good()) { // file doesnt exist
        const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fOutDir+"/Fits/","NPRanking" + fSaveSuf + "_");
        Common::MergeTxTFiles(inPaths, fileToRead);
    }

    RankingManager manager{};
    manager.SetOutputPath(fileToRead);
    manager.SetAtlasLabel(fFitList[0]->fAtlasLabel);
    manager.SetLumiLabel(fFitList[0]->fLumiLabel);
    manager.SetCmeLabel(fFitList[0]->fCmeLabel);
    manager.SetUseHEPDataFormat(fHEPDataFormat);
    manager.SetName(fOutDir);
    manager.SetSuffix(fSaveSuf);
    manager.SetMaxNPPlot(fFitList[0]->fRankingMaxNP);
    manager.SetRankingPOIName(fFitList[0]->fRankingPOIName);
    manager.SetRankingCanvasSize(fFitList[0]->fNPRankingCanvasSize);

    manager.PlotRanking(GetFitRegions(), flagSysts, flagGammas);
}

//__________________________________________________________________________________
//
void MultiFit::GetLikelihoodScan( RooWorkspace *ws, const std::string& varName, RooDataSet* data,bool recreate) const{
    WriteInfoStatus("MultiFit::GetLikelihoodScan", "Running likelihood scan for the parameter = " + varName);
    static const TString LHDir("LHoodPlots/");

    double minVal(-3);
    double maxVal(3);

    std::unique_ptr<TGraph> graph;
    if(recreate){
        LikelihoodScanManager manager{};
        manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps);
        manager.SetNCPU(fCPU);
        manager.SetOffSet(true);
        
        auto scanResult = manager.Run1DScan(ws, varName, data);

        const std::vector<double> x = scanResult.first;
        const std::vector<double> y = scanResult.second;

        if (x.empty() || y.empty()) {
            WriteWarningStatus("MultiFit::GetLikelihoodScan", "Skipping LHscan");
            return;
        }

        graph = std::make_unique<TGraph>(fLHscanSteps, &x[0], &y[0]);
        minVal = x[0];
        maxVal = x.back();

    }
    else{
        TFile* f = TFile::Open(fName+"/"+LHDir+"NLLscan_"+varName+"_curve.root", "READ");
        graph.reset(dynamic_cast<TGraph*>(f->Get("LHscan")));
        graph->SetLineColor(kRed);
        graph->SetLineWidth(3);
        minVal = graph->GetXaxis()->GetXmin();
        maxVal = graph->GetXaxis()->GetXmax();
    }

    TCanvas can("NLLscan");
    can.SetTopMargin(0.1);
    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(varName);

    can.SetTitle(cname);
    can.SetName(cname);
    can.cd();
    graph->Draw("ALP");

    std::size_t nValid(0);
    for (const auto& ifit : fFitList) {
        if (ifit->fUseInFit) ++nValid;
    }

    // take the LH curves also for other fits
    std::vector<TGraph*> curve_fit;
    TLegend leg(0.5,0.85-0.06*(nValid+1),0.75,0.85);
    leg.SetFillColor(kWhite);
    leg.SetBorderSize(0);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    if(fCompare){
        for(auto fit : fFitList){
            if (!fit->fUseInFit) continue;
            TFile* f = nullptr;
            if(fit->fFitResultsFile!=""){
                std::vector<std::string> v = Vectorize(fit->fFitResultsFile,'/');
                f = TFile::Open(v[0]+"/"+LHDir+"NLLscan_"+varName+"_curve.root");
            }
            else{
                f = TFile::Open(fit->fName+"/"+LHDir+"NLLscan_"+varName+"_curve.root");
            }
            if(f!=nullptr) {
                curve_fit.push_back(dynamic_cast<TGraph*>(f->Get("LHscan")));
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
            crv->Draw("LP same");
            leg.AddEntry(crv,fFitList[idx]->fLabel.c_str(),"l");
            idx++;
        }
        leg.AddEntry(graph.get(),"Combined","l");
    }

    graph->GetXaxis()->SetRangeUser(minVal,maxVal);

    // y axis
    graph->GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") graph->GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") graph->GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());


    TLatex tex{};
    tex.SetTextColor(kGray+2);

    TLine l1s(minVal,0.5,maxVal,0.5);
    l1s.SetLineStyle(kDashed);
    l1s.SetLineColor(kGray);
    l1s.SetLineWidth(2);
    if(graph->GetMaximum()>2){
        l1s.Draw();
        tex.DrawLatex(maxVal,0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }

    if(graph->GetMaximum()>2){
        TLine l2s(minVal,2,maxVal,2);
        l2s.SetLineStyle(kDashed);
        l2s.SetLineColor(kGray);
        l2s.SetLineWidth(2);
        l2s.Draw();
        tex.DrawLatex(maxVal,2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
    }
    //
    if(graph->GetMaximum()>4.5){
        TLine l3s(minVal,4.5,maxVal,4.5);
        l3s.SetLineStyle(kDashed);
        l3s.SetLineColor(kGray);
        l3s.SetLineWidth(2);
        l3s.Draw();
        tex.DrawLatex(maxVal,4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
    }
    //
    TLine lv0(0,graph->GetMinimum(),0,graph->GetMaximum());
    lv0.Draw();
    //
    TLine lh0(minVal,0,maxVal,0);
    lh0.Draw();

    system(TString("mkdir -vp ")+fName+"/"+LHDir);

    if(fCompare){
        leg.Draw();
        if (fFitList[0]->fAtlasLabel != "none") ATLASLabel(0.15,0.93,fFitList[0]->fAtlasLabel.c_str(),kBlack);
        myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
        if(fLabel!="") myText(0.2,0.85,kBlack,Form("#kern[-1]{%s}",fLabel.c_str()));
    }

    can.RedrawAxis();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        can.SaveAs( fName+"/"+LHDir+"NLLscan_"+varName+"."+format);
    }

    if(recreate){
        // write it to a ROOT file as well
        std::unique_ptr<TFile> f(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varName+"_curve.root","UPDATE"));
        f->cd();
        graph->Write("LHscan",TObject::kOverwrite);
        f->Close();
    }
}

//____________________________________________________________________________________
//
void MultiFit::Get2DLikelihoodScan( RooWorkspace *ws, const std::vector<std::string>& varNames, RooDataSet* data) const{
    if (varNames.size() != 2){
        WriteErrorStatus("MultiFit::Get2DLikelihoodScan", "Wrong number of parameters provided for 2D likelihood scan, returning");
        return;
    }

    LikelihoodScanManager manager{};
    manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps);
    manager.SetScanParamsY(fLHscanMinY, fLHscanMaxY, fLHscanStepsY);
    manager.SetNCPU(fCPU);
    manager.SetOffSet(!fParal2D);
    if (fParal2D) {
        manager.SetParalel2Dstep(fParal2Dstep);
    }

    const std::pair<std::string, std::string> names = std::make_pair(varNames.at(0), varNames.at(1));
    auto scanResult = manager.Run2DScan(ws, names, data);

    const std::vector<double> x = scanResult.x;
    const std::vector<double> y = scanResult.y;
    std::vector<std::vector<double> > z = scanResult.z;

    if (x.empty() || y.empty() || z.empty()) {
        WriteWarningStatus("MultiFit::Get2DLikelihoodScan", "Skipping LHscan");
        return;
    }

    // find minimium in z
    auto findMin = [](const std::vector<std::vector<double> >& vec) {
        double min(9999999);
        for (const auto& i : vec) {
            for (const auto& j : i) {
                if (j < min) min = j;
            }
        }
        return min;
    };

    const double zmin = findMin(z);

    // make plots
    TCanvas can("2D_NLLscan");
    can.cd();

    TGraph2D graph(fLHscanSteps * fLHscanStepsY);

    TH2D h_nll("NLL", "NLL", fLHscanSteps, x[0], x.back(), fLHscanStepsY, y[0], y.back());
    unsigned int i=0;
    for (int ipoint = 0; ipoint < fLHscanSteps; ++ipoint) {
        if (fParal2D && ipoint!=fParal2Dstep) continue;
        for (int jpoint = 0; jpoint < fLHscanStepsY; ++jpoint) {
            if (!fParal2D) { // if you are paralellizing, no knowledge of the absolute minimum in each job
                // shift the likelihood values to zero
                z[ipoint][jpoint] -= zmin;
            }
            h_nll.SetBinContent(ipoint+1, jpoint+1, z[ipoint][jpoint]);
            i = ipoint * fLHscanStepsY + jpoint;
            graph.SetPoint(i,x[ipoint],y[jpoint],z[ipoint][jpoint]);
        }
    }

    TString LHDir("LHoodPlots/");
    system(TString("mkdir -vp ")+fName+"/"+LHDir);

    if (!fParal2D) { // Only draw and save graph when not running parallel
        gStyle->SetPalette(57); // Reset Palette to default (Pruning or Correlation matrinx changes this)
        graph.Draw("colz");
        graph.GetXaxis()->SetRangeUser(x[0],x.back());
        graph.GetYaxis()->SetRangeUser(y[0],y.back());

        // y axis
        graph.GetXaxis()->SetTitle(varNames.at(0).c_str());
        graph.GetYaxis()->SetTitle(varNames.at(1).c_str());

        // Print the canvas
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            can.SaveAs(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"."+format);
        }

        // write it to a ROOT file as well
        std::unique_ptr<TFile> f(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"_curve.root","UPDATE"));
        f->cd();
        graph.Write(("LHscan_2D_"+varNames.at(0)+"_"+varNames.at(1)).c_str(),TObject::kOverwrite);
        f->Close();
    }

    // Write histogram to Root file as well
    if (fParal2D) {
        std::ostringstream step_os;
        step_os << fParal2Dstep;
        std::string paral2Dstep_str=step_os.str();
        std::unique_ptr<TFile> f2(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"_step"+paral2Dstep_str+"_histo.root","UPDATE"));
        h_nll.Write("NLL",TObject::kOverwrite);
        f2->Close();
    } else {
        std::unique_ptr<TFile> f2(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"_histo.root","UPDATE"));
        h_nll.Write("NLL",TObject::kOverwrite);
        f2->Close();
    }
}

//____________________________________________________________________________________
//
void MultiFit::PlotSummarySoverB() const {
    WriteInfoStatus("MultiFit::PlotSummarySoverB", "....................................");
    WriteInfoStatus("MultiFit::PlotSummarySoverB", "Producing S/B plot...");

    bool includeBonly = false;
    if(fBonlySuffix!="") includeBonly = true;

    fFitList[0]->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".txt");
    double muFit = fFitList[0]->fFitResults->GetNuisParValue(fPOI);
    if (Common::HistFromFile( fOutDir+"/Limits/"+fName+fSaveSuf+".root/limit" ) == nullptr) {
        WriteWarningStatus("MultiFit::PlotSummarySoverB", "Histo pointer is nullptr, skipping plotting.");
        return;
    }
    double muLimit = Common::HistFromFile( fOutDir+"/Limits/"+fName+fSaveSuf+".root/limit" )->GetBinContent(1);

    std::vector<string> fileNames;
    std::vector<string> fileNamesBonly;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(unsigned int i_reg=0;i_reg<ifit->fRegions.size();i_reg++){
            if(ifit->fRegions[i_reg]->fRegionType==Region::VALIDATION) continue;
            fileNames.push_back(ifit->fName+"/Histograms/"+ifit->fRegions[i_reg]->fName+"_postFit.root");
            if(includeBonly)
                fileNamesBonly.push_back(ifit->fName+"/Histograms/"+ifit->fRegions[i_reg]->fName+fBonlySuffix+"_postFit.root");
        }
    }
    int Nhist = (int)fileNames.size();

    //
    // create a list of all the samples
    std::vector<string> sigList;
    std::vector<string> bkgList;
    std::vector<string> dataList;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(unsigned int i_smp=0;i_smp<ifit->fSamples.size();i_smp++){
            if(ifit->fSamples[i_smp]->fType==Sample::SIGNAL && Common::FindInStringVector(sigList,ifit->fSamples[i_smp]->fName)<0)
                sigList.push_back(ifit->fSamples[i_smp]->fName);
            if(ifit->fSamples[i_smp]->fType==Sample::BACKGROUND && Common::FindInStringVector(bkgList,ifit->fSamples[i_smp]->fName)<0)
                bkgList.push_back(ifit->fSamples[i_smp]->fName);
            if(ifit->fSamples[i_smp]->fType==Sample::DATA && Common::FindInStringVector(dataList,ifit->fSamples[i_smp]->fName)<0)
                dataList.push_back(ifit->fSamples[i_smp]->fName);
        }
    }
    //
    // create a list of all the systematics
    std::vector<string> systList; systList.clear();
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        // actual systematics
        for(unsigned int i_syst=0;i_syst<ifit->fSystematics.size();i_syst++){
            if(Common::FindInStringVector(systList,ifit->fSystematics[i_syst]->fName)<0)
                systList.push_back(ifit->fSystematics[i_syst]->fName);
        }
        // norm factors
        for(unsigned int i_norm=0;i_norm<ifit->fNormFactors.size();i_norm++){
            if(Common::FindInStringVector(systList,ifit->fNormFactors[i_norm]->fName)<0)
                systList.push_back(ifit->fNormFactors[i_norm]->fName);
        }
    }

    unsigned int Nsyst = systList.size();

    std::vector<std::unique_ptr<TFile> > file;
    std::vector<std::unique_ptr<TFile> > fileBonly;
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
        file.emplace_back(std::move(TFile::Open(fileNames[i_hist].c_str())));
        if(includeBonly){
            WriteDebugStatus("MultiFit::PlotSummarySoverB", "Opening file " + fileNamesBonly[i_hist]);
            fileBonly.emplace_back(std::move(TFile::Open(fileNamesBonly[i_hist].c_str())));
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
            h_tmp = dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+sigList[i_sig]+"_postFit").c_str()));
            if(h_tmp!=nullptr){
                WriteDebugStatus("MultiFit::PlotSummarySoverB", " ... FOUND");
                if(h_sig[i_hist]==nullptr) h_sig[i_hist] = h_tmp;
                else                   h_sig[i_hist]->Add(h_tmp);
            }
        }
        for(unsigned int i_bkg=0;i_bkg<bkgList.size();i_bkg++){
            WriteDebugStatus("MultiFit::PlotSummarySoverB", "  Getting histogram h_"+bkgList[i_bkg]+"_postFit");
            h_tmp = dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_postFit").c_str()));
            if(includeBonly) h_tmpBonly = dynamic_cast<TH1D*>(fileBonly[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_postFit").c_str()));
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
                h_tmp = dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_"+systList[i_syst]+"_Up_postFit").c_str()));
                if(h_tmp!=nullptr){
                    if(h_syst_up[i_syst][i_hist]==nullptr)   h_syst_up[i_syst][i_hist] = h_tmp;
                    else                                 h_syst_up[i_syst][i_hist]->Add(h_tmp);
                }
                // down
                if(h_tmp!=nullptr){
                    h_tmp = dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_"+systList[i_syst]+"_Down_postFit").c_str()));
                    if(h_syst_down[i_syst][i_hist]==nullptr) h_syst_down[i_syst][i_hist] = h_tmp;
                    else                                 h_syst_down[i_syst][i_hist]->Add(h_tmp);
                }
            }
        }
        for(unsigned int i_data=0;i_data<dataList.size();i_data++){
            h_tmp = dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+dataList[i_data]).c_str()));
            if(h_tmp!=nullptr){
                if(h_data[i_hist]==nullptr) h_data[i_hist] = h_tmp;
                else                    h_data[i_hist]->Add(h_tmp);
            }
        }

        if(TRExFitter::PREFITONPOSTFIT) {
            h_tot_bkg_prefit[i_hist] = dynamic_cast<TH1D*>(file[i_hist]->Get("h_tot_bkg_prefit"));
        }

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

    for (auto& ifile : file) {
        ifile->Close();
    }
    for (auto& ifile : fileBonly) {
        ifile->Close();
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

    std::vector<double> SoverSqrtB;

    for(int i_bin=1;i_bin<=h_bkg_comb->GetNbinsX();i_bin++){
        const double sig = h_sig_comb->GetBinContent(i_bin);
        const double bkg = h_bkg_comb->GetBinContent(i_bin);
        SoverSqrtB.emplace_back(sig/bkg);
    }

    std::unique_ptr<TH1D> h_bkg_ord(Rebin(h_bkg_comb,SoverSqrtB,false));
    std::unique_ptr<TH1D> h_bkgBonly_ord(nullptr);
    if(includeBonly) h_bkgBonly_ord = std::unique_ptr<TH1D> (Rebin(h_bkgBonly_comb,SoverSqrtB,false));
    std::unique_ptr<TH1D> h_sig_ord(Rebin(h_sig_comb,SoverSqrtB,false));
    std::unique_ptr<TH1D> h_data_ord(Rebin(h_data_comb,SoverSqrtB));
    std::unique_ptr<TH1D> h_tot_bkg_prefit_ord(nullptr);
    if(TRExFitter::PREFITONPOSTFIT) h_tot_bkg_prefit_ord = std::unique_ptr<TH1D>(Rebin(h_tot_bkg_prefit_comb,SoverSqrtB,false));

    std::vector<std::unique_ptr<TH1D> > h_syst_up_ord  (Nsyst);
    std::vector<std::unique_ptr<TH1D> > h_syst_down_ord(Nsyst);
    for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
        h_syst_up_ord  [i_syst] = std::move(std::unique_ptr<TH1D>(Rebin(static_cast<TH1D*>(h_syst_up_comb  [i_syst]),SoverSqrtB,false)));
        h_syst_down_ord[i_syst] = std::move(std::unique_ptr<TH1D>(Rebin(static_cast<TH1D*>(h_syst_down_comb[i_syst]),SoverSqrtB,false)));
    }

    for(int i_bin=0;i_bin<h_bkg_ord->GetNbinsX()+2;i_bin++){
        double err_tot = h_bkg_ord->GetBinError(i_bin); // this should be the stat unc
        double errUp(0);
        double errDown(0);
        double err(0);
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            for(unsigned int j_syst=0;j_syst<systList.size();j_syst++){
                const double corr = fFitList[0]->fFitResults->fCorrMatrix->GetCorrelation( systList[i_syst],systList[j_syst] );
                errUp   += corr * h_syst_up_ord  [i_syst]->GetBinContent(i_bin) * h_syst_up_ord  [j_syst]->GetBinContent(i_bin);
                errDown += corr * h_syst_down_ord[i_syst]->GetBinContent(i_bin) * h_syst_down_ord[j_syst]->GetBinContent(i_bin);
            }
        }
        errUp   = std::sqrt(errUp);
        errDown = std::sqrt(errDown);
        err = std::abs(errUp+errDown)/2.;
        err_tot = std::hypot(err_tot, err);
        h_bkg_ord->SetBinError(i_bin,err_tot);
    }

    TCanvas c("c","c",600,600);

    TPad pad0("pad0","pad0",0,0.28,1,1,0,0,0);
    pad0.SetTicks(1,1);
    pad0.SetTopMargin(0.05);
    pad0.SetBottomMargin(0);
    pad0.SetLeftMargin(0.14);
    pad0.SetRightMargin(0.05);
    pad0.SetFrameBorderMode(0);
    //
    TPad pad1("pad1","pad1",0,0,1,0.28,0,0,0);
    pad1.SetTicks(1,1);
    pad1.SetTopMargin(0.0);
    pad1.SetBottomMargin(0.37);
    pad1.SetLeftMargin(0.14);
    pad1.SetRightMargin(0.05);
    pad1.SetFrameBorderMode(0);

    pad1.Draw();
    pad0.Draw();
    pad0.cd();

    h_sig_ord->SetLineColor(kRed);
    h_sig_ord->SetFillColor(kRed);

    std::unique_ptr<TH1D> h_sig_ord_lim(static_cast<TH1D*>(h_sig_ord->Clone("h_sig_ord_lim")));
    h_sig_ord_lim->Scale(muLimit/muFit);
    h_sig_ord_lim->SetFillColor(kOrange);
    h_sig_ord_lim->SetLineColor(kOrange);
    std::unique_ptr<TH1D> h_sig_ord_lim_diff(static_cast<TH1D*>(h_sig_ord_lim->Clone("h_sig_ord_lim_diff")));
    h_sig_ord_lim_diff->Add(h_sig_ord.get(),-1);

    THStack h_s{};
    h_s.Add(h_bkg_ord.get());
    h_s.Add(h_sig_ord.get());
    h_s.Add(h_sig_ord_lim_diff.get());
    h_data_ord->Draw("EX0");

    h_s.Draw("HISTsame ][");
    TH1D* h_err = static_cast<TH1D*>(h_bkg_ord->Clone("h_err"));
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

    std::unique_ptr<TLegend> leg(nullptr);
    if(includeBonly) leg = std::make_unique<TLegend>(0.6,0.50,0.90,0.92);
    else             leg = std::make_unique<TLegend>(0.6,0.57,0.90,0.92);
    leg->SetFillStyle(0);
    leg->SetMargin(0.2);
    leg->SetBorderSize(0);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->AddEntry(h_data_ord.get(),"Data","lep");
    leg->AddEntry(h_sig_ord_lim.get(),Form(("%s ("+fPOIName+"_{95%% excl.}=%.1f)").c_str(),fSignalTitle.c_str(),muLimit),"f");
    leg->AddEntry(h_sig_ord.get(),    Form(("%s ("+fPOIName+"_{fit}=%.1f)"       ).c_str(),fSignalTitle.c_str(),muFit),"f");
    leg->AddEntry(h_bkg_ord.get(),"Background","f");
    leg->AddEntry(h_err,"Bkgd. Unc.","f");
    if(includeBonly) leg->AddEntry(h_bkgBonly_ord.get(),("Bkgd. ("+fPOIName+"=0)").c_str(),"l");
    if(TRExFitter::PREFITONPOSTFIT) leg->AddEntry(h_tot_bkg_prefit_comb,"Pre-Fit Bkgd.","l");
    leg->Draw();

    if (fFitList[0]->fAtlasLabel != "none") ATLASLabelNew(0.17,0.87, (char*)fFitList[0]->fAtlasLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(0.17,0.80,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()) );
    if(fLabel!="") myText(0.17,0.18,kBlack,Form("%s Combined",fLabel.c_str()) );
    else           myText(0.17,0.18,kBlack,"Combined");
    std::string channels = "";
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        if (!fFitList.at(i_fit)->fUseInFit) continue;
        if(i_fit!=0){
            if(i_fit==fFitList.size()-1) channels += " and ";
            else                         channels += ", ";
        }
        channels += fFitList[i_fit]->fLabel;
    }
    myText(0.17,0.13,kBlack,channels.c_str());
    myText(0.17,0.05,kBlack,"Post-Fit");

    pad0.RedrawAxis();
    pad0.SetLogy();

    pad1.cd();
    pad1.GetFrame()->SetY1(2);
    std::unique_ptr<TH1D> h_ratio(static_cast<TH1D*>(h_data_ord->Clone("h_ratio")));
    std::unique_ptr<TH1D> h_den(static_cast<TH1D*>(h_bkg_ord->Clone("h_den")));
    for(int i_bin=0;i_bin<h_den->GetNbinsX()+2;i_bin++){
        h_den->SetBinError(i_bin,0);
    }

    std::unique_ptr<TH1D> h_ratioBonly(nullptr);
    if(includeBonly){
        h_ratioBonly = std::unique_ptr<TH1D>(static_cast<TH1D*>(h_bkgBonly_ord->Clone("h_ratioBonly")));
        h_ratioBonly->Divide(h_den.get());
        h_ratioBonly->SetLineStyle(kDashed);
        h_ratioBonly->SetLineColor(kBlack);
    }

    std::unique_ptr<TH1D> h_stackSig(static_cast<TH1D*>(h_sig_ord ->Clone("h_sig_ratio")));
    h_stackSig->Add(h_bkg_ord.get());
    h_stackSig->Divide(h_den.get());
    h_stackSig->SetFillColor(0);
    h_stackSig->SetFillStyle(0);
    h_stackSig->SetLineColor(kRed);

    std::unique_ptr<TH1D> h_stackSigLim (static_cast<TH1D*>(h_sig_ord_lim ->Clone("h_sig_lim_ratio")));
    h_stackSigLim->Add(h_bkg_ord.get());
    h_stackSigLim->Divide(h_den.get());
    h_stackSigLim->SetFillColor(0);
    h_stackSigLim->SetFillStyle(0);
    h_stackSigLim->SetLineStyle(kDashed);
    h_stackSigLim->SetLineColor(kOrange+1);

    std::unique_ptr<TH1D> h_ratio2(static_cast<TH1D*>(h_err->Clone("h_ratio2")));
    h_ratio2->SetMarkerSize(0);
    h_ratio->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetTitle("Data / Bkgd.");
    h_ratio->GetYaxis()->SetTitleSize(20);
    h_ratio->GetYaxis()->SetTitleOffset(2.);
    h_ratio->GetYaxis()->SetLabelSize(20); // 0.04
    h_ratio ->Divide(h_den.get());
    h_ratio2->Divide(h_den.get());
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
    TLine hline(h_ratio->GetXaxis()->GetXmin(),1,h_ratio->GetXaxis()->GetXmax(),1);
    h_ratio->Draw("E1");
    h_ratio2->Draw("same E2");
    hline.Draw();
    h_stackSig->Draw("same HIST ][");
    h_stackSigLim->Draw("same HIST ][");
    h_ratio->Draw("same E1");

    if(includeBonly) h_ratioBonly->Draw("same HIST ][");

    TLegend leg2(0.17,0.64,0.75,0.98);
    leg2.SetFillStyle(0);
    leg2.SetMargin(0.1);
    leg2.SetBorderSize(0);
    leg2.SetTextSize(gStyle->GetTextSize());
    leg2.AddEntry(h_stackSigLim.get(),Form(("%s ("+fPOIName+"_{95%% excl.}=%.1f) + Bkgd.").c_str(),fSignalTitle.c_str(),muLimit),"l");
    leg2.AddEntry(h_stackSig.get(),   Form(("%s ("+fPOIName+"_{fit}=%.1f) + Bkgd."       ).c_str(),fSignalTitle.c_str(),muFit)  ,"l");
    leg2.Draw();

    std::unique_ptr<TLegend> leg3(nullptr);
    if(includeBonly){
        leg3 = std::make_unique<TLegend>(0.17,0.4,0.75,0.5);
        leg3->SetFillStyle(0);
        leg3->SetMargin(0.1);
        leg3->SetBorderSize(0);
        leg3->SetTextSize(gStyle->GetTextSize());
        leg3->AddEntry(h_ratioBonly.get(),"Bkgd. (from Bkgd-only fit)","l");
        leg3->Draw();
    }


    pad0.RedrawAxis();
    pad1.RedrawAxis();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs( (fOutDir+"/SoverB_postFit."+format).c_str() );
    }
}

//____________________________________________________________________________________
//
TH1D* MultiFit::Combine(vector<TH1D*> h) const{
    int Nbins = 0;
    int Nhist = h.size();
    for(int i_hist=0;i_hist<Nhist;i_hist++){
        if(h[i_hist]==nullptr) WriteWarningStatus("MultiFit::Combine", "empty histogram " + std::to_string(i_hist));
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
// merge bins in bins of SoverSqrtB
TH1D* MultiFit::Rebin(TH1D* h, const vector<double>& vec, bool isData) const{
    TH1D* h_new = new TH1D(Form("%s_rebin",h->GetName()),Form("%s_rebin",h->GetTitle()),17,-3.8,-0.5);
    h_new->Sumw2();
    // new way
    for(int j_bin=1;j_bin<=h->GetNbinsX();j_bin++){
        double value=std::log10(vec[j_bin-1]);
        if ( value<h_new->GetXaxis()->GetXmin() ) value=0.9999*h_new->GetXaxis()->GetXmin();
        if ( value>h_new->GetXaxis()->GetXmax() ) {
            double tmpvalue=1.0001*h_new->GetXaxis()->GetXmax();
            WriteDebugStatus("MultiFit::Rebin", "turning: " + std::to_string(value) + " in: " + std::to_string(tmpvalue));
            value=tmpvalue;
        }
        int i_bin=h_new->FindBin(value);
        h_new->SetBinContent(i_bin,h_new->GetBinContent(i_bin)+h->GetBinContent(j_bin));
        if (!isData) {
            h_new->SetBinError(i_bin,std::hypot(h_new->GetBinError(i_bin), h->GetBinError(j_bin)));
        }
    }
    if (isData) {
        for(int j_bin=1;j_bin<=h_new->GetNbinsX();j_bin++){
            h_new->SetBinError(j_bin, std::sqrt(h_new->GetBinContent(j_bin) ) );
        }
    }
    h_new->SetMinimum(1);
    return h_new;
}

//____________________________________________________________________________________
// combine individual results from grouped impact evaluation into one table
void MultiFit::BuildGroupedImpactTable() const{
    WriteInfoStatus("MultiFit::BuildGroupedImpactTable", "merging grouped impact evaluations");
    const std::string targetName = fOutDir+"/Fits/GroupedImpact"+fSaveSuf+".txt";

    if(std::ifstream(targetName).good()){
        WriteWarningStatus("MultiFit::BuildGroupedImpactTable","file " + targetName + " already exists, will not overwrite");
    }
    else{
        const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fOutDir+"/Fits/","GroupedImpact"+fSaveSuf+"_");
        Common::MergeTxTFiles(inPaths, targetName);
    }
}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<NormFactor> > MultiFit::GetFitNormFactors() const {
    std::vector<std::shared_ptr<NormFactor> > result;
    
    std::vector<std::string> names; 
    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& nf : ifit->fNormFactors) {
            // only add the unique NFs
            if (std::find(names.begin(), names.end(), nf->fName) != names.end()) continue;
            names.emplace_back(nf->fName);

            result.emplace_back(nf);
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::map<std::string, double> MultiFit::GetFixedNPs() const {

    std::map<std::string, double> result;
    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for (const auto& ifixed : ifit->fFitFixedNPs) {
            auto it = result.find(ifixed.first);
            if (it != result.end()) continue;

            result[ifixed.first] = ifixed.second;
        }
    }
    return result;
}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<Systematic> > MultiFit::GetFitSystematics() const {
    
    std::vector<std::shared_ptr<Systematic> > result;

    std::vector< std::string > names;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& isyst : ifit->fSystematics) {
            const std::string systName = isyst->fNuisanceParameter;
            if(Common::FindInStringVector(names,systName) < 0) {
                names.push_back(systName);
                result.emplace_back(isyst);
            }
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<Region* > MultiFit::GetFitRegions() const {
    std::vector<Region*> result;
        
    std::vector< std::string > names;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for (const auto& ireg : ifit->fRegions) {
            const std::string regName = ireg->fName;
            if(Common::FindInStringVector(names,regName) < 0) {
                names.emplace_back(regName);
                result.emplace_back(ireg);
            }
        }
    }

    return result;
}
