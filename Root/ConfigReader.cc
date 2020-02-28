// Class include
#include "TRExFitter/ConfigReader.h"

// Framework inclused
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/TRExFit.h"
#include "TRExFitter/TruthSample.h"
#include "TRExFitter/UnfoldingSample.h"
#include "TRExFitter/UnfoldingSystematic.h"

// ROOT includes
#include "TColor.h"
#include "TSystem.h"

// c++ includes
#include <algorithm>
#include <iostream>

#define _STRINGIZE(x) #x
#define STRINGIZE(x) _STRINGIZE(x)

//__________________________________________________________________________________
//
ConfigReader::ConfigReader(TRExFit *fitter) :
    fFitter(fitter),
    fParser(new ConfigParser()),
    fAllowWrongRegionSample(false),
    fNonGhostIsSet(false),
    fOnlyLHscan(""),
    fHasAtLeastOneValidRegion(false),
    fHasAtLeastOneValidSample(false) {
    WriteInfoStatus("ConfigReader::ConfigReader", "Started reading the config");
}

//__________________________________________________________________________________
// Read the full config file
int ConfigReader::ReadFullConfig(const std::string& fileName, const std::string& opt, const std::string& option){
    // initialize ConfigParser for the actual config
    fParser->ReadFile(fileName);

    // initialize checker COnfigParser to cross check the input
    ConfigParser refConfig;
    std::string homeArea("$TREXFITTER_HOME");
#ifdef TREXFITTER_HOME
    homeArea = std::string(STRINGIZE(TREXFITTER_HOME));
#endif
    refConfig.ReadFile(gSystem->ExpandPathName((homeArea+"/jobSchema.config").c_str()));
    int sc = fParser->CheckSyntax(&refConfig);

    if (sc != 0) return sc;

    // syntax of the config is ok
    // read different types of settings
    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc+= ReadJobOptions();

    sc+= ReadGeneralOptions();

    sc+= ReadFitOptions();

    sc+= ReadLimitOptions();

    sc+= ReadSignificanceOptions();

    sc+= ReadRegionOptions(opt);

    sc+= ReadSampleOptions();

    sc+= ReadNormFactorOptions();

    sc+= ReadShapeFactorOptions();

    sc+= ReadSystOptions();

    sc+= ReadUnfoldingOptions();

    sc+= ReadTruthSamples();

    sc+= ReadUnfoldingSampleOptions();

    sc+= ReadUnfoldingSystematicOptions();

    sc+= UnfoldingCorrections();

    sc+= PostConfig(opt);

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadCommandLineOptions(const std::string& option){
    std::vector< std::string > optVec = Vectorize(option,':');
    std::map< std::string,std::string > optMap;

    for(const std::string& iopt : optVec){
        std::vector< std::string > optPair;
        optPair = Vectorize(iopt,'=');
        if (optPair.size() < 2){
            WriteErrorStatus("ConfigReader::ReadCommandLineOptions", "Cannot read your command line option, please check this!");
            return 1;
        }
        optMap[optPair[0]] = optPair[1];
    }
    if(optMap["Regions"]!=""){
        fOnlyRegions = Vectorize(optMap["Regions"],',');
    }
    if(optMap["Samples"]!=""){
        fOnlySamples = Vectorize(optMap["Samples"],',');
    }
    if(optMap["Systematics"]!=""){
        fOnlySystematics = Vectorize(optMap["Systematics"],',');
    }
    if(optMap["Exclude"]!=""){
        fToExclude = Vectorize(optMap["Exclude"],',');
    }
    if(optMap["Suffix"]!=""){
        fFitter->fSuffix = optMap["Suffix"]; // used for input & output  plots, txt files & workspaces - NOT for histograms file
    }
    if(optMap["SaveSuffix"]!=""){
        fFitter->fSaveSuffix = optMap["SaveSuffix"]; // ... and this one for histograms file
    }
    if(optMap["Update"]!="" && optMap["Update"]!="FALSE"){
        fFitter->fUpdate = true;
    }
    if(optMap["StatOnly"]!="" && optMap["StatOnly"]!="FALSE"){
        fFitter->fStatOnly = true;
    }
    if(optMap["StatOnlyFit"]!="" && optMap["StatOnlyFit"]!="FALSE"){
        fFitter->fStatOnlyFit = true;
    }
    if(optMap["Ranking"]!=""){
        fFitter->fRankingOnly = optMap["Ranking"];
    }
    if(optMap["Signal"]!=""){
        fOnlySignals.clear();
        fOnlySignals.emplace_back(optMap["Signal"]);
    }
    if(optMap["Signals"]!=""){
        fOnlySignals = Vectorize(optMap["Signals"],',');
    }
    if(optMap["FitResults"]!=""){
        fFitter->fFitResultsFile = optMap["FitResults"];
    }
    if(optMap["FitType"]!=""){
        if(optMap["FitType"]=="SPLUSB") fFitter->SetFitType(TRExFit::SPLUSB);
        if(optMap["FitType"]=="BONLY")  fFitter->SetFitType(TRExFit::BONLY);
        if(optMap["FitType"]=="UNFOLDING")  fFitter->SetFitType(TRExFit::UNFOLDING);
    }
    if(optMap["LumiScale"]!=""){
        fFitter->fLumiScale = atof(optMap["LumiScale"].c_str());
    }
    if(optMap["BootstrapIdx"]!=""){
        fFitter->fBootstrapIdx = atoi(optMap["BootstrapIdx"].c_str());
    }
    if(optMap["BootstrapSyst"]!=""){
        fFitter->fBootstrapSyst = optMap["BootstrapSyst"];
    }
    if(optMap["GroupedImpact"]!=""){
        fFitter->fGroupedImpactCategory = optMap["GroupedImpact"];
    }
    if(optMap["OutputDir"]!=""){
        fFitter->fDir = RemoveQuotes(optMap["OutputDir"]);
        if(fFitter->fDir.back() != '/') fFitter->fDir += '/';
        gSystem->mkdir(fFitter->fDir.c_str());
    }
    if(optMap["Job"]!=""){
        fFitter->fName = RemoveQuotes(optMap["Job"]);
    }
    if(optMap["LimitParamValue"]!=""){
        fFitter->fLimitParamValue = atof(optMap["LimitParamValue"].c_str());
    }
    if(optMap["LHscan"]!=""){
        fOnlyLHscan = optMap["LHscan"];
    }
    if(optMap["Parallel2Dscan"]!="" && optMap["Parallel2Dscan"]!="FALSE"){
        fFitter->fParal2D = true;
    }
    if(optMap["Parallel2DscanStep"]!=""){
        fFitter->fParal2Dstep = atoi(optMap["Parallel2DscanStep"].c_str());
    }
    if(optMap["FitBlind"]!=""){
        fFitter->fFitIsBlind = true;
    }
    if(optMap["BlindedParameters"]!=""){
        fFitter->fBlindedParameters = Vectorize(optMap["BlindedParameters"],',');
    }
    //
    WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "-------------------------------------------");
    WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "Running options: ");
    if(fOnlyRegions.size()>0){
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "  Only these Regions: ");
        for(const std::string& ireg : fOnlyRegions){
            WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "    " + ireg);
        }
    }
    if(fOnlySamples.size()>0){
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "  Only these Samples: ");
        for(const std::string& isamp : fOnlySamples){
            WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "    " + isamp);
        }
    }
    if(fOnlySystematics.size()>0){
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "  Only these Systematics: ");
        for(const std::string& isyst : fOnlySystematics){
            WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "    " + isyst);
        }
    }
    if(fToExclude.size()>0){
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "  Exclude: ");
        for(const std::string& iexcl : fToExclude){
            WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "    " + iexcl);
        }
    }
    if(fOnlySignals.size()>0){
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "  Only Signals: ");
        std::string toPrint = "";
        for(const auto& s : fOnlySignals) toPrint += " " + s;
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "    " + toPrint);
    }
    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadJobOptions(){
    std::string param = ""; // helper string

    ConfigSet *confSet = fParser->GetConfigSet("Job");
    if (confSet == nullptr){
        WriteErrorStatus("ConfigReader::ReadJobOptions", "You need to provide JOB settings!");
        return 1;
    }

    if (fFitter->fDir == "") {
        // default
        if (fFitter->fName == "MyMeasurement") fFitter->fName = CheckName(confSet->GetValue());
    } else {
        if (fFitter->fName == "MyMeasurement") fFitter->fName = fFitter->fDir + CheckName(confSet->GetValue());
        else fFitter->fName = fFitter->fDir + fFitter->fName;
    }
    fFitter->fInputName = CheckName(confSet->GetValue());

    //Set DebugLevel
    param = confSet->Get("DebugLevel");
    if( param != "")  TRExFitter::SetDebugLevel( atoi(param.c_str()) );

    param = confSet->Get("AllowWrongRegionSample");
    if( param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") fAllowWrongRegionSample = true;
        else if (param == "FALSE") fAllowWrongRegionSample = false;
    }

    // Set outputDir
    param = confSet->Get("OutputDir");
    if(param != ""){
        if (fFitter->fDir.size() == 0){
            fFitter->fDir = RemoveQuotes(param);
            if(fFitter->fDir.back() != '/') fFitter->fDir += '/';
            fFitter->fName = fFitter->fDir + fFitter->fName;
            gSystem->mkdir((fFitter->fName).c_str(), true);
        }
    }

    // Set Label
    param = confSet->Get("Label");
    if(param!="") fFitter->fLabel = RemoveQuotes(param);
    else          fFitter->fLabel = fFitter->fName;

    // Set POI
    fFitter->SetPOI(CheckName(confSet->Get("POI")));

    // POI unit (if any)
    param = confSet->Get("POIUnit");
    if(param!="") fFitter->fPOIunit = RemoveQuotes(param);

    // Set reading option
    param = confSet->Get("ReadFrom");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if(      param=="HIST" || param=="HISTOGRAMS")  fFitter->fInputType = 0;
    else if( param=="NTUP" || param=="NTUPLES" )    fFitter->fInputType = 1;
    else{
        WriteErrorStatus("ConfigReader::ReadJobOptions", "Invalid \"ReadFrom\" argument. Options: \"HIST\", \"NTUP\"");
        return 1;
    }

    // set default MERGEUNDEROVERFLOW
    if(fFitter->fInputType==0)      TRExFitter::MERGEUNDEROVERFLOW = false;
    else if(fFitter->fInputType==1) TRExFitter::MERGEUNDEROVERFLOW = true;

    // Set MergeUnderOverFlow from config
    param = confSet->Get("MergeUnderOverFlow");
    if(param!=""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  TRExFitter::MERGEUNDEROVERFLOW = true;
        else if( param == "FALSE" ) TRExFitter::MERGEUNDEROVERFLOW = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'MergeUnderOverFlow' option but you did not provide valid setting. Using default (FALSE)");
            TRExFitter::MERGEUNDEROVERFLOW = false;
        }
    }

    // Set paths for the unfolding preprocessing
    param = confSet->Get("ResponseMatrixName");
    if (param != "") {
        fFitter->fResponseMatrixNames.clear();
        fFitter->fResponseMatrixNames.emplace_back(CheckName(param));
    }

    param = confSet->Get("ResponseMatrixNames");
    if (param != "") {
        fFitter->fResponseMatrixNames = Vectorize(param, ',');
    }

    param = confSet->Get("ResponseMatrixFile");
    if (param != "") {
        fFitter->fResponseMatrixFiles.clear();
        fFitter->fResponseMatrixFiles.emplace_back(CheckName(param));
    }

    param = confSet->Get("ResponseMatrixFiles");
    if (param != "") {
        fFitter->fResponseMatrixFiles = Vectorize(param, ',');
    }

    param = confSet->Get("ResponseMatrixPath");
    if (param != "") {
        fFitter->fResponseMatrixPaths.clear();
        fFitter->fResponseMatrixPaths.emplace_back(CheckName(param));
    }

    param = confSet->Get("ResponseMatrixPaths");
    if (param != "") {
        fFitter->fResponseMatrixPaths = Vectorize(param, ',');
    }

    param = confSet->Get("ResponseMatrixNameNominal");
    if(param!=""){
      fFitter->fResponseMatrixNamesNominal.clear();
      fFitter->fResponseMatrixNamesNominal.emplace_back(CheckName(param));
    }

    param = confSet->Get("AcceptanceName");
    if (param != "") {
        fFitter->fAcceptanceNames.clear();
        fFitter->fAcceptanceNames.emplace_back(CheckName(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceNames");
    if (param != "") {
        fFitter->fAcceptanceNames = Vectorize(param, ',');
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceFile");
    if (param != "") {
        fFitter->fAcceptanceFiles.clear();
        fFitter->fAcceptanceFiles.emplace_back(CheckName(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceFiles");
    if (param != "") {
        fFitter->fAcceptanceFiles = Vectorize(param, ',');
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptancePath");
    if (param != "") {
        fFitter->fAcceptancePaths.clear();
        fFitter->fAcceptancePaths.emplace_back(CheckName(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptancePaths");
    if (param != "") {
        fFitter->fAcceptancePaths = Vectorize(param, ',');
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceNameNominal");
    if(param!=""){
        fFitter->fAcceptanceNamesNominal.clear();
        fFitter->fAcceptanceNamesNominal.emplace_back(CheckName(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("SelectionEffName");
    if (param != "") {
        fFitter->fSelectionEffNames.clear();
        fFitter->fSelectionEffNames.emplace_back(CheckName(param));
    }

    param = confSet->Get("SelectionEffNames");
    if (param != "") {
        fFitter->fSelectionEffNames = Vectorize(param, ',');
    }

    param = confSet->Get("SelectionEffFile");
    if (param != "") {
        fFitter->fSelectionEffFiles.clear();
        fFitter->fSelectionEffFiles.emplace_back(CheckName(param));
    }

    param = confSet->Get("SelectionEffFiles");
    if (param != "") {
        fFitter->fSelectionEffFiles = Vectorize(param, ',');
    }

    param = confSet->Get("SelectionEffPath");
    if (param != "") {
        fFitter->fSelectionEffPaths.clear();
        fFitter->fSelectionEffPaths.emplace_back(CheckName(param));
    }

    param = confSet->Get("SelectionEffPaths");
    if (param != "") {
        fFitter->fSelectionEffPaths = Vectorize(param, ',');
    }

    param = confSet->Get("SelectionEffNameNominal");
    if(param!=""){
      fFitter->fSelectionEffNamesNominal.clear();
      fFitter->fSelectionEffNamesNominal.emplace_back(CheckName(param));
    }

    param = confSet->Get("MigrationName");
    if (param != "") {
        fFitter->fMigrationNames.clear();
        fFitter->fMigrationNames.emplace_back(CheckName(param));
    }

    param = confSet->Get("MigrationNames");
    if (param != "") {
        fFitter->fMigrationNames = Vectorize(param, ',');
    }

    param = confSet->Get("MigrationFile");
    if (param != "") {
        fFitter->fMigrationFiles.clear();
        fFitter->fMigrationFiles.emplace_back(CheckName(param));
    }

    param = confSet->Get("MigrationFiles");
    if (param != "") {
        fFitter->fMigrationFiles = Vectorize(param, ',');
    }

    param = confSet->Get("MigrationPath");
    if (param != "") {
        fFitter->fMigrationPaths.clear();
        fFitter->fMigrationPaths.emplace_back(CheckName(param));
    }

    param = confSet->Get("MigrationPaths");
    if (param != "") {
        fFitter->fMigrationPaths = Vectorize(param, ',');
    }

    param = confSet->Get("MigrationNameNominal");
    if(param!=""){
      fFitter->fMigrationNamesNominal.clear();
      fFitter->fMigrationNamesNominal.emplace_back(CheckName(param));
    }

    // Set paths
    // HIST option only
    if(fFitter->fInputType==0){
        if (confSet->Get("MCweight") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'MCweight:' option, ignoring.");
        }
        if (confSet->Get("Selection") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'Selection:' option, ignoring.");
        }
        if (confSet->Get("NtupleName") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtupleName:' option, ignoring.");
        }
        if (confSet->Get("NtupleFile") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtupleFile:' option, ignoring.");
        }
        if (confSet->Get("NtuplePath") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtuplePath:' option, ignoring.");
        }
        if (confSet->Get("NtuplePaths") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtuplePaths:' option, ignoring.");
        }
        //
        param = confSet->Get("HistoName");
        if(param!=""){
            fFitter->fHistoNames.clear();
            fFitter->fHistoNames.push_back( CheckName(param) );
        }
        param = confSet->Get("HistoNames");
        if(param!=""){
            fFitter->fHistoNames = Vectorize( param,',' );
        }
        param = confSet->Get("HistoNameNominal");
        if(param!=""){
          fFitter->fHistoNamesNominal.clear();
          fFitter->fHistoNamesNominal.push_back( CheckName(param) );
        }
        param = confSet->Get("HistoFile");
        if(param!=""){
            fFitter->fHistoFiles.clear();
            fFitter->fHistoFiles.push_back( CheckName(param) );
        }
        param = confSet->Get("HistoFiles");
        if(param!=""){
            fFitter->fHistoFiles = Vectorize( param,',' );
        }
        param = confSet->Get("HistoPath");
        if(param!=""){
            fFitter->AddHistoPath( CheckName(param) );
        }
        param = confSet->Get("HistoPaths");
        if(param!=""){
            fFitter->fHistoPaths = Vectorize( param,',' );
        }
    }
    // Setting for NTUP only
    if(fFitter->fInputType==1){
        if (confSet->Get("HistoName") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NTUP option but you provided 'HistoName:' option, ignoring.");
        }
        if (confSet->Get("HistoFile") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NTUP option but you provided 'HistoFile:' option, ignoring.");
        }
        if (confSet->Get("HistoPath") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NTUP option but you provided 'HistoPath:' option, ignoring.");
        }
        if (confSet->Get("HistoPaths") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NTUP option but you provided 'HistoPaths:' option, ignoring.");
        }
        //
        param = confSet->Get("NtupleName");
        if(param!=""){
            fFitter->SetNtupleName( CheckName(param) );
        }
        param = confSet->Get("NtupleNames");
        if(param!=""){
            fFitter->fNtupleNames = Vectorize( param,',' );
        }
        param = confSet->Get("NtupleFile");
        if( param != "" ){
            fFitter->SetNtupleFile( CheckName(param) );
        }
        param = confSet->Get("NtupleFiles");
        if(param!=""){
            fFitter->fNtupleFiles = Vectorize( param,',' );
        }
        param = confSet->Get("NtuplePath");
        if( param != "" ) {
            fFitter->AddNtuplePath( CheckName(param) );
        }
        param = confSet->Get("NtuplePaths");
        if( param != "" ){
            std::vector<std::string> paths = Vectorize( param,',' );
            for(const std::string& ipath : paths){
                fFitter->AddNtuplePath( ipath );
            }
        }
        param = confSet->Get("MCweight");
        if(param!=""){
            fFitter->SetMCweight( RemoveQuotes(param) );
        }
        param = confSet->Get("Selection");
        if(param!=""){
            fFitter->SetSelection( RemoveQuotes(param) );
        }
    }

    // Set lumi
    param = confSet->Get("Lumi");
    if( param != "" ) fFitter->SetLumi( atof(param.c_str()) );

    // Set LumiScale
    param = confSet->Get("LumiScale");
    if( param != "" ){
        WriteWarningStatus("ConfigReader::ReadJobOptions", "\"LumiScale\" is only done for quick tests since it is inefficient.");
        WriteWarningStatus("ConfigReader::ReadJobOptions", "To normalize all the samples to the luminosity, use \"Lumi\" instead.");
        fFitter->fLumiScale = atof(param.c_str());
    }

    // Set Smoothing option
    param = confSet->Get("SmoothingOption");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "MAXVARIATION" ) fFitter->fSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
        else if (param == "TTBARRESONANCE") fFitter->fSmoothOption = HistoTools::SmoothOption::TTBARRESONANCE;
        else if (param == "COMMONTOOLSMOOTHMONOTONIC") fFitter->fSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHMONOTONIC;
        else if (param == "COMMONTOOLSMOOTHPARABOLIC") fFitter->fSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHPARABOLIC;
        else if (param == "KERNELRATIOUNIFORM") fFitter->fSmoothOption = HistoTools::SmoothOption::KERNELRATIOUNIFORM;
        else if (param == "KERNELDELTAGAUSS") fFitter->fSmoothOption = HistoTools::SmoothOption::KERNELDELTAGAUSS;
        else if (param == "KERNELRATIOGAUSS") fFitter->fSmoothOption = HistoTools::SmoothOption::KERNELRATIOGAUSS;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'SmoothingOption' option but you did not provide valid input. Using default (MAXVARIATION)");
            fFitter->fSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
        }
    }

    // Set SystPruningShape
    param = confSet->Get("SystPruningShape");
    if( param != "") fFitter->fThresholdSystPruning_Shape = atof(param.c_str());

    // Set SystPruningNorm
    param = confSet->Get("SystPruningNorm");
    if( param != "")  fFitter->fThresholdSystPruning_Normalisation = atof(param.c_str());

    // Set SystLarge
    param = confSet->Get("SystLarge");
    if( param != "")  fFitter->fThresholdSystLarge = atof(param.c_str());

    // Set IntCodeOverall
    param = confSet->Get("IntCodeOverall");
    if( param != "")  fFitter->fIntCode_overall = atoi(param.c_str());

    // Set IntCodeShape
    param = confSet->Get("IntCodeShape");
    if( param != "")  fFitter->fIntCode_shape = atoi(param.c_str());

    // Set MCstatThreshold
    param = confSet->Get("MCstatThreshold");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param=="NONE")  fFitter->SetStatErrorConfig( false, 0. );
        else{
            fFitter->SetStatErrorConfig( true, atof(param.c_str()));
        }
    }
    else{
        fFitter->SetStatErrorConfig( true, 0. );
    }

    //Set MCstatConstraint
    param = confSet->Get("MCstatConstraint");
    if( param != "") fFitter->fStatErrCons = RemoveQuotes(param);

    // Set UseGammaPulls
    param = confSet->Get("UseGammaPulls");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") fFitter->fUseGammaPulls = true;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'UseGammaPulls' option but you did not set it to TRUE. Using default (FALSE)");
            fFitter->fUseGammaPulls = false;
        }
    }

    // plotting options are in special function
    if (SetJobPlot(confSet) != 0) return 1;

    // Set TableOptions
    param = confSet->Get("TableOptions");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fFitter->fTableOptions = RemoveQuotes(param);
    }

    // Set CorrelationThreshold
    param = confSet->Get("CorrelationThreshold");
    if( param != ""){
        // set it if it was not set already
        if (TRExFitter::CORRELATIONTHRESHOLD < 0){
            TRExFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
        }
    }

    // Set HistoChecks
    param = confSet->Get("HistoChecks");
    if(param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "NOCRASH" ){
            TRExFitter::HISTOCHECKCRASH = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'HistoChecks' option but you did not set it to NOCRASH.");
        }
    }

    // Set LumiLabel
    param = confSet->Get("LumiLabel");
    if( param != "") fFitter->fLumiLabel = RemoveQuotes(param);

    // Set CmeLabel
    param = confSet->Get("CmeLabel");
    if( param != "") fFitter->fCmeLabel = RemoveQuotes(param);

    // Set SplitHistoFiles
    param = confSet->Get("SplitHistoFiles");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TRExFitter::SPLITHISTOFILES = true;
        } else if (param == "FALSE") {
            TRExFitter::SPLITHISTOFILES = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'SplitHistoFiles' option but you did not provide valid setting. Using default (false)");
            TRExFitter::SPLITHISTOFILES = false;
        }
    }

    // Set BlindingThreshold"
    param = confSet->Get("BlindingThreshold");
    if( param != ""){
        fFitter->fBlindingThreshold = atof(param.c_str());
    }

    // Set BlindingType
    param = confSet->Get("BlindingType");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param=="SOVERB" )                fFitter->fBlindingType = Common::SOVERB;
        else if ( param=="SOVERSPLUSB" )     fFitter->fBlindingType = Common::SOVERSPLUSB;
        else if ( param=="SOVERSQRTB" )      fFitter->fBlindingType = Common::SOVERSQRTB;
        else if ( param=="SOVERSQRTSPLUSB" ) fFitter->fBlindingType = Common::SOVERSQRTSPLUSB;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'BlindingType' option but did not provide a valid setting. Using default (SOVERB)");
            fFitter->fBlindingType = Common::SOVERB;
        }
    }

    // Set RankingMaxNP
    param = confSet->Get("RankingMaxNP");
    if( param != ""){
        fFitter->fRankingMaxNP = atoi(param.c_str());
    }

    // Set ImageFormat
    param = confSet->Get("ImageFormat");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param,',');
        if (tmp.size() > 0) fFitter->fImageFormat = tmp.at(0);
        else {
            WriteErrorStatus("ConfigReader::ReadJobOptions", "You specified 'ImageFormat' option but we cannot split the setting. Please check");
        }
        TRExFitter::IMAGEFORMAT = tmp;
    }

    // Set StatOnly
    param = confSet->Get("StatOnly");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fStatOnly = true;
        } else if (param == "FALSE"){
            fFitter->fStatOnly = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'StatOnly' option but you did not provide valid setting. Using default (false)");
            fFitter->fStatOnly = false;
        }
    }

    // Set FixNPforStatOnly
    param = confSet->Get("FixNPforStatOnly");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fFixNPforStatOnlyFit = true;
        } else if (param == "FALSE"){
            fFitter->fFixNPforStatOnlyFit = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'FixNPforStatOnly' option but you did not provide valid setting. Using default (false)");
            fFitter->fFixNPforStatOnlyFit = false;
        }
    }

    // Set InputFolder
    param = confSet->Get("InputFolder");
    if( param != "" ){
        fFitter->fInputFolder = RemoveQuotes(param);
    }

    // Set InputName
    param = confSet->Get("InputName");
    if( param != "" ){
        fFitter->fInputName = RemoveQuotes(param);
    }

    // Set WorkspaceFileName
    param = confSet->Get("WorkspaceFileName");
    if( param != "" ){
        fFitter->fWorkspaceFileName = RemoveQuotes(param);
    }

    // Set KeepPruning
    param = confSet->Get("KeepPruning");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fKeepPruning = true;
        else if (param == "FALSE") fFitter->fKeepPruning = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'KeepPruning' option but you did not provide valid setting. Using default (false)");
            fFitter->fKeepPruning = false;
        }
    }

    // Set AtlasLabel
    param = confSet->Get("AtlasLabel");
    if( param != "" ){
        fFitter->fAtlasLabel = RemoveQuotes(param);
    }

    // Set CleanTables
    param = confSet->Get("CleanTables");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fCleanTables = true;
        else if (param == "FALSE") fFitter->fCleanTables = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'CleanTables' option but you did not provide valid setting. Using default (false)");
            fFitter->fCleanTables = false;
        }
    }

    // Set SystCategoryTables
    param = confSet->Get("SystCategoryTables");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fSystCategoryTables = true;
        else if (param == "FALSE") fFitter->fSystCategoryTables = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'SystCategoryTables' option but you did not provide valid setting. Using default (false)");
            fFitter->fSystCategoryTables = false;
        }
    }

    // Set Suffix
    param = confSet->Get("Suffix");
    if( param != "" ){
        fFitter->fSuffix = RemoveQuotes(param);
    }

    // Set SaveSuffix
    param = confSet->Get("SaveSuffix");
    if( param != "" ){
        fFitter->fSaveSuffix = RemoveQuotes(param);
    }

    // Set HideNP
    param = confSet->Get("HideNP");
    if( param != "" ){
        fFitter->fVarNameHide = Vectorize(param,',');
    }

    // Set RegionGroups
    param = confSet->Get("RegionGroups");
    if( param != "" ) {
        std::vector<std::string> groups = Vectorize(param,',');
        for(const std::string& igroup : groups) fFitter->fRegionGroups.push_back(igroup);
    }

    // Set KeepPrefitBlindedBins
    param = confSet->Get("KeepPrefitBlindedBins");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fKeepPrefitBlindedBins = true;
        else if (param == "FALSE") fFitter->fKeepPrefitBlindedBins = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'KeepPrefitBlindedBins' option but you did not provide valid setting. Using default (false)");
            fFitter->fKeepPrefitBlindedBins = false;
        }
    }

    // Set CustomAsimov
    param = confSet->Get("CustomAsimov");
    if( param != "" ){
        fFitter->fCustomAsimov = RemoveQuotes(param);
    }

    // Set GetChi2
    param = confSet->Get("GetChi2");
    if( param != "" ){ // can be TRUE, SYST+STAT, STAT-ONLY... (if it contains STAT and no SYST => stat-only, ptherwise stat+syst)
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fGetChi2 = 2;
        }
        else if( param.find("SYST")!=std::string::npos ){
            fFitter->fGetChi2 = 2;
        }
        else if( param.find("STAT")!=std::string::npos ){
            fFitter->fGetChi2 = 1;
        } else {
            WriteErrorStatus("ConfigReader::ReadJobOptions", "You specified 'GetChi2' option but you did not provide valid option. Check this!");
            return 1;
        }
    }

    // Set DoTables
    param = confSet->Get("DoTables");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoTables = true;
        else if( param == "FALSE" ) fFitter->fDoTables = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'DoTables' option but you did not provide valid setting. Using default (false)");
            fFitter->fDoTables = false;
        }
    }

    // Set CustomFunctions
    param = confSet->Get("CustomFunctions");
    if( param != "" ) {
        fFitter->fCustomFunctions = Vectorize(param,',');
    }

    // Set Bootstrap
    param = confSet->Get("Bootstrap");
    if( param != "" ){
        fFitter->fBootstrap = RemoveQuotes(param);
    }

    // Set Bootstrap
    param = confSet->Get("BootstrapSyst");
    if( param != "" ){
        fFitter->fBootstrapSyst = RemoveQuotes(param);
    }

    // Set DecorrSuff
    param = confSet->Get("DecorrSuff");
    if( param != ""){
        fFitter->fDecorrSuff = RemoveQuotes(param);
    }

    // Set DecorrSysts
    param = confSet->Get("DecorrSysts");
    if( param != ""){
        fFitter->fDecorrSysts = Vectorize(param,',');
    }

    // Set SmoothMorphingTemplates
    param = confSet->Get("SmoothMorphingTemplates");
    if ( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE"){
            fFitter->fSmoothMorphingTemplates = "TRUE";
        }
        else if (param == "FALSE"){
            fFitter->fSmoothMorphingTemplates = "";
        } else {
            fFitter->fSmoothMorphingTemplates = param;
        }
    }

    // Set POIPrecision
    param = confSet->Get("POIPrecision");
    if( param != ""){
        fFitter->fPOIPrecision = stoi(param);
        if (fFitter->fPOIPrecision < 1 || fFitter->fPOIPrecision > 5){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "Parameter POIPrecision has value smaller than 1 or larger than 5. Using default (2).");
            fFitter->fPOIPrecision = 2;
        }
    }

    // Set UseGammasForCorr
    param = confSet->Get("UseGammasForCorr");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE"){
            fFitter->fuseGammasForCorr = true;
        } else if (param == "FALSE") {
            fFitter->fuseGammasForCorr = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified UseGammasForCorr option but did not provide valid parameter. Using default (false)");
            fFitter->fuseGammasForCorr = false;
        }
    }

    // Set PropagateSystsForMorphing
    param = confSet->Get("PropagateSystsForMorphing");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE"){
            fFitter->fPropagateSystsForMorphing = true;
        } else if (param == "FALSE") {
            fFitter->fPropagateSystsForMorphing = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified PropagateSystsForMorphing option but did not provide valid parameter. Using default (false)");
            fFitter->fPropagateSystsForMorphing = false;
        }
    }

    // Set UseATLASRounding
    param = confSet->Get("UseATLASRounding");
    if ( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fUseATLASRoundingTxt = true;
            fFitter->fUseATLASRoundingTex = true;
        } else if (param == "FALSE") {
            fFitter->fUseATLASRoundingTxt = false;
            fFitter->fUseATLASRoundingTex = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified UseATLASRounding option but did not provide valid parameter. Using default (false)");
            fFitter->fUseATLASRoundingTxt = false;
            fFitter->fUseATLASRoundingTex = false;
        }
    }
    param = confSet->Get("UseATLASRoundingTxt");
    if ( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fUseATLASRoundingTxt = true;
        } else if (param == "FALSE") {
            fFitter->fUseATLASRoundingTxt = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified UseATLASRoundingTxt option but did not provide valid parameter. Using default (false)");
            fFitter->fUseATLASRoundingTxt = false;
        }
    }
    param = confSet->Get("UseATLASRoundingTex");
    if ( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fUseATLASRoundingTex = true;
        } else if (param == "FALSE") {
            fFitter->fUseATLASRoundingTex = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified UseATLASRoundingTex option but did not provide valid parameter. Using default (false)");
            fFitter->fUseATLASRoundingTex = false;
        }
    }

    // Set RankingPOIName
    param = confSet->Get("RankingPOIName");
    if( param != ""){
        fFitter->fRankingPOIName = RemoveQuotes(param);
    }

    // Set PruningType
    param = confSet->Get("PruningType");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "SEPARATESAMPLE") {
            fFitter->fPruningType = TRExFit::SEPARATESAMPLE;
        } else if (param == "BACKGROUNDREFERENCE") {
            fFitter->fPruningType = TRExFit::BACKGROUNDREFERENCE;
        } else if (param == "COMBINEDREFERENCE") {
            fFitter->fPruningType = TRExFit::COMBINEDREFERENCE;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified PruningType option but did not provide valid parameter. Using default (SEPARATESAMPLE)");
            fFitter->fPruningType = TRExFit::SEPARATESAMPLE;
        }
    }

    param = confSet->Get("DoSystNormalizationPlots");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fDoSystNormalizationPlots = true;
        } else if (param == "FALSE") {
            fFitter->fDoSystNormalizationPlots = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified DoSystNormalizationPlots option but did not provide valid parameter. Using default (TRUE)");
            fFitter->fDoSystNormalizationPlots = true;
        }
    }

    // Set PrePostFitCanvasSize
    param = confSet->Get("PrePostFitCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param, ',');
        if (tmp.size() != 2){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified PrePostFitCanvasSize option but did not provide 2 parameters. Ignoring.");
        }
        const int& x = std::stoi(tmp.at(0));
        const int& y = std::stoi(tmp.at(1));

        if (x <= 100 || y <= 100){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified PrePostFitCanvasSize option but at least one parameter is <= 100. Ignoring.");
        }
        fFitter->fPrePostFitCanvasSize.emplace_back(x);
        fFitter->fPrePostFitCanvasSize.emplace_back(y);
    }

    // Set SummaryCanvasSize
    param = confSet->Get("SummaryCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param, ',');
        if (tmp.size() != 2){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified SummaryCanvasSize option but did not provide 2 parameters. Ignoring.");
        }
        const int& x = std::stoi(tmp.at(0));
        const int& y = std::stoi(tmp.at(1));

        if (x <= 100 || y <= 100){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified SummaryCanvasSize option but at least one parameter is <= 100. Ignoring.");
        }
        fFitter->fSummaryCanvasSize.emplace_back(x);
        fFitter->fSummaryCanvasSize.emplace_back(y);
    }

    // Set MergeCanvasSize
    param = confSet->Get("MergeCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param, ',');
        if (tmp.size() != 2){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified MergeCanvasSize option but did not provide 2 parameters. Ignoring.");
        }
        const int& x = std::stoi(tmp.at(0));
        const int& y = std::stoi(tmp.at(1));

        if (x <= 100 || y <= 100){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified MergeCanvasSize option but at least one parameter is <= 100. Ignoring.");
        }
        fFitter->fMergeCanvasSize.emplace_back(x);
        fFitter->fMergeCanvasSize.emplace_back(y);
    }

    // Set PieChartCanvasSize
    param = confSet->Get("PieChartCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param, ',');
        if (tmp.size() != 2){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified PieChartCanvasSize option but did not provide 2 parameters. Ignoring.");
        }
        const int& x = std::stoi(tmp.at(0));
        const int& y = std::stoi(tmp.at(1));

        if (x <= 100 || y <= 100){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified PieChartCanvasSize option but at least one parameter is <= 100. Ignoring.");
        }
        fFitter->fPieChartCanvasSize.emplace_back(x);
        fFitter->fPieChartCanvasSize.emplace_back(y);
    }

    // Set NPRankingCanvasSize
    param = confSet->Get("NPRankingCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param, ',');
        if (tmp.size() != 2){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NPRankingCanvasSize option but did not provide 2 parameters. Ignoring.");
        }
        const int& x = std::stoi(tmp.at(0));
        const int& y = std::stoi(tmp.at(1));

        if (x <= 100 || y <= 100){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NPRankingCanvasSize option but at least one parameter is <= 100. Ignoring.");
        }
        fFitter->fNPRankingCanvasSize.emplace_back(x);
        fFitter->fNPRankingCanvasSize.emplace_back(y);
    }

    // Set ReorderNPs
    param = confSet->Get("ReorderNPs");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fReorderNPs = true;
        else if( param == "FALSE" ) fFitter->fReorderNPs = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'ReorderNPs' option but you did not provide valid setting. Using default (false)");
            fFitter->fReorderNPs = false;
        }
    }

    // success
    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetJobPlot(ConfigSet *confSet){

    // Plot option
    std::string param = confSet->Get("PlotOptions");
    std::vector<std::string> vec;
    if( param != ""){
        vec = Vectorize(RemoveQuotes(param),',');
        if( std::find(vec.begin(), vec.end(), "YIELDS") !=vec.end() )  TRExFitter::SHOWYIELDS     = true;
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TRExFitter::SHOWSTACKSIG   = false;
        if( std::find(vec.begin(), vec.end(), "SKIPSIG")!=vec.end() )  TRExFitter::ADDSTACKSIG    = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TRExFitter::SHOWNORMSIG    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TRExFitter::SHOWOVERLAYSIG = true;
        if( std::find(vec.begin(), vec.end(), "LEFT")   !=vec.end() )  TRExFitter::LEGENDLEFT     = true; // forces leg entry to be left-aligned when adding yields to legend
        if( std::find(vec.begin(), vec.end(), "RIGHT")  !=vec.end() )  TRExFitter::LEGENDRIGHT    = true; // forces leg entry to be right-aligned even when no yields in legend
        if( std::find(vec.begin(), vec.end(), "CHI2")   !=vec.end() )  TRExFitter::SHOWCHI2       = true;
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TRExFitter::PREFITONPOSTFIT= true;
        if( std::find(vec.begin(), vec.end(), "POISSONIZE")        !=vec.end() )  TRExFitter::POISSONIZE     = true;
        if( std::find(vec.begin(), vec.end(), "NOXERR") !=vec.end() )  TRExFitter::REMOVEXERRORS  = true;
        if( std::find(vec.begin(), vec.end(), "OPRATIO") !=vec.end() ) TRExFitter::OPRATIO        = true;
        if( std::find(vec.begin(), vec.end(), "NORATIO") !=vec.end() ) TRExFitter::NORATIO        = true;
    }

    // Set PlotOptionsSummary
    param = confSet->Get("PlotOptionsSummary");
    if( param != ""){
        vec = Vectorize(RemoveQuotes(param),',');
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TRExFitter::SHOWSTACKSIG_SUMMARY   = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TRExFitter::SHOWNORMSIG_SUMMARY    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TRExFitter::SHOWOVERLAYSIG_SUMMARY = true;
    }
    else{
        WriteDebugStatus("ConfigReader::SetJobPlot", "PlotOptionsSummary not specified setting Summary values to 'PlotOptions'");
        TRExFitter::SHOWSTACKSIG_SUMMARY   = TRExFitter::SHOWSTACKSIG    ;
        TRExFitter::SHOWNORMSIG_SUMMARY    = TRExFitter::SHOWNORMSIG     ;
        TRExFitter::SHOWOVERLAYSIG_SUMMARY = TRExFitter::SHOWOVERLAYSIG  ;
    }

    // Set SystControlPlots
    param = confSet->Get("SystControlPlots");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TRExFitter::SYSTCONTROLPLOTS = true;
        } else if (param == "FALSE"){
            TRExFitter::SYSTCONTROLPLOTS = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'SystControlPlots' option but you did not provide valid setting. Using default (TRUE)");
            TRExFitter::SYSTCONTROLPLOTS = true;
        }
    }

    // Set SystDataPlots
    param = confSet->Get("SystDataPlots");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TRExFitter::SYSTDATAPLOT = true;
            fFitter->fSystDataPlot_upFrame=true;
        } else if( param == "FALSE" ){
            TRExFitter::SYSTDATAPLOT = false;
            fFitter->fSystDataPlot_upFrame= false;
        } else if( param == "FILLUPFRAME" ){
            TRExFitter::SYSTDATAPLOT = true;
            fFitter->fSystDataPlot_upFrame=true;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'SystDataPlots' option but the value is not 'TRUE' nor 'FILLUPFRAME'. Setting SystDataPlot and SystDataPlot_upFrame to false");
            TRExFitter::SYSTDATAPLOT = false;
            fFitter->fSystDataPlot_upFrame=false;
        }
    }

    // Set SystErrorBars
    param = confSet->Get("SystErrorBars");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE" ){
            TRExFitter::SYSTERRORBARS = true;
        } else if (param == "FALSE"){
            TRExFitter::SYSTERRORBARS = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'SystErrorBars' option but you did not provide valid setting. Using default (TRUE)");
            TRExFitter::SYSTERRORBARS = true;
        }
    }

    // Set GuessMCStatEmptyBins
    param = confSet->Get("GuessMCStatEmptyBins");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TRExFitter::GUESSMCSTATERROR = true;
        } else if (param == "FALSE") {
            TRExFitter::GUESSMCSTATERROR = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'GuessMCStatEmptyBins' option but you did not provide valid setting. Using default (TRUE)");
            TRExFitter::GUESSMCSTATERROR = true;
        }
    }

    // Set CorrectNormForNegativeIntegral
    param = confSet->Get("CorrectNormForNegativeIntegral");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL = true;
        } else if (param == "FALSE") {
            TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'CorrectNormForNegativeIntegral' option but you did not provide valid setting. Using default (FALSE)");
            TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL = false;
        }
    }

    // Set SuppressNegativeBinWarnings
    param = confSet->Get("SuppressNegativeBinWarnings");
     if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if ( param == "TRUE" ){
            fFitter->fSuppressNegativeBinWarnings = true;
        } else if (param == "FALSE"){
            fFitter->fSuppressNegativeBinWarnings = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified SuppressNegativeBinWarnings option but did not provide valid parameter. Using default (false)");
            fFitter->fSuppressNegativeBinWarnings = false;
        }
    }

    // Set SignalRegionsPlot
    param = confSet->Get("SignalRegionsPlot");
    if(param != ""){
        fFitter->fRegionsToPlot = Vectorize(param,',');
    }

    // Set SummaryPlotRegions
    param = confSet->Get("SummaryPlotRegions");
    if(param != ""){
        fFitter->fSummaryPlotRegions = Vectorize(param,',');
    }

    // Set SummaryPlotLabels
    param = confSet->Get("SummaryPlotLabels");
    if(param != ""){
        fFitter->fSummaryPlotLabels = Vectorize(param,',');
    }

    // Set SummaryPlotValidationRegions
    param = confSet->Get("SummaryPlotValidationRegions");
    if(param != ""){
        fFitter->fSummaryPlotValidationRegions = Vectorize(param,',');
    }

    // Set SummaryPlotValidationLabels
    param = confSet->Get("SummaryPlotValidationLabels");
    if(param != ""){
        fFitter->fSummaryPlotValidationLabels = Vectorize(param,',');
    }

    // Set SummaryPlotYmin
    param = confSet->Get("SummaryPlotYmin");
    if(param != "") fFitter->fYmin = atof(param.c_str());

    // Set SummaryPlotYmax
    param = confSet->Get("SummaryPlotYmax");
    if(param != "") fFitter->fYmax = atof(param.c_str());

    // Set SummaryLogY
    param = confSet->Get("SummaryLogY");
    if(param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fSummaryLogY = true;
        } else if (param == "FALSE") {
            fFitter->fSummaryLogY = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified SummaryLogY option but did not provide valid parameter. Using default (TRUE)");
        }
    }

    // Set RatioYmin
    param = confSet->Get("RatioYmin");
    if(param != "") {
        fFitter->fRatioYmin = atof(param.c_str());
        fFitter->fRatioYminPostFit = fFitter->fRatioYmin;
    }

    // Set RatioYmax
    param = confSet->Get("RatioYmax");
    if(param != ""){
        fFitter->fRatioYmax = atof(param.c_str());
        fFitter->fRatioYmaxPostFit = fFitter->fRatioYmax;
    }

    // Set RatioYminPostFit
    param = confSet->Get("RatioYminPostFit");
    if(param != "") fFitter->fRatioYminPostFit = atof(param.c_str());

    // Set RatioYmaxPostFit
    param = confSet->Get("RatioYmaxPostFit");
    if(param != "") fFitter->fRatioYmaxPostFit = atof(param.c_str());

    // Set RatioTitle
    param = confSet->Get("RatioYtitle");
    if(param != "") fFitter->fRatioYtitle = RemoveQuotes(param.c_str());

    // Set RatioTitle
    param = confSet->Get("RatioType");
    if(param != "") fFitter->fRatioType = RemoveQuotes(param.c_str());

    // Set Label and Legend position
    param = confSet->Get("LabelX");
    if(param != "") fFitter->fLabelX = atof(param.c_str());
    param = confSet->Get("LabelY");
    if(param != "") fFitter->fLabelY = atof(param.c_str());
    param = confSet->Get("LegendX1");
    if(param != "") fFitter->fLegendX1 = atof(param.c_str());
    param = confSet->Get("LegendX2");
    if(param != "") fFitter->fLegendX2 = atof(param.c_str());
    param = confSet->Get("LegendY");
    if(param != "") fFitter->fLegendY = atof(param.c_str());

    // Set Label and Legend position for Summary
    param = confSet->Get("LabelXSummary");
    if(param != "") fFitter->fLabelXSummary = atof(param.c_str());
    param = confSet->Get("LabelYSummary");
    if(param != "") fFitter->fLabelYSummary = atof(param.c_str());
    param = confSet->Get("LegendX1Summary");
    if(param != "") fFitter->fLegendX1Summary = atof(param.c_str());
    param = confSet->Get("LegendX2Summary");
    if(param != "") fFitter->fLegendX2Summary = atof(param.c_str());
    param = confSet->Get("LegendYSummary");
    if(param != "") fFitter->fLegendYSummary = atof(param.c_str());

    // Set Label and Legend position for Merge
    param = confSet->Get("LabelXMerge");
    if(param != "") fFitter->fLabelXMerge = atof(param.c_str());
    param = confSet->Get("LabelYMerge");
    if(param != "") fFitter->fLabelYMerge = atof(param.c_str());
    param = confSet->Get("LegendX1Merge");
    if(param != "") fFitter->fLegendX1Merge = atof(param.c_str());
    param = confSet->Get("LegendX2Merge");
    if(param != "") fFitter->fLegendX2Merge = atof(param.c_str());
    param = confSet->Get("LegendYMerge");
    if(param != "") fFitter->fLegendYMerge = atof(param.c_str());

    // Set LegendNColumns
    param = confSet->Get("LegendNColumns");
    if(param != "") fFitter->fLegendNColumns = atoi(param.c_str());
    param = confSet->Get("LegendNColumnsMerge");
    if(param != "") fFitter->fLegendNColumnsMerge = atoi(param.c_str());
    param = confSet->Get("LegendNColumnsSummary");
    if(param != "") fFitter->fLegendNColumnsSummary = atoi(param.c_str());

    // Set DoSummaryPlot
    param = confSet->Get("DoSummaryPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoSummaryPlot = true;
        else if( param == "FALSE" ) fFitter->fDoSummaryPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoSummaryPlot' option but did not provide valid parameter. Using default (false)");
            fFitter->fDoSummaryPlot = false;
        }
    }

    // Set DoMergedPlot
    param = confSet->Get("DoMergedPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoMergedPlot = true;
        else if( param == "FALSE" ) fFitter->fDoMergedPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoMergedPlot' option but did not provide valid parameter. Using default (false)");
            fFitter->fDoMergedPlot = false;
        }
    }

    // Set DoSignalRegionsPlot
    param = confSet->Get("DoSignalRegionsPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoSignalRegionsPlot = true;
        else if( param == "FALSE" ) fFitter->fDoSignalRegionsPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoSignalRegionsPlot' option but did not provide valid parameter. Using default (false)");
            fFitter->fDoSignalRegionsPlot = false;
        }
    }
    param = confSet->Get("DoPieChartPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoPieChartPlot = true;
        else if( param == "FALSE" ) fFitter->fDoPieChartPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoPieChartPlot' option but did not provide valid parameter. Using default (false)");
            fFitter->fDoPieChartPlot = false;
        }
    }

    // Set RankingPlot
    param = confSet->Get("RankingPlot");
    if( param != ""){
        fFitter->fRankingPlot = RemoveQuotes(param);
    }

    // exclude a template from morphing
    param = confSet->Get("ExcludeFromMorphing");
    if( param != ""){
        fFitter->fExcludeFromMorphing = CheckName(param);
    }

    param = confSet->Get("ScaleSamplesToData");
    if( param != ""){
        fFitter->fScaleSamplesToData = Vectorize(param,',');
    }

    param = confSet->Get("MaxNtupleEvents");
    if(param != "") fFitter->fDebugNev = atoi(param.c_str());

    param = confSet->Get("PruningShapeOption");
    if(param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "MAXBIN") {
            fFitter->fPruningShapeOption = PruningUtil::SHAPEOPTION::MAXBIN;
        } else if (param == "KSTEST") {
            fFitter->fPruningShapeOption = PruningUtil::SHAPEOPTION::KSTEST;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'PruningShapeOption' option but did not provide valid parameter. Using default (MAXBIN)");
            fFitter->fPruningShapeOption = PruningUtil::SHAPEOPTION::MAXBIN;
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadGeneralOptions(){
    ConfigSet* confSet = fParser->GetConfigSet("Options");
    if (confSet != nullptr){
        for(int i=0; i < confSet->GetN(); i++){
            if(confSet->GetConfigValue(i) != ""){
                TRExFitter::OPTION[confSet->GetConfigName(i)] = atof(confSet->GetConfigValue(i).c_str());
            }
        }
    } else {
        WriteDebugStatus("ConfigReader::ReadGeneralOptions", "You do not have 'Options' option in the config. It is ok, we just want to let you know.");
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadFitOptions(){
    std::string param = "";

    ConfigSet *confSet = fParser->GetConfigSet("Fit");
    if (confSet == nullptr){
        WriteInfoStatus("ConfigReader::ReadFitOptions", "You do not have Fit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set FitType
    param = confSet->Get("FitType");
    if( param != "" && fFitter->fFitType == TRExFit::UNDEFINED ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "SPLUSB" ){
            fFitter->SetFitType(TRExFit::SPLUSB);
        }
        else if( param == "BONLY" ){
            fFitter->SetFitType(TRExFit::BONLY);
        }
        else if( param == "UNFOLDING" ){
            fFitter->SetFitType(TRExFit::UNFOLDING);
        }
        else{
            WriteErrorStatus("ConfigReader::ReadFitOptions", "Unknown FitType argument : " + confSet->Get("FitType"));
            return 1;
        }
    }
    else if( fFitter->fFitType == TRExFit::UNDEFINED ){
        WriteWarningStatus("ConfigReader::ReadFitOptions","Setting default fit Type SPLUSB");
        fFitter->SetFitType(TRExFit::SPLUSB);
    }

    // Set FitRegion
    param = confSet->Get("FitRegion");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if( param != "" ){
        if( param == "CRONLY" ){
            fFitter->SetFitRegion(TRExFit::CRONLY);
        }
        else if( param == "CRSR" ){
            fFitter->SetFitRegion(TRExFit::CRSR);
        }
        else{
            fFitter->SetFitRegion(TRExFit::USERSPECIFIC);
            fFitter->fFitRegionsToFit = Vectorize(param,',');
            if(fFitter->fFitRegionsToFit.size()==0){
                WriteErrorStatus("ConfigReader::ReadFitOptions", "Unknown FitRegion argument : " + confSet->Get("FitRegion"));
                return 1;
            }
        }
    }

    // Set FitBlind
    param = confSet->Get("FitBlind");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fFitIsBlind = true;
        } else if ( param == "FALSE" ){
            fFitter->fFitIsBlind = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'FitBlind' option but did not provide valid parameter. Using default (false)");
            fFitter->fFitIsBlind = false;
        }
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
         fFitter->fFitPOIAsimov = atof(param.c_str());
    }

    // Set NPValues
    param = confSet->Get("NPValues");
    if( param != "" ){
        std::vector < std::string > temp_vec = Vectorize(param,',',false);
        for(std::string iNP : temp_vec){
            std::vector < std::string > np_value = Vectorize(iNP,':');
            if(np_value.size()==2){
                fFitter->fFitNPValues.insert( std::pair < std::string, double >( np_value[0], atof(np_value[1].c_str()) ) );
            } else {
                WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'NPValues' option but did not provide 2 parameters for each NP which is expected. Ignoring");
            }
        }
    }

    // Set FixNPs
    param = confSet->Get("FixNPs");
    if( param != "" ){
        std::vector < std::string > temp_fixedNPs = Vectorize(param,',',false);
        for(std::string iNP : temp_fixedNPs){
            std::vector < std::string > fixed_nps = Vectorize(iNP,':');
            if(fixed_nps.size()==2){
                fFitter->fFitFixedNPs.insert( std::pair < std::string, double >( fixed_nps[0], atof(fixed_nps[1].c_str()) ) );
            } else {
                WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'FixNPs' option but did not provide 2 parameters for each NP which is expected. Ignoring");
            }
        }
    }

    // Set doLHscan
    param = confSet->Get("doLHscan");
    if( param != "" ){
        if (fOnlyLHscan==""){
            fFitter->fVarNameLH = Vectorize(param,',');
        } else {
            fFitter->fVarNameLH.emplace_back(fOnlyLHscan);
        }
    }

    // Set do2DLHscan
    param = confSet->Get("do2DLHscan");
    if( param != "" ){
        const std::vector<std::string> tmp = Vectorize(param,':');
        for (const auto& ivec : tmp){
            const std::vector<std::string> v = Vectorize(ivec,',');
            if (v.size() == 2){
                fFitter->fVarName2DLH.emplace_back(v);
            } else {
                WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'do2DLHscan' option but did not provide correct input. Ignoring");
            }
        }
    }

    // Set LHscanMin
    param = confSet->Get("LHscanMin");
    if ( param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'LHscanMin' option but did not set doLHscan. Ignoring");
        } else {
            fFitter->fLHscanMin = std::stof(param);
        }
    }

    // Set LHscanMax
    param = confSet->Get("LHscanMax");
    if ( param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'LHscanMax' option but did not set doLHscan. Ignoring");
        } else {
            fFitter->fLHscanMax = std::stof(param);
        }
    }

    // Set LHscanSteps
    param = confSet->Get("LHscanSteps");
    if ( param != "" ) {
        fFitter->fLHscanSteps = std::stoi(param);
        if(fFitter->fLHscanSteps < 3 || fFitter->fLHscanSteps > 100){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "LHscanSteps is smaller than 3 or larger than 100, setting to default (30)");
            fFitter->fLHscanSteps = 30;
        }
    }
    // Set LHscanMin for second variable
    param = confSet->Get("LHscanMinY");
    if ( param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'LHscanMinY' option but did not set doLHscan. Ignoring");
        } else {
            fFitter->fLHscanMinY = std::stof(param);
        }
    }

    // Set LHscanMax for second variable
    param = confSet->Get("LHscanMaxY");
    if ( param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'LHscanMaxY' option but did not set doLHscan. Ignoring");
        } else {
            fFitter->fLHscanMaxY = std::stof(param);
        }
    }

    // Set LHscanSteps for second variable
    param = confSet->Get("LHscanStepsY");
    if ( param != "" ) {
        if (fFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'LHscanStepsY' option but did not set doLHscan. Ignoring");
        } else {
            fFitter->fLHscanStepsY = std::stoi(param);
            if(fFitter->fLHscanStepsY < 3 || fFitter->fLHscanStepsY > 100){
                WriteWarningStatus("ConfigReader::ReadFitOptions", "LHscanSteps is smaller than 3 or larger than 100, setting to default (30)");
                fFitter->fLHscanStepsY = fFitter->fLHscanSteps;
            }
        }
    }
    else {
        fFitter->fLHscanStepsY = fFitter->fLHscanSteps;
    }

    // Set Paral2D
    param = confSet->Get("Parallel2Dscan");
    if ( param != "" ) {
        if( param == "TRUE" ){
            fFitter->fParal2D = true;
        } else if (param == "FALSE") {
            fFitter->fParal2D = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'Parallel2Dscan' option but did not provide valid parameter. Using default (false)");
            fFitter->fParal2D = false;
        }
    }

    // Set Paral2Dstep
    param = confSet->Get("Parallel2DscanStep");
    if ( param != "" ) {
        fFitter->fParal2Dstep = std::atoi( param.c_str());
        if (fFitter->fParal2Dstep < 1 || fFitter->fParal2Dstep>=fFitter->fLHscanSteps ){
            WriteErrorStatus("ConfigReader::ReadFitOptions", "You specified a step for 2D LHscan outside the allowed range.");
            return 1;
        }
    }

    // Set UseMinos
    param = confSet->Get("UseMinos");
    if( param != "" ){
        fFitter->fVarNameMinos = Vectorize(param,',');
    }

    // Set SetRandomInitialNPval
    param = confSet->Get("SetRandomInitialNPval");
    if( param != ""){
        fFitter->fUseRnd = true;
        fFitter->fRndRange = std::atof(param.c_str());
    }

    // Set SetRandomInitialNPvalSeed
    param = confSet->Get("SetRandomInitialNPvalSeed");
    if( param != ""){
        fFitter->fRndSeed = std::atol(param.c_str());
    }

    // Set NumCPU
    param = confSet->Get("NumCPU");
    if( param != "" ){
        TRExFitter::NCPU = std::atoi( param.c_str());
    }

    // Set StatOnlyFit
    param = confSet->Get("StatOnlyFit");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fStatOnlyFit = true;
        } else if (param == "FALSE") {
            fFitter->fStatOnlyFit = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'StatOnlyFit' option but did not provide valid parameter. Using default (false)");
            fFitter->fStatOnlyFit = false;
        }
    }

    // Set GetGoodnessOfFit
    param = confSet->Get("GetGoodnessOfFit");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fGetGoodnessOfFit = true;
            fFitter->fSaturatedModel = true;
        } else if (param == "FALSE"){
            fFitter->fGetGoodnessOfFit = false;
            fFitter->fSaturatedModel = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'GetGoodnessOfFit' option but did not provide valid parameter. Using default (false)");
            fFitter->fGetGoodnessOfFit = false;
            fFitter->fSaturatedModel = false;
        }
    }

    // Set DoNonProfileFit
    param = confSet->Get("DoNonProfileFit");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fDoNonProfileFit = true;
        } else if (param == "FALSE"){
            fFitter->fDoNonProfileFit = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'DoNonProfileFit' option but did not provide valid parameter. Using default (false)");
            fFitter->fDoNonProfileFit = false;
        }
    }

    // Set NonProfileFitSystThreshold
    param = confSet->Get("NonProfileFitSystThreshold");
    if( param != "" ){
        fFitter->fNonProfileFitSystThreshold = std::atof(param.c_str());
    }

    // Set FitToys
    param = confSet->Get("FitToys");
    if( param != "" ){
        fFitter->fFitToys = std::atoi( param.c_str());
    }

    // Set ToysHistoMin
    param = confSet->Get("ToysHistoMin");
    if( param != "" ){
        fFitter->fToysHistoMin = std::stof( param.c_str());
    }

    // Set ToysHistoMax
    param = confSet->Get("ToysHistoMax");
    if( param != "" ){
        fFitter->fToysHistoMax = std::stof( param.c_str());

        if (fFitter->fToysHistoMin > fFitter->fToysHistoMax){
            WriteErrorStatus("ConfigReader::ReadFitOptions", "Minimum for toys is larger than maximum for toys");
            return 1;
        }
    }

    // Set ToysHistoNbins
    param = confSet->Get("ToysHistoNbins");
    if( param != "" ){
        fFitter->fToysHistoNbins = std::atoi( param.c_str());
        if (fFitter->fToysHistoNbins < 2){
            WriteErrorStatus("ConfigReader::ReadFitOptions", "Number of bins for toys is < 2");
            return 1;
        }
    }

    // Set ToysPseudodataNP
    param = confSet->Get("ToysPseudodataNP");
    if( param != "" ){
        fFitter->fToysPseudodataNP = param.c_str();
    }

    // Set ToysPseudodataNPShift
    param = confSet->Get("ToysPseudodataNPShift");
    if( param != "" ){
        fFitter->fToysPseudodataNPShift = std::atoi( param.c_str() );
        if (fFitter->fToysPseudodataNPShift < - 3. || fFitter->fToysPseudodataNPShift > 3.) {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "ToysPseudodataNPShift smaller than -3 or larger than +3. This is probably wrong... Setting to default (+1).");
        }
    }

    // Set TemplateInterpolationOption
    param = confSet->Get("TemplateInterpolationOption");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "LINEAR"){
            fFitter->fTemplateInterpolationOption = TRExFit::LINEAR;
        } else if (param == "SMOOTHLINEAR"){
            fFitter->fTemplateInterpolationOption = TRExFit::SMOOTHLINEAR;
        } else if (param == "SQUAREROOT"){
            fFitter->fTemplateInterpolationOption = TRExFit::SQUAREROOT;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'TemplateInterpolationOption' option but did not provide valid parameter. Using default (LINEAR)");
            fFitter->fTemplateInterpolationOption = TRExFit::LINEAR;
        }
    }

    // Set BlindedParameters
    param = confSet->Get("BlindedParameters");
    if( param != "" ){
        std::vector<std::string> tmp = Vectorize( param.c_str(),',');
        for (auto itmp : tmp){
            itmp = RemoveQuotes(itmp);
        }
        fFitter->fBlindedParameters = tmp;
        if (fFitter->fBlindedParameters.size() ==  0){
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'BlindedParameters' option but you did not provide valid parameters. Ignoring");
            fFitter->fBlindedParameters.clear();
        }
    }

    // Saturated Model fit for goodness of fit
    param = confSet->Get("SaturatedModel");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fSaturatedModel = true;
        } else if (param == "FALSE"){
            fFitter->fSaturatedModel = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'SaturatedModel' option but did not provide valid parameter. Using default (false)");
            fFitter->fSaturatedModel = false;
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadLimitOptions(){
    std::string param = "";

    ConfigSet* confSet = fParser->GetConfigSet("Limit");
    if (confSet == nullptr){
        WriteDebugStatus("ConfigReader::ReadLimitOptions", "You do not have Limit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set LimitType
    param = confSet->Get("LimitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "ASYMPTOTIC" ){
            fFitter->SetLimitType(TRExFit::ASYMPTOTIC);
        }
        else if( param == "TOYS" ){
            fFitter->SetLimitType(TRExFit::TOYS);
        }
        else{
            WriteErrorStatus("ConfigReader::ReadLimitOptions", "Unknown LimitType argument : " + confSet->Get("LimitType"));
            return 1;
        }
    }

    // Set LimitBlind
    param = confSet->Get("LimitBlind");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fLimitIsBlind = true;
        } else if ( param == "FALSE" ){
            fFitter->fLimitIsBlind = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadLimitOptions", "You specified 'LimitBlind' option but did not provide valid parameter. Using default (false)");
            fFitter->fLimitIsBlind = false;
        }
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
        fFitter->fLimitPOIAsimov = atof(param.c_str());
    }

    // Set SignalInjection
    param = confSet->Get("SignalInjection");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fSignalInjection = true;
        } else if ( param == "FALSE" ){
            fFitter->fSignalInjection = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadLimitOptions", "You specified 'SignalInjection' option but did not provide valid parameter. Using default (false)");
            fFitter->fSignalInjection = false;
        }
    }

    // Set SignalInjectionValue
    param = confSet->Get("SignalInjectionValue");
    if( param != "" ){
        fFitter->fSignalInjectionValue = std::stof(param);
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fFitter->fLimitParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fFitter->fLimitParamValue = std::stof(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fFitter->fLimitOutputPrefixName = param;
    }

    param = confSet->Get("ConfidenceLevel");
    if( param != "" ){
        double conf = std::stod(param);
        if (conf <= 0 || conf >= 1){
            WriteWarningStatus("ConfigReader::ReadLimitOptions", "Confidence level is <= 0 or >=1. Setting to default 0.95");
            conf = 0.95;
        }
        fFitter->fLimitsConfidence = conf;
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadSignificanceOptions(){
    std::string param = "";

    ConfigSet* confSet = fParser->GetConfigSet("Significance");
    if (confSet == nullptr){
        WriteDebugStatus("ConfigReader::ReadSignificanceOptions", "You do not have Significance option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set LimitBlind
    param = confSet->Get("SignificanceBlind");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fSignificanceIsBlind = true;
        } else if ( param == "FALSE" ){
            fFitter->fSignificanceIsBlind = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadSignificanceOptions", "You specified 'SignificanceBlind' option but did not provide valid parameter. Using default (false)");
            fFitter->fSignificanceIsBlind = false;
        }
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
        fFitter->fSignificancePOIAsimov = atof(param.c_str());
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fFitter->fSignificanceParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fFitter->fSignificanceParamValue = std::stof(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fFitter->fSignificanceOutputPrefixName = param;
    }


    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadRegionOptions(const std::string& opt){

    fAvailableRegions = GetAvailableRegions();

    // Check of the regions from commands like exist
    if (fOnlyRegions.size() > 0){
        if (!CheckPresence(fOnlyRegions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "You set regions that do not exist in your command line options");
            } else {
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You set regions that do not exist in your command line options");
                return 1;
            }
        }
    }

    fHasAtLeastOneValidRegion = false;

    int nReg = 0;
    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("Region",nReg);
        if (confSet == nullptr) break;

        nReg++;
        if(fOnlyRegions.size()>0 && Common::FindInStringVector(fOnlyRegions,RemoveQuotes(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,RemoveQuotes(confSet->GetValue()))>=0) continue;
        fRegNames.push_back( CheckName(confSet->GetValue()) );
        fRegions.emplace_back( CheckName(confSet->GetValue()) );
        Region *reg;
        reg = fFitter->NewRegion(CheckName(confSet->GetValue()));
        reg->fGetChi2 = fFitter->fGetChi2;
        reg->SetVariableTitle(RemoveQuotes(confSet->Get("VariableTitle")));
        reg->SetLabel(RemoveQuotes(confSet->Get("Label")),RemoveQuotes(confSet->Get("ShortLabel")));

        std::string param = "";
        // Set axisTitle
        param = confSet->Get("YaxisTitle");
        if( param != "") reg->fYTitle = RemoveQuotes(param);

        // Set YmaxScale
        param = confSet->Get("YmaxScale");
        if(param!="") reg->fYmaxScale = atof(param.c_str());

        // Set Ymin
        param = confSet->Get("Ymin");
        if(param!="") reg->fYmin = atof(param.c_str());

        // Set Ymax
        param = confSet->Get("Ymax");
        if(param!="") reg->fYmax = atof(param.c_str());

        // Set RatioYmin
        param = confSet->Get("RatioYmin");
        if(param!="") {
            reg->fRatioYmin = atof(param.c_str());
            reg->fRatioYminPostFit = reg->fRatioYmin;
        }

        // Set RatioYmax
        param = confSet->Get("RatioYmax");
        if(param!="") {
            reg->fRatioYmax = atof(param.c_str());
            reg->fRatioYmaxPostFit = reg->fRatioYmax;
        }

        // Set RatioYminPostFit
        param = confSet->Get("RatioYminPostFit");
        if(param!="") reg->fRatioYminPostFit = atof(param.c_str());

        // Set RatioYmaxPostFit
        param = confSet->Get("RatioYmaxPostFit");
        if(param!="") reg->fRatioYmaxPostFit = atof(param.c_str());

        // Set TexLabel
        param = confSet->Get("TexLabel");
        if( param != "") reg->fTexLabel = RemoveQuotes(param);

        // Set LumiLabel
        param = confSet->Get("LumiLabel");
        if( param != "") reg->fLumiLabel = RemoveQuotes(param);
        else reg->fLumiLabel = fFitter->fLumiLabel;

        // Set CmeLabel
        param = confSet->Get("CmeLabel");
        if( param != "") reg->fCmeLabel = RemoveQuotes(param);
        else reg->fCmeLabel = fFitter->fCmeLabel;

        // Set LogScale
        param = confSet->Get("LogScale");
        if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") reg->fLogScale = true;
            else if(param=="FALSE") reg->fLogScale = false;
            else {
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "You specified 'LogScale' option but did not provide valid parameter. Using default (false)");
                reg->fLogScale = false;
            }
        }

        // Set Group
        param = confSet->Get("Group");
        if( param != "") reg->fGroup = RemoveQuotes(param);

        // Paths for the unfolding code
        param = confSet->Get("ResponseMatrixFile");
        if (param != "") {
            reg->fResponseMatrixFiles.clear();
            reg->fResponseMatrixFiles.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixFiles");
        if (param != "") {
            reg->fResponseMatrixFiles = Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("ResponseMatrixName");
        if (param != "") {
            reg->fResponseMatrixNames.clear();
            reg->fResponseMatrixNames.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNames");
        if (param != "") {
            reg->fResponseMatrixNames = Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("ResponseMatrixPath");
        if (param != "") {
            reg->fResponseMatrixPaths.clear();
            reg->fResponseMatrixPaths.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixPaths");
        if (param != "") {
            reg->fResponseMatrixPaths = Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixFileSuff");
        if (param != "") {
            reg->fResponseMatrixFileSuffs.clear();
            reg->fResponseMatrixFileSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixFileSuffs");
        if (param != "") {
            reg->fResponseMatrixFileSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixNameSuff");
        if (param != "") {
            reg->fResponseMatrixNameSuffs.clear();
            reg->fResponseMatrixNameSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNameSuffs");
        if (param != "") {
            reg->fResponseMatrixNameSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixPathSuff");
        if (param != "") {
            reg->fResponseMatrixPathSuffs.clear();
            reg->fResponseMatrixPathSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixPathSuffs");
        if (param != "") {
            reg->fResponseMatrixPathSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceName");
        if (param != "") {
            reg->fAcceptanceNames.clear();
            reg->fAcceptanceNames.emplace_back(RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceNames");
        if (param != "") {
            reg->fAcceptanceNames = Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        // Paths for the unfolding code
        param = confSet->Get("AcceptancePath");
        if (param != "") {
            reg->fAcceptancePaths.clear();
            reg->fAcceptancePaths.emplace_back(RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptancePaths");
        if (param != "") {
            reg->fAcceptancePaths = Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceFileSuff");
        if (param != "") {
            reg->fAcceptanceFileSuffs.clear();
            reg->fAcceptanceFileSuffs.emplace_back(RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceFileSuffs");
        if (param != "") {
            reg->fAcceptanceFileSuffs = Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceNameSuff");
        if (param != "") {
            reg->fAcceptanceNameSuffs.clear();
            reg->fAcceptanceNameSuffs.emplace_back(RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceNameSuffs");
        if (param != "") {
            reg->fAcceptanceNameSuffs = Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptancePathSuff");
        if (param != "") {
            reg->fAcceptancePathSuffs.clear();
            reg->fAcceptancePathSuffs.emplace_back(RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptancePathSuffs");
        if (param != "") {
            reg->fAcceptancePathSuffs = Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("SelectionEffName");
        if (param != "") {
            reg->fSelectionEffNames.clear();
            reg->fSelectionEffNames.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNames");
        if (param != "") {
            reg->fSelectionEffNames = Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("SelectionEffPath");
        if (param != "") {
            reg->fSelectionEffPaths.clear();
            reg->fSelectionEffPaths.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPaths");
        if (param != "") {
            reg->fSelectionEffPaths = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFileSuff");
        if (param != "") {
            reg->fSelectionEffFileSuffs.clear();
            reg->fSelectionEffFileSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFileSuffs");
        if (param != "") {
            reg->fSelectionEffFileSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffNameSuff");
        if (param != "") {
            reg->fSelectionEffNameSuffs.clear();
            reg->fSelectionEffNameSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNameSuffs");
        if (param != "") {
            reg->fSelectionEffNameSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPathSuff");
        if (param != "") {
            reg->fSelectionEffPathSuffs.clear();
            reg->fSelectionEffPathSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathSuffs");
        if (param != "") {
            reg->fSelectionEffPathSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationName");
        if (param != "") {
            reg->fMigrationNames.clear();
            reg->fMigrationNames.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNames");
        if (param != "") {
            reg->fMigrationNames = Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("MigrationPath");
        if (param != "") {
            reg->fMigrationPaths.clear();
            reg->fMigrationPaths.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPaths");
        if (param != "") {
            reg->fMigrationPaths = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFileSuff");
        if (param != "") {
            reg->fMigrationFileSuffs.clear();
            reg->fMigrationFileSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFileSuffs");
        if (param != "") {
            reg->fMigrationFileSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationNameSuff");
        if (param != "") {
            reg->fMigrationNameSuffs.clear();
            reg->fMigrationNameSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNameSuffs");
        if (param != "") {
            reg->fMigrationNameSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPathSuff");
        if (param != "") {
            reg->fMigrationPathSuffs.clear();
            reg->fMigrationPathSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathSuffs");
        if (param != "") {
            reg->fMigrationPathSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("NumberOfRecoBins");
        if (param != "") {
            const int bins = std::stoi(param);
            if (bins < 2 || bins > 100) {
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "Number of reco bins is < 2 or > 100. This does not seem correct");
                return 1;
            }
            reg->fNumberUnfoldingRecoBins = bins;
        } else {
            if (fFitter->fFitType == TRExFit::UNFOLDING) {
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You need to provide the number of reco bins (NumberOfRecoBins) in each Region.");
                return 1;
            }
        }

        param = confSet->Get("NormalizeMigrationMatrix");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "TRUE") {
                reg->fNormalizeMigrationMatrix = true;
            } else if (param == "FALSE") {
                reg->fNormalizeMigrationMatrix = false;
            } else {
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `NormalizeMigrationMatrix` option, but you did not provide any reasonable option. Using the default (TRUE)!");
                reg->fNormalizeMigrationMatrix = true;
            }
        }


        // Paths for the unfolding code
        // Paths for the unfolding code
        // Paths for the unfolding code
        // Setting based on input type
        if (fFitter->fInputType == 0){
            if (SetRegionHIST(reg, confSet) != 0) return 1;
        } else if (fFitter->fInputType == 1){
            if (SetRegionNTUP(reg, confSet) != 0) return 1;
        } else {
            WriteErrorStatus("ConfigReader::ReadRegionOptions", "Unknown input type: " +std::to_string(fFitter->fInputType));
            return 1;
        }

        // Set Rebin
        param = confSet->Get("Rebin");
        if(param != "") reg->Rebin(atoi(param.c_str()));

        // Set post-process rebin ("Rebinning")
        param = confSet->Get("Rebinning");
        if(param != ""){
            std::vector < std::string > vec_bins = Vectorize(param, ',');
            if (vec_bins.size() == 0){
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Rebinning` option, but you did not provide any reasonable option. Check this!");
                return 1;
            }
            // eventually add auto-binning
            const unsigned int nBounds = vec_bins.size();
            double *bins = new double[nBounds];
            for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
                bins[iBound] = atof(vec_bins[iBound].c_str());
            }
            reg -> SetRebinning(nBounds-1,bins);
        }

        // Set Binning
        param = confSet->Get("Binning");
        if(param != "" && param !="-"){
            std::vector < std::string > vec_bins = Vectorize(param, ',');
            if (vec_bins.size() == 0){
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option, but you did not provide any reasonable option. Check this!");
                return 1;
            }
            if(vec_bins[0]=="AutoBin"){
                if (vec_bins.size() < 2){
                    WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option with Autobin, but you did not provide any reasonable option. Check this!");
                    return 1;
                }
                reg -> fBinTransfo = vec_bins[1];
                if(vec_bins[1]=="TransfoD"){
                    if (vec_bins.size() < 4){
                        WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option with TransfoD, but you did not provide any reasonable option. Check this!");
                        return 1;
                    }
                    reg -> fTransfoDzSig=Common::convertStoD(vec_bins[2]);
                    reg -> fTransfoDzBkg=Common::convertStoD(vec_bins[3]);
                    if(vec_bins.size()>4){
                        for(const std::string& ibkg : vec_bins){
                            reg -> fAutoBinBkgsInSig.push_back(ibkg);
                        }
                    }
                }
                else if(vec_bins[1]=="TransfoF"){
                    if (vec_bins.size() < 4){
                        WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option with TransfoF, but you did not provide any reasonable option. Check this!");
                        return 1;
                    }
                    reg -> fTransfoFzSig=Common::convertStoD(vec_bins[2]);
                    reg -> fTransfoFzBkg=Common::convertStoD(vec_bins[3]);
                    if(vec_bins.size()>4){
                        for(const std::string& ibkg : vec_bins){
                            reg -> fAutoBinBkgsInSig.push_back(ibkg);
                        }
                    }
                }
                else if(vec_bins[1]=="TransfoJ"){
                    if(vec_bins.size() > 2) reg -> fTransfoJpar1=Common::convertStoD(vec_bins[2]);
                    else reg -> fTransfoJpar1 = 5.;
                    if(vec_bins.size() > 3) reg -> fTransfoJpar2=Common::convertStoD(vec_bins[3]);
                    else reg -> fTransfoJpar2 = 1.;
                    if(vec_bins.size() > 4) reg -> fTransfoJpar3=Common::convertStoD(vec_bins[4]);
                    else reg -> fTransfoJpar3 = 5.;
                    if(vec_bins.size()>5){
                        for(const std::string& ibkg : vec_bins){
                            reg -> fAutoBinBkgsInSig.push_back(ibkg);
                        }
                    }
                }
                else{
                    WriteErrorStatus("ConfigReader::ReadRegionOptions", "Unknown transformation: " + vec_bins[1] + ", try again");
                    return 1;
                }
            }
            else{
                const unsigned int nBounds = vec_bins.size();
                double *bins = new double[nBounds];
                for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
                    bins[iBound] = atof(vec_bins[iBound].c_str());
                }
                reg -> SetBinning(nBounds-1,bins);
            }
        }

        // Set BinWidth
        param = confSet->Get("BinWidth");
        if(param != "") reg->fBinWidth = atof(param.c_str());

        // Set Type
        param = confSet->Get("Type");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="CONTROL" )     reg -> SetRegionType(Region::CONTROL);
            else if( param=="VALIDATION" )  reg -> SetRegionType(Region::VALIDATION);
            else if( param=="SIGNAL" )      reg -> SetRegionType(Region::SIGNAL);
            else {
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified 'Type' option in region but did not provide valid parameter. Please check this!");
                return 1;
            }
        }
        if(reg -> fRegionType != Region::VALIDATION){
            reg->fUseGammaPulls = fFitter->fUseGammaPulls;
            fHasAtLeastOneValidRegion = true;
        }

        // Set DataType
        param = confSet->Get("DataType");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="DATA" )     reg -> SetRegionDataType(Region::REALDATA);
            else if( param=="ASIMOV" )  reg -> SetRegionDataType(Region::ASIMOVDATA);
            else{
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "DataType is not recognised: " + param);
                return 1;
            }
        }

        // Set SkipSmoothing
        param = confSet->Get("SkipSmoothing");
        if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE")  reg->fSkipSmoothing = true;
            else if(param=="FALSE")  reg->fSkipSmoothing = false;
            else {
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "You specified 'SkipSmoothing' option in region but did not provide valid parameter. Using default (FALSE)");
                reg->fSkipSmoothing = false;
            }
        }

        // Set AutomaticDropBins
        param = confSet->Get("AutomaticDropBins");
        if (param != "") {
            bool isOK(true);
            if (reg->fDropBins.size() > 0) {
                WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'AutomaticDropBins' option but you previously set DropBins for region " + reg->fName + " ignoring the automatic option");
                isOK = false;
            }
            if (isOK) {
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if (param == "TRUE") {
                    reg->SetAutomaticDropBins(true);
                } else if (param == "FALSE") {
                    reg->SetAutomaticDropBins(false);
                } else {
                    WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'AutomaticDropBins' option but did not provide a valid setting. Using default (FALSE)");
                    reg->SetAutomaticDropBins(false);
                }
            }
        }

        // Set DropBins
        param = confSet->Get("DropBins");
        if( param != "" ){
            if (reg->GetAutomaticDropBins()) {
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "You specified set `AutomaticDropBins` to TRUE, but using DropBins will disable it for region " + reg->fName + "!.");
            }
            reg->SetAutomaticDropBins(false);
            reg->fDropBins.clear();
            const std::vector<std::string>& s = Vectorize( param,',' );
            for(const std::string& is : s){
                reg->fDropBins.push_back(atoi(is.c_str()));
            }
        }

        // Set BinLabels
        param = confSet->Get("BinLabels");
        if( param != "" ){
            std::vector<std::string> vec_string = Vectorize( param,',' );
            reg->fBinLabels = vec_string;
        }

        // Set XaxisRange
        param = confSet->Get("XaxisRange");
        if( param != "" ){
            std::vector<std::string> vec_string = Vectorize( param,',' );
            if (vec_string.size() != 2){
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "Setting 'XaxisRange' needs exactly two parameters (floats). Ignoring.");
            }
            const double min = std::stod(vec_string.at(0));
            const double max = std::stod(vec_string.at(1));

            std::vector<double> range{};
            if (min < max){
                range.emplace_back(min);
                range.emplace_back(max);
            } else {
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "Setting 'XaxisRange' needs the first parameter to be smaller than the second parameter. Ignoring.");
            }
            reg->fXaxisRange = range;
        }

        // Inter-region smoothing
        param = confSet->Get("IsBinOfRegion");
        if( param != "" ){
            std::vector<std::string> vec_string = Vectorize( param,':' );
            if (vec_string.size() != 2){
                WriteWarningStatus("ConfigReader::IsBinOfRegion", "Setting 'IsBinOfRegion' needs exactly two parameters (in the form string:int). Ignoring.");
            }
            else{
                reg->fIsBinOfRegion[vec_string.at(0)] = std::stof(vec_string.at(1));
            }
        }

    }

    if (!fHasAtLeastOneValidRegion && Common::OptionRunsFit(opt)){
        WriteErrorStatus("ConfigReader::ReadRegionOptions","You need to provide at least one region that is not Validation otherwise the fit will crash.");
        return 1;
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetRegionHIST(Region* reg, ConfigSet *confSet){
    std::string param = "";

    // Set HistoFile
    param = confSet->Get("HistoFile");
    if(param!=""){
        reg->fHistoFiles.clear();
        reg->fHistoFiles.push_back( RemoveQuotes(param) );
    }
    // Set HistoFiles
    param = confSet->Get("HistoFiles");
    if(param!="") reg->fHistoFiles = Vectorize( param,',' );

    // Set HistoName
    param = confSet->Get("HistoName");
    if(param!=""){
        reg->fHistoNames.clear();
        reg->fHistoNames.push_back( RemoveQuotes(param) );
    }
    // Set HistoNames
    param = confSet->Get("HistoNames");
    if(param!="") reg->fHistoNames = Vectorize( param,',' );

    // Set HistoPath
    param = confSet->Get("HistoPath");
    if(param!=""){
        reg->fHistoPaths.clear();
        reg->fHistoPaths.push_back( RemoveQuotes(param) );
    }
    // Set HistoFiles
    param = confSet->Get("HistoPaths");
    if(param!="") reg->fHistoPaths = Vectorize( param,',' );

    // Set HistoFileSuff
    param = confSet->Get("HistoFileSuff");
    if(param !=""){
        reg->fHistoFileSuffs.clear();
        reg->fHistoFileSuffs.push_back( RemoveQuotes(param) );
    }
    // Set HistoPathSuffs
    param = confSet->Get("HistoFileSuffs");
    if(param!=""){
        reg->fHistoFileSuffs.clear();
        std::vector<std::string> paths = Vectorize( param,',' );
        for(std::string ipath : paths){
            reg->fHistoFileSuffs.push_back( RemoveQuotes(ipath) );
        }
    }

    // Set HistoNameSuff
    param = confSet->Get("HistoNameSuff");
    if(param !=""){
        reg->fHistoNameSuffs.clear();
        reg->fHistoNameSuffs.push_back( RemoveQuotes(param) );
    }
    // Set HistoNameSuffs
    param = confSet->Get("HistoNameSuffs");
    if(param!=""){
        reg->fHistoNameSuffs.clear();
        std::vector<std::string> paths = Vectorize( param,',' );
        for(std::string ipath : paths){
            reg->fHistoNameSuffs.push_back( RemoveQuotes(ipath) );
        }
    }

    // Set HistoPathSuff
    param = confSet->Get("HistoPathSuff");
    if(param !=""){
        reg->fHistoPathSuffs.clear();
        reg->fHistoPathSuffs.push_back( RemoveQuotes(param) );
    }
    // Set HistoPathSuffs
    param = confSet->Get("HistoPathSuffs");
    if(param!=""){
        reg->fHistoPathSuffs.clear();
        std::vector<std::string> paths = Vectorize( param,',' );
        for(std::string ipath : paths){
            reg->fHistoPathSuffs.push_back( RemoveQuotes(ipath) );
        }
    }

    // Check for NTUP inputs
    if (ConfigHasNTUP(confSet)){
        WriteWarningStatus("ConfigReader::SetRegionHIST", "Found some NTUP settings in Region, while input option is HIST. Ignoring them.");
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetRegionNTUP(Region* reg, ConfigSet *confSet){
    std::string param = "";

    // Set Variable
    std::vector<std::string> variable = Vectorize(confSet->Get("Variable"),',');
    if (variable.size() == 0){
        WriteErrorStatus("ConfigReader::SetRegionNTUP", "Variable option is required but not present. Please check it!");
        return 1;

    }
    // fix variable vector if special functions are used
    if(variable[0].find("Alt$")!=std::string::npos || variable[0].find("MaxIf$")!=std::string::npos ||variable[0].find("MinIf$")!=std::string::npos ){
        if (variable.size() > 1){
            do{
            variable[0]+=","+variable[1];
            variable.erase(variable.begin()+1);
            }
            while(variable[1].find("Alt$")!=std::string::npos);
            variable[0]+=","+variable[1];
            variable.erase(variable.begin()+1);
        } else {
            WriteErrorStatus("ConfigReader::SetRegionNTUP", "Variable option has weird input. Please check it!");
            return 1;
        }
    }

    std::vector<std::string> corrVar  = Vectorize(variable[0],'|');
    if(corrVar.size()==2){
        if (variable.size() < 3){
            WriteErrorStatus("ConfigReader::SetRegionNTUP", "Corr size == 2 but variable size < 3. Check this!");
            return 1;
        }
        WriteDebugStatus("ConfigReader::SetRegionNTUP", "Have a correlation variable in reg " + fRegNames.back() + " : ");
        WriteDebugStatus("ConfigReader::SetRegionNTUP", corrVar[0] + " and " + corrVar[1]);
        reg->SetVariable(  "corr_"+corrVar[0]+"_"+corrVar[1], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()), corrVar[0].c_str(), corrVar[1].c_str() );
    }
    else {
        if (variable.size() < 4){
            WriteErrorStatus("ConfigReader::SetRegionNTUP", "Corr size != 2 but variable size < 4. Check this!");
            return 1;
        }
        WriteDebugStatus("ConfigReader::SetRegionNTUP", "Have a usual variable in reg " + fRegNames.back() + " : ");
        WriteDebugStatus("ConfigReader::SetRegionNTUP", variable[0] + " and size of corrVar=" + std::to_string(corrVar.size()));
        reg->SetVariable(  variable[0], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()) );
    }

    // Set VariableForSample
    param = confSet->Get("VariableForSample");
    if( param != "" ){
        std::vector < std::string > temp_samplesAndVars = Vectorize(param,',',false);
        for(std::string ivar : temp_samplesAndVars){
          std::vector < std::string > vars = Vectorize(ivar,':');
            if(vars.size()==2){
                reg->SetAlternativeVariable(vars[1], vars[0]);
            }
        }
    }

    // Set Selection
    param = confSet->Get("Selection");
    if(param != "") reg->AddSelection( RemoveQuotes(param) );

    // Set SelectionForSample
    param = confSet->Get("SelectionForSample");
    if( param != "" ){
        std::vector < std::string > temp_samplesAndSels = Vectorize(param,',',false);
        for(std::string ivar : temp_samplesAndSels){
          std::vector < std::string > vars = Vectorize(ivar,':');
            if(vars.size()==2){
                reg->SetAlternativeSelection(vars[1], vars[0]);
            }
        }
    }

    // Set MCweight
    param = confSet->Get("MCweight");
    if (param != "") reg->fMCweight = RemoveQuotes(param); // this will override the global MCweight, if any

    // Set NtupleFile
    param = confSet->Get("NtupleFile");
    if(param!=""){
        reg->fNtupleFiles.clear();
        reg->fNtupleFiles.push_back( RemoveQuotes(param) );
    }
    // Set NtupleFiles
    param = confSet->Get("NtupleFiles");
    if(param!="") reg->fNtupleFiles = Vectorize( param,',' );

    // Set NtupleFileSuff
    param = confSet->Get("NtupleFileSuff");
    if(param!="") {
        reg->fNtupleFileSuffs.clear();
        reg->fNtupleFileSuffs.push_back( RemoveQuotes(param) );
    }
    // Set NtupleFileSuffs
    param = confSet->Get("NtupleFileSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Vectorize( param,',' );
        reg->fNtupleFileSuffs = paths;
    }

    // Set NtupleName
    param = confSet->Get("NtupleName");
    if(param!="") {
        reg->fNtupleNames.clear();
        reg->fNtupleNames.push_back( RemoveQuotes(param) );
    }
    // Set NtupleNames
    param = confSet->Get("NtupleNames");
    if(param!="") reg->fNtupleNames = Vectorize( param,',' );

    // Set NtupleNameSuff
    param = confSet->Get("NtupleNameSuff");
    if(param!="") {
        reg->fNtupleNameSuffs.clear();
        reg->fNtupleNameSuffs.push_back( RemoveQuotes(param) );
    }
    // Set NtupleNameSuffs
    param = confSet->Get("NtupleNameSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Vectorize( param,',' );
        reg->fNtupleNameSuffs = paths;
    }

    // Set NtuplePath
    param = confSet->Get("NtuplePath");
    if(param!="") {
        reg->fNtuplePaths.clear();
        reg->fNtuplePaths.push_back( RemoveQuotes(param) );
    }
    // Set NtuplePaths
    param = confSet->Get("NtuplePaths");
    if(param!="") reg->fNtuplePaths = Vectorize( param,',' );

    // Set NtuplePathSuff
    param = confSet->Get("NtuplePathSuff");
    if(param != "") {
        reg->fNtuplePathSuffs.clear();
        reg->fNtuplePathSuffs.push_back( RemoveQuotes(param) );
    }
    // Set NtuplePathSuffs
    param = confSet->Get("NtuplePathSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Vectorize( param,',' );
        reg->fNtuplePathSuffs = paths;
    }

    // Check for NTUP inputs
    if (ConfigHasHIST(confSet)){
        WriteWarningStatus("ConfigReader::SetRegionNTUP", "Found some HIST settings in Region, while input option is NTUP. Ignoring them.");
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadSampleOptions() {

    fAvailableSamples = GetAvailableSamples();

    // Check of the Samples from commands like exist
    if (fOnlyRegions.size() > 0){
        if (!CheckPresence(fOnlySamples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You set samples that do not exist in your command line options");
            } else {
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "You set samples that do not exist in your command line options");
                return 1;
            }
        }
    }

    fHasAtLeastOneValidSample = false;

    int nSmp = 0;
    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("Sample",nSmp);
        if (confSet == nullptr) break;
        nSmp++;

        Sample *sample = nullptr;
        NormFactor *nfactor = nullptr;
        ShapeFactor *sfactor = nullptr;
        int type = 0;
        std::string param = "";

        if(fOnlySamples.size()>0 && Common::FindInStringVector(fOnlySamples,RemoveQuotes(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,RemoveQuotes(confSet->GetValue()))>=0) continue;
        type = Sample::BACKGROUND;

        // Set Type
        param = confSet->Get("Type");
        if (param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "SIGNAL"){
                type = Sample::SIGNAL;
                fNonGhostIsSet = true;
                fHasAtLeastOneValidSample = true;
            }
            else if(param == "DATA"){
                type = Sample::DATA;
                fNonGhostIsSet = true;
            }
            else if(param == "GHOST"){
                if (fNonGhostIsSet){
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Please define GHOST samples first and then other samples");
                    return 1;
                }
                type = Sample::GHOST;
            }
            else if(param == "BACKGROUND"){
                type = Sample::BACKGROUND;
                fNonGhostIsSet = true;
                fHasAtLeastOneValidSample = true;
            }
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'Type' option in sample but did not provide valid parameter. Using default (BACKGROUND)");
            }
            if(fOnlySignals.size()>0){
                if(type==Sample::SIGNAL){
                    bool skip = true;
                    for(const auto& s : fOnlySignals){
                        if(CheckName(confSet->GetValue())==s) skip = false;
                    }
                    if(skip) continue;
                }
                else if(type==Sample::GHOST){
                    for(const auto& s : fOnlySignals){
                        if(CheckName(confSet->GetValue())==s) type = Sample::SIGNAL;
                    }
                }
            }
        }
        sample = fFitter->NewSample(CheckName(confSet->GetValue()),type);
        fSamples.emplace_back(CheckName(confSet->GetValue()));

        // Set Title
        param = confSet->Get("Title");
        if (param != "") sample->SetTitle(RemoveQuotes(param));

        // Set TexTitle
        param = confSet->Get("TexTitle");
        if(param!="") sample->fTexTitle = RemoveQuotes(param);

        // Set Group
        param = confSet->Get("Group");
        if(param!="") sample->fGroup = RemoveQuotes(param);

        // HIST input
        if (fFitter->fInputType == 0){
            // Set HistoFile
            param = confSet->Get("HistoFile");
            if(param!="") sample->fHistoFiles.push_back( RemoveQuotes(param) );

            // Set HistoFiles
            param = confSet->Get("HistoFiles");
            if(param!="") sample->fHistoFiles = Vectorize( param, ',' );

            // Set HistoFileSuff
            param = confSet->Get("HistoFileSuff");
            if(param!="") sample->fHistoFileSuffs.push_back( RemoveQuotes(param) );

            // Set HistoFileSuffs
            param = confSet->Get("HistoFileSuffs");
            if(param!="") sample->fHistoFileSuffs = Vectorize( param, ',' );

            // Set HistoName
            param = confSet->Get("HistoName");
            if(param!="") sample->fHistoNames.push_back( RemoveQuotes(param) );

            // Set HistoNames
            param = confSet->Get("HistoNames");
            if(param!="") sample->fHistoNames = Vectorize( param, ',' );

            // Set HistoNameSuff
            param = confSet->Get("HistoNameSuff");
            if(param!="") sample->fHistoNameSuffs.push_back( RemoveQuotes(param) );

            // Set HistoNameSuffs
            param = confSet->Get("HistoNameSuffs");
            if(param!="") sample->fHistoNameSuffs = Vectorize( param, ',' );

            // Set HistoPath
            param = confSet->Get("HistoPath");
            if(param!="") sample->fHistoPaths.push_back( RemoveQuotes(param) );

            // Set HistoPaths
            param = confSet->Get("HistoPaths");
            if(param!="") sample->fHistoPaths = Vectorize( param, ',' );

            // Set HistoPathSuff
            param = confSet->Get("HistoPathSuff");
            if(param!="") sample->fHistoPathSuffs.push_back( RemoveQuotes(param) );

            // Set HistoPathSuffs
            param = confSet->Get("HistoPathSuffs");
            if(param!="") sample->fHistoPathSuffs = Vectorize( param, ',' );

            if (ConfigHasNTUP(confSet)){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You provided some NTUP options but your input type is HIST. Options will be ignored");
            }
        } else if (fFitter->fInputType == 1){ // NTUP input
            // Set NtupleFile
            param = confSet->Get("NtupleFile");
            if(param!="") sample->fNtupleFiles.push_back( RemoveQuotes(param) );

            // Set NtupleFiles
            param = confSet->Get("NtupleFiles");
            if(param!="") sample->fNtupleFiles = Vectorize( param ,',' );

            param = confSet->Get("NtupleFileSuff");
            if(param!="") sample->fNtupleFileSuffs.push_back( RemoveQuotes(param) );

            // Set NtupleFileSuffs
            param = confSet->Get("NtupleFileSuffs");
            if(param!="") sample->fNtupleFileSuffs = Vectorize( param ,',' );

            // Set NtupleName
            param = confSet->Get("NtupleName");
            if(param!="") sample->fNtupleNames.push_back( RemoveQuotes(param) );

            // Set NtupleNames
            param = confSet->Get("NtupleNames");
            if(param!="") sample->fNtupleNames = Vectorize( param ,',' );

            // Set NtupleNameSuff
            param = confSet->Get("NtupleNameSuff");
            if(param!="") sample->fNtupleNameSuffs.push_back( RemoveQuotes(param) );

            // Set NtupleNameSuffs
            param = confSet->Get("NtupleNameSuffs");
            if( param != "" ) sample->fNtupleNameSuffs = Vectorize( param,',' );

            // Set NtuplePath
            param = confSet->Get("NtuplePath");
            if(param!="") sample->fNtuplePaths.push_back( RemoveQuotes(param) );

            // Set NtuplePaths
            param = confSet->Get("NtuplePaths");
            if(param != "") sample->fNtuplePaths = Vectorize( param ,',' );

            // Set NtuplePathSuff
            param = confSet->Get("NtuplePathSuff");
            if(param!="") sample->fNtuplePathSuffs.push_back( RemoveQuotes(param) );

            // Set NtuplePathSuffs
            param = confSet->Get("NtuplePathSuffs");
            if(param != "") sample->fNtuplePathSuffs = Vectorize( param ,',' );

            if (ConfigHasHIST(confSet)){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You provided some HIST options but your input type is NTUP. Options will be ignored");
            }
        } else {
            WriteErrorStatus("ConfigReader::ReadSampleOptions", "No valid input type provided. Please check this!");
            return 1;
        }

        // Set FillColor
        param = confSet->Get("FillColor");
        if(param != "") sample->SetFillColor(atoi(param.c_str()));

        // Set LineColor
        param = confSet->Get("LineColor");
        if(param != "") sample->SetLineColor(atoi(param.c_str()));

        // Convert a string to an RGB value
        auto convert_str_to_RGB = [] (const std::string& str) {
          int num = std::stoi(str);
          if (num < 0) throw std::invalid_argument("RGB value out of range [0, 255]");
          if (num > 255) throw std::invalid_argument("RGB value out of range [0, 255]");
          return num;
        };

        // Convert a vector of strings to a vector of RGB values
        auto create_RGB_array = [&convert_str_to_RGB] (const std::vector<std::string>& str_vec) {
          std::array<int, 3> int_arr{{-1}};
          for (auto itr = str_vec.begin(); itr != str_vec.end(); ++itr) {
            int_arr.at(itr - str_vec.begin()) = convert_str_to_RGB(*itr);
          }
          return int_arr;
        };

        // Set FillColor from RGB values if given
        param = confSet->Get("FillColorRGB");
        if (param != "") {
          std::vector<std::string> rgb_strings = Vectorize(param, ',');
          if (rgb_strings.size() != 3) {
            WriteErrorStatus("ConfigReader::ReadSampleOptions", "No valid input for 'FillColorRGB' provided. Please check this!");
            return 1;
          } else {
            auto col_arr = create_RGB_array(rgb_strings);
            TColor* tcol = new TColor{TColor::GetFreeColorIndex(),
                                      (float)col_arr[0]/255,
                                      (float)col_arr[1]/255,
                                      (float)col_arr[2]/255};
            sample->fFillColor = tcol->GetNumber();
          }
        }

        // Set LineColor from RGB values if given
        param = confSet->Get("LineColorRGB");
        if (param != "") {
          std::vector<std::string> rgb_strings = Vectorize(param, ',');
          if (rgb_strings.size() != 3) {
            WriteErrorStatus("ConfigReader::ReadSampleOptions", "No valid input for 'LineColorRGB' provided. Please check this!");
            return 1;
          } else {
            auto col_arr = create_RGB_array(rgb_strings);
            TColor* tcol = new TColor{TColor::GetFreeColorIndex(),
                                      (float)col_arr[0]/255,
                                      (float)col_arr[1]/255,
                                      (float)col_arr[2]/255};
            sample->fLineColor = tcol->GetNumber();
          }
        }

        // Set NormalizedByTheory
        param = confSet->Get("NormalizedByTheory");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") sample->NormalizedByTheory(false);
            else if(param=="TRUE") sample->NormalizedByTheory(true);
            else{
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'NormalizedByTheory' option but did not provide valid parameter. Using default (true)");
                sample->NormalizedByTheory(true);
            }
        }

        // Set MCweight and Selection
        if(fFitter->fInputType==1){
            param = confSet->Get("MCweight");
            if(param != "") sample->SetMCweight( RemoveQuotes(param) );

            param = confSet->Get("Selection");
            if(param!="") sample->SetSelection( RemoveQuotes(param) );
        }

        // to specify only certain regions
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');
        sample->fRegions.clear();

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
                return 1;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has regions to exclude set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has regions to exclude set up that do not exist");
                return 1;
            }
        }

        for(int i_reg=0;i_reg<fFitter->fNRegions;i_reg++){
            std::string regName = fFitter->fRegions[i_reg]->fName;
            if( (regions_str=="" || regions_str=="all" || Common::FindInStringVector(regions,regName)>=0)
                && Common::FindInStringVector(exclude,regName)<0 ){
                sample->fRegions.push_back( fFitter->fRegions[i_reg]->fName );
            }
        }

        // Set NormFactor
        param = confSet->Get("NormFactor");
        if(param!=""){
            // check if the normfactor is called just with the name or with full definition
            const unsigned int sz = Vectorize(param,',').size();
            if (sz != 1 && sz != 4 && sz != 5){
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "No valid input for 'NormFactor' provided. Please check this!");
                return 1;
            }
            if( sz > 1 ){
                bool isConst = false;
                if( Vectorize(param,',').size()>4){
                    std::string tmp = Vectorize(param,',')[4];
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                    if (tmp == "TRUE"){
                        isConst = true;
                    }
                }
                if (sz > 3)
                    nfactor = sample->AddNormFactor(
                    Vectorize(param,',')[0],
                    atof(Vectorize(param,',')[1].c_str()),
                    atof(Vectorize(param,',')[2].c_str()),
                    atof(Vectorize(param,',')[3].c_str()),
                    isConst
                );
            }
            else{
                nfactor = sample->AddNormFactor( Vectorize(param,',')[0] );
            }
            nfactor->fRegions = sample->fRegions;
            if( Common::FindInStringVector(fFitter->fNormFactorNames,nfactor->fName)<0 ){
                fFitter->fNormFactors.push_back( nfactor );
                fFitter->fNormFactorNames.push_back( nfactor->fName );
                fFitter->fNNorm++;
            }
        }

        // Set ShapeFactor
        param = confSet->Get("ShapeFactor");
        if(param!=""){
            // check if the normfactor is called just with the name or with full definition
            const unsigned int sz = Vectorize(param,',').size();
            if (sz != 1 && sz != 4 && sz != 5){
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "No valid input for 'ShapeFactor' provided. Please check this!");
                return 1;
            }
            if( sz > 1 ){
                bool isConst = false;
                if( Vectorize(param,',').size()>4){
                    std::string tmp = Vectorize(param,',')[4];
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                    if (tmp == "TRUE"){
                        isConst = true;
                    }
                }
                if (sz > 3)
                    sfactor = sample->AddShapeFactor(
                    Vectorize(param,',')[0],
                    atof(Vectorize(param,',')[1].c_str()),
                    atof(Vectorize(param,',')[2].c_str()),
                    atof(Vectorize(param,',')[3].c_str()),
                    isConst
                );
            }
            else{
                sfactor = sample->AddShapeFactor( Vectorize(param,',')[0] );
            }
            sfactor->fRegions = sample->fRegions;
            if( Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
                fFitter->fShapeFactors.push_back( sfactor );
                fFitter->fShapeFactorNames.push_back( sfactor->fName );
                fFitter->fNShape++;
            }
        }

        // Set LumiScale
        param = confSet->Get("LumiScale");
        if(param!="") sample->fLumiScales.push_back( atof(param.c_str()) );

        // Set LumiScales
        param = confSet->Get("LumiScales");
        if(param!=""){
            std::vector<std::string> lumiScales_str = Vectorize( param ,',' );
            for(const std::string& ilumiscale : lumiScales_str){
                sample->fLumiScales.push_back( atof(ilumiscale.c_str()) );
            }
        }

        // Set IgnoreSelection
        // to skip global & region selection for this sample
        param = confSet->Get("IgnoreSelection");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE") sample->fIgnoreSelection = "TRUE";
            else if(param == "FALSE") sample->fIgnoreSelection = "FALSE";
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'IgnoreSelection' option but did not provide valid parameter. Using default (false)");
                sample->fIgnoreSelection = "FALSE";
            }
        }

        // Set IgnoreWeights
        // to skip global & region weights for this sample
        param = confSet->Get("IgnoreWeight");
        if(param!=""){
            sample->fIgnoreWeight = RemoveQuotes(param);
        }

        // Set UseMCstat
        // to skip MC stat uncertainty for this sample
        param = confSet->Get("UseMCstat");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") sample->fUseMCStat = false;
            else if (param == "TRUE") sample->fUseMCStat = true;
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'UseMCstat' option but did not provide valid parameter. Using default (true)");
                sample->fUseMCStat = true;
            }
        }

        // Set UseSystematics
        // to skip MC systematics for this sample
        param = confSet->Get("UseSystematics");
        // set it to false for ghost samples and data and true for other samples
        if(type == Sample::GHOST || type == Sample::DATA) sample->fUseSystematics = false;
        else                                              sample->fUseSystematics = true;
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") sample->fUseSystematics = false;
            else if(param == "TRUE" ) sample->fUseSystematics = true;
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'UseSystematics' option but did not provide valid parameter. Using default (true)");
                sample->fUseSystematics = true;
            }
        }

        // Set DivideBy
        param = confSet->Get("DivideBy");;
        if (param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for DivideBy that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for DivideBy that do not exist");
                    return 1;
                }
            }
            sample->fDivideBy = RemoveQuotes(param);
        }

        // Set MultiplyBy
        param = confSet->Get("MultiplyBy");
        if (param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for MultiplyBy that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for MultiplyBy that do not exist");
                    return 1;
                }
            }
            sample->fMultiplyBy = RemoveQuotes(param);
        }

        // Set SubtractSample
        param = confSet->Get("SubtractSample");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for SubtractSample that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for SubtractSample that do not exist");
                    return 1;
                }
            }
            sample->fSubtractSamples.push_back( RemoveQuotes(param) );
        }

        // Set SubtractSamples
        param = confSet->Get("SubtractSamples");
        if(param != ""){
            std::vector<std::string> tmp = Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for SubtractSamples that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for SubtractSamples that do not exist");
                    return 1;
                }
            }
            sample->fSubtractSamples = Vectorize(param,',');
        }

        // Set AddSample
        param = confSet->Get("AddSample");
        if(param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for AddSample that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for AddSample that do not exist");
                    return 1;
                }
            }
            sample->fAddSamples.push_back( RemoveQuotes(param) );
        }

        // Set AddSamples
        param = confSet->Get("AddSamples");
        if(param!=""){
            std::vector<std::string> tmp = Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for AddSamples that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for AddSamples that do not exist");
                    return 1;
                }
            }
            sample->fAddSamples = Vectorize(param,',');
        }

        // Set NormToSample
        param = confSet->Get("NormToSample");
        if(param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for NormToSample that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + CheckName(confSet->GetValue()) + " has samples set up for NormToSample that do not exist");
                    return 1;
                }
            }
            sample->fNormToSample = RemoveQuotes(param);
        }

        // Set BuildPullTable
        // enable pull tables
        param = confSet->Get("BuildPullTable");
        if( param != "" ){ // can be TRUE, NORM-ONLY, NORM+SHAPE (TRUE is equal to NORM-ONLY)
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param == "TRUE" ){
                sample->fBuildPullTable = 1;
                fFitter->fWithPullTables = true;
            }
            if( param.find("NORM-ONLY")!=std::string::npos ){
                sample->fBuildPullTable = 1;
                fFitter->fWithPullTables = true;
            }
            else if( param.find("NORM+SHAPE")!=std::string::npos ){
                sample->fBuildPullTable = 2;
                fFitter->fWithPullTables = true;
            }
        }

        // Set Smooth
        // allow smoothing of nominal histogram?
        param = confSet->Get("Smooth");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE") sample->fSmooth = true;
            else if(param == "FALSE") sample->fSmooth = false;
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'Smooth' option but did not provide valid parameter. Using default (false)");
                sample->fSmooth = false;
            }
        }

        // Set AsimovReplacementFor
        // AsimovReplacementFor
        param = confSet->Get("AsimovReplacementFor");
        if(RemoveQuotes(param) != ""){
            if(Vectorize(param,',').size() == 2){
                sample->fAsimovReplacementFor.first  = Vectorize(param,',')[0];
                std::string tmp = Vectorize(param,',')[1];
                sample->fAsimovReplacementFor.second = tmp;
            } else {
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "You specified 'AsimovReplacementFor' option but did not provide 2 parameters. Please check this");
                return 1;
            }
        }

        // Set SeparateGammas
        // separate gammas
        param = confSet->Get("SeparateGammas");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE"){
                sample->fSeparateGammas = true;
                if(confSet->Get("UseMCstat") == "") sample->fUseMCStat = false; // remove the usual gammas for this sample (only if no UseMCstat is specified!!)
            } else if (param == "FALSE") sample->fSeparateGammas = false;
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'SeparateGammas' option but did not provide valid parameter. Using default (false)");
                sample->fSeparateGammas = false;
            }
        }

        // Scale up/down MC stat size
        param = confSet->Get("MCstatScale");
        if(param != ""){
            sample->fMCstatScale = atof(param.c_str());
        }

        // Set CorrelateGammasInRegions
        // in the form    CorrelateGammasInRegions: SR1:SR2,CR1:CR2:CR3
        param = confSet->Get("CorrelateGammasInRegions");
        if(param != ""){
            std::vector<std::string> sets = Vectorize(param,',',false);
            for(std::string set : sets){
                std::vector<std::string> regions_corr = Vectorize(set,':');
                WriteDebugStatus("ConfigReader::ReadSampleOptions", "Correlating gammas for this sample in regions " + set);
                sample->fCorrelateGammasInRegions.push_back(regions_corr);
            }
        }

        // Set CorrelateGammasWithSample
        param = confSet->Get("CorrelateGammasWithSample");
        if(param != ""){
            sample->fCorrelateGammasWithSample = RemoveQuotes(param);
        }

        // Set SystFromSample
        param = confSet->Get("SystFromSample");
        if(param != ""){
            sample->fSystFromSample = RemoveQuotes(param);
        }

        // Set Morphing
        param = confSet->Get("Morphing");
        if(param != ""){
            if(fFitter->fExcludeFromMorphing!=sample->fName){
                std::vector<std::string> morph_par = Vectorize(param,',');
                if (morph_par.size() != 2){
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Morphing requires exactly 2 parameters, but " + std::to_string(morph_par.size()) + " provided");
                    return 1;
                }
                std::string name      = morph_par.at(0);
                double value = std::stod(morph_par.at(1));
                WriteDebugStatus("ConfigReader::ReadSampleOptions", "Morphing: Adding " + name + ", with value: " + std::to_string(value));
                if (!fFitter->MorphIsAlreadyPresent(name, value)) fFitter->AddTemplateWeight(name, value);
                // set proper normalization
                std::string morphName = "morph_"+name+"_"+Common::ReplaceString(std::to_string(value),"-","m");
                NormFactor *nf = sample->AddNormFactor(morphName, 1, 0, 10, false);
                fFitter->fNormFactors.push_back( nf );
                fFitter->fNormFactorNames.push_back( nf->fName );
                fFitter->fNNorm++;
                sample->fIsMorph[name] = true;
                sample->fMorphValue[name] = value;
                if(Common::FindInStringVector(fFitter->fMorphParams,name)<0) fFitter->fMorphParams.push_back( name );
            }
        }

    }

    // build new samples if AsimovReplacementFor are specified
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        if(fFitter->fSamples[i_smp]->fAsimovReplacementFor.first!=""){
            if (std::find(fSamples.begin(), fSamples.end(),fFitter->fSamples[i_smp]->fAsimovReplacementFor.second ) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSampleOptions", "Sample: " + fSamples[i_smp] + " has sample set up for AsimovReplacementFor that does not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSampleOptions", "Sample: " + fSamples[i_smp] + " has sample set up for AsimovReplacementFor that does not exist");
                    return 1;
                }
            }
            WriteDebugStatus("ConfigReader::ReadSampleOptions", "Creating sample " + fFitter->fSamples[i_smp]->fAsimovReplacementFor.first);
            Sample *ca = fFitter->NewSample("customAsimov_"+fFitter->fSamples[i_smp]->fAsimovReplacementFor.first,Sample::GHOST);
            ca->SetTitle("Pseudo-Data ("+fFitter->fSamples[i_smp]->fAsimovReplacementFor.first+")");
            ca->fUseSystematics = false;
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadNormFactorOptions(){
    std::string param = "";

    int nNorm = 0;
    NormFactor *nfactor = nullptr;
    Sample *sample = nullptr;

    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("NormFactor", nNorm);
        if (confSet == nullptr) break;
        nNorm++;

        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,CheckName(confSet->GetValue()))>=0) continue;

        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Vectorize(samples_str,',');
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadNormFactorOptions", "NormFactor: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "NormFactor: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
                return 1;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadNormFactorOptions", "NormFactor: " + CheckName(confSet->GetValue()) + " has regions set up for excluding that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "NormFactor: " + CheckName(confSet->GetValue()) + " has regions set up for excluding that do not exist");
                return 1;
            }
        }

        if (samples.size() > 0 && !CheckPresence(samples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadNormFactorOptions", "NormFactor: " + CheckName(confSet->GetValue()) + " has samples set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "NormFactor: " + CheckName(confSet->GetValue()) + " has samples set up that do not exist");
                return 1;
            }
        }

        nfactor = new NormFactor(CheckName(confSet->GetValue()));

        TRExFitter::SYSTMAP[nfactor->fName] = nfactor->fName;
        if( Common::FindInStringVector(fFitter->fNormFactorNames,nfactor->fName)<0 ){
            fFitter->fNormFactors.push_back( nfactor );
            fFitter->fNormFactorNames.push_back( nfactor->fName );
            fFitter->fNNorm++;
        }
        else{
            nfactor = fFitter->fNormFactors[ Common::FindInStringVector(fFitter->fNormFactorNames,nfactor->fName) ];
        }

        // Set NuisanceParameter = Name
        nfactor->fNuisanceParameter = nfactor->fName;
        TRExFitter::NPMAP[nfactor->fName] = nfactor->fName;

        if (SystHasProblematicName(nfactor->fNuisanceParameter)) return 1;

        // Set Constant
        param = confSet->Get("Constant");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") nfactor->fConst = true;
            else if(param=="FALSE") nfactor->fConst = false;
            else {
                WriteWarningStatus("ConfigReader::ReadNormFactorOptions", "You specified 'Constant' option but did not provide valid parameter. Using default (false)");
                nfactor->fConst = false;
            }
        }

        // Set Category
        param = confSet->Get("Category");
        if(param != "") nfactor->fCategory = RemoveQuotes(param);

        // Set SubCategory
        param = confSet->Get("SubCategory");
        if(param != "") nfactor->fSubCategory = RemoveQuotes(param);

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            nfactor->fTitle = RemoveQuotes(param);
            TRExFitter::SYSTMAP[nfactor->fName] = nfactor->fTitle;
        }

        // Set TexTitle
        param = confSet->Get("TexTitle");
        if(param != "") TRExFitter::SYSTTEX[nfactor->fName] = RemoveQuotes(param);

        // Set Min
        param = confSet->Get("Min");
        if(param!="") nfactor->fMin = atof(param.c_str());

        // Set Max
        param = confSet->Get("Max");
        if(param!="") nfactor->fMax = atof(param.c_str());

        // Set Nominal
        param = confSet->Get("Nominal");
        if(param!="") nfactor->fNominal = atof(param.c_str());

        // Set Expression
        param = confSet->Get("Expression");
        if(param!=""){
            std::vector<std::string> v = Vectorize(param,':');
            if (v.size() < 2){
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "You specified 'Expression' option but did not provide 2 parameters. Please check this");
                return 1;
            }
            nfactor->fExpression = std::make_pair(v[0],v[1]);
            // title will contain the expression FIXME
            nfactor->fTitle = v[0];
            TRExFitter::SYSTMAP[nfactor->fName] = v[0];
            // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
            nfactor->fNuisanceParameter = v[1];
            TRExFitter::NPMAP[nfactor->fName] = v[1];
            // set nominal, min and max according to the norm factor the expression depends on FIXME
            for(NormFactor *nf : fFitter->fNormFactors){
                if(nf->fNuisanceParameter == v[1]){
                    nfactor->fNominal = nf->fNominal;
                    nfactor->fMin = nf->fMin;
                    nfactor->fMax = nf->fMax;
                }
            }
        }

        // Set tau
        param = confSet->Get("Tau");
        if (param != "") {
            nfactor->fTau = std::stof(param);
        }

        // save list of
        if (regions.size() == 0 || exclude.size() == 0){
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "Region or exclude region size is equal to zero. Please check this");
                return 1;
        }
        if(regions[0] == "all") {
            nfactor->fRegions = GetAvailableRegions();
        } else {
            nfactor->fRegions = regions;
        }
        if(exclude[0] != "")    nfactor->fExclude = exclude;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sample = fFitter->fSamples[i_smp];
            if(sample->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || Common::FindInStringVector(samples, sample->fName)>=0 )
               && (exclude[0]==""    || Common::FindInStringVector(exclude, sample->fName)<0 ) ){
                sample->AddNormFactor(nfactor);
            }
        }
    }
    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadShapeFactorOptions(){
    std::string param = "";
    int nShape = 0;
    ShapeFactor *sfactor = nullptr;
    Sample *sample = nullptr;

    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("ShapeFactor",nShape);
        if (confSet == nullptr) break;
        nShape++;

        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,CheckName(confSet->GetValue()))>=0) continue;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Vectorize(samples_str,',');
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadShapeFactorOptions", "ShapeFactor: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadShapeFactorOptions", "ShapeFactor: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
                return 1;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadShapeFactorOptions", "ShapeFactor: " + CheckName(confSet->GetValue()) + " has regions set up for excluding that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadShapeFactorOptions", "ShapeFactor: " + CheckName(confSet->GetValue()) + " has regions set up for excluding that do not exist");
                return 1;
            }
        }

        if (samples.size() > 0 && !CheckPresence(samples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadShapeFactorOptions", "ShapeFactor: " + CheckName(confSet->GetValue()) + " has samples set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadShapeFactorOptions", "ShapeFactor: " + CheckName(confSet->GetValue()) + " has samples set up that do not exist");
                return 1;
            }
        }

        sfactor = new ShapeFactor(CheckName(confSet->GetValue()));
        if( Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
            fFitter->fShapeFactors.push_back( sfactor );
            fFitter->fShapeFactorNames.push_back( sfactor->fName );
            fFitter->fNShape++;
        }
        else{
            sfactor = fFitter->fShapeFactors[ Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName) ];
        }

        // Set NuisanceParameter = Name
        sfactor->fNuisanceParameter = sfactor->fName;
        TRExFitter::NPMAP[sfactor->fName] = sfactor->fName;

        if (SystHasProblematicName(sfactor->fNuisanceParameter)) return 1;

        // Set Constant
        param = confSet->Get("Constant");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") sfactor->fConst = true;
            else if(param=="FALSE") sfactor->fConst = false;
            else {
                WriteWarningStatus("ConfigReader::ReadShapeFactorOptions", "You specified 'Constant' option but did not provide valid parameter. Using default (false)");
                sfactor->fConst = false;
            }
        }

        // Set Category
        param = confSet->Get("Category");
        if(param!="") sfactor->fCategory = RemoveQuotes(param);

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            sfactor->fTitle = RemoveQuotes(param);
            TRExFitter::SYSTMAP[sfactor->fName] = sfactor->fTitle;
        }
        param = confSet->Get("TexTitle");
        if(param != "") TRExFitter::SYSTTEX[sfactor->fName] = RemoveQuotes(param);

        // Set Min
        param = confSet->Get("Min");
        if(param!="") sfactor->fMin = atof(param.c_str());

        // Set Max
        param = confSet->Get("Max");
        if(param!="") sfactor->fMax = atof(param.c_str());

        // Set Nominal
        param = confSet->Get("Nominal");
        if(param!="") sfactor->fNominal = atof(param.c_str());

        if (regions.size() == 0 || exclude.size() == 0){
            WriteErrorStatus("ConfigReader::ReadShapeFactorOptions", "Region or exclude region size is equal to zero. Please check this");
            return 1;
        }
        // save list of
        if(regions[0] == "all") {
            sfactor->fRegions = GetAvailableRegions();
        } else {
            sfactor->fRegions = regions;
        }
        if(exclude[0]!="")    sfactor->fExclude = exclude;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sample = fFitter->fSamples[i_smp];
            if(sample->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || Common::FindInStringVector(samples, sample->fName)>=0 )
               && (exclude[0]==""    || Common::FindInStringVector(exclude, sample->fName)<0 ) ){
                sample->AddShapeFactor(sfactor);
            }
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadSystOptions(){

    if (fOnlySystematics.size() > 0){
        std::vector<std::string> availableSysts = GetAvailableSysts();
        if (!CheckPresence(fOnlySystematics, availableSysts)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You set systematics that do not exist in your command line options");
            } else {
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "You set systematics that do not exist in your command line options");
                return 1;
            }
        }
    }

    int nSys = 0;

    Sample *sample = nullptr;
    //Addition for StatOnly fit: dummy systematic for the significance computation and limit setting
    if (fFitter->fStatOnly) {
        int typed = Systematic::OVERALL;
        Systematic *sysd = new Systematic("Dummy",typed);
        sysd->fOverallUp   = 0.;
        sysd->fOverallDown = -0.;
        sysd->fScaleUp   = 1.;
        sysd->fScaleDown   = 1.;
        fFitter->fSystematics.push_back( sysd );
        TRExFitter::SYSTMAP[sysd->fName] = "Dummy";
        fFitter->fNSyst++;
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sample = fFitter->fSamples[i_smp];
            if(sample->fType == Sample::SIGNAL ) {
                sample->AddSystematic(sysd);
            }
        }
    }

    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("Systematic",nSys);
        if (confSet == nullptr) break;
        nSys++;

        std::string param = "";
        Systematic *sys = nullptr;
        if(fOnlySystematics.size()>0 && Common::FindInStringVector(fOnlySystematics,CheckName(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,CheckName(confSet->GetValue()))>=0) continue;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        std::string excludeRegionSample_str = confSet->Get("ExcludeRegionSample");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Vectorize(samples_str,',');
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');

        if(samples_str!="all" && confSet->Get("DummyForSamples")!=""){
            std::vector<std::string> addSamples = Vectorize(confSet->Get("DummyForSamples"),',');
            for(const auto& addSmp : addSamples){
                if(Common::FindInStringVector(samples,addSmp)<0){
                    samples.emplace_back(addSmp);
                }
            }
        }

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up that do not exist");
                return 1;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples/regions set up for excluding that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples/regions set up for excluding that do not exist");
                return 1;
            }
        }

        if (samples.size() > 0 && !CheckPresence(samples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up that do not exist");
            } else {
                WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up that do not exist");
                return 1;
            }
        }

        fExcludeRegionSample = Vectorize(excludeRegionSample_str,',');
        int type = Systematic::HISTO;

        // Set type
        param = confSet->Get("Type");
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);

        if(param == "OVERALL") type = Systematic::OVERALL;
        else if(param == "SHAPE") type = Systematic::SHAPE;
        else if (param == "STAT") type = Systematic::STAT;

        std::string decorrelate = confSet->Get("Decorrelate");

        sys = new Systematic(CheckName(confSet->GetValue()),type);
        TRExFitter::SYSTMAP[sys->fName] = sys->fTitle;
        if(param == "OVERALL") sys->fIsNormOnly=true;

        sys->fSamples = samples;

        // SetCategory
        param = confSet->Get("Category");
        if(param != ""){
            sys->fCategory = RemoveQuotes(param);
            sys->fSubCategory = RemoveQuotes(param); //SubCategory defaults to the Category setting, if the Category is explicitly set
        }

        // CombineName
        param = confSet->Get("CombineName");
        if(param != ""){
            sys->fCombineName = RemoveQuotes(param);
        }

        // CombineType
        param = confSet->Get("CombineType");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            const std::string tmp = RemoveQuotes(param);

            if (tmp == "STANDARDDEVIATION") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::STANDARDDEVIATION;
            } else if (tmp == "ENVELOPE") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::ENVELOPE;
            } else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'CombineType' option but did not provide valid parameter. Using default (ENVELOPE)");
                sys->fCombineType = Systematic::COMBINATIONTYPE::ENVELOPE;
            }
        }

        // SetSubCategory
        param = confSet->Get("SubCategory");
        if (param != ""){
            sys->fSubCategory = RemoveQuotes(param); // note this needs to happen after Category was set, in order to overwrite the default if required
        }

        // Set IsFreeParameter
        // Experimental
        param = confSet->Get("IsFreeParameter");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "TRUE") sys->fIsFreeParameter = true;
            else if (param == "FALSE") sys->fIsFreeParameter = false;
            else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'IsFreeParameter' option but did not provide valid parameter. Using default (false)");
                sys->fIsFreeParameter = false;
            }
        }

        // Set StoredName
        // New: name to use when writing / reading the Histograms file
        param = confSet->Get("StoredName");
        if(param != "") sys->fStoredName = RemoveQuotes(param);

        if(type==Systematic::HISTO || type==Systematic::SHAPE){
            bool hasUp(false);
            bool hasDown(false);

            if(fFitter->fInputType==0){ // HIST input
                param = confSet->Get("HistoPathUp");
                if(param!=""){
                    sys->fHistoPathsUp.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathDown");
                if(param!=""){
                    sys->fHistoPathsDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoPathSufUp");
                if(param!=""){
                    sys->fHistoPathSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathSufDown");
                if(param!=""){
                    sys->fHistoPathSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoFileUp");
                if(param!=""){
                    sys->fHistoFilesUp.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoFileDown");
                if(param!=""){
                    sys->fHistoFilesDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoFilesUp");
                if(param!=""){
                    sys->fHistoFilesUp = Vectorize(param, ',');
                    hasUp   = true;
                }
                param = confSet->Get("HistoFilesDown");
                if(param!=""){
                    sys->fHistoFilesDown = Vectorize(param, ',');
                    hasDown = true;
                }
                param = confSet->Get("HistoFileSufUp");
                if(param!=""){
                    sys->fHistoFileSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoFileSufDown");
                if(param!=""){
                    sys->fHistoFileSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoNameUp");
                if(param!=""){
                    sys->fHistoNamesUp.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameDown");
                if(param!=""){
                    sys->fHistoNamesDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoNameSufUp");
                if(param!=""){
                    sys->fHistoNameSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameSufDown");
                if(param!=""){
                    sys->fHistoNameSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                // For reference file when using systematics on it - like JER on data
                param = confSet->Get("HistoPathUpRefSample");
                if(param!=""){
                    sys->fHistoPathsUpRefSample.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathDownRefSample");
                if(param!=""){
                    sys->fHistoPathsDownRefSample.push_back(RemoveQuotes(param));
                    hasDown   = true;
                }
                param = confSet->Get("HistoPathsUpRefSample");
                if(param!=""){
                    sys->fHistoPathsUpRefSample = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathsDownRefSample");
                if(param!=""){
                    sys->fHistoPathsDownRefSample = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("HistoPathSufUpRefSample");
                if(param!=""){
                    sys->fHistoPathSufUpRefSample = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathSufDownRefSample");
                if(param!=""){
                    sys->fHistoPathSufDownRefSample = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoFileUpRefSample");
                if(param!=""){
                    sys->fHistoFilesUpRefSample.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoFileDownRefSample");
                if(param!=""){
                    sys->fHistoFilesDownRefSample.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoFilesUpRefSample");
                if(param!=""){
                    sys->fHistoFilesUpRefSample = Vectorize(param, ',');
                    hasUp   = true;
                }
                param = confSet->Get("HistoFilesDownRefSample");
                if(param!=""){
                    sys->fHistoFilesDownRefSample = Vectorize(param, ',');
                    hasDown = true;
                }
                param = confSet->Get("HistoFilesUpRefSample");
                if(param!=""){
                    sys->fHistoFilesUpRefSample = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("HistoFilesDownRefSample");
                if(param!=""){
                    sys->fHistoFilesDownRefSample = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("HistoFileSufUpRefSample");
                if(param!=""){
                    sys->fHistoFileSufUpRefSample = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoFileSufDownRefSample");
                if(param!=""){
                    sys->fHistoFileSufDownRefSample = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoNameUpRefSample");
                if(param!=""){
                    sys->fHistoNamesUpRefSample.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameDownRefSample");
                if(param!=""){
                    sys->fHistoNamesDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoNamesUpRefSample");
                if(param!=""){
                    sys->fHistoNamesUpRefSample = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("HistoNamesDownRefSample");
                if(param!=""){
                    sys->fHistoNamesDownRefSample = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("HistoNameSufUpRefSample");
                if(param!=""){
                    sys->fHistoNameSufUpRefSample = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameSufDownRefSample");
                if(param!=""){
                    sys->fHistoNameSufDownRefSample = RemoveQuotes(param);
                    hasDown = true;
                }
            }
            else if(fFitter->fInputType==1){ // NTUP option
                param = confSet->Get("NtuplePathUp");
                if(param!=""){
                    sys->fNtuplePathsUp.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathDown");
                if(param!=""){
                    sys->fNtuplePathsDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathsUp");
                if(param!=""){
                    sys->fNtuplePathsUp = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathsDown");
                if(param!=""){
                    sys->fNtuplePathsDown = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathSufUp");
                if(param!=""){
                    sys->fNtuplePathSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathSufDown");
                if(param!=""){
                    sys->fNtuplePathSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileUp");
                if(param!=""){
                    sys->fNtupleFilesUp.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFileDown");
                if(param!=""){
                    sys->fNtupleFilesDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("NtupleFilesUp");
                if(param!=""){
                    sys->fNtupleFilesUp = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFilesDown");
                if(param!=""){
                    sys->fNtupleFilesDown = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileSufUp");
                if(param!=""){
                    sys->fNtupleFileSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFileSufDown");
                if(param!=""){
                    sys->fNtupleFileSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameUp");
                if(param!=""){
                    sys->fNtupleNamesUp.push_back(RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNameDown");
                if(param!=""){
                    sys->fNtupleNamesDown.push_back(RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("NtupleNamesUp");
                if(param!=""){
                    sys->fNtupleNamesUp = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNamesDown");
                if(param!=""){
                    sys->fNtupleNamesDown = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameSufUp");
                if(param!=""){
                    sys->fNtupleNameSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNameSufDown");
                if(param!=""){
                    sys->fNtupleNameSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("WeightUp");
                if(param!=""){
                    sys->fWeightUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("WeightDown");
                if(param!=""){
                    sys->fWeightDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("WeightSufUp");
                if(param!=""){
                    sys->fWeightSufUp = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("WeightSufDown");
                if(param!=""){
                    sys->fWeightSufDown = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("IgnoreWeight");
                if(param!=""){
                    sys->fIgnoreWeight = RemoveQuotes(param);
                }
                // For reference file when using systematics on it - like JER on data
                param = confSet->Get("NtuplePathUpRefSample");
                if(param!=""){
                    sys->fNtuplePathsUpRefSample.push_back( RemoveQuotes(param) );
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathDownRefSample");
                if(param!=""){
                    sys->fNtuplePathsDownRefSample.push_back( RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathsUpRefSample");
                if(param!=""){
                    sys->fNtuplePathsUpRefSample = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathsDownRefSample");
                if(param!=""){
                    sys->fNtuplePathsDownRefSample = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathSufUpRefSample");
                if(param!=""){
                    sys->fNtuplePathSufUpRefSample = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathSufDownRefSample");
                if(param!=""){
                    sys->fNtuplePathSufDownRefSample = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileUpRefSample");
                if(param!=""){
                    sys->fNtupleFilesUpRefSample.push_back( RemoveQuotes(param) );
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFileDownRefSample");
                if(param!=""){
                    sys->fNtupleFilesDownRefSample.push_back( RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("NtupleFilesUpRefSample");
                if(param!=""){
                    sys->fNtupleFilesUpRefSample = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFilesDownRefSample");
                if(param!=""){
                    sys->fNtupleFilesDownRefSample = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileSufUpRefSample");
                if(param!=""){
                    sys->fNtupleFileSufUpRefSample = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFileSufDownRefSample");
                if(param!=""){
                    sys->fNtupleFileSufDownRefSample = RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameUpRefSample");
                if(param!=""){
                    sys->fNtupleNamesUpRefSample.push_back( RemoveQuotes(param) );
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNameDownRefSample");
                if(param!=""){
                    sys->fNtupleNamesDownRefSample.push_back( RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("NtupleNamesUpRefSample");
                if(param!=""){
                    sys->fNtupleNamesUpRefSample = Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNamesDownRefSample");
                if(param!=""){
                    sys->fNtupleNamesDownRefSample = Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameSufUpRefSample");
                if(param!=""){
                    sys->fNtupleNameSufUpRefSample = RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNameSufDownRefSample");
                if(param!=""){
                    sys->fNtupleNameSufDownRefSample = RemoveQuotes(param);
                    hasDown = true;
                }
            }
            sys->fHasUpVariation   = hasUp  ;
            sys->fHasDownVariation = hasDown;

            // Set Symmetrisation
            param = confSet->Get("Symmetrisation");
            if(param != ""){
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if(param == "ONESIDED"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEONESIDED;
                }
                else if(param == "TWOSIDED"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZETWOSIDED;
                }
                else if(param == "ABSMEAN"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEABSMEAN;
                }
                else if(param == "MAXIMUM"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEMAXIMUM;
                }
                else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Symetrisation scheme is not recognized ... ");
                    return 1;
                }
            }

            // Set Smoothing
            param = confSet->Get("Smoothing");
            if(param != "") sys->fSmoothType = atoi(param.c_str());

            param = confSet->Get("PreSmoothing");
            if(param != ""){
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if(param == "TRUE") sys->fPreSmoothing = true;
                else if(param == "FALSE") sys->fPreSmoothing = false;
                else {
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'PreSmoothing' option but did not provide valid parameter. Using default (false)");
                    sys->fPreSmoothing = false;
                }
            }

            param = confSet->Get("SmoothingOption");
            if(param != ""){
                sys->fSampleSmoothing = true;
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if( param == "MAXVARIATION" ) sys->fSampleSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
                else if (param == "TTBARRESONANCE") sys->fSampleSmoothOption = HistoTools::SmoothOption::TTBARRESONANCE;
                else if (param == "COMMONTOOLSMOOTHMONOTONIC") sys->fSampleSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHMONOTONIC;
                else if (param == "COMMONTOOLSMOOTHPARABOLIC") sys->fSampleSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHPARABOLIC;
                else if (param == "KERNELRATIOUNIFORM") sys->fSampleSmoothOption = HistoTools::SmoothOption::KERNELRATIOUNIFORM;
                else if (param == "KERNELDELTAGAUSS") sys->fSampleSmoothOption = HistoTools::SmoothOption::KERNELDELTAGAUSS;
                else if (param == "KERNELRATIOGAUSS") sys->fSampleSmoothOption = HistoTools::SmoothOption::KERNELRATIOGAUSS;
                else {
                    WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'SmoothingOption' option but you didn't provide valid input. Using default from job block");
                    sys->fSampleSmoothing = false;
                }
            }
        } // end of if(type==Systematic::HISTO || type==Systematic::SHAPE)
        else if(type==Systematic::OVERALL){
            // Set OverallUp
            param = confSet->Get("OverallUp");
            if (param != "") sys->fOverallUp = atof( param.c_str() );

            // Set OverallUDown
            param = confSet->Get("OverallDown");
            if (param != "") sys->fOverallDown = atof( param.c_str() );
        }

        // Set ScaleUp
        param = confSet->Get("ScaleUp");
        if(param!=""){
            std::vector < std::string > temp_vec = Vectorize(param,',',false);
            if(temp_vec.size()==1 && Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleUp = atof(param.c_str());
            }
            else{
                for(std::string ireg : temp_vec){
                    std::vector < std::string > reg_value = Vectorize(ireg,':');
                    if(reg_value.size()==2){
                        sys->fScaleUpRegions.insert( std::pair < std::string, double >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }

        // Set ScaleDown
        param = confSet->Get("ScaleDown");
        if(param!=""){
            std::vector < std::string > temp_vec = Vectorize(param,',',false);
            if(temp_vec.size()==1 && Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleDown = atof(param.c_str());
            }
            else{
                for(std::string ireg : temp_vec){
                    std::vector < std::string > reg_value = Vectorize(ireg,':');
                    if(reg_value.size()==2){
                        sys->fScaleDownRegions.insert( std::pair < std::string, double >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }

        // Set SampleUp
        // a new way of defining systematics, just comparing directly with GHOST samples
        // --> this can be used only if this systematic is applied to a single sample
        param = confSet->Get("SampleUp");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in SampleUp that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in SampleUp that do not exist");
                    return 1;
                }
            }
            sys->fSampleUp = RemoveQuotes(param);
        }

        // Set SampleDown
        param = confSet->Get("SampleDown");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in SampleDown that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in SampleDown that do not exist");
                    return 1;
                }
            }
            sys->fSampleDown = RemoveQuotes(param);
        }

        // Set ReferenceSample
        // this to obtain syst variation relatively to given sample
        param = confSet->Get("ReferenceSample");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fAvailableSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in ReferenceSample that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in ReferenceSample that do not exist");
                    return 1;
                }
            }
            sys->fReferenceSample = RemoveQuotes(param);
        }

        // Set KeepReferenceOverallVar
        param = confSet->Get("KeepReferenceOverallVar");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") sys->fKeepReferenceOverallVar = false;
            else if(param == "TRUE" ) sys->fKeepReferenceOverallVar = true;
            else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'KeepReferenceOverallVar' option but did not provide valid parameter. Using default (true)");
                sys->fKeepReferenceOverallVar = true;
            }
        }

        // Set ReferenceSmoothing
        // this is to obtain syst variation relatively to given sample
        param = confSet->Get("ReferenceSmoothing");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fAvailableSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in ReferenceSmoothing that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in ReferenceSmoothing that do not exist");
                    return 1;
                }
            }
            if (std::find(samples.begin(), samples.end(), RemoveQuotes(param)) == samples.end()){
                WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " requires that the ReferenceSample appears in Samples for this systematic");
                return 1;
            }
            sys->fReferenceSmoothing = RemoveQuotes(param);
        }

        // Set ReferencePruning
        // this allows to prune wrt to a sample and then apply the result to all samples
        param = confSet->Get("ReferencePruning");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), RemoveQuotes(param)) == fAvailableSamples.end()){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in ReferencePruning that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in ReferencePruning that do not exist");
                    return 1;
                }
            }
            if (std::find(samples.begin(), samples.end(), RemoveQuotes(param)) == samples.end()){
                WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " requires that the ReferencePruning appears in Samples for this systematic");
                return 1;
            }
            sys->fReferencePruning = RemoveQuotes(param);
        }

        // Set DropShapeIn
        param = confSet->Get("DropShapeIn");
        if(param!="") {
            std::vector<std::string> tmp = Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions)){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up in DropShapeIn that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up in DropShapeIn that do not exist");
                    return 1;
                }
            }
            sys->fDropShapeIn = tmp;
        }

        // Set DropNorm
        param = confSet->Get("NoPruning");
        if (param!="") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "TRUE") {
                sys->fNoPruning = true;
            } else if (param == "FALSE") {
                sys->fNoPruning = false;
            } else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'NoPruning' option but did not provide valid parameter. Using default (FALSE)");
                sys->fNoPruning = false;
            }
        }

        // Set DropNorm
        param = confSet->Get("DropNorm");
        if(param!=""){
            std::vector<std::string> tmp = Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up in DropNorm that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up in DropNorm that do not exist");
                    return 1;
                }
            }
            sys->fDropNormIn = tmp;
        }

        // Set DropNormSpecial
        param = confSet->Get("DropNormSpecial");
        if(param!=""){
            std::vector<std::string> tmp = Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up in DropNorm that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has regions set up in DropNorm that do not exist");
                    return 1;
                }
            }
            sys->fDropNormSpecialIn = tmp;
        }

        // Set KeepNormForSamples
        param = confSet->Get("KeepNormForSamples");
        if(param!="") {
            std::vector<std::string> tmp = Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in KeepNormForSamples that do not exist");
                } else {
                    WriteErrorStatus("ConfigReader::ReadSystOptions", "Systematic: " + CheckName(confSet->GetValue()) + " has samples set up in KeepNormForSamples that do not exist");
                    return 1;
                }
            }
            sys->fKeepNormForSamples = tmp;
        }

        if (regions.size() == 0 || exclude.size() == 0){
            WriteErrorStatus("ConfigReader::ReadSystOptions", "Region or exclude region size is equal to zero. Please check this");
            return 1;
        }

        // Set DummyForSamples
        param = confSet->Get("DummyForSamples");
        if(param!="") {
            std::vector<std::string> tmp = Vectorize(param,',');
            sys->fDummyForSamples = tmp;
        }

        // Set ForceShape
        param = confSet->Get("ForceShape");
        if(param!="") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "NOSHAPE") {
                sys->fForceShape = HistoTools::FORCESHAPETYPE::NOSHAPE;
            } else if (param == "LINEAR") {
                sys->fForceShape = HistoTools::FORCESHAPETYPE::LINEAR;
            } else if (param == "TRIANGULAR") {
                sys->fForceShape = HistoTools::FORCESHAPETYPE::TRIANGULAR;
            } else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'ForceShape' option but did not provide valid parameter. Using default (NOSHAPE)");
                sys->fForceShape = HistoTools::FORCESHAPETYPE::NOSHAPE;
            }
        }

        // Set SubtractRefSampleVar
        // New: for systematics which also vary Data (e.g. JER with Full NPs)
        // This will subtract linearly the relative variation on Data from each relative variation on MC
        param = confSet->Get("SubtractRefSampleVar");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE" ) sys->fSubtractRefSampleVar = true;// default is false
            else if(param == "FALSE" ) sys->fSubtractRefSampleVar = false;// default is false
            else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'PreSmoothing' option but did not provide valid parameter. Using default (false)");
                sys->fSubtractRefSampleVar = false;
            }
        }

        if(regions[0]!="all") sys->fRegions = regions;
        if(exclude[0]!="")    sys->fExclude = exclude;

        for (std::string ier : fExcludeRegionSample){
            std::vector<std::string> pair_ERS = Vectorize(ier,':');
            if (pair_ERS.size()==2){
                sys->fExcludeRegionSample.push_back(pair_ERS);
            }
        }
        if ( decorrelate == "" && type != Systematic::STAT) {
            if (SetSystNoDecorelate(confSet, sys, samples, exclude) != 0) return 1;
        }
        else if (decorrelate == "REGION" || type == Systematic::STAT)  {
            if (SetSystRegionDecorelate(confSet, sys, samples, exclude, regions, type) != 0) return 1;
        }
        else if (decorrelate == "SAMPLE")  {
            if (SetSystSampleDecorelate(confSet, sys, samples, exclude) != 0) return 1;
        }
        else if (decorrelate == "SHAPEACC")  {
            if (SetSystShapeDecorelate(confSet, sys, samples, exclude) != 0) return 1;
        }
        else {
            WriteErrorStatus("ConfigReader::ReadSystOptions", "decorrelate option: " + decorrelate  + "  not supported ...");
            WriteErrorStatus("ConfigReader::ReadSystOptions", "       PLEASE USE ONLY: REGION, SAMPLE, SHAPEACC");
            return 1;
        }


    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystNoDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude){
    Sample *sam = nullptr;

    fFitter->fSystematics.push_back( sys );
    fFitter->fNSyst++;

    // Set NuisanceParameter
    std::string param = confSet->Get("NuisanceParameter");
    if(param!=""){
        sys->fNuisanceParameter = RemoveQuotes(param);
        TRExFitter::NPMAP[sys->fName] = sys->fNuisanceParameter;
    }
    else{
        sys->fNuisanceParameter = sys->fName;
        TRExFitter::NPMAP[sys->fName] = sys->fName;
    }

    if (SystHasProblematicName(sys->fNuisanceParameter)) return 1;

    // Set Title
    param = confSet->Get("Title");
    if(param != ""){
        sys->fTitle = RemoveQuotes(param);
        TRExFitter::SYSTMAP[sys->fName] = sys->fTitle;
    }

    // Set TexTitle
    param = confSet->Get("TexTitle");
    if(param!="") TRExFitter::SYSTTEX[sys->fName] = RemoveQuotes(param);

    // attach the syst to the proper samples
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        sam = fFitter->fSamples[i_smp];
        // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
        if(sam->fType == Sample::DATA){
            if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
              sam->AddSystematic(sys);
            }
            else continue;
        }
        if(!sam->fUseSystematics) continue;
        if((samples[0]=="all" || Common::FindInStringVector(samples, sam->fName)>=0 )
           && (exclude[0]=="" || Common::FindInStringVector(exclude, sam->fName)<0 ) ){
            sam->AddSystematic(sys);
        }
    }
    if(Common::FindInStringVector(fFitter->fDecorrSysts,sys->fNuisanceParameter)>=0){
        WriteInfoStatus("ConfigReader::ReadSystOptions","Decorrelating systematic with NP = " + sys->fNuisanceParameter);
        sys->fNuisanceParameter += fFitter->fDecorrSuff;
    }

    FixReferenceSamples(sys);

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystRegionDecorelate(ConfigSet *confSet,
                                          Systematic *sys,
                                          const std::vector<std::string>& samples,
                                          const std::vector<std::string>& exclude,
                                          const std::vector<std::string>& regions,
                                          int type) {
    Sample *sam = nullptr;
    std::string param = "";

    for(const auto& ireg : fRegNames) {
        bool keepReg=false;
        if ( regions[0]=="all" ) keepReg=true;
        else {
            for ( const std::string& iGoodReg : regions) {
                if ( ireg == iGoodReg ) keepReg=true;
            }
        }
        for ( const std::string& iBadReg : exclude) {
            if ( iBadReg == ireg ) keepReg=false;
        }

        if (!keepReg) {
            WriteInfoStatus("ConfigReader::SetSystRegionDecorelate", "IGNORING REGION: " + ireg);
            continue;
        }
        WriteInfoStatus("ConfigReader::SetSystRegionDecorelate", "--> KEEPING IT!!! " + ireg);

        Region* reg = fFitter->GetRegion(ireg);

        if (type == Systematic::STAT) {
            unsigned int nbins = reg->fHistoNBinsRebin>0 ? reg->fHistoNBinsRebin : reg->fNbins;
            WriteInfoStatus("ConfigReader::SetSystRegionDecorelate", ireg + " " + std::to_string(nbins));
            // decorrelate by bin
            for (unsigned int i_bin = 0; i_bin < nbins; i_bin++) {
                Systematic* mySys= new Systematic(*sys);
                mySys->fName=(mySys->fName)+"_"+ireg+"_bin"+std::to_string(i_bin);
                std::vector<std::string> tmpReg;
                tmpReg.push_back( ireg );
                mySys->fRegions = tmpReg;
                std::vector<int> tmpBin;
                tmpBin.push_back( i_bin );
                mySys->fBins = tmpBin;
                fFitter->fSystematics.push_back( mySys );

                // Set NuisanceParameter
                param = confSet->Get("NuisanceParameter");
                if(param != ""){
                    mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+ireg+"_bin"+std::to_string(i_bin);
                    TRExFitter::NPMAP[mySys->fName] = mySys->fNuisanceParameter;
                } else {
                    mySys->fNuisanceParameter = mySys->fName;
                    TRExFitter::NPMAP[mySys->fName] = mySys->fName;
                }

                if (SystHasProblematicName(mySys->fNuisanceParameter)) return 1;

                // Set Title
                param = confSet->Get("Title");
                if(param != ""){
                    mySys->fTitle = (sys->fTitle)+"_"+ireg+"_bin"+std::to_string(i_bin);
                    TRExFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
                }
                fFitter->fNSyst++;
                for (int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
                    sam = fFitter->fSamples[i_smp];
                    if(sam->fType == Sample::DATA) continue;
                    if(!sam->fUseSystematics) continue;
                    if(   (samples[0]=="all" || Common::FindInStringVector(samples, sam->fName)>=0 )
                       && (exclude[0]==""    || Common::FindInStringVector(exclude, sam->fName)<0 ) ){
                        sam->AddSystematic(mySys);
                    }
                }

                FixReferenceSamples(mySys);
            }
        } else {
            //
            // cloning the sys for each region
            Systematic* mySys= new Systematic(*sys);
            mySys->fName=(mySys->fName)+"_"+ireg;
            std::vector<std::string> tmpReg;
            tmpReg.push_back( ireg );
            mySys->fRegions = tmpReg;
            fFitter->fSystematics.push_back( mySys );

            // Set NuisanceParameter
            param = confSet->Get("NuisanceParameter");
            if(param != ""){
                mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+ireg;
                TRExFitter::NPMAP[mySys->fName] = mySys->fNuisanceParameter;
            }
            else{
                mySys->fNuisanceParameter = mySys->fName;
                TRExFitter::NPMAP[mySys->fName] = mySys->fName;
            }

            if (SystHasProblematicName(mySys->fNuisanceParameter)) return 1;

            // Set Title
            param = confSet->Get("Title");
            if(param != ""){
                mySys->fTitle = (sys->fTitle)+" ("+reg->fLabel+")";
                TRExFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
            }
            fFitter->fNSyst++;
            //
            for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
                sam = fFitter->fSamples[i_smp];
                // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
                if(sam->fType == Sample::DATA){
                    if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
                        sam->AddSystematic(mySys);
                    }
                    else continue;
                }
                if(!sam->fUseSystematics) continue;
                if(   (samples[0]=="all" || Common::FindInStringVector(samples, sam->fName)>=0 )
                    && (exclude[0]==""    || Common::FindInStringVector(exclude, sam->fName)<0 ) ){
                    sam->AddSystematic(mySys);
                }
            }
            FixReferenceSamples(mySys);
        }
    }

    delete sys;

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystSampleDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string> &samples, const std::vector<std::string> &exclude){
    Sample *sam = nullptr;
    std::string param = "";

    // (this is really messy)
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        sam = fFitter->fSamples[i_smp];
        // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
        if(sam->fType == Sample::DATA){
          if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
            sam->AddSystematic(sys);
          }
          else continue;
        }
        if(sam->fType == Sample::GHOST) continue;
        bool keepSam=false;
        if ( samples[0]=="all" ) keepSam=true;
        else {
            if ( find(samples.begin(), samples.end(), sam->fName)!=samples.end() ) keepSam=true;
        }
        if ( find(exclude.begin(), exclude.end(), sam->fName)!=exclude.end() ) keepSam=false;
        if (!keepSam) {
            WriteInfoStatus("ConfigReader::SetSystSampleDecorelate", " IGNORING SAMPLE: " + sam->fName);
            continue;
        }
        WriteInfoStatus("ConfigReader::SetSystSampleDecorelate", " --> KEEPING SAMPLE: " + sam->fName);
        //
        // cloning the sys for each sample
        Systematic* mySys= new Systematic(*sys);
        mySys->fName=(mySys->fName)+"_"+sam->fName;
        fFitter->fSystematics.push_back( mySys );

        // Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+sam->fName;
            TRExFitter::NPMAP[mySys->fName] = mySys->fNuisanceParameter;
        }
        else{
            mySys->fNuisanceParameter = mySys->fName;
            TRExFitter::NPMAP[mySys->fName] = mySys->fName;
        }

        if (SystHasProblematicName(mySys->fNuisanceParameter)) return 1;

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            mySys->fTitle = (sys->fTitle)+" "+sam->fTitle;
            TRExFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
        }
        fFitter->fNSyst++;

        // Add sample/syst cross-reference
        sam->AddSystematic(mySys);
        mySys->fSamples = { sam->fName };

        FixReferenceSamples(mySys);
    }

    delete sys;

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystShapeDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string> &samples, const std::vector<std::string> &exclude){
    Sample *sam = nullptr;
    std::string param = "";

    // cloning the sys
    Systematic* mySys1= new Systematic(*sys);
    mySys1->fName=(mySys1->fName)+"_Acc";
    fFitter->fSystematics.push_back( mySys1 );
    mySys1->fIsNormOnly = true;
    mySys1->fIsShapeOnly = false;

    // Set NuisanceParameter
    param = confSet->Get("NuisanceParameter");
    if(param != ""){
        mySys1->fNuisanceParameter = (sys->fNuisanceParameter)+"_Acc";
        TRExFitter::NPMAP[mySys1->fName] = mySys1->fNuisanceParameter;
    }
    else{
        mySys1->fNuisanceParameter = mySys1->fName;
        TRExFitter::NPMAP[mySys1->fName] = mySys1->fName;
    }

    if (SystHasProblematicName(mySys1->fNuisanceParameter)) return 1;

    // Set Title
    param = confSet->Get("Title");
    if(param != ""){
        mySys1->fTitle = (sys->fTitle)+" Acc";
        TRExFitter::SYSTMAP[mySys1->fName] = mySys1->fTitle;
    }
    fFitter->fNSyst++;

    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        sam = fFitter->fSamples[i_smp];
        // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
        if(sam->fType == Sample::DATA){
            if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
                sam->AddSystematic(mySys1);
            }
            else continue;
        }
        if(!sam->fUseSystematics) continue;
        if(   (samples[0]=="all" || Common::FindInStringVector(samples, sam->fName)>=0 )
           && (exclude[0]==""    || Common::FindInStringVector(exclude, sam->fName)<0 ) ){
            sam->AddSystematic(mySys1);
        }
    }
    FixReferenceSamples(mySys1);

    if ( sys->fType!=Systematic::OVERALL ) {
        // cloning the sys
        Systematic* mySys2= new Systematic(*sys);
        mySys2->fName=(mySys2->fName)+"_Shape";
        mySys2->fIsNormOnly=false;
        mySys2->fIsShapeOnly=true;

        fFitter->fSystematics.push_back( mySys2 );

        //Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            mySys2->fNuisanceParameter = (sys->fNuisanceParameter)+"_Shape";
            TRExFitter::NPMAP[mySys2->fName] = mySys2->fNuisanceParameter;
        }
        else{
            mySys2->fNuisanceParameter = mySys2->fName;
            TRExFitter::NPMAP[mySys2->fName] = mySys2->fName;
        }

        if (SystHasProblematicName(mySys2->fNuisanceParameter)) return 1;

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            mySys2->fTitle = (sys->fTitle)+" Shape";
            TRExFitter::SYSTMAP[mySys2->fName] = mySys2->fTitle;
        }
        fFitter->fNSyst++;

        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sam = fFitter->fSamples[i_smp];
            if(sam->fType == Sample::DATA) continue;
            if(!sam->fUseSystematics) continue;
            if(   (samples[0]=="all" || Common::FindInStringVector(samples, sam->fName)>=0 )
               && (exclude[0]==""    || Common::FindInStringVector(exclude, sam->fName)<0 ) ){
                sam->AddSystematic(mySys2);
            }
        }
        FixReferenceSamples(mySys2);
    }

    delete sys;

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::PostConfig(const std::string& opt){

    if (!fHasAtLeastOneValidSample && Common::OptionRunsFit(opt)){
        WriteErrorStatus("ConfigReader::ReadSampleOptions","You need to provide at least one sample that is either SIGNAL or BACKGROUND, otherwise the fit will crash.");
        return 1;
    }

    // if StatOnly, also sets to OFF the MC stat
    if(fFitter->fStatOnly){
        WriteInfoStatus("ConfigReader::PostConfig","StatOnly option is setting to OFF the MC-stat (gammas) as well.");
        WriteInfoStatus("ConfigReader::PostConfig","To keep them on use the command line option 'Systematics=NONE'");
        WriteInfoStatus("ConfigReader::PostConfig","or comment out all Systematics in config file.");
        fFitter->SetStatErrorConfig( false, 0. );
    }

    // add nuisance parameter - systematic title correspondence
    for(auto syst : fFitter->fSystematics){
        if(syst->fNuisanceParameter!=syst->fName) TRExFitter::SYSTMAP[syst->fNuisanceParameter] = syst->fTitle;
    }
    // add nuisance parameter - norm-factor title correspondence & fix nuisance parameter
    for(auto norm : fFitter->fNormFactors){
        if(TRExFitter::NPMAP[norm->fName]=="") TRExFitter::NPMAP[norm->fName] = norm->fName;
        if(norm->fNuisanceParameter!=norm->fName) TRExFitter::SYSTMAP[norm->fNuisanceParameter] = norm->fTitle;
    }

    // morphing
    if (fFitter->fMorphParams.size()!=0){
        // template fitting stuff
        fFitter->fTemplateWeightVec = fFitter->GetTemplateWeightVec(fFitter->fTemplateInterpolationOption);
        for(const TRExFit::TemplateWeight& itemp : fFitter->fTemplateWeightVec){
            std::string normName = "morph_"+itemp.name+"_"+Common::ReplaceString(std::to_string(itemp.value),"-","m");
            TRExFitter::SYSTMAP[normName] = itemp.function;
            TRExFitter::NPMAP[normName]   = itemp.name;
            // get the norm factor corresponding to each template
            for(auto norm : fFitter->fNormFactors){
                if(norm->fName == normName){
                    // find a norm factor in the config corresponding to the morphing parameter
                    // NB: it should be there in the config, otherwise an error message is shown and the code crashe
                    bool found = false;
                    for(auto norm2 : fFitter->fNormFactors){
                        if(norm2->fName == itemp.name){
                            norm->fNominal = norm2->fNominal;
                            norm->fMin = norm2->fMin;
                            norm->fMax = norm2->fMax;
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        WriteErrorStatus("ConfigReader::PostConfig", "No NormFactor with name " + itemp.name + " found (needed for morphing");
                        WriteErrorStatus("ConfigReader::PostConfig", "Please add to the config something like:");
                        WriteErrorStatus("ConfigReader::PostConfig", "  NormFactor: " + itemp.name);
                        WriteErrorStatus("ConfigReader::PostConfig", "    Min: <the min value for which you have provided template>");
                        WriteErrorStatus("ConfigReader::PostConfig", "    Max: <the min value for which you have provided template>");
                        WriteErrorStatus("ConfigReader::PostConfig", "    Samples: none");
                        return 1;
                    }
                }
            }
        }
    }

    // Check and fix shape factor number of bins
    for(auto sfactor : fFitter->fShapeFactors){
        WriteInfoStatus("ConfigReader::PostConfig","Checking consistency of shape factor " + sfactor->fName + " across regions...");
        int nbins = 0;
        for(const auto& regName : sfactor->fRegions){
            Region *reg = fFitter->GetRegion(regName);
            if(reg==nullptr) continue; // skip non-existing regions
            if(nbins==0) nbins = reg->fNbins;
            else{
                if(nbins!=reg->fNbins){
                    WriteErrorStatus("ConfigReader::PostConfig","Shape factor " + sfactor->fName + " assigned to regions with different number of bins. Please check.");
                    return 1;
                }
            }
        }
        sfactor->fNbins = nbins;
    }

    // if SaturatedModel was set, add proper shape factors
    if(fFitter->fSaturatedModel){
        for(auto reg : fFitter->fRegions){
            std::string sfactorName = "saturated_model_sf_" + reg->fName;
            ShapeFactor *sfactor = new ShapeFactor(sfactorName);
            if( Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
                fFitter->fShapeFactors.push_back( sfactor );
                fFitter->fShapeFactorNames.push_back( sfactor->fName );
                fFitter->fNShape++;
            }
            else{
                sfactor = fFitter->fShapeFactors[ Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName) ];
            }
            // Set NuisanceParameter = Name, Title = Name
            sfactor->fNuisanceParameter = sfactor->fName;
            TRExFitter::NPMAP[sfactor->fName] = sfactor->fName;
            sfactor->fTitle = sfactor->fName;
            TRExFitter::SYSTMAP[sfactor->fName] = sfactor->fTitle;
            if (SystHasProblematicName(sfactor->fNuisanceParameter)) return 1;
            // Set it to constant by default (will be made non-constant later when fitting)
//             sfactor->fConst = true;
            // needed? FIXME
            sfactor->fMin = 0;
            sfactor->fMax = 10;
            sfactor->fNominal = 1;
            // set nbins
            sfactor->fNbins = reg->fNbins;
            // save list of
            sfactor->fRegions.push_back(reg->fName);
            // attach the syst to all non-data samples
            for(auto sample : fFitter->fSamples){
                if(sample->fType == Sample::DATA) continue;
                sample->AddShapeFactor(sfactor);
            }
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadUnfoldingOptions() {
    ConfigSet *confSet = fParser->GetConfigSet("Unfolding");

    if (fFitter->fFitType == TRExFit::UNFOLDING && !confSet) {
        WriteErrorStatus("ConfigReader::ReadUnfoldingOptions", "You set FitType == UNFOLDING, but didnt provide Unfolding block!");
        return 1;
    }

    if (!confSet) {
        return 0;
    }

    std::string param = confSet->Get("MatrixOrientation");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUTHONHORIZONTAL") {
            fFitter->fMatrixOrientation = FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS;
        } else if (param == "TRUTHONVERTICAL") {
            fFitter->fMatrixOrientation = FoldingManager::MATRIXORIENTATION::TRUTHONVERTICALAXIS;
        } else {
            WriteWarningStatus("ConfigReader::ReadUnfoldingOptions", "You specified 'MatrixOrientation' option, but you didnt provide a valid config. Setting to TRUTHONHORIZONTAL.");
            fFitter->fMatrixOrientation = FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS;
        }
    }

    param = confSet->Get("TruthDistributionPath");
    if (param != "") {
        fFitter->fTruthDistributionPath = RemoveQuotes(param);
    }

    param = confSet->Get("TruthDistributionFile");
    if (param != "") {
        fFitter->fTruthDistributionFile = RemoveQuotes(param);
    }

    param = confSet->Get("TruthDistributionName");
    if (param != "") {
        fFitter->fTruthDistributionName = RemoveQuotes(param);
    }

    param = confSet->Get("NumberOfTruthBins");
    if (param == "") {
        WriteErrorStatus("ConfigReader::ReadUnfoldingOptions", "You need to set the number of truth bins!");
        return 1;
    } else {
        const int bins = std::stoi(param);
        if (bins < 2 || bins > 100) {
            WriteErrorStatus("ConfigReader::ReadUnfoldingOptions", "You set the number of truth bins which is < 2 or > 100, that loooks wrong");
            return 1;
        }
        fFitter->fNumberUnfoldingTruthBins = bins;
    }

    param = confSet->Get("NumberOfRecoBins");
    if (param != "") {
        const int bins = std::stoi(param);
        if (bins < 2 || bins > 100) {
            WriteErrorStatus("ConfigReader::ReadUnfoldingOptions", "You set the number of reco bins which is < 2 or > 100, that loooks wrong");
            return 1;
        }
        fFitter->fNumberUnfoldingRecoBins = bins;
    }

    param = confSet->Get("Tau");
    if (param != "") {
        const std::vector<std::string>& tmp = Vectorize(param, ',');
        if (tmp.size()==1) {
            fTaus.emplace_back(-1,std::stof(param));
        }
        else {
            for (auto& i : tmp) {
                const std::vector<std::string>& oneTau = Vectorize(i, ':');
                if (oneTau.size() != 2) {
                    WriteErrorStatus("ConfigReader::ReadUnfoldingOptions", "Wrong format for Tau!");
                    return 1;
                }

                int bin = std::stoi(oneTau.at(0));
                double value = std::stof(oneTau.at(1));
                fTaus.emplace_back(bin, value);
            }
        }
    }

    param = confSet->Get("UnfoldingResultMin");
    if (param != "") {
        fFitter->fUnfoldingResultMin = std::stod(param);
    }

    param = confSet->Get("UnfoldingResultMax");
    if (param != "") {
        fFitter->fUnfoldingResultMax = std::stod(param);
    }

    param = confSet->Get("TitleX");
    if (param != "") {
        fFitter->fUnfoldingTitleX = RemoveQuotes(param);
    }

    param = confSet->Get("TitleY");
    if (param != "") {
        fFitter->fUnfoldingTitleY = RemoveQuotes(param);
    }

    param = confSet->Get("TitleOffsetX");
    if (param != "") {
        fFitter->fUnfoldingTitleOffsetX = std::stod(param);
    }

    param = confSet->Get("TitleOffsetY");
    if (param != "") {
        fFitter->fUnfoldingTitleOffsetY = std::stod(param);
    }

    param = confSet->Get("RatioYmax");
    if (param != "") {
        fFitter->fUnfoldingRatioYmax = std::stod(param);
    }

    param = confSet->Get("RatioYmin");
    if (param != "") {
        fFitter->fUnfoldingRatioYmin = std::stod(param);
    }

    param = confSet->Get("LogX");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fUnfoldingLogX = true;
        } else {
            fFitter->fUnfoldingLogX = false;
        }
    }

    param = confSet->Get("LogY");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fUnfoldingLogY = true;
        } else {
            fFitter->fUnfoldingLogY = false;
        }
    }

    param = confSet->Get("MigrationTitleX");
    if (param != "") {
        fFitter->fMigrationTitleX = RemoveQuotes(param);
    }

    param = confSet->Get("MigrationTitleY");
    if (param != "") {
        fFitter->fMigrationTitleY = RemoveQuotes(param);
    }

    param = confSet->Get("MigrationLogX");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fMigrationLogX = true;
        } else {
            fFitter->fMigrationLogX = false;
        }
    }

    param = confSet->Get("MigrationLogY");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fMigrationLogY = true;
        } else {
            fFitter->fMigrationLogY = false;
        }
    }

    param = confSet->Get("MigrationTitleOffsetX");
    if (param != "") {
        fFitter->fMigrationTitleOffsetX = std::stod(param);
    }

    param = confSet->Get("MigrationTitleOffsetY");
    if (param != "") {
        fFitter->fMigrationTitleOffsetY = std::stod(param);
    }

    param = confSet->Get("MigrationZmin");
    if (param != "") {
        fFitter->fMigrationZmin = std::stod(param);
    }

    param = confSet->Get("MigrationZmax");
    if (param != "") {
        fFitter->fMigrationZmax = std::stod(param);
    }

    param = confSet->Get("ResponseZmin");
    if (param != "") {
        fFitter->fResponseZmin = std::stod(param);
    }

    param = confSet->Get("ResponseZmax");
    if (param != "") {
        fFitter->fResponseZmax = std::stod(param);
    }

    param = confSet->Get("PlotSystematicMigrations");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fPlotSystematicMigrations = true;
        } else {
            fFitter->fPlotSystematicMigrations = false;
        }
    }

    param = confSet->Get("MigrationText");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fFitter->fMigrationText = true;
        } else {
            fFitter->fMigrationText = false;
        }
    }

    param = confSet->Get("NominalTruthSample");
    if (param == "") {
        WriteErrorStatus("ConfigReader::ReadUnfoldingOptions", "You need to set NominalTruthSample option!");
        return 1;
    } else {
        fFitter->fNominalTruthSample = RemoveQuotes(param);
    }

    param = confSet->Get("AlternativeAsimovTruthSample");
    if (param != "") {
        fFitter->fAlternativeAsimovTruthSample = RemoveQuotes(param);
        fFitter->fFitIsBlind = false;
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadTruthSamples() {

    int isample(0);

    bool found(false);
    bool foundAlternative(false);

    std::vector<std::string> names;

    while(true) {
        ConfigSet *confSet = fParser->GetConfigSet("TruthSample",isample);
        if (!confSet) break;
        ++isample;

        auto sample = std::make_unique<TruthSample>(RemoveQuotes(confSet->GetValue()));
        if (sample->GetName() == fFitter->fNominalTruthSample) {
            found = true;
        }
        if (sample->GetName() == fFitter->fAlternativeAsimovTruthSample) {
            foundAlternative = true;
        }

        if (std::find(names.begin(), names.end(), RemoveQuotes(confSet->GetValue())) == names.end()) {
            names.emplace_back(RemoveQuotes(confSet->GetValue()));
        } else {
            WriteErrorStatus("ConfigReader::ReadTruthSamples", "Multiply defined TruthSample: " + RemoveQuotes(confSet->GetValue()));
            return 1;
        }

        std::string param = confSet->Get("Title");
        if (param != "") {
            sample->SetTitle(RemoveQuotes(param));
        }

        param = confSet->Get("LineStyle");
        if (param != "") {
            sample->SetLineStyle(std::stoi(param));
        }

        param = confSet->Get("LineColor");
        if (param != "") {
            sample->SetLineColor(std::stoi(param));
        }

        param = confSet->Get("UseForPlotting");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "TRUE") {
                sample->SetUseForPlotting(true);
            } else if (param == "FALSE") {
                sample->SetUseForPlotting(false);
            }
        }

        param = confSet->Get("TruthDistributionPath");
        if (param != "") {
            sample->SetTruthDistributionPath(RemoveQuotes(param));
        }

        param = confSet->Get("TruthDistributionFile");
        if (param != "") {
            sample->SetTruthDistributionFile(RemoveQuotes(param));
        }

        param = confSet->Get("TruthDistributionName");
        if (param != "") {
            sample->SetTruthDistributionName(RemoveQuotes(param));
        }

        fFitter->fTruthSamples.emplace_back(std::move(sample));
    }

    if (isample != 0 && !found) {
        WriteErrorStatus("ConfigReader::ReadTruthSamples", "The NominalTruthSample not found in any of the TruthSample");
        return 1;
    }

    if (fFitter->fAlternativeAsimovTruthSample != "" && !foundAlternative) {
        WriteErrorStatus("ConfigReader::ReadTruthSamples", "The AlternativeAsimovTruthSample is set but doesnt match any of the TruthSamples");
        return 1;
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadUnfoldingSampleOptions() {

    int isample = 0;
    while(true) {
        ConfigSet *confSet = fParser->GetConfigSet("UnfoldingSample",isample);
        if (confSet == nullptr) break;
        fHasAtLeastOneValidSample = true;
        ++isample;

        auto sample = std::make_unique<UnfoldingSample>();
        sample->SetName(RemoveQuotes(confSet->GetValue()));

        std::string param = confSet->Get("Title");
        if (param != "") {
            sample->SetTitle(RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixFile");
        if (param != "") {
            sample->fResponseMatrixFiles.clear();
            sample->fResponseMatrixFiles.emplace_back(RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixFiles");
        if (param != "") {
            sample->fResponseMatrixFiles = Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixName");
        if (param != "") {
            sample->fResponseMatrixNames.clear();
            sample->fResponseMatrixNames.emplace_back(RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixNames");
        if (param != "") {
            sample->fResponseMatrixNames = Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPath");
        if (param != "") {
            sample->fResponseMatrixPaths.clear();
            sample->fResponseMatrixPaths.emplace_back(RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPaths");
        if (param != "") {
            sample->fResponseMatrixPaths = Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixFileSuff");
        if (param != "") {
            sample->fResponseMatrixFileSuffs.clear();
            sample->fResponseMatrixFileSuffs.emplace_back(RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixFileSuffs");
        if (param != "") {
            sample->fResponseMatrixFileSuffs = Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixNameSuff");
        if (param != "") {
            sample->fResponseMatrixNameSuffs.clear();
            sample->fResponseMatrixNameSuffs.emplace_back(RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixNameSuffs");
        if (param != "") {
            sample->fResponseMatrixNameSuffs = Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPathSuff");
        if (param != "") {
            sample->fResponseMatrixPathSuffs.clear();
            sample->fResponseMatrixPathSuffs.emplace_back(RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPathSuffs");
        if (param != "") {
            sample->fResponseMatrixPathSuffs = Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("AcceptanceFile");
        if (param != "") {
            sample->fAcceptanceFiles.clear();
            sample->fAcceptanceFiles.emplace_back(RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFiles");
        if (param != "") {
            sample->fAcceptanceFiles = Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceName");
        if (param != "") {
            sample->fAcceptanceNames.clear();
            sample->fAcceptanceNames.emplace_back(RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNames");
        if (param != "") {
            sample->fAcceptanceNames = Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePath");
        if (param != "") {
            sample->fAcceptancePaths.clear();
            sample->fAcceptancePaths.emplace_back(RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePaths");
        if (param != "") {
            sample->fAcceptancePaths = Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuff");
        if (param != "") {
            sample->fAcceptanceFileSuffs.clear();
            sample->fAcceptanceFileSuffs.emplace_back(RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffs");
        if (param != "") {
            sample->fAcceptanceFileSuffs = Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuff");
        if (param != "") {
            sample->fAcceptanceNameSuffs.clear();
            sample->fAcceptanceNameSuffs.emplace_back(RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffs");
        if (param != "") {
            sample->fAcceptanceNameSuffs = Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuff");
        if (param != "") {
            sample->fAcceptancePathSuffs.clear();
            sample->fAcceptancePathSuffs.emplace_back(RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffs");
        if (param != "") {
            sample->fAcceptancePathSuffs = Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("SelectionEffFile");
        if (param != "") {
            sample->fSelectionEffFiles.clear();
            sample->fSelectionEffFiles.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFiles");
        if (param != "") {
            sample->fSelectionEffFiles = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffName");
        if (param != "") {
            sample->fSelectionEffNames.clear();
            sample->fSelectionEffNames.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNames");
        if (param != "") {
            sample->fSelectionEffNames = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPath");
        if (param != "") {
            sample->fSelectionEffPaths.clear();
            sample->fSelectionEffPaths.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPaths");
        if (param != "") {
            sample->fSelectionEffPaths = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFileSuff");
        if (param != "") {
            sample->fSelectionEffFileSuffs.clear();
            sample->fSelectionEffFileSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFileSuffs");
        if (param != "") {
            sample->fSelectionEffFileSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffNameSuff");
        if (param != "") {
            sample->fSelectionEffNameSuffs.clear();
            sample->fSelectionEffNameSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNameSuffs");
        if (param != "") {
            sample->fSelectionEffNameSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPathSuff");
        if (param != "") {
            sample->fSelectionEffPathSuffs.clear();
            sample->fSelectionEffPathSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathSuffs");
        if (param != "") {
            sample->fSelectionEffPathSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFile");
        if (param != "") {
            sample->fMigrationFiles.clear();
            sample->fMigrationFiles.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFiles");
        if (param != "") {
            sample->fMigrationFiles = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationName");
        if (param != "") {
            sample->fMigrationNames.clear();
            sample->fMigrationNames.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNames");
        if (param != "") {
            sample->fMigrationNames = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPath");
        if (param != "") {
            sample->fMigrationPaths.clear();
            sample->fMigrationPaths.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPaths");
        if (param != "") {
            sample->fMigrationPaths = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFileSuff");
        if (param != "") {
            sample->fMigrationFileSuffs.clear();
            sample->fMigrationFileSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFileSuffs");
        if (param != "") {
            sample->fMigrationFileSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationNameSuff");
        if (param != "") {
            sample->fMigrationNameSuffs.clear();
            sample->fMigrationNameSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNameSuffs");
        if (param != "") {
            sample->fMigrationNameSuffs = Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPathSuff");
        if (param != "") {
            sample->fMigrationPathSuffs.clear();
            sample->fMigrationPathSuffs.emplace_back(RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathSuffs");
        if (param != "") {
            sample->fMigrationPathSuffs = Vectorize(param, ',');
        }

        // Set FillColor
        param = confSet->Get("FillColor");
        if(param != "") sample->SetFillColor(std::atoi(param.c_str()));

        // Set LineColor
        param = confSet->Get("LineColor");
        if(param != "") sample->SetLineColor(std::atoi(param.c_str()));

        // Set Regions
        param = confSet->Get("Regions");
        if(param == "") {
            sample->fRegions.emplace_back("all");
        } else {
            sample->fRegions = Vectorize(param, ',');
        }

        fFitter->fUnfoldingSamples.emplace_back(std::move(sample));

    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadUnfoldingSystematicOptions() {

    int isyst(0);

    while(true) {
        ConfigSet *confSet = fParser->GetConfigSet("UnfoldingSystematic", isyst);
        if (confSet == nullptr) break;

        ++isyst;

        auto syst = std::make_unique<UnfoldingSystematic>();
        syst->SetName(RemoveQuotes(confSet->GetValue()));

        // Set NuisanceParameter
        std::string param = confSet->Get("NuisanceParameter");
        if(param != ""){
            syst->fNuisanceParameter = RemoveQuotes(param);
            TRExFitter::NPMAP[syst->GetName()] = syst->fNuisanceParameter;
        } else {
            syst->fNuisanceParameter = syst->GetName();
            TRExFitter::NPMAP[syst->GetName()] = syst->GetName();
        }

        param = confSet->Get("Type");
        if (param == "") {
            syst->SetType(Systematic::HISTO);
        } else {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);

            if(param == "OVERALL") syst->SetType(Systematic::OVERALL);
            else if(param == "SHAPE") syst->SetType(Systematic::SHAPE);
            else if (param == "STAT") syst->SetType(Systematic::STAT);
            else if (param == "HISTO") syst->SetType(Systematic::HISTO);
        }

        param = confSet->Get("Samples");
        if (param == "") {
            syst->fSamples.emplace_back("all");
        } else {
            syst->fSamples = Vectorize(param, ',');
        }

        param = confSet->Get("Regions");
        if (param == "") {
            syst->fRegions.emplace_back("all");
        } else {
            syst->fRegions = Vectorize(param, ',');
        }

        param = confSet->Get("Title");
        if (param == "") {
            syst->SetTitle(RemoveQuotes(confSet->GetValue()));
        } else {
            syst->SetTitle(RemoveQuotes(param));
        }

        // SetCategory
        param = confSet->Get("Category");
        if(param != ""){
            syst->fCategory = RemoveQuotes(param);
            syst->fCategory = RemoveQuotes(param); //SubCategory defaults to the Category setting, if the Category is explicitly set
        }

        // SetSubCategory
        param = confSet->Get("SubCategory");
        if (param != ""){
            syst->fSubCategory = RemoveQuotes(param); // note this needs to happen after Category was set, in order to overwrite the default if required
        }

        bool hasUp(false);
        bool hasDown(false);

        param = confSet->Get("ResponseMatrixPathUp");
        if (param != "") {
            syst->fResponseMatrixPathsUp.clear();
            syst->fResponseMatrixPathsUp.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathsUp");
        if (param != "") {
            syst->fResponseMatrixPathsUp = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathDown");
        if (param != "") {
            syst->fResponseMatrixPathsDown.clear();
            syst->fResponseMatrixPathsDown.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathsDown");
        if (param != "") {
            syst->fResponseMatrixPathsDown = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffUp");
        if (param != "") {
            syst->fResponseMatrixPathSuffsUp.clear();
            syst->fResponseMatrixPathSuffsUp.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffsUp");
        if (param != "") {
            syst->fResponseMatrixPathSuffsUp = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffDown");
        if (param != "") {
            syst->fResponseMatrixPathSuffsDown.clear();
            syst->fResponseMatrixPathSuffsDown.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffsDown");
        if (param != "") {
            syst->fResponseMatrixPathSuffsDown = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNameUp");
        if (param != "") {
            syst->fResponseMatrixNamesUp.clear();
            syst->fResponseMatrixNamesUp.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNamesUp");
        if (param != "") {
            syst->fResponseMatrixNamesUp = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNameDown");
        if (param != "") {
            syst->fResponseMatrixNamesDown.clear();
            syst->fResponseMatrixNamesDown.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNamesDown");
        if (param != "") {
            syst->fResponseMatrixNamesDown = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffUp");
        if (param != "") {
            syst->fResponseMatrixNameSuffsUp.clear();
            syst->fResponseMatrixNameSuffsUp.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffsUp");
        if (param != "") {
            syst->fResponseMatrixNameSuffsUp = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffDown");
        if (param != "") {
            syst->fResponseMatrixNameSuffsDown.clear();
            syst->fResponseMatrixNameSuffsDown.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffsDown");
        if (param != "") {
            syst->fResponseMatrixNameSuffsDown = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFileUp");
        if (param != "") {
            syst->fResponseMatrixFilesUp.clear();
            syst->fResponseMatrixFilesUp.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFilesUp");
        if (param != "") {
            syst->fResponseMatrixFilesUp = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFileDown");
        if (param != "") {
            syst->fResponseMatrixFilesDown.clear();
            syst->fResponseMatrixFilesDown.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFilesDown");
        if (param != "") {
            syst->fResponseMatrixFilesDown = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffUp");
        if (param != "") {
            syst->fResponseMatrixFileSuffsUp.clear();
            syst->fResponseMatrixFileSuffsUp.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffsUp");
        if (param != "") {
            syst->fResponseMatrixFileSuffsUp = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffDown");
        if (param != "") {
            syst->fResponseMatrixFileSuffsDown.clear();
            syst->fResponseMatrixFileSuffsDown.emplace_back(RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffsDown");
        if (param != "") {
            syst->fResponseMatrixFileSuffsDown = Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("AcceptancePathUp");
        if (param != "") {
            syst->fAcceptancePathsUp.clear();
            syst->fAcceptancePathsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathsUp");
        if (param != "") {
            syst->fAcceptancePathsUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathDown");
        if (param != "") {
            syst->fAcceptancePathsDown.clear();
            syst->fAcceptancePathsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathsDown");
        if (param != "") {
            syst->fAcceptancePathsDown = Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffUp");
        if (param != "") {
            syst->fAcceptancePathSuffsUp.clear();
            syst->fAcceptancePathSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffsUp");
        if (param != "") {
            syst->fAcceptancePathSuffsUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffDown");
        if (param != "") {
            syst->fAcceptancePathSuffsDown.clear();
            syst->fAcceptancePathSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffsDown");
        if (param != "") {
            syst->fAcceptancePathSuffsDown = Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameUp");
        if (param != "") {
            syst->fAcceptanceNamesUp.clear();
            syst->fAcceptanceNamesUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNamesUp");
        if (param != "") {
            syst->fAcceptanceNamesUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameDown");
        if (param != "") {
            syst->fAcceptanceNamesDown.clear();
            syst->fAcceptanceNamesDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNamesDown");
        if (param != "") {
            syst->fAcceptanceNamesDown = Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffUp");
        if (param != "") {
            syst->fAcceptanceNameSuffsUp.clear();
            syst->fAcceptanceNameSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffsUp");
        if (param != "") {
            syst->fAcceptanceNameSuffsUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffDown");
        if (param != "") {
            syst->fAcceptanceNameSuffsDown.clear();
            syst->fAcceptanceNameSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffsDown");
        if (param != "") {
            syst->fAcceptanceNameSuffsDown = Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileUp");
        if (param != "") {
            syst->fAcceptanceFilesUp.clear();
            syst->fAcceptanceFilesUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFilesUp");
        if (param != "") {
            syst->fAcceptanceFilesUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileDown");
        if (param != "") {
            syst->fAcceptanceFilesDown.clear();
            syst->fAcceptanceFilesDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFilesDown");
        if (param != "") {
            syst->fAcceptanceFilesDown = Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffUp");
        if (param != "") {
            syst->fAcceptanceFileSuffsUp.clear();
            syst->fAcceptanceFileSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffsUp");
        if (param != "") {
            syst->fAcceptanceFileSuffsUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffDown");
        if (param != "") {
            syst->fAcceptanceFileSuffsDown.clear();
            syst->fAcceptanceFileSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffsDown");
        if (param != "") {
            syst->fAcceptanceFileSuffsDown = Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }
        param = confSet->Get("SelectionEffPathUp");
        if (param != "") {
            syst->fSelectionEffPathsUp.clear();
            syst->fSelectionEffPathsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("SelectionEffPathsUp");
        if (param != "") {
            syst->fSelectionEffPathsUp = Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("SelectionEffPathDown");
        if (param != "") {
            syst->fSelectionEffPathsDown.clear();
            syst->fSelectionEffPathsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffPathsDown");
        if (param != "") {
            syst->fSelectionEffPathsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffPathSuffUp");
        if (param != "") {
            syst->fSelectionEffPathSuffsUp.clear();
            syst->fSelectionEffPathSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffPathSuffsUp");
        if (param != "") {
            syst->fSelectionEffPathSuffsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffPathSuffDown");
        if (param != "") {
            syst->fSelectionEffPathSuffsDown.clear();
            syst->fSelectionEffPathSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffPathSuffsDown");
        if (param != "") {
            syst->fSelectionEffPathSuffsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNameUp");
        if (param != "") {
            syst->fSelectionEffNamesUp.clear();
            syst->fSelectionEffNamesUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNamesUp");
        if (param != "") {
            syst->fSelectionEffNamesUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNameDown");
        if (param != "") {
            syst->fSelectionEffNamesDown.clear();
            syst->fSelectionEffNamesDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNamesDown");
        if (param != "") {
            syst->fSelectionEffNamesDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNameSuffUp");
        if (param != "") {
            syst->fSelectionEffNameSuffsUp.clear();
            syst->fSelectionEffNameSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNameSuffsUp");
        if (param != "") {
            syst->fSelectionEffNameSuffsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNameSuffDown");
        if (param != "") {
            syst->fSelectionEffNameSuffsDown.clear();
            syst->fSelectionEffNameSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNameSuffsDown");
        if (param != "") {
            syst->fSelectionEffNameSuffsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFileUp");
        if (param != "") {
            syst->fSelectionEffFilesUp.clear();
            syst->fSelectionEffFilesUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFilesUp");
        if (param != "") {
            syst->fSelectionEffFilesUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFileDown");
        if (param != "") {
            syst->fSelectionEffFilesDown.clear();
            syst->fSelectionEffFilesDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFilesDown");
        if (param != "") {
            syst->fSelectionEffFilesDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFileSuffUp");
        if (param != "") {
            syst->fSelectionEffFileSuffsUp.clear();
            syst->fSelectionEffFileSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFileSuffsUp");
        if (param != "") {
            syst->fSelectionEffFileSuffsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFileSuffDown");
        if (param != "") {
            syst->fSelectionEffFileSuffsDown.clear();
            syst->fSelectionEffFileSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFileSuffsDown");
        if (param != "") {
            syst->fSelectionEffFileSuffsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationPathUp");
        if (param != "") {
            syst->fMigrationPathsUp.clear();
            syst->fMigrationPathsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationPathsUp");
        if (param != "") {
            syst->fMigrationPathsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationPathDown");
        if (param != "") {
            syst->fMigrationPathsDown.clear();
            syst->fMigrationPathsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationPathsDown");
        if (param != "") {
            syst->fMigrationPathsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationPathSuffUp");
        if (param != "") {
            syst->fMigrationPathSuffsUp.clear();
            syst->fMigrationPathSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationPathSuffsUp");
        if (param != "") {
            syst->fMigrationPathSuffsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationPathSuffDown");
        if (param != "") {
            syst->fMigrationPathSuffsDown.clear();
            syst->fMigrationPathSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationPathSuffsDown");
        if (param != "") {
            syst->fMigrationPathSuffsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationNameUp");
        if (param != "") {
            syst->fMigrationNamesUp.clear();
            syst->fMigrationNamesUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationNamesUp");
        if (param != "") {
            syst->fMigrationNamesUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationNameDown");
        if (param != "") {
            syst->fMigrationNamesDown.clear();
            syst->fMigrationNamesDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationNamesDown");
        if (param != "") {
            syst->fMigrationNamesDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationNameSuffUp");
        if (param != "") {
            syst->fMigrationNameSuffsUp.clear();
            syst->fMigrationNameSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationNameSuffsUp");
        if (param != "") {
            syst->fMigrationNameSuffsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationNameSuffDown");
        if (param != "") {
            syst->fMigrationNameSuffsDown.clear();
            syst->fMigrationNameSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationNameSuffsDown");
        if (param != "") {
            syst->fMigrationNameSuffsDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationFileUp");
        if (param != "") {
            syst->fMigrationFilesUp.clear();
            syst->fMigrationFilesUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationFilesUp");
        if (param != "") {
            syst->fMigrationFilesUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationFileDown");
        if (param != "") {
            syst->fMigrationFilesDown.clear();
            syst->fMigrationFilesDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationFilesDown");
        if (param != "") {
            syst->fMigrationFilesDown = Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationFileSuffUp");
        if (param != "") {
            syst->fMigrationFileSuffsUp.clear();
            syst->fMigrationFileSuffsUp.emplace_back(RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationFileSuffsUp");
        if (param != "") {
            syst->fMigrationFileSuffsUp = Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationFileSuffDown");
        if (param != "") {
            syst->fMigrationFileSuffsDown.clear();
            syst->fMigrationFileSuffsDown.emplace_back(RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationFileSuffsDown");
        if (param != "") {
            syst->fMigrationFileSuffsDown = Vectorize(param,',');
            hasDown = true;
        }

        syst->fHasUpVariation = hasUp;
        syst->fHasDownVariation = hasDown;

        // Set Symmetrisation
        param = confSet->Get("Symmetrisation");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "ONESIDED"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEONESIDED);
            }
            else if(param == "TWOSIDED"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZETWOSIDED);
            }
            else if(param == "ABSMEAN"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEABSMEAN);
            }
            else if(param == "MAXIMUM"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEMAXIMUM);
            }
            else {
                WriteErrorStatus("ConfigReader::ReadUnfoldingSystematicOptions", "Symetrisation scheme is not recognized ... ");
                return 1;
            }
        }

        param = confSet->Get("SmoothingOption");
        if(param != ""){
            syst->fSampleSmoothing = true;
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param == "MAXVARIATION" ) syst->SetSmoothOption(HistoTools::SmoothOption::MAXVARIATION);
            else if (param == "TTBARRESONANCE") syst->SetSmoothOption(HistoTools::SmoothOption::TTBARRESONANCE);
            else if (param == "COMMONTOOLSMOOTHMONOTONIC") syst->SetSmoothOption(HistoTools::SmoothOption::COMMONTOOLSMOOTHMONOTONIC);
            else if (param == "COMMONTOOLSMOOTHPARABOLIC") syst->SetSmoothOption(HistoTools::SmoothOption::COMMONTOOLSMOOTHPARABOLIC);
            else if (param == "KERNELRATIOUNIFORM") syst->SetSmoothOption(HistoTools::SmoothOption::KERNELRATIOUNIFORM);
            else if (param == "KERNELDELTAGAUSS") syst->SetSmoothOption(HistoTools::SmoothOption::KERNELDELTAGAUSS);
            else if (param == "KERNELRATIOGAUSS") syst->SetSmoothOption(HistoTools::SmoothOption::KERNELRATIOGAUSS);
            else {
                WriteWarningStatus("ConfigReader::ReadUnfoldingSystematicOptions", "You specified 'SmoothingOption' option but you didn't provide valid input. Using default from job block");
                syst->fSampleSmoothing = false;
            }
        }

        fFitter->fUnfoldingSystematics.emplace_back(std::move(syst));
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::UnfoldingCorrections() {
    if (fFitter->fFitType != TRExFit::FitType::UNFOLDING) return 0;

    int sc(0);

    // First process Samples
    sc += ProcessUnfoldingSamples();

    // Then process Systematics
    sc += ProcessUnfoldingSystematics();

    // Add norm factors
    sc += AddUnfoldingNormFactors();

    return sc;
}

//__________________________________________________________________________________
//
std::string ConfigReader::CheckName( std::string name ){
    name = RemoveQuotes(name);
    if( std::isdigit( name.at(0) ) ){
        WriteErrorStatus("ConfigReader::CheckName", "Failed to browse name: " + name + ". A number has been detected at the first position of the name.");
        WriteErrorStatus("ConfigReader::CheckName", "           This can lead to unexpected behaviour in HistFactory. Please change the name. ");
        exit(EXIT_FAILURE);
    } else {
        return RemoveQuotes(name);
    }
}

//__________________________________________________________________________________
//
bool ConfigReader::ConfigHasNTUP(ConfigSet* confSet){
    if (confSet->Get("Variable") != "" || confSet->Get("VariableForSample") != "" || confSet->Get("Selection") != "" || confSet->Get("NtupleName") != "" || confSet->Get("NtupleNameSuff") != "" || confSet->Get("MCweight") != "" || confSet->Get("NtuplePathSuff") != "" || confSet->Get("NtupleFile") != "" || confSet->Get("NtupleFiles") != "" || confSet->Get("NtupleNames") != "" || confSet->Get("NtuplePath") != "" || confSet->Get("NtuplePaths") != "" || confSet->Get("NtuplePathSuffs") != "") return true;
    else return false;
}

//__________________________________________________________________________________
//
bool ConfigReader::ConfigHasHIST(ConfigSet* confSet){
    if (confSet->Get("HistoFile") != "" || confSet->Get("HistoName") != "" || confSet->Get("HistoPathSuff") != "" || confSet->Get("HistoPathSuffs") != "" || confSet->Get("HistoPath") != "" ) return true;
    else return false;
}


//__________________________________________________________________________________
//
bool ConfigReader::CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2){
    for (const auto& i : v1){
        if (i == "") continue;
        std::string s = i;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "NONE") continue;
        if (s == "ALL") continue;
        if (Common::FindInStringVector(v2, i) < 0){
            return false;
        }
    }

    return true;
}

//__________________________________________________________________________________
//
bool ConfigReader::CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2, const std::vector<std::string> &v3){
    for (const auto& i : v1){
        if (i == "") continue;
        std::string s = i;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "NONE") continue;
        if (s == "ALL") continue;
        if (Common::FindInStringVector(v2, i) < 0){
            if (Common::FindInStringVector(v3, i) < 0){
                return false;
            }
        }
    }

    return true;
}


//__________________________________________________________________________________
//
std::vector<std::string> ConfigReader::GetAvailableRegions(){
    std::vector<std::string> availableRegions;

    int nReg = 0;
    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("Region",nReg);
        if (confSet == nullptr) break;

        nReg++;
        std::string tmp = RemoveQuotes(confSet->GetValue());
        availableRegions.emplace_back(CheckName(tmp));
    }
    return availableRegions;
}

//__________________________________________________________________________________
//
std::vector<std::string> ConfigReader::GetAvailableSamples(){
    std::vector<std::string> availableSamples;
    int nSample = 0;
    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("Sample",nSample);
        if (confSet == nullptr) break;

        nSample++;
        std::string tmp = RemoveQuotes(confSet->GetValue());
        availableSamples.emplace_back(CheckName(tmp));
    }
    return availableSamples;
}

//__________________________________________________________________________________
//
std::vector<std::string> ConfigReader::GetAvailableSysts(){
    std::vector<std::string> availableSysts;
    int nSyst = 0;
    while(true){
        ConfigSet *confSet = fParser->GetConfigSet("Systematic",nSyst);
        if (confSet == nullptr) break;

        nSyst++;
        std::string tmp = RemoveQuotes(confSet->GetValue());
        availableSysts.emplace_back(CheckName(tmp));
    }
    return availableSysts;
}

//__________________________________________________________________________________
//
bool ConfigReader::SystHasProblematicName(const std::string& name){
    if ((name.find("gamma") != std::string::npos) || (name.find("alpha") != std::string::npos)){
        WriteErrorStatus("ConfigReader::SystHasProblematicName", "NP " + name + " has a problematic name, please change it.");
        WriteErrorStatus("ConfigReader::SystHasProblematicName", "You should not be using names with: \"gamma\", \"alpha\" as these are used internally and can cause problems.");
        return true;
    }

    return false;
}

//__________________________________________________________________________________
//
int ConfigReader::ProcessUnfoldingSamples() {

    for (const auto& ireg : fFitter->fRegions) {
        if (ireg->fRegionType != Region::RegionType::SIGNAL) continue;

        for (const auto& isample : fFitter->fUnfoldingSamples) {
            if(isample->fRegions[0] != "all" &&
                Common::FindInStringVector(isample->fRegions, ireg->fName) < 0) continue;

            // Convert the UnfoldingSample to sample and adjust the paths
            const std::vector<Sample*> samples = isample->ConvertToSample(ireg, fFitter->fNumberUnfoldingTruthBins, fFitter->fName);
            fFitter->fSamples.insert(fFitter->fSamples.end(), samples.begin(), samples.end());
            fFitter->fNSamples += fFitter->fNumberUnfoldingTruthBins;
        }

        // add custom asimov sample
        if (fFitter->fAlternativeAsimovTruthSample != "") {
            Sample* sample = new Sample("AlternativeSignal_"+ireg->fName, Sample::SampleType::GHOST);
            sample->fRegions = Common::ToVec(ireg->fName);
            sample->fHistoPaths = Common::ToVec(fFitter->fName+"/UnfoldingHistograms");
            sample->fHistoFiles = Common::ToVec("FoldedHistograms");
            const std::string histoName = ireg->fName + "_AlternativeAsimov";
            sample->fHistoNames = Common::ToVec(histoName);
            fFitter->fSamples.emplace_back(std::move(sample));
            ++fFitter->fNSamples;
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::ProcessUnfoldingSystematics() {
    for (const auto& ireg : fFitter->fRegions) {
        if (ireg->fRegionType != Region::RegionType::SIGNAL) continue;

        for (const auto& isyst : fFitter->fUnfoldingSystematics) {
            if(isyst->fRegions[0] != "all" &&
                Common::FindInStringVector(isyst->fRegions, ireg->fName) < 0) continue;

            if (fFitter->fUnfoldingSamples.size() == 0) {
                WriteErrorStatus("ConfigReader::ProcessUnfoldingSystematics", "No UnfoldingSamples set!");
                return 1;
            }
            std::string unfoldingSampleName("");
            for (const auto& isample : fFitter->fUnfoldingSamples) {
                if(isample->fRegions[0] != "all" &&
                    Common::FindInStringVector(isample->fRegions, ireg->fName) < 0) continue;

                unfoldingSampleName = isample->GetName();
            }

            if (unfoldingSampleName == "") {
                WriteErrorStatus("ConfigReader::ProcessUnfoldingSystematics", "Sample name not set!");
                exit(EXIT_FAILURE);
            }

            const std::vector<Systematic*> systs = isyst->ConvertToSystematic(ireg,
                                                                              fFitter->fNumberUnfoldingTruthBins,
                                                                              fFitter->fName,
                                                                              unfoldingSampleName,
                                                                              fFitter->fSamples);
            fFitter->fSystematics.insert(fFitter->fSystematics.end(), systs.begin(), systs.end());
            fFitter->fNSyst += fFitter->fNumberUnfoldingTruthBins;
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReader::AddUnfoldingNormFactors() {

    for (int i = 0; i < fFitter->fNumberUnfoldingTruthBins; ++i) {
        const std::string name = "Bin_" + std::to_string(i+1);
        NormFactor* nf = new NormFactor(name);

        TRExFitter::SYSTMAP[nf->fName] = nf->fName;

        if(Common::FindInStringVector(fFitter->fNormFactorNames, nf->fName) < 0) {
           fFitter->fNormFactors.emplace_back(nf);
           fFitter->fNormFactorNames.emplace_back(nf->fName);
           fFitter->fNNorm++;
        }

        nf->fNuisanceParameter = nf->fName;
        TRExFitter::NPMAP[nf->fName] = nf->fName;
        nf->fCategory = "TruthBins";
        nf->fTitle = "UnfoldedTruthBin_"+std::to_string(i);
        nf->fMin = fFitter->fUnfoldingResultMin;
        nf->fMax = fFitter->fUnfoldingResultMax;
        nf->fNominal = 1;
        nf->fRegions = GetAvailableRegions();
        std::vector<std::string> sampleNames;
        for (const auto& ireg : fFitter->fRegions) {
            if (ireg->fRegionType == Region::RegionType::SIGNAL) {
                const std::string sampleName = ireg->fName + "_Truth_bin_" + std::to_string(i+1);
                sampleNames.emplace_back(sampleName);
            }
        }
        for(auto& isample : fFitter->fSamples) {
            if (std::find(sampleNames.begin(), sampleNames.end(), isample->fName) == sampleNames.end()) continue;
            isample->AddNormFactor(nf);
        }

        // check if the bin is in the list specifies taus
        bool hasTau = false;
        double tauValue = 0.;
        if (fTaus.size()==1) {
            if (fTaus[0].first==-1){
                hasTau = true;
                tauValue = fTaus[0].second;
            }
        }
        else {
            auto it = std::find_if(fTaus.begin(), fTaus.end(),
            [&i](const std::pair<int, double>& element){ return element.first == i+1;});
            if (it != fTaus.end()) {
                hasTau = true;
                tauValue = it->second;
            }
        }
        if (hasTau) {
            WriteInfoStatus("ConfigReader::AddUnfoldingNormFactors", "Setting truth bin: " + std::to_string(i+1) + " to use tau = " + std::to_string(tauValue));
            nf->fTau = tauValue;
        }
    }

    return 0;
}

//__________________________________________________________________________________
//
void ConfigReader::FixReferenceSamples(Systematic* sys) const {

    if (sys->fHistoPathsUpRefSample.size()     == 0) sys->fHistoPathsUpRefSample     = sys->fHistoPathsUp;
    if (sys->fHistoPathsDownRefSample.size()   == 0) sys->fHistoPathsDownRefSample   = sys->fHistoPathsDown;
    if (sys->fHistoPathSufUpRefSample.size()   == 0) sys->fHistoPathSufUpRefSample   = sys->fHistoPathSufUp;
    if (sys->fHistoPathSufDownRefSample.size() == 0) sys->fHistoPathSufDownRefSample = sys->fHistoPathSufDown;
    if (sys->fHistoFilesUpRefSample.size()     == 0) sys->fHistoFilesUpRefSample     = sys->fHistoFilesUp;
    if (sys->fHistoFilesDownRefSample.size()   == 0) sys->fHistoFilesDownRefSample   = sys->fHistoFilesDown;
    if (sys->fHistoFileSufUpRefSample.size()   == 0) sys->fHistoFileSufUpRefSample   = sys->fHistoFileSufUp;
    if (sys->fHistoFileSufDownRefSample.size() == 0) sys->fHistoFileSufDownRefSample = sys->fHistoFileSufDown;
    if (sys->fHistoNamesUpRefSample.size()     == 0) sys->fHistoNamesUpRefSample     = sys->fHistoNamesUp;
    if (sys->fHistoNamesDownRefSample.size()   == 0) sys->fHistoNamesDownRefSample   = sys->fHistoNamesDown;
    if (sys->fHistoNameSufUpRefSample.size()   == 0) sys->fHistoNameSufUpRefSample   = sys->fHistoNameSufUp;
    if (sys->fHistoNameSufDownRefSample.size() == 0) sys->fHistoNameSufDownRefSample = sys->fHistoNameSufDown;

    if (sys->fNtuplePathsUpRefSample.size()    == 0) sys->fNtuplePathsUpRefSample     = sys->fNtuplePathsUp;
    if (sys->fNtuplePathsDownRefSample.size()  == 0) sys->fNtuplePathsDownRefSample   = sys->fNtuplePathsDown;
    if (sys->fNtuplePathSufUpRefSample.size()  == 0) sys->fNtuplePathSufUpRefSample   = sys->fNtuplePathSufUp;
    if (sys->fNtuplePathSufDownRefSample.size()== 0) sys->fNtuplePathSufDownRefSample = sys->fNtuplePathSufDown;
    if (sys->fNtupleFilesUpRefSample.size()    == 0) sys->fNtupleFilesUpRefSample     = sys->fNtupleFilesUp;
    if (sys->fNtupleFilesDownRefSample.size()  == 0) sys->fNtupleFilesDownRefSample   = sys->fNtupleFilesDown;
    if (sys->fNtupleFileSufUpRefSample.size()  == 0) sys->fNtupleFileSufUpRefSample   = sys->fNtupleFileSufUp;
    if (sys->fNtupleFileSufDownRefSample.size()== 0) sys->fNtupleFileSufDownRefSample = sys->fNtupleFileSufDown;
    if (sys->fNtupleNamesUpRefSample.size()    == 0) sys->fNtupleNamesUpRefSample     = sys->fNtupleNamesUp;
    if (sys->fNtupleNamesDownRefSample.size()  == 0) sys->fNtupleNamesDownRefSample   = sys->fNtupleNamesDown;
    if (sys->fNtupleNameSufUpRefSample.size()  == 0) sys->fNtupleNameSufUpRefSample   = sys->fNtupleNameSufUp;
    if (sys->fNtupleNameSufDownRefSample.size()== 0) sys->fNtupleNameSufDownRefSample = sys->fNtupleNameSufDown;
}

