#include "TtHFitter/ConfigReader.h"

#include "TtHFitter/TtHFit.h"
#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/Common.h"
#include "TtHFitter/Region.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/NormFactor.h"
#include "TtHFitter/ShapeFactor.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/HistoTools.h"


ConfigReader::ConfigReader(TtHFit *fitter){
    fFitter = fitter;
    WriteInfoStatus("ConfigReader::ConfigReader", "Started reading the config");
}

ConfigReader::~ConfigReader(){
}


int ConfigReader::ReadFullConfig(const std::string& fileName, const std::string& option){
    // initialize ConfigParser for the actual config
    fParser.ReadFile(fileName);

    // initialize checker COnfigParser to cross check the input
    //ConfigParser refConfig;
    //refConfig.ReadFile("jobSchema.config");
    //int sc = fParser.CheckSyntax(&refConfig);
    int sc = 0;

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

    sc+= ReadRegionOptions();
    
    sc+= ReadSampleOptions();
    
    sc+= ReadNormFactorOptions();
    
    sc+= ReadShapeFactorOptions();
    
    sc+= ReadSystOptions();
    
    sc+= PostConfig();

    return sc;
}

int ConfigReader::ReadCommandLineOptions(std::string option){
    std::vector< std::string > optVec = Vectorize(option,':');
    std::map< std::string,std::string > optMap;
    
    for(const std::string& iopt : optVec){
        std::vector< std::string > optPair;
        optPair = Vectorize(iopt,'=');
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
        fOnlySignal = optMap["Signal"];
    }
    if(optMap["FitResults"]!=""){
        fFitter->fFitResultsFile = optMap["FitResults"];
    }
    if(optMap["FitType"]!=""){
        if(optMap["FitType"]=="SPLUSB") fFitter->SetFitType(TtHFit::SPLUSB);
        if(optMap["FitType"]=="BONLY")  fFitter->SetFitType(TtHFit::BONLY);
    }
    if(optMap["LumiScale"]!=""){
        fFitter->fLumiScale = atof(optMap["LumiScale"].c_str());
    }
    if(optMap["BootstrapIdx"]!=""){
        fFitter->fBootstrapIdx = atoi(optMap["BootstrapIdx"].c_str());
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
    if(fOnlySignal!=""){
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "  Only Signal: ");
        WriteInfoStatus("ConfigReader::ReadCommandLineOptions", "    " + fOnlySignal);
    }
    return 0;
}

int ConfigReader::ReadJobOptions(){
    std::string param = ""; // helper string

    ConfigSet *confSet = fParser.GetConfigSet("Job");
    if (confSet == nullptr){
        WriteErrorStatus("ConfigReader::ReadJobOptions", "You need to provide JOB settings!");
        return 1;
    }
    
    fFitter->fName = CheckName(confSet->GetValue());
    fFitter->fInputName = fFitter->fName;
    
    //Set DebugLevel
    param = confSet->Get("DebugLevel");
    if( param != "")  TtHFitter::SetDebugLevel( atoi(param.c_str()) );
    
    // Set outputDir
    param = confSet->Get("OutputDir");
    if(param != ""){
      fFitter->fDir = param;
      if(fFitter->fDir.back() != '/') fFitter->fDir += '/';
      fFitter->fName = fFitter->fDir + fFitter->fName;
      gSystem->mkdir((fFitter->fName).c_str(), true);
    }

    // Set Label
    param = confSet->Get("Label");
    if(param!="") fFitter->fLabel = param;
    else          fFitter->fLabel = fFitter->fName;

    // Set POI
    fFitter->SetPOI(CheckName(confSet->Get("POI")));

    //Set reading option
    param = confSet->Get("ReadFrom");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if(      param=="HIST" || param=="HISTOGRAMS")  fFitter->fInputType = 0;
    else if( param=="NTUP" || param=="NTUPLES" )    fFitter->fInputType = 1;
    else{
        WriteErrorStatus("ConfigReader::ReadJobOptions", "Invalid \"ReadFrom\" argument. Options: \"HIST\", \"NTUP\"");
        return 1;
    }
    
    // set default MERGEUNDEROVERFLOW
    if(fFitter->fInputType==0)      TtHFitter::MERGEUNDEROVERFLOW = false;
    else if(fFitter->fInputType==1) TtHFitter::MERGEUNDEROVERFLOW = true;
   
    // Set MergeUnderOverFlow from config 
    param = confSet->Get("MergeUnderOverFlow");
    if(param!=""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  TtHFitter::MERGEUNDEROVERFLOW = true;
        else if( param == "FALSE" ) TtHFitter::MERGEUNDEROVERFLOW = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'MergeUnderOverFlow' option but you didn't provide valid setting. Using default (FALSE)");
            TtHFitter::MERGEUNDEROVERFLOW = false;
        }
    }
    
    //Set paths
    // HIST option only
    if(fFitter->fInputType==0){
        fFitter->AddHistoPath( confSet->Get("HistoPath") );
        if (confSet->Get("NtuplePath") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtuplePath:' option, ignoring.");
        }
        if (confSet->Get("NtuplePaths") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtuplePaths:' option, ignoring.");
        }
        if (confSet->Get("MCweight") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'MCweight:' option, ignoring.");
        }
        if (confSet->Get("Selection") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'Selection:' option, ignoring.");
        }
        if (confSet->Get("NtupleName") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtupleName:' option, ignoring.");
        }
    }
    // Setting for NTUP only
    if(fFitter->fInputType==1){
        if (confSet->Get("HistoPath") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NTUP option but you provided 'HistoPath:' option, ignoring.");
        }
        fFitter->SetNtupleFile( confSet->Get("NtupleFile") );
        if(confSet->Get("NtuplePath")!="") {
            fFitter->AddNtuplePath( confSet->Get("NtuplePath") ); 
        }
        param = confSet->Get("NtuplePaths");
        if( param != "" ){
            std::vector<std::string> paths = Vectorize( param,',' );
            for(const std::string& ipath : paths){
                fFitter->AddNtuplePath( ipath );
            }
        }
        param = confSet->Get("MCweight");
        if(param!="") fFitter->SetMCweight(param);

        param = confSet->Get("Selection");
        if(param!="") fFitter->SetSelection(param);
        fFitter->SetNtupleName( confSet->Get("NtupleName") );
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

    // Set TtresSmoothing
    param = confSet->Get("TtresSmoothing");
    if( param != ""){ 
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fTtresSmoothing = true;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'TtresSmoothing' option but you didn't set it to TRUE. Using default (FALSE)");
            fFitter->fTtresSmoothing = false;
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
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'MCstatThreshold' option but you didn't set it to NONE. Using default (TRUE, 0)");
            fFitter->SetStatErrorConfig( true, 0.);
        }
    }
    else{
        fFitter->SetStatErrorConfig( true, 0. );
    }

    //Set MCstatConstraint
    param = confSet->Get("MCstatConstraint");
    if( param != "")  fFitter->fStatErrCons = param;

    // Set UseGammaPulls
    param = confSet->Get("UseGammaPulls");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") fFitter->fUseGammaPulls = true;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'UseGammaPulls' option but you didn't set it to TRUE. Using default (FALSE)");
            fFitter->fUseGammaPulls = false;
        }
    }

    // plotting options are in special function
    if (SetJobPlot(confSet) != 0) return 1;

    // Set TableOptions 
    param = confSet->Get("TableOptions");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fFitter->fTableOptions = param;
    }

    // Set CorrelationThreshold
    param = confSet->Get("CorrelationThreshold");
    if( param != ""){
        TtHFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }

    // Set HistoChecks
    param = confSet->Get("HistoChecks");
    if(param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "NOCRASH" ){
            TtHFitter::HISTOCHECKCRASH = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'HistoChecks' option but you didn't set it to NOCRASH.");
        }
    }

    // Set LumiLabel
    param = confSet->Get("LumiLabel");
    if( param != "") fFitter->fLumiLabel = param;

    // Set CmeLabel
    param = confSet->Get("CmeLabel");
    if( param != "") fFitter->fCmeLabel = param;

    // Set SplitHistoFiles
    param = confSet->Get("SplitHistoFiles");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TtHFitter::SPLITHISTOFILES = true;
        } else if (param == "FALSE") {
            TtHFitter::SPLITHISTOFILES = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'SplitHistoFiles' option but you didn't provide valid setting. Using default (false)");
            TtHFitter::SPLITHISTOFILES = false;
        }
    }

    // Set BlindingThreshold"
    param = confSet->Get("BlindingThreshold");
    if( param != ""){
        fFitter->fBlindingThreshold = atof(param.c_str());
    }

    // Set RankingMaxNP
    param = confSet->Get("RankingMaxNP");
    if( param != ""){
        fFitter->fRankingMaxNP = atoi(param.c_str());
    }

    // Set ReduceNPforRanking
    param = confSet->Get("ReduceNPforRanking");
    if( param != ""){
        fFitter->fReduceNPforRanking = atof(param.c_str());
    }

    // Set ImageFormat
    param = confSet->Get("ImageFormat");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param,',');
        if (tmp.size() > 0) fFitter->fImageFormat = tmp.at(0);
        else {
            WriteErrorStatus("ConfigReader::ReadJobOptions", "You specified 'ImageFormat' option but we cannot split the setting. Please check");
        }
        TtHFitter::IMAGEFORMAT = tmp;
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
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'StatOnly' option but you didn't provide valid setting. Using default (false)");
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
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'FixNPforStatOnly' option but you didn't provide valid setting. Using default (false)");
            fFitter->fFixNPforStatOnlyFit = false;
        }
    }

    // Set InputFolder
    param = confSet->Get("InputFolder");
    if( param != "" ){
        fFitter->fInputFolder = param;
    }

    // Set InputName
    param = confSet->Get("InputName");
    if( param != "" ){
        fFitter->fInputName = param;
    }

    // Set WorkspaceFileName
    param = confSet->Get("WorkspaceFileName");
    if( param != "" ){
        fFitter->fWorkspaceFileName = param;
    }

    // Set KeepPruning
    param = confSet->Get("KeepPruning");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fKeepPruning = true;
        else if (param == "FALSE") fFitter->fKeepPruning = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'KeepPruning' option but you didn't provide valid setting. Using default (false)");
            fFitter->fKeepPruning = false;
        }
    }


    // Set AtlasLabel
    param = confSet->Get("AtlasLabel");
    if( param != "" ){
        fFitter->fAtlasLabel = param;
    }
        
    // Set CleanTables
    param = confSet->Get("CleanTables");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fCleanTables = true;
        else if (param == "FALSE") fFitter->fCleanTables = false;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'CleanTables' option but you didn't provide valid setting. Using default (false)");
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
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'SystCategoryTables' option but you didn't provide valid setting. Using default (false)");
            fFitter->fSystCategoryTables = false;
        }
    }

    // Set Suffix
    param = confSet->Get("Suffix");
    if( param != "" ){
        fFitter->fSuffix = param;
    }
    
    // Set SaveSuffix
    param = confSet->Get("SaveSuffix");
    if( param != "" ){
        fFitter->fSaveSuffix = param;
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
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'KeepPrefitBlindedBins' option but you didn't provide valid setting. Using default (false)");
            fFitter->fKeepPrefitBlindedBins = false;
        }
    }

    // Set CustomAsimov
    param = confSet->Get("CustomAsimov");
    if( param != "" ){
        fFitter->fCustomAsimov = param;
    }

    // Set RandomPOISeed
    param = confSet->Get("RandomPOISeed");
    if( param != "" ){
        int seed = atoi(param.c_str());
        if(seed>=0){
             fFitter->fRandomPOISeed = seed;
        }
        else {
            WriteErrorStatus("ConfigReader::ReadJobOptions", "You specified 'RandomPOISeed' option but the provided value is < 0. Check this!");
            return 1;
        }
    }

    // Set GetChi2
    param = confSet->Get("GetChi2");
    if( param != "" ){ // can be TRUE, SYST+STAT, STAT-ONLY... (if it contains STAT and no SYST => stat-only, ptherwise stat+syst)
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fGetChi2 = 2;
        }
        if( param.find("SYST")!=std::string::npos ){
            fFitter->fGetChi2 = 2;
        }
        else if( param.find("STAT")!=std::string::npos ){
            fFitter->fGetChi2 = 1;
        } else {
            WriteErrorStatus("ConfigReader::ReadJobOptions", "You specified 'GetChi2' option but you didn;t provide valid option. Check this!");
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
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'DoTables' option but you didn't provide valid setting. Using default (false)");
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
        fFitter->fBootstrap = param;
    }

    // Set RunROOTMacros
    param = confSet->Get("RunROOTMacros");
    if ( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE"){
            fFitter->fRunROOTMacros = true; 
        }
        else if (param == "FALSE"){
            fFitter->fRunROOTMacros = false; 
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified RunROOTMacros option but didnt provide valid parameter. Using default (false)");
            fFitter->fRunROOTMacros = false; 
        }
    }

    // Set DecorrSuff
    param = confSet->Get("DecorrSuff");
    if( param != ""){
        fFitter->fDecorrSuff = param;
    }

    // success
    return 0;
}

int ConfigReader::SetJobPlot(ConfigSet *confSet){

    // Plot option
    std::string param = confSet->Get("PlotOptions");
    std::vector<std::string> vec;
    if( param != ""){
        vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "YIELDS") !=vec.end() )  TtHFitter::SHOWYIELDS     = true;
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TtHFitter::SHOWSTACKSIG   = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TtHFitter::SHOWNORMSIG    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TtHFitter::SHOWOVERLAYSIG = true;
        if( std::find(vec.begin(), vec.end(), "LEFT")   !=vec.end() )  TtHFitter::LEGENDLEFT     = true;
        if( std::find(vec.begin(), vec.end(), "CHI2")   !=vec.end() )  TtHFitter::SHOWCHI2       = true;
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TtHFitter::PREFITONPOSTFIT= true;
        if( std::find(vec.begin(), vec.end(), "POISSONIZE")        !=vec.end() )  TtHFitter::POISSONIZE     = true;
        if( std::find(vec.begin(), vec.end(), "NOXERR") !=vec.end() )  TtHFitter::REMOVEXERRORS  = true;
        if( std::find(vec.begin(), vec.end(), "NOENDERR") !=vec.end() )TtHFitter::NOENDERR       = true;
    }

    // Set PlotOptionsSummary
    param = confSet->Get("PlotOptionsSummary");
    if( param != ""){
        vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TtHFitter::SHOWSTACKSIG_SUMMARY   = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TtHFitter::SHOWNORMSIG_SUMMARY    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TtHFitter::SHOWOVERLAYSIG_SUMMARY = true;
    }
    else{
        WriteDebugStatus("ConfigReader::SetJobPlot", "PlotOptionsSummary not specified setting Summary values to 'PlotOptions'");
        TtHFitter::SHOWSTACKSIG_SUMMARY   = TtHFitter::SHOWSTACKSIG    ;
        TtHFitter::SHOWNORMSIG_SUMMARY    = TtHFitter::SHOWNORMSIG     ;
        TtHFitter::SHOWOVERLAYSIG_SUMMARY = TtHFitter::SHOWOVERLAYSIG  ;
    }
    
    // Set SystControlPlots
    param = confSet->Get("SystControlPlots");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TtHFitter::SYSTCONTROLPLOTS = true;
        } else if (param == "FALSE"){
            TtHFitter::SYSTCONTROLPLOTS = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'SystControlPlots' option but you didn't provide valid setting. Using default (FALSE)");
            TtHFitter::SYSTCONTROLPLOTS = false;
        }
    }

    // Set SystDataPlots
    param = confSet->Get("SystDataPlots");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TtHFitter::SYSTDATAPLOT = true;
            fFitter->fSystDataPlot_upFrame=false;
        } else if( param == "FILLUPFRAME" ){
            TtHFitter::SYSTDATAPLOT = true;
            fFitter->fSystDataPlot_upFrame=true;
        } else {
            WriteDebugStatus("ConfigReader::SetJobPlot", "You specified 'SystDataPlots' option but the value is not 'TRUE' nor 'FILLUPFRAME'. Setting SystDataPlot and SystDataPlot_upFrame to false");
            TtHFitter::SYSTDATAPLOT = false;
            fFitter->fSystDataPlot_upFrame=false;
        }
    }

    // Set SystErrorBars
    param = confSet->Get("SystErrorBars");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE" ){
            TtHFitter::SYSTERRORBARS = true;
        } else if (param == "FALSE"){
            TtHFitter::SYSTERRORBARS = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'SystErrorBars' option but you didn't provide valid setting. Using default (FALSE)");
            TtHFitter::SYSTERRORBARS = false;
        }
    }

    // Set GuessMCStatEmptyBins
    param = confSet->Get("GuessMCStatEmptyBins");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            TtHFitter::GUESSMCSTATERROR = true;
        } else if (param == "FALSE") {
            TtHFitter::GUESSMCSTATERROR = false;
        } else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'GuessMCStatEmptyBins' option but you didn't provide valid setting. Using default (FALSE)");
            TtHFitter::GUESSMCSTATERROR = false;
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
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified SuppressNegativeBinWarnings option but didnt provide valid parameter. Using default (false)");
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

    // Set DoSummaryPlot
    param = confSet->Get("DoSummaryPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoSummaryPlot = true;
        else if( param == "FALSE" ) fFitter->fDoSummaryPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoSummaryPlot' option but didnt provide valid parameter. Using default (false)");
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
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoMergedPlot' option but didnt provide valid parameter. Using default (false)");
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
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoSignalRegionsPlot' option but didnt provide valid parameter. Using default (false)");
            fFitter->fDoSignalRegionsPlot = false;
        }
    }
    param = confSet->Get("DoPieChartPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoPieChartPlot = true;
        else if( param == "FALSE" ) fFitter->fDoPieChartPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoPieChartPlot' option but didnt provide valid parameter. Using default (false)");
            fFitter->fDoPieChartPlot = false;
        }
    }

    // Set RankingPlot
    param = confSet->Get("RankingPlot");
    if( param != ""){
        fFitter->fRankingPlot = param;
    }

    return 0;
}

int ConfigReader::ReadGeneralOptions(){
    ConfigSet* confSet = fParser.GetConfigSet("Options");        
    if (confSet != nullptr){
        for(int i=0; i < confSet->GetN(); i++){
            if(confSet->GetConfigValue(i) != ""){
                TtHFitter::OPTION[confSet->GetConfigName(i)] = atof(confSet->GetConfigValue(i).c_str());
            }
        }
    } else {
        WriteDebugStatus("ConfigReader::ReadGeneralOptions", "You do not have 'Options' option in the config. It is ok, we just want to let you know.");
    }

    return 0;
}

int ConfigReader::ReadFitOptions(){
    std::string param = "";

    ConfigSet *confSet = fParser.GetConfigSet("Fit");
    if (confSet == nullptr){
        WriteInfoStatus("ConfigReader::ReadFitOptions", "You do not have Fit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up 
    }

    //Set FitType
    param = confSet->Get("FitType");
    if( param != "" && fFitter->fFitType == TtHFit::UNDEFINED ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "SPLUSB" ){
            fFitter->SetFitType(TtHFit::SPLUSB);
        }
        else if( param == "BONLY" ){
            fFitter->SetFitType(TtHFit::BONLY);
        }
        else{
            WriteErrorStatus("ConfigReader::ReadFitOptions", "Unknown FitType argument : " + confSet->Get("FitType"));
            return 1;
        }
    }
    else if( fFitter->fFitType == TtHFit::UNDEFINED ){
        WriteInfoStatus("ConfigReader::ReadFitOptions","Setting default fit Type SPLUSB");
        fFitter->SetFitType(TtHFit::SPLUSB);
    }

    // Set FitRegion
    param = confSet->Get("FitRegion");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if( param != "" ){
        if( param == "CRONLY" ){
            fFitter->SetFitRegion(TtHFit::CRONLY);
        }
        else if( param == "CRSR" ){
            fFitter->SetFitRegion(TtHFit::CRSR);
        }
        else{
            fFitter->SetFitRegion(TtHFit::USERSPECIFIC);
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
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'FitBlind' option but didnt provide valid parameter. Using default (false)");
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
        std::vector < std::string > temp_vec = Vectorize(param,',');
        for(std::string iNP : temp_vec){
            std::vector < std::string > np_value = Vectorize(iNP,':');
            if(np_value.size()==2){
                fFitter->fFitNPValues.insert( std::pair < std::string, double >( np_value[0], atof(np_value[1].c_str()) ) );
            } else {
                WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'NPValues' option but didnt provide 2 parameters for each NP which is expected. Ignoring");
            }
        }
    }
    
    // Set FixNPs
    param = confSet->Get("FixNPs");
    if( param != "" ){
        std::vector < std::string > temp_fixedNPs = Vectorize(param,',');
        for(std::string iNP : temp_fixedNPs){
            std::vector < std::string > fixed_nps = Vectorize(iNP,':');
            if(fixed_nps.size()==2){
                fFitter->fFitFixedNPs.insert( std::pair < std::string, double >( fixed_nps[0], atof(fixed_nps[1].c_str()) ) );
            } else {
                WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'FixNPs' option but didnt provide 2 parameters for each NP which is expected. Ignoring");
            }
        }
    }

    // Set doLHscan
    param = confSet->Get("doLHscan");
    if( param != "" ){
        fFitter->fVarNameLH = Vectorize(param,',');
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
        TtHFitter::NCPU = std::atoi( param.c_str());
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
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'StatOnlyFit' option but didnt provide valid parameter. Using default (false)");
            fFitter->fStatOnlyFit = false;
        }
    }

    // Set GetGoodnessOfFit
    param = confSet->Get("GetGoodnessOfFit");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fGetGoodnessOfFit = true;
        } else if (param == "FALSE"){
            fFitter->fGetGoodnessOfFit = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadFitOptions", "You specified 'GetGoodnessOfFit' option but didnt provide valid parameter. Using default (false)");
            fFitter->fGetGoodnessOfFit = false;
        }
    }

    return 0;
}

int ConfigReader::ReadLimitOptions(){
    std::string param = "";

    ConfigSet* confSet = fParser.GetConfigSet("Limit");
    if (confSet == nullptr){
        WriteInfoStatus("ConfigReader::ReadLimitOptions", "You do not have Limit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up 
    }

    // Set LimitType    
    param = confSet->Get("LimitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "ASYMPTOTIC" ){
            fFitter->SetLimitType(TtHFit::ASYMPTOTIC);
        }
        else if( param == "TOYS" ){
            fFitter->SetLimitType(TtHFit::TOYS);
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
            WriteWarningStatus("ConfigReader::ReadLimitOptions", "You specified 'LimitBlind' option but didnt provide valid parameter. Using default (false)");
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
            WriteWarningStatus("ConfigReader::ReadLimitOptions", "You specified 'SignalInjection' option but didnt provide valid parameter. Using default (false)");
            fFitter->fSignalInjection = false;
        }
    }

    return 0; 
}

int ConfigReader::ReadRegionOptions(){
    std::string param = "";
    int nReg = 0;

    Region *reg;

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("Region",nReg);
        if (confSet == nullptr) break;

        nReg++;
        if(fOnlyRegions.size()>0 && FindInStringVector(fOnlyRegions,confSet->GetValue())<0) continue;
        if(fToExclude.size()>0 && FindInStringVector(fToExclude,confSet->GetValue())>=0) continue;
        fRegNames.push_back( CheckName(confSet->GetValue()) ); //why the CheckName is needed?? A: cs->GetValue() might have leading/trailing spaces...
        
        reg = fFitter->NewRegion((confSet->GetValue()));
        reg->fGetChi2 = fFitter->fGetChi2;
        reg->SetVariableTitle(confSet->Get("VariableTitle"));
        reg->SetLabel(confSet->Get("Label"),confSet->Get("ShortLabel"));

        // Set axisTitle
        param = confSet->Get("YaxisTitle");
        if( param != "") reg->fYTitle = param;

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
        if( param != "") reg->fTexLabel = param;

        // Set LumiLabel
        param = confSet->Get("LumiLabel");
        if( param != "") reg->fLumiLabel = param;
        else reg->fLumiLabel = fFitter->fLumiLabel;

        // Set CmeLabel
        param = confSet->Get("CmeLabel");
        if( param != "") reg->fCmeLabel = param;
        else reg->fCmeLabel = fFitter->fCmeLabel;

        // Set LogScale
        param = confSet->Get("LogScale");
        if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") reg->fLogScale = true;
            else if(param=="FALSE") reg->fLogScale = false;
            else {
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "You specified 'LogScale' option but didnt provide valid parameter. Using default (false)");
                reg->fLogScale = false;
            }
        }

        // Set Group
        param = confSet->Get("Group");
        if( param != "") reg->fGroup = param;

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

        // Set Binning
        param = confSet->Get("Binning");
        if(param != "" && param !="-"){
            std::vector < std::string > vec_bins = Vectorize(confSet->Get("Binning"), ',');
            if (vec_bins.size() == 0){
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option, but you didn't provide any reasonable option. Check this!");
                return 1;
            }
            if(vec_bins[0]=="AutoBin"){
                if (vec_bins.size() < 2){
                    WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option with Autobin, but you didn't provide any reasonable option. Check this!");
                    return 1;
                }
                reg -> fBinTransfo = vec_bins[1];
                if(vec_bins[1]=="TransfoD"){
                    if (vec_bins.size() < 4){
                        WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option with TransfoD, but you didn't provide any reasonable option. Check this!");
                        return 1;
                    }
                    reg -> fTransfoDzSig=convertStoD(vec_bins[2]);
                    reg -> fTransfoDzBkg=convertStoD(vec_bins[3]);
	                if(vec_bins.size()>4){
	                    for(const std::string& ibkg : vec_bins){
	    	                reg -> fAutoBinBkgsInSig.push_back(ibkg);
	                    }
	                }
                }
                else if(vec_bins[1]=="TransfoF"){
                    if (vec_bins.size() < 4){
                        WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified `Binning` option with TransfoF, but you didn't provide any reasonable option. Check this!");
                        return 1;
                    }
                    reg -> fTransfoFzSig=convertStoD(vec_bins[2]);
                    reg -> fTransfoFzBkg=convertStoD(vec_bins[3]);
	                if(vec_bins.size()>4){
	                    for(const std::string& ibkg : vec_bins){
	    	                reg -> fAutoBinBkgsInSig.push_back(ibkg);
	                    }
	                }
                }
                else if(vec_bins[1]=="TransfoJ"){
                    if(vec_bins.size() > 2) reg -> fTransfoJpar1=convertStoD(vec_bins[2]);
                    else reg -> fTransfoJpar1 = 5.;
                    if(vec_bins.size() > 3) reg -> fTransfoJpar2=convertStoD(vec_bins[3]);
                    else reg -> fTransfoJpar2 = 1.;
                    if(vec_bins.size() > 4) reg -> fTransfoJpar3=convertStoD(vec_bins[4]);
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
        param == confSet->Get("Type");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="CONTROL" )     reg -> SetRegionType(Region::CONTROL);
            else if( param=="VALIDATION" )  reg -> SetRegionType(Region::VALIDATION);
            else if( param=="SIGNAL" )      reg -> SetRegionType(Region::SIGNAL);
            else {
                WriteErrorStatus("ConfigReader::ReadRegionOptions", "You specified 'Type' option in region but didnt provide valid parameter. Please check this!");
                return 1;
            }
        }
        if(reg -> fRegionType != Region::VALIDATION) reg->fUseGammaPulls = fFitter->fUseGammaPulls;

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
                WriteWarningStatus("ConfigReader::ReadRegionOptions", "You specified 'SkipSmoothing' option in region but didnt provide valid parameter. Using default (FALSE)");
                reg->fSkipSmoothing = false;
            }
        }

        // Set DropBins
        param = confSet->Get("DropBins");
        if( param != "" ){
            reg->fDropBins.clear();
            std::vector<std::string> s = Vectorize( param,',' );
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

    }
    return 0;
}

int ConfigReader::SetRegionHIST(Region* reg, ConfigSet *confSet){
    std::string param = "";

    // Set HistoFile
    param = confSet->Get("HistoFile");
    if(param!="") reg->fHistoFiles.push_back( param );

    // Set HistoName
    param = confSet->Get("HistoName"); 
    if(param!="") reg->SetHistoName( param );

    // Set HistoPathSuff
    param = confSet->Get("HistoPathSuff");
    if(param !="") {
        reg->fHistoPathSuffs.clear();
        reg->fHistoPathSuffs.push_back( Fix(param) );
    }

    // Set HistoPathSuffs
    param = confSet->Get("HistoPathSuffs");
    if( param != "" ){
        reg->fHistoPathSuffs.clear();
        std::vector<std::string> paths = Vectorize( param,',' );
        for(std::string ipath : paths){
            reg->fHistoPathSuffs.push_back( Fix(ipath) );
        }
    }
   
    // Check for NTUP inputs
    if (ConfigHasNTUP(confSet)){
        WriteWarningStatus("ConfigReader::SetRegionHIST", "Found some NTUP settings in Region, while input option is HIST. Ignoring them.");
    }
    
    return 0;
}

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
        std::vector < std::string > temp_samplesAndVars = Vectorize(param,',');
        for(std::string ivar : temp_samplesAndVars){
          std::vector < std::string > vars = Vectorize(ivar,':');
            if(vars.size()==2){
                reg->SetAlternativeVariable(vars[1], vars[0]);
            }
        }
    }

    // Set Selection
    param = confSet->Get("Selection");
    if(param != "")
    reg->AddSelection( param );

    // Set NtupleName
    param = confSet->Get("NtupleName");
    if(param!="") { 
        reg->fNtupleNames.clear();
        reg->fNtupleNames.push_back(param);
    }

    // Set NtupleNameSuff
    param = confSet->Get("NtupleNameSuff");
    if(param!="") {
        reg->fNtupleNameSuffs.clear();
        reg->fNtupleNameSuffs.push_back( param );
    }

    // Set NtupleNameSuffs
    param = confSet->Get("NtupleNameSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Vectorize( param,',' );
        reg->fNtupleNameSuffs = paths;
    }
    
    // Set MCweight
    param = confSet->Get("MCweight");
    if (param != "") reg->fMCweight = param; // this will override the global MCweight, if any

    // Set NtuplePathSuff
    param = confSet->Get("NtuplePathSuff");
    if(param != "") {
        reg->fNtuplePathSuffs.clear();
        reg->fNtuplePathSuffs.push_back( param );
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

int ConfigReader::ReadSampleOptions(){
    int nSmp = 0;
    Sample *sample = nullptr;
    NormFactor *nfactor = nullptr;
    ShapeFactor *sfactor = nullptr;

    int type = 0;
    std::string param = "";

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("Sample",nSmp);
        if (confSet == nullptr) break;
        nSmp++;
        
        if(fOnlySamples.size()>0 && FindInStringVector(fOnlySamples,confSet->GetValue())<0) continue;
        if(fToExclude.size()>0 && FindInStringVector(fToExclude,confSet->GetValue())>=0) continue;
        type = Sample::BACKGROUND;

        // Set Type
        param = confSet->Get("Type");
        if (param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "SIGNAL") type = Sample::SIGNAL;
            else if(param == "DATA") type = Sample::DATA;
            else if(param == "GHOST") type = Sample::GHOST;
            else if(param == "BACKGROUND") type = Sample::BACKGROUND;
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'Type' option in sample but didnt provide valid parameter. Using default (BACKGROUND)");
            }
            if(fOnlySignal != "" && type==Sample::SIGNAL && confSet->GetValue()!=fOnlySignal) continue;
        } 
        sample = fFitter->NewSample((confSet->GetValue()),type);
    
        // Set Title
        param = confSet->Get("Title");
        if (param != "") sample->SetTitle(param);

        // Set TexTitle
        param = confSet->Get("TexTitle");
        if(param!="") sample->fTexTitle = param;

        // Set Group
        param = confSet->Get("Group");
        if(param!="") sample->fGroup = param;

        // HIST input
        if (fFitter->fInputType == 0){
            // Set HistoFile
            param = confSet->Get("HistoFile");
            if(param!="") sample->AddHistoFile( param );

            // Set HistoName
            param = confSet->Get("HistoName");
            if(param!="") sample->fHistoNames.push_back( param );

            // Set HistoPath
            param = confSet->Get("HistoPath");
            if(param!="") sample->AddHistoPath( param );
        
            if (ConfigHasNTUP(confSet)){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You provided some NTUP options but your input type is HIST. Options will be ingored");
            }
        } else if (fFitter->fInputType == 1){ // NTUP input
            // Set NtupleFile
            param = confSet->Get("NtupleFile");
            if(param!="") sample->AddNtupleFile( param );

            // Set NtupleFiles
            param = confSet->Get("NtupleFiles");
            if(param!="") sample->fNtupleFiles = Vectorize( param ,',' );

            // Set NtupleName
            param = confSet->Get("NtupleName");
            if(param!="") sample->AddNtupleName( param );

            // Set NtupleNames
            param = confSet->Get("NtupleNames");
            if(param!="") sample->fNtupleNames = Vectorize( param ,',' );
            
            // Set NtuplePath
            param = confSet->Get("NtuplePath");
            if(param!="") sample->AddNtuplePath( param );

            // Set NtuplePaths
            param = confSet->Get("NtuplePaths");
            if(param != "") sample->fNtuplePaths = Vectorize( param ,',' );

            // Set NtupleNameSuff
            param = confSet->Get("NtupleNameSuff");
            if(param != "") {
                sample->fNtupleNameSuffs.clear();
                sample->fNtupleNameSuffs.push_back( param );
            }
            
            // Set NtupleNameSuffs
            param = confSet->Get("NtupleNameSuffs");
            if( param != "" ){
                std::vector<std::string> paths = Vectorize( param,',' );
                sample->fNtupleNameSuffs = paths;
            }

            if (ConfigHasHIST(confSet)){
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You provided some HIST options but your input type is NTUP. Options will be ingored");
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
            if( FindInStringVector(fFitter->fNormFactorNames,nfactor->fName)<0 ){
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
            if( FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
                fFitter->fShapeFactors.push_back( sfactor );
                fFitter->fShapeFactorNames.push_back( sfactor->fName );
                fFitter->fNShape++;
            }
        }

        // Set NormalizedByTheory
        param = confSet->Get("NormalizedByTheory");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") sample->NormalizedByTheory(false);
            else if(param=="TRUE") sample->NormalizedByTheory(true);
            else{
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'NormalizedByTheory' option but didnt provide valid parameter. Using default (true)");
                sample->NormalizedByTheory(true);
            }
        }

        // Set MCweight and Selection
        if(fFitter->fInputType==1){
            param = confSet->Get("MCweight");
            if(param != "") sample->SetMCweight( param );

            param = confSet->Get("Selection");
            if(param!="") sample->SetSelection( param );
        }


        // to specify only certain regions
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');
        sample->fRegions.clear();

        for(int i_reg=0;i_reg<fFitter->fNRegions;i_reg++){
            std::string regName = fFitter->fRegions[i_reg]->fName;
            if( (regions_str=="" || regions_str=="all" || FindInStringVector(regions,regName)>=0)
                && FindInStringVector(exclude,regName)<0 ){
                sample->fRegions.push_back( fFitter->fRegions[i_reg]->fName );
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
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'IgnoreSelection' option but didnt provide valid parameter. Using default (false)");
                sample->fIgnoreSelection = "FALSE";
            }
        }
        
        // Set UseMCstat
        // to skip MC stat uncertainty for this sample
        param = confSet->Get("UseMCstat");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") sample->fUseMCStat = false;
            else if (param == "TRUE") sample->fUseMCStat = true;
            else {
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'UseMCstat' option but didnt provide valid parameter. Using default (true)");
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
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'UseSystematics' option but didnt provide valid parameter. Using default (true)");
                sample->fUseSystematics = true;
            }
        }

        // Set DivideBy
        param = confSet->Get("DivideBy");;
        if (param != "") sample->fDivideBy = param;

        // Set MultiplyBy
        param = confSet->Get("MultiplyBy");
        sample->fMultiplyBy = param;

        // Set SubtractSample
        param = confSet->Get("SubtractSample");
        if(param!="") sample->fSubtractSamples.push_back(param);

        // Set SubtractSamples
        param = confSet->Get("SubtractSamples");
        if(param != "") sample->fSubtractSamples = Vectorize(param,',');

        // Set AddSample
        param = confSet->Get("AddSample");
        if(param != "") sample->fAddSamples.push_back(param);

        // Set AddSamples
        param = confSet->Get("AddSamples");
        if(param!="") sample->fAddSamples = Vectorize(param,',');
        
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
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'Smooth' option but didnt provide valid parameter. Using default (false)");
                sample->fSmooth = false;
            }
        }

        // Set AsimovReplacementFor
        // AsimovReplacementFor
        param = confSet->Get("AsimovReplacementFor");
        if(param != ""){
            if(Vectorize(param,',').size() == 2){
                sample->fAsimovReplacementFor.first  = Vectorize(param,',')[0];
                sample->fAsimovReplacementFor.second = Vectorize(param,',')[1];
            } else {
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "You specified 'AsimovReplacementFor' option but didnt provide 2 parameters. Please check this");
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
                WriteWarningStatus("ConfigReader::ReadSampleOptions", "You specified 'SeparateGammas' option but didnt provide valid parameter. Using default (false)");
                sample->fSeparateGammas = false;
            }
        }

        // Set CorrelateGammasInRegions
        // in the form    CorrelateGammasInRegions: SR1:SR2,CR1:CR2:CR3
        param = confSet->Get("CorrelateGammasInRegions");
        if(param != ""){
            std::vector<std::string> sets = Vectorize(param,',');
            for(std::string set : sets){
                std::vector<std::string> regions = Vectorize(set,':');
                WriteDebugStatus("ConfigReader::ReadSampleOptions", "Correlating gammas for this sample in regions " + set);
                sample->fCorrelateGammasInRegions.push_back(regions);
            }
        }

        // Set Morphing
        param = confSet->Get("Morphing");
        if(param != ""){
            std::vector<std::string> morph_par = Vectorize(param,',');
            if (morph_par.size() != 2){
                WriteErrorStatus("ConfigReader::ReadSampleOptions", "Morphing requires exactly 2 parameters, but " + std::to_string(morph_par.size()) + " provided");
                return 1;
            }
            fFitter->fRunMorphing = true;
            std::string name      = morph_par.at(0);
            float value = std::stof(morph_par.at(1));
            WriteDebugStatus("ConfigReader::ReadSampleOptions", "Morphing: Adding " + name + ", with value: " + std::to_string(value));
            if (!fFitter->MorphIsAlreadyPresent(name, value)) fFitter->AddTemplateWeight(name, value);
            // set proper normalization
            std::string morphName = "morph_"+name+"_"+ReplaceString(std::to_string(value),"-","m");
            NormFactor *nf = sample->AddNormFactor(morphName, 1, 0, 10, false);
            fFitter->fNormFactors.push_back( nf );
            fFitter->fNormFactorNames.push_back( nf->fName );
            fFitter->fNNorm++;

            sample->fIsMorph = true;
        }
        
    }
    
    // build new samples if AsimovReplacementFor are specified
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        if(fFitter->fSamples[i_smp]->fAsimovReplacementFor.first!=""){
            WriteDebugStatus("ConfigReader::ReadSampleOptions", "Creating sample " + fFitter->fSamples[i_smp]->fAsimovReplacementFor.first);
            Sample *ca = fFitter->NewSample("customAsimov_"+fFitter->fSamples[i_smp]->fAsimovReplacementFor.first,Sample::GHOST);
            ca->SetTitle("Pseudo-Data ("+fFitter->fSamples[i_smp]->fAsimovReplacementFor.first+")");
            ca->fUseSystematics = false;
        }
    }

    if (nSmp == 0){
        WriteErrorStatus("ConfigReader::ReadSampleOptions", "No 'Sample' provided. You need to provide at least one 'Sample' object. Check this!");
        return 1;
    }

    return 0;
}

int ConfigReader::ReadNormFactorOptions(){
    std::string param = "";

    int nNorm = 0;
    NormFactor *nfactor = nullptr;
    Sample *sample = nullptr;

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("NormFactor", nNorm);
        if (confSet == nullptr) break;
        nNorm++;
        
        if(fToExclude.size()>0 && FindInStringVector(fToExclude,confSet->GetValue())>=0) continue;

        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Vectorize(samples_str,',');
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');

        nfactor = new NormFactor(CheckName(confSet->GetValue()));
        
        TtHFitter::SYSTMAP[nfactor->fName] = nfactor->fName;
        if( FindInStringVector(fFitter->fNormFactorNames,nfactor->fName)<0 ){
            fFitter->fNormFactors.push_back( nfactor );
            fFitter->fNormFactorNames.push_back( nfactor->fName );
            fFitter->fNNorm++;
        }
        else{
            nfactor = fFitter->fNormFactors[ FindInStringVector(fFitter->fNormFactorNames,nfactor->fName) ];
        }

        // Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            nfactor->fNuisanceParameter = param;
            TtHFitter::NPMAP[nfactor->fName] = nfactor->fNuisanceParameter;
        }
        else{
            nfactor->fNuisanceParameter = nfactor->fName;
            TtHFitter::NPMAP[nfactor->fName] = nfactor->fName;
        }

        // Set Constant
        param = confSet->Get("Constant");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") nfactor->fConst = true;
            else if(param=="FALSE") nfactor->fConst = false;
            else {
                WriteWarningStatus("ConfigReader::ReadNormFactorOptions", "You specified 'Constant' option but didnt provide valid parameter. Using default (false)");
                nfactor->fConst = false;
            }
        }

        // Set Category
        param = confSet->Get("Category");
        if(param != "") nfactor->fCategory = param;

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            nfactor->fTitle = param;
            TtHFitter::SYSTMAP[nfactor->fName] = nfactor->fTitle;
        }
        
        // Set TexTitle
        param = confSet->Get("TexTitle");
        if(param != "") TtHFitter::SYSTTEX[nfactor->fName] = param;
        
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
            std::vector<std::string> v = Vectorize(param,',');
            if (v.size() != 2){
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "You specified 'Expression' option but didnt provide 2 parameters. Please check this");
                return 1;
            }
            nfactor->fExpression = std::make_pair(v[0],v[1]);
            // title will contain the expression FIXME
            nfactor->fTitle = v[0];
            TtHFitter::SYSTMAP[nfactor->fName] = v[0];
            // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
            nfactor->fNuisanceParameter = v[1];
            TtHFitter::NPMAP[nfactor->fName] = v[1];
            // set nominal, min and max according to the norm factor the expression depends on FIXME
            for(NormFactor *nf : fFitter->fNormFactors){
                if(nf->fNuisanceParameter == v[1]){
                    nfactor->fNominal = nf->fNominal;
                    nfactor->fMin = nf->fMin;
                    nfactor->fMax = nf->fMax;
                }
            }
        }
        
        // save list of
        if (regions.size() == 0 || exclude.size() == 0){
                WriteErrorStatus("ConfigReader::ReadNormFactorOptions", "Region or excude region size is equal to zero. Please check this");
                return 1;
            
        }
        if(regions[0] != "all") nfactor->fRegions = regions;
        if(exclude[0] != "")    nfactor->fExclude = exclude;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sample = fFitter->fSamples[i_smp];
            if(sample->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || FindInStringVector(samples, sample->fName)>=0 )
               && (exclude[0]==""    || FindInStringVector(exclude, sample->fName)<0 ) ){
                sample->AddNormFactor(nfactor);
            }
        }
    }
    return 0;
}

int ConfigReader::ReadShapeFactorOptions(){
    std::string param = "";
    int nShape = 0;
    ShapeFactor *sfactor = nullptr;
    Sample *sample = nullptr;

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("ShapeFactor",nShape);
        if (confSet == nullptr) break;
        nShape++;

        if(fToExclude.size()>0 && FindInStringVector(fToExclude,confSet->GetValue())>=0) continue;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Vectorize(samples_str,',');
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');
        sfactor = new ShapeFactor(CheckName(confSet->GetValue()));
        if( FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
            fFitter->fShapeFactors.push_back( sfactor );
            fFitter->fShapeFactorNames.push_back( sfactor->fName );
            fFitter->fNShape++;
        }
        else{
            sfactor = fFitter->fShapeFactors[ FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName) ];
        }

        // Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            sfactor->fNuisanceParameter = param;
            TtHFitter::NPMAP[sfactor->fName] = sfactor->fNuisanceParameter;
        }
        else{
            sfactor->fNuisanceParameter = sfactor->fName;
            TtHFitter::NPMAP[sfactor->fName] = sfactor->fName;
        }

        // Set Constant
        param = confSet->Get("Constant");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") sfactor->fConst = true;
            else if(param=="FALSE") sfactor->fConst = false;
            else {
                WriteWarningStatus("ConfigReader::ReadShapeFactorOptions", "You specified 'Constant' option but didnt provide valid parameter. Using default (false)");
                sfactor->fConst = false;
            }
        }

        // Set Category
        param = confSet->Get("Category");
        if(param!="") sfactor->fCategory = param;
    
        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            sfactor->fTitle = param;
            TtHFitter::SYSTMAP[sfactor->fName] = sfactor->fTitle;
        }
        param = confSet->Get("TexTitle");
        if(param != "") TtHFitter::SYSTTEX[sfactor->fName] = param;
        
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
            WriteErrorStatus("ConfigReader::ReadShapeFactorOptions", "Region or excude region size is equal to zero. Please check this");
            return 1;
        }
        // save list of
        if(regions[0]!="all") sfactor->fRegions = regions;
        if(exclude[0]!="")    sfactor->fExclude = exclude;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sample = fFitter->fSamples[i_smp];
            if(sample->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || FindInStringVector(samples, sample->fName)>=0 )
               && (exclude[0]==""    || FindInStringVector(exclude, sample->fName)<0 ) ){
                sample->AddShapeFactor(sfactor);
            }
        }
    }

    return 0;
}

int ConfigReader::ReadSystOptions(){
    std::string param = "";
    int nSys = 0;
    Systematic *sys = nullptr;
    Sample *sample = nullptr;
    int type = 0;

    int typed = 0;
    //Addition for StatOnly fit: dummy systematic for the significance computation and limit setting
    Systematic *sysd;
    if (fFitter->fStatOnly) {
        typed = Systematic::OVERALL;
        sysd = new Systematic("Dummy",typed);
        sysd->fOverallUp   = 0.;
        sysd->fOverallDown = -0.;
        sysd->fScaleUp   = 1.;
        sysd->fScaleDown   = 1.;
        fFitter->fSystematics.push_back( sysd );
        TtHFitter::SYSTMAP[sysd->fName] = "Dummy";
        fFitter->fNSyst++;
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sample = fFitter->fSamples[i_smp];
            if(sample->fType == Sample::SIGNAL ) {
                sample->AddSystematic(sysd);
            }
        }
    }

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("Systematic",nSys);
        if (confSet == nullptr) break;
        nSys++;
        
        if(fOnlySystematics.size()>0 && FindInStringVector(fOnlySystematics,confSet->GetValue())<0) continue;
        if(fToExclude.size()>0 && FindInStringVector(fToExclude,confSet->GetValue())>=0) continue;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        std::string excludeRegionSample_str = confSet->Get("ExcludeRegionSample");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Vectorize(samples_str,',');
        std::vector<std::string> regions = Vectorize(regions_str,',');
        std::vector<std::string> exclude = Vectorize(exclude_str,',');
        fExcludeRegionSample = Vectorize(excludeRegionSample_str,',');
        type = Systematic::HISTO;

        // Set type
        param = confSet->Get("Type");
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);

        if(param == "OVERALL") type = Systematic::OVERALL;
        else if(param == "SHAPE") type = Systematic::SHAPE;
        else if (param == "STAT") type = Systematic::STAT;

        std::string decorrelate = confSet->Get("Decorrelate");
        
        sys = new Systematic(CheckName(confSet->GetValue()),type);
        TtHFitter::SYSTMAP[sys->fName] = sys->fTitle;
        if(param == "OVERALL") sys->fIsNormOnly=true;
        
        // SetCategory
        param = confSet->Get("Category");
        if(param != "") sys->fCategory = param;

        // Set IsFreeParameter
        // Experimental
        param = confSet->Get("IsFreeParameter");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param != "TRUE") sys->fIsFreeParameter = true;
            else if (param != "TRUE") sys->fIsFreeParameter = false;
            else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'IsFreeParameter' option but didnt provide valid parameter. Using default (false)");
                sys->fIsFreeParameter = false;
            }
        }

        // Set StoredName
        // New: name to use when writing / reading the Histograms file
        param = confSet->Get("StoredName");
        if(param != "") sys->fStoredName = param;
        
        bool hasUp   = false;
        bool hasDown = false;
        if(type==Systematic::HISTO || type==Systematic::SHAPE){
            if(fFitter->fInputType==0){ // HIST input
                if(confSet->Get("HistoPathUp")!=""){
                    sys->fHistoPathsUp.push_back(confSet->Get("HistoPathUp"));
                    hasUp   = true;
                }
                if(confSet->Get("HistoPathDown")!=""){
                    sys->fHistoPathsDown.push_back(confSet->Get("HistoPathDown"));
                    hasDown = true;
                }
                if(confSet->Get("HistoPathSufUp")!=""){
                    sys->fHistoPathSufUp = confSet->Get("HistoPathSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("HistoPathSufDown")!=""){
                    sys->fHistoPathSufDown = confSet->Get("HistoPathSufDown");
                    hasDown = true;
                }
                if(confSet->Get("HistoFileUp")!=""){
                    sys->fHistoFilesUp.push_back(confSet->Get("HistoFileUp"));
                    hasUp   = true;
                }
                if(confSet->Get("HistoFileDown")!=""){
                    sys->fHistoFilesDown.push_back(confSet->Get("HistoFileDown"));
                    hasDown = true;
                }
                if(confSet->Get("HistoFileSufUp")!=""){
                    sys->fHistoFileSufUp = confSet->Get("HistoFileSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("HistoFileSufDown")!=""){
                    sys->fHistoFileSufDown = confSet->Get("HistoFileSufDown");
                    hasDown = true;
                }
                if(confSet->Get("HistoNameUp")!=""){
                    sys->fHistoNamesUp.push_back(confSet->Get("HistoNameUp"));
                    hasUp   = true;
                }
                if(confSet->Get("HistoNameDown")!=""){
                    sys->fHistoNamesDown.push_back(confSet->Get("HistoNameDown"));
                    hasDown = true;
                }
                if(confSet->Get("HistoNameSufUp")!=""){
                    sys->fHistoNameSufUp = confSet->Get("HistoNameSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("HistoNameSufDown")!=""){
                    sys->fHistoNameSufDown = confSet->Get("HistoNameSufDown");
                    hasDown = true;
                }
            }
            else if(fFitter->fInputType==1){ // NTUP option
                if(confSet->Get("NtuplePathUp")!=""){
                    sys->fNtuplePathsUp.push_back(confSet->Get("NtuplePathUp"));
                    hasUp   = true;
                }
                if(confSet->Get("NtuplePathDown")!=""){
                    sys->fNtuplePathsDown.push_back(confSet->Get("NtuplePathDown"));
                    hasDown = true;
                }
                if(confSet->Get("NtuplePathsUp")!=""){
                    sys->fNtuplePathsUp = Vectorize(confSet->Get("NtuplePathsUp"),',');
                    hasUp   = true;
                }
                if(confSet->Get("NtuplePathsDown")!=""){
                    sys->fNtuplePathsDown = Vectorize(confSet->Get("NtuplePathsDown"),',');
                    hasDown = true;
                }
                if(confSet->Get("NtuplePathSufUp")!=""){
                    sys->fNtuplePathSufUp = confSet->Get("NtuplePathSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("NtuplePathSufDown")!=""){
                    sys->fNtuplePathSufDown = confSet->Get("NtuplePathSufDown");
                    hasDown = true;
                }
                if(confSet->Get("NtupleFileUp")!=""){
                    sys->fNtupleFilesUp.push_back(confSet->Get("NtupleFileUp"));
                    hasUp   = true;
                }
                if(confSet->Get("NtupleFileDown")!=""){
                    sys->fNtupleFilesDown .push_back(confSet->Get("NtupleFileDown"));
                    hasDown = true;
                }
                if(confSet->Get("NtupleFilesUp")!=""){
                    sys->fNtupleFilesUp = Vectorize(confSet->Get("NtupleFilesUp"),',');
                    hasUp   = true;
                }
                if(confSet->Get("NtupleFilesDown")!=""){
                    sys->fNtupleFilesDown = Vectorize(confSet->Get("NtupleFilesDown"),',');
                    hasDown = true;
                }
                if(confSet->Get("NtupleFileSufUp")!=""){
                    sys->fNtupleFileSufUp = confSet->Get("NtupleFileSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("NtupleFileSufDown")!=""){
                    sys->fNtupleFileSufDown = confSet->Get("NtupleFileSufDown");
                    hasDown = true;
                }
                if(confSet->Get("NtupleNameUp")!=""){
                    sys->fNtupleNamesUp.push_back(confSet->Get("NtupleNameUp"));
                    hasUp   = true;
                }
                if(confSet->Get("NtupleNameDown")!=""){
                    sys->fNtupleNamesDown.push_back( confSet->Get("NtupleNameDown"));
                    hasDown = true;
                }
                if(confSet->Get("NtupleNamesUp")!=""){
                    sys->fNtupleNamesUp = Vectorize(confSet->Get("NtupleNamesUp"),',');
                    hasUp   = true;
                }
                if(confSet->Get("NtupleNamesDown")!=""){
                    sys->fNtupleNamesDown = Vectorize(confSet->Get("NtupleNamesDown"),',');
                    hasDown = true;
                }
                if(confSet->Get("NtupleNameSufUp")!=""){
                    sys->fNtupleNameSufUp = confSet->Get("NtupleNameSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("NtupleNameSufDown")!=""){
                    sys->fNtupleNameSufDown = confSet->Get("NtupleNameSufDown");
                    hasDown = true;
                }
                if(confSet->Get("WeightUp")!=""){
                    sys->fWeightUp = confSet->Get("WeightUp");
                    hasUp   = true;
                }
                if(confSet->Get("WeightDown")!=""){
                    sys->fWeightDown = confSet->Get("WeightDown");
                    hasDown = true;
                }
                if(confSet->Get("WeightSufUp")!=""){
                    sys->fWeightSufUp = confSet->Get("WeightSufUp");
                    hasUp   = true;
                }
                if(confSet->Get("WeightSufDown")!=""){
                    sys->fWeightSufDown = confSet->Get("WeightSufDown");
                    hasDown = true;
                }
                if(confSet->Get("IgnoreWeight")!="") sys->fIgnoreWeight = confSet->Get("IgnoreWeight");
            }
            sys->fHasUpVariation   = hasUp  ;
            sys->fHasDownVariation = hasDown;

            // Set Symmetrisation
            param = confSet->Get("Symmetrisation");
            if(param != ""){
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if(param == "ONESIDED") sys->fSymmetrisationType = HistoTools::SYMMETRIZEONESIDED;
                else if(param == "TWOSIDED") sys->fSymmetrisationType = HistoTools::SYMMETRIZETWOSIDED;
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
                    WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'PreSmoothing' option but didnt provide valid parameter. Using default (false)");
                    sys->fPreSmoothing = false;
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
            std::vector < std::string > temp_vec = Vectorize(param,',');
            if(temp_vec.size()==1 && Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleUp = atof(param.c_str());
            }
            else{
                for(std::string ireg : temp_vec){
                    std::vector < std::string > reg_value = Vectorize(ireg,':');
                    if(reg_value.size()==2){
                        sys->fScaleUpRegions.insert( std::pair < std::string, float >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }

        // Set ScaleDown
        param = confSet->Get("ScaleDown");
        if(param!=""){
            std::vector < std::string > temp_vec = Vectorize(param,',');
            if(temp_vec.size()==1 && Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleDown = atof(param.c_str());
            }
            else{
                for(std::string ireg : temp_vec){
                    std::vector < std::string > reg_value = Vectorize(ireg,':');
                    if(reg_value.size()==2){
                        sys->fScaleDownRegions.insert( std::pair < std::string, float >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }

        // Set SampleUp
        // a new way of defining systematics, just comparing directly with GHOST samples
        // --> this can be used only if this systematic is applied to a single sample
        param = confSet->Get("SampleUp");
        if(param!="") sys->fSampleUp = param;
        
        // Set SampleDown
        param = confSet->Get("SampleDown");
        if(param!="") sys->fSampleDown = param;
        
        // Set ReferenceSample    
        // this to obtain syst variation relatively to given sample
        param = confSet->Get("ReferenceSample");
        if(param!="") sys->fReferenceSample = param;

        // Set KeepReferenceOverallVar
        param = confSet->Get("KeepReferenceOverallVar");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") sys->fKeepReferenceOverallVar = false;
            else if(param == "TRUE" ) sys->fKeepReferenceOverallVar = true;
            else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'KeepReferenceOverallVar' option but didnt provide valid parameter. Using default (true)");
                sys->fKeepReferenceOverallVar = true;
            }
        }

        // Set DropShapeIn
        param = confSet->Get("DropShapeIn");
        if(param!="") sys->fDropShapeIn = Vectorize(param,',');
        
        // Set DropNorm
        param = confSet->Get("DropNorm");
        if(param!="") sys->fDropNormIn = Vectorize(param,',');

        // Set KeepNormForSamples
        param = confSet->Get("KeepNormForSamples");
        if(param!="") sys->fKeepNormForSamples = Vectorize(param,',');
       
        if (regions.size() == 0 || exclude.size() == 0){
            WriteErrorStatus("ConfigReader::ReadSystOptions", "Region or excude region size is equal to zero. Please check this");
            return 1;
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

        // Set SubtractRefSampleVar
        // New: for systeamtics which also vary Data (e.g. JER with Full NPs)
        // This will subtract linearly the relative variation on Data from each relative variation on MC
        param = confSet->Get("SubtractRefSampleVar");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE" ) sys->fSubtractRefSampleVar = true;// default is false
            else if(param == "FALSE" ) sys->fSubtractRefSampleVar = false;// default is false
            else {
                WriteWarningStatus("ConfigReader::ReadSystOptions", "You specified 'PreSmoothing' option but didnt provide valid parameter. Using default (false)");
                sys->fSubtractRefSampleVar = false;
            }
        }
        
        if(FindInStringVector(fFitter->fDecorrSysts,sys->fNuisanceParameter)>=0){
            WriteInfoStatus("ConfigReader::ReadSystOptions","Decorrelating systematic with NP = " + sys->fNuisanceParameter);
            sys->fNuisanceParameter += fFitter->fDecorrSuff;
        }
    
    }
    
    return 0;
}

int ConfigReader::SetSystNoDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude){
    Sample *sam = nullptr;

    fFitter->fSystematics.push_back( sys );
    fFitter->fNSyst++;

    // Set NuisanceParameter
    std::string param = confSet->Get("NuisanceParameter");
    if(param!=""){
        sys->fNuisanceParameter = param;
        TtHFitter::NPMAP[sys->fName] = sys->fNuisanceParameter;
    }
    else{
        sys->fNuisanceParameter = sys->fName;
        TtHFitter::NPMAP[sys->fName] = sys->fName;
    }

    // Set Title
    param = confSet->Get("Title");
    if(param != ""){
        sys->fTitle = param;
        TtHFitter::SYSTMAP[sys->fName] = sys->fTitle;
    }

    // Set TexTitle
    param = confSet->Get("TexTitle");
    if(param!="") TtHFitter::SYSTTEX[sys->fName] = param;
    
    // attach the syst to the proper samples
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        sam = fFitter->fSamples[i_smp];
        if(sam->fType == Sample::DATA) continue;
        // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
        if(sam->fType == Sample::DATA){
            if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
              sam->AddSystematic(sys);
            }
            else continue;
        }
        if(!sam->fUseSystematics) continue;
        if((samples[0] == "all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
           && (exclude[0] == "" || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
        ){
            if((samples[0]=="all" || FindInStringVector(samples, sam->fName)>=0 )
               && (exclude[0]=="" || FindInStringVector(exclude, sam->fName)<0 ) ){
                sam->AddSystematic(sys);
            }
        }
    }
    
    return 0;
}

int ConfigReader::SetSystRegionDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude, const std::vector<std::string> regions, int type){
    Sample *sam = nullptr;
    std::string param = "";
    
    for(std::string ireg : fRegNames) {
        bool keepReg=false;
        if ( regions[0]=="all" ) keepReg=true;
        else {
            for ( const std::string& iGoodReg : regions) {
                if ( ireg == iGoodReg ) keepReg=true;
            }
        }
        for ( const std::string iBadReg : exclude) {
            if ( iBadReg == ireg ) keepReg=false;
        }

        if (!keepReg) {
            WriteInfoStatus("ConfigReader::SetSystRegionDecorelate", "IGNORING REGION: " + ireg);
            continue;
        }
        WriteInfoStatus("ConfigReader::SetSystRegionDecorelate", "--> KEEPING IT!!! " + ireg);

        if (type == Systematic::STAT) {
          Region* reg = fFitter->GetRegion(ireg);
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
                    TtHFitter::NPMAP[mySys->fName] = sys->fNuisanceParameter;
                } else {
                    mySys->fNuisanceParameter = mySys->fName;
                    TtHFitter::NPMAP[mySys->fName] = mySys->fName;
                }
        
                // Set Title
                param = confSet->Get("Title");
                if(param != ""){
                    mySys->fTitle = (sys->fTitle)+"_"+ireg+"_bin"+std::to_string(i_bin);
                    TtHFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
                }
                fFitter->fNSyst++;
                for (int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
                    sam = fFitter->fSamples[i_smp];
                    if(sam->fType == Sample::DATA) continue;
                    if(!sam->fUseSystematics) continue;
                    if(   (samples[0]=="all" || FindInStringVector(samples, sam->fName)>=0 )
                       && (exclude[0]==""    || FindInStringVector(exclude, sam->fName)<0 ) ){
                        sam->AddSystematic(mySys);
                    }
                }
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
            param == confSet->Get("NuisanceParameter");
            if(param != ""){
                mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+ireg;
                TtHFitter::NPMAP[mySys->fName] = sys->fNuisanceParameter;
            }
            else{
                mySys->fNuisanceParameter = mySys->fName;
                TtHFitter::NPMAP[mySys->fName] = mySys->fName;
            }

            // Set Title
            param = confSet->Get("Title");
            if(param != ""){
                mySys->fTitle = (sys->fTitle)+"_"+ireg;
                TtHFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
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
                if(   (samples[0]=="all" || FindInStringVector(samples, sam->fName)>=0 )
                    && (exclude[0]==""    || FindInStringVector(exclude, sam->fName)<0 ) ){
                    sam->AddSystematic(mySys);
                }
            }
        }
    }
    delete sys;

    return 0;
}

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
        // cloning the sys for each region
        Systematic* mySys= new Systematic(*sys);
        mySys->fName=(mySys->fName)+"_"+sam->fName;
        fFitter->fSystematics.push_back( mySys );

        // Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+sam->fName;
            TtHFitter::NPMAP[mySys->fName] = sys->fNuisanceParameter;
        }
        else{
            mySys->fNuisanceParameter = mySys->fName;
            TtHFitter::NPMAP[mySys->fName] = mySys->fName;
        }

        // Set Title
        param == confSet->Get("Title");
        if(param != ""){
            mySys->fTitle = (sys->fTitle)+"_"+sam->fName;
            TtHFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
        }
        fFitter->fNSyst++;
        sam->AddSystematic(mySys);
    }
    delete sys;
    
    return 0;
}

int ConfigReader::SetSystShapeDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string> &samples, const std::vector<std::string> &exclude){
    Sample *sam = nullptr;
    std::string param = "";

    // cloning the sys
    Systematic* mySys1= new Systematic(*sys);
    mySys1->fName=(mySys1->fName)+"_Acc";
    fFitter->fSystematics.push_back( mySys1 );
    mySys1->fIsNormOnly = true;

    // Set NuisanceParameter
    param = confSet->Get("NuisanceParameter");
    if(param != ""){
        mySys1->fNuisanceParameter = (sys->fNuisanceParameter)+"_Acc";
        TtHFitter::NPMAP[mySys1->fName] = sys->fNuisanceParameter;
    }
    else{
        mySys1->fNuisanceParameter = mySys1->fName;
        TtHFitter::NPMAP[mySys1->fName] = mySys1->fName;
    }

    // Set Title
    param = confSet->Get("Title");
    if(param != ""){
        mySys1->fTitle = (sys->fTitle)+"_Acc";
        TtHFitter::SYSTMAP[mySys1->fName] = mySys1->fTitle;
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
        if(   (samples[0]=="all" || FindInStringVector(samples, sam->fName)>=0 )
           && (exclude[0]==""    || FindInStringVector(exclude, sam->fName)<0 ) ){
            sam->AddSystematic(mySys1);
        }
    }

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
            TtHFitter::NPMAP[mySys2->fName] = sys->fNuisanceParameter;
        }
        else{
            mySys2->fNuisanceParameter = mySys2->fName;
            TtHFitter::NPMAP[mySys2->fName] = mySys2->fName;
        }

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            mySys2->fTitle = (sys->fTitle)+"_Shape";
            TtHFitter::SYSTMAP[mySys2->fName] = mySys2->fTitle;
        }
        fFitter->fNSyst++;

        
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            sam = fFitter->fSamples[i_smp];
            if(sam->fType == Sample::DATA) continue;
            if(!sam->fUseSystematics) continue;
            if(   (samples[0]=="all" || FindInStringVector(samples, sam->fName)>=0 )
               && (exclude[0]==""    || FindInStringVector(exclude, sam->fName)<0 ) ){
                sam->AddSystematic(mySys2);
            }
        }
    }
    delete sys;
    
    return 0;
}

int ConfigReader::PostConfig(){
    // if StatOnly, also sets to OFF the MC stat
    if(fFitter->fStatOnly){
        WriteInfoStatus("ConfigReader::PostConfig","StatOnly option is setting to OFF the MC-stat (gammas) as well.");
        WriteInfoStatus("ConfigReader::PostConfig","To keep them on use the command line option 'Systematics=NONE'");
        WriteInfoStatus("ConfigReader::PostConfig","or comment out all Systematics in config file.");
        fFitter->SetStatErrorConfig( false, 0. );
    }

    // add nuisance parameter - systematic title correspondence
    for(auto syst : fFitter->fSystematics){
        if(syst->fNuisanceParameter!=syst->fName) TtHFitter::SYSTMAP[syst->fNuisanceParameter] = syst->fTitle;
    }
    // add nuisance parameter - norm-factor title correspondence & fix nuisance parameter
    for(auto norm : fFitter->fNormFactors){
        if(TtHFitter::NPMAP[norm->fName]=="") TtHFitter::NPMAP[norm->fName] = norm->fName;
        if(norm->fNuisanceParameter!=norm->fName) TtHFitter::SYSTMAP[norm->fNuisanceParameter] = norm->fTitle;
    }

    // morphing
    if (fFitter->fRunMorphing){
        // template fitting stuff
        fFitter->fTemplateWeightVec = fFitter->GetTemplateWeightVec(fFitter->fTemplateInterpolationOption);
        for(const TtHFit::TemplateWeight& itemp : fFitter->fTemplateWeightVec){
            std::string normName = "morph_"+itemp.name+"_"+ReplaceString(std::to_string(itemp.value),"-","m");
            TtHFitter::SYSTMAP[normName] = itemp.function;
            TtHFitter::NPMAP[normName]   = itemp.name;
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

    return 0;
}

std::string ConfigReader::CheckName( const std::string &name ){
    if( std::isdigit( name.at(0) ) ){
        WriteErrorStatus("ConfigReader::CheckName", "Failed to browse name: " + name + ". A number has been detected at the first position of the name.");
        WriteErrorStatus("ConfigReader::CheckName", "           This can lead to unexpected behaviours in HistFactory. Please change the name. ");
        WriteErrorStatus("ConfigReader::CheckName", "           The code is about to crash.");
        std::abort();
    } else {
        return name;
    }
}

bool ConfigReader::ConfigHasNTUP(ConfigSet* confSet){
    if (confSet->Get("Variable") != "" || confSet->Get("VariableForSample") != "" || confSet->Get("Selection") != "" || confSet->Get("NtupleName") != "" || confSet->Get("NtupleNameSuff") != "" || confSet->Get("MCweight") != "" || confSet->Get("NtuplePathSuff") != "" || confSet->Get("NtupleFile") != "" || confSet->Get("NtupleFiles") != "" || confSet->Get("NtupleNames") != "" || confSet->Get("NtuplePath") != "" || confSet->Get("NtuplePaths") != "" || confSet->Get("NtuplePathSuffs") != "") return true;
    else return false;
}

bool ConfigReader::ConfigHasHIST(ConfigSet* confSet){
    if (confSet->Get("HistoFile") != "" || confSet->Get("HistoName") != "" || confSet->Get("HistoPathSuff") != "" || confSet->Get("HistoPathSuffs") != "" || confSet->Get("HistoPath") != "" ) return true;
    else return false;
}
