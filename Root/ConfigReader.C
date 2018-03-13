#include "TtHFitter/ConfigReader.h"

#include "TtHFitter/TtHFit.h"
#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/Common.h"
#include "TtHFitter/Region.h"


ConfigReader::ConfigReader(TtHFit *fitter){
    fFitter = fitter;
    WriteInfoStatus("ConfigReader::ConfigReader", "Started reading the config");
}

ConfigReader::~ConfigReader(){
    WriteInfoStatus("ConfigReader::~ConfigReader", "Finished with the config reading");
}


int ConfigReader::ReadFullConfig(const std::string& fileName, const std::string& option){
    // initialize ConfigParser for the actual config
    fParser.ReadFile(fileName);

    // initialize checker COnfigParser to cross check the input
    ConfigParser refConfig;
    refConfig.ReadFile("jobSchema.config");
    int sc = fParser.CheckSyntax(&refConfig);

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

    return 0;
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
    std::vector<std::string> fRegNames;

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

std::string ConfigReader::CheckName( const std::string &name ){
    if( std::isdigit( name.at(0) ) ){
        WriteErrorStatus("ConfigReader::CheckName", "ERROR in browsing name: " + name + ". A number has been detected at the first position of the name.");
        WriteErrorStatus("ConfigReader::CheckName", "           This can lead to unexpected behaviours in HistFactory. Please change the name. ");
        WriteErrorStatus("ConfigReader::CheckName", "           The code is about to crash.");
        std::abort();
    } else {
        return name;
    }
}

bool ConfigReader::ConfigHasNTUP(ConfigSet* confSet){
    if (confSet->Get("Variable") != "" || confSet->Get("VariableForSample") != "" || confSet->Get("Selection") != "" || confSet->Get("NtupleName") != "" || confSet->Get("NtupleNameSuff") != "" || confSet->Get("MCweight") != "" || confSet->Get("NtuplePathSuff") != "") return true;
    else return false;
}

bool ConfigReader::ConfigHasHIST(ConfigSet* confSet){
    if (confSet->Get("HistoFile") != "" || confSet->Get("HistoName") != "" || confSet->Get("HistoPathSuff") != "" || confSet->Get("HistoPathSuffs") != "" ) return true;
    else return false;
}
