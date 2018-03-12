#include "TtHFitter/ConfigReader.h"

#include "TtHFitter/TtHFit.h"
#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/Common.h"


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

    fConfSet = fParser.GetConfigSet("Job");
    if (fConfSet == nullptr){
        WriteErrorStatus("ConfigReader::ReadJobOptions", "You need to provide JOB settings!");
        exit(EXIT_FAILURE);
    }
    
    fFitter->fName = CheckName(fConfSet->GetValue());
    fFitter->fInputName = fFitter->fName;
    
    //Set DebugLevel
    param = fConfSet->Get("DebugLevel");
    if( param != "")  TtHFitter::SetDebugLevel( atoi(param.c_str()) );
    
    // Set outputDir
    param = fConfSet->Get("OutputDir");
    if(param != ""){
      fFitter->fDir = param;
      if(fFitter->fDir.back() != '/') fFitter->fDir += '/';
      fFitter->fName = fFitter->fDir + fFitter->fName;
      gSystem->mkdir((fFitter->fName).c_str(), true);
    }

    // Set Label
    param = fConfSet->Get("Label");
    if(param!="") fFitter->fLabel = param;
    else          fFitter->fLabel = fFitter->fName;

    // Set POI
    fFitter->SetPOI(CheckName(fConfSet->Get("POI")));

    //Set reading option
    param = fConfSet->Get("ReadFrom");
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
    param = fConfSet->Get("MergeUnderOverFlow");
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
        fFitter->AddHistoPath( fConfSet->Get("HistoPath") );
        if (fConfSet->Get("NtuplePath") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtuplePath:' option, ignoring.");
        }
        if (fConfSet->Get("NtuplePaths") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtuplePaths:' option, ignoring.");
        }
        if (fConfSet->Get("MCweight") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'MCweight:' option, ignoring.");
        }
        if (fConfSet->Get("Selection") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'Selection:' option, ignoring.");
        }
        if (fConfSet->Get("NtupleName") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified HIST option but you provided 'NtupleName:' option, ignoring.");
        }
    }
    // Setting for NTUP only
    if(fFitter->fInputType==1){
        if (fConfSet->Get("HistoPath") != ""){
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified NTUP option but you provided 'HistoPath:' option, ignoring.");
        }
        fFitter->SetNtupleFile( fConfSet->Get("NtupleFile") );
        if(fConfSet->Get("NtuplePath")!="") {
            fFitter->AddNtuplePath( fConfSet->Get("NtuplePath") ); 
        }
        param = fConfSet->Get("NtuplePaths");
        if( param != "" ){
            std::vector<std::string> paths = Vectorize( param,',' );
            for(const std::string& ipath : paths){
                fFitter->AddNtuplePath( ipath );
            }
        }
        param = fConfSet->Get("MCweight");
        if(param!="") fFitter->SetMCweight(param);

        param = fConfSet->Get("Selection");
        if(param!="") fFitter->SetSelection(param);
        fFitter->SetNtupleName( fConfSet->Get("NtupleName") );
    }
    
    // Set lumi
    param = fConfSet->Get("Lumi");
    if( param != "" ) fFitter->SetLumi( atof(param.c_str()) );

    // Set LumiScale
    param = fConfSet->Get("LumiScale");
    if( param != "" ){
        WriteWarningStatus("ConfigReader::ReadJobOptions", "\"LumiScale\" is only done for quick tests since it is inefficient.");
        WriteWarningStatus("ConfigReader::ReadJobOptions", "To normalize all the samples to the luminosity, use \"Lumi\" instead.");
        fFitter->fLumiScale = atof(param.c_str());
    }

    // Set TtresSmoothing
    param = fConfSet->Get("TtresSmoothing");
    if( param != ""){ 
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fFitter->fTtresSmoothing = true;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'TtresSmoothing' option but you didn't set it to TRUE. Using default (FALSE)");
            fFitter->fTtresSmoothing = false;
        }
    }

    // Set SystPruningShape
    param = fConfSet->Get("SystPruningShape");
    if( param != "") fFitter->fThresholdSystPruning_Shape = atof(param.c_str());

    // Set SystPruningNorm
    param = fConfSet->Get("SystPruningNorm");
    if( param != "")  fFitter->fThresholdSystPruning_Normalisation = atof(param.c_str());

    // Set SystLarge
    param = fConfSet->Get("SystLarge");
    if( param != "")  fFitter->fThresholdSystLarge = atof(param.c_str());

    // Set IntCodeOverall
    param = fConfSet->Get("IntCodeOverall");
    if( param != "")  fFitter->fIntCode_overall = atoi(param.c_str());

    // Set IntCodeShape
    param = fConfSet->Get("IntCodeShape");
    if( param != "")  fFitter->fIntCode_shape = atoi(param.c_str());

    // Set MCstatThreshold
    param = fConfSet->Get("MCstatThreshold");
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
    param = fConfSet->Get("MCstatConstraint");
    if( param != "")  fFitter->fStatErrCons = param;

    // Set UseGammaPulls
    param = fConfSet->Get("UseGammaPulls");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") fFitter->fUseGammaPulls = true;
        else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'UseGammaPulls' option but you didn't set it to TRUE. Using default (FALSE)");
            fFitter->fUseGammaPulls = false;
        }
    }

    // plotting options are in special function
    if (SetJobPlot() != 0) return 1;

    // Set TableOptions 
    param = fConfSet->Get("TableOptions");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fFitter->fTableOptions = param;
    }

    // Set CorrelationThreshold
    param = fConfSet->Get("CorrelationThreshold");
    if( param != ""){
        TtHFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }

    // Set HistoChecks
    param = fConfSet->Get("HistoChecks");
    if(param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "NOCRASH" ){
            TtHFitter::HISTOCHECKCRASH = false;
        } else {
            WriteWarningStatus("ConfigReader::ReadJobOptions", "You specified 'HistoChecks' option but you didn't set it to NOCRASH.");
        }
    }

    // Set LumiLabel
    param = fConfSet->Get("LumiLabel");
    if( param != "") fFitter->fLumiLabel = param;

    // Set CmeLabel
    param = fConfSet->Get("CmeLabel");
    if( param != "") fFitter->fCmeLabel = param;

    // Set SplitHistoFiles
    param = fConfSet->Get("SplitHistoFiles");
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
    param = fConfSet->Get("BlindingThreshold");
    if( param != ""){
        fFitter->fBlindingThreshold = atof(param.c_str());
    }

    // Set RankingMaxNP
    param = fConfSet->Get("RankingMaxNP");
    if( param != ""){
        fFitter->fRankingMaxNP = atoi(param.c_str());
    }

    // Set ReduceNPforRanking
    param = fConfSet->Get("ReduceNPforRanking");
    if( param != ""){
        fFitter->fReduceNPforRanking = atof(param.c_str());
    }

    // Set ImageFormat
    param = fConfSet->Get("ImageFormat");
    if( param != ""){
        std::vector<std::string> tmp = Vectorize(param,',');
        if (tmp.size() > 0) fFitter->fImageFormat = tmp.at(0);
        else {
            WriteErrorStatus("ConfigReader::ReadJobOptions", "You specified 'ImageFormat' option but we cannot split the setting. Please check");
        }
        TtHFitter::IMAGEFORMAT = tmp;
    }

    // Set StatOnly
    param = fConfSet->Get("StatOnly");
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
    param = fConfSet->Get("FixNPforStatOnly");
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
    param = fConfSet->Get("InputFolder");
    if( param != "" ){
        fFitter->fInputFolder = param;
    }

    // Set InputName
    param = fConfSet->Get("InputName");
    if( param != "" ){
        fFitter->fInputName = param;
    }

    // Set WorkspaceFileName
    param = fConfSet->Get("WorkspaceFileName");
    if( param != "" ){
        fFitter->fWorkspaceFileName = param;
    }

    // Set KeepPruning
    param = fConfSet->Get("KeepPruning");
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
    param = fConfSet->Get("AtlasLabel");
    if( param != "" ){
        fFitter->fAtlasLabel = param;
    }
        
    // Set CleanTables
    param = fConfSet->Get("CleanTables");
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
    param = fConfSet->Get("SystCategoryTables");
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
    param = fConfSet->Get("Suffix");
    if( param != "" ){
        fFitter->fSuffix = param;
    }
    
    // Set SaveSuffix
    param = fConfSet->Get("SaveSuffix");
    if( param != "" ){
        fFitter->fSaveSuffix = param;
    }

    // Set HideNP
    param = fConfSet->Get("HideNP");
    if( param != "" ){
        fFitter->fVarNameHide = Vectorize(param,',');
    }
        
    // Set RegionGroups
    param = fConfSet->Get("RegionGroups");
    if( param != "" ) {
        std::vector<std::string> groups = Vectorize(param,',');
        for(const std::string& igroup : groups) fFitter->fRegionGroups.push_back(igroup);
    }

    // Set KeepPrefitBlindedBins
    param = fConfSet->Get("KeepPrefitBlindedBins");
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
    param = fConfSet->Get("CustomAsimov");
    if( param != "" ){
        fFitter->fCustomAsimov = param;
    }

    // Set RandomPOISeed
    param = fConfSet->Get("RandomPOISeed");
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
    param = fConfSet->Get("GetChi2");
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
    param = fConfSet->Get("DoTables");
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
    param = fConfSet->Get("CustomFunctions");
    if( param != "" ) {
        fFitter->fCustomFunctions = Vectorize(param,',');
    }

    // Set Bootstrap
    param = fConfSet->Get("Bootstrap");
    if( param != "" ){
        fFitter->fBootstrap = param;
    }

    // Set RunROOTMacros
    param = fConfSet->Get("RunROOTMacros");
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
    param = fConfSet->Get("DecorrSuff");
    if( param != ""){
        fFitter->fDecorrSuff = param;
    }

    // success
    return 0;
}

int ConfigReader::SetJobPlot(){

    // Plot option
    std::string param = fConfSet->Get("PlotOptions");
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
    param = fConfSet->Get("PlotOptionsSummary");
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
    param = fConfSet->Get("SystControlPlots");
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
    param = fConfSet->Get("SystDataPlots");
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
    param = fConfSet->Get("SystErrorBars");
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
    param = fConfSet->Get("GuessMCStatEmptyBins");
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
    param = fConfSet->Get("SuppressNegativeBinWarnings"); 
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
    param = fConfSet->Get("SignalRegionsPlot");
    if(param != ""){
        fFitter->fRegionsToPlot = Vectorize(param,',');
    }

    // Set SummaryPlotRegions
    param = fConfSet->Get("SummaryPlotRegions");
    if(param != ""){
        fFitter->fSummaryPlotRegions = Vectorize(param,',');
    }

    // Set SummaryPlotLabels
    param = fConfSet->Get("SummaryPlotLabels");
    if(param != ""){
        fFitter->fSummaryPlotLabels = Vectorize(param,',');
    }

    // Set SummaryPlotValidationRegions
    param = fConfSet->Get("SummaryPlotValidationRegions");  
    if(param != ""){
        fFitter->fSummaryPlotValidationRegions = Vectorize(param,',');
    }

    // Set SummaryPlotValidationLabels
    param = fConfSet->Get("SummaryPlotValidationLabels");
    if(param != ""){
        fFitter->fSummaryPlotValidationLabels = Vectorize(param,',');
    }

    // Set SummaryPlotYmin
    param = fConfSet->Get("SummaryPlotYmin");
    if(param != "") fFitter->fYmin = atof(param.c_str());

    // Set SummaryPlotYmax
    param = fConfSet->Get("SummaryPlotYmax");
    if(param != "") fFitter->fYmax = atof(param.c_str());

    // Set RatioYmin
    param = fConfSet->Get("RatioYmin");
    if(param != "") {
        fFitter->fRatioYmin = atof(param.c_str());
        fFitter->fRatioYminPostFit = fFitter->fRatioYmin;
    }
    
    // Set RatioYmax
    param = fConfSet->Get("RatioYmax");
    if(param != ""){
        fFitter->fRatioYmax = atof(param.c_str());
        fFitter->fRatioYmaxPostFit = fFitter->fRatioYmax; 
    }  

    // Set RatioYminPostFit
    param = fConfSet->Get("RatioYminPostFit");
    if(param != "") fFitter->fRatioYminPostFit = atof(param.c_str());

    // Set RatioYmaxPostFit
    param = fConfSet->Get("RatioYmaxPostFit");
    if(param != "") fFitter->fRatioYmaxPostFit = atof(param.c_str());    

    // Set DoSummaryPlot
    param = fConfSet->Get("DoSummaryPlot");
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
    param = fConfSet->Get("DoMergedPlot");
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
    param = fConfSet->Get("DoSignalRegionsPlot");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param == "TRUE" )  fFitter->fDoSignalRegionsPlot = true;
        else if( param == "FALSE" ) fFitter->fDoSignalRegionsPlot = false;
        else {
            WriteWarningStatus("ConfigReader::SetJobPlot", "You specified 'DoSignalRegionsPlot' option but didnt provide valid parameter. Using default (false)");
            fFitter->fDoSignalRegionsPlot = false;
        }
    }
    param = fConfSet->Get("DoPieChartPlot");
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
    param = fConfSet->Get("RankingPlot");
    if( param != ""){
        fFitter->fRankingPlot = param;
    }

    return 0;
}

int ConfigReader::ReadGeneralOptions(){
    fConfSet = fParser.GetConfigSet("Options");        
    if (fConfSet != nullptr){
        for(int i=0; i < fConfSet->GetN(); i++){
            if(fConfSet->GetConfigValue(i) != ""){
                TtHFitter::OPTION[fConfSet->GetConfigName(i)] = atof(fConfSet->GetConfigValue(i).c_str());
            }
        }
    } else {
        WriteDebugStatus("ConfigReader::ReadGeneralOptions", "You do not have 'Options' option in the config. It is ok, we just want to let you know.");
    }

    return 0;
}

int ConfigReader::ReadFitOptions(){
    std::string param = "";

    fConfSet = fParser.GetConfigSet("Fit");
    if (fConfSet == nullptr){
        WriteInfoStatus("ConfigReader::ReadFitOptions", "You do not have Fit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up 
    }

    //Set FitType
    param = fConfSet->Get("FitType");
    if( param != "" && fFitter->fFitType == TtHFit::UNDEFINED ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "SPLUSB" ){
            fFitter->SetFitType(TtHFit::SPLUSB);
        }
        else if( param == "BONLY" ){
            fFitter->SetFitType(TtHFit::BONLY);
        }
        else{
            WriteErrorStatus("ConfigReader::ReadFitOptions", "Unknown FitType argument : " + fConfSet->Get("FitType"));
            return 1;
        }
    }
    else if( fFitter->fFitType == TtHFit::UNDEFINED ){
        WriteInfoStatus("ConfigReader::ReadFitOptions","Setting default fit Type SPLUSB");
        fFitter->SetFitType(TtHFit::SPLUSB);
    }

    // Set FitRegion
    param = fConfSet->Get("FitRegion");
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
                WriteErrorStatus("ConfigReader::ReadFitOptions", "Unknown FitRegion argument : " + fConfSet->Get("FitRegion"));
                return 1;
            }
        }
    }
    
    // Set FitBlind
    param = fConfSet->Get("FitBlind");
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
    param = fConfSet->Get("POIAsimov");
    if( param != "" ){
         fFitter->fFitPOIAsimov = atof(param.c_str());
    }

    // Set NPValues
    param = fConfSet->Get("NPValues");
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
    param = fConfSet->Get("FixNPs");
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
    param = fConfSet->Get("doLHscan");
    if( param != "" ){
        fFitter->fVarNameLH = Vectorize(param,',');
    }
    
    // Set UseMinos
    param = fConfSet->Get("UseMinos");
    if( param != "" ){
        fFitter->fVarNameMinos = Vectorize(param,',');
    }

    // Set SetRandomInitialNPval
    param = fConfSet->Get("SetRandomInitialNPval");
    if( param != ""){
        fFitter->fUseRnd = true;
        fFitter->fRndRange = std::atof(param.c_str());
    }

    // Set SetRandomInitialNPvalSeed
    param = fConfSet->Get("SetRandomInitialNPvalSeed");
    if( param != ""){
        fFitter->fRndSeed = std::atol(param.c_str());
    }

    // Set NumCPU
    param = fConfSet->Get("NumCPU");
    if( param != "" ){
        TtHFitter::NCPU = std::atoi( param.c_str());
    }

    // Set StatOnlyFit
    param = fConfSet->Get("StatOnlyFit");
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
    param = fConfSet->Get("GetGoodnessOfFit");
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

