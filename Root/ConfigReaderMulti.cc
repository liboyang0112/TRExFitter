// Class include
#include "TRExFitter/MultiFit.h"

// Framework includes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReaderMulti.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/TRExFit.h"

// ROOT includes
#include "TSystem.h"

// c++ includes
#include <algorithm>

#define _STRINGIZE(x) #x
#define STRINGIZE(x) _STRINGIZE(x)

//_______________________________________________________________________________________
//
ConfigReaderMulti::ConfigReaderMulti(MultiFit *multiFitter) :
    fMultiFitter(multiFitter),
    fGlobalSuffix(""),
    fOnlyLHscan("") {
    WriteInfoStatus("ConfigReaderMulti::ConfigReaderMulti", "Started reading the config for multifit");
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadFullConfig(const std::string& fileName, const std::string& opt, const std::string& option){
    // initialize ConfigParser for the actual config
    fParser.ReadFile(fileName);

    // initialize checker COnfigParser to cross check the input
    ConfigParser refConfig;
    std::string homeArea("$TREXFITTER_HOME");
#ifdef TREXFITTER_HOME
    homeArea = std::string(STRINGIZE(TREXFITTER_HOME));
#endif
    refConfig.ReadFile(gSystem->ExpandPathName((homeArea+"/multiFitSchema.config").c_str()));
    int sc = fParser.CheckSyntax(&refConfig);

    if (sc != 0) return sc;

    // syntax of the config is ok
    // read different types of settings
    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc+= ReadJobOptions();

    sc+= ReadLimitOptions();

    sc+= ReadSignificanceOptions();

    sc+= ReadFitOptions(opt, option);

    // make directory
    gSystem->mkdir(fMultiFitter->fOutDir.c_str());

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadCommandLineOptions(const std::string &option){
    // Read options (to skip stuff, or include only some regions, samples, systs...)
    // Syntax: .. .. Regions=ge4jge2b:Exclude=singleTop,wjets
    std::map< std::string,std::string > optMap;

    int sc(0);

    std::vector< std::string > optVec = Common::Vectorize(option,':');
    for(const std::string& iopt : optVec){
        std::vector< std::string > optPair;
        optPair = Common::Vectorize(iopt,'=');
        if (optPair.size() < 2){
            WriteErrorStatus("ConfigReaderMulti::ReadCommandLineOptions", "Cannot read your command line option, please check this!");
            ++sc;
        } else {
            optMap[optPair[0]] = optPair[1];
        }
    }

    if(optMap["Ranking"]!=""){
        fMultiFitter->fRankingOnly = optMap["Ranking"];
    }
    if(optMap["GroupedImpact"]!=""){
        fMultiFitter->fGroupedImpactCategory = optMap["GroupedImpact"];
    }

    if(optMap["Suffix"]!=""){
        fGlobalSuffix = optMap["Suffix"];
    }
    if(optMap["LHscan"]!=""){
        fOnlyLHscan = optMap["LHscan"];
    }
    if(optMap["Parallel2Dscan"]!="" && optMap["Parallel2Dscan"]!="FALSE"){
        fMultiFitter->fParal2D = true;
    }
    if(optMap["Parallel2DscanStep"]!=""){
        fMultiFitter->fParal2Dstep = atoi(optMap["Parallel2DscanStep"].c_str());
    }

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadJobOptions() {

    int sc(0);

    std::string param = "";
    ConfigSet *confSet = fParser.GetConfigSet("MultiFit");
    if (confSet == nullptr){
        WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "Cannot find 'MultiFit' in your config which is required. Please check this!");
        ++sc;
    }

    fMultiFitter->fName = Common::CheckName(confSet->GetValue());

    // Set OutputDir
    param = confSet->Get("OutputDir");
    if(param !=""){
        fMultiFitter->fDir = Common::CheckName(param);
        if(fMultiFitter->fDir.back() != '/') fMultiFitter->fDir += '/';
        fMultiFitter->fOutDir = fMultiFitter->fDir + fMultiFitter->fName;
        gSystem->mkdir((fMultiFitter->fOutDir).c_str(), true);
    }
    else{
        fMultiFitter->fOutDir = "./" + fMultiFitter->fName;
    }

    // Set Label
    param = confSet->Get("Label");
    if( param != "") fMultiFitter->fLabel = Common::RemoveQuotes(param);

    // Set LumiLabel
    param = confSet->Get("LumiLabel");
    if( param != "") fMultiFitter->fLumiLabel = Common::RemoveQuotes(param);

    // Set CmeLabel
    param = confSet->Get("CmeLabel");
    if( param != "") fMultiFitter->fCmeLabel = Common::RemoveQuotes(param);

    // Set Label of combination
    param = confSet->Get("CombiLabel");
    if( param != "") fMultiFitter->fCombiLabel = Common::RemoveQuotes(param);

    // Set SaveSuf
    param = confSet->Get("SaveSuf");
    if( param != "") fMultiFitter->fSaveSuf = Common::RemoveQuotes(param);
    else fMultiFitter->fSaveSuf             = fGlobalSuffix;

    // Set ShowObserved
    param = confSet->Get("ShowObserved");
    if( param != ""){
        fMultiFitter->fShowObserved = Common::StringToBoolean(param); 
    }

    // Set LimitTitle
    param = confSet->Get("LimitTitle");
    if( param != "") fMultiFitter->fLimitTitle = Common::RemoveQuotes(param);
    // --- these lines should be removed at some point, since now we protect the text inside ""
    if(fMultiFitter->fLimitTitle.find("95CL")!=std::string::npos){
         fMultiFitter->fLimitTitle.replace(fMultiFitter->fLimitTitle.find("95CL"),4,"95% CL");
    }
    // ---

    // Ser POITitle
    param = confSet->Get("POITitle");
    if( param != "") fMultiFitter->fPOITitle = Common::RemoveQuotes(param);

    // Set CompareLimits
    param = confSet->Get("CompareLimits");
    if( param != ""){
        fMultiFitter->fCompareLimits = Common::StringToBoolean(param);
    }

    // Ser ComparePOI
    param = confSet->Get("ComparePOI");
    if( param != ""){
        fMultiFitter->fComparePOI = Common::StringToBoolean(param);
    }

    // Set ComparePulls
    param = confSet->Get("ComparePulls");
    if( param != "" ){
        fMultiFitter->fComparePulls = Common::StringToBoolean(param);
    }

    // Set PlotCombCorrMatrix
    param = confSet->Get("PlotCombCorrMatrix");
    if (param != ""){
        fMultiFitter->fPlotCombCorrMatrix = Common::StringToBoolean(param);
    }

    // Set CorrelationThreshold
    param = confSet->Get("CorrelationThreshold");
    if( param != ""){
        TRExFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }

    // Set UseGammasForCorr
    param = confSet->Get("UseGammasForCorr");
    if( param != ""){
        fMultiFitter->fuseGammasForCorr = Common::StringToBoolean(param);
    }

    // Set Combine
    param = confSet->Get("Combine");
    if( param != ""){
        fMultiFitter->fCombine = Common::StringToBoolean(param);
    }

    // Set Compare
    param = confSet->Get("Compare");
    if( param != ""){
        fMultiFitter->fCompare = Common::StringToBoolean(param);
    }

    // Set StatOnly
    param = confSet->Get("StatOnly");
    if( param != "" ){
        fMultiFitter->fStatOnly = Common::StringToBoolean(param);
    }

    // Set IncludeStatOnly
    param = confSet->Get("IncludeStatOnly");
    if( param != ""){
        fMultiFitter->fIncludeStatOnly = Common::StringToBoolean(param);
    }

    // Set POIName
    param = confSet->Get("POIName");
    if( param != "" ) {
        const auto tmp =  Common::Vectorize(param, ',');
        for (const auto& i : tmp) {
            fMultiFitter->AddPOI(i);
        }
    }

    // Set POILabel
    param = confSet->Get("POILabel");
    if( param != "" ) fMultiFitter->fPOIName = Common::RemoveQuotes(param);

    // Set POINominal
    param = confSet->Get("POINominal");
    if( param != "" ) fMultiFitter->fPOINominal = std::stof(param);

    // Set POIRange
    param = confSet->Get("POIRange");
    if( param != ""){
        const auto tmp = Common::Vectorize(param, ',');   
        if (tmp.size() != fMultiFitter->fPOIs.size()) {
            WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'POIRange' option but you didn't pass the same number of parametrs as the number of POIs. Please check this!");
            ++sc;
        }
        for (const auto& irange : tmp) {
            const auto vec = Common::Vectorize(irange,':');
            if (vec.size()==2 ) {
                fMultiFitter->fPOIMin.emplace_back(atof( vec[0].c_str() ));
                fMultiFitter->fPOIMax.emplace_back(atof( vec[1].c_str() ));
            } else {
                WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'POIRange' option but you didn't provide valid setting. Please check this!");
                ++sc;
            }
        }
    }

    // Set LimitMax
    param = confSet->Get("LimitMax");
    if( param != "" ) {
        fMultiFitter->fLimitMax = atof( param.c_str() );
    }

    // Set POIPrecision
    param = confSet->Get("POIPrecision");
    if( param != "" ) {
        const auto tmp = Common::Vectorize(param, ',');
        if (tmp.size() != fMultiFitter->fPOIs.size()) {
            WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'POIPrecision' option but you didn't pass the same number of parametrs as the number of POIs. Please check this!");
            ++sc;
        }
        for (const auto& i : tmp) { 
            fMultiFitter->fPOIPrecision.emplace_back(Common::RemoveQuotes(i).c_str());
        }
    }

    //Set DataName
    param = confSet->Get("DataName");
    if( param != "" ) fMultiFitter->fDataName = Common::RemoveQuotes(param);

    // Set FitType
    param = confSet->Get("FitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param=="SPLUSB")      fMultiFitter->fFitType = 1;
        else if(param=="BONLY")  fMultiFitter->fFitType = 2;
        else if(param=="UNFOLDING")  fMultiFitter->fFitType = 3;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'FitType' option but you didn't provide valid setting. Using default (SPLUSB)");
            fMultiFitter->fFitType = 1;
        }
    }

    // Set CombineChByCh
    param = confSet->Get("CombineChByCh");
    if( param != ""){
        fMultiFitter->fCombineChByCh = Common::StringToBoolean(param);
    }

    // Set NPCategories
    param = confSet->Get("NPCategories");
    if( param != "" ) {
        std::vector<std::string> categ = Common::Vectorize(param,',');
        for(const std::string &icat : categ)
            fMultiFitter->fNPCategories.push_back(icat);
    }

    // Set SetRandomInitialNPval
    param = confSet->Get("SetRandomInitialNPval");
    if( param != ""){
        fMultiFitter->fUseRnd = true;
        fMultiFitter->fRndRange = atof(param.c_str());
    }

    // Set SetRandomInitialNPvalSeed
    param = confSet->Get("SetRandomInitialNPvalSeed");
    if( param != ""){
        fMultiFitter->fRndSeed = atol(param.c_str());
    }

    // Set NumCPU
    param = confSet->Get("NumCPU");
    if( param != "" ){
        fMultiFitter->fCPU = atoi( param.c_str());
    }

    // Set FastFit
    param = confSet->Get("FastFit");
    if (param != ""){
        fMultiFitter->fFastFit = Common::StringToBoolean(param);
    }

    // Set FastFitForRanking
    param = confSet->Get("FastFitForRanking");
    if (param != ""){
        fMultiFitter->fFastFitForRanking = Common::StringToBoolean(param);
    }

    // Set NuisParListFile
    param = confSet->Get("NuisParListFile");
    if( param != "" ) fMultiFitter->fNuisParListFile = Common::RemoveQuotes(param);

    // Set PlotSoverB
    param = confSet->Get("PlotSoverB");
    if (param != ""){
        fMultiFitter->fPlotSoverB = Common::StringToBoolean(param);
    }

    // Set SignalTitle
    param = confSet->Get("SignalTitle");
    if( param != "" ) fMultiFitter->fSignalTitle = Common::RemoveQuotes(param);

    // Set FitResultsFile
    param = confSet->Get("FitResultsFile");
    if( param != "" ) fMultiFitter->fFitResultsFile = Common::RemoveQuotes(param);

    // Set LimitsFile
    param = confSet->Get("LimitsFile");
    if( param != "" ) fMultiFitter->fLimitsFile = Common::RemoveQuotes(param);

    // Set BonlySuffix
    param = confSet->Get("BonlySuffix");
    if( param != "" ) fMultiFitter->fBonlySuffix = Common::RemoveQuotes(param);

    // Set ShowSystForPOI
    param = confSet->Get("ShowSystForPOI");
    if( param != ""){
        fMultiFitter->fShowSystForPOI = Common::StringToBoolean(param);
    }

    // Set GetGoodnessOfFit
    param = confSet->Get("GetGoodnessOfFit");
    if( param != "" ){
        fMultiFitter->fGetGoodnessOfFit = Common::StringToBoolean(param);
    }

    // Set doLHscan
    param = confSet->Get("doLHscan");
    if( param != "" ){
        if (fOnlyLHscan==""){
            fMultiFitter->fVarNameLH = Common::Vectorize(param,',');
        } else {
            fMultiFitter->fVarNameLH.emplace_back(fOnlyLHscan);
        }
    }

    // Set do2DLHscan
    param = confSet->Get("do2DLHscan");
    if( param != "" ){
        const std::vector<std::string> tmp = Common::Vectorize(param,':');
        for (const auto& ivec : tmp){
            const std::vector<std::string> v = Common::Vectorize(ivec,',');
            if (v.size() == 2){
                fMultiFitter->fVarName2DLH.emplace_back(v);
            } else {
                WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "You specified 'do2DLHscan' option but did not provide correct input. Ignoring");
            }
        }
    }

    // Set LHscanMin
    param = confSet->Get("LHscanMin");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'LHscanMin' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanMin = std::stof(param);
        }
    }

    // Set LHscanMax
    param = confSet->Get("LHscanMax");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'LHscanMax' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanMax = std::stof(param);
        }
    }
    // Set LHscanMin for second variable
    param = confSet->Get("LHscanMinY");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'LHscanMinY' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanMinY = std::stof(param);
        }
    }

    // Set LHscanMax for second variable
    param = confSet->Get("LHscanMaxY");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'LHscanMaxY' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanMaxY = std::stof(param);
        }
    }

    // Set LHscanSteps for second variable
    param = confSet->Get("LHscanStepsY");
    if ( param != "" ) {
        if (fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'LHscanStepsY' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanStepsY = std::stoi(param);
            if(fMultiFitter->fLHscanStepsY < 3 || fMultiFitter->fLHscanStepsY > 100){
                WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "LHscanSteps is smaller than 3 or larger than 100, setting to defaut (30)");
                fMultiFitter->fLHscanStepsY = fMultiFitter->fLHscanSteps;
            }
        }
    }
    else {
        fMultiFitter->fLHscanStepsY = fMultiFitter->fLHscanSteps;
    }

    // Set Paral2D
    param = confSet->Get("Parallel2Dscan");
    if ( param != "" ) {
        fMultiFitter->fParal2D = Common::StringToBoolean(param);
    }

    // Set Paral2Dstep
    param = confSet->Get("Parallel2DscanStep");
    if ( param != "" ) {
        fMultiFitter->fParal2Dstep = std::atoi( param.c_str());
        if (fMultiFitter->fParal2Dstep < 1 || fMultiFitter->fParal2Dstep>=fMultiFitter->fLHscanSteps ){
            WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "You specified a step for 2D LHscan outside the allowed range.");
            ++sc;
        }
    }

    // Set ShowTotalOnly
    param = confSet->Get("ShowTotalOnly");
    if ( param != "" ) {
        fMultiFitter->fShowTotalOnly = Common::StringToBoolean(param);
    }

    // Set LHscanSteps
    param = confSet->Get("LHscanSteps");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'LHscanSteps' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanSteps = std::stoi(param);
            if(fMultiFitter->fLHscanSteps < 3 || fMultiFitter->fLHscanSteps > 100){
                WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "LHscanSteps is smaller than 3 or larger than 100, setting to defaut (30)");
                fMultiFitter->fLHscanSteps = 30;
            }
        }
    }

    // Set PlotOptions
    param = confSet->Get("PlotOptions");
    if( param != ""){
        std::vector<std::string> vec = Common::Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TRExFitter::PREFITONPOSTFIT= true;
    }

    // Set POIAsimovl
    param = confSet->Get("POIAsimov");
    if( param != ""){
        fMultiFitter->fPOIAsimov = std::stof(param);
    }

    // Set POIInitial
    param = confSet->Get("POIInitial");
    if( param != ""){
        fMultiFitter->fPOIInitial = std::stof(param);
    }
    
    // Set HEPDataFormat
    param = confSet->Get("HEPDataFormat");
    if( param != ""){
        fMultiFitter->fHEPDataFormat = Common::StringToBoolean(param);
    }
    
    // Set FitStrategy
    param = confSet->Get("FitStrategy");
    if (param != "") {
        fMultiFitter->fFitStrategy = stoi(param);
        if (fMultiFitter->fFitStrategy > 2) {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "FitStrategy > 2, setting to default (-1)");
            fMultiFitter->fFitStrategy = -1;
        }
    }

    // Set BinnedLikelihoodOptimization
    param = confSet->Get("BinnedLikelihoodOptimization");
    if (param != "") {
        fMultiFitter->fBinnedLikelihood = Common::StringToBoolean(param);
    }
    
    // Set UsePOISinRanking
    param = confSet->Get("UsePOISinRanking");
    if (param != "") {
        fMultiFitter->fUsePOISinRanking = Common::StringToBoolean(param);
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReaderMulti::ReadLimitOptions(){

    int sc(0);

    std::string param = "";

    ConfigSet* confSet = fParser.GetConfigSet("Limit");
    if (confSet == nullptr){
        WriteDebugStatus("ConfigReader::ReadLimitOptions", "You do not have Limit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set POI
    param = confSet->Get("POI");
    if( param != "" ){
        fMultiFitter->fPOIforLimit = param;
        fMultiFitter->AddPOI(param);
    }
    
    // Set LimitBlind
    param = confSet->Get("LimitBlind");
    if( param != "" ){
        fMultiFitter->fLimitIsBlind = Common::StringToBoolean(param);
    }

    // Set SignalInjection
    param = confSet->Get("SignalInjection");
    if( param != "" ){
        fMultiFitter->fSignalInjection = Common::StringToBoolean(param);
    }

    // Set SignalInjectionValue
    param = confSet->Get("SignalInjectionValue");
    if( param != "" ){
        fMultiFitter->fSignalInjectionValue = std::stof(param);
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fMultiFitter->fLimitParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fMultiFitter->fLimitParamValue = std::stof(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fMultiFitter->fLimitOutputPrefixName = param;
    }

    param = confSet->Get("ConfidenceLevel");
    if( param != "" ){
        double conf = std::stod(param);
        if (conf <= 0 || conf >= 1){
            WriteWarningStatus("ConfigReaderMulti::ReadLimitOptions", "Confidence level is <= 0 or >=1. Setting to default 0.95");
            conf = 0.95;
        }
        fMultiFitter->fLimitsConfidence = conf;
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReaderMulti::ReadSignificanceOptions(){

    int sc(0);

    std::string param = "";

    ConfigSet* confSet = fParser.GetConfigSet("Significance");
    if (confSet == nullptr){
        WriteDebugStatus("ConfigReaderMulti::ReadSignificanceOptions", "You do not have Significance option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set POI
    param = confSet->Get("POI");
    if( param != "" ){
        fMultiFitter->fPOIforSig = param;
        fMultiFitter->AddPOI(param);
    }
    
    // Set LimitBlind
    param = confSet->Get("SignificanceBlind");
    if( param != "" ){
        fMultiFitter->fSignificanceIsBlind = Common::StringToBoolean(param);
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
        fMultiFitter->fSignificancePOIAsimov = atof(param.c_str());
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fMultiFitter->fSignificanceParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fMultiFitter->fSignificanceParamValue = std::stof(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fMultiFitter->fSignificanceOutputPrefixName = param;
    }

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadFitOptions(const std::string& opt, const std::string& options){

    int sc(0);

    int nFit = 0;

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("Fit", nFit);
        if (confSet == nullptr) break;
        nFit++;

        // Set Options
        std::string fullOptions;
        std::string param = confSet->Get("Options");
        if(param!="" && options!="") fullOptions = options+";"+Common::RemoveQuotes(param);
        else if(param!="") fullOptions = Common::RemoveQuotes(param);
        else fullOptions = options;

        // name
        fMultiFitter->fFitNames.push_back(Common::CheckName(confSet->GetValue()));

        // Set Label
        param = confSet->Get("Label");
        std::string label = Common::CheckName(confSet->GetValue());
        if(param!="") label = Common::RemoveQuotes(param);

        // Set suf
        param = confSet->Get("LoadSuf");
        std::string loadSuf = "";
        if(param!="") loadSuf = Common::RemoveQuotes(param);
        else          loadSuf = fGlobalSuffix;

        // config file
        std::string confFile = "";
        param = confSet->Get("ConfigFile");
        if(param!="") confFile = Common::RemoveQuotes(param);

        // workspace
        std::string wsFile = "";
        param = confSet->Get("Workspace");
        if(param!="") wsFile = Common::RemoveQuotes(param);

        bool useInFit(true);
        param = confSet->Get("UseInFit");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") {
                useInFit = false;
            }
        }

        bool useInComparison(true);
        param = confSet->Get("UseInComparison");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") {
                useInComparison = false;
            }
        }

        fMultiFitter->AddFitFromConfig(confFile, opt, fullOptions, label, loadSuf, wsFile, useInFit, useInComparison);

        // Set FitResultsFile
        param = confSet->Get("FitResultsFile");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fFitResultsFile = Common::RemoveQuotes(param);
        fMultiFitter->fLimitsFiles.push_back("");

        // Set LimitsFile
        param = confSet->Get("LimitsFile");
        if( param != "" ) fMultiFitter->fLimitsFiles[fMultiFitter->fFitList.size()-1] = Common::RemoveQuotes(param);

        // Set POIName
        param = confSet->Get("POIName");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fPOIs[0] = Common::RemoveQuotes(param);

        // Set Directory
        param = confSet->Get("Directory");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fName = Common::RemoveQuotes(param);

        // Set InputName
        param = confSet->Get("InputName");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fInputName = Common::RemoveQuotes(param);
    }

    if (nFit == 0){
        WriteErrorStatus("ConfigReaderMulti::ReadFitOptions", "You need to provide at least one 'Fit' option. Please check this!");
        ++sc;
    }

    return sc;
}
