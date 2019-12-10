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
    std::vector< std::string > optVec;

    optVec = Vectorize(option,':');
    for(const std::string& iopt : optVec){
        std::vector< std::string > optPair;
        optPair = Vectorize(iopt,'=');
        if (optPair.size() < 2){
            WriteErrorStatus("ConfigReaderMulti::ReadCommandLineOptions", "Cannot read your command line option, please check this!");
            return 1;
        }
        optMap[optPair[0]] = optPair[1];
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

    return 0;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadJobOptions(){
    std::string param = "";
    ConfigSet *confSet = fParser.GetConfigSet("MultiFit");
    if (confSet == nullptr){
        WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "Cannot find 'MultiFit' in your config which is required. Please check this!");
        return 1;
    }

    fMultiFitter->fName = CheckName(confSet->GetValue());

    // Set OutputDir
    param = confSet->Get("OutputDir");
    if(param !=""){
        fMultiFitter->fDir = CheckName(param);
        if(fMultiFitter->fDir.back() != '/') fMultiFitter->fDir += '/';
        fMultiFitter->fOutDir = fMultiFitter->fDir + fMultiFitter->fName;
        gSystem->mkdir((fMultiFitter->fOutDir).c_str(), true);
    }
    else{
        fMultiFitter->fOutDir = "./" + fMultiFitter->fName;
    }

    // Set Label
    param = confSet->Get("Label");
    if( param != "") fMultiFitter->fLabel = RemoveQuotes(param);

    // Set LumiLabel
    param = confSet->Get("LumiLabel");
    if( param != "") fMultiFitter->fLumiLabel = RemoveQuotes(param);

    // Set CmeLabel
    param = confSet->Get("CmeLabel");
    if( param != "") fMultiFitter->fCmeLabel = RemoveQuotes(param);

    // Set SaveSuf
    param = confSet->Get("SaveSuf");
    if( param != "") fMultiFitter->fSaveSuf = RemoveQuotes(param);
    else fMultiFitter->fSaveSuf             = fGlobalSuffix;

    // Set ShowObserved
    param = confSet->Get("ShowObserved");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") fMultiFitter->fShowObserved = true;
        else if (param == "FALSE") fMultiFitter->fShowObserved = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'ShowObserved' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fShowObserved = false;
        }
    }

    // Set LimitTitle
    param = confSet->Get("LimitTitle");
    if( param != "") fMultiFitter->fLimitTitle = RemoveQuotes(param);
    // --- these lines should be removed at some point, since now we protect the text inside ""
    if(fMultiFitter->fLimitTitle.find("95CL")!=std::string::npos){
         fMultiFitter->fLimitTitle.replace(fMultiFitter->fLimitTitle.find("95CL"),4,"95% CL");
    }
    // ---

    // Ser POITitle
    param = confSet->Get("POITitle");
    if( param != "") fMultiFitter->fPOITitle = RemoveQuotes(param);

    // Set CompareLimits
    param = confSet->Get("CompareLimits");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "FALSE") fMultiFitter->fCompareLimits = false;
        else if (param == "TRUE") fMultiFitter->fCompareLimits = true;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'CompareLimits' option but you didn't provide valid setting. Using default (TRUE)");
            fMultiFitter->fCompareLimits = true;
        }
    }

    // Ser ComparePOI
    param = confSet->Get("ComparePOI");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "FALSE") fMultiFitter->fComparePOI = false;
        else if (param == "TRUE") fMultiFitter->fComparePOI = true;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'ComparePOI' option but you didn't provide valid setting. Using default (TRUE)");
            fMultiFitter->fComparePOI = true;
        }
    }

    // Set ComparePulls
    param = confSet->Get("ComparePulls");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "FALSE") fMultiFitter->fComparePulls  = false;
        else if(param == "TRUE") fMultiFitter->fComparePulls  = true;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'ComparePulls' option but you didn't provide valid setting. Using default (TRUE)");
            fMultiFitter->fComparePulls  = true;
        }
    }

    // Set PlotCombCorrMatrix
    param = confSet->Get("PlotCombCorrMatrix");
    if (param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fPlotCombCorrMatrix  = true;
        else if(param == "FALSE") fMultiFitter->fPlotCombCorrMatrix  = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'PlotCombCorrMatrix' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fPlotCombCorrMatrix  = false;
        }
    }

    // Set CorrelationThreshold
    param = confSet->Get("CorrelationThreshold");
    if( param != ""){
        TRExFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }

    // Set UseGammasForCorr
    param = confSet->Get("UseGammasForCorr");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE"){
            fMultiFitter->fuseGammasForCorr = true;
        } else if (param == "FALSE") {
            fMultiFitter->fuseGammasForCorr = false;
        } else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified UseGammasForCorr option but did not provide valid parameter. Using default (false)");
            fMultiFitter->fuseGammasForCorr = false;
        }
    }

    // Set Combine
    param = confSet->Get("Combine");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fCombine = true;
        else if(param == "FALSE") fMultiFitter->fCombine = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'Combine' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fCombine = false;
        }
    }

    // Set Compare
    param = confSet->Get("Compare");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") fMultiFitter->fCompare = true;
        else if (param == "FALSE") fMultiFitter->fCompare = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'Compare' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fCompare = false;
        }
    }

    // Set StatOnly
    param = confSet->Get("StatOnly");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fStatOnly = true;
        else if(param == "FALSE") fMultiFitter->fStatOnly = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'StatOnly' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fStatOnly = false;
        }
    }

    // Set IncludeStatOnly
    param = confSet->Get("IncludeStatOnly");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fIncludeStatOnly = true;
        else if(param == "FALSE") fMultiFitter->fIncludeStatOnly = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'IncludeStatOnly' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fIncludeStatOnly = false;
        }
    }

    // Set POIName
    param = confSet->Get("POIName");
    if( param != "" ) fMultiFitter->fPOI = RemoveQuotes(param);

    // Set POILabel
    param = confSet->Get("POILabel");
    if( param != "" ) fMultiFitter->fPOIName = RemoveQuotes(param);

    // Set POINominal
    param = confSet->Get("POINominal");
    if( param != "" ) fMultiFitter->fPOINominal = std::stof(param);

    // Set POIRange
    // Set POIRange
    param = confSet->Get("POIRange");
    if( param != ""){
        if (Vectorize(param,',').size()==2 ) {
            fMultiFitter->fPOIMin = atof( Vectorize(param,',')[0].c_str() );
            fMultiFitter->fPOIMax = atof( Vectorize(param,',')[1].c_str() );
        } else {
            WriteErrorStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'POIRange' option but you didn't provide valid setting. Please check this!");
            return 1;
        }
    }

    // Set LimitMax
    param = confSet->Get("LimitMax");
    if( param != "" ) {
        fMultiFitter->fLimitMax = atof( param.c_str() );
    }

    // Set POIVal
    param = confSet->Get("POIVal");
    if( param != "" ) fMultiFitter->fPOIVal = atof(param.c_str());

    // Set POIPrecision
    param = confSet->Get("POIPrecision");
    if( param != "" ) fMultiFitter->fPOIPrecision = RemoveQuotes(param).c_str();

    //Set DataName
    param = confSet->Get("DataName");
    if( param != "" ) fMultiFitter->fDataName = RemoveQuotes(param);

    // Set FitType
    param = confSet->Get("FitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param=="SPLUSB")      fMultiFitter->fFitType = 1;
        else if(param=="BONLY")  fMultiFitter->fFitType = 2;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'FitType' option but you didn't provide valid setting. Using default (SPLUSB)");
            fMultiFitter->fFitType = 1;
        }
    }

    // Set CombineChByCh
    param = confSet->Get("CombineChByCh");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fCombineChByCh = true;
        else if(param == "FALSE") fMultiFitter->fCombineChByCh = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'CombineChByCh' option but you didn't provide valid setting. Using default (TRUE)");
            fMultiFitter->fCombineChByCh = true;
        }
    }

    // Set NPCategories
    param = confSet->Get("NPCategories");
    if( param != "" ) {
        std::vector<std::string> categ = Vectorize(param,',');
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
        TRExFitter::NCPU = atoi( param.c_str());
    }

    // Set FastFit
    param = confSet->Get("FastFit");
    if (param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" )  fMultiFitter->fFastFit = true;
        else if( param == "FALSE" )  fMultiFitter->fFastFit = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'FastFit' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fFastFit = false;
        }
    }

    // Set FastFitForRanking
    param = confSet->Get("FastFitForRanking");
    if (param != ""){
        if( param == "TRUE" )  fMultiFitter->fFastFitForRanking = true;
        else if( param == "FALSE" )  fMultiFitter->fFastFitForRanking = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'FastFitForRanking' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fFastFitForRanking = false;
        }
    }

    // Set NuisParListFile
    param = confSet->Get("NuisParListFile");
    if( param != "" ) fMultiFitter->fNuisParListFile = RemoveQuotes(param);

    // Set PlotSoverB
    param = confSet->Get("PlotSoverB");
    if (param != ""){
        if( param == "TRUE" ) fMultiFitter->fPlotSoverB = true;
        else if( param == "FALSE" ) fMultiFitter->fPlotSoverB = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'PlotSoverB' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fPlotSoverB = false;
        }
    }

    // Set SignalTitle
    param = confSet->Get("SignalTitle");
    if( param != "" ) fMultiFitter->fSignalTitle = RemoveQuotes(param);

    // Set FitResultsFile
    param = confSet->Get("FitResultsFile");
    if( param != "" ) fMultiFitter->fFitResultsFile = RemoveQuotes(param);

    // Set LimitsFile
    param = confSet->Get("LimitsFile");
    if( param != "" ) fMultiFitter->fLimitsFile = RemoveQuotes(param);

    // Set BonlySuffix
    param = confSet->Get("BonlySuffix");
    if( param != "" ) fMultiFitter->fBonlySuffix = RemoveQuotes(param);

    // Set ShowSystForPOI
    param = confSet->Get("ShowSystForPOI");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fShowSystForPOI = true;
        else if(param == "FALSE") fMultiFitter->fShowSystForPOI = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'ShowSystForPOI' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fShowSystForPOI = false;
        }
    }

    // Set GetGoodnessOfFit
    param = confSet->Get("GetGoodnessOfFit");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "TRUE") fMultiFitter->fGetGoodnessOfFit = true;
        else if(param == "FALSE") fMultiFitter->fGetGoodnessOfFit = false;
        else {
            WriteWarningStatus("ConfigReaderMulti::ReadJobOptions", "You specified 'GetGoodnessOfFit' option but you didn't provide valid setting. Using default (FALSE)");
            fMultiFitter->fGetGoodnessOfFit = false;
        }
    }

    // Set doLHscan
    param = confSet->Get("doLHscan");
    if( param != "" ){
        if (fOnlyLHscan==""){
            fMultiFitter->fVarNameLH = Vectorize(param,',');
        } else {
            fMultiFitter->fVarNameLH.emplace_back(fOnlyLHscan);
        }
    }

    // Set do2DLHscan
    param = confSet->Get("do2DLHscan");
    if( param != "" ){
        const std::vector<std::string> tmp = Vectorize(param,':');
        for (const auto& ivec : tmp){
            const std::vector<std::string> v = Vectorize(ivec,',');
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
            WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "You specified 'LHscanMinY' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanMinY = std::stof(param);
        }
    }

    // Set LHscanMax for second variable
    param = confSet->Get("LHscanMaxY");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "You specified 'LHscanMaxY' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanMaxY = std::stof(param);
        }
    }

    // Set LHscanSteps for second variable
    param = confSet->Get("LHscanStepsY");
    if ( param != "" ) {
        if (fMultiFitter->fVarName2DLH.size() == 0){
            WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "You specified 'LHscanStepsY' option but did not set doLHscan. Ignoring");
        } else {
            fMultiFitter->fLHscanStepsY = std::stoi(param);
            if(fMultiFitter->fLHscanStepsY < 3 || fMultiFitter->fLHscanStepsY > 100){
                WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "LHscanSteps is smaller than 3 or larger than 100, setting to defaut (30)");
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
        if( param == "TRUE" ){
            fMultiFitter->fParal2D = true;
        } else if (param == "FALSE") {
            fMultiFitter->fParal2D = false;
        } else {
            WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "You specified 'Parallel2Dscan' option but did not provide valid parameter. Using default (false)");
            fMultiFitter->fParal2D = false;
        }
    }

    // Set Paral2Dstep
    param = confSet->Get("Parallel2DscanStep");
    if ( param != "" ) {
        fMultiFitter->fParal2Dstep = std::atoi( param.c_str());
        if (fMultiFitter->fParal2Dstep < 1 || fMultiFitter->fParal2Dstep>=fMultiFitter->fLHscanSteps ){
            WriteErrorStatus("ConfigReaderMulti::ReadFitOptions", "You specified a step for 2D LHscan outside the allowed range.");
            return 1;
        }
    }

    // Set ShowTotalOnly
    param = confSet->Get("ShowTotalOnly");
    if ( param != "" ) {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "TRUE") {
            fMultiFitter->fShowTotalOnly = true;
        } else if (param == "FALSE"){
            fMultiFitter->fShowTotalOnly = false;
        }
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
        std::vector<std::string> vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TRExFitter::PREFITONPOSTFIT= true;
    }

    // Set POIInitial
    param = confSet->Get("POIInitia;");
    if( param != ""){
        fMultiFitter->fPOIInitial = std::stof(param);
    }

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReaderMulti::ReadLimitOptions(){
    std::string param = "";

    ConfigSet* confSet = fParser.GetConfigSet("Limit");
    if (confSet == nullptr){
        WriteDebugStatus("ConfigReader::ReadLimitOptions", "You do not have Limit option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set LimitBlind
    param = confSet->Get("LimitBlind");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fMultiFitter->fLimitIsBlind = true;
        } else if ( param == "FALSE" ){
            fMultiFitter->fLimitIsBlind = false;
        } else {
            WriteWarningStatus("ConfigReaderMulti::ReadLimitOptions", "You specified 'LimitBlind' option but did not provide valid parameter. Using default (false)");
            fMultiFitter->fLimitIsBlind = false;
        }
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
        fMultiFitter->fLimitPOIAsimov = atof(param.c_str());
    }

    // Set SignalInjection
    param = confSet->Get("SignalInjection");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fMultiFitter->fSignalInjection = true;
        } else if ( param == "FALSE" ){
            fMultiFitter->fSignalInjection = false;
        } else {
            WriteWarningStatus("ConfigReaderMulti::ReadLimitOptions", "You specified 'SignalInjection' option but did not provide valid parameter. Using default (false)");
            fMultiFitter->fSignalInjection = false;
        }
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

    return 0;
}

//__________________________________________________________________________________
//
int ConfigReaderMulti::ReadSignificanceOptions(){
    std::string param = "";

    ConfigSet* confSet = fParser.GetConfigSet("Significance");
    if (confSet == nullptr){
        WriteDebugStatus("ConfigReaderMulti::ReadSignificanceOptions", "You do not have Significance option in the config. It is ok, we just want to let you know.");
        return 0; // it is ok to not have Fit set up
    }

    // Set LimitBlind
    param = confSet->Get("SignificanceBlind");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fMultiFitter->fSignificanceIsBlind = true;
        } else if ( param == "FALSE" ){
            fMultiFitter->fSignificanceIsBlind = false;
        } else {
            WriteWarningStatus("ConfigReaderMulti::ReadSignificanceOptions", "You specified 'SignificanceBlind' option but did not provide valid parameter. Using default (false)");
            fMultiFitter->fSignificanceIsBlind = false;
        }
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

    return 0;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadFitOptions(const std::string& opt, const std::string& options){
    int nFit = 0;

    while(true){
        ConfigSet *confSet = fParser.GetConfigSet("Fit", nFit);
        if (confSet == nullptr) break;
        nFit++;

        // Set Options
        std::string fullOptions;
        std::string param = confSet->Get("Options");
        if(param!="" && options!="") fullOptions = options+";"+RemoveQuotes(param);
        else if(param!="") fullOptions = RemoveQuotes(param);
        else fullOptions = options;

        // name
        fMultiFitter->fFitNames.push_back(CheckName(confSet->GetValue()));

        // Set Label
        param = confSet->Get("Label");
        std::string label = CheckName(confSet->GetValue());
        if(param!="") label = RemoveQuotes(param);

        // Set suf
        param = confSet->Get("LoadSuf");
        std::string loadSuf = "";
        if(param!="") loadSuf = RemoveQuotes(param);
        else          loadSuf = fGlobalSuffix;

        // config file
        std::string confFile = "";
        param = confSet->Get("ConfigFile");
        if(param!="") confFile = RemoveQuotes(param);

        // workspace
        std::string wsFile = "";
        param = confSet->Get("Workspace");
        if(param!="") wsFile = RemoveQuotes(param);

        // show obs
        param = confSet->Get("ShowObserved");
        if (param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") fMultiFitter->fFitShowObserved.push_back(false);
            else if (param == "TRUE") fMultiFitter->fFitShowObserved.push_back(true);
            else {
                WriteWarningStatus("ConfigReaderMulti::ReadFitOptions", "You specified 'ShowObserved' option but you didn't provide valid setting. Using default (TRUE)");
                fMultiFitter->fFitShowObserved.push_back(true);
            }
        }

        fMultiFitter->AddFitFromConfig(confFile,opt,fullOptions,label,loadSuf,wsFile);

        // Set FitResultsFile
        param = confSet->Get("FitResultsFile");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fFitResultsFile = RemoveQuotes(param);
        fMultiFitter->fLimitsFiles.push_back("");

        // Set LimitsFile
        param = confSet->Get("LimitsFile");
        if( param != "" ) fMultiFitter->fLimitsFiles[fMultiFitter->fFitList.size()-1] = RemoveQuotes(param);

        // Set POIName
        param = confSet->Get("POIName");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fPOI = RemoveQuotes(param);

        // Set Directory
        param = confSet->Get("Directory");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fName = RemoveQuotes(param);

        // Set InputName
        param = confSet->Get("InputName");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fInputName = RemoveQuotes(param);
    }

    if (nFit == 0){
        WriteErrorStatus("ConfigReaderMulti::ReadFitOptions", "You need to provide at least one 'Fit' option. Please check this!");
        return 1;
    }

    return 0;
}

//__________________________________________________________________________________
//
std::string ConfigReaderMulti::CheckName( const std::string &name ){
    if( std::isdigit( name.at(0) ) ){
        WriteErrorStatus("ConfigReaderMulti::CheckName", "Failed to browse name: " + name + ". A number has been detected at the first position of the name.");
        WriteErrorStatus("ConfigReaderMulti::CheckName", "           This can lead to unexpected behaviour in HistFactory. Please change the name. ");
        WriteErrorStatus("ConfigReaderMulti::CheckName", "           The code is about to crash.");
        std::abort();
    } else {
        return RemoveQuotes(name);
    }
}
