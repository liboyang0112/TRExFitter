#include "TtHFitter/ConfigReader.h"

#include "TtHFitter/TtHFit.h"
#include "TtHFitter/StatusLogbook.h"


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
    /*std::string param; // helper string
    std::vector< string > vec;
    int type;
    //
    // Read options (to skip stuff, or include only some regions, samples, systs...)
    // Syntax: .. .. Regions=ge4jge2b:Exclude=singleTop,wjets
    // */

    // read different types of settings
    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc+= ReadJobOptions();

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
    fConfSet = *fParser.GetConfigSet("Job");
    
    return 0;
}
