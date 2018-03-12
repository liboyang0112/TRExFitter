#include "TtHFitter/ConfigReader.h"
#include "TtHFitter/TtHFit.h"


ConfigReader::ConfigReader(TtHFit *fitter){
    fFitter = fitter;
}

ConfigReader::~ConfigReader(){
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
    std::map< std::string,std::string > optMap;
    std::vector< std::string > optVec;
    std::vector< std::string > onlyRegions;
    std::vector< std::string > onlySamples;
    std::vector< std::string > onlySystematics;
    std::vector< std::string > toExclude;
    std::string onlySignal = "";*/

    // read different types of settings
    if (option != ""){
        sc+= ReadCommandLineOptions();
    }

    sc+= ReadJobOptions();


    return 0;
}

int ConfigReader::ReadCommandLineOptions(){
    return 0;
}

int ConfigReader::ReadJobOptions(){
    fConfSet = *fParser.GetConfigSet("Job");
    
    return 0;
}
