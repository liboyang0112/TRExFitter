#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include "TtHFitter/ConfigParser.h"

#include <string>

class TtHFit;

class ConfigReader {
    public:
       ConfigReader(TtHFit * fitter);

       ~ConfigReader();

       int ReadFullConfig(const std::string& fileName, const std::string& option);

    private:

        // private methods:
        int ReadCommandLineOptions();
        
        int ReadJobOptions();

        // private members
        TtHFit *fFitter;
    
        ConfigSet fConfSet;

        ConfigParser fParser;
};

#endif
