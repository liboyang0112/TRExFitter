#ifndef CONFIGREADER_H
#define CONFIGREADER_H

// TtHFitter class includes
#include "TtHFitter/ConfigParser.h"

// standard c++ includes
#include <string>

//Forward class declaration
class TtHFit;
class Region;

/**
 * \class ConfigReader
 * \brief Class that reads config parameters
 * \author Tomas Dado
 */

class ConfigReader {
    public:

        /**
          * The default constructor.
          * @param fitter A pointer to TtHFit class
          */
        ConfigReader(TtHFit * fitter);

        /**
          * The default destructor
          */ 
        ~ConfigReader();

        /**
          * Reads the config and passes parameters to TtHFit
          * @param string Config path 
          * @param string Additional options
          * @return int status code
          */ 
        int ReadFullConfig(const std::string& fileName, const std::string& option);

    private:

        /**
          * Helper function to read Command line settings
          * @param string config options
          * @return int status code
          */
        int ReadCommandLineOptions(std::string option);
        
        /**
          * Helper function to read JOB settings
          * @return int status code
          */
        int ReadJobOptions();

        /**
          * Helper function to read plotting settings in JOB
          * @param ConfigSet A pointer needed to parse the config
          * @return int status code
          */
        int SetJobPlot(ConfigSet *confSet);

        /**
          * Helper function to read plotting settings in Options
          * @return int status code
          */
        int ReadGeneralOptions();

        /**
          * Helper function to read Fit settings
          * @return int status code
          */
        int ReadFitOptions();

        /**
          * Helper function to read Limit settings
          * @return int status code
          */
        int ReadLimitOptions();

        /**
          * Helper function to read Region settings
          * @return int status code
          */
        int ReadRegionOptions();

        /**
          * Helper function to read Region settings based on input type
          * @param Region A pointer needed to parse the input
          * @param ConfigSet A pointer neede to parse the input
          * @return int status code
          */
        int SetRegionHIST(Region* reg, ConfigSet *confSet);

        /**
          * Helper function to read Region settings based on input type
          * @param Region A pointer needed to parse the input
          * @param ConfigSet A pointer neede to parse the input
          * @return int status code
          */
        int SetRegionNTUP(Region* reg, ConfigSet *confSet);

        /**
          * Helper function to check if config has settings for NTUP
          * @param ConfigSet A pointer needed to parse the input
          * @return bool 
          */ 
        bool ConfigHasNTUP(ConfigSet* confSet);

        /**
          * Helper function to check if config has settings for HIST
          * @param ConfigSet A pointer needed to parse the input
          * @return bool 
          */ 
        bool ConfigHasHIST(ConfigSet* confSet);

        /**
          * Helper function to read Sample settings
          * @return int status code
          */
        int ReadSampleOptions();

        /**
          * Helper function to check the consistency of the input
          * @param string Input parameter
          * @return string Corrected parameter
          */
        std::string CheckName(const std::string &name);

        /**
          * Pointer to TtHFit class, set during initialization
          */
        TtHFit *fFitter;
    
        /**
          * Instance of ConfigParser used to parse the text
          */
        ConfigParser fParser;
   
        /**
          * vector of strings, one for each region
          */ 
        std::vector< std::string > fOnlyRegions;
        
        /**
          * vector of strings, one for each sample
          */ 
        std::vector< std::string > fOnlySamples;
        
        /**
          * vector of strings, one for each systematics
          */ 
        std::vector< std::string > fOnlySystematics;
        
        /**
          * vector of strings, one for each exclude region
          */ 
        std::vector< std::string > fToExclude;
       
        /**
          * vector of names, one for each region
          */ 
        std::vector<std::string> fRegNames; 

        /**
          *  string for signal only
          */ 
        std::string fOnlySignal = "";
};

#endif
