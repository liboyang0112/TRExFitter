#ifndef CONFIGREADER_H
#define CONFIGREADER_H

// TtHFitter class includes
#include "TtHFitter/ConfigParser.h"

// standard c++ includes
#include <string>
#include <iostream>
#include <fstream>

//Forward class declaration
class TtHFit;
class Region;
class Systematic;

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
          * Helper function to read NormFactor settings
          * @return int status code
          */
        int ReadNormFactorOptions();

        /**
          * Helper function to read ShapeFactor settings
          * @return int status code
          */
        int ReadShapeFactorOptions();

        /**
          * Helper function to read Systematic settings
          * @return int status code
          */
        int ReadSystOptions();

        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @param vector of strings Needed for this setting
          * @param vector of strings Needed for this setting
          * @return int status code
          */
        int SetSystNoDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude);

        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @param vector of strings Needed for this setting
          * @param vector of strings Needed for this setting
          * @param vector of strings Needed for this setting
          * @param int Flag needed for for this setting
          * @return int status code
          */
        int SetSystRegionDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude, const std::vector<std::string> regions, int type);
    
        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @return int status code
          */
        int SetSystSampleDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude);
    
        /**
          * Helper function to read Part of Syst config
          * @param COnfigSet A pointer needed for reading
          * @param Systematic A pointer to syst that is being set
          * @return int status code
          */
        int SetSystShapeDecorelate(ConfigSet *confSet, Systematic *sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude);

        /**
          * Helper function that is run after config is read
          * @return int status code
          */
        int PostConfig();
    
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
          * vector of strings, one for each exclude region sample
          */ 
        std::vector< std::string > fExcludeRegionSample; 
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
