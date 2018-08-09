#ifndef CONFIGREADER_H
#define CONFIGREADER_H

// Framework includes
#include "TRExFitter/ConfigParser.h"

/// c++ includes
#include <string>
#include <vector>

///Forward class declaration
class ConfigSet;
class TRExFit;
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
          * @param fitter A pointer to TRExFit class
          */
        ConfigReader(TRExFit * fitter);

        /**
          * The default destructor
          */ 
        ~ConfigReader();

        /**
          * Reads the config and passes parameters to TRExFit
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
        int ReadCommandLineOptions(const std::string& option);
        
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
          * Helper function to read Significance settings
          * @return int status code
          */
        int ReadSignificanceOptions();

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
          * Helper function to check if elements of one vector are present in another
          * @param vector of parameters to check
          * @param vector of paramaeters to check to
          * @return True if all exist, False if at least one does not exist
          */            
        bool CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2);

        /**
          * Helper function to check if elements of one vector are present in another
          * @param vector of parameters to check
          * @param vector of paramaeters to check to
          * @param vector of paramaeters to check to
          * @return True if all exist, False if at least one does not exist
          */            
        bool CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2, const std::vector<std::string> &v3);

        /**
          * Pointer to TRExFit class, set during initialization
          */
        TRExFit *fFitter;
    
        /**
          * Instance of ConfigParser used to parse the text
          */
        ConfigParser fParser;

        /**
          * flag to control if wrong samples/regions are ok
          */
        bool fAllowWrongRegionSample;

        /**
          * flag to control if other than ghost sampels have been set already
          */
        bool fNonGhostIsSet; 

        /**
          * vector of strings, one for each sample, needed for cross-checks
          */
        std::vector< std::string > fSamples; 
   
        /**
          * vector of strings, one for each region, needed for cross-checks
          */
        std::vector< std::string > fRegions; 
   
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