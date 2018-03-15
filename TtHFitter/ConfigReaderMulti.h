#ifndef CONFIGREADERMULTI_H
#define CONFIGREADERMULTI_H

// TtHFitter class includes
#include "TtHFitter/ConfigParser.h"

// standard c++ includes
#include <string>

//Forward class declaration
class MultiFit;

/**
 * \class ConfigReaderMulti
 * \brief Class that reads config parameters
 * \author Tomas Dado
 */

class ConfigReaderMulti {
    public:

        /**
          * The default constructor.
          * @param MultiFIt A pointer to MultiFit class
          */
        ConfigReaderMulti(MultiFit * multiFitter);

        /**
          * The default destructor
          */ 
        ~ConfigReaderMulti();

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
        int ReadCommandLineOptions(const std::string &option);
        
        /**
          * Helper function to read JOB settings
          * @return int status code
          */
        int ReadJobOptions();

        /**
          * Helper function to read Fit settings
          * @return int status code
          */
        int ReadFitOptions();

        /**
          * Pointer to MultiFit class, set during initialization
          */
        MultiFit *fMultiFitter;
    
        /**
          * Instance of ConfigParser used to parse the text
          */
        ConfigParser fParser;

        /**
          * String Used for global settings
          */
        std::string fGlobalSuffix;
};

#endif
