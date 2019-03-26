#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

/// c++ includes
#include <string>
#include <vector>

const static int MAXconfig = 2000;

// Functions
std::string RemoveSpaces(const std::string& s);
std::string RemoveComments(const std::string& s);
std::string RemoveQuotes(const std::string& s);
std::string Fix(const std::string& s);
std::vector<std::string> Vectorize(const std::string& s,char c,bool removeQuotes=true);
std::string First(const std::string& s);
std::string Second(const std::string& s);
std::string ReadValueFromConfig(const std::string& fileName,const std::string& option);

// classes
class Config {
public:
    Config();
    ~Config();

    std::string fName;
    std::string fValue;
};

class ConfigSet {
public:
    ConfigSet();
    ~ConfigSet();

    void SetConfig(const std::string& name,const std::string& value);
    std::string Get(const std::string& name) const;
    std::string operator[](const std::string& name) const;
    Config GetConfig(int i) const;
    std::string GetConfigName(int i) const;
    std::string GetConfigValue(int i) const;
    int GetN() const;
    int size() const;
    void Set(const std::string& name,const std::string& value);
    std::string GetName() const;
    std::string GetValue() const;

    int fN;
    std::string fName;
    std::string fValue;
    Config fConfig[MAXconfig];
};


class ConfigParser {
public:
    ConfigParser();
    ~ConfigParser();

    std::vector<ConfigSet*> fConfSets;
    void ReadFile(const std::string& fileName);
    ConfigSet *GetConfigSet(int i=0);
    ConfigSet *GetConfigSet(const std::string& name,int i=0);
    int CheckSyntax(ConfigParser *refConfigParser);

    /**
      * Method that checks provided setting with reference config file
      * @param ConfigSet Pointer to actual configuration setting type 
      * @param ConfigSet Pointer to referece config parser 
      * @param string Name of the setting set
      * @param string Name of the setting
      * @return int Status code
      */
    int SettingIsValid(ConfigSet *cs, ConfigParser *refConfigParser, const std::string &setting_set, const std::string &setting) const;

    /**
      * Helper method that checks single setting with reference file
      * @param ConfigSet Pointer to actual configuration setting type 
      * @param ConfigSet Pointer to referece config parser 
      * @param string Name of the setting set
      * @param string Name of the setting
      * @return int Status code
      */
    int CheckSingleSetting(ConfigSet *cs, ConfigSet *cs_ref, const std::string &setting_set, const std::string &setting) const;

    /**
      * Helper method that checks validity of provided parameters for given setting with reference file
      * @param string Setting in users config
      * @param string vector List of possible settings
      * @param string Name of the setting set
      * @param string Name of the setting
      * @return int Status code
      */
    int CheckParameters(std::string current, const std::vector<std::string> &possible_settings, const std::string& setting_set, const std::string &setting) const;

    /**
      * @param string Name of the setting set
      * @param string Name of the setting
      * @param string Setting value in users config
      * @param string Possible setting from reference file
      * @param char Delimiter whn multiple parameters are provided
      * @return bool Is valid setting
      */
    bool SettingMultipleParamIsOK(const std::string &setting_set, const std::string& setting, const std::string& current, const std::string& possible, const char delimiter = ',') const;

// private:
    int fN;
};

#endif
