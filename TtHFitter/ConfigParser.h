#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

const int MAXconfig = 2000;

// Functions
std::string Fix(std::string s);
std::vector<std::string> Vectorize(std::string s,char c);
std::string First(std::string s);
std::string Second(std::string s);

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

    void SetConfig(std::string name,std::string value);
    std::string Get(std::string name);
    std::string operator[](std::string name);
    Config GetConfig(int i);
    std::string GetConfigName(int i);
    std::string GetConfigValue(int i);
    int GetN();
    int size();
    void Set(std::string name,std::string value);
    std::string GetName();
    std::string GetValue();

private:
    int fN;
    std::string fName;
    std::string fValue;
    Config fConfig[MAXconfig];
};


class ConfigParser {
public:
    ConfigParser();
    ~ConfigParser();

    ConfigSet *fConfSets[MAXconfig];
    void ReadFile(std::string fileName);
    ConfigSet *GetConfigSet(int i=0);
    ConfigSet *GetConfigSet(std::string name,int i=0);

private:
    int fN;
};

#endif
