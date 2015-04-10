#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef _CONFIGPARSER
#define _CONFIGPARSER

using namespace std;

const int MAXconfig = 200;

// Functions
string Fix(string s);
vector<string> Vectorize(string s,char c);
string First(string s);
string Second(string s);

// classes
class Config {
public:
  Config();
  ~Config();
  
  string fName;
  string fValue;
};

class ConfigSet {
public:
  ConfigSet();
  ~ConfigSet();
  
  void SetConfig(string name,string value);
  string Get(string name);
  string operator[](string name);
  Config GetConfig(int i);
  string GetConfigName(int i);
  string GetConfigValue(int i);
  int GetN();
  int size();
  void Set(string name,string value);
  string GetName();
  string GetValue();
  
private:
  int fN;
  string fName;
  string fValue;
  Config fConfig[MAXconfig];
};


class ConfigParser {
public:
  ConfigParser();
  ~ConfigParser();
  
  ConfigSet *fConfSets[MAXconfig];
  void ReadFile(string fileName);
  ConfigSet *GetConfigSet(int i=0);
  ConfigSet *GetConfigSet(string name,int i=0);
  
private:
  int fN;
};

#endif