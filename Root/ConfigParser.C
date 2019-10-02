// Class include
#include "TRExFitter/ConfigParser.h"

// Framework includes
#include "TRExFitter/StatusLogbook.h"

// c++ includes
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>

// using namespace std;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// -- Functions --
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
// Removes leading and trailing white spaces
std::string RemoveSpaces(const std::string& s){
    if(s=="") return "";
    std::string ss = s;
    if (ss.find_first_not_of(' ')>=std::string::npos){
        ss = "";
    }
    else if (ss.find_first_not_of(' ')>0){
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' '));
    }
    else{
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1);
    }
    return ss;
}

//__________________________________________________________________________________
// Removes everything after '%' or '#', but only if not inside quotation marks!!
std::string RemoveComments(const std::string& s){
    if(s=="") return "";
    std::string ss = "";
    bool insideQuotes = false;
    for(unsigned long i=0;i<s.size();i++){
        if(s[i]=='"'){
            if(!insideQuotes) insideQuotes = true;
            else              insideQuotes = false;
        }
        if((s[i]=='%' || s[i]=='#') && !insideQuotes) break;
        ss += s[i];
    }
    return RemoveSpaces(ss);
}

//__________________________________________________________________________________
// Removes '"'
std::string RemoveQuotes(const std::string& s){
    if(s=="") return "";
    std::string ss = s;
    replace( ss.begin(), ss.end(), '"', ' ');
    return RemoveSpaces(ss);
}

//__________________________________________________________________________________
// Removes leading and trailing white spaces, removes everything after "%", removes '"'
std::string Fix(const std::string& s){
    if(s=="") return "";
    std::string ss = s.substr( 0, s.find_first_of('%') );
    replace( ss.begin(), ss.end(), '"', ' ');
    if (ss.find_first_not_of(' ')>=std::string::npos){
        ss = "";
    }
    else if (ss.find_first_not_of(' ')>0){
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' '));
    }
    else{
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1);
    }
    return ss;
}

//__________________________________________________________________________________
//
std::vector<std::string> Vectorize(const std::string& s,char c,bool removeQuotes){
    std::vector<std::string> v;
    v.clear();
    std::string ss = RemoveComments(s);
    if(ss==""){
        v.push_back("");
        return v;
    }
    std::string t;
    bool insideQuotes = false;
    for(unsigned long i=0;i<ss.size();i++){
        if(!insideQuotes && ss[i]==c){
            if(removeQuotes) v.push_back(RemoveQuotes(t));
            else             v.push_back(RemoveSpaces(t));
            t = "";
        }
        else if(!insideQuotes && ss[i]=='"'){
            insideQuotes = true;
            t += ss[i];
        }
        else if(insideQuotes && ss[i]=='"'){
            insideQuotes = false;
            t += ss[i];
        }
        else{
            t += ss[i];
        }
    }
    if(RemoveQuotes(t)!=""){
        if(removeQuotes) v.push_back(RemoveQuotes(t));
        else             v.push_back(RemoveSpaces(t));
    }
    return v;
}

//__________________________________________________________________________________
//
std::string First(const std::string& s){
    std::string first;
    first = s.substr( 0, s.find_first_of(':') );
    return Fix(first);
}

//__________________________________________________________________________________
//
std::string Second(const std::string& s){
    std::string second;
    second = s.substr( s.find_first_of(':')+1,std::string::npos );
    second = RemoveComments(second);
    if(second==""){
        WriteErrorStatus("ConfigParser","No value set for parameter "+First(s)+" in the config.");
        exit(EXIT_FAILURE);
    }
    return second;
}

//_______________________________________________________________________________________
//
void ReplaceStringInPlace(std::string& subject, const std::string& search,
                                                const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
}

//_______________________________________________________________________________________
// used to pre-read a config file to check single option values
std::string ReadValueFromConfig(const std::string& fileName,const std::string& option){
    std::string value = "";
    std::ifstream file(fileName.c_str());
    if(!file.is_open()) return value;
    std::string str;
    while(true){
        file >> str;
        if(!file.good()) break;
        if(str==(option+":")){
            file >> value;
            break;
        }
    }
    file.close();
    file.clear();
    return value;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// -- Classes --
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
// Config
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
Config::Config(){}

//__________________________________________________________________________________
//
Config::~Config(){}

//----------------------------------------------------------------------------------
// ConfigSet
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
ConfigSet::ConfigSet(){
    fN = 0;
}

//__________________________________________________________________________________
//
ConfigSet::~ConfigSet(){}

//__________________________________________________________________________________
//
void ConfigSet::SetConfig(const std::string& name,const std::string& value){
    // check if is name is there
    bool isThere = false;
    int index = -1;
    for(int i=0;i<fN;i++){
        if(fConfig[i].fName == name){
            isThere = true;
            index = i;
        }
    }
    if(!isThere){
        fConfig[fN].fName = name;
        fConfig[fN].fValue = value;
        fN++;
    }
    else{
        fConfig[index].fName = name;
        fConfig[index].fValue = value;
    }
}

//__________________________________________________________________________________
//
std::string ConfigSet::Get(const std::string& name) const{
    for(int i=0;i<fN;i++){
        if(fConfig[i].fName == name){
            return fConfig[i].fValue;
        }
    }
    return "";
}

//__________________________________________________________________________________
//
std::string ConfigSet::operator[](const std::string& name) const{
    return Get(name);
}

//__________________________________________________________________________________
//
Config ConfigSet::GetConfig(int i) const{
    return fConfig[i];
}

//__________________________________________________________________________________
//
std::string ConfigSet::GetConfigName(int i) const{
    return fConfig[i].fName;
}

//__________________________________________________________________________________
//
std::string ConfigSet::GetConfigValue(int i) const{
    return fConfig[i].fValue;
}

//__________________________________________________________________________________
//
int ConfigSet::GetN() const{
    return fN;
}

//__________________________________________________________________________________
//
int ConfigSet::size() const{
    return fN;
}

//__________________________________________________________________________________
//
void ConfigSet::Set(const std::string& name,const std::string& value){
    fName = name;
    fValue = value;
}

//__________________________________________________________________________________
//
std::string ConfigSet::GetName() const{
    return fName;
}

//__________________________________________________________________________________
//
std::string ConfigSet::GetValue() const{
    return fValue;
}

//----------------------------------------------------------------------------------
// ConfigParser
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
ConfigParser::ConfigParser(){
    fN = 0;
    fConfSets.reserve(MAXconfig);
}

//__________________________________________________________________________________
//
ConfigParser::~ConfigParser(){
}

//__________________________________________________________________________________
//
void ConfigParser::ReadFile(const std::string& fileName){
    WriteInfoStatus("ConfigParser::ReadFile", "Reading config file: " + fileName);
    std::ifstream file(fileName.c_str());
    if(!file.is_open()){
        WriteErrorStatus("ConfigParser::ReadFile", "The config file cannot be opened!");
        exit(-1);
    }
    std::string str, val;
    //
    std::map<std::string,std::string> fReplacement;
    std::string replacementFileName="";
    // loop over the file once to automatically find the replacement fileName
    while (getline(file, str)){
        replace( str.begin(), str.end(), '\n', ' ');
        replace( str.begin(), str.end(), '\r', ' ');
        replace( str.begin(), str.end(), '\t', ' ');
        if(str[str.find_first_not_of(' ')]=='%') continue;
        if(str[str.find_first_not_of(' ')]=='#') continue;
        if ( str.find("ReplacementFile")==std::string::npos ) continue;
        replacementFileName=Second(str);
        break;
    }
    file.close();
    file.clear();
    if (replacementFileName!="" && (fileName.find("jobSchema.config") == std::string::npos)) {
        WriteInfoStatus("ConfigParser::ReadFile", "Opening replacement file: " + replacementFileName + " to fill the map");
        std::ifstream fileR(replacementFileName.c_str());
        if(!fileR.is_open()){
            WriteErrorStatus("ConfigParser::ReadFile", "The replacement file: " + replacementFileName + " cannot be opened!");
            exit(-1);
        }
        while (getline(fileR, str)){
            replace( str.begin(), str.end(), '\n', ' ');
            replace( str.begin(), str.end(), '\r', ' ');
            replace( str.begin(), str.end(), '\t', ' ');
            if ( str.find("%")!=std::string::npos || str.find("#")!=std::string::npos) continue;
            WriteInfoStatus("ConfigParser::ReadFile", "putting in the map: [" + First(str) + "]=" + Second(str));
            fReplacement[First(str)]=Second(str);
        }
        fileR.close();
    }
    //
    file.open(fileName.c_str());
    std::vector<std::string> valVec;
    bool reading = false;
    int n = 1;
    int k = 0;
    while (getline(file, str)){
        replace( str.begin(), str.end(), '\n', ' ');
        replace( str.begin(), str.end(), '\r', ' ');
        replace( str.begin(), str.end(), '\t', ' ');
        if(str[str.find_first_not_of(' ')]=='%') continue;
        if(str[str.find_first_not_of(' ')]=='#') continue;
        //
        if (str.find("XXX")!=std::string::npos) {
            WriteInfoStatus("ConfigParser::ReadFile", "BEFORE replacement: " + str);
            std::map<std::string,std::string>::iterator itr=fReplacement.begin();
            for ( ;itr!=fReplacement.end();++itr) {
                std::string oldV=itr->first;
                std::string newV=itr->second;
                ReplaceStringInPlace(str, oldV, newV);
            }
            WriteInfoStatus("ConfigParser::ReadFile", " AFTER replacement: " + str);
        }
        //
        if(str.find_first_not_of(' ')==std::string::npos){
            reading = false;
        }
        else{
            valVec = Vectorize(Second(str),';',false);
            if(!reading){
                n = valVec.size();
                for(k=0;k<n;k++){
                    fConfSets.emplace_back(std::make_unique<ConfigSet>());
                    fConfSets[fN]->Set( First(str),RemoveSpaces(valVec[k]) );
                    fN++;
                }
                reading = true;
            }
            else{
                for(k=0;k<n;k++){
                    if(k>=(int)valVec.size()) val = valVec[valVec.size()-1];
                    else                      val = valVec[k];
                    fConfSets[fN-n+k]->SetConfig(First(str),RemoveSpaces(val) );
                }
            }
        }
    }
    file.close();
}

//__________________________________________________________________________________
// Returns the i-th configSet
ConfigSet *ConfigParser::GetConfigSet(int i){
    return fConfSets[i].get();
}

//__________________________________________________________________________________
// Returns the i-th configSet with given name
ConfigSet *ConfigParser::GetConfigSet(const std::string& name,int i){
    int k = 0;
    for(int j=0;j<fN;j++){
        if(fConfSets[j]->GetName() == name){
            if(i==k) return fConfSets[j].get();
            k++;
        }
    }
    return nullptr;
}

//__________________________________________________________________________________
//
int ConfigParser::CheckSyntax(ConfigParser *refConfigParser){
    int exitStatus = 0;
    // loop on all the confic sets
    for(int i_cs = 0; i_cs<fN; i_cs++){
        ConfigSet *cs = fConfSets[i_cs].get();
        // check if the same exists in the reference
        for (int i_c = 0; i_c < cs->fN; i_c++){
            Config c = cs->fConfig[i_c];
            exitStatus+= SettingIsValid(cs, refConfigParser, cs->fName, c.fName);
        }
    }
    return exitStatus;
}

//_______________________________________________________________________________________
//
int ConfigParser::SettingIsValid(ConfigSet *cs, ConfigParser *refConfigParser, const std::string &setting_set, const std::string &setting) const{
    if (refConfigParser == nullptr){
        WriteErrorStatus("ConfigParser::SettingIsPresentAndValid", "Invalid pointer to the reference ConfigParser. Please check this!");
        exit(EXIT_FAILURE);
    }
    ConfigSet *cs_ref = nullptr;
    bool refIsFound = false;
    // check if the setting type exists in the reference
    for (int i_cs =0; i_cs < refConfigParser->fN; i_cs++){
        cs_ref = refConfigParser->fConfSets[i_cs].get();
        if (cs_ref->fName == setting_set){
            refIsFound = true;
            break;
        }
    }
    //
    // config set is not present in the reference
    if (!refIsFound){
        WriteErrorStatus("ConfigParser::SettingIsValid", "Cannot find config set '" + setting_set + "' in reference config. Please check this!");
        return 1;
    }
    //
    // check the validity of the single setting
    if(CheckSingleSetting(cs, cs_ref, setting_set, setting)) return 1;
    return 0;
}

//_______________________________________________________________________________________
//
int ConfigParser::CheckSingleSetting(ConfigSet *cs, ConfigSet *cs_ref, const std::string &setting_set, const std::string &setting) const {
    std::string param = cs->Get(setting);
    std::string ref_param = cs_ref->Get(setting);
    //
    // this option used to exist, write a message informing the user of the new syntax
    if (setting == "TtresSmoothing"){
        WriteErrorStatus("ConfigParser::CheckSingleSetting", "The TtresSmoothing option is deprecated, use SmoothingOption: TTBARRESONANCE instead.");
        return 1;
    }
    else if (ref_param == ""){
        WriteErrorStatus("ConfigParser::CheckSingleSetting", "Cannot find setting '" + setting + "' for setting set " + setting_set +  " in reference config. Please check this!");
        return 1;
    }
    //
    // there is nothing to check if the reference setting is simply 'std::string'
    if (ref_param == "string") return 0;
    //
    // need to check the consistency of the provided settings
    // first check the number of provided parameters
    std::vector<std::string> current_settings = Vectorize(param,',');
    std::vector<std::string> possible_settings = Vectorize(ref_param,'/');
    std::vector<unsigned int> possible_sizes;
    //
    for (const std::string &ioption : possible_settings){
        possible_sizes.push_back( Vectorize(ioption,',').size() );
    }
    unsigned int current_setting_size = current_settings.size();
    //
    // search the vector for the possible sizes
    auto it = std::find(possible_sizes.begin(), possible_sizes.end(), current_setting_size);
    if (it == possible_sizes.end()){
        WriteErrorStatus("ConfigParser::CheckSingleSetting", "You provided " + std::to_string(current_setting_size) +" parameters, for setting set '" + setting_set + "' and setting '" + setting );
        std::string tmp = "Possible setting sizes are: ";
        for (const unsigned int& i : possible_sizes){
            tmp+= std::to_string(i) + " ";
        }
        tmp+= ". Please check this!";
        WriteErrorStatus("ConfigParser::CheckSingleSetting", tmp );
        return 1;
    }
    //
    // sizes are correct, not we need to check the actual input
    if(CheckParameters(param, possible_settings, setting_set, setting)) return 1;
    return 0;
}

//_______________________________________________________________________________________
//
int ConfigParser::CheckParameters(std::string current, const std::vector<std::string> &possible_settings, const std::string &setting_set, const std::string &setting) const{
    if (possible_settings.size() == 1){
        if(!SettingMultipleParamIsOK(setting_set, setting, current, possible_settings.at(0))) return 1;
    }
    else {
        // multiple settings are possible
        bool isFound = false;
        std::transform(current.begin(), current.end(), current.begin(), ::toupper);
        for (const std::string& isetting : possible_settings){
            // for std::string
            if(SettingMultipleParamIsOK(setting_set, setting, current, isetting)){
                isFound = true;
                break;
            }
        }
        if (!isFound){
            WriteErrorStatus("ConfigParser::CheckParameters", "Parameter " + current +" is not valid, for setting set '" + setting_set + "' and setting '" + setting );
            std::string tmp = "Possible values: ";
            for (const std::string &i : possible_settings){
                tmp+= i + " ";
            }
            tmp+= ". Please check this!";
            WriteErrorStatus("ConfigParser::CheckParameters", tmp );
            return 1;
        }
    }
    return 0;
}

//_______________________________________________________________________________________
//
bool ConfigParser::SettingMultipleParamIsOK(const std::string& setting_set, const std::string& setting, const std::string& current, const std::string& possible, const char delimiter) const{
    std::vector<std::string> current_vec = Vectorize(current, delimiter);
    std::vector<std::string> possible_vec = Vectorize(possible, delimiter);
    if (current_vec.size() != possible_vec.size()) return false;
    //
    // check setting by setting
    for (unsigned int iparam = 0; iparam < possible_vec.size(); iparam++){
        if (possible_vec.at(iparam) == "string"){
            continue; // nothing to check
        } else if (possible_vec.at(iparam) == "int"){
            try {
                std::stoi(current_vec.at(iparam));
            } catch (std::exception &e){
                WriteErrorStatus("ConfigParser::SettingMultipleParamIsOK", "Parameter " + current_vec.at(iparam) + " is not valid, for setting set '" + setting_set + "' and setting '" + setting + ", for parameter number " + std::to_string(iparam+1) + ". Please check this!" );
                return false;
            }
        } else if (possible_vec.at(iparam) == "float"){
            try {
                std::stod(current_vec.at(iparam));
            } catch (std::exception &e){
                WriteErrorStatus("ConfigParser::SettingMultipleParamIsOK", "Parameter " + current_vec.at(iparam) + " is not valid, for setting set '" + setting_set + "' and setting '" + setting + ", for parameter number " + std::to_string(iparam+1) + ". Please check this!" );
                return false;
            }
        } else {
            if (current_vec.at(iparam) != possible_vec.at(iparam)){
                return false;
            }
        }
    }
    return true;
}
