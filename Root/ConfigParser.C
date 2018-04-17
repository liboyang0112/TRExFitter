#include "TtHFitter/ConfigParser.h"
#include "TtHFitter/StatusLogbook.h"
#include <map>
#include <exception>

using namespace std;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// -- Functions --
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
string Fix(string s){ // removes leading and trailing white spaces, remove everything after "%", remove '"'...
    string ss = s.substr( 0, s.find_first_of('%') );
    replace( ss.begin(), ss.end(), '"', ' ');
//   if(ss.find("#PERCENT#")!=string::npos) ss.replace( ss.find("#PERCENT#"), 9, "%");
//   ss = ss.substr( ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1 );
    // Fix by Arthur
    if (ss.find_first_not_of(' ')>0){
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' '));
    }
    else{
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1);
    }
    return ss;
}

//__________________________________________________________________________________
//
vector<string> Vectorize(string s,char c){
    vector<string> v;
    v.clear();
    if(s==""){
        v.push_back("");
        return v;
    }
    string t;
    while(true){
        t = s.substr(0,s.find_first_of(c));
        v.push_back(Fix(t));
        if(t==s) break;
        s = s.substr(s.find_first_of(c)+1,string::npos);
    }
    return v;
}

//__________________________________________________________________________________
//
string First(string s){
    string first;
    first = s.substr( 0, s.find_first_of(':') );
    return Fix(first);
}

//__________________________________________________________________________________
//
string Second(string s){
    string second;
    second = s.substr( s.find_first_of(':')+1,string::npos );
    return Fix(second);
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
std::string ReadValueFromConfig(std::string fileName,std::string option){
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
void ConfigSet::SetConfig(string name,string value){
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
string ConfigSet::Get(string name){
    for(int i=0;i<fN;i++){
        if(fConfig[i].fName == name){
            return fConfig[i].fValue;
        }
    }
    return "";
}

//__________________________________________________________________________________
//
string ConfigSet::operator[](string name){
    return Get(name);
}

//__________________________________________________________________________________
//
Config ConfigSet::GetConfig(int i){
    return fConfig[i];
}

//__________________________________________________________________________________
//
string ConfigSet::GetConfigName(int i){
    return fConfig[i].fName;
}

//__________________________________________________________________________________
//
string ConfigSet::GetConfigValue(int i){
    return fConfig[i].fValue;
}

//__________________________________________________________________________________
//
int ConfigSet::GetN(){
    return fN;
}

//__________________________________________________________________________________
//
int ConfigSet::size(){
    return fN;
}

//__________________________________________________________________________________
//
void ConfigSet::Set(string name,string value){
    fName = name;
    fValue = value;
}

//__________________________________________________________________________________
//
string ConfigSet::GetName(){
    return fName;
}

//__________________________________________________________________________________
//
string ConfigSet::GetValue(){
    return fValue;
}

//----------------------------------------------------------------------------------
// ConfigParser
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
ConfigParser::ConfigParser(){
    fN = 0;
}

//__________________________________________________________________________________
//
ConfigParser::~ConfigParser(){
}

//__________________________________________________________________________________
//
void ConfigParser::ReadFile(string fileName){
    WriteInfoStatus("ConfigParser::ReadFile", "Reading config file: " + fileName);
    ifstream file(fileName.c_str());
    if(!file.is_open()){
        WriteErrorStatus("ConfigParser::ReadFile", "The config file cannot be opened!");
        exit(-1);
    }
    string str, val;

    map<string,string> fReplacement;
    string replacementFileName="";
    //loop over the file once to automatically find the replacement fileName
    while (getline(file, str)){
        replace( str.begin(), str.end(), '\n', ' ');
        replace( str.begin(), str.end(), '\r', ' ');
        replace( str.begin(), str.end(), '\t', ' ');
        if ( str.find("%")!=string::npos || str.find("#")!=string::npos) continue;
        if ( str.find("ReplacementFile")==string::npos ) continue;
        replacementFileName=Second(str);
        break;
    }
    file.close();
    file.clear();
    if (replacementFileName!="" && (fileName.find("jobSchema.config") == string::npos)) {
        WriteInfoStatus("ConfigParser::ReadFile", "Opening replacement file: " + replacementFileName + " to fill the map");
        ////replacementFileName="Common_XS_unc_Replacement.txt";
        ifstream fileR(replacementFileName.c_str());
        if(!fileR.is_open()){
            WriteErrorStatus("ConfigParser::ReadFile", "The replacement file: " + replacementFileName + " cannot be opend!");
            exit(-1);
        }
        while (getline(fileR, str)){
            replace( str.begin(), str.end(), '\n', ' ');
            replace( str.begin(), str.end(), '\r', ' ');
            replace( str.begin(), str.end(), '\t', ' ');
            if ( str.find("%")!=string::npos || str.find("#")!=string::npos) continue;
            WriteInfoStatus("ConfigParser::ReadFile", "putting in the map: [" + First(str) + "]=" + Second(str));
            fReplacement[First(str)]=Second(str);
        }
        fileR.close();
    }

    file.open(fileName.c_str());
    vector<string> valVec;
    bool reading = false;
    int n = 1;
    int k = 0;
    while (getline(file, str)){
        replace( str.begin(), str.end(), '\n', ' ');
        replace( str.begin(), str.end(), '\r', ' ');
//       if(str[0]=='%') continue;
        if(str[str.find_first_not_of(' ')]=='%') continue;

        if (str.find("XXX")!=string::npos) {
            WriteInfoStatus("ConfigParser::ReadFile", "BEFORE replacement: " + str);
            map<string,string>::iterator itr=fReplacement.begin();
            for ( ;itr!=fReplacement.end();++itr) {
                string oldV=itr->first;
                string newV=itr->second;
                //std::cout << "VALERIO: " << oldV << " , " << newV << std::endl;
                //replace( str.begin(), str.end(), oldV, newV);
                ReplaceStringInPlace(str, oldV, newV);
            }
            WriteInfoStatus("ConfigParser::ReadFile", " AFTER replacement: " + str);
        }

        if(str.find_first_not_of(' ')==string::npos){
            reading = false;
        }
        else{
            valVec = Vectorize(Second(str),';');
            if(!reading){
                n = valVec.size();
                for(k=0;k<n;k++){
                    fConfSets[fN] = new ConfigSet();
                    fConfSets[fN]->Set( First(str),Fix(valVec[k]) );
                    fN++;
                }
                reading = true;
            }
            else{
                for(k=0;k<n;k++){
                    if(k>=(int)valVec.size()) val = valVec[valVec.size()-1];
                    else                                            val = valVec[k];
                    fConfSets[fN-n+k]->SetConfig(First(str),Fix(val) );
                }
            }
        }
    }
    file.close();
}

//__________________________________________________________________________________
//
ConfigSet *ConfigParser::GetConfigSet(int i){ // returns the i-th configSet
    return fConfSets[i];
}

//__________________________________________________________________________________
//
ConfigSet *ConfigParser::GetConfigSet(string name,int i){ // returns the i-th configSet with given name
    int k = 0;
    for(int j=0;j<fN;j++){
        if(fConfSets[j]->GetName() == name){
            if(i==k) return fConfSets[j];
            k++;
        }
    }
    return 0x0;
}

//__________________________________________________________________________________
//
int ConfigParser::CheckSyntax(ConfigParser *refConfigParser){
    int exitStatus = 0;
    // loop on all the confic sets
    for(int i_cs = 0; i_cs<fN; i_cs++){
        ConfigSet *cs = fConfSets[i_cs];
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
        cs_ref = refConfigParser->fConfSets[i_cs];
        if (cs_ref->fName == setting_set){
            refIsFound = true;
            break;
        }
    }

    // config set is not present in the reference
    if (!refIsFound){
        WriteErrorStatus("ConfigParser::SettingIsValid", "Cannot find config set '" + setting_set + "' in reference config. Please check this!");
        return 1;
    }

    // check the validity of the single setting
    if(CheckSingleSetting(cs, cs_ref, setting_set, setting)) return 1;

    return 0;
}

//_______________________________________________________________________________________
//
int ConfigParser::CheckSingleSetting(ConfigSet *cs, ConfigSet *cs_ref, const std::string &setting_set, const std::string &setting) const {
    std::string param = cs->Get(setting);
    std::string ref_param = cs_ref->Get(setting);

    if (ref_param == ""){
        WriteErrorStatus("ConfigParser::CheckSingleSetting", "Cannot find setting '" + setting + "' for setting set " + setting_set +  " in reference config. Please check this!");
        return 1;
    }

    // there is nothing to check if the reference setting is simply 'string'
    if (ref_param == "string") return 0;

    // need to check the consistency of the provided settings
    // first check the number of provided parameters
    std::vector<std::string> current_settings = Vectorize(param,',');

    std::vector<std::string> possible_settings = Vectorize(ref_param,'/');
    std::vector<unsigned int> possible_sizes;

    for (const std::string &ioption : possible_settings){
        possible_sizes.push_back( Vectorize(ioption,',').size() );
    }

    unsigned int current_setting_size = current_settings.size();
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

    // sizes are correct, not we need to check the actual input
    if(CheckParameters(param, possible_settings, setting_set, setting)) return 1;

    return 0;
}

//_______________________________________________________________________________________
//
int ConfigParser::CheckParameters(std::string current, const std::vector<std::string> &possible_settings, const std::string &setting_set, const std::string &setting) const{
    if (possible_settings.size() == 1){
        if(!SettingMultipleParamIsOK(setting_set, current, possible_settings.at(0))) return 1;
    }
    else {
        // multiple settings are possible
        bool isFound = false;
        std::transform(current.begin(), current.end(), current.begin(), ::toupper);
        for (const std::string& isetting : possible_settings){
            // for string
            if(SettingMultipleParamIsOK(setting_set, current, isetting)){
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
bool ConfigParser::SettingMultipleParamIsOK(const std::string& setting_set, const std::string& current, const std::string& possible, const char delimiter) const{

    std::vector<std::string> current_vec = Vectorize(current, delimiter);
    std::vector<std::string> possible_vec = Vectorize(possible, delimiter);
    
    if (current_vec.size() != possible_vec.size()) return false;

    // check setting by setting
    for (unsigned int iparam = 0; iparam < possible_vec.size(); iparam++){
        if (possible_vec.at(iparam) == "string"){
            continue; // nothing to check
        } else if (possible_vec.at(iparam) == "int"){
            try {
                std::stoi(current_vec.at(iparam));
            } catch (std::exception &e){
                WriteErrorStatus("ConfigParser::SettingMultipleParamIsOK", "Parameter " + current_vec.at(iparam) + " is not valid, for setting set '" + setting_set + "' and setting '" + current + ", for parameter number " + std::to_string(iparam+1) + ". Please check this!" );
                return false;
            }
        } else if (possible_vec.at(iparam) == "float"){
            try {
                std::stof(current_vec.at(iparam));
            } catch (std::exception &e){
                WriteErrorStatus("ConfigParser::SettingMultipleParamIsOK", "Parameter " + current_vec.at(iparam) + " is not valid, for setting set '" + setting_set + "' and setting '" + current + ", for parameter number " + std::to_string(iparam+1) + ". Please check this!" );
                return false;
            }
        }
    }

    return true;
}
