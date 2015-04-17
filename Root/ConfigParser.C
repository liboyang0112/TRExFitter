#include "TtHFitter/ConfigParser.h" 

// -- Functions --

string Fix(string s){ // removes leading and trailing white spaces, remove everything after "%", remove '"'...
  string ss = s.substr( 0, s.find_first_of('%') );
  replace( ss.begin(), ss.end(), '"', ' ');
  ss = ss.substr( ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1 );
//   cout << "\"" << ss << "\"" << endl;
  return ss;
}

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

string First(string s){
  string first;
  first = s.substr( 0, s.find_first_of(':') );
  return Fix(first);
}

string Second(string s){
  string second;
  second = s.substr( s.find_first_of(':')+1,string::npos );
  return Fix(second);
}


// -- Classes -- 

// Config

Config::Config(){};
Config::~Config(){};

// ConfigSet

ConfigSet::ConfigSet(){
  fN = 0;
}

ConfigSet::~ConfigSet(){}
  
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

string ConfigSet::Get(string name){
  for(int i=0;i<fN;i++){
    if(fConfig[i].fName == name){
      return fConfig[i].fValue;
    }
  }
  return "";
}

string ConfigSet::operator[](string name){
  return Get(name);
}

Config ConfigSet::GetConfig(int i){
  return fConfig[i];
}

string ConfigSet::GetConfigName(int i){
  return fConfig[i].fName;
}

string ConfigSet::GetConfigValue(int i){
  return fConfig[i].fValue;
}

int ConfigSet::GetN(){
  return fN;
}

int ConfigSet::size(){
  return fN;
}

void ConfigSet::Set(string name,string value){
  fName = name;
  fValue = value;
}

string ConfigSet::GetName(){
  return fName;
}

string ConfigSet::GetValue(){
  return fValue;
}


// ConfigParser

ConfigParser::ConfigParser(){
  fN = 0;
}
ConfigParser::~ConfigParser(){
}

void ConfigParser::ReadFile(string fileName){
  cout << "Reading config file " << fileName << endl;
  ifstream file(fileName.c_str());
  string str;
  bool reading = false;
  int n = 1;
  int k = 0;
  vector<string> valVec;
  while (getline(file, str)){
    replace( str.begin(), str.end(), '\n', ' ');
    replace( str.begin(), str.end(), '\r', ' ');
    if(str[0]=='%') continue;
    if(str.find_first_not_of(' ')==string::npos){
      reading = false;
    }
    else{
      valVec = Vectorize(Second(str),';');
      if(!reading){
        n = valVec.size();
        for(k=0;k<n;k++){
          fConfSets[fN] = new ConfigSet();
          fConfSets[fN]->Set( First(str),valVec[k] );
          fN++;
        }
        reading = true;
      }
      else{
        for(k=0;k<n;k++){
          fConfSets[fN-n+k]->SetConfig(First(str),valVec[k]);
        }
      }
    }
  }
}

ConfigSet *ConfigParser::GetConfigSet(int i){ // returns the i-th configSet
  return fConfSets[i];
}

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
