#include "TtHFitter/FitResults.h"

FitResults::FitResults(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
    fCorrMatrix = 0;
    fNuisPar.clear();
}
FitResults::~FitResults(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
    if(fCorrMatrix) delete fCorrMatrix;
    
    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        if(fNuisPar[i]) delete fNuisPar[i];
    }
    fNuisPar.clear();
}

void FitResults::AddNuisPar(NuisParameter *par){
  fNuisPar.push_back(par);
  string p = par->fName;
  fNuisParIdx[p] = (int)fNuisParNames.size();
  fNuisParNames.push_back(p);
  fNuisParIsThere[p] = true;
}

float FitResults::GetNuisParValue(string p){
  int idx = -1;
  if(fNuisParIsThere[p]){
    idx = fNuisParIdx[p];
  }
  else{
    cout << "  WARNING: NP " << p << " not found... Returning 0." << endl;
    return 0.;
  }
  return fNuisPar[idx]->fFitValue;
}

float FitResults::GetNuisParErrUp(string p){
  int idx = -1;
  if(fNuisParIsThere[p]){
    idx = fNuisParIdx[p];
  }
  else{
    cout << "  WARNING: NP " << p << " not found... Returning error = 1." << endl;
    return 1.;
  }
  return fNuisPar[idx]->fPostFitUp;
}

float FitResults::GetNuisParErrDown(string p){
  int idx = -1;
  if(fNuisParIsThere[p]){
    idx = fNuisParIdx[p];
  }
  else{
    cout << "  WARNING: NP " << p << " not found... Returning error = 1." << endl;
    return 1.;
  }
  return fNuisPar[idx]->fPostFitDown;
}

void FitResults::ReadFromTXT(string fileName){
  bool includeCorrelations = true;
  bool invertedCorrMatrix = true;
  bool print = true;
  //
  CorrelationMatrix* matrix = new CorrelationMatrix();
  NuisParameter *np;
  //
  // get fitted NP's
  std::ifstream in;
  in.open(fileName.c_str());
  string input;
  string line;
  bool readingNP = false;
  bool readingCM = false;
  int i = 0;
  int j = 0;
  int Nsyst_corr = 0;
  //
  string name;
  float value, up, down;
  float corr;
  //
  // read file line by line
  while(std::getline(in, line)){
    if(line=="") continue;
    if(line=="NUISANCE_PARAMETERS"){
      cout << "--------------------" << endl;
      cout << "Reading Nuisance Parameters..." << endl;
      cout << "--------------------" << endl;
      readingNP = true;
      continue;
    }
    else if(line=="CORRELATION_MATRIX"){
      cout << "--------------------" << endl;
      cout << "Reading Correlation Matrix..." << endl;
      cout << "--------------------" << endl;
      readingNP = false;
      readingCM = true;
      std::getline(in, line); // skip 1 line
      Nsyst_corr = atof(line.substr(0,line.find(" ")).c_str());
      continue;
    }
    std::istringstream iss(line);
    if(readingNP){
      iss >> input;
      if(input=="" || input=="CORRELATION_MATRIX"){
        readingNP = false;
      }
      while(input.find("\\")!=string::npos) input = input.replace(input.find("\\"),1,"");
      name = input;
      // clean the syst name...
      name = ReplaceString(name,"alpha_","");
      name = ReplaceString(name,"gamma_","");
      AddNuisPar(new NuisParameter(name));
      iss >> value >> up >> down;
      np = fNuisPar[fNuisParIdx[name]];
      np->fFitValue = value;
      np->fPostFitUp = up;
      np->fPostFitDown = down;
      if(print) cout << name << ": " << value << " +" << up << " " << down << endl;
      i++;
    }
      
    if(readingCM){
      if(!includeCorrelations) break;
      for(int i_sys=0;i_sys<Nsyst_corr;i_sys++){
        iss >> corr;
          if(invertedCorrMatrix){
              matrix->SetCorrelation(fNuisParNames[Nsyst_corr-i_sys-1],fNuisParNames[j],corr);
          }
        else matrix->SetCorrelation(fNuisParNames[i_sys],fNuisParNames[j],corr);
      }
      j++;
    }
  }
    if(includeCorrelations){
        if(print){
            for(int i_sys=0;i_sys<Nsyst_corr;i_sys++){
                for(int j_sys=0;j_sys<Nsyst_corr;j_sys++){
                    cout << Form("\t%.4f",matrix->fMatrix[i_sys][j_sys]);
                }
                cout << endl;
            }
        }
    }
  fCorrMatrix = matrix;
  //
  int TOTsyst = fNuisParNames.size();
  cout << "Found " << TOTsyst << " systematics." << endl;
  if(TOTsyst<=0) cout << "WARNING: No systematics found in fit result file..." << endl;
    
}
