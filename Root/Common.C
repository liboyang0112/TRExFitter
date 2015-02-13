#include "TtHFitter/Common.h"

// -------------------------------------------------------------------------------------------------
// FUNCTIONS

TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight){
  TH1F* h = new TH1F("h","h",nbin,xmin,xmax);
  cout << "  Extracting histogram from  " << ntuple << "  ..." << endl;
  TChain *t = new TChain();
  t->Add(ntuple.c_str());
  h->Sumw2();
  t->Draw( Form("%s>>h",variable.c_str()), Form("(%s)*(%s)",weight.c_str(),selection.c_str()), "goff");
  t->~TChain();
  return h;
}

TH1* HistFromFile(string fileName,string histoName){
  if(fileName=="") return 0x0;
  if(histoName=="") return 0x0;
  TH1* h;
  TDirectory *dir = gDirectory;
  TFile *f = new TFile(fileName.c_str());
  h = (TH1*)f->Get(histoName.c_str())->Clone();
  dir->cd();
  return h;
}

void WriteHistToFile(TH1* h,string fileName,string option){
  TDirectory *dir = gDirectory;
  TFile *f = new TFile(fileName.c_str(),option.c_str());
  h->Write("",TObject::kOverwrite);
  f->~TFile();
  dir->cd();
}

vector<string> CreatePathsList( vector<string> paths, vector<string> pathSufs, 
                                vector<string> files, vector<string> fileSufs, 
                                vector<string> names, vector<string> nameSufs){
  // turn the empty vectors into vectors containing one "" entry
  if(paths.size()==0) paths.push_back("");
  if(pathSufs.size()==0) pathSufs.push_back("");
  if(files.size()==0) files.push_back("");
  if(fileSufs.size()==0) fileSufs.push_back("");
  if(names.size()==0) names.push_back("");
  if(nameSufs.size()==0) nameSufs.push_back("");
  //
  vector<string> output;
  string fullPath;
  output.clear();
  for(int i_path=0;i_path<(int)paths.size();i_path++){
    for(int i_pathSuf=0;i_pathSuf<(int)pathSufs.size();i_pathSuf++){
      for(int i_file=0;i_file<(int)files.size();i_file++){
        for(int i_fileSuf=0;i_fileSuf<(int)fileSufs.size();i_fileSuf++){
          for(int i_name=0;i_name<(int)names.size();i_name++){
            for(int i_nameSuf=0;i_nameSuf<(int)nameSufs.size();i_nameSuf++){
              fullPath  = paths[i_path];
              fullPath += pathSufs[i_pathSuf];
              fullPath += "/";
              fullPath += files[i_file];
              fullPath += fileSufs[i_fileSuf];
              fullPath += ".root/";
              fullPath += names[i_name];
              fullPath += nameSufs[i_nameSuf];
//               cout << "   Adding " << fullPath << endl;
              output.push_back( fullPath );
            }
          }
        }
      }
    }
  }
  return output;
}

vector<string> CombinePathSufs( vector<string> pathSufs, vector<string> newPathSufs ){
  vector<string> output; output.clear();
  if(pathSufs.size()==0) pathSufs.push_back("");
  if(newPathSufs.size()==0) newPathSufs.push_back("");
  for(int i=0;i<(int)pathSufs.size();i++){
    for(int j=0;j<(int)newPathSufs.size();j++){
      output.push_back(pathSufs[i]+"/"+newPathSufs[j]);
    }
  }
  return output;
}

vector<string> ToVec(string s){
  vector<string> output;
  output.clear();
  output.push_back(s);
  return output;
}


void TtHFitter::SetDebugLevel(int level){
  DEBUGLEVEL = level;
}


// string RemovePrefix(string s,string prefix){
//   string output = s;
//   s = 
// }
std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return subject;
}