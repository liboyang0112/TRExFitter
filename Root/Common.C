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

TH1* HistFromFile(string fullName){
  string fileName  = fullName.substr(0,fullName.find_last_of(".")+5);
  string histoName = fullName.substr(fullName.find_last_of(".")+6,string::npos);
  return HistFromFile(fileName,histoName);
}

TH1* HistFromFile(string fileName,string histoName){
  if(fileName=="") return 0x0;
  if(histoName=="") return 0x0;
  cout << "  Extracting histogram  " << histoName << "  from file  " << fileName << "  ..." << endl;
  TH1* h;
  TFile *f = new TFile(fileName.c_str());
  h = (TH1*)f->Get(histoName.c_str())->Clone();
  h->SetDirectory(0);
  f->Close();
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
              fullPath += ".root";
              if(names[i_name]!="" || nameSufs[i_nameSuf]!=""){
                fullPath += "/";
                fullPath += names[i_name];
                fullPath += nameSufs[i_nameSuf];
              }
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

void SmoothSystHistos(TH1* h_nominal,TH1* h_syst_up,TH1* h_syst_down){
  TH1* h_const = (TH1*)h_nominal->Clone("h_const");
  for(int i_bin=1;i_bin<=h_const->GetNbinsX();i_bin++){
    h_const->SetBinContent(i_bin,1);
  }
  float yield_nominal = h_nominal->Integral();
  float yield_syst_up = h_syst_up->Integral();
  float yield_syst_down = h_syst_down->Integral();
  //
  // --- --- --- //
  // Smoothing implementation: just a simple example here...
  // 1) Symmetrize if needed (if up/down variation = nominal)
  if(Separation(h_nominal,h_syst_up) == 0 && Separation(h_nominal,h_syst_down))
    return;
  if(Separation(h_nominal,h_syst_up) == 0){
    h_syst_up->Add(h_syst_down,-1);
    h_syst_up->Add(h_nominal,1);
  }
  if(Separation(h_nominal,h_syst_down) == 0){
    h_syst_down->Add(h_syst_up,-1);
    h_syst_down->Add(h_nominal,1);
  }
  // 2) Smooth
  // turn the syst histos into relative differeneces
  h_syst_up->Add(h_nominal,-1);
  h_syst_down->Add(h_nominal,-1);
  h_syst_up->Divide(h_nominal);
  h_syst_down->Divide(h_nominal);
  // add an offset
  h_syst_up->Add(h_const,100);
  h_syst_down->Add(h_const,100);
  // --- smooth ---
  h_syst_up->Smooth(2);
  h_syst_down->Smooth(2);
  // ---        ---
  // go back to syst histo
  h_syst_up->Add(h_const,-100);
  h_syst_down->Add(h_const,-100);
  h_syst_up->Multiply(h_nominal);
  h_syst_down->Multiply(h_nominal);
  h_syst_up->Add(h_nominal,1);
  h_syst_down->Add(h_nominal,1);
  //
  delete h_const;
}

float Separation(TH1* h1,TH1* h2){
  float sep = 0;
  for(int i_bin=1;i_bin<=h1->GetNbinsX();i_bin++){
    sep += TMath::Abs( h1->GetBinContent(i_bin) - h2->GetBinContent(i_bin) );
  }
  return sep;
}
