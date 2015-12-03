#include "TtHFitter/Common.h"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// VARIABLES
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int TtHFitter::DEBUGLEVEL = 1;
bool TtHFitter::SHOWYIELDS = false;
bool TtHFitter::SHOWNORMSIG = false;
bool TtHFitter::SYSTCONTROLPLOTS = false;
bool TtHFitter::SPLITHISTOFILES = false;
bool TtHFitter::HISTOCHECKCRASH = true;
float TtHFitter::CORRELATIONTHRESHOLD = -1;
std::map <string,string> TtHFitter::SYSTMAP; //(TtHFitter::SYSTMAP).clear();

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight){
    TH1F* h = new TH1F("h","h",nbin,xmin,xmax);
    if(TtHFitter::DEBUGLEVEL>0) cout << "  Extracting histogram from  " << ntuple << "  ..." << endl;
    TChain *t = new TChain();
    t->Add(ntuple.c_str());
    h->Sumw2();
    t->Draw( Form("%s>>h",variable.c_str()), Form("(%s)*(%s)",weight.c_str(),selection.c_str()), "goff");
    MergeUnderOverFlow(h);
    t->~TChain();
    return h;
}

//__________________________________________________________________________________
//
TH1* HistFromFile(string fullName){
    string fileName  = fullName.substr(0,fullName.find_last_of(".")+5);
    string histoName = fullName.substr(fullName.find_last_of(".")+6,string::npos);
    return HistFromFile(fileName,histoName);
}

//__________________________________________________________________________________
//
TH1* HistFromFile(string fileName,string histoName){
    if(fileName=="") return 0x0;
    if(histoName=="") return 0x0;
    if(TtHFitter::DEBUGLEVEL>0) cout << "  Extracting histogram  " << histoName << "  from file  " << fileName << "  ..." << endl;
    TFile *f = new TFile(fileName.c_str());
    TH1 *h = (TH1*)f->Get(histoName.c_str())->Clone();
    h->SetDirectory(0);
    f->Close();
    delete f;
    return h;
}

//__________________________________________________________________________________
//
void WriteHistToFile(TH1* h,string fileName,string option){
    TDirectory *dir = gDirectory;
    TFile *f = new TFile(fileName.c_str(),option.c_str());
    h->Write("",TObject::kOverwrite);
    delete f;
    dir->cd();
}

//__________________________________________________________________________________
//
void MergeUnderOverFlow(TH1* h){
    int nbins = h->GetNbinsX();
    h->AddBinContent( 1, h->GetBinContent(0) ); // merge first bin with underflow bin
    h->SetBinError(   1, sqrt( pow(h->GetBinError(1),2)+pow(h->GetBinError(0),2)) ); // increase the stat uncertainty as well
    h->AddBinContent( nbins, h->GetBinContent(nbins+1) ); // merge first bin with overflow bin
    h->SetBinError(   nbins, sqrt( pow(h->GetBinError(nbins),2)+pow(h->GetBinError(nbins+1),2)) ); // increase the stat uncertainty as well
    // set under/overflow bins to 0
    h->SetBinContent( 0, 0. );
    h->SetBinContent( nbins+1, 0. );
}

//__________________________________________________________________________________
//
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

//__________________________________________________________________________________
//
vector<string> CombinePathSufs( vector<string> pathSufs, vector<string> newPathSufs ){
    vector<string> output; output.clear();
    if(pathSufs.size()==0) pathSufs.push_back("");
    if(newPathSufs.size()==0) newPathSufs.push_back("");
    for(int i=0;i<(int)pathSufs.size();i++){
        for(int j=0;j<(int)newPathSufs.size();j++){
            output.push_back(pathSufs[i]+newPathSufs[j]);
        }
    }
    return output;
}

//__________________________________________________________________________________
//
vector<string> ToVec(string s){
    vector<string> output;
    output.clear();
    output.push_back(s);
    return output;
}

//__________________________________________________________________________________
//
void TtHFitter::SetDebugLevel(int level){
    DEBUGLEVEL = level;
}

//__________________________________________________________________________________
//
string ReplaceString(string subject, const string& search,
                          const string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return subject;
}

//__________________________________________________________________________________
//
int FindInStringVector(std::vector< string > v, string s){
    int idx = -1;
    string s1;
    string s2;
    for(unsigned int i=0;i<v.size();i++){
        s1 = v[i];
        if(s1==s){
            idx = (int)i;
            break;
        }
        // if last character is "*"...
        if(s1[s1.size()-1]=='*'){
            s2 = s1.substr(0,s1.size()-1);
            if(s.find(s2)!=string::npos){
                idx = (int)i;
                break;
            }
        }
    }
    return idx;
}

//__________________________________________________________________________________
//
double GetSeparation( TH1F* S1, TH1F* B1 ) {
  // taken from TMVA!!!
  TH1F* S=new TH1F(*S1);
  TH1F* B=new TH1F(*B1);
  Double_t separation = 0;
  if ((S->GetNbinsX() != B->GetNbinsX()) || (S->GetNbinsX() <= 0)) {
    cout << "<GetSeparation> signal and background"
         << " histograms have different number of bins: "
         << S->GetNbinsX() << " : " << B->GetNbinsX() << endl;
  }
  if (S->GetXaxis()->GetXmin() != B->GetXaxis()->GetXmin() ||
      S->GetXaxis()->GetXmax() != B->GetXaxis()->GetXmax() ||
      S->GetXaxis()->GetXmax() <= S->GetXaxis()->GetXmin()) {
    cout << S->GetXaxis()->GetXmin() << " " << B->GetXaxis()->GetXmin()
         << " " << S->GetXaxis()->GetXmax() << " " << B->GetXaxis()->GetXmax()
         << " " << S->GetXaxis()->GetXmax() << " " << S->GetXaxis()->GetXmin() << endl;
    cout << "<GetSeparation> signal and background"
         << " histograms have different or invalid dimensions:" << endl;
  }
  Int_t    nstep  = S->GetNbinsX();
  Double_t intBin = (S->GetXaxis()->GetXmax() - S->GetXaxis()->GetXmin())/nstep;
  Double_t nS     = S->GetSumOfWeights()*intBin;
  Double_t nB     = B->GetSumOfWeights()*intBin;
  if (nS > 0 && nB > 0) {
    for (Int_t bin=0; bin <= nstep + 1; bin++) {
      Double_t s = S->GetBinContent( bin )/Double_t(nS);
      Double_t b = B->GetBinContent( bin )/Double_t(nB);
 if (s + b > 0) separation += 0.5*(s - b)*(s - b)/(s + b);
    }
    separation *= intBin;
  }
  else {
    cout << "<GetSeparation> histograms with zero entries: "
         << nS << " : " << nB << " cannot compute separation"
         << endl;
    separation = 0;
  }
  return separation;
}

//__________________________________________________________________________________
// Code to blind bins with S/B > threshold
// - the code kills this kind of bins in data
// - in addition a histogram is returned, with bin content 0 or 1 depending on the bin beeing blinded or not
TH1F* BlindDataHisto( TH1* h_data, TH1* h_bkg, TH1* h_sig, float threshold ) {
  TH1F* h_blind = (TH1F*)h_data->Clone();
  for(int i_bin=1;i_bin<h_data->GetNbinsX()+1;i_bin++){
    if( h_sig->GetBinContent(i_bin) / h_bkg->GetBinContent(i_bin) > threshold ){
      std::cout << "Common::INFO: Blinding bin n." << i_bin << std::endl;
      h_data->SetBinContent(i_bin,0.);
      h_blind->SetBinContent(i_bin,1.);
    }
    else{
      h_blind->SetBinContent(i_bin,0.);
    }
  }
  return h_blind;
}
