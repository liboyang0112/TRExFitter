// Header include
#include "TtHFitter/Common.h"

// Framework includes
#include "TtHFitter/HistoTools.h"
#include "TtHFitter/StatusLogbook.h"

// ATLAS stuff
#include "AtlasUtils/AtlasStyle.h"
#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

// ROOT stuff
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObject.h"
#include "TString.h"

#include <sstream>

using namespace std;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// VARIABLES
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int TtHFitter::DEBUGLEVEL = 1;
bool TtHFitter::SHOWYIELDS = false;
bool TtHFitter::SHOWSTACKSIG = true;
bool TtHFitter::SHOWNORMSIG = false;
bool TtHFitter::SHOWOVERLAYSIG = false;
bool TtHFitter::SHOWCHI2 = false;
bool TtHFitter::SHOWSTACKSIG_SUMMARY = true;
bool TtHFitter::SHOWNORMSIG_SUMMARY = false;
bool TtHFitter::SHOWOVERLAYSIG_SUMMARY = false;
bool TtHFitter::LEGENDLEFT = false;
bool TtHFitter::PREFITONPOSTFIT = false;
bool TtHFitter::POISSONIZE = false;
bool TtHFitter::SYSTCONTROLPLOTS = false;
bool TtHFitter::SYSTERRORBARS = true;
bool TtHFitter::SYSTDATAPLOT = false;
bool TtHFitter::SPLITHISTOFILES = false;
bool TtHFitter::HISTOCHECKCRASH = true;
bool TtHFitter::GUESSMCSTATERROR = true;
bool TtHFitter::REMOVEXERRORS = false;
bool TtHFitter::NOENDERR = false;
float TtHFitter::CORRELATIONTHRESHOLD = -1;
bool TtHFitter::MERGEUNDEROVERFLOW = false;
std::map <string,string> TtHFitter::SYSTMAP;
std::map <string,string> TtHFitter::SYSTTEX;
std::map <string,string> TtHFitter::NPMAP;
std::vector <string> TtHFitter::IMAGEFORMAT;
int TtHFitter::NCPU = 1;
//
std::map <string,float> TtHFitter::OPTION;
std::map<string,TFile*> TtHFitter::TFILEMAP;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight){
    TH1F* h = new TH1F("h","h",nbin,xmin,xmax);
    WriteVerboseStatus("Common::HistFromNtuple", "    Extracting histogram " + variable + " from  " + ntuple + "  ...");
    WriteVerboseStatus("Common::HistFromNtuple", "        with weight  (" + weight + ")*("+selection+")  ...");
    TChain *t = new TChain();
    t->Add(ntuple.c_str());
    h->Sumw2();
    TString drawVariable = Form("%s>>h",variable.c_str()), drawWeight = Form("(%s)*(%s)",weight.c_str(),selection.c_str());
    t->Draw(drawVariable, drawWeight, "goff");
    if(TtHFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(h);
    delete t;
    return h;
}

//__________________________________________________________________________________
//
TH1F* HistFromNtupleBinArr(string ntuple, string variable, int nbin, double *bins, string selection, string weight){
    TH1F* h = new TH1F("h","h",nbin,bins);
    WriteVerboseStatus("Common::HistFromNtupleBinArr", "  Extracting histogram " + variable + " from  " + ntuple + "  ...");
    WriteVerboseStatus("Common::HistFromNtupleBinArr", "      with weight  (" + weight + ")*("+selection+")  ...");
    TChain *t = new TChain();
    t->Add(ntuple.c_str());
    h->Sumw2();
    TString drawVariable = Form("%s>>h",variable.c_str()), drawWeight = Form("(%s)*(%s)",weight.c_str(),selection.c_str());
    t->Draw(drawVariable, drawWeight, "goff");
    if(TtHFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(h);
    delete t;
    return h;
}

//__________________________________________________________________________________
//
TFile* GetFile(string fileName){
    auto it = TtHFitter::TFILEMAP.find(fileName);
    if(it != TtHFitter::TFILEMAP.end()) return it->second;
    else {
       TFile *f = new TFile(fileName.c_str());
       TtHFitter::TFILEMAP.insert(std::pair<string,TFile*>(fileName,f));
       return f;
    }
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
    bool hasCustomAsimov = false;
    if (fileName.find("customAsimov") != std::string::npos) hasCustomAsimov = true;
    WriteVerboseStatus("Common::HistFromFile", "  Extracting histogram    " + histoName + "  from file    " + fileName + "    ...");
    TH1 *h = 0x0;
    TFile *f = GetFile(fileName);
    if(f == 0x0){
            WriteErrorStatus("Common::HistFromFile", "cannot find input file '" + fileName + "'");
            return h;
    }
    h = static_cast<TH1*>(f->Get(histoName.c_str()));
    if(h == 0x0){
            if (!hasCustomAsimov) WriteErrorStatus("Common::HistFromFile", "cannot find histogram '" + histoName + "' from input file '" + fileName + "'");
            else WriteDebugStatus("Common::HistFromFile", "cannot find histogram '" + histoName + "' from input file '" + fileName + "', but its customAsimov histogram so this should not be a problem");
            return h;
    }
    h = static_cast<TH1*>(h->Clone());
    if(h!=0x0) h->SetDirectory(0);
    if(TtHFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(h);
    return h;
}

//__________________________________________________________________________________
//
void WriteHistToFile(TH1* h,string fileName,string option){
    TDirectory *dir = gDirectory;
    TFile *f = new TFile(fileName.c_str(),option.c_str());
    h->Write("",TObject::kOverwrite);
    h->SetDirectory(0);
//     h->AddDirectory(0);
    delete f;
    dir->cd();
}

//__________________________________________________________________________________
//
void WriteHistToFile(TH1* h,TFile *f){
    TDirectory *dir = gDirectory;
    f->cd();
    h->Write("",TObject::kOverwrite);
    h->SetDirectory(0);
//     h->AddDirectory(0);
    dir->cd();
}

//__________________________________________________________________________________
//
void MergeUnderOverFlow(TH1* h){
    int nbins = h->GetNbinsX();
    h->AddBinContent( 1, h->GetBinContent(0) ); // merge first bin with underflow bin
    h->SetBinError(     1, sqrt( pow(h->GetBinError(1),2)+pow(h->GetBinError(0),2)) ); // increase the stat uncertainty as well
    h->AddBinContent( nbins, h->GetBinContent(nbins+1) ); // merge first bin with overflow bin
    h->SetBinError(     nbins, sqrt( pow(h->GetBinError(nbins),2)+pow(h->GetBinError(nbins+1),2)) ); // increase the stat uncertainty as well
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
                            fullPath    = paths[i_path];
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
vector< pair < string,vector<double> > > processString(string target) {
  size_t pos = 0;
  vector<pair <string,vector<double> > > output;
  while((pos = target.find("[",pos)) !=string::npos) {
    pair <string, vector<double> > onePair;
    vector<double> values;
    double oneValue;
    int length = target.find("]",pos) - pos;
    std::stringstream ss(target.substr(pos+1,length-1));
    //std::cout << ss.str() << std::endl;
    while (ss>>oneValue){
      values.push_back(oneValue);
      if (ss.peek() == ','){
	ss.ignore();
      }
    }
    onePair.first = target.substr(0,pos);
    onePair.second = values;
    output.push_back(onePair);
    target.erase(0,pos+length+2);
    pos = 0;
  }
  return output; 
}

//__________________________________________________________________________________
// taking into account wildcards on both
bool StringsMatch(std::string s1,std::string s2){
    if(wildcmp(s1.c_str(),s2.c_str())>0 || wildcmp(s2.c_str(),s1.c_str())>0) return true;
    return false;
}

//__________________________________________________________________________________
// taking into account wildcards on first argument
int wildcmp(const char *wild, const char *string) {
    // Written by Jack Handy - <A href="mailto:jakkhandy@hotmail.com">jakkhandy@hotmail.com</A>
    const char *cp = NULL, *mp = NULL;
    while ((*string) && (*wild != '*')) {
        if ((*wild != *string) && (*wild != '?')) {
            return 0;
        }
        wild++;
        string++;
    }
    while (*string) {
        if (*wild == '*') {
          if (!*++wild) {
              return 1;
          }
          mp = wild;
          cp = string+1;
        } else if ((*wild == *string) || (*wild == '?')) {
            wild++;
            string++;
        } else {
            wild = mp;
            string = cp++;
        }
    }
    while (*wild == '*') {
        wild++;
    }
    return !*wild;
}

//__________________________________________________________________________________
//
int FindInStringVector(std::vector< string > v, string s){
    int idx = -1;
    string s1;
//     string s2;
    for(unsigned int i=0;i<v.size();i++){
        s1 = v[i];
        if(StringsMatch(s1,s)){
            idx = (int)i;
            break;
        }
/*        
        if(s1==s){
            idx = (int)i;
            break;
        }
//       // if there is a "*"...
//       int wildcard_pos = s.find("*");
//       if(wildcard_pos!=string::npos){
//               int foo = s1.find(s.substr(0,wildcard_pos)), bar = s1.find(s.substr(wildcard_pos+1));
//               if(foo!=string::npos && bar!=string::npos)
//                       idx = (int)i;
//                       break;
//       }
        // if both first and last character are "*"...
        if(s1[0]=='*' && s1[s1.size()-1]=='*'){
            s2 = s1.substr(1,s1.size()-1);
                if(s.find(s2)!=string::npos){
                    idx = (int)i;
                    break;
                }
        }
        // if last character is "*"...
        else if(s1[s1.size()-1]=='*'){
            s2 = s1.substr(0,s1.size()-1);
            if(s.find(s2)!=string::npos && s2[0]==s[0]){
                idx = (int)i;
                break;
            }
        }
        // if first character is "*"...
        else if(s1[0]=='*'){
            s2 = s1.substr(1,s1.size());
            if(s.find(s2)!=string::npos && s2[s2.size()-1]==s[s.size()-1]){
                idx = (int)i;
                break;
            }
        }*/
    }
    return idx;
}

//__________________________________________________________________________________
//
int FindInStringVectorOfVectors(std::vector< std::vector<string> > v, string s, string ss){
    int idx = -1;
    string s1;
    string s11;
    string s2;
    string s21;
    for(unsigned int i=0;i<v.size();i++){
        s1 = v[i][0];
        s2 = v[i][1];
//         if(s1==s && s2==ss){
        if(StringsMatch(s1,s) && StringsMatch(s2,ss)){
            idx = (int)i;
            break;
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
        WriteErrorStatus("Common::GetSeparation", "signal and background histograms have different number of bins: " + std::to_string(S->GetNbinsX()) + " : " + std::to_string(B->GetNbinsX()));
    }
    if (S->GetXaxis()->GetXmin() != B->GetXaxis()->GetXmin() ||
            S->GetXaxis()->GetXmax() != B->GetXaxis()->GetXmax() ||
            S->GetXaxis()->GetXmax() <= S->GetXaxis()->GetXmin()) {
        WriteErrorStatus("Common::GetSeparation", "signal and background histograms have different or invalid dimensions:");
        WriteErrorStatus("Common::GetSeparation", "Signal Xmin: " + std::to_string(S->GetXaxis()->GetXmin()) + ", background Xmin " + std::to_string(B->GetXaxis()->GetXmin()));
        WriteErrorStatus("Common::GetSeparation", "Signal Xmax: " + std::to_string(S->GetXaxis()->GetXmax()) + ", background Xmax " + std::to_string(B->GetXaxis()->GetXmax()));
    }
    Int_t        nstep  = S->GetNbinsX();
    Double_t intBin = (S->GetXaxis()->GetXmax() - S->GetXaxis()->GetXmin())/nstep;
    Double_t nS         = S->GetSumOfWeights()*intBin;
    Double_t nB         = B->GetSumOfWeights()*intBin;
    if (nS > 0 && nB > 0) {
        for (Int_t bin=0; bin <= nstep + 1; bin++) {
            Double_t s = S->GetBinContent( bin )/Double_t(nS);
            Double_t b = B->GetBinContent( bin )/Double_t(nB);
    if (s + b > 0) separation += 0.5*(s - b)*(s - b)/(s + b);
        }
        separation *= intBin;
    }
    else {
        WriteErrorStatus("Common::GetSeparation", "histograms with zero entries: signal: " + std::to_string(nS) + " : background : " + std::to_string(nB) + " cannot compute separation");
        separation = 0;
    }
    return separation;
}

//__________________________________________________________________________________
// Code to blind bins with S/B > threshold
// - the code kills this kind of bins in data
// - in addition a histogram is returned, with bin content 0 or 1 depending on the bin beeing blinded or not
TH1F* BlindDataHisto( TH1* h_data, TH1* h_bkg, TH1* h_sig, float threshold ) {
    TH1F* h_blind = (TH1F*)h_data->Clone("h_blind");
    for(int i_bin=1;i_bin<h_data->GetNbinsX()+1;i_bin++){
        if( h_sig->GetBinContent(i_bin) / h_bkg->GetBinContent(i_bin) > threshold ){
            WriteDebugStatus("Common::BlindDataHisto", "Blinding bin n." + std::to_string(i_bin));
            h_data->SetBinContent(i_bin,0.);
            h_blind->SetBinContent(i_bin,1.);
        }
        else{
            h_blind->SetBinContent(i_bin,0.);
        }
    }
    return h_blind;
}

//__________________________________________________________________________________
// This one to blind according to a given histogram containing already info on bins to blind
void BlindDataHisto( TH1* h_data, TH1* h_blind ) {
    for(int i_bin=1;i_bin<h_data->GetNbinsX()+1;i_bin++){
        if(h_blind->GetBinContent(i_bin)!=0){
            WriteDebugStatus("Common::BlindDataHisto", "Blinding bin n." + std::to_string(i_bin));
            h_data->SetBinContent(i_bin,0.);
        }
    }
}

//__________________________________________________________________________________
//
double convertStoD(string toConvert){
    double converted;
    std::string::size_type pos;
    try{
        converted = std::stod(toConvert, &pos);
        if(pos != toConvert.size()){
            WriteErrorStatus("Common::BlindDataHisto", "Convert string -> double, partially converted object: " +  toConvert);
            exit(1);
        }
    }
    catch(const std::exception& err){
        WriteErrorStatus("Common::BlindDataHisto", "Convert string -> double, exception catched: " + toConvert +    " " + err.what());
        exit(1);
    }
    return converted;
}

//__________________________________________________________________________________
//
struct BinNom {
    double N;
    double dN2;
    double edge;
    BinNom(double _N, double _dN2, double _edge) { N = _N; dN2 = _dN2; edge = _edge; }
};

//__________________________________________________________________________________
//
bool systFluctuationNominal(std::vector<BinNom> &hist) {
    auto dM = [](const BinNom &b) { return sqrt(b.dN2); };
    //auto dMoverN = [dM](const BinNom &b) {
    //    double N = b.N;
    //    if (N == 0) N = 1e-16;
    //    return dM(b)/N;
    //};
    auto N = [](const BinNom &b) {
        return b.N;
    };
    int Nbins = hist.size();
    for (int k = 1; k < Nbins; ++k) {
        double variation_prev = fabs(N(hist[k]) - N(hist[k-1]));
        double err = max(dM(hist[k]), dM(hist[k-1]));
        if (variation_prev < err) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
void SmoothHistogramTtres( TH1* h) {
    double origIntegral = h->Integral();

    h->Smooth(2);

    if(h->Integral()!=0){
        h->Scale(origIntegral/h->Integral());
    }
}

//__________________________________________________________________________________
// to smooth a nominal histogram, taking into account the statistical uncertinaty on each bin (note: no empty bins, please!!)
bool SmoothHistogram( TH1* h, float nsigma ){
    int nbinsx = h->GetNbinsX();
    double error;
    float integral = h->IntegralAndError(1,h->GetNbinsX(),error);
    //
    // if not flat, go on with the smoothing
    int Nmax = 5;
    for(int i=0;i<Nmax;i++){
        TH1* h0 = (TH1*)h->Clone("h0");
        h->Smooth();
        bool changesApplied = false;
        for(int i_bin=1;i_bin<=nbinsx;i_bin++){
            if( TMath::Abs(h->GetBinContent(i_bin) - h0->GetBinContent(i_bin)) > nsigma*h0->GetBinError(i_bin) ){
                h->SetBinContent(i_bin,h0->GetBinContent(i_bin));
            }
            else{
                changesApplied = true;
            }
            // bring bins < 1e-6 to 1e-06
            if(h->GetBinContent(i_bin)<1e-06) h->SetBinContent(i_bin,1e-06);
        }
        if(!changesApplied) break;
        delete h0;
    }

    //
    // try to see if it's consistent with being flat
    //
    // make sure you didn't change the integral
    if(h->Integral()>0){
        h->Scale(integral/h->Integral());
    }
    //
    // fix stat error so that the total stat error is unchanged, and it's distributed among all bins
    for(int i_bin=1;i_bin<=nbinsx;i_bin++){
        float N = integral;
        float E = error;
        float n = h->GetBinContent(i_bin);
        h->SetBinError(i_bin,E*sqrt(n)/sqrt(N));
    }
    //
    return false; // this is actual behaviour that was implemented previously
}

//__________________________________________________________________________________
//
void DropBins(TH1* h,const std::vector<int> &v){
    TH1* h_new = (TH1*)h->Clone(h->GetName());
    for(int i_bin=1;i_bin<=h_new->GetNbinsX();i_bin++){
        if(find(v.begin(),v.end(),i_bin-1)!=v.end()){
            h->SetBinContent(i_bin,-1.);
            h->SetBinError(i_bin,0.);
        }
    }
}

//__________________________________________________________________________________
//
float CorrectIntegral(TH1* h,float *err){
    float integral = 0.;
    float error = 0.;
    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
        if(h->GetBinContent(i_bin)<0) continue;
        integral+=h->GetBinContent(i_bin);
        if(h->GetBinError(i_bin)<=0) continue;
        error += pow(h->GetBinError(i_bin),1);
    }
    if(err!=0) *err = sqrt(error);
    return integral;
}

//__________________________________________________________________________________
//
void CloseFiles( const std::set < std::string> &files_names ){
    for( const auto &fullName : files_names ){
        std::string file = fullName.substr(0,fullName.find_last_of(".")+5);
        auto it = TtHFitter::TFILEMAP.find(file);
        if(it != TtHFitter::TFILEMAP.end()){
            //the file exists. Let's close it, and delete the pointer
            it->second->Close();
            TtHFitter::TFILEMAP.erase(file);
        }
    }
}

//__________________________________________________________________________________
//
TH1F* MergeHistograms(vector<TH1*> hVec){
    if(hVec.size()==0) return 0x0;
    if(hVec[0]==0x0) return 0x0;
    // get total number of bins
    int Nbins = 0;
    for(auto h : hVec){
        Nbins += h->GetNbinsX();
    }
    // build array of bin edges
    float *bins = new float[Nbins];
    // coutner
    int k_bin = 0;
    // first edge from first histogram
    bins[0] = hVec[0]->GetXaxis()->GetBinLowEdge(1);
    k_bin ++;
    // define the offset, which will be increased by the last bin UpEdge of a histogram at the end of the loop on its bins
    float offset = 0;
    //
    for(auto h : hVec){
        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
            bins[k_bin] = h->GetXaxis()->GetBinUpEdge(i_bin) + offset;
            if(i_bin==h->GetNbinsX()) offset += h->GetXaxis()->GetBinUpEdge(i_bin)-h->GetXaxis()->GetBinLowEdge(1);
            k_bin ++;
        }
    }
    // create the new histogram
    TH1F *hOut = new TH1F("h_merge","h_merge",Nbins,bins);
    hOut->SetTitle(hVec[0]->GetTitle());
    hOut->SetLineColor(hVec[0]->GetLineColor());
    hOut->SetLineStyle(hVec[0]->GetLineStyle());
    hOut->SetLineWidth(hVec[0]->GetLineWidth());
    hOut->SetFillColor(hVec[0]->GetFillColor());
    hOut->SetFillStyle(hVec[0]->GetFillStyle());
    // fill it
    k_bin = 1;
    for(auto h : hVec){
        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
            hOut->SetBinContent(k_bin,h->GetBinContent(i_bin));
            hOut->SetBinError(k_bin,h->GetBinError(i_bin));
            k_bin ++;
        }
    }
    delete [] bins;
    // return
    return hOut;
}
