// Header include
#include "TRExFitter/Common.h"

// Framework includes
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/StatusLogbook.h"

// ATLAS stuff
#include "AtlasUtils/AtlasStyle.h"
#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

// ROOT stuff
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

// c++ stuff
#include <iostream>
#include <iomanip>
#include <numeric>

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// VARIABLES
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int TRExFitter::DEBUGLEVEL = 1;
bool TRExFitter::SHOWYIELDS = false;
bool TRExFitter::SHOWSTACKSIG = true;
bool TRExFitter::ADDSTACKSIG = true;
bool TRExFitter::SHOWNORMSIG = false;
bool TRExFitter::SHOWOVERLAYSIG = false;
bool TRExFitter::SHOWCHI2 = false;
bool TRExFitter::SHOWSTACKSIG_SUMMARY = true;
bool TRExFitter::SHOWNORMSIG_SUMMARY = false;
bool TRExFitter::SHOWOVERLAYSIG_SUMMARY = false;
bool TRExFitter::LEGENDLEFT = false;
bool TRExFitter::LEGENDRIGHT = false;
bool TRExFitter::PREFITONPOSTFIT = false;
bool TRExFitter::POISSONIZE = false;
bool TRExFitter::SYSTCONTROLPLOTS = false;
bool TRExFitter::SYSTERRORBARS = true;
bool TRExFitter::SYSTDATAPLOT = false;
bool TRExFitter::SPLITHISTOFILES = false;
bool TRExFitter::HISTOCHECKCRASH = true;
bool TRExFitter::GUESSMCSTATERROR = true;
bool TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL = false;
bool TRExFitter::REMOVEXERRORS = false;
double TRExFitter::CORRELATIONTHRESHOLD = -1.;
bool TRExFitter::MERGEUNDEROVERFLOW = false;
bool TRExFitter::OPRATIO = false;
bool TRExFitter::NORATIO = false;
std::map <std::string,std::string> TRExFitter::SYSTMAP;
std::map <std::string,std::string> TRExFitter::SYSTTEX;
std::map <std::string,std::string> TRExFitter::NPMAP;
std::vector <std::string> TRExFitter::IMAGEFORMAT;
int TRExFitter::NCPU = 1;
//
std::map<std::string,double> TRExFitter::OPTION;
std::map<std::string,TFile*> TRExFitter::TFILEMAP;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
TH1D* HistFromNtuple(const std::string& ntuple, const std::string& variable, int nbin, double xmin,
                     double xmax, const std::string& selection, const std::string& weight, int Nev){
    TH1D* h = new TH1D("h","h",nbin,xmin,xmax);
    WriteVerboseStatus("Common::HistFromNtuple", "    Extracting histogram " + variable + " from  " + ntuple + "  ...");
    WriteVerboseStatus("Common::HistFromNtuple", "        with weight  (" + weight + ")*("+selection+")  ...");

    bool hasWildcard = false;
    // check whether file actually exists, AccessPathName() returns FALSE if file can be accessed
    // see https://root.cern.ch/root/html602/TSystem.html#TSystem:AccessPathName
    const std::string fileName = ntuple.substr(0,ntuple.find_last_of("/")); // remove tree name from string to obtain path to file
    if (fileName.find('*') != std::string::npos) hasWildcard = true;
    if (gSystem->AccessPathName(fileName.c_str()) == kTRUE && !hasWildcard ){
        if (TRExFitter::HISTOCHECKCRASH) {
            WriteErrorStatus("Common::HistFromNtuple", "Cannot find input file in: " + fileName);
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("Common::HistFromNtuple", "Cannot find input file in: " + fileName);
        }
    }

    TChain *t = new TChain();
    if (t->Add(ntuple.c_str()) == 0 && hasWildcard){
      WriteWarningStatus("Common::HistFromNtuple", "You used wildcards, but added zero files from " + fileName);
    }
    h->Sumw2();
    TString drawVariable = Form("%s>>h",variable.c_str()), drawWeight = Form("(%s)*(%s)",weight.c_str(),selection.c_str());
    if(Nev>=0) t->Draw(drawVariable, drawWeight, "goff", Nev);
    else       t->Draw(drawVariable, drawWeight, "goff");
    if(TRExFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(h);
    delete t;
    return h;
}

//__________________________________________________________________________________
//
TH1D* HistFromNtupleBinArr(const std::string& ntuple, const std::string& variable, int nbin, double *bins,
                           const std::string& selection, const std::string& weight, int Nev){
    TH1D* h = new TH1D("h","h",nbin,bins);
    WriteVerboseStatus("Common::HistFromNtupleBinArr", "  Extracting histogram " + variable + " from  " + ntuple + "  ...");
    WriteVerboseStatus("Common::HistFromNtupleBinArr", "      with weight  (" + weight + ")*("+selection+")  ...");
    TChain *t = new TChain();

    bool hasWildcard = false;
    // check whether file actually exists, AccessPathName() returns FALSE if file can be accessed
    // see https://root.cern.ch/root/html602/TSystem.html#TSystem:AccessPathName
    const std::string fileName = ntuple.substr(0,ntuple.find_last_of("/")); // remove tree name from string to obtain path to file
    if (fileName.find('*') != std::string::npos) hasWildcard = true;
    if (gSystem->AccessPathName(fileName.c_str()) == kTRUE && !hasWildcard ){
        if (TRExFitter::HISTOCHECKCRASH) {
            WriteErrorStatus("Common::HistFromNtupleBinArr", "Cannot find input file in: " + fileName);
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("Common::HistFromNtupleBinArr", "Cannot find input file in: " + fileName);
        }
    }
    if (t->Add(ntuple.c_str()) == 0 && hasWildcard) {
      WriteWarningStatus("Common::HistFromNtupleBinArr", "You used wildcards, but added zero files from " + fileName);
    }
    h->Sumw2();
    TString drawVariable = Form("%s>>h",variable.c_str()), drawWeight = Form("(%s)*(%s)",weight.c_str(),selection.c_str());
    if(Nev>=0) t->Draw(drawVariable, drawWeight, "goff", Nev);
    else       t->Draw(drawVariable, drawWeight, "goff");
    if(TRExFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(h);
    delete t;
    return h;
}

//__________________________________________________________________________________
//
TFile* GetFile(const std::string& fileName){
    auto it = TRExFitter::TFILEMAP.find(fileName);
    if(it != TRExFitter::TFILEMAP.end()) return it->second;
    else {
       TFile *f = new TFile(fileName.c_str());
       TRExFitter::TFILEMAP.insert(std::pair<std::string,TFile*>(fileName,f));
       return f;
    }
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1> HistFromFile(const std::string& fullName){
    std::string fileName  = fullName.substr(0,fullName.find_last_of(".")+5);
    std::string histoName = fullName.substr(fullName.find_last_of(".")+6,std::string::npos);
    return HistFromFile(fileName,histoName);
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1> HistFromFile(const std::string& fileName, const std::string& histoName){
    if(fileName=="") return nullptr;
    if(histoName=="") return nullptr;
    bool hasCustomAsimov = false;
    if (fileName.find("customAsimov") != std::string::npos) hasCustomAsimov = true;
    WriteVerboseStatus("Common::HistFromFile", "  Extracting histogram    " + histoName + "  from file    " + fileName + "    ...");
    std::unique_ptr<TH1> h = nullptr;
    TFile *f = GetFile(fileName);
    if(!f){
            WriteErrorStatus("Common::HistFromFile", "cannot find input file '" + fileName + "'");
            return nullptr;
    }
    h = std::unique_ptr<TH1>(static_cast<TH1*>(f->Get(histoName.c_str())));
    if(!h){
            if (!hasCustomAsimov) WriteErrorStatus("Common::HistFromFile", "cannot find histogram '" + histoName + "' from input file '" + fileName + "'");
            else WriteDebugStatus("Common::HistFromFile", "cannot find histogram '" + histoName + "' from input file '" + fileName + "', but its customAsimov histogram so this should not be a problem");
            return nullptr;
    }
    if(h!=nullptr) h->SetDirectory(0);
    if(TRExFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(h.get());
    return h;
}

//__________________________________________________________________________________
//
void WriteHistToFile(TH1* h, const std::string& fileName, std::string option){
    TDirectory *dir = gDirectory;
    TFile *f = new TFile(fileName.c_str(),option.c_str());
    h->Write("",TObject::kOverwrite);
    h->SetDirectory(0);
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
    // set under/overflow bins and its errors to 0
    h->SetBinContent( 0, 0. );
    h->SetBinContent( nbins+1, 0. );
    h->SetBinError( 0, 0. );
    h->SetBinError( nbins+1, 0. );
}

//__________________________________________________________________________________
//
std::vector<std::string> CreatePathsList( std::vector<std::string> paths, std::vector<std::string> pathSufs,
                                std::vector<std::string> files, std::vector<std::string> fileSufs,
                                std::vector<std::string> names, std::vector<std::string> nameSufs){
    // turn the empty vectors into vectors containing one "" entry
    if(paths.size()==0) paths.push_back("");
    if(pathSufs.size()==0) pathSufs.push_back("");
    if(files.size()==0) files.push_back("");
    if(fileSufs.size()==0) fileSufs.push_back("");
    if(names.size()==0) names.push_back("");
    if(nameSufs.size()==0) nameSufs.push_back("");
    //
    std::vector<std::string> output;
    std::string fullPath;
    for (const auto& ipath : paths) {
        for(const auto& ipathSuf : pathSufs){
            for(const auto& ifile : files){
                for(const auto& ifileSuf : fileSufs){
                    for(const auto& iname : names){
                        for(const auto& inameSuf : nameSufs){
                            fullPath    = ipath;
                            fullPath += ipathSuf;
                            fullPath += "/";
                            fullPath += ifile;
                            fullPath += ifileSuf;
                            fullPath += ".root";
                            if(iname !="" || inameSuf!=""){
                                fullPath += "/";
                                fullPath += iname;
                                fullPath += inameSuf;
                            }
                            output.emplace_back( fullPath );
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
std::vector<std::string> CombinePathSufs( std::vector<std::string> pathSufs,
                                          std::vector<std::string> newPathSufs ){
    std::vector<std::string> output;
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
std::vector<std::string> ToVec(const std::string& s){
    std::vector<std::string> output;
    output.clear();
    output.push_back(s);
    return output;
}

//__________________________________________________________________________________
//
void TRExFitter::SetDebugLevel(int level){
    DEBUGLEVEL = level;
}

//__________________________________________________________________________________
//
std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

//__________________________________________________________________________________
//
std::vector< std::pair < std::string,std::vector<double> > > processString(std::string target) {
    size_t pos = 0;
    std::vector<std::pair <std::string,std::vector<double> > > output;
    while((pos = target.find("[",pos)) !=std::string::npos) {
        std::pair <std::string, std::vector<double> > onePair;
        std::vector<double> values;
        double oneValue;
        int length = target.find("]",pos) - pos;
        std::stringstream ss(target.substr(pos+1,length-1));
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
bool StringsMatch(const std::string& s1, const std::string& s2){
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
int FindInStringVector(const std::vector<std::string>& v, const std::string& s){
    int idx = -1;
    std::string s1;
    for(unsigned int i=0;i<v.size();i++){
        s1 = v[i];
        if(StringsMatch(s1,s)){
            idx = (int)i;
            break;
        }
    }
    return idx;
}

//__________________________________________________________________________________
//
int FindInStringVectorOfVectors(const std::vector< std::vector<std::string> >& v, const std::string& s, const std::string& ss){
    int idx = -1;
    std::string s1;
    std::string s11;
    std::string s2;
    std::string s21;
    for(unsigned int i=0;i<v.size();i++){
        s1 = v[i][0];
        s2 = v[i][1];
        if(StringsMatch(s1,s) && StringsMatch(s2,ss)){
            idx = (int)i;
            break;
        }
    }
    return idx;
}

//__________________________________________________________________________________
//
double GetSeparation( TH1D* S1, TH1D* B1 ) {
    // taken from TMVA!!!
    std::unique_ptr<TH1> S = std::make_unique<TH1D>(*S1);
    std::unique_ptr<TH1> B = std::make_unique<TH1D>(*B1);
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
    Int_t nstep     = S->GetNbinsX();
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
        WriteErrorStatus("Common::GetSeparation", "histograms with zero entries: signal: " + std::to_string(nS) + " : background : " + std::to_string(nB) + " cannot compute separation");
        separation = 0;
    }
    return separation;
}

//__________________________________________________________________________________
// Code to blind bins with (h_data yield) / (h_bkg yield) > threshold
// - the code kills this kind of bins in data
// - also set uncertainties in blinded bins to zero
// - in addition a histogram is returned, with bin content 0 or 1 depending on the bin beeing blinded or not
// when takeSqrt is true, take the sqrt of the denominator when evaluating the blinding
TH1D* BlindDataHisto( TH1* h_data, TH1* h_bkg, TH1* h_sig, double threshold, bool takeSqrt) {
    TH1D* h_blind = (TH1D*)h_data->Clone("h_blind");
    for(int i_bin=1;i_bin<h_data->GetNbinsX()+1;i_bin++){
        double tmpDenominator = h_bkg->GetBinContent(i_bin);
        if(takeSqrt) tmpDenominator = sqrt(tmpDenominator); // for calculating S/sqrt(B) and S/sqrt(S+B)
        if( h_sig->GetBinContent(i_bin) / tmpDenominator > threshold ){
            WriteDebugStatus("Common::BlindDataHisto", "Blinding bin n." + std::to_string(i_bin));
            h_data->SetBinContent(i_bin,0.);
            h_data->SetBinError(i_bin,0.);
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
            h_data->SetBinError(i_bin,0.);
        }
    }
}

//__________________________________________________________________________________
//
double convertStoD(std::string toConvert){
    double converted;
    std::string::size_type pos;
    try{
        converted = std::stod(toConvert, &pos);
        if(pos != toConvert.size()){
            WriteErrorStatus("Common::BlindDataHisto", "Convert string -> double, partially converted object: " +  toConvert);
            exit(EXIT_FAILURE);
        }
    }
    catch(const std::exception& err){
        WriteErrorStatus("Common::BlindDataHisto", "Convert string -> double, exception caught: " + toConvert +    " " + err.what());
        exit(EXIT_FAILURE);
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
    auto dM = [](const BinNom &b) {
        return sqrt(b.dN2);
    };
    auto N = [](const BinNom &b) {
        return b.N;
    };
    int Nbins = hist.size();
    for (int k = 1; k < Nbins; ++k) {
        double variation_prev = std::fabs(N(hist[k]) - N(hist[k-1]));
        double err = std::max(dM(hist[k]), dM(hist[k-1]));
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
bool SmoothHistogram( TH1* h, double nsigma ){
    int nbinsx = h->GetNbinsX();
    double error = 0.;
    double integral = h->IntegralAndError(1,h->GetNbinsX(),error);
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
        double N = integral;
        double E = error;
        double n = h->GetBinContent(i_bin);
        h->SetBinError(i_bin,E*sqrt(n)/sqrt(N));
    }
    //
    return false; // this is actual behaviour that was implemented previously
}

//__________________________________________________________________________________
//
void DropBins(TH1* h,const std::vector<int> &v){
    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
        if(find(v.begin(),v.end(),i_bin-1)!=v.end()){
            h->SetBinContent(i_bin,-1.);
            h->SetBinError(i_bin,0.);
        }
    }
}

//__________________________________________________________________________________
//
double CorrectIntegral(TH1* h, double * err){
    double integral = 0.;
    double error = 0.;
    for( int i_bin=1; i_bin <= h->GetNbinsX(); i_bin++){
        if(h->GetBinContent(i_bin)<0) continue;
        integral += h->GetBinContent(i_bin);
        if(h->GetBinError(i_bin)<=0) continue;
        error += pow(h->GetBinError(i_bin), 2);
    }
    if(err!=0) *err = sqrt(error);
    return integral;
}

//__________________________________________________________________________________
//
void CloseFiles( const std::set < std::string> &files_names ){
    for( const auto &fullName : files_names ){
        std::string file = fullName.substr(0,fullName.find_last_of(".")+5);
        auto it = TRExFitter::TFILEMAP.find(file);
        if(it != TRExFitter::TFILEMAP.end()){
            //the file exists. Let's close it, and delete the pointer
            it->second->Close();
            TRExFitter::TFILEMAP.erase(file);
        }
    }
}

//__________________________________________________________________________________
//
TH1D* MergeHistograms(const std::vector<TH1*>& hVec){
    if(hVec.size()==0) return nullptr;
    if(hVec[0]==nullptr) return nullptr;
    // build vector of bin edges
    std::vector<double> binVec;
    binVec.push_back( hVec[0]->GetXaxis()->GetBinLowEdge(1) );
    // define the offset, which will be increased by the last bin UpEdge of a histogram at the end of the loop on its bins
    double offset = 0;
    //
    for(unsigned int i_h=0;i_h<hVec.size();i_h++){
        TH1* h = hVec[i_h];
        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
            if(i_h==0) binVec.push_back( h->GetXaxis()->GetBinUpEdge(i_bin) + offset );
            else       binVec.push_back( h->GetXaxis()->GetBinUpEdge(i_bin) - h->GetXaxis()->GetBinLowEdge(1) + offset );
            if(i_bin==h->GetNbinsX()){
                if(i_h==0) offset += h->GetXaxis()->GetBinUpEdge(i_bin);
                else       offset += h->GetXaxis()->GetBinUpEdge(i_bin) - h->GetXaxis()->GetBinLowEdge(1);
            }
        }
    }
    int Nbins = binVec.size()-1;
    // create the new histogram
    TH1D* hOut = new TH1D("h_merge","h_merge",Nbins,&binVec[0]);
    hOut->SetTitle(hVec[0]->GetTitle());
    hOut->SetLineColor(hVec[0]->GetLineColor());
    hOut->SetLineStyle(hVec[0]->GetLineStyle());
    hOut->SetLineWidth(hVec[0]->GetLineWidth());
    hOut->SetFillColor(hVec[0]->GetFillColor());
    hOut->SetFillStyle(hVec[0]->GetFillStyle());
    // fill it
    int k_bin = 1;
    for(const auto& h : hVec){
        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
            hOut->SetBinContent(k_bin,h->GetBinContent(i_bin));
            hOut->SetBinError(k_bin,h->GetBinError(i_bin));
            k_bin ++;
        }
    }
    // return
    return hOut;
}

//___________________________________________________________
//
int ApplyATLASrounding(double &mean, double &error){
    if (error < 0 ){
        WriteWarningStatus("Common::ApplyATLASrounding", "Error value is < 0. Not applying rounding.");
        return -1;
    }

    int sig = 0;
    int iterations = ApplyErrorRounding(error,sig);
    if (iterations > 100) { // something went wrong
        WriteWarningStatus("Common::ApplyATLASrounding", "Problem with applying PDG rounding rules to error.");
        return -1;
    }

    // now apply the correct rounding for nominal value
    RoundToSig(mean, iterations);

    // return the number of decimal digits (for later printing avoiding exponent...)
    int decPlaces = iterations;
    if(iterations<0) decPlaces = 0;
    return decPlaces;
}

//___________________________________________________________
// FIXME : still to fix the 100
int ApplyErrorRounding(double& error,int& sig){
    int iterations = 0;

    if (error == 0) {
        WriteWarningStatus("Common::ApplyErrorRounding", "Error is zero, you should have a look at this.");
        return 0;
    }

    while (error < 100) {
        error*= 10;
        iterations++;
        if (iterations > 15){
            WriteWarningStatus("Common::ApplyErrorRounding", "Too many iterations in determination of decimal places. Not applying rounding");
            return 999;
        }
    }

    while (error >= 1000) {
        error/= 10;
        iterations--;
        if (iterations < -15){
            WriteWarningStatus("Common::ApplyErrorRounding", "Too many iterations in determination of decimal places. Not applying rounding");
            return 999;
        }
    }

    // PDG rounding rules
    if (error >= 100 && error < 355){
        sig = 2;
    } else if (error >= 355 && error < 950) {
        sig = 1;
    } else if (error >= 950 && error < 1000) {
        error = 1000;
        sig = 1;
    } else {
        WriteWarningStatus("Common::ApplyErrorRounding", "3 significant digit are < 100 or > 999. This should not happen.");
        return 999;
    }

    // have three significant digits, now round
    // according to the number of decimal places
    error/= std::pow(10, (3-sig));
    error = std::round(error);
    error*= std::pow(10, (3-sig));

    // now we need to get back to original value
    // this is not optimal but should be optimized by compiler
    if(iterations>0) error/= std::pow(10, std::abs(iterations));
    if(iterations<0) error*= std::pow(10, std::abs(iterations));

    // return number of iterations needed minus 2
    // this will be used to match precision of mean to precision
    // of rounded error
    return (iterations - 3 + sig);
}

//___________________________________________________________
//
void RoundToSig(double& value, const int& n){
    if (n == 0) {
        value = std::round(value);
        return;
    }

    if (n > 0) { // will multiply
        value*= std::pow(10,n);
        value = std::round(value);
        value/= std::pow(10,n);
    } else if (n < 0) { // will divide
        value/= std::pow(10,std::abs(n));
        value = std::round(value);
        value*= std::pow(10,std::abs(n));
    }
}

//___________________________________________________________
//
unsigned int NCharactersInString(const std::string& s,const char c){
    unsigned int N = 0;
    for(unsigned int i_c=0;i_c<s.size();i_c++){
        if(s[i_c]==c) N++;
    }
    return N;
}

//___________________________________________________________
// for the moment just checks the number of parenthesis, but can be expanded
bool CheckExpression(const std::string& s){
    if(s.find("Alt$")!=std::string::npos){
        return true;
    }
    int nParOpen = NCharactersInString(s,'(');
    int nParClose = NCharactersInString(s,')');
    if(nParOpen!=nParClose) return false;
    // ...
    return true;
}

//----------------------------------------------------------------------------------
//
std::string FloatToPseudoHex(const float value){
    std::string s = std::to_string(value);
    std::string first = s.substr(0,s.find('.'));
    std::string second = s.substr(s.find('.')+1, s.length());

    //Count the number of "0" after the comma
    int count = 0;
    for (unsigned int i = 0; i < second.size(); i++) {
      if (second[i] != '0')
        break;
      count++;
    }

    int value1 = std::stoi(first);
    const int value2 = std::stoi(second);

    // add 1234 to the first digit so it is not easily readable, we will subtract it in the decoding
    value1+=1234;
    // add 5678 to the number of '0'
    count+=5678;

    std::stringstream ss;
    ss << std::hex << value1 << "." << std::hex << count  << "." << std::hex << value2;

    return ss.str();
}


//----------------------------------------------------------------------------------
//
std::string DoubleToPseudoHex(const double value){
    std::string s = std::to_string(value);
    std::string first = s.substr(0,s.find('.'));
    std::string second = s.substr(s.find('.')+1, s.length());

    //Count the number of "0" after the comma
    int count = 0;
    for (unsigned int i = 0; i < second.size(); i++) {
      if (second[i] != '0')
        break;
      count++;
    }

    int value1 = std::stoi(first);
    const int value2 = std::stoi(second);

    // add 1234 to the first digit so it is not easily readable, we will subtract it in the decoding
    value1+=1234;
    // add 5678 to the number of '0'
    count+=5678;

    std::stringstream ss;
    ss << std::hex << value1 << "." << std::hex << count  << "." << std::hex << value2;

    return ss.str();
}

//----------------------------------------------------------------------------------
//
float HexToFloat(const std::string& s){
    std::string first = s.substr(0,s.find('.'));
    std::string rest = s.substr(s.find('.')+1, s.length());
    std::string zeros = rest.substr(0,rest.find('.'));
    std::string second = rest.substr(rest.find('.')+1, rest.length());

    unsigned int i1, i2, n0;

    std::stringstream ss;
    ss << std::hex << first;
    ss >> i1;

    std::stringstream ss1;
    ss1 << std::hex << second;
    ss1 >> i2;

    std::stringstream ss2;
    ss2 << std::hex << zeros;
    ss2 >> n0;

    int signed1 = static_cast<int>(i1);
    // need to subtract the 1234 we added
    signed1-= 1234;
    // need to substract the 5678
    n0-= 5678;

    std::string result = std::to_string(signed1)+".";

    for (unsigned int i = 0; i < n0; i++)
      result += "0";

    result += std::to_string(i2);

    return std::stof(result);
}

//----------------------------------------------------------------------------------
//
double HexToDouble(const std::string& s){
    std::string first = s.substr(0,s.find('.'));
    std::string rest = s.substr(s.find('.')+1, s.length());
    std::string zeros = rest.substr(0,rest.find('.'));
    std::string second = rest.substr(rest.find('.')+1, rest.length());

    unsigned int i1, i2, n0;

    std::stringstream ss;
    ss << std::hex << first;
    ss >> i1;

    std::stringstream ss1;
    ss1 << std::hex << second;
    ss1 >> i2;

    std::stringstream ss2;
    ss2 << std::hex << zeros;
    ss2 >> n0;

    int signed1 = static_cast<int>(i1);
    // need to subtract the 1234 we added
    signed1-= 1234;
    // need to substract the 5678
    n0-= 5678;

    std::string result = std::to_string(signed1)+".";

    for (unsigned int i = 0; i < n0; i++)
      result += "0";

    result += std::to_string(i2);

    return std::stod(result);
}


//___________________________________________________________
//
void ScaleNominal(const SampleHist* const sig, TH1* hist){
    for(size_t i_nf=0; i_nf<sig->fSample->fNormFactors.size(); ++i_nf){
        NormFactor *nf = sig->fSample->fNormFactors[i_nf];
        // if this norm factor is a morphing one
        if(nf->fName.find("morph_")!=std::string::npos || nf->fExpression.first!=""){
            std::string formula = TRExFitter::SYSTMAP[nf->fName];
            std::string name = TRExFitter::NPMAP[nf->fName];
            formula = ReplaceString(formula,name,"x");
            auto f_morph = std::unique_ptr<TF1>(new TF1("f_morph",formula.c_str(),nf->fMin,nf->fMax));
            const double& scale = f_morph->Eval(nf->fNominal);
            hist->Scale(scale);
            WriteDebugStatus("Common::ScaleNominal", nf->fName + " => Scaling " + sig->fSample->fName + " by " + std::to_string(scale));
        }
        else{
            hist->Scale(nf->fNominal);
            WriteDebugStatus("Common::ScaleNominal", nf->fName + " => Scaling " + sig->fSample->fName + " by " + std::to_string(sig->fSample->fNormFactors[i_nf]->fNominal));
        }
    }
}

//___________________________________________________________
//
std::size_t GetSampleIndexFromList(const std::vector<Sample*>& list, const std::string name){
    for (std::size_t i = 0; i < list.size(); ++i){
        if (list.at(i)->fName == name) return i;
    }

    return 9999;
}


//____________________________________________________________________________________
//
double GetNominalMorphScale(const SampleHist* const sh){
    double scale = 1.;
    if (!sh) return 1.;
    if (!(sh->fSample)) return 1.;
    for (unsigned int i_nf = 0; i_nf < sh->fSample->fNormFactors.size(); i_nf++){
        NormFactor *nf = sh->fSample->fNormFactors[i_nf];
        if (!nf) continue;
        std::string nfName = nf->fName;

        if(nfName.find("morph_")!=std::string::npos || nf->fExpression.first!=""){
            std::string formula = TRExFitter::SYSTMAP[nfName];
            std::string name = TRExFitter::NPMAP[nfName];
            WriteDebugStatus("Common::GetNominalMorphScale", "formula: " +formula);
            WriteDebugStatus("Common::GetNominalMorphScale", "name: " +name);
            std::vector < std::pair < std::string,std::vector<double> > > nameS;
            if(nfName.find("morph_")!=std::string::npos) {
                nameS.push_back(std::make_pair(name,std::vector<double>{double(nf->fNominal),double(nf->fMin),double(nf->fMax)}));
            } else {
              nameS = processString(name);
            }
            std::vector <double> nfNominalvec;
            for (unsigned int j = 0; j<nameS.size(); j++){
                formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                nfNominalvec.push_back(nameS[j].second[0]);
            }
            double *nfNominal = nfNominalvec.data();
            WriteDebugStatus("Common::GetNominalMorphScale", "formula: " +formula);
            for(unsigned int j = 0; j<nameS.size(); ++j) {
                WriteDebugStatus("Common::GetNominalMorphScale", "nfNominal["+std::to_string(j)+"]: "+std::to_string(nfNominal[j]));
            }
            TFormula f_morph("f_morph",formula.c_str());
            scale *= f_morph.EvalPar(nfNominal,nullptr);
        } else {
            scale *= sh->fSample->fNormFactors[i_nf]->fNominal;
        }
    }

    return scale;
}

//___________________________________________________________
//
bool OptionRunsFit(const std::string& opt){
    if (opt.find("w")!=std::string::npos) return true;
    if (opt.find("f")!=std::string::npos) return true;
    if (opt.find("l")!=std::string::npos) return true;
    if (opt.find("s")!=std::string::npos) return true;
    if (opt.find("r")!=std::string::npos) return true;
    if (opt.find("i")!=std::string::npos) return true;
    if (opt.find("x")!=std::string::npos) return true;
    return false;
}

//___________________________________________________________
//
std::unique_ptr<TH1> GetHistCopyNoError(const TH1* const hist){
    if (hist == nullptr) return nullptr;
    std::unique_ptr<TH1> result(static_cast<TH1*>(hist->Clone()));

    for (int ibin = 0; ibin <= hist->GetNbinsX(); ++ibin){
        result->SetBinError(ibin, 0.);
    }

    return result;
}

// BW helper functions to pad bin numbers for gamma plots
// replaces them with zero padded versions.  "Gamma Bin 1" -> "Gamma Bin 0001"

std::vector<std::string> mysplit(const std::string & s, const char delimiter) {
    std::vector<std::string> answer;
    std::string token;

    // this converts a single char into a string
    // the constructor std:string( n, char )
    // produces a string of n copies of char
    std::string localDelim( 1, delimiter );
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        //keep the delimiter for easy reconstruction
        answer.push_back(token+localDelim);
    }

    //remove the trailing delim from the last token
    std::string last = answer.back();
    if( !last.empty() ) last.pop_back();

    answer.pop_back();
    answer.push_back( last );

    return answer;
}

//___________________________________________________________
//
std::string addpad( const std::string & input, const char filler, const unsigned width ) {
    std::stringstream mySS;

    mySS.fill(filler);
    mySS.width(width);

    mySS << input;

    return mySS.str();

}

//___________________________________________________________
//
std::string pad_trail( const std::string & input ) {

    std::vector<std::string> words = mysplit( input, ' ' );

    std::string paddedValue = addpad( words.back(), '0', 4 );

    words.pop_back();
    words.push_back( paddedValue );

    std::string answer = std::accumulate( words.begin(), words.end(), std::string("") );

    return answer;
}

//___________________________________________________________
// Helper functions to drop norm or shape part from systematic variations 
void DropNorm(TH1* hUp,TH1* hDown,TH1* hNom){
    const double intNom = hNom->Integral();
    if(hUp!=nullptr){
        const double intUp = hUp->Integral();
        if(intUp!=0) hUp->Scale(intNom/intUp);
        else WriteErrorStatus("Common::DropNorm","Integral of up variation = 0. Cannot drop normalization.");
    }
    if(hDown!=nullptr){
        const double intDown = hDown->Integral();
        if(intDown!=0) hDown->Scale(intNom/intDown);
        else WriteErrorStatus("Common::DropNorm","Integral of down variation = 0. Cannot drop normalization.");
    }
}
void DropShape(TH1* hUp,TH1* hDown,TH1* hNom){
    const double intNom = hNom->Integral();
    if(intNom==0){
        WriteErrorStatus("Common::DropShape","Integral of nominal histogram = 0. Cannot drop shape of syst variations.");
        return;
    }
    if(hUp!=nullptr){
        const double ratioUp = hUp->Integral()/intNom;
        delete hUp;
        hUp = static_cast<TH1*>(hNom->Clone());
        hUp->Scale(ratioUp);
    }
    if(hDown!=nullptr){
        const double ratioDown = hDown->Integral()/intNom;
        delete hDown;
        hDown = static_cast<TH1*>(hNom->Clone());
        hDown->Scale(ratioDown);
    }
}

//___________________________________________________________
//
void ScaleMCstatInHist(TH1* hist, const double scale) {
    if (std::fabs(scale-1) < 1e-6) return; // basically scale == 1 but floating precision

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
        hist->SetBinError(ibin, scale * hist->GetBinError(ibin));
    }
}
