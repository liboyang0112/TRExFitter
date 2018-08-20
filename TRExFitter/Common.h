#ifndef COMMON_H
#define COMMON_H

/// c++ stuff
#include <string>
#include <vector>
#include <map>
#include <set>

/// Forward class declaration
class TFile;
class TH1;
class TH1D;

namespace TRExFitter{
    extern int DEBUGLEVEL;
    void SetDebugLevel(int level=0);
    extern bool SHOWYIELDS; // flag to show or not yields in plots
    extern bool SHOWSTACKSIG;  // flag to show signal or not
    extern bool SHOWNORMSIG;  // flag to show normalized signal or not
    extern bool SHOWOVERLAYSIG;  // flag to show overlayed signal or not
    extern bool SHOWCHI2;
    extern bool SHOWSTACKSIG_SUMMARY;  // flag to show signal or not in Summary Plot
    extern bool SHOWNORMSIG_SUMMARY;  // flag to show normalized signal or not in Summary Plot
    extern bool SHOWOVERLAYSIG_SUMMARY;  // flag to show overlayed signal or not in Summary Plot
    extern bool LEGENDLEFT;  // flag to show sample names on left aligned in the legend
    extern bool PREFITONPOSTFIT;  // flag to show prefit background as dashed line on postfit plots
    extern bool POISSONIZE;
    extern bool SYSTCONTROLPLOTS;
    extern bool SYSTDATAPLOT;
    extern bool SYSTERRORBARS;
    extern bool SPLITHISTOFILES;
    extern bool HISTOCHECKCRASH;
    extern bool REMOVEXERRORS;
    extern bool NOENDERR;
    extern float CORRELATIONTHRESHOLD;
    extern bool MERGEUNDEROVERFLOW;
    extern std::map< std::string,std::string > SYSTMAP;
    extern std::map< std::string,std::string > SYSTTEX;
    extern std::map< std::string,std::string > NPMAP;
    extern std::vector< std::string > IMAGEFORMAT;
    extern int NCPU;
    //
    extern std::map< std::string, float > OPTION;
    extern std::map<std::string,TFile*> TFILEMAP;
    extern bool GUESSMCSTATERROR;
}

const int MAXregions = 100;
const int MAXsamples = 100;
const int MAXsyst = 500;
const int MAXnorm = 10;

TFile* GetFile(const std::string& fileName);
TH1D* HistFromNtuple(const std::string& ntuple, const std::string& variable, int nbin, float xmin, float xmax, const std::string& selection, const std::string& weight);
TH1D* HistFromNtupleBinArr(const std::string& ntuple, const std::string& variable, int nbin, double *bins, const std::string& selection, const std::string& weight);
TH1* HistFromFile(const std::string& fullName);
TH1* HistFromFile(const std::string& fileName, const std::string& histoName);
void WriteHistToFile(TH1* h, const std::string& fileName, std::string option="UPDATE");
void WriteHistToFile(TH1* h, TFile *f);
void MergeUnderOverFlow(TH1* h);
std::vector<std::string> CreatePathsList(std::vector<std::string> paths, std::vector<std::string> pathSufs,
                                         std::vector<std::string> files, std::vector<std::string> fileSufs,
                                         std::vector<std::string> names, std::vector<std::string> nameSufs);
std::vector<std::string> CombinePathSufs(std::vector<std::string> pathSufs, std::vector<std::string> newPathSufs );
std::vector<std::string> ToVec(const std::string& s);
std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace);

bool StringsMatch(const std::string& s1, const std::string& s2);
int wildcmp(const char *wild, const char *string);

int FindInStringVector(const std::vector<std::string>& v, const std::string& s);
int FindInStringVectorOfVectors(const std::vector<std::vector<std::string> >& v, const std::string& s, const std::string& ss);
double GetSeparation( TH1D* S1, TH1D* B1 );

TH1D* BlindDataHisto( TH1* h_data, TH1* h_bkg, TH1* h_sig, float threshold=0.02 );
void BlindDataHisto( TH1* h_data, TH1* h_blind );
double convertStoD(std::string toConvert);

bool SmoothHistogram( TH1* h, float nsigma=2. ); // forceFlat: 0 force no flat, 1 force flat, -1 keep it free
void SmoothHistogramTtres( TH1* h);

void DropBins(TH1* h, const std::vector<int> &v);

double CorrectIntegral(TH1* h, double *err=0);

void CloseFiles( const std::set<std::string> &set);

TH1D* MergeHistograms(std::vector<TH1*> hVec);

#endif
