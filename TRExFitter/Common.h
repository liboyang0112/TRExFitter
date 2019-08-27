#ifndef COMMON_H
#define COMMON_H

/// c++ stuff
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <set>

/// TRExFitter stuff
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/NormFactor.h"

// ROOT stuff
#include "TF1.h"

/// Forward class declaration
class TFile;
class TH1;
class TH1D;

namespace TRExFitter{
    extern int DEBUGLEVEL;
    void SetDebugLevel(int level=0);
    extern bool SHOWYIELDS; // flag to show or not yields in plots
    extern bool SHOWSTACKSIG;  // flag to show signal or not
    extern bool ADDSTACKSIG;  // flag to add signal to total or not
    extern bool SHOWNORMSIG;  // flag to show normalized signal or not
    extern bool SHOWOVERLAYSIG;  // flag to show overlayed signal or not
    extern bool SHOWCHI2;
    extern bool SHOWSTACKSIG_SUMMARY;  // flag to show signal or not in Summary Plot
    extern bool SHOWNORMSIG_SUMMARY;  // flag to show normalized signal or not in Summary Plot
    extern bool SHOWOVERLAYSIG_SUMMARY;  // flag to show overlayed signal or not in Summary Plot
    extern bool LEGENDLEFT;  // flag to show sample names on left aligned in the legend
    extern bool LEGENDRIGHT;  // flag to show sample names on right aligned in the legend
    extern bool PREFITONPOSTFIT;  // flag to show prefit background as dashed line on postfit plots
    extern bool POISSONIZE;
    extern bool SYSTCONTROLPLOTS;
    extern bool SYSTDATAPLOT;
    extern bool SYSTERRORBARS;
    extern bool SPLITHISTOFILES;
    extern bool HISTOCHECKCRASH;
    extern bool REMOVEXERRORS;
    extern bool OPRATIO;
    extern bool NORATIO; // flag to hide ratio pad
    extern double CORRELATIONTHRESHOLD;
    extern bool MERGEUNDEROVERFLOW;
    extern std::map< std::string,std::string > SYSTMAP;
    extern std::map< std::string,std::string > SYSTTEX;
    extern std::map< std::string,std::string > NPMAP;
    extern std::vector< std::string > IMAGEFORMAT;
    extern int NCPU;
    //
    extern std::map< std::string, double > OPTION;
    extern std::map<std::string,TFile*> TFILEMAP;
    extern bool GUESSMCSTATERROR;
    extern bool CORRECTNORMFORNEGATIVEINTEGRAL;
}

const int MAXregions = 100;
const int MAXsamples = 100;
const int MAXsyst = 500;
const int MAXnorm = 10;

TFile* GetFile(const std::string& fileName);
TH1D* HistFromNtuple(const std::string& ntuple, const std::string& variable, int nbin, double xmin, double xmax, const std::string& selection, const std::string& weight, int Nev=-1);
TH1D* HistFromNtupleBinArr(const std::string& ntuple, const std::string& variable, int nbin, double *bins, const std::string& selection, const std::string& weight, int Nev=-1);
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
std::vector< std::pair < std::string,std::vector<double> > > processString(std::string target);

bool StringsMatch(const std::string& s1, const std::string& s2);
int wildcmp(const char *wild, const char *string);

int FindInStringVector(const std::vector<std::string>& v, const std::string& s);
int FindInStringVectorOfVectors(const std::vector<std::vector<std::string> >& v, const std::string& s, const std::string& ss);
double GetSeparation( TH1D* S1, TH1D* B1 );

TH1D* BlindDataHisto( TH1* h_data, TH1* h_bkg, TH1* h_sig, double threshold=0.02, bool takeSqrt=false );
void BlindDataHisto( TH1* h_data, TH1* h_blind );
double convertStoD(std::string toConvert);

bool SmoothHistogram( TH1* h, double nsigma=2. ); // forceFlat: 0 force no flat, 1 force flat, -1 keep it free
void SmoothHistogramTtres( TH1* h);

void DropBins(TH1* h, const std::vector<int> &v);

double CorrectIntegral(TH1* h, double *err=0);

void CloseFiles( const std::set<std::string> &set);

TH1D* MergeHistograms(std::vector<TH1*> hVec);

/**
  * A function to apply ATLAS/PDG rounding rules to values
  * @param A reference to mean value
  * @param A reference to uncertainty
  */
int ApplyATLASrounding(double& mean, double& error);

/**
  * A helper function to round error according to PDG rules
  * @param The value of error that will be rounded
  * @return number of iterations of multiplication/division by 10 needed to reach the same precision for nominal value
  */
int ApplyErrorRounding(double& error, int& sig);

/**
  * A helper function to round value to n decimal palces
  * @param A value that needs to be rounded
  * @param Number of multiplications/divisions by 10 needed to get the value that can be rounded
  */
void RoundToSig(double& value, const int& n);

std::string FloatToPseudoHex(const float value);
std::string DoubleToPseudoHex(const double value);

float HexToFloat(const std::string& s);
double HexToDouble(const std::string& s);

/**
    * A helper function to scale samples (signal) to nominakl SFs
    * @param SampleHist
    * @param Histogram that will be scaled
    */ 
void ScaleNominal(const SampleHist* const sig, TH1* hist);

TH1* CloneNoError(TH1* h,const char* name="");

unsigned int NCharactersInString(const std::string& s,const char c);

bool CheckExpression(const std::string& s);

std::size_t GetSampleIndexFromList(const std::vector<Sample*>& list, const std::string name);

/**
    * Helper function to calculate nominal scale factor, for morphed samples as well
    * @param pointer to SampleHist for which we need to calculate the scale factor
    * @return scale factor
    */
double GetNominalMorphScale(const SampleHist* const sh);

/**
 * Helper function to parsee the string to indetify if the chosen option needs to run the fit
 * @return true if needs to run the fit
 */
bool OptionRunsFit(const std::string& opt);

/**
 * Helper function to make a copy of histogram with no errors in bins
 * This is useful when doing some scaling operations like Add/Divide 
 * without modifying the original uncertainty in bins
 * @param histogram to be copied
 * @return histogramw with no errors
 */
std::unique_ptr<TH1> GetHistCopyNoError(const TH1* const hist);

/// BW added functions to help pad bin numbers in gamma NP plot

std::vector<std::string> mysplit(const std::string & s, char delimiter);
std::string addpad( const std::string & input, const char filler, const unsigned width );
std::string pad_trail( const std::string & input );


#endif
