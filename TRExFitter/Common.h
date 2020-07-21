#ifndef COMMON_H
#define COMMON_H

/// TRExFitter stuff
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/NormFactor.h"

// ROOT stuff
#include "TF1.h"

/// c++ stuff
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class Region;
class TFile;
class TGraphAsymmErrors;
class TH1;
class TH2;
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
    //
    extern std::map< std::string, double > OPTION;
    extern std::map<std::string, std::shared_ptr<TFile> > TFILEMAP;
    extern bool GUESSMCSTATERROR;
    extern bool CORRECTNORMFORNEGATIVEINTEGRAL;
}

namespace Common {

enum BlindingType {
  SOVERB = 1,
  SOVERSPLUSB = 2,
  SOVERSQRTB = 3,
  SOVERSQRTSPLUSB = 4
};

std::shared_ptr<TFile> GetFile(const std::string& fileName);
TH1D* HistFromNtuple(const std::string& ntuple, const std::string& variable, int nbin, double xmin, double xmax, const std::string& selection, const std::string& weight, int Nev=-1);
TH1D* HistFromNtupleBinArr(const std::string& ntuple, const std::string& variable, int nbin, double *bins, const std::string& selection, const std::string& weight, int Nev=-1);
std::unique_ptr<TH1> HistFromFile(const std::string& fullName);
std::unique_ptr<TH1> HistFromFile(const std::string& fileName, const std::string& histoName);
std::unique_ptr<TH2> Hist2DFromFile(const std::string& fullName);
std::unique_ptr<TH2> Hist2DFromFile(const std::string& fileName, const std::string& histoName);
void WriteHistToFile(TH1* h, const std::string& fileName, const std::string& option="UPDATE");
void WriteHistToFile(TH1* h, std::shared_ptr<TFile> f);
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

/**
  * Function to blind data and retrieve the blinding histogram
  * @param data histogram
  * @param indices ob blinded bins
  * @return Histogram with non-zero bins on postions to be blinded
  */ 
std::unique_ptr<TH1D> BlindDataHisto(TH1* h_data, const std::vector<int>& blindedBins);

/**
  * Function to blind data histogram based on a blinding histogram
  * @param data histogram
  * @param blinding histogram
  */ 
void BlindDataHisto(TH1* h_data, TH1* h_blind);

double convertStoD(std::string toConvert);

bool SmoothHistogram( TH1* h, double nsigma=2. ); // forceFlat: 0 force no flat, 1 force flat, -1 keep it free
void SmoothHistogramTtres( TH1* h);

void DropBins(TH1* h, const std::vector<int> &v);

double CorrectIntegral(TH1* h, double *err=0);

void CloseFiles( const std::set<std::string> &set);

TH1D* MergeHistograms(const std::vector<TH1*>& hVec);

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

std::string KeepSignificantDigits(double value, const int n);

std::string DoubleToPseudoHex(const double value);

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

void ScaleMCstatInHist(TH1* hist, const double scale);

/// BW added functions to help pad bin numbers in gamma NP plot

std::vector<std::string> mysplit(const std::string & s, char delimiter);
std::string addpad( const std::string & input, const char filler, const unsigned width );
std::string pad_trail( const std::string & input );

// Helper functions to drop norm or shape part from systematic variations 
void DropNorm(TH1* hUp,TH1* hDown,TH1* hNom);
void DropShape(TH1* hUp,TH1* hDown,TH1* hNom);

void SetHistoBinsFromOtherHist(TH1* toSet, const TH1* other);

/**
 * Helper function to get the integral of a histogram only considering positive bins
 * @param pinter to histogram
 * @return the effective integral
 */
double EffIntegral(const TH1* const h);
    
/**
  * A helper function that gets the indices of the blinded bins in a region
  * @param the given Region
  * @oaram blinding type
  * @param blinding threshold
  * @return indices of the bins (ROOT index convention)
  */ 
std::vector<int> GetBlindedBins(const Region* reg,
                                const BlindingType type,
                                const double threshold);

/**
  * A helper function to retrive the blinded bins from histograms
  * @param signal histogram
  * @param background histogram
  * @oaram blinding type
  * @param blinding threshold
  * @return blidned bins
  */
std::vector<int> ComputeBlindedBins(const TH1* signal,
                                    const TH1* bkg,
                                    const BlindingType type,
                                    const double threshold);

/**
  * A helper function to combine histograms from a vector of full paths
  * @param A vector where each element represents the full path
  * @return a combined histogram
  */ 
std::unique_ptr<TH1> CombineHistosFromFullPaths(const std::vector<std::string>& paths);

/**
  * A helper function to combine 2D histograms from a vector of full paths
  * @param A vector where each element represents the full path
  * @return a combined 2D histogram
  */ 
std::unique_ptr<TH2> CombineHistos2DFromFullPaths(const std::vector<std::string>& paths);

/**
  * A helper function to calculate error band on ratio
  * @param Graph with total uncertainties
  * @param data
  * @return ratio graph
  */  
std::unique_ptr<TGraphAsymmErrors> GetRatioBand(const TGraphAsymmErrors* total, const TH1D* data);

/**
 * A helepr function to divide bin content by bin width
 * @param histogram
 */
void ScaleByBinWidth(TH1* h);

/**
 * A helepr function to divide bin content by bin width
 * @param raph
 */
void ScaleByBinWidth(TGraphAsymmErrors* g);

/**
  * A helper function to transform an integer into a string, with leading zeros
  * @param input int
  * @param number of places (default = 3)
  * @return string
  */  
std::string IntToFixLenStr(int i,int n=3);

/**
  * Case insensitive string to boolean
  * @param string
  * @return conversion result
  */  
bool StringToBoolean(std::string param);

/**
  * Get apths to files containing a key from a folder
  * @param folder path
  * @param key
  * @param key2
  * @return vector of paths
  */  
std::vector<std::string> GetFilesMatchingString(const std::string& folder, const std::string& key, const std::string& key2);

/**
  * Merge txt files into another one
  * @param paths to the input files
  * @param path to the output file
  * @return outputfile path
  */  
void MergeTxTFiles(const std::vector<std::string>& input, const std::string& out);

/**
  * A helper function to check the validity of a string
  * @param
  * return name
  */
std::string CheckName(const std::string& name);

/**
  * A helper function to remove quotes
  * @param
  * return name
  */
std::string RemoveQuotes(const std::string& name);

/**
  * A helper function to remove spaces
  * @param
  * return name
  */
std::string RemoveSpaces(const std::string& name);

/**
  * A helper function to remove comments
  * @param
  * return name
  */
std::string RemoveComments(const std::string& s);

/**
  * A helper function to vectorise string
  * @param
  * @param splitter
  * @param remove Quotes
  * return splitted string
  */
std::vector<std::string> Vectorize(const std::string& s,char c,bool removeQuotes=true);
}

#endif
