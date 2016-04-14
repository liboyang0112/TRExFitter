#include "TtHFitter/Common.h"

#include "TtHFitter/TtHFit.h"

#ifndef __MultiFit__
#define __MultiFit__

class RooWorkspace;

class MultiFit {
public:
    
    MultiFit(string name="MyMultiFit");
    ~MultiFit();

    void ReadConfigFile(string configFile,string options);
    void AddFitFromConfig(string configFile,string options,string label,string loadSuf="");
    RooWorkspace* CombineWS();
    void SaveCombinedWS();
    std::map < std::string, double > FitCombinedWS( int fitType=1, string inputData="" );
    void GetCombinedLimit(string inputData="obsData"); // or asimovData    
    
    void ComparePOI(string POI);
    void CompareLimit();
    void ComparePulls();

    std::vector< string > fFitNames;
    std::vector< TtHFit* > fFitList;
    std::vector< string > fFitLabels;
    std::vector< string > fFitSuffs;
    
    bool fCombine;
    bool fCompare;
    
    bool fCompareLimits;
    bool fComparePOI;
    bool fComparePulls;
    
    string fName;
    string fLabel;
    bool fShowObserved;
    string fLimitTitle;
    string fPOITitle;
    
    string fPOI;
    float fPOIMin;
    float fPOIMax;
    
    string fLumiLabel;
    string fCmeLabel;
    
    ConfigParser *fConfig;
    
    string fSaveSuf;
    std::vector< bool > fFitShowObserved;
    
    string fDataName;
    int fFitType;
};
    
#endif
