#include "TtHFitter/Common.h"

#include "TtHFitter/TtHFit.h"

#ifndef __MultiFit__
#define __MultiFit__

class MultiFit {
public:
    
    MultiFit(string name="MyMultiFit");
    ~MultiFit();

    void ReadConfigFile(string configFile,string options);
    void AddFitFromConfig(string configFile,string options,string label,string loadSuf="");
    void ComparePOI(string POI);
    void CompareLimit();
    void ComparePulls();

    std::vector< TtHFit* > fFitList;
    std::vector< string > fFitLabels;
    std::vector< string > fFitSuffs;
    
    bool fCompareLimits;
    bool fComparePOI;
    bool fComparePulls;
    
    string fName;
    string fLabel;
    bool fShowObserved;
    string fLimitTitle;
    string fPOITitle;
    
    string fLumiLabel;
    string fCmeLabel;
    
    ConfigParser *fConfig;
    
    string fSaveSuf;
    std::vector< bool > fFitShowObserved;
};
    
#endif
