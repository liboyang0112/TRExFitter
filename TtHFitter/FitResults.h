#include "TtHFitter/Common.h"
#include "TtHFitter/NuisParameter.h"
#include "TtHFitter/CorrelationMatrix.h"

#ifndef __FitResults__
#define __FitResults__

class FitResults {
public:
    FitResults();
    ~FitResults();
    
    //
    // Functions
    //
    void AddNuisPar(NuisParameter *par);
    float GetNuisParValue(string p);
    float GetNuisParErrUp(string p);
    float GetNuisParErrDown(string p);
    void ReadFromTXT(string fileName);
    void DrawPulls(string path);
    void DrawCorrelationMatrix( const string &path, const double corrMin = -1. );
    
    //
    // Data members
    //
    vector<string> fNuisParNames;
    map<string,int> fNuisParIdx;
    map<string,bool> fNuisParIsThere;
    
    std::vector < NuisParameter* > fNuisPar;
    CorrelationMatrix *fCorrMatrix;
    
};

#endif
