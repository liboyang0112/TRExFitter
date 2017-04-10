#include "TtHFitter/Common.h"
#include "TtHFitter/NuisParameter.h"
#include "TtHFitter/NormFactor.h"
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
    void DrawNPPulls(const string &path, const string &category, const std::vector < NormFactor* > &normFactors);
    void DrawNormFactors(const string &path, const std::vector < NormFactor* > &normFactor );
    void DrawGammaPulls(const string &path );
    void DrawCorrelationMatrix(string path, const double corrMin = -1. );

    //
    // Data members
    //
    vector<string> fNuisParNames;
    map<string,int> fNuisParIdx;
    map<string,bool> fNuisParIsThere;

    vector<string> fNuisParToHide; // NPs to hide

    std::vector < NuisParameter* > fNuisPar;
    CorrelationMatrix *fCorrMatrix;

};

#endif
