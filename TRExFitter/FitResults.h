#ifndef FITRESULTS_H
#define FITRESULTS_H

/// c++ includes
#include <map>
#include <string>
#include <vector>

/// Forward class declaration
class CorrelationMatrix;
class NormFactor;
class NuisParameter;

class FitResults {
public:
    FitResults();
    ~FitResults();

    //
    // Functions
    //
    void AddNuisPar(NuisParameter *par);
    float GetNuisParValue(std::string p);
    float GetNuisParErrUp(std::string p);
    float GetNuisParErrDown(std::string p);
    void ReadFromTXT(std::string fileName);
    void DrawNPPulls(const std::string &path, const std::string &category, const std::vector < NormFactor* > &normFactors);
    void DrawNormFactors(const std::string &path, const std::vector < NormFactor* > &normFactor );
    void DrawGammaPulls(const std::string &path );
    void DrawCorrelationMatrix(std::string path, const double corrMin = -1. );
    void SetPOIPrecision(const int& precision){fPOIPrecision = precision;}

    //
    // Data members
    //
    std::vector<std::string> fNuisParNames;
    std::map<std::string,int> fNuisParIdx;
    std::map<std::string,bool> fNuisParIsThere;

    std::vector<std::string> fNuisParToHide; // NPs to hide

    std::vector < NuisParameter* > fNuisPar;
    CorrelationMatrix *fCorrMatrix;

    int fPOIPrecision;

};

#endif
