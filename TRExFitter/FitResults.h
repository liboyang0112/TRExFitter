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
    float GetNuisParValue(const std::string& p) const;
    float GetNuisParErrUp(const std::string& p) const;
    float GetNuisParErrDown(const std::string& p) const;
    void ReadFromTXT(const std::string& fileName);
    void DrawNPPulls(const std::string &path, const std::string &category, const std::vector < NormFactor* > &normFactors) const;
    void DrawNormFactors(const std::string &path, const std::vector < NormFactor* > &normFactor ) const;
    void DrawGammaPulls(const std::string &path ) const;
    void DrawCorrelationMatrix(const std::string& path, const double corrMin = -1. ) const;
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
