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
    float GetNuisParValue(const std::string& p);
    float GetNuisParErrUp(const std::string& p);
    float GetNuisParErrDown(const std::string& p);
    void ReadFromTXT(const std::string& fileName, const std::vector<std::string>& blinded);
    void DrawNPPulls(const std::string &path, const std::string &category, const std::vector < NormFactor* > &normFactors, const std::vector<std::string>& blinded) const;
    void DrawNormFactors(const std::string &path, const std::vector < NormFactor* > &normFactor, const std::vector<std::string>& blinded ) const;
    void DrawGammaPulls(const std::string &path, const std::vector<std::string>& blinded ) const;

    /**
      * Function to draw correlation matrix 
      * @param Path to the output file
      * @param Flag to include gammas on the matrix
      * @param Minimum correlation considered for plotting
      */
    void DrawCorrelationMatrix(const std::string& path, const bool& useGammas, const double corrMin = -1. );
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
