#ifndef FITRESULTS_H
#define FITRESULTS_H

/// c++ includes
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/// Forward class declaration
class CorrelationMatrix;
class NormFactor;
class NuisParameter;

class FitResults {
public:
    explicit FitResults();

    ~FitResults();

    FitResults(const FitResults& f) = delete;
    FitResults(FitResults&& f) = delete;
    FitResults& operator=(const FitResults& f) = delete;
    FitResults& operator=(FitResults&& f) = delete;

    //
    // Functions
    //
    void AddNuisPar(NuisParameter *par);
    double GetNuisParValue(const std::string& p);
    double GetNuisParErrUp(const std::string& p);
    double GetNuisParErrDown(const std::string& p);
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
    std::unordered_map<std::string,std::size_t> fNuisParIdx;
    std::unordered_map<std::string,bool> fNuisParIsThere;

    std::vector<std::string> fNuisParToHide; // NPs to hide

    std::vector < std::unique_ptr<NuisParameter> > fNuisPar;
    std::unique_ptr<CorrelationMatrix> fCorrMatrix;

    int fPOIPrecision;

    double fNLL;

};

#endif
