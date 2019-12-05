#ifndef NORMFACTOR_H
#define NORMFACTOR_H

/// c++ includes
#include <string>
#include <vector>

class NormFactor{
public:
    explicit NormFactor();
    
    explicit NormFactor(const std::string& name, double nominal=1., double min=0., double max=10., bool isConst=false, std::string subCategory="NormFactors");
    
    ~NormFactor() = default;
    NormFactor(const NormFactor& n) = delete;
    NormFactor(NormFactor&& n) = delete;
    NormFactor& operator=(const NormFactor& n) = delete;
    NormFactor& operator=(NormFactor&& n) = delete;

    void Print() const;

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;
    std::string fSubCategory;

    double fNominal;
    double fMin;
    double fMax;
    bool fConst;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;

    std::pair<std::string,std::string> fExpression;
    
    double fTau; // for Tikhonov regularization
};

#endif
