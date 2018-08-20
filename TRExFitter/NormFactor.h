#ifndef NORMFACTOR_H
#define NORMFACTOR_H

/// c++ includes
#include <string>
#include <vector>

class NormFactor{
public:
    NormFactor();
    NormFactor(const std::string& name, float nominal=1, float min=0, float max=10, bool isConst=false, std::string subCategory="NormFactors");
    ~NormFactor();

    void Print();

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;
    std::string fSubCategory;

    float fNominal;
    float fMin;
    float fMax;
    bool fConst;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;
    
    std::pair<std::string,std::string> fExpression;
};

#endif
