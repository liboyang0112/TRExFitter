#ifndef NORMFACTOR_H
#define NORMFACTOR_H

#include "TtHFitter/Common.h"

class NormFactor{
public:
    NormFactor();
    NormFactor(std::string name, float nominal=1, float min=0, float max=10, bool isConst=false);
    ~NormFactor();
    void Set(std::string name, float nominal=1, float min=0, float max=10, bool isConst=false);

    void Print();

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;

    float fNominal;
    float fMin;
    float fMax;
    bool fConst;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;
    
    std::pair<std::string,std::string> fExpression;
};

#endif
