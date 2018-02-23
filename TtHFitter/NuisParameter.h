#ifndef NUISPARAMETER_H
#define NUISPARAMETER_H

#include "TtHFitter/Common.h"

class NuisParameter {
public:
    NuisParameter(std::string name);
    ~NuisParameter();

    std::string fName;
    std::string fTitle;
    std::string fCategory;
    float fStartValue;
    float fFitValue;
    float fPostFitUp; // this should be like +0.8... So alpha+deltaAlpha = fFitValue + fPostFitUp
    float fPostFitDown; // this like -0.7...
    int fConstrainType;
    int fInterpCode;
};

#endif
