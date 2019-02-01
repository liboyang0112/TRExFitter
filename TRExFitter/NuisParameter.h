#ifndef NUISPARAMETER_H
#define NUISPARAMETER_H

/// c++ includes
#include <string>

class NuisParameter {
public:
    NuisParameter(const std::string& name);
    ~NuisParameter();

    std::string fName;
    std::string fTitle;
    std::string fCategory;
    double fStartValue;
    double fFitValue;
    double fPostFitUp; // this should be like +0.8... So alpha+deltaAlpha = fFitValue + fPostFitUp
    double fPostFitDown; // this like -0.7...
    int fConstrainType;
    int fInterpCode;
};

#endif
