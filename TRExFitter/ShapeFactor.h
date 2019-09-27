#ifndef SHAPEFACTOR_H
#define SHAPEFACTOR_H

/// c++ includes
#include <string>
#include <vector>

class ShapeFactor{
public:
    ShapeFactor();
    ShapeFactor(const std::string& name, double nominal=1., double min=0., double max=10., bool isConst=false);
    ~ShapeFactor();
    void Set(const std::string& name, double nominal=1., double min=0., double max=10., bool isConst=false);

    void Print() const;

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;

    double fNominal;
    double fMin;
    double fMax;
    bool fConst;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;

    int fNbins;
};

#endif
