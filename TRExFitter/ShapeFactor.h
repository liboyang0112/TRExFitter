#ifndef SHAPEFACTOR_H
#define SHAPEFACTOR_H

/// c++ includes
#include <string>
#include <vector>

class ShapeFactor{
public:
    ShapeFactor();
    ShapeFactor(std::string name, float nominal=1, float min=0, float max=10, bool isConst=false);
    ~ShapeFactor();
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
};

#endif
