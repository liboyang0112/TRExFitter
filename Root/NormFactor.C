// Class include
#include "TRExFitter/NormFactor.h"

// Framework includes
#include "TRExFitter/StatusLogbook.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

//__________________________________________________________________________________
//
NormFactor::NormFactor() : 
    fName(""),
    fNuisanceParameter(""),
    fTitle(""),
    fCategory(""),
    fSubCategory(""),
    fNominal(0),
    fMin(0),
    fMax(0),
    fConst(false),
    fTau(0) {
}

//__________________________________________________________________________________
//
NormFactor::NormFactor(const std::string& name, double nominal, double min, double max, bool isConst, std::string subCategory) : 
    fName(name),
    fNuisanceParameter(name),
    fTitle(name),
    fCategory(""),
    fSubCategory(subCategory),
    fNominal(nominal),
    fMin(min),
    fMax(max),
    fConst(isConst),
    fTau(0) {
}

//__________________________________________________________________________________
//
void NormFactor::Print() const{
    if (fConst) WriteInfoStatus("NormFactor::Print", fName + "\t" + std::to_string(fNominal) + ", " + std::to_string(fMin) + ", " + std::to_string(fMax) + "  (CONSTANT)");
    else WriteInfoStatus("NormFactor::Print", fName + "\t" + std::to_string(fNominal) + ", " + std::to_string(fMin) + ", " + std::to_string(fMax));
}
