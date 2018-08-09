// Class include
#include "TRExFitter/NormFactor.h"

// Framework includes
#include "TRExFitter/StatusLogbook.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

//__________________________________________________________________________________
//
NormFactor::NormFactor():fName(""),fNominal(0),fMin(0),fMax(0),fConst(false){}

//__________________________________________________________________________________
//
NormFactor::NormFactor(std::string name, float nominal, float min, float max, bool isConst, std::string subCategory){
    fName = name;
    fNominal = nominal;
    fMin = min;
    fMax = max;
    fConst = isConst;
    //
    fNuisanceParameter = name;
    fTitle = name;
    //
    fSubCategory = subCategory;
    fExpression = std::make_pair("","");
}

//__________________________________________________________________________________
//
NormFactor::~NormFactor(){}

//__________________________________________________________________________________
//
void NormFactor::Print(){
    if (fConst) WriteInfoStatus("NormFactor::Print", fName + "\t" + std::to_string(fNominal) + ", " + std::to_string(fMin) + ", " + std::to_string(fMax) + "  (CONSTANT)");
    else WriteInfoStatus("NormFactor::Print", fName + "\t" + std::to_string(fNominal) + ", " + std::to_string(fMin) + ", " + std::to_string(fMax));
}
