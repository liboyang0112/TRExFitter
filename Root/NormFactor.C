#include "TtHFitter/NormFactor.h"
#include "TtHFitter/StatusLogbook.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

//__________________________________________________________________________________
//
NormFactor::NormFactor():fName(""),fNominal(0),fMin(0),fMax(0),fConst(false){}

//__________________________________________________________________________________
//
NormFactor::NormFactor(std::string name, float nominal, float min, float max, bool isConst){
    Set(name,nominal,min,max,isConst);
}

//__________________________________________________________________________________
//
NormFactor::~NormFactor(){}

//__________________________________________________________________________________
//
void NormFactor::Set(std::string name, float nominal, float min, float max, bool isConst){
    fName = name;
    fNominal = nominal;
    fMin = min;
    fMax = max;
    fConst = isConst;
    //
    fNuisanceParameter = name;
    fTitle = name;
    //
    fExpression = std::make_pair("","");
}

//__________________________________________________________________________________
//
void NormFactor::Print(){
    if (fConst) WriteInfoStatus("NormFactor::Print", fName + "\t" + std::to_string(fNominal) + ", " + std::to_string(fMin) + ", " + std::to_string(fMax) + "  (CONSTANT)");
    else WriteInfoStatus("NormFactor::Print", fName + "\t" + std::to_string(fNominal) + ", " + std::to_string(fMin) + ", " + std::to_string(fMax));
}
