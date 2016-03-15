#include "TtHFitter/NormFactor.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

//__________________________________________________________________________________
//
NormFactor::NormFactor():fName(""),fNominal(0),fMin(0),fMax(0),fConst(false){}

//__________________________________________________________________________________
//
NormFactor::NormFactor(string name, float nominal, float min, float max, bool isConst){
    Set(name,nominal,min,max,isConst);
}

//__________________________________________________________________________________
//
NormFactor::~NormFactor(){}

//__________________________________________________________________________________
//
void NormFactor::Set(string name, float nominal, float min, float max, bool isConst){
    fName = name;
    fNominal = nominal;
    fMin = min;
    fMax = max;
    fConst = isConst;
    //
    fNuisanceParameter = name;
    fTitle = name;
}

//__________________________________________________________________________________
//
void NormFactor::Print(){
    cout << "        NormFactor: " << fName << "\t" << fNominal << ", " << fMin << ", " << fMax;
    if(fConst) cout << " (CONSTANT)"; 
    cout << endl;
}