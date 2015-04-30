#include "TtHFitter/NormFactor.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

//__________________________________________________________________________________
//
NormFactor::NormFactor():fName(""),fNominal(0),fMin(0),fMax(0){}

//__________________________________________________________________________________
//
NormFactor::NormFactor(string name, float nominal, float min, float max){
    Set(name,nominal,min,max);
}

//__________________________________________________________________________________
//
NormFactor::~NormFactor(){}

//__________________________________________________________________________________
//
void NormFactor::Set(string name, float nominal, float min, float max){
    fName = name;
    fNominal = nominal;
    fMin = min;
    fMax = max;
}

//__________________________________________________________________________________
//
void NormFactor::Print(){
    cout << "        NormFactor: " << fName << "\t" << fNominal << ", " << fMin << ", " << fMax << endl;
}