#include "TtHFitter/ShapeFactor.h"

// -------------------------------------------------------------------------------------------------
// class ShapeFactor

//__________________________________________________________________________________
//
ShapeFactor::ShapeFactor():fName(""),fNominal(0),fMin(0),fMax(0),fConst(false){}

//__________________________________________________________________________________
//
ShapeFactor::ShapeFactor(string name, float nominal, float min, float max, bool isConst){
    Set(name,nominal,min,max,isConst);
}

//__________________________________________________________________________________
//
ShapeFactor::~ShapeFactor(){}

//__________________________________________________________________________________
//
void ShapeFactor::Set(string name, float nominal, float min, float max, bool isConst){
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
void ShapeFactor::Print(){
    cout << "        ShapeFactor: " << fName << "\t" << fNominal << ", " << fMin << ", " << fMax;
    if(fConst) cout << " (CONSTANT)"; 
    cout << endl;
}
