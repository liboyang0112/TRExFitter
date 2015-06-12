#include "TtHFitter/NuisParameter.h"

//__________________________________________________________________________________
//
NuisParameter::NuisParameter(string name){
    fName = name;
    fTitle = name;
    fStartValue = 0;
    fFitValue = 0;
    fPostFitUp = 1;
    fPostFitDown = -1;
    fConstrainType = 0; // 0 -> Gaussian
    fInterpCode = 4; // 4 -> Pol interp. + Exp extrap.
}

//__________________________________________________________________________________
//
NuisParameter::~NuisParameter(){
}
