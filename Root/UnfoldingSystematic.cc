#include "TRExFitter/UnfoldingSystematic.h"

UnfoldingSystematic::UnfoldingSystematic() :
    fCategory(""),
    fSubCategory(""),
    fHasUpVariation(false),
    fHasDownVariation(false),
    fSampleSmoothing(false),
    fNuisanceParameter(""),
    fName(""),
    fTitle(""),
    fType(0),
    fSymmetrisationType(HistoTools::SymmetrizationType::NOSYMMETRIZATION),
    fSampleSmoothingOption(HistoTools::SmoothOption::MAXVARIATION) 
{
}


