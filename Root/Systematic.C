#include "TtHFitter/Systematic.h"

// -------------------------------------------------------------------------------------------------
// Systematic

//_____________________________________________________________________________
//
Systematic::Systematic(string name,int type,float up,float down){
    fName = name;
    fTitle = name;
    fNuisanceParameter = name;
    fType = type;
    fCategory = "";
    
    fSmoothType = 0;
    fSymmetrisationType = 0;
    //
    fOverallUp = up;
    fOverallDown = down;
    //
    fHasUpVariation = true;
    fHasDownVariation = true;
    //
    fIsFreeParameter = false;
    //
    fReferenceSample = "";
    fKeepReferenceOverallVar = true;
    //
    fWeightUp = "";
    fWeightSufUp = "";
    fNtuplePathsUp.clear();
    fNtuplePathSufUp = "";
    fNtupleFilesUp.clear();
    fNtupleFileSufUp = "";
    fNtupleNamesUp.clear();
    fNtupleNameSufUp = "";
    //
    fWeightDown = "";
    fWeightSufDown = "";  
    fNtuplePathsDown.clear();
    fNtuplePathSufDown = "";
    fNtupleFilesDown.clear();
    fNtupleFileSufDown = "";
    fNtupleNamesDown.clear();
    fNtupleNameSufDown = "";
    //
    fIgnoreWeight = "";
    //
    fHistoPathsUp.clear();
    fHistoPathSufUp = "";
    fHistoFilesUp.clear();
    fHistoFileSufUp = "";
    fHistoNamesUp.clear();
    fHistoNameSufUp = "";
    //
    fHistoPathsDown.clear();
    fHistoPathSufDown = "";
    fHistoFilesDown.clear();
    fHistoFileSufDown = "";
    fHistoNamesDown.clear();
    fHistoNameSufDown = "";
    //
    fRegions.clear();
    fExclude.clear();
    fDropShapeIn.clear();
}

//_____________________________________________________________________________
//
Systematic::~Systematic(){
    fNtuplePathsUp.clear();
    fNtupleFilesUp.clear();
    fNtupleNamesUp.clear();
    fNtuplePathsDown.clear();
    fNtupleFilesDown.clear();
    fNtupleNamesDown.clear();
    fHistoPathsUp.clear();
    fHistoFilesUp.clear();
    fHistoNamesUp.clear(); 
    fHistoPathsDown.clear();
    fHistoFilesDown.clear();
    fHistoNamesDown.clear();
}
