#include "TtHFitter/Systematic.h"

using namespace std;

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
    fSubCategory = "Uncategorised";
    fStoredName = name;

    fSmoothType = 0;
    fSymmetrisationType = 0;
    fPreSmoothing = false;
    //
    fOverallUp   = up;
    fOverallDown = down;
    //
    fScaleUp   = 1.;
    fScaleDown = 1.;
    //
    fScaleUpRegions  .clear();
    fScaleDownRegions.clear();
    //
    fHasUpVariation = true;
    fHasDownVariation = true;
    //
    fIsFreeParameter = false;
    fIsShapeOnly = false;
    fIsNormOnly = false;
    //
    fReferenceSample = "";
    fKeepReferenceOverallVar = true;
    //
    fSubtractRefSampleVar = false;
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
    fHistoPathsUpData.clear();
    fHistoPathSufUpData = "";
    fHistoFilesUpData.clear();
    fHistoFileSufUpData = "";
    fHistoNamesUpData.clear();
    fHistoNameSufUpData = "";
    //
    fHistoPathsDown.clear();
    fHistoPathSufDown = "";
    fHistoFilesDown.clear();
    fHistoFileSufDown = "";
    fHistoNamesDown.clear();
    fHistoNameSufDown = "";
    //
    fHistoPathsDownData.clear();
    fHistoPathSufDownData = "";
    fHistoFilesDownData.clear();
    fHistoFileSufDownData = "";
    fHistoNamesDownData.clear();
    fHistoNameSufDownData = "";
    //
    fRegions.clear();
    fExclude.clear();
    fExcludeRegionSample.clear();
    fDropShapeIn.clear();
    fDropNormIn.clear();
    //
    fSampleUp = "";
    fSampleDown = "";
}


Systematic::Systematic(Systematic &sys) {
    fName = sys.fName;
    fTitle = sys.fTitle;
    fNuisanceParameter = sys.fNuisanceParameter;
    fType = sys.fType;
    fCategory = sys.fCategory;
    fSubCategory = sys.fSubCategory;
    fStoredName = sys.fStoredName;

    fSmoothType = sys.fSmoothType;
    fSymmetrisationType = sys.fSymmetrisationType;
    //
    fOverallUp = sys.fOverallUp;
    fOverallDown = sys.fOverallDown;
    //
    fHasUpVariation = sys.fHasUpVariation;
    fHasDownVariation = sys.fHasDownVariation;
    //
    fIsFreeParameter = sys.fIsFreeParameter;
    fIsShapeOnly = sys.fIsShapeOnly;
    fIsNormOnly = sys.fIsNormOnly;
    //
    fReferenceSample = sys.fReferenceSample;
    fKeepReferenceOverallVar = sys.fKeepReferenceOverallVar;
    //
    fSubtractRefSampleVar = sys.fSubtractRefSampleVar;
    //
    fWeightUp = sys.fWeightUp;
    fWeightSufUp = sys.fWeightSufUp;
    fNtuplePathsUp = sys.fNtuplePathsUp;
    fNtuplePathSufUp = sys.fNtuplePathSufUp;
    fNtupleFilesUp = sys.fNtupleFilesUp;
    fNtupleFileSufUp = sys.fNtupleFileSufUp;
    fNtupleNamesUp = sys.fNtupleNamesUp;
    fNtupleNameSufUp = sys.fNtupleNameSufUp;
    //
    fWeightDown = sys.fWeightDown;
    fWeightSufDown = sys.fWeightSufDown;
    fNtuplePathsDown = sys.fNtuplePathsDown;
    fNtuplePathSufDown = sys.fNtuplePathSufDown;
    fNtupleFilesDown = sys.fNtupleFilesDown;
    fNtupleFileSufDown = sys.fNtupleFileSufDown;
    fNtupleNamesDown = sys.fNtupleNamesDown;
    fNtupleNameSufDown = sys.fNtupleNameSufDown;
    //
    fIgnoreWeight = sys.fIgnoreWeight;
    //
    fHistoPathsUp = sys.fHistoPathsUp;
    fHistoPathSufUp = sys.fHistoPathSufUp;
    fHistoFilesUp = sys.fHistoFilesUp;
    fHistoFileSufUp = sys.fHistoFileSufUp;
    fHistoNamesUp = sys.fHistoNamesUp;
    fHistoNameSufUp = sys.fHistoNameSufUp;
    //
    fHistoPathsUpData = sys.fHistoPathsUpData;
    fHistoPathSufUpData = sys.fHistoPathSufUpData;
    fHistoFilesUpData = sys.fHistoFilesUpData;
    fHistoFileSufUpData = sys.fHistoFileSufUpData;
    fHistoNamesUpData = sys.fHistoNamesUpData;
    fHistoNameSufUpData = sys.fHistoNameSufUpData;
    //
    fHistoPathsDown = sys.fHistoPathsDown;
    fHistoPathSufDown = sys.fHistoPathSufDown;
    fHistoFilesDown = sys.fHistoFilesDown;
    fHistoFileSufDown = sys.fHistoFileSufDown;
    fHistoNamesDown = sys.fHistoNamesDown;
    fHistoNameSufDown = sys.fHistoNameSufDown;
    //
    fHistoPathsDownData = sys.fHistoPathsDownData;
    fHistoPathSufDownData = sys.fHistoPathSufDownData;
    fHistoFilesDownData = sys.fHistoFilesDownData;
    fHistoFileSufDownData = sys.fHistoFileSufDownData;
    fHistoNamesDownData = sys.fHistoNamesDownData;
    fHistoNameSufDownData = sys.fHistoNameSufDownData;
    //
    fRegions = sys.fRegions;
    fExclude = sys.fExclude;
    fExcludeRegionSample = sys.fExcludeRegionSample;
    fDropShapeIn = sys.fDropShapeIn;
    fDropNormIn = sys.fDropNormIn;
    //
    fSampleUp = "";
    fSampleDown = "";
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
    fHistoPathsUpData.clear();
    fHistoFilesUpData.clear();
    fHistoNamesUpData.clear();
    fHistoPathsDown.clear();
    fHistoFilesDown.clear();
    fHistoNamesDown.clear();
    fHistoPathsDownData.clear();
    fHistoFilesDownData.clear();
    fHistoNamesDownData.clear();
}
