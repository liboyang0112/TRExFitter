// Class include
#include "TRExFitter/Systematic.h"

// -------------------------------------------------------------------------------------------------
// Systematic

//_____________________________________________________________________________
//
Systematic::Systematic(const std::string& name,int type,double up,double down){
    fName = name;
    fTitle = name;
    fNuisanceParameter = name;
    fType = type;
    fCategory = "";
    fSubCategory = "Uncategorised";
    fStoredName = name;

    fSmoothType = 0;
    fSymmetrisationType = HistoTools::NOSYMMETRIZATION;
    fPreSmoothing = false;
    fSampleSmoothing = false;
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
    fReferenceSmoothing = "";
    fReferencePruning = "";
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
    fNtuplePathsUpRefSample.clear();
    fNtuplePathSufUpRefSample = "";
    fNtupleFilesUpRefSample.clear();
    fNtupleFileSufUpRefSample = "";
    fNtupleNamesUpRefSample.clear();
    fNtupleNameSufUpRefSample = "";
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
    fNtuplePathsDownRefSample.clear();
    fNtuplePathSufDownRefSample = "";
    fNtupleFilesDownRefSample.clear();
    fNtupleFileSufDownRefSample = "";
    fNtupleNamesDownRefSample.clear();
    fNtupleNameSufDownRefSample = "";
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
    fHistoPathsUpRefSample.clear();
    fHistoPathSufUpRefSample = "";
    fHistoFilesUpRefSample.clear();
    fHistoFileSufUpRefSample = "";
    fHistoNamesUpRefSample.clear();
    fHistoNameSufUpRefSample = "";
    //
    fHistoPathsDown.clear();
    fHistoPathSufDown = "";
    fHistoFilesDown.clear();
    fHistoFileSufDown = "";
    fHistoNamesDown.clear();
    fHistoNameSufDown = "";
    //
    fHistoPathsDownRefSample.clear();
    fHistoPathSufDownRefSample = "";
    fHistoFilesDownRefSample.clear();
    fHistoFileSufDownRefSample = "";
    fHistoNamesDownRefSample.clear();
    fHistoNameSufDownRefSample = "";
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


Systematic::Systematic(const Systematic &sys) {
    fName = sys.fName;
    fTitle = sys.fTitle;
    fNuisanceParameter = sys.fNuisanceParameter;
    fType = sys.fType;
    fCategory = sys.fCategory;
    fSubCategory = sys.fSubCategory;
    fStoredName = sys.fStoredName;

    fSmoothType = sys.fSmoothType;
    fSymmetrisationType = sys.fSymmetrisationType;
    fPreSmoothing = sys.fPreSmoothing;
    fSampleSmoothing = sys.fSampleSmoothing;
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
    fReferenceSmoothing = sys.fReferenceSmoothing;
    fReferencePruning = sys.fReferencePruning;
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
    fNtuplePathsUpRefSample = sys.fNtuplePathsUpRefSample;
    fNtuplePathSufUpRefSample = sys.fNtuplePathSufUpRefSample;
    fNtupleFilesUpRefSample = sys.fNtupleFilesUpRefSample;
    fNtupleFileSufUpRefSample = sys.fNtupleFileSufUpRefSample;
    fNtupleNamesUpRefSample = sys.fNtupleNamesUpRefSample;
    fNtupleNameSufUpRefSample = sys.fNtupleNameSufUpRefSample;
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
    fNtuplePathsDownRefSample = sys.fNtuplePathsDownRefSample;
    fNtuplePathSufDownRefSample = sys.fNtuplePathSufDownRefSample;
    fNtupleFilesDownRefSample = sys.fNtupleFilesDownRefSample;
    fNtupleFileSufDownRefSample = sys.fNtupleFileSufDownRefSample;
    fNtupleNamesDownRefSample = sys.fNtupleNamesDownRefSample;
    fNtupleNameSufDownRefSample = sys.fNtupleNameSufDownRefSample;
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
    fHistoPathsUpRefSample = sys.fHistoPathsUpRefSample;
    fHistoPathSufUpRefSample = sys.fHistoPathSufUpRefSample;
    fHistoFilesUpRefSample = sys.fHistoFilesUpRefSample;
    fHistoFileSufUpRefSample = sys.fHistoFileSufUpRefSample;
    fHistoNamesUpRefSample = sys.fHistoNamesUpRefSample;
    fHistoNameSufUpRefSample = sys.fHistoNameSufUpRefSample;
    //
    fHistoPathsDown = sys.fHistoPathsDown;
    fHistoPathSufDown = sys.fHistoPathSufDown;
    fHistoFilesDown = sys.fHistoFilesDown;
    fHistoFileSufDown = sys.fHistoFileSufDown;
    fHistoNamesDown = sys.fHistoNamesDown;
    fHistoNameSufDown = sys.fHistoNameSufDown;
    //
    fHistoPathsDownRefSample = sys.fHistoPathsDownRefSample;
    fHistoPathSufDownRefSample = sys.fHistoPathSufDownRefSample;
    fHistoFilesDownRefSample = sys.fHistoFilesDownRefSample;
    fHistoFileSufDownRefSample = sys.fHistoFileSufDownRefSample;
    fHistoNamesDownRefSample = sys.fHistoNamesDownRefSample;
    fHistoNameSufDownRefSample = sys.fHistoNameSufDownRefSample;
    //
    fRegions = sys.fRegions;
    fExclude = sys.fExclude;
    fExcludeRegionSample = sys.fExcludeRegionSample;
    fDropShapeIn = sys.fDropShapeIn;
    fDropNormIn = sys.fDropNormIn;
    //
    fSampleUp = "";
    fSampleDown = "";

    fCombineName = "";
    fCombineType = COMBINATIONTYPE::ENVELOPE;
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
    fNtuplePathsUpRefSample.clear();
    fNtupleFilesUpRefSample.clear();
    fNtupleNamesUpRefSample.clear();
    fNtuplePathsDownRefSample.clear();
    fNtupleFilesDownRefSample.clear();
    fNtupleNamesDownRefSample.clear();
    fHistoPathsUp.clear();
    fHistoFilesUp.clear();
    fHistoNamesUp.clear();
    fHistoPathsUpRefSample.clear();
    fHistoFilesUpRefSample.clear();
    fHistoNamesUpRefSample.clear();
    fHistoPathsDown.clear();
    fHistoFilesDown.clear();
    fHistoNamesDown.clear();
    fHistoPathsDownRefSample.clear();
    fHistoFilesDownRefSample.clear();
    fHistoNamesDownRefSample.clear();
}
