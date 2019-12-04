// Class include
#include "TRExFitter/Systematic.h"

// -------------------------------------------------------------------------------------------------
// Systematic

//_____________________________________________________________________________
//
Systematic::Systematic(const std::string& name,int type,double up,double down) :
    fName(name),
    fNuisanceParameter(name),
    fTitle(name),
    fCategory(""),
    fSubCategory("Uncategorised"),
    fStoredName(name),
    fType(type),
    fSmoothType(0),
    fPreSmoothing(false),
    fSampleSmoothing(false),
    fSymmetrisationType(HistoTools::NOSYMMETRIZATION),
    fReferenceSample(""),
    fKeepReferenceOverallVar(true),
    fReferenceSmoothing(""),
    fReferencePruning(""),
    fSubtractRefSampleVar(false),
    fOverallUp(up),
    fOverallDown(down),
    fScaleUp(1.),
    fScaleDown(1.),
    fHasUpVariation(true),
    fHasDownVariation(true),
    fIsFreeParameter(false),
    fIsShapeOnly(false),
    fIsNormOnly(false),
    fWeightUp(""),
    fWeightSufUp(""),
    fNtuplePathSufUp(""),
    fNtupleFileSufUp(""),
    fNtupleNameSufUp(""),
    fNtuplePathSufUpRefSample(""),
    fNtupleFileSufUpRefSample(""),
    fNtupleNameSufUpRefSample(""),
    fWeightDown(""),
    fWeightSufDown(""),
    fNtuplePathSufDown(""),
    fNtupleFileSufDown(""),
    fNtupleNameSufDown(""),
    fNtuplePathSufDownRefSample(""),
    fNtupleFileSufDownRefSample(""),
    fNtupleNameSufDownRefSample(""),
    fIgnoreWeight(""),
    fHistoPathSufUp(""),
    fHistoFileSufUp(""),
    fHistoNameSufUp(""),
    fHistoPathSufUpRefSample(""),
    fHistoFileSufUpRefSample(""),
    fHistoNameSufUpRefSample(""),
    fHistoPathSufDown(""),
    fHistoFileSufDown(""),
    fHistoNameSufDown(""),
    fHistoPathSufDownRefSample(""),
    fHistoFileSufDownRefSample(""),
    fHistoNameSufDownRefSample(""),
    fSampleUp(""),
    fSampleDown("") {
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
    fScaleUp = sys.fScaleUp;
    fScaleDown = sys.fScaleDown;
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
    fSamples = sys.fSamples;
    fExcludeRegionSample = sys.fExcludeRegionSample;
    fDropShapeIn = sys.fDropShapeIn;
    fDropNormIn = sys.fDropNormIn;
    //
    fSampleUp = "";
    fSampleDown = "";

    fCombineName = "";
    fCombineType = COMBINATIONTYPE::ENVELOPE;
}
