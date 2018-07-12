#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#include "TtHFitter/Common.h"

class Systematic {
public:

    enum SystType{
        OVERALL, // 0
        SHAPE, // 1
        HISTO, // 2
        STAT // 3
    };

    Systematic(std::string name,int type=0,float up=0,float down=0);
    Systematic( Systematic &sys);  // copy constructor
    ~Systematic();

    // -------
    // Members
    // -------

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;
    std::string fSubCategory;
    std::string fStoredName;
    int fType;
    int fSmoothType;
    bool fPreSmoothing;
    int fSymmetrisationType;
    std::string fReferenceSample;
    bool fKeepReferenceOverallVar;

    bool fSubtractRefSampleVar;

    float fOverallUp;
    float fOverallDown;

    float fScaleUp;
    float fScaleDown;

    std::map<std::string,float> fScaleUpRegions;
    std::map<std::string,float> fScaleDownRegions;


    bool fHasUpVariation;
    bool fHasDownVariation;

    bool fIsFreeParameter;
    bool fIsShapeOnly;
    bool fIsNormOnly;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;
    std::vector<std::vector<std::string> > fExcludeRegionSample;
    std::vector<std::string> fDropShapeIn;
    std::vector<std::string> fDropNormIn;
    std::vector<std::string> fKeepNormForSamples;
    std::vector<int> fBins;

    // from ntuples - up
    std::string fWeightUp;
    std::string fWeightSufUp;
    std::vector<std::string> fNtuplePathsUp;
    std::string fNtuplePathSufUp;
    std::vector<std::string> fNtupleFilesUp;
    std::string fNtupleFileSufUp;
    std::vector<std::string> fNtupleNamesUp;
    std::string fNtupleNameSufUp;

    // from ntuples - down
    std::string fWeightDown;
    std::string fWeightSufDown;
    std::vector<std::string> fNtuplePathsDown;
    std::string fNtuplePathSufDown;
    std::vector<std::string> fNtupleFilesDown;
    std::string fNtupleFileSufDown;
    std::vector<std::string> fNtupleNamesDown;
    std::string fNtupleNameSufDown;

    std::string fIgnoreWeight;

    // from histos - up
    std::vector<std::string> fHistoPathsUp;
    std::string fHistoPathSufUp;
    std::vector<std::string> fHistoFilesUp;
    std::string fHistoFileSufUp;
    std::vector<std::string> fHistoNamesUp;
    std::string fHistoNameSufUp;

    // needed for systematics on data - like JER
    // up variation
    std::vector<std::string> fHistoPathsUpData;
    std::string fHistoPathSufUpData;
    std::vector<std::string> fHistoFilesUpData;
    std::string fHistoFileSufUpData;
    std::vector<std::string> fHistoNamesUpData;
    std::string fHistoNameSufUpData;

    // from histos - down
    std::vector<std::string> fHistoPathsDown;
    std::string fHistoPathSufDown;
    std::vector<std::string> fHistoFilesDown;
    std::string fHistoFileSufDown;
    std::vector<std::string> fHistoNamesDown;
    std::string fHistoNameSufDown;
    
    // needed for systematics on data - like JER
    // Down variation
    std::vector<std::string> fHistoPathsDownData;
    std::string fHistoPathSufDownData;
    std::vector<std::string> fHistoFilesDownData;
    std::string fHistoFileSufDownData;
    std::vector<std::string> fHistoNamesDownData;
    std::string fHistoNameSufDownData;

    //
    std::string fSampleUp;
    std::string fSampleDown;

};

#endif
