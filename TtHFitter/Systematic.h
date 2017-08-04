#include "TtHFitter/Common.h"

#ifndef __Systematic__
#define __Systematic__

class Systematic {
public:

    enum SystType{
        OVERALL, // 0
        SHAPE, // 1
        HISTO, // 2
        STAT // 3
    };
  
    Systematic(string name,int type=0,float up=0,float down=0);
    Systematic( Systematic &sys);  // copy constructor
    ~Systematic();

    // -------
    // Members
    // -------
  
    string fName;
    string fNuisanceParameter;
    string fTitle;
    string fCategory;
    string fStoredName;
    int fType;
    int fSmoothType;
    bool fPreSmoothing;
    int fSymmetrisationType;
    string fReferenceSample;
    bool fKeepReferenceOverallVar;

    bool fSubtractRefSampleVar;
      
    float fOverallUp;
    float fOverallDown;
    
    float fScaleUp;
    float fScaleDown;

    bool fHasUpVariation;
    bool fHasDownVariation;
    
    bool fIsFreeParameter;
    bool fIsShapeOnly;    
    bool fIsNormOnly;

    std::vector<string> fRegions;
    std::vector<string> fExclude;
    std::vector<std::vector<std::string> > fExcludeRegionSample;
    std::vector<string> fDropShapeIn;
    std::vector<string> fDropNormIn;
    std::vector<string> fKeepNormForSamples;
    std::vector<int> fBins;
    
    // from ntuples - up
    string fWeightUp;
    string fWeightSufUp;  
    std::vector<string> fNtuplePathsUp;
    string fNtuplePathSufUp;
    std::vector<string> fNtupleFilesUp;
    string fNtupleFileSufUp;
    std::vector<string> fNtupleNamesUp;
    string fNtupleNameSufUp;

    // from ntuples - down
    string fWeightDown;
    string fWeightSufDown;  
    std::vector<string> fNtuplePathsDown;
    string fNtuplePathSufDown;
    std::vector<string> fNtupleFilesDown;
    string fNtupleFileSufDown;
    std::vector<string> fNtupleNamesDown;
    string fNtupleNameSufDown;

    string fIgnoreWeight;
    
    // from histos - up
    std::vector<string> fHistoPathsUp;
    string fHistoPathSufUp;
    std::vector<string> fHistoFilesUp;
    string fHistoFileSufUp;
    std::vector<string> fHistoNamesUp;
    string fHistoNameSufUp;
    
    // from histos - down
    std::vector<string> fHistoPathsDown;
    string fHistoPathSufDown;
    std::vector<string> fHistoFilesDown;
    string fHistoFileSufDown;
    std::vector<string> fHistoNamesDown;
    string fHistoNameSufDown;
};

#endif
