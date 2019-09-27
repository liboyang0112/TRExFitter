#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

/// Framework includes
#include "TRExFitter/HistoTools.h"

/// c++ includes
#include <map>
#include <string>
#include <vector>

class Systematic {
public:

    enum SystType{
         OVERALL, // 0
         SHAPE, // 1
         HISTO, // 2
         STAT // 3
    };

    Systematic(const std::string& name,int type=0,double up=0.,double down=0.);
    Systematic(const Systematic &sys);  // copy constructor
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
    HistoTools::SymmetrizationType fSymmetrisationType;
    std::string fReferenceSample;
    bool fKeepReferenceOverallVar;
    std::string fReferenceSmoothing;
    std::string fReferencePruning;

    bool fSubtractRefSampleVar;

    double fOverallUp;
    double fOverallDown;

    double fScaleUp;
    double fScaleDown;

    std::map<std::string,double> fScaleUpRegions;
    std::map<std::string,double> fScaleDownRegions;


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
    std::vector<std::string> fDummyForSamples;
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

    // needed for systematics on data - like JER
    // up variation
    std::vector<std::string> fNtuplePathsUpRefSample;
    std::string fNtuplePathSufUpRefSample;
    std::vector<std::string> fNtupleFilesUpRefSample;
    std::string fNtupleFileSufUpRefSample;
    std::vector<std::string> fNtupleNamesUpRefSample;
    std::string fNtupleNameSufUpRefSample;

    // from ntuples - down
    std::string fWeightDown;
    std::string fWeightSufDown;
    std::vector<std::string> fNtuplePathsDown;
    std::string fNtuplePathSufDown;
    std::vector<std::string> fNtupleFilesDown;
    std::string fNtupleFileSufDown;
    std::vector<std::string> fNtupleNamesDown;
    std::string fNtupleNameSufDown;

    // needed for systematics on data - like JER
    // Down variation
    std::vector<std::string> fNtuplePathsDownRefSample;
    std::string fNtuplePathSufDownRefSample;
    std::vector<std::string> fNtupleFilesDownRefSample;
    std::string fNtupleFileSufDownRefSample;
    std::vector<std::string> fNtupleNamesDownRefSample;
    std::string fNtupleNameSufDownRefSample;

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
    std::vector<std::string> fHistoPathsUpRefSample;
    std::string fHistoPathSufUpRefSample;
    std::vector<std::string> fHistoFilesUpRefSample;
    std::string fHistoFileSufUpRefSample;
    std::vector<std::string> fHistoNamesUpRefSample;
    std::string fHistoNameSufUpRefSample;

    // from histos - down
    std::vector<std::string> fHistoPathsDown;
    std::string fHistoPathSufDown;
    std::vector<std::string> fHistoFilesDown;
    std::string fHistoFileSufDown;
    std::vector<std::string> fHistoNamesDown;
    std::string fHistoNameSufDown;

    // needed for systematics on data - like JER
    // Down variation
    std::vector<std::string> fHistoPathsDownRefSample;
    std::string fHistoPathSufDownRefSample;
    std::vector<std::string> fHistoFilesDownRefSample;
    std::string fHistoFileSufDownRefSample;
    std::vector<std::string> fHistoNamesDownRefSample;
    std::string fHistoNameSufDownRefSample;

    //
    std::string fSampleUp;
    std::string fSampleDown;

    std::vector<std::string> fSamples;

};

#endif
