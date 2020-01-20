#include "TRExFitter/UnfoldingSystematic.h"

#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/Systematic.h"

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

    std::vector<std::string> fResponseMatrixPathsUp;
    std::vector<std::string> fResponseMatrixPathsDown;
    std::vector<std::string> fResponseMatrixNamesUp;
    std::vector<std::string> fResponseMatrixNamesDown;
    std::vector<std::string> fResponseMatrixFilesUp;
    std::vector<std::string> fResponseMatrixFilesDown;
    std::vector<std::string> fResponseMatrixPathSuffsUp;
    std::vector<std::string> fResponseMatrixPathSuffsDown;
    std::vector<std::string> fResponseMatrixNameSuffsUp;
    std::vector<std::string> fResponseMatrixNameSuffsDown;
    std::vector<std::string> fResponseMatrixFileSuffsUp;
    std::vector<std::string> fResponseMatrixFileSuffsDown;

    bool fSampleSmoothing;
    HistoTools::SymmetrizationType fSymmetrisationType;
    HistoTools::SmoothOption fSampleSmoothingOption;


std::vector<Systematic*> UnfoldingSystematic::ConvertToSystematic(const Region* reg,
                                                                  const int bins,
                                                                  const std::string& name,
                                                                  std::vector<Sample*>& samples) const {
    std::vector<Systematic*> result;
    for (int ibin = 0; ibin < bins; ++ibin) {
        const std::string sampleName = "Truth_bin_" + std::to_string(ibin+1);

        Systematic* syst = new Systematic(fName, fType);
        syst->fNuisanceParameter = fNuisanceParameter;
        syst->fTitle = fTitle;
        syst->fHasUpVariation = fHasUpVariation;
        syst->fHasDownVariation = fHasDownVariation;
        syst->fRegions = Common::ToVec(reg->fName);
        syst->fSamples = Common::ToVec(sampleName);
        syst->fCategory = fCategory;
        syst->fSubCategory = fSubCategory;
        syst->fSampleSmoothing = fSampleSmoothing;
        syst->fSymmetrisationType = fSymmetrisationType;
        syst->fSampleSmoothOption = fSampleSmoothingOption;
        
        TRExFitter::SYSTMAP[syst->fName] = syst->fTitle;

        // Paths
        if (fHasUpVariation) {
            syst->fHistoPathsUp = Common::ToVec(name + "/UnfoldingHistograms");
            syst->fHistoFilesUp = Common::ToVec("FoldedHistograms");
            const std::string histoName = fName + "_Up/" + reg->fName + "_" + fName + "_bin_" + std::to_string(ibin);
            syst->fHistoNamesUp = Common::ToVec(histoName);
        }
        if (fHasDownVariation) {
            syst->fHistoPathsDown = Common::ToVec(name + "/UnfoldingHistograms");
            syst->fHistoFilesDown = Common::ToVec("FoldedHistograms");
            const std::string histoName = fName + "_Down/" + reg->fName + "_" + fName + "_bin_" + std::to_string(ibin);
            syst->fHistoNamesDown = Common::ToVec(histoName);
        }
        
        if (syst->fHistoPathsUpRefSample.size()     == 0) syst->fHistoPathsUpRefSample     = syst->fHistoPathsUp;
        if (syst->fHistoPathsDownRefSample.size()   == 0) syst->fHistoPathsDownRefSample   = syst->fHistoPathsDown;
        if (syst->fHistoPathSufUpRefSample.size()   == 0) syst->fHistoPathSufUpRefSample   = syst->fHistoPathSufUp;
        if (syst->fHistoPathSufDownRefSample.size() == 0) syst->fHistoPathSufDownRefSample = syst->fHistoPathSufDown;
        if (syst->fHistoFilesUpRefSample.size()     == 0) syst->fHistoFilesUpRefSample     = syst->fHistoFilesUp;
        if (syst->fHistoFilesDownRefSample.size()   == 0) syst->fHistoFilesDownRefSample   = syst->fHistoFilesDown;
        if (syst->fHistoFileSufUpRefSample.size()   == 0) syst->fHistoFileSufUpRefSample   = syst->fHistoFileSufUp;
        if (syst->fHistoFileSufDownRefSample.size() == 0) syst->fHistoFileSufDownRefSample = syst->fHistoFileSufDown;
        if (syst->fHistoNamesUpRefSample.size()     == 0) syst->fHistoNamesUpRefSample     = syst->fHistoNamesUp;
        if (syst->fHistoNamesDownRefSample.size()   == 0) syst->fHistoNamesDownRefSample   = syst->fHistoNamesDown;
        if (syst->fHistoNameSufUpRefSample.size()   == 0) syst->fHistoNameSufUpRefSample   = syst->fHistoNameSufUp;
        if (syst->fHistoNameSufDownRefSample.size() == 0) syst->fHistoNameSufDownRefSample = syst->fHistoNameSufDown;


        for(auto isample : samples) {
            if(!isample->fUseSystematics) continue;
            if(Common::FindInStringVector(syst->fSamples, isample->fName) < 0) continue;
            
            isample->AddSystematic(syst);
        }

        result.emplace_back(syst);
    }

    return result;
}
