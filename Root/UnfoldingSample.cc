#include "TRExFitter/UnfoldingSample.h"

#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"

UnfoldingSample::UnfoldingSample() :
    fName(""),
    fTitle(""),
    fFillColor(0),
    fLineColor(0),
    fHasResponse(false),
    fHasAcceptance(false)
{
}

std::vector<std::shared_ptr<Sample> > UnfoldingSample::ConvertToSample(const Region* reg,
                                                                       const int bins,
                                                                       const std::string& name) const {

    std::vector<std::shared_ptr<Sample> > result;
    for (int ibin = 0; ibin < bins; ++ibin) {
        const std::string sampleName = reg->fName + "_Truth_bin_" + std::to_string(ibin+1);

        std::shared_ptr<Sample> sample = std::make_shared<Sample>(sampleName, Sample::SampleType::SIGNAL);
        sample->SetTitle(fTitle); 
        sample->SetFillColor(fFillColor); 
        sample->SetLineColor(fLineColor); 
        sample->fRegions = Common::ToVec(reg->fName);

        // paths
        sample->fHistoPaths = Common::ToVec(name + "/UnfoldingHistograms");
        sample->fHistoFiles = Common::ToVec("FoldedHistograms");
        const std::string histoName = "nominal/" + reg->fName + "_" + fName + "_bin_" + std::to_string(ibin);
        sample->fHistoNames = Common::ToVec(histoName);
    
        result.emplace_back(sample);
    }

    return result;
}
