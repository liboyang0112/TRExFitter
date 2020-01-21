#include "TRExFitter/UnfoldingSample.h"

#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"

UnfoldingSample::UnfoldingSample() :
    fName(""),
    fTitle(""),
    fFillColor(0),
    fLineColor(0)
{
}

std::vector<Sample*> UnfoldingSample::ConvertToSample(const Region* reg,
                                                      const int bins,
                                                      const std::string& name) const {

    std::vector<Sample*> result;
    for (int ibin = 0; ibin < bins; ++ibin) {
        const std::string sampleName = "Truth_bin_" + std::to_string(ibin+1);

        Sample* sample = new Sample(sampleName, Sample::SampleType::SIGNAL);
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
