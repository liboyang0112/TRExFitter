#include "TRExFitter/TruthSample.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/TRExFit.h"

#include "TH1.h"

TruthSample::TruthSample(const std::string& name) :
    fName(name),
    fTitle(""),
    fLineStyle(0),
    fLineColor(1),
    fTruthDistributionPath(""),
    fTruthDistributionFile(""),
    fTruthDistributionName("")
{
}

std::unique_ptr<TH1> TruthSample::GetHisto(const TRExFit* fitter) const {

    const std::string path = fTruthDistributionPath == "" ? fitter->fTruthDistributionPath : fTruthDistributionPath;
    const std::string file = fTruthDistributionFile == "" ? fitter->fTruthDistributionFile : fTruthDistributionFile;
    const std::string name = fTruthDistributionName == "" ? fitter->fTruthDistributionName : fTruthDistributionName;

    const std::string full = path+"/"+file+".root/"+name;
    
    return Common::HistFromFile(full);
}
