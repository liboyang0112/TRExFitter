#include "TRExFitter/TruthSample.h"

TruthSample::TruthSample(const std::string& name) :
    fName(name),
    fTitle(""),
    fFillColor(1),
    fLineColor(1),
    fTruthDistributionPath(""),
    fTruthDistributionFile(""),
    fTruthDistributionName("")
{
}
