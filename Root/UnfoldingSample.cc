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

Sample* UnfoldingSample::ConvertToSample(const Region* reg) const {
    return nullptr;
}
