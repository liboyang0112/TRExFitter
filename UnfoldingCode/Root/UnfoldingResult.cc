#include "UnfoldingCode/UnfoldingResult.h"

#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"

#include <exception>

UnfoldingResult::UnfoldingResult() :
    fTruthDistribution(nullptr)
{
}

UnfoldingResult::~UnfoldingResult() {
}

void UnfoldingResult::SetTruthDistribution(const TH1* truth) {
    if (!truth) {
        throw std::runtime_error{"UnfoldingResult::SetTruthDistribution: Passing nullptr"};
    }

    fTruthDistribution.reset(static_cast<TH1D*>(truth->Clone()));
}

std::unique_ptr<TGraphAsymmErrors> UnfoldingResult::GetUnfoldedResult() const {

    if (!fTruthDistribution) {
        throw std::runtime_error{"UnfoldingResult::GetUnfoldedResult: No truth distribution is set"};
    }

    if (static_cast<int>(fFitValues.size()) != fTruthDistribution->GetNbinsX()) {
        throw std::runtime_error{"UnfoldingResult::GetUnfoldedResult: Size of the passed FitValues doesnt math the number of bins of the truth distribution"};
    }

    // Use cosntructor from TH1 to get the correct binninb and horizontal errors
    auto result = std::make_unique<TGraphAsymmErrors>(fTruthDistribution.get());

    for (int ibin = 1; ibin < fTruthDistribution->GetNbinsX(); ++ibin) {
        const double error_x_low = result->GetErrorXlow(ibin);
        const double error_x_high = result->GetErrorXhigh(ibin);
        double x;
        double y;
        result->GetPoint(ibin, x, y);

        const double mean = fFitValues.at(ibin-1).nominal * fTruthDistribution->GetBinContent(ibin);
        const double up   = fFitValues.at(ibin-1).up      * fTruthDistribution->GetBinContent(ibin);
        const double down = fFitValues.at(ibin-1).down    * fTruthDistribution->GetBinContent(ibin);

        result->SetPoint(ibin, x, mean);
        result->SetPointError(ibin, error_x_low, error_x_high, down, up);
    }

    return result;
}
