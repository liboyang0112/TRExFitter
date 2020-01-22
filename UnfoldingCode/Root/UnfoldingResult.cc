#include "UnfoldingCode/UnfoldingResult.h"

#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"

#include <exception>
#include <fstream>
#include <iomanip>

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
    fTruthDistribution->SetDirectory(nullptr);
}

std::unique_ptr<TGraphAsymmErrors> UnfoldingResult::GetUnfoldedResultErrorBand() const {

    if (!fTruthDistribution) {
        throw std::runtime_error{"UnfoldingResult::GetUnfoldedResultErrorband: No truth distribution is set"};
    }

    if (static_cast<int>(fFitValues.size()) != fTruthDistribution->GetNbinsX()) {
        throw std::runtime_error{"UnfoldingResult::GetUnfoldedResultErrorBand: Size of the passed FitValues doesnt math the number of bins of the truth distribution"};
    }

    // Use cosntructor from TH1 to get the correct binninb and horizontal errors
    auto result = std::make_unique<TGraphAsymmErrors>(fTruthDistribution.get());

    for (int ibin = 1; ibin <= fTruthDistribution->GetNbinsX(); ++ibin) {
        const double error_x_low = result->GetErrorXlow(ibin-1);
        const double error_x_high = result->GetErrorXhigh(ibin-1);
        double x;
        double y;
        result->GetPoint(ibin-1, x, y);

        const double mean = fFitValues.at(ibin-1).nominal * fTruthDistribution->GetBinContent(ibin);
        const double up   = std::fabs((fFitValues.at(ibin-1).up      * fTruthDistribution->GetBinContent(ibin)) - mean);
        const double down = std::fabs((fFitValues.at(ibin-1).down    * fTruthDistribution->GetBinContent(ibin)) - mean);

        result->SetPoint(ibin-1, x, mean);
        result->SetPointError(ibin-1, error_x_low, error_x_high, down, up);
    }

    return result;
}
    
void UnfoldingResult::AddFitValue(const double nominal, const double up, const double down) {
    FitValue value;
    value.nominal = nominal;
    value.up = up;
    value.down = down;

    fFitValues.emplace_back(value);
}

std::unique_ptr<TH1D> UnfoldingResult::GetUnfoldedResult() const {
    if (!fTruthDistribution) {
        throw std::runtime_error{"UnfoldingResult::GetUnfoldedResult: No truth distribution is set"};
    }

    if (static_cast<int>(fFitValues.size()) != fTruthDistribution->GetNbinsX()) {
        throw std::runtime_error{"UnfoldingResult::GetUnfoldedResult: Size of the passed FitValues doesnt math the number of bins of the truth distribution"};
    }

    std::unique_ptr<TH1D> result(static_cast<TH1D*>(fTruthDistribution->Clone()));
    result->SetDirectory(nullptr);

    for (int ibin = 1; ibin <= fTruthDistribution->GetNbinsX(); ++ibin) {
        const double mean = fFitValues.at(ibin-1).nominal * fTruthDistribution->GetBinContent(ibin);
        result->SetBinContent(ibin, mean);
        result->SetBinError(ibin, 0);
    }
    
    return result;
}

void UnfoldingResult::DumpResults(std::ofstream* stream) const {
    if (!stream) {
        throw std::runtime_error{"UnfoldingResult::DumpResults : Nullptr passed"};
    }

    if (!stream->is_open() || !stream->good()) {
        throw std::runtime_error{"UnfoldingResult::DumpResults : Problematic stream"};
    }

    if (static_cast<int>(fFitValues.size()) != fTruthDistribution->GetNbinsX()) {
        throw std::runtime_error{"UnfoldingResult::DumpResults: Size of the passed FitValues doesnt math the number of bins of the truth distribution"};
    }

    *stream << std::fixed << std::setprecision(4);
    *stream << "Yields in bins for unfolded data\n";
    *stream << "Bin number nominal Up uncertainty Down uncertainty\n";
    for (int ibin = 1; ibin <= fTruthDistribution->GetNbinsX(); ++ibin) {
        const double mean = fFitValues.at(ibin-1).nominal * fTruthDistribution->GetBinContent(ibin);
        const double up   = fFitValues.at(ibin-1).up * fTruthDistribution->GetBinContent(ibin) - mean;
        const double down = fFitValues.at(ibin-1).down * fTruthDistribution->GetBinContent(ibin) - mean;
        *stream << ibin << " " << mean << " " << up << " " << down << "\n";
    }
}
