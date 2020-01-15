#include "UnfoldingCode/UnfoldingTools.h"

#include "TH2.h"

#include <exception>

void UnfoldingTools::NormalizeMatrix(TH2* matrix, const bool byRow) {
    if (!matrix) {
        throw std::runtime_error{"UnfoldingTools::NormalizeMatrix: Nullptr passed!"};
    }

    if (byRow) {
        for (int ibiny = 1; ibiny <= matrix->GetNbinsY(); ++ibiny) {
            double sum(0.);
            for (int ibinx = 1; ibinx <= matrix->GetNbinsX(); ++ibinx) {
                sum += matrix->GetBinContent(ibinx, ibiny);
            }
            if (sum < 1e-9) continue;
            for (int ibinx = 1; ibinx <= matrix->GetNbinsX(); ++ibinx) {
                const double content = matrix->GetBinContent(ibinx, ibiny);
                matrix->SetBinContent(ibinx, ibiny, content/sum);
            }
        }
    } else {
        for (int ibinx = 1; ibinx <= matrix->GetNbinsX(); ++ibinx) {
            double sum(0.);
            for (int ibiny = 1; ibiny <= matrix->GetNbinsY(); ++ibiny) {
                sum += matrix->GetBinContent(ibinx, ibiny);
            }
            if (sum < 1e-9) continue;
            for (int ibiny = 1; ibiny <= matrix->GetNbinsY(); ++ibiny) {
                const double content = matrix->GetBinContent(ibinx, ibiny);
                matrix->SetBinContent(ibinx, ibiny, content/sum);
            }
        }
    }
}
