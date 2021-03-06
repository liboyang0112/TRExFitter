#include "UnfoldingCode/FoldingManager.h"

#include "UnfoldingCode/UnfoldingTools.h"

#include "TDirectory.h"
#include "TFile.h"

#include <exception>

FoldingManager::FoldingManager() : 
    fAcceptance(nullptr),
    fSelectionEfficiency(nullptr),
    fMigrationMatrix(nullptr),
    fResponseMatrix(nullptr),
    fTruthDistribution(nullptr),
    fMatrixOrientation(FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS)
{}

//__________________________________________________________________________________
//
FoldingManager::~FoldingManager() {
}

//__________________________________________________________________________________
//
void FoldingManager::SetResponseMatrix(const TH2* matrix) {
    if (!matrix) {
        throw std::runtime_error{"FoldingManager::SetResponseMatrix: Passed nullptr as the response matrix"};
    }

    fResponseMatrix.reset(static_cast<TH2D*>(matrix->Clone()));
    UnfoldingTools::Correct2DMatrix(fResponseMatrix.get());
}

//__________________________________________________________________________________
//
const TH2D* FoldingManager::GetResponseMatrix() const {
    return fResponseMatrix.get();
}

//__________________________________________________________________________________
//
void FoldingManager::SetTruthDistribution(const TH1* truth) {
    if (!truth) {
        throw std::runtime_error{"FoldingManager::SetTruthDistribution: Passed nullptr as the truth distribution"};
    }

    fTruthDistribution.reset(static_cast<TH1D*>(truth->Clone()));
}

//__________________________________________________________________________________
//
const TH1D* FoldingManager::GetTruthDistribution() const {
    return fTruthDistribution.get();
}

//__________________________________________________________________________________
//
void FoldingManager::SetSelectionEfficiency(const TH1* eff) {
    if (!eff) {
        throw std::runtime_error{"FoldingManager::SetSelectionEfficiency: Passed nullptr as the selection efficiency"};
    }
    
    fSelectionEfficiency.reset(static_cast<TH1D*>(eff->Clone()));
}

//__________________________________________________________________________________
//
const TH1D* FoldingManager::GetSelectionEfficiency() const {
    return fSelectionEfficiency.get();
}

//__________________________________________________________________________________
//
void FoldingManager::SetAcceptance(const TH1* acc) {
    if (!acc) {
        throw std::runtime_error{"FoldingManager::SetAcceptance: Passed nullptr as the acceptance"};
    }
    
    fAcceptance.reset(static_cast<TH1D*>(acc->Clone()));
}

//__________________________________________________________________________________
//
const TH1D* FoldingManager::GetAcceptance() const {
    return fAcceptance.get();
}

//__________________________________________________________________________________
//
void FoldingManager::SetMigrationMatrix(const TH2* matrix, const bool normalize) {
    if (!matrix) {
        throw std::runtime_error{"FoldingManager::SetMigrationMatrix: Passed nullptr as the migration matrix"};
    }

    fMigrationMatrix.reset(static_cast<TH2D*>(matrix->Clone()));
    UnfoldingTools::Correct2DMatrix(fMigrationMatrix.get());
    if (normalize) {
        const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
        if (horizontal) {
            UnfoldingTools::NormalizeMatrix(fMigrationMatrix.get(), false);
        } else {
            UnfoldingTools::NormalizeMatrix(fMigrationMatrix.get(), true);
        }
    }
}

//__________________________________________________________________________________
//
const TH2D* FoldingManager::GetMigrationMatrix() const {
    return fResponseMatrix.get();
}

//__________________________________________________________________________________
//
void FoldingManager::SetMatrixOrientation(const FoldingManager::MATRIXORIENTATION opt) {
    fMatrixOrientation = opt;
}

//__________________________________________________________________________________
//
FoldingManager::MATRIXORIENTATION FoldingManager::GetMatrixOrientation() const {
    return fMatrixOrientation;
}

//__________________________________________________________________________________
//
void FoldingManager::CalculateResponseMatrix(bool forceRecalculate) {
    if (!forceRecalculate && fResponseMatrix) {
        throw std::runtime_error{"FoldingManager::CalculateResponsematrix: Response matrix exists and you didnt want to recalculate it!"};
    }

    if (!CheckConsistencyForResponse()) {
        throw std::runtime_error{"FoldingManager::CalculateResponsematrix: Migration and selection efficiensy is not consistent"};
    }

    fResponseMatrix = MultiplyAcceptanceEfficiencyAndMigration(fAcceptance.get(), fSelectionEfficiency.get(), fMigrationMatrix.get());
}

//__________________________________________________________________________________
//
void FoldingManager::FoldTruth() {
    if (!CheckConsistencyForFolding()) {
        throw std::runtime_error{"FoldingManager::FoldTruth: Inconsistent truth distributions and the response matrix"};
    }

    PrepareFoldedDistributions(fTruthDistribution.get(), fResponseMatrix.get());
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1> FoldingManager::TotalFold(const TH1* const truth) const {

    if (!truth) {
        throw std::runtime_error{"FoldingManager::TotalFold: Passed nullptr"};
    }
    if (!fResponseMatrix) {
        throw std::runtime_error{"FoldingManager::TotalFold: Response matrix is nullptr"};
    }

    std::unique_ptr<TH1> result(nullptr);
    
    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    const int nRecoBins   = horizontal ? fResponseMatrix->GetNbinsY() : fResponseMatrix->GetNbinsX();
    const int nTruthBins  = horizontal ? fResponseMatrix->GetNbinsX() : fResponseMatrix->GetNbinsY();

    horizontal ? result.reset(fResponseMatrix->ProjectionY()) :
                 result.reset(fResponseMatrix->ProjectionX());
    result->Reset();

    for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
        double content(0);
        double error(0);
        for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
            content += horizontal ?
                       (truth->GetBinContent(itruth) * fResponseMatrix->GetBinContent(itruth, ireco)) :
                       (truth->GetBinContent(itruth) * fResponseMatrix->GetBinContent(ireco, itruth));

            error = horizontal ?
                    std::hypot(error, (truth->GetBinContent(itruth) * fResponseMatrix->GetBinError(itruth, ireco))):
                    std::hypot(error, (truth->GetBinContent(itruth) * fResponseMatrix->GetBinError(ireco, itruth)));
        }

        result->SetBinContent(ireco, content);
        result->SetBinError(  ireco, error);
    }

    return result;
}

//__________________________________________________________________________________
//
bool FoldingManager::CheckConsistencyForResponse() const {

    if (fResponseMatrix) return true;

    if (!fSelectionEfficiency || !fMigrationMatrix) return false;

    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    if (horizontal  && fSelectionEfficiency->GetNbinsX() != fMigrationMatrix->GetNbinsX()) return false;
    if (!horizontal && fSelectionEfficiency->GetNbinsX() != fMigrationMatrix->GetNbinsY()) return false;

    if (fAcceptance) {
        if (horizontal  && fAcceptance->GetNbinsX() != fMigrationMatrix->GetNbinsY()) return false;
        if (!horizontal && fAcceptance->GetNbinsX() != fMigrationMatrix->GetNbinsX()) return false;
    }

    return true;
}

//__________________________________________________________________________________
//
bool FoldingManager::CheckConsistencyForFolding() const {
    if (!fResponseMatrix) return false;

    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    if (horizontal  && fResponseMatrix->GetNbinsX() != fTruthDistribution->GetNbinsX()) return false;
    if (!horizontal && fResponseMatrix->GetNbinsY() != fTruthDistribution->GetNbinsX()) return false;

    return true;
}

//__________________________________________________________________________________
//
std::unique_ptr<TH2D> FoldingManager::MultiplyAcceptanceEfficiencyAndMigration(const TH1D* acc, const TH1D* sel, const TH2D* mig) const {
    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    const int nRecoBins   = horizontal ? mig->GetNbinsY() : mig->GetNbinsX();
    const int nTruthBins  = horizontal ? mig->GetNbinsX() : mig->GetNbinsY();
    
    std::unique_ptr<TH2D> result(static_cast<TH2D*>(mig->Clone()));
    for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
        const double eff = sel->GetBinContent(itruth);
        if (eff > 1) {
            throw std::runtime_error{"FoldingManager::MultiplyAcceptanceEfficiencyAndMigration: Selection efficiency is > 1"};
        }
        for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
            double acceptance(1.);
            if (fAcceptance) {
                if (acc->GetBinContent(ireco) > 1) {
                    throw std::runtime_error{"FoldingManager::MultiplyAcceptanceEfficiencyAndMigration: The acceptance in bin" + std::to_string(itruth) + " is > 1"};
                }
                acceptance = acc->GetBinContent(ireco);
            }
            const double content = horizontal ? mig->GetBinContent(itruth, ireco) : mig->GetBinContent(ireco, itruth);
            const double error = horizontal ? mig->GetBinError(itruth, ireco) : mig->GetBinError(ireco, itruth);

            if (horizontal) {
                result->SetBinContent(itruth, ireco, content*eff/acceptance);
                /// this assumes no error on truth
                result->SetBinError(itruth, ireco, error*eff/acceptance);
            } else {
                result->SetBinContent(ireco, itruth, content*eff/acceptance);
                /// this assumes no error on truth
                result->SetBinError(ireco, itruth, error*eff/acceptance);
            }
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void FoldingManager::PrepareFoldedDistributions(const TH1D* truth, const TH2D* response) {
    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    const int nRecoBins   = horizontal ? response->GetNbinsY() : response->GetNbinsX();
    const int nTruthBins  = horizontal ? response->GetNbinsX() : response->GetNbinsY();

    fFoldedDistributions.clear();

    for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
        // doing projection to get the bin edges correctly
        std::unique_ptr<TH1D> h = horizontal ? std::unique_ptr<TH1D>(response->ProjectionY()) : std::unique_ptr<TH1D>(response->ProjectionX());
        h->Reset();
        /// what to do with the uncertainties??
        fFoldedDistributions.emplace_back(*h);
        for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
            double content = horizontal ?  
                             (truth->GetBinContent(itruth) * response->GetBinContent(itruth, ireco)) :
                             (truth->GetBinContent(itruth) * response->GetBinContent(ireco, itruth));

            if (content < 0) content = 1e-6;
            const double error   = horizontal ?  
                                   (truth->GetBinContent(itruth) * response->GetBinError(itruth, ireco)) :
                                   (truth->GetBinContent(itruth) * response->GetBinError(ireco, itruth));

            fFoldedDistributions.back().SetBinContent(ireco, content);
            /// this assumes no error on truth
            fFoldedDistributions.back().SetBinError  (ireco, error);
        }
    }
}
    
//__________________________________________________________________________________
//
const std::vector<TH1D>& FoldingManager::GetFoldedDistributions() const {
    return fFoldedDistributions;
}

//__________________________________________________________________________________
//
void FoldingManager::WriteFoldedToHisto(TFile* f,
                                        const std::string& dir,
                                        const std::string& path) const {
    if (!f) {
        throw std::runtime_error{"FoldingManager::WriteFoldedToHisto: File is nullptr"};
    }

    if (fFoldedDistributions.size() == 0) {
        throw std::runtime_error{"FoldingManager::WriteFoldedToHisto: Size of the folded distributions is 0. Did you forget to run FoldingManager::FoldTruth?"};
    }

    if (dir != "") {
        const TDirectory* directory = dynamic_cast<TDirectory*>(f->Get(dir.c_str()));
        if(!directory) {
            throw std::runtime_error{"FoldingManager::WriteFoldedToHisto: Cannot cd to directory: " + dir};
        }
    }

    for (std::size_t ihist = 0; ihist < fFoldedDistributions.size(); ++ihist) {
        if (dir != "") f->cd(dir.c_str());
        else f->cd();

        const std::string fullName = path+"_bin_"+std::to_string(ihist);
        fFoldedDistributions.at(ihist).Write(fullName.c_str());
    }
}

//__________________________________________________________________________________
//
void FoldingManager::WriteTruthToHisto(TFile* f,
                                       const std::string& dir,
                                       const std::string& path) const {
    if (!f) {
        throw std::runtime_error{"FoldingManager::WriteTruthToHisto: File is nullptr"};
    }

    if (!fTruthDistribution) {
        throw std::runtime_error{"FoldingManager::WriteTruthToHisto: Truth distribution is a nullpt!"};
    }

    if (dir != "") {
        const TDirectory* directory = dynamic_cast<TDirectory*>(f->Get(dir.c_str()));
        if(!directory) {
            throw std::runtime_error{"FoldingManager::WriteTruthToHisto: Cannot cd to directory: " + dir};
        }
    }

    if (dir != "") f->cd(dir.c_str());
    else f->cd();
    fTruthDistribution->Write(path.c_str());
}

//__________________________________________________________________________________
//
void FoldingManager::WriteMigrationToHisto(TFile* f,
                                           const std::string& dir,
                                           const std::string& path) const {
    if (!f) {
        throw std::runtime_error{"FoldingManager::WriteMigrationToHisto: File is nullptr"};
    }

    if (!fMigrationMatrix) {
        throw std::runtime_error{"FoldingManager::WriteMigrationToHisto: Migration matrix is a nullpt!"};
    }

    if (dir != "") {
        const TDirectory* directory = dynamic_cast<TDirectory*>(f->Get(dir.c_str()));
        if(!directory) {
            throw std::runtime_error{"FoldingManager::WriteMigrationToHisto: Cannot cd to directory: " + dir};
        }

    }
    if (dir != "") f->cd(dir.c_str());
    else f->cd();
    fMigrationMatrix->Write(path.c_str());
}

//__________________________________________________________________________________
//
void FoldingManager::WriteSelectionEffToHisto(TFile* f,
                                              const std::string& dir,
                                              const std::string& path) const {
    if (!f) {
        throw std::runtime_error{"FoldingManager::WriteSelectionEffToHisto: File is nullptr"};
    }

    if (!fSelectionEfficiency) {
        throw std::runtime_error{"FoldingManager::WriteSelectionEffToHisto: Selection efficiency is a nullpt!"};
    }

    if (dir != "") {
        const TDirectory* directory = dynamic_cast<TDirectory*>(f->Get(dir.c_str()));
        if(!directory) {
            throw std::runtime_error{"FoldingManager::WriteSelectionEffToHisto: Cannot cd to directory: " + dir};
        }
    }

    if (dir != "") f->cd(dir.c_str());
    else f->cd();
    fSelectionEfficiency->Write(path.c_str());
}
