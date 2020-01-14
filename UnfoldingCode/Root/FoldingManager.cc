#include "UnfoldingCode/FoldingManager.h"

#include "TFile.h"

#include <exception>

FoldingManager::FoldingManager() : 
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
}

//__________________________________________________________________________________
//
void FoldingManager::SetResponseMatrix(TFile* file, const std::string& path) {
    if (!file) {
        throw std::runtime_error{"FoldingManager::SetResponseMatrix: Passed nullptr as the input file"};
    }

    const TH2D* tmp = dynamic_cast<TH2D*>(file->Get(path.c_str()));
    if (!tmp) {
        throw std::runtime_error{"FoldingManager::SetResponseMatrix: Cannot read response matrix from: "+path};
    }

    fResponseMatrix.reset(static_cast<TH2D*>(tmp->Clone()));
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
void FoldingManager::SetTruthDistribution(TFile* file, const std::string& path) {
    if (!file) {
        throw std::runtime_error{"FoldingManager::SetTruthDistribution: Passed nullptr as the input file"};
    }

    const TH1D* tmp = dynamic_cast<TH1D*>(file->Get(path.c_str()));
    if (!tmp) {
        throw std::runtime_error{"FoldingManager::SetTruthDistribution: Cannot read truth distribution from: "+path};
    }

    fTruthDistribution.reset(static_cast<TH1D*>(tmp->Clone()));
}

//__________________________________________________________________________________
//
const TH1D* FoldingManager::GetTruthDistribution() const {
    return fTruthDistribution.get();
}

void FoldingManager::SetMigrationMatrix(const TH2* matrix) {
    if (!matrix) {
        throw std::runtime_error{"FoldingManager::SetMigrationMatrix: Passed nullptr as the migration matrix"};
    }

    fMigrationMatrix.reset(static_cast<TH2D*>(matrix->Clone()));
}

//__________________________________________________________________________________
//
void FoldingManager::SetMigrationMatrix(TFile* file, const std::string& path) {
    if (!file) {
        throw std::runtime_error{"FoldingManager::SetMigrationMatrix: Passed nullptr as the input file"};
    }

    const TH2D* tmp = dynamic_cast<TH2D*>(file->Get(path.c_str()));
    if (!tmp) {
        throw std::runtime_error{"FoldingManager::SetMigrationMatrix: Cannot read migration matrix from: "+path};
    }

    fMigrationMatrix.reset(static_cast<TH2D*>(tmp->Clone()));
}

//__________________________________________________________________________________
//
const TH2D* FoldingManager::GetMigrationMatrix() const {
    return fResponseMatrix.get();
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
void FoldingManager::SetSelectionEfficiency(TFile* file, const std::string& path) {
    if (!file) {
        throw std::runtime_error{"FoldingManager::SetSelectionEfficiency: Passed nullptr as the input file"};
    }

    const TH1D* tmp = dynamic_cast<TH1D*>(file->Get(path.c_str()));
    if (!tmp) {
        throw std::runtime_error{"FoldingManager::SetSelectionEfficiency: Cannot read selection efficiency from: "+path};
    }

    fSelectionEfficiency.reset(static_cast<TH1D*>(tmp->Clone()));
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

    fResponseMatrix = MultiplyEfficiencyAndMigration(fSelectionEfficiency.get(), fMigrationMatrix.get());
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
bool FoldingManager::CheckConsistencyForResponse() const {

    if (fResponseMatrix) return true;

    if (!fSelectionEfficiency || !fMigrationMatrix) return false;

    if (fSelectionEfficiency->GetNbinsX() != fMigrationMatrix->GetNbinsX()) return false;

    return true;
}

//__________________________________________________________________________________
//
bool FoldingManager::CheckConsistencyForFolding() const {
    if (!fResponseMatrix) return false;

    if (fResponseMatrix->GetNbinsX() != fTruthDistribution->GetNbinsX()) return false;

    return true;
}

//__________________________________________________________________________________
//
std::unique_ptr<TH2D> FoldingManager::MultiplyEfficiencyAndMigration(const TH1D* sel, const TH2D* mig) const {
    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    const int nRecoBins   = horizontal ? mig->GetNbinsY() : mig->GetNbinsX();
    const int nTruthBins  = horizontal ? mig->GetNbinsX() : mig->GetNbinsY();
    const double recoMin  = horizontal ? mig->GetYaxis()->GetXmin() : mig->GetXaxis()->GetXmin();
    const double recoMax  = horizontal ? mig->GetYaxis()->GetXmax() : mig->GetXaxis()->GetXmax();
    const double truthMin = horizontal ? mig->GetXaxis()->GetXmin() : mig->GetYaxis()->GetXmin();
    const double truthMax = horizontal ? mig->GetXaxis()->GetXmax() : mig->GetYaxis()->GetXmax();
    
    std::unique_ptr<TH2D> result =  horizontal ? std::make_unique<TH2D>("", "", nTruthBins, truthMin, truthMax, nRecoBins, recoMin, recoMax) :
                                                 std::make_unique<TH2D>("", "", nRecoBins, recoMin, recoMax, nTruthBins, truthMin, truthMax);
    for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
        for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
            const double content = horizontal ? mig->GetBinContent(itruth, ireco) : mig->GetBinContent(ireco, itruth);
            result->SetBinContent(itruth, ireco, content*sel->GetBinContent(itruth));
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
    const double recoMin  = horizontal ? response->GetYaxis()->GetXmin() : response->GetXaxis()->GetXmin();
    const double recoMax  = horizontal ? response->GetYaxis()->GetXmax() : response->GetXaxis()->GetXmax();

    fFoldedDistributions.clear();

    for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
        fFoldedDistributions.emplace_back("","",nRecoBins, recoMin, recoMax);
        for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
            const double content = horizontal ?  
                                   (truth->GetBinContent(itruth) * response->GetBinContent(itruth, ireco)) :
                                   (truth->GetBinContent(itruth) * response->GetBinContent(ireco, itruth));
            fFoldedDistributions.back().SetBinContent(ireco, content);
        }
    }
}
    
//__________________________________________________________________________________
//
const std::vector<TH1D>& FoldingManager::GetFoldedDistributions() const {
    return fFoldedDistributions;
}
