#include "UnfoldingCode/FoldingManager.h"

#include "TFile.h"

#include <exception>

FoldingManager::FoldingManager() : 
    fSelectionEfficiency(nullptr),
    fMigrationMatrix(nullptr),
    fResponseMatrix(nullptr),
    fTruthDistribution(nullptr)
{}

FoldingManager::~FoldingManager() {
}

void FoldingManager::SetResponseMatrix(const TH2* matrix) {
    if (!matrix) {
        throw std::runtime_error{"FoldingManager::SetResponseMatrix: Passed nullptr as the response matrix"};
    }

    fResponseMatrix.reset(static_cast<TH2D*>(matrix->Clone()));
}

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

const TH2D* FoldingManager::GetResponseMatrix() const {
    return fResponseMatrix.get();
}

void FoldingManager::SetTruthDistribution(const TH1* truth) {
    if (!truth) {
        throw std::runtime_error{"FoldingManager::SetTruthDistribution: Passed nullptr as the truth distribution"};
    }

    fTruthDistribution.reset(static_cast<TH1D*>(truth->Clone()));
}

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

const TH1D* FoldingManager::GetTruthDistribution() const {
    return fTruthDistribution.get();
}

void FoldingManager::SetMigrationMatrix(const TH2* matrix) {
    if (!matrix) {
        throw std::runtime_error{"FoldingManager::SetMigrationMatrix: Passed nullptr as the migration matrix"};
    }

    fMigrationMatrix.reset(static_cast<TH2D*>(matrix->Clone()));
}

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

const TH2D* FoldingManager::GetMigrationMatrix() const {
    return fResponseMatrix.get();
}

void FoldingManager::SetSelectionEfficiency(const TH1* eff) {
    if (!eff) {
        throw std::runtime_error{"FoldingManager::SetSelectionEfficiency: Passed nullptr as the selection efficiency"};
    }

    fSelectionEfficiency.reset(static_cast<TH1D*>(eff->Clone()));
}

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

void FoldingManager::CalculateResponseMatrix(bool forceRecalculate) {
    if (!forceRecalculate && fResponseMatrix) {
        throw std::runtime_error{"FoldingManager::CalculateResponsematrix: Response matrix exists and you didnt want to recalculate it!"};
    }

    if (!CheckConsistencyForResponse()) {
        throw std::runtime_error{"FoldingManager::CalculateResponsematrix: Migration and selection efficiensy is not consistent"};
    }

    fResponseMatrix = MultiplyEfficiencyAndMigration(fSelectionEfficiency.get(), fMigrationMatrix.get());
}

void FoldingManager::FoldTruth() {
    if (!CheckConsistencyForFolding()) {
        throw std::runtime_error{"FoldingManager::FoldTruth: Inconsistent truth distributions and the response matrix"};
    }

    GetFoldedDistributions(fTruthDistribution.get(), fResponseMatrix.get());
}

bool FoldingManager::CheckConsistencyForResponse() const {

    if (fResponseMatrix) return true;

    if (!fSelectionEfficiency || !fMigrationMatrix) return false;

    if (fSelectionEfficiency->GetNbinsX() != fMigrationMatrix->GetNbinsX()) return false;

    return true;
}

bool FoldingManager::CheckConsistencyForFolding() const {
    if (!fResponseMatrix) return false;

    if (fResponseMatrix->GetNbinsX() != fTruthDistribution->GetNbinsX()) return false;

    return true;
}

std::unique_ptr<TH2D> FoldingManager::MultiplyEfficiencyAndMigration(const TH1D* sel, const TH2D* mig) const {
    const int nRecoBins   = mig->GetNbinsY();
    const int nTruthBins  = mig->GetNbinsX();
    const double recoMin  = mig->GetYaxis()->GetXmin();
    const double recoMax  = mig->GetYaxis()->GetXmax();
    const double truthMin = mig->GetXaxis()->GetXmin();
    const double truthMax = mig->GetXaxis()->GetXmax();
    
    std::unique_ptr<TH2D> result = std::make_unique<TH2D>("", "", nTruthBins, truthMin, truthMax, nRecoBins, recoMin, recoMax);
    for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
        for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
            const double content = mig->GetBinContent(itruth, ireco);
            result->SetBinContent(itruth, ireco, content*sel->GetBinContent(itruth));
        }
    }

    return result;
}

void FoldingManager::GetFoldedDistributions(const TH1D* truth, const TH2D* response) {
    const int nRecoBins   = response->GetNbinsY();
    const int nTruthBins  = response->GetNbinsX();
    const double recoMin  = response->GetYaxis()->GetXmin();
    const double recoMax  = response->GetYaxis()->GetXmax();

    fFoldedDistributions.clear();

    for (int itruth = 1; itruth <= nTruthBins; ++itruth) {
        fFoldedDistributions.emplace_back("","",nRecoBins, recoMin, recoMax);
        for (int ireco = 1; ireco <= nRecoBins; ++ireco) {
            const double content = truth->GetBinContent(itruth) * response->GetBinContent(itruth, ireco);
            fFoldedDistributions.back().SetBinContent(ireco, content);
        }
    }
}
    
const std::vector<TH1D>& FoldingManager::GetFoldedDistributions() const {
    return fFoldedDistributions;
}
