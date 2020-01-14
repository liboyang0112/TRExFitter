#ifndef FOLDINGMANAGER_H_
#define FOLDINGMANAGER_H_

#include "TH1.h"
#include "TH2D.h"

#include <memory>
#include <string>
#include <vector>

class TFile;

class FoldingManager {
public:
    explicit FoldingManager();

    ~FoldingManager();

    FoldingManager(const FoldingManager& rhs) = delete;
    FoldingManager(FoldingManager&& rhs) = delete;
    FoldingManager& operator=(const FoldingManager& rhs) = delete;
    FoldingManager& operator=(FoldingManager&& rhs) = delete;

    void SetResponseMatrix(const TH2* matrix);
    
    void SetResponseMatrix(TFile* file, const std::string& path);

    const TH2D* GetResponseMatrix() const;

    void SetTruthDistribution(const TH1* truth);
    
    void SetTruthDistribution(TFile* file, const std::string& path);

    const TH1D* GetTruthDistribution() const;

    void SetSelectionEfficiency(const TH1* eff);
    
    void SetSelectionEfficiency(TFile* file, const std::string& path);

    const TH1D* GetSelectionEfficiency() const;

    void SetMigrationMatrix(const TH2* matrix);
    
    void SetMigrationMatrix(TFile* file, const std::string& path);

    const TH2D* GetMigrationMatrix() const;

    void CalculateResponseMatrix(bool forceRecalculate);

    void FoldTruth();

    const std::vector<TH1D>& GetFoldedDistributions() const;

private:
    std::unique_ptr<TH1D> fSelectionEfficiency;
    std::unique_ptr<TH2D> fMigrationMatrix;
    std::unique_ptr<TH2D> fResponseMatrix;
    std::unique_ptr<TH1D> fTruthDistribution;
    std::vector<TH1D> fFoldedDistributions;

    bool CheckConsistencyForResponse() const;
    
    bool CheckConsistencyForFolding() const;

    std::unique_ptr<TH2D> MultiplyEfficiencyAndMigration(const TH1D* sel, const TH2D* mig) const;

    void PrepareFoldedDistributions(const TH1D* truth, const TH2D* response);
};

#endif
