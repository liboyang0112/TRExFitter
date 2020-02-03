// Class include
#include "TRExFitter/PruningUtil.h"
#include "TRExFitter/Common.h"

// C++ includes
#include <memory>

// -------------------------------------------------------------------------------------------------
// class PruningUtil

//__________________________________________________________________________________
//
PruningUtil::PruningUtil() :
    fStrategy(0),
    fShapeOption(PruningUtil::SHAPEOPTION::MAXBIN),
    fThresholdNorm(-1),
    fThresholdShape(-1),
    fThresholdIsLarge(-1) {
}

//__________________________________________________________________________________
//
void PruningUtil::SetStrategy(const int strat) {
    fStrategy = strat;
}

//__________________________________________________________________________________
//
void PruningUtil::SetShapeOption(const PruningUtil::SHAPEOPTION opt) {
    fShapeOption = opt;
}

//__________________________________________________________________________________
//
void PruningUtil::SetThresholdNorm(const double thres) {
    fThresholdNorm = thres;
}

//__________________________________________________________________________________
//
void PruningUtil::SetThresholdShape(const double thres) {
    fThresholdShape = thres;
}

//__________________________________________________________________________________
//
void PruningUtil::SetThresholdIsLarge(const double thres) {
    fThresholdIsLarge = thres;
}

//__________________________________________________________________________________
//
int PruningUtil::CheckSystPruning(const TH1* const hUp,
                                  const TH1* const hDown,
                                  const TH1* const hNom,
                                  const TH1* hTot) {

    if(fStrategy!=0 && hTot==nullptr){
        std::cout << "PruningUtil::ERROR: hTot set to 0 while asking for relative pruning... Reverting to sample-by-sample pruning." << std::endl;
        fStrategy = 0;
    }
    std::unique_ptr<TH1> hRef = nullptr;
    if(fStrategy==0) hRef = std::unique_ptr<TH1>(static_cast<TH1*>(hNom->Clone()));
    else hRef = std::unique_ptr<TH1>(static_cast<TH1*>(hTot->Clone()));

    int res = 0;

    // create shape-only syst variations
    std::unique_ptr<TH1> hShapeUp        = nullptr;
    if(hUp) hShapeUp     = std::unique_ptr<TH1>(static_cast<TH1*>(hUp  ->Clone(Form("%s_shape",hUp  ->GetName()))));
    if(hShapeUp) hShapeUp->Scale( Common::EffIntegral(hNom)/Common::EffIntegral(hShapeUp.get()) );
    std::unique_ptr<TH1> hShapeDown      = nullptr;
    if(hDown) hShapeDown = std::unique_ptr<TH1>(static_cast<TH1*>(hDown->Clone(Form("%s_shape",hDown->GetName()))));
    if(hShapeDown) hShapeDown->Scale( Common::EffIntegral(hNom)/Common::EffIntegral(hShapeDown.get()) );

    // get norm effects
    const double normUp   = std::fabs((Common::EffIntegral(hUp  )-Common::EffIntegral(hNom))/Common::EffIntegral(hRef.get()));
    const double normDown = std::fabs((Common::EffIntegral(hDown)-Common::EffIntegral(hNom))/Common::EffIntegral(hRef.get()));

    // check if systematic has no shape --> 1
    bool hasShape(true);
    if(fThresholdShape>=0) {
        if (fShapeOption == PruningUtil::SHAPEOPTION::MAXBIN) {
            hasShape = HasShapeRelative(hNom,hShapeUp.get(),hShapeDown.get(),hRef.get(),fThresholdShape);
        }
        if (fShapeOption == PruningUtil::SHAPEOPTION::KSTEST) {
            hasShape = HasShapeKS(hNom,hShapeUp.get(),hShapeDown.get(),fThresholdShape);
        }
    }

    // check if systematic norm effect is under threshold
    bool hasNorm = true;
    if(fThresholdNorm>=0) hasNorm = ((normUp >= fThresholdNorm) || (normDown >= fThresholdNorm));

    // now check for crazy systematics
    bool hasGoodShape = true;
    if(fThresholdIsLarge>=0) {
        if (fShapeOption == PruningUtil::SHAPEOPTION::MAXBIN) {
            hasGoodShape = !HasShapeRelative(hNom,hShapeUp.get(),hShapeDown.get(),hRef.get(),fThresholdIsLarge);
        }
        if (fShapeOption == PruningUtil::SHAPEOPTION::KSTEST) {
            hasGoodShape = !HasShapeKS(hNom,hShapeUp.get(),hShapeDown.get(),fThresholdIsLarge);
        }
    }
    bool hasGoodNorm = true;
    if(fThresholdIsLarge>=0) hasGoodNorm = ((normUp <= fThresholdIsLarge) && (normDown <= fThresholdIsLarge));

    if(!hasGoodShape && !hasGoodNorm) res = -4;
    else if(!hasGoodShape) res = -3;
    else if(!hasGoodNorm) res = -2;
    else if(!hasShape && !hasNorm) res = 3;
    else if(!hasShape) res = 1;
    else if(!hasNorm) res = 2;

    return res;
}

//_________________________________________________________________________
//
bool PruningUtil::HasShapeRelative(const TH1* const hNom,
                                   const TH1* const hUp,
                                   const TH1* const hDown,
                                   const TH1* const combined,
                                   const double threshold) const {
    if (!hNom || !hUp || !hDown || !combined) return false;

    if (hUp->GetNbinsX() == 1) return false;

    const double& integralUp = hUp->Integral();
    const double& integralDown = hDown->Integral();
    const double& integralCombined = combined->Integral();

    if ((integralUp != integralUp) || integralUp == 0) return false;
    if ((integralDown != integralDown) || integralDown == 0) return false;
    if ((integralCombined != integralCombined) || integralCombined == 0) return false;

    bool hasShape = false;

    for (int ibin = 1; ibin <= hUp->GetNbinsX(); ++ibin){
        const double& nominal  = hNom->GetBinContent(ibin);
        if(nominal<0) continue;
        const double& comb     = combined->GetBinContent(ibin);
        const double& up       = hUp->GetBinContent(ibin);
        const double& down     = hDown->GetBinContent(ibin);
        const double& up_err   = std::fabs((up-nominal)/comb);
        const double& down_err = std::fabs((down-nominal)/comb);
        if(up_err>=threshold || down_err>=threshold){
            hasShape = true;
            break;
        }
    }

    return hasShape;
}

//_________________________________________________________________________
//
bool PruningUtil::HasShapeKS(const TH1* const hNom,
                             const TH1* const hUp,
                             const TH1* const hDown,
                             const double threshold) const {
    if (!hNom || !hUp || !hDown) return false;

    if (hUp->GetNbinsX() == 1) return false;

    bool hasShape(false);

    const double& upProb   = hUp->KolmogorovTest(hNom, "X");
    const double& downProb = hDown->KolmogorovTest(hNom, "X");
    const double probThreshold = 1 - threshold;

    if ((upProb <= probThreshold) || (downProb <= probThreshold)) {
        hasShape = true;
    }

    return hasShape;
}
