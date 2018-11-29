#ifndef HISTOTOOLS_H
#define HISTOTOOLS_H

/// c++ includes
#include <vector>
#include <string>

/// Foward class declaration
class TH1;
class TH1D;
class SystematicHist;

namespace HistoTools {

    /// Enum that stores different symmetrization options
    enum SymmetrizationType {
        NOSYMMETRIZATION = 0, /// no symmetrization applied
        SYMMETRIZEONESIDED = 1, ///symmetrize one-sided systematics (e.g. JER)
        SYMMETRIZETWOSIDED = 2, /// symmetrize two-sided systematics (protects from statistical fluctuations)
        SYMMETRIZEABSMEAN = 3, /// symmetrize two-sided systematics by taking mean of the absolute values from up and down shift (protects from statistical fluctuations)
        SYMMETRIZEMAXIMUM = 4, /// symmetrize two-sided systematics by taking the larger variation from up and down and symmetrizing (protects from statistical fluctuations)
    };

    /// Enum that stores different smoothing level options
    enum SmoothLevel {
        SMOOTHDEPENDENT = 10,
        SMOOTHINDEPENDENT = 100,
        UNKNOWN = 1000
    };

    /// Enum that stores different smoothing options
    enum SmoothOption {
        MAXVARIATION = 0,
        TTBARRESONANCE = 1,
        COMMONTOOLSMOOTHMONOTONIC = 2,
        COMMONTOOLSMOOTHPARABOLIC = 3,
        KERNELRATIOUNIFORM = 4,
        KERNELDELTAGAUSS = 5,
        KERNELRATIOGAUSS = 6
    };

    /**
    * In RooStats, input histogram variable binning is not supported => convert to a constant binning
    * by creating an histogram with the same number of bins but with constant binning between 0 and 1
    * - now in case some bins are < 0 (due to bin drop functionality), they are ignored for the regBin histos
    * @param Original histogram
    * @return Trasformed histogram
    */
    TH1D* TranformHistogramBinning(TH1* originalHist);

    /**
     * A helper function to smooth and symmetrize histograms
     * @param level of smoothing
     * @param type of symmetrization
     * @param nominal histogram
     * @param up variation histogram
     * @param down variation histogram
     * @param symmetrized/smoothed up variation
     * @param symmetrized/smoothed down variation
     * @param scale for the up variation
     * @param scale for the down varaition
     * @param smoothing option
     * @param apply ttbar resonance smoothing
     */
    void ManageHistograms(int smoothingLevel, const SymmetrizationType& symType, TH1* hNom, TH1* originUp, TH1* originDown,
        TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown, const SmoothOption &smoothOpt, bool TtresSmoothing = false);

    /**
     * A helper function to symmetrize histograms
     * @param symmetrization type
     * @param nominal histogram
     * @param up variation histogram
     * @param down variation histogram
     * @param symmetrized/smoothed up variation
     * @param symmetrized/smoothed down variation
     * @param scale for the up variation
     * @param scale for the down varaition
     */
    void SymmetrizeHistograms(const SymmetrizationType& symType,  TH1* hNom, TH1* originUp, TH1* originDown,
        TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown);

    /**
     * A helper function to smooth histograms
     * @param level of smoothing
     * @param symmetrized/smoothed up variation
     * @param symmetrized/smoothed down variation
     * @param smooth option
     * @param apply ttbar resonance smoothing
     */
    void SmoothHistograms(int smoothingLevels,  TH1* hNom, TH1* &modifiedUp, TH1* &modifiedDown,
        const SmoothOption &smoothOpt, bool TtresSmoothing = false);

    /**
     * A helper function to Symmetrize systematics using one sided method
     * @param Nominal histogram
     * @param Systematic histogram
     * @param Flag for up/down variation
     * @return Modified histogram
     */
    TH1D* SymmetrizeOneSided(const TH1* const h_nominal, const TH1* const h_syst, bool &isUp );

    /**
     * A helper function to invert one sided shift
     * @param Systematic histogram
     * @param Nominal histogram
     * @return Modified histogram
     */
    TH1D* InvertShift(const TH1* const h_syst, const TH1* const h_nominal);

    /**
     * A helper function to calculate separation between histograms
     * @param First histogram
     * @param Second histogram
     * @return Separation value
     */
    float Separation(const TH1* const h1, const TH1* const h2);

    /**
     * A helper function to calculate symmetrization using two sided variation
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    TH1D* SymmetrizeTwoSided(const TH1* const var1, const TH1* const var2, const TH1* const hnom);

    /**
     * A helper function to calculate symmetrization using mean of the absolute values
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    TH1D* SymmetrizeAbsMean(const TH1* const var1, const TH1* const var2, const TH1* const hnom);

    /**
     * A helper function to calculate symmetrization using maximum variation
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    TH1D* SymmetrizeMaximum(const TH1* const var1, const TH1* const var2, const TH1* const hnom);

    /**
     * A helper function to print warning if the two variations have the same shift and TwoSided
     * symmetrization is used.
     * @param Original variation up
     * @param Original variation down
     * @param Nominal histogram
     * @param Modified syst variation up
     * @param Modifiec syst variation down
     */
    void CheckSameShift(const TH1* const var1, const TH1* const var2, const TH1* const hnom,
    const TH1D* const tmp1, const TH1D* const tmp2);

    /**
     * A helper function to scale systeamtic histograms
     * @param Systematic histogram that will be scaled
     * @param Npminal histogram
     * @param scale factor
     */
    void Scale(TH1* h_syst, TH1* h_nominal, float factor);

    /**
     * A helper function to check if the syst variation has shape effect
     * @param Nominal histogram
     * @param Systematic histogram
     * @param Threshold used for the check
     * @return Has shape
     */
    bool HasShape(TH1* nom, SystematicHist* sh, float threshold);

    /**
     * A helper function to check for various weird features
     * @param Nominal histogram
     * @param Systematic histogram
     * @param Check for nullpointer
     * @param Crash when histo is nullptr?
     * @return Histo is OK
     */
    bool CheckHistograms(TH1* nom, SystematicHist* sh, bool checkNull = true, bool causeCrash = false);

    bool HasShapeRelative(const TH1* const hNom, const TH1* const up, const TH1* const down, const TH1* const combined, float threshold);
}
#endif
