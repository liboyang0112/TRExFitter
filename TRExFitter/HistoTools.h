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

    enum HistogramOperations {
        //Symmetrisation operations are units
        SYMMETRIZEONESIDED = 1, //symmetrize one-sided systematics (e.g. JER)
        SYMMETRIZETWOSIDED = 2, // symmetrize two-sided systematics (protects from statistical fluctuations)
        SYMMETRIZEABSMEAN = 3, // symmetrize two-sided systematics by taking mean ob the aboslute values from up and down shift (protects from statistical fluctuations)
        SYMMETRIZEMAXIMUM = 4, // symmetrize two-sided systematics by taking the larger variation from up and down and symmetrizing (protects from statistical fluctuations)

        //Smoothing operations are 10th
        SMOOTH = 10,

        SMOOTH_INDEPENDENT = 100,

        //Other (possible) functionnalities
        UNKNOWN = 1000
    };

    enum SmoothOption {
        MAXVARIATION = 0,
        TTBARRESONANCE = 1,
        COMMONTOOLSMOOTHMONOTONIC = 2,
        COMMONTOOLSMOOTHPARABOLIC = 3,
        KERNELRATIOUNIFORM = 4,
        KERNELDELTAGAUSS = 5,
        KERNELRATIOGAUSS = 6
    };

    TH1D* TranformHistogramBinning(TH1* originalHist);

    void ManageHistograms(int histOps,  TH1* hNom, TH1* originUp, TH1* originDown, TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown, const SmoothOption &smoothOpt, bool TtresSmoothing = false);
    void SymmetrizeHistograms(int histOps,  TH1* hNom, TH1* originUp, TH1* originDown, TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown);
    void SmoothHistograms(int histOps,  TH1* hNom, TH1* &modifiedUp, TH1* &modifiedDown, const SmoothOption &smoothOpt, bool TtresSmoothing = false);

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
     * A helper function to prnt warning if the two variations have the same shift and TwoSided
     * symmetrization is used.
     * @param Original variation up
     * @param Original variation down
     * @param Nominal histogram
     * @param Modified syst variation up
     * @param Modifiec syst variation down
     */
    void CheckSameShift(const TH1* const var1, const TH1* const var2, const TH1* const hnom,
    const TH1D* const tmp1, const TH1D* const tmp2);

    void Scale(TH1* h_syst, TH1* h_nominal, float factor);


    struct Bin {
        double N;
        double S;
        double dN2;
        double dS2;
        double edge;
        Bin(double _N, double _S, double _dN2, double _dS2, double _edge) { N = _N; S = _S; dN2 = _dN2; dS2 = _dS2; edge = _edge; }
    };
    double avgError(std::vector<Bin> &hist, bool independentVar);
    bool systSmallerThanStat(std::vector<Bin> &hist, bool independentVar, double avgError);

    //Has systematic
    bool HasShape(TH1* nom, SystematicHist* sh, float threshold);

    //Histograms checker
    bool CheckHistograms(TH1* nom, SystematicHist* sh, bool checkNull = true, bool causeCrash = false);
}
#endif
