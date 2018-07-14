#ifndef HISTOTOOLS_H
#define HISTOTOOLS_H

#include <vector>
#include <string>

//foward declarations
class TH1;
class TH1F;
class SystematicHist;

namespace HistoTools {

    enum HistogramOperations {
        //Symmetrisation operations are units
        SYMMETRIZEONESIDED = 1, //symmetrize one-sided systematics (e.g. JER)
        SYMMETRIZETWOSIDED = 2, // symmetrize two-sided systematics (protects from statistical fluctuations)

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

    TH1F* TranformHistogramBinning(TH1* originalHist);

    void ManageHistograms(int histOps,  TH1* hNom, TH1* originUp, TH1* originDown, TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown, const SmoothOption &smoothOpt, bool TtresSmoothing = false);
    void SymmetrizeHistograms(int histOps,  TH1* hNom, TH1* originUp, TH1* originDown, TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown);
    void SmoothHistograms(int histOps,  TH1* hNom, TH1* originUp, TH1* originDown, TH1* &modifiedUp, TH1* &modifiedDown, const SmoothOption &smoothOpt, bool TtresSmoothing = false);

    //Symmetrisation functions
    TH1F* SymmetrizeOneSided( TH1* h_nominal, TH1* h_syst, bool &isUp );
    TH1F* InvertShift(TH1* h_syst, TH1* h_nominal);
    float Separation(TH1* h1,TH1* h2);
    TH1F* SymmetrizeTwoSided(TH1* var1, TH1* var2, TH1* hnom);

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
