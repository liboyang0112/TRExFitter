#ifndef __HistoTools__
#define __HistoTools__

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
        
        //Other (possible) functionnalities
        UNKNOWN = 100
    };
    
    TH1F* TranformHistogramBinning(TH1* originalHist);
    
    void ManageHistograms(int histOps,  TH1* hNom, TH1* originUp, TH1* originDown, TH1* &modifiedUp, TH1* &modifiedDown);
    
    //Symmetrisation functions
    TH1F* SymmetrizeOneSided( TH1* h_nominal, TH1* h_syst, bool &isUp );
    TH1F* InvertShift(TH1* h_syst, TH1* h_nominal);
    float Separation(TH1* h1,TH1* h2);
    TH1F* SymmetrizeTwoSided(TH1* var1, TH1* var2, TH1* hnom);
    
    //Smoothing uilities
    int rebin_getMaxVar(TH1* hnom,TH1* hsyst, double tolerance);
    int getBinWidth(TH1 *ratio);
    void Smooth_maxVariations(TH1* hsyst,TH1* hnom, int nbins);
    int get_nVar(TH1* hratio);
    
    //Has systematic
    bool HasShape(TH1* nom, SystematicHist* sh, float threshold);
}
#endif
