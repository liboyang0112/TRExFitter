#ifndef SAMPLEHIST_H
#define SAMPLEHIST_H

/// Framework includes
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/PruningUtil.h"

/// ROOT includes
#include "Rtypes.h"

/// c++ includes
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class TFile;
class TH1;
class TPad;
class Sample;
class NormFactor;
class ShapeFactor;
class SystematicHist;

class SampleHist {
public:
    explicit SampleHist();
    explicit SampleHist(Sample *sample,TH1 *hist);
    explicit SampleHist(Sample *sample, const std::string& histoName, const std::string& fileName);
    ~SampleHist();
    SampleHist(const SampleHist& s) = delete;
    SampleHist(SampleHist&& s) = delete;
    SampleHist& operator=(const SampleHist& s) = delete;
    SampleHist& operator=(SampleHist&& s) = delete;

    TH1* GetHist() const;
    const Sample* GetSample() const {return fSample;}
    SystematicHist* AddOverallSyst(const std::string& name,const std::string& storedName,double up,double down);
    SystematicHist* AddStatSyst(const std::string& name,const std::string& storedName,int i_bin);
    SystematicHist* AddHistoSyst(const std::string& name,const std::string& storedName,TH1* h_up,TH1* h_down);
    SystematicHist* AddHistoSyst(const std::string& name,const std::string& storedName, const std::string& histoName_up,
                                 const std::string& fileName_up, const std::string& histoName_down,
                                 const std::string& fileName_down, int pruned=0);
    SystematicHist* GetSystematic(const std::string& systName) const;
    SystematicHist* GetSystFromNP(const std::string& NuisParName) const;
    NormFactor* AddNormFactor(const std::string& name,double nominal, double min, double max);
    NormFactor* AddNormFactor(NormFactor *normFactor);
    NormFactor* GetNormFactor(const std::string& name) const;
    ShapeFactor* AddShapeFactor(const std::string& name,double nominal, double min, double max);
    ShapeFactor* AddShapeFactor(ShapeFactor *shapeFactor);
    ShapeFactor* GetShapeFactor(const std::string& name) const;

    bool HasSyst(const std::string& name) const;
    bool HasNorm(const std::string& name) const;
    bool HasShapeFactor(const std::string& name) const;

    void WriteToFile(TFile *f=0x0,bool reWriteOrig=true);
    void ReadFromFile();

    void FixEmptyBins(const bool suppress);
    void NegativeTotalYieldWarning(TH1* hist, double yield) const;

    void Print() const;

    void Rebin(int ngroup = 2, const Double_t* xbins = 0);
    void DrawSystPlot( const std::string &syst="all", TH1* h_data=0x0,
                       bool SumAndData=false, bool bothPanels=false ) const;
    void SmoothSyst(const HistoTools::SmoothOption &opt, std::string syst="all", bool force=false);

    void Divide(  SampleHist* sh);
    void Multiply(SampleHist* sh);
    void Add(     SampleHist* sh,double scale=1.);
    void Scale(double scale);

    void SampleHistAdd(SampleHist* h, double scale = 1.);
    void SampleHistAddNominal(SampleHist* h, double scale);
    void CloneSampleHist(SampleHist* h, const std::set<std::string>& names, double scale = 1.);
    void SystPruning(PruningUtil *pu,TH1* hTot=nullptr);

    void DrawSystPlotUpper(TPad* pad0,
                           TH1* nominal,
                           TH1* nominal_orig,
                           TH1* syst_up,
                           TH1* syst_up_orig,
                           TH1* syst_down,
                           TH1* syst_down_orig,
                           TH1* data,
                           TH1* tmp,
                           bool SumAndData,
                           bool bothPanels) const;

    void DrawSystPlotRatio(TPad* pad1,
                           TH1* nominal,
                           TH1* nominal_orig,
                           TH1* syst_up,
                           TH1* syst_up_orig,
                           TH1* syst_down,
                           TH1* syst_down_orig,
                           TH1* data,
                           TH1* tmp,
                           bool SumAndData) const;

    std::string fName;
    Sample *fSample;
    std::unique_ptr<TH1> fHist;
    std::unique_ptr<TH1> fHist_orig;
    std::unique_ptr<TH1> fHist_regBin;
    std::unique_ptr<TH1> fHist_preSmooth; // new - to use only for syst plots
    std::unique_ptr<TH1> fHist_postFit;
    std::string fFileName;
    std::string fHistoName;
    bool fIsData;
    bool fIsSig;
    std::map<std::string,bool> fIsMorph;

    int fNSyst;
    std::vector < std::unique_ptr<SystematicHist> > fSyst;
    int fNNorm;
    std::vector < std::unique_ptr<NormFactor> > fNormFactors;
    int fNShape;
    std::vector < std::unique_ptr<ShapeFactor> > fShapeFactors;

    // other useful info
    std::string fFitName;
    std::string fRegionName;
    std::string fRegionLabel;
    std::string fVariableTitle;
    bool fSystSmoothed;
};

#endif

