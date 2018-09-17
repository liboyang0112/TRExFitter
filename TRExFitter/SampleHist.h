#ifndef SAMPLEHIST_H
#define SAMPLEHIST_H

/// Framework includes
#include "TRExFitter/HistoTools.h"

/// ROOT includes
#include "Rtypes.h"

/// c++ includes
#include <map>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class TFile;
class TH1;
class Sample;
class NormFactor;
class ShapeFactor;
class SystematicHist;

class SampleHist {
public:
    SampleHist();
    SampleHist(Sample *sample,TH1 *hist);
    SampleHist(Sample *sample, const std::string& histoName, const std::string& fileName);
    ~SampleHist();

    TH1* GetHist() const;
    Sample* GetSample() const;
    SystematicHist* AddOverallSyst(const std::string& name,float up,float down);
    SystematicHist* AddStatSyst(const std::string& name,int i_bin);
    SystematicHist* AddHistoSyst(const std::string& name,TH1* h_up,TH1* h_down);
    SystematicHist* AddHistoSyst(const std::string& name, const std::string& histoName_up,
                                 const std::string& fileName_up, const std::string& histoName_down,
                                 const std::string& fileName_down, int pruned=0);
    SystematicHist* GetSystematic(const std::string& systName) const;
    SystematicHist* GetSystFromNP(const std::string& NuisParName) const;
    NormFactor* AddNormFactor(const std::string& name,float nominal, float min, float max);
    NormFactor* AddNormFactor(NormFactor *normFactor);
    NormFactor* GetNormFactor(const std::string& name) const;
    ShapeFactor* AddShapeFactor(const std::string& name,float nominal, float min, float max);
    ShapeFactor* AddShapeFactor(ShapeFactor *shapeFactor);
    ShapeFactor* GetShapeFactor(const std::string& name) const;

    bool HasSyst(const std::string& name) const;
    bool HasNorm(const std::string& name) const;
    bool HasShapeFactor(const std::string& name) const;

    void WriteToFile(TFile *f=0x0,bool reWriteOrig=true);
    void ReadFromFile();

    void FixEmptyBins(const bool suppress);

    void Print() const;

    void Rebin(int ngroup = 2, const Double_t* xbins = 0);
    void DrawSystPlot( const std::string &syst="all", TH1* h_data=0x0,
                       bool SumAndData=false, bool bothPanels=false ) const;
    void SmoothSyst(const HistoTools::SmoothOption &opt, std::string syst="all",bool force=false, bool TtresSmoothing = false);

    void Divide(  SampleHist* sh);
    void Multiply(SampleHist* sh);
    void Add(     SampleHist* sh,float scale=1.);
    void Scale(float scale);

    void SampleHistAdd(SampleHist* h, float scale = 1.);
    void CloneSampleHist(SampleHist* h, const std::set<std::string>& names, float scale = 1.);

    std::string fName;
    Sample *fSample;
    TH1 *fHist;
    TH1 *fHist_orig;
    TH1 *fHist_regBin;
    TH1 *fHist_preSmooth; // new - to use only for syst plots
    TH1 *fHist_postFit;
    std::string fFileName;
    std::string fHistoName;
    bool fIsData;
    bool fIsSig;
    std::map<std::string,bool> fIsMorph;

    int fNSyst;
    std::vector < SystematicHist* > fSyst;
    int fNNorm;
    std::vector < NormFactor* > fNormFactors;
    int fNShape;
    std::vector < ShapeFactor* > fShapeFactors;

    // other useful info
    std::string fFitName;
    std::string fRegionName;
    std::string fRegionLabel;
    std::string fVariableTitle;
    bool fSystSmoothed;
};

#endif

