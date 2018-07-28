#ifndef SAMPLEHIST_H
#define SAMPLEHIST_H

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

#include "TtHFitter/Common.h"

#include "TtHFitter/SystematicHist.h"
#include "TtHFitter/Sample.h"
#include "TtHFitter/HistoTools.h"

class SampleHist {
public:
    SampleHist();
    SampleHist(Sample *sample,TH1 *hist);
    SampleHist(Sample *sample, std::string histoName, std::string fileName);
    ~SampleHist();

    TH1* GetHist();
    Sample* GetSample();
    SystematicHist* AddOverallSyst(std::string name,float up,float down);
    SystematicHist* AddStatSyst(std::string name,int i_bin);
    SystematicHist* AddHistoSyst(std::string name,TH1* h_up,TH1* h_down);
    SystematicHist* AddHistoSyst(std::string name,std::string histoName_up, std::string fileName_up,std::string histoName_down, std::string fileName_down, int pruned=0);
    SystematicHist* GetSystematic(std::string systName);
    SystematicHist* GetSystFromNP(std::string NuisParName);
    NormFactor* AddNormFactor(std::string name,float nominal, float min, float max);
    NormFactor* AddNormFactor(NormFactor *normFactor);
    NormFactor* GetNormFactor(std::string name);
    ShapeFactor* AddShapeFactor(std::string name,float nominal, float min, float max);
    ShapeFactor* AddShapeFactor(ShapeFactor *shapeFactor);
    ShapeFactor* GetShapeFactor(std::string name);

    bool HasSyst(std::string name);
    bool HasNorm(std::string name);
    bool HasShapeFactor(std::string name);

    void WriteToFile(TFile *f=0x0);
    void ReadFromFile();

    void FixEmptyBins(const bool suppress);

    void Print();

    void Rebin(int ngroup = 2, const Double_t* xbins = 0);
    void DrawSystPlot( const std::string &syst="all", TH1* h_data=0x0, bool SumAndData=false, bool bothPanels=false );
    void SmoothSyst(const HistoTools::SmoothOption &opt, std::string syst="all",bool force=false, bool TtresSmoothing = false);

    void Divide(  SampleHist* sh);
    void Multiply(SampleHist* sh);
    void Add(     SampleHist* sh,float scale=1.);
    void Scale(float scale);

    void SampleHistAdd(SampleHist* h, float scale = 1.);
    void CloneSampleHist(SampleHist* h, std::set<std::string> names, float scale = 1.);

    std::string fName;
    Sample *fSample;
    TH1 *fHist;
    TH1 *fHist_orig;  // new
    TH1 *fHist_regBin;  // new
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

