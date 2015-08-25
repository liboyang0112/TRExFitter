#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

#include "TtHFitter/Common.h"
#include "TtHFitter/Systematic.h"
#include "TtHFitter/NormFactor.h"

#ifndef __Sample__
#define __Sample__

class Sample {
public:
  
    enum SampleType{
        BACKGROUND, // 0
        SIGNAL, // 1
        DATA // 2
    };
    
    Sample(string name,int type=0);
    ~Sample();
    
    // -------
    // Methods
    // -------

    // comsetics
    void SetTitle(string title);
    void SetFillColor(int color);
    void SetLineColor(int color);
    void NormalizedByTheory(const bool norm);
    
    // read from ntupes
    void AddNtuplePath(string path);
    void AddNtupleFile(string file);
    void AddNtupleName(string name);
    void SetMCweight(string weight);
    void SetSelection(string selection);

    // read from histos
    void AddHistoPath(string path);
    void AddHistoFile(string file);
    void AddHistoName(string name);
    
    // norm factors and systs
    void AddNormFactor(NormFactor *factor);
    void AddSystematic(Systematic *syst);
    NormFactor* AddNormFactor(string name,float nominal,float min,float max);
    Systematic* AddSystematic(string name,int type=0,float up=0,float down=0);

    // -------
    // Members
    // -------

    string fName;
    int fType;
    string fFitName;
    string fTitle;
    string fGroup;
    int fFillColor;
    int fLineColor;
    bool fNormalizedByTheory;
    
    // to read from ntuples
    string fSelection;
    string fMCweight;
    std::vector<string> fNtuplePaths;
    std::vector<string> fNtuplePathSuffs;
    std::vector<string> fNtupleFiles;
    std::vector<string> fNtupleFileSuffs;
    std::vector<string> fNtupleNames;
    std::vector<string> fNtupleNameSuffs;
    
    // to read from histograms
    // <path>/<file>.root/<name>
    std::vector<string> fHistoPaths;
    std::vector<string> fHistoPathSuffs;
    std::vector<string> fHistoFiles;
    std::vector<string> fHistoFileSuffs;
    std::vector<string> fHistoNames;
    std::vector<string> fHistoNameSuffs;

    // systematics & norm.factors
    int fNSyst;
    std::vector < Systematic* > fSystematics;
    int fNNorm;
    std::vector < NormFactor* > fNormFactors;
};

#endif

