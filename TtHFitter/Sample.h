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
#include "TtHFitter/ShapeFactor.h"

#ifndef __Sample__
#define __Sample__

class Sample {
public:
  
    enum SampleType{
        BACKGROUND, // 0
        SIGNAL, // 1
        DATA, // 2
        GHOST // 3
    };
    
    Sample(std::string name,int type=0);
    ~Sample();
    
    // -------
    // Methods
    // -------

    // comsetics
    void SetTitle(std::string title);
    void SetFillColor(int color);
    void SetLineColor(int color);
    void NormalizedByTheory(const bool norm);
    
    // read from ntupes
    void AddNtuplePath(std::string path);
    void AddNtupleFile(std::string file);
    void AddNtupleName(std::string name);
    void SetMCweight(std::string weight);
    void SetSelection(std::string selection);

    // read from histos
    void AddHistoPath(std::string path);
    void AddHistoFile(std::string file);
    void AddHistoName(std::string name);
    
    // norm factors and systs
    void AddNormFactor(NormFactor *factor);
    void AddShapeFactor(ShapeFactor *factor);
    void AddSystematic(Systematic *syst);
    NormFactor* AddNormFactor(std::string name,float nominal=1,float min=0,float max=10,bool isConst=false);
    ShapeFactor* AddShapeFactor(string name,float nominal=1,float min=0,float max=10,bool isConst=false);
    Systematic* AddSystematic(std::string name,int type=0,float up=0,float down=0);
    bool HasNormFactor(std::string name);
    bool HasSystematic(std::string name);
    
    // -------
    // Members
    // -------

    std::string fName;
    int fType;
    std::string fFitName;
    std::string fTitle;
    std::string fTexTitle;
    std::string fGroup;
    int fFillColor;
    int fLineColor;
    bool fNormalizedByTheory;
    std::vector<std::string> fRegions;
    std::vector<float> fLumiScales;
    std::string fIgnoreSelection;
    bool fUseMCStat;
    bool fUseSystematics;
    std::string fDivideBy;
    std::string fMultiplyBy;
    std::vector<std::string> fSubtractSamples;
    std::vector<std::string> fAddSamples;
    bool fSmooth;
    
    // to read from ntuples
    std::string fSelection;
    std::string fMCweight;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtuplePathSuffs;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleFileSuffs;
    std::vector<std::string> fNtupleNames;
    std::vector<std::string> fNtupleNameSuffs;
    
    // to read from histograms
    // <path>/<file>.root/<name>
    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoPathSuffs;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoFileSuffs;
    std::vector<std::string> fHistoNames;
    std::vector<std::string> fHistoNameSuffs;

    // systematics & norm.factors
    int fNSyst;
    std::vector < Systematic* > fSystematics;
    int fNNorm;
    std::vector < NormFactor* > fNormFactors;
    int fNShape;
    std::vector < ShapeFactor* > fShapeFactors; 
   
    std::pair<std::string,std::string> fAsimovReplacementFor;
    
    bool fSeparateGammas;
};

#endif

