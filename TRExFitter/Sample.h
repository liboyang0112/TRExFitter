#ifndef SAMPLE_H
#define SAMPLE_H

/// c++ includes
#include <map>
#include <string>
#include <vector>

/// Forward class declaration
class NormFactor;
class ShapeFactor;
class Systematic;

class Sample {
public:

    enum SampleType{
        BACKGROUND, // 0
        SIGNAL, // 1
        DATA, // 2
        GHOST // 3
    };

    Sample(const std::string& name,int type=0);
    ~Sample();

    // -------
    // Methods
    // -------

    // comsetics
    void SetTitle(const std::string& title);
    void SetFillColor(int color);
    void SetLineColor(int color);
    void SetFillColorRGB(const std::array<int, 3>& rgb);
    void SetLineColorRGB(const std::array<int, 3>& rgb);
    void NormalizedByTheory(const bool norm);

    // read from ntupes
    void AddNtuplePath(const std::string& path);
    void AddNtupleFile(const std::string& file);
    void AddNtupleName(const std::string& name);
    void SetMCweight(const std::string& weight);
    void SetSelection(const std::string& selection);

    // read from histos
    void AddHistoPath(const std::string& path);
    void AddHistoFile(const std::string& file);
    void AddHistoName(const std::string& name);

    // norm factors and systs
    void AddNormFactor(NormFactor *factor);
    void AddShapeFactor(ShapeFactor *factor);
    void AddSystematic(Systematic *syst);
    NormFactor* AddNormFactor(const std::string& name,float nominal=1,float min=0,float max=10,bool isConst=false);
    ShapeFactor* AddShapeFactor(const std::string& name,float nominal=1,float min=0,float max=10,bool isConst=false);
    Systematic* AddSystematic(const std::string& name,int type=0,float up=0,float down=0);
    bool HasNormFactor(const std::string& name) const;
    bool HasSystematic(const std::string& name) const;

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
    std::array<int, 3> fFillColorRGB;
    std::array<int, 3> fLineColorRGB;
    bool fNormalizedByTheory;
    std::vector<std::string> fRegions;
    std::vector<float> fLumiScales;
    std::string fIgnoreSelection;
    std::string fIgnoreWeight;
    bool fUseMCStat;
    bool fUseSystematics;
    std::string fDivideBy;
    std::string fMultiplyBy;
    std::vector<std::string> fSubtractSamples;
    std::vector<std::string> fAddSamples;
    std::string fNormToSample;
    bool fSmooth;
    int fBuildPullTable;
    std::map<std::string,bool> fIsMorph;
    std::map<std::string,float> fMorphValue;

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
    std::vector<std::vector<std::string>> fCorrelateGammasInRegions;
};

#endif

