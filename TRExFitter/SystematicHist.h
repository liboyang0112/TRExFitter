#ifndef SYSTEMATICHIST_H
#define SYSTEMATICHIST_H

/// c++ includes 
#include <string>

/// Forwards class declaration
class TFile;
class TH1;
class Systematic;

class SystematicHist {
public:
    SystematicHist(const std::string& name);
    ~SystematicHist();

    void WriteToFile(TFile *f=nullptr,bool reWriteOrig=true) const;
    void ReadFromFile();
    bool IsShape() const;

    void Print() const;

    void Divide(TH1* h);
    void Divide(SystematicHist *syh);
    void Multiply(TH1* h);
    void Multiply(SystematicHist *syh);
    void Add(TH1* h,float scale=1.);
    void Add(SystematicHist *syh,float scale=1.);

    std::string fName;
    Systematic *fSystematic;

    bool fIsOverall;
    bool fIsShape;
    int fSmoothType;
    int fSymmetrisationType;

    bool fShapePruned;
    bool fNormPruned;
    bool fBadShape;
    bool fBadNorm;

    TH1* fHistUp;
    TH1* fHistUp_orig;
    TH1* fHistUp_preSmooth;
    TH1* fHistShapeUp;
    float fNormUp;
    std::string fFileNameUp;
    std::string fHistoNameUp;
    std::string fFileNameShapeUp;
    std::string fHistoNameShapeUp;
    TH1* fHistUp_original;
    TH1* fHistUp_postFit;

    TH1* fHistDown;
    TH1* fHistDown_orig;
    TH1* fHistDown_preSmooth;
    TH1* fHistShapeDown;
    float fNormDown;
    std::string fFileNameDown;
    std::string fHistoNameDown;
    std::string fFileNameShapeDown;
    std::string fHistoNameShapeDown;
    TH1* fHistDown_original;
    TH1* fHistDown_postFit;

    float fScaleUp;
    float fScaleDown;
};

#endif
