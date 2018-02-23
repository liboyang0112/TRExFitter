#ifndef SYSTEMATICHIST_H
#define SYSTEMATICHIST_H

#include "TtHFitter/Common.h"
#include "TtHFitter/Systematic.h"

class SystematicHist {
public:
    SystematicHist(std::string name);
    ~SystematicHist();

    void WriteToFile(TFile *f=0x0);
    void ReadFromFile();
    bool IsShape();

    void Print();

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
