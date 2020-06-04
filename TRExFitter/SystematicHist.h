#ifndef SYSTEMATICHIST_H
#define SYSTEMATICHIST_H

/// Framework includes
#include "TRExFitter/HistoTools.h"

/// c++ includes
#include <memory>
#include <string>

/// Forwards class declaration
class TFile;
class TH1;
class Systematic;

class SystematicHist {
public:
    explicit SystematicHist(const std::string& name);

    ~SystematicHist();

    SystematicHist(const SystematicHist& s) = delete;
    SystematicHist(SystematicHist&& s) = delete;
    SystematicHist& operator=(const SystematicHist& s) = delete;
    SystematicHist& operator=(SystematicHist&& s) = delete;

    void WriteToFile(TFile *f=nullptr,bool reWriteOrig=true) const;
    void ReadFromFile();
    bool IsShape() const;

    void Print() const;

    void Divide(TH1* h);
    void Divide(SystematicHist *syh);
    void Multiply(TH1* h);
    void Multiply(SystematicHist *syh);
    void Add(TH1* h,double scale=1.);
    void Add(SystematicHist *syh,double scale=1.);

    std::string fName;
    std::shared_ptr<Systematic> fSystematic;

    bool fIsOverall;
    bool fIsShape;
    int fSmoothType;
    HistoTools::SymmetrizationType fSymmetrisationType;

    bool fShapePruned;
    bool fNormPruned;
    bool fBadShape;
    bool fBadNorm;

    std::unique_ptr<TH1> fHistUp;
    std::unique_ptr<TH1> fHistUp_orig;
    std::unique_ptr<TH1> fHistUp_preSmooth;
    std::unique_ptr<TH1> fHistShapeUp;
    double fNormUp;
    std::string fFileNameUp;
    std::string fHistoNameUp;
    std::string fFileNameShapeUp;
    std::string fHistoNameShapeUp;
    std::unique_ptr<TH1> fHistUp_postFit;

    std::unique_ptr<TH1> fHistDown;
    std::unique_ptr<TH1> fHistDown_orig;
    std::unique_ptr<TH1> fHistDown_preSmooth;
    std::unique_ptr<TH1> fHistShapeDown;
    double fNormDown;
    std::string fFileNameDown;
    std::string fHistoNameDown;
    std::string fFileNameShapeDown;
    std::string fHistoNameShapeDown;
    std::unique_ptr<TH1> fHistDown_postFit;

    double fScaleUp;
    double fScaleDown;
};

#endif
