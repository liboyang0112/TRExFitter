#ifndef TRUTHSAMPLE_H_
#define TRUTHSAMPLE_H_

#include <memory>
#include <string>

class TH1;
class TRExFit;

class TruthSample {

public:

    explicit TruthSample(const std::string& name);
    ~TruthSample() = default;

    TruthSample(const TruthSample& s) = default;
    TruthSample& operator=(const TruthSample& s) = default;
    TruthSample(TruthSample&& s) = default;
    TruthSample& operator=(TruthSample&& s) = default;

    inline const std::string& GetName() const {return fName;};
    inline void SetTitle(const std::string& s) {fTitle = s;}
    inline const std::string& GetTitle() const {return fTitle;}
    inline void SetFillColor(const int c) {fFillColor = c;}
    inline int GetFillColor() const {return fFillColor;}
    inline void SetLineColor(const int c) {fLineColor = c;}
    inline int GetLineColor() const {return fLineColor;}
    inline void SetTruthDistributionPath(const std::string& s) {fTruthDistributionPath = s;}
    inline void SetTruthDistributionFile(const std::string& s) {fTruthDistributionFile = s;}
    inline void SetTruthDistributionName(const std::string& s) {fTruthDistributionName = s;}

    std::unique_ptr<TH1> GetHisto(const TRExFit* fitter) const;

private:

    std::string fName;
    std::string fTitle;
    int fFillColor;
    int fLineColor;
    std::string fTruthDistributionPath;
    std::string fTruthDistributionFile;
    std::string fTruthDistributionName;
};

#endif
