#ifndef TRUTHSAMPLE_H_
#define TRUTHSAMPLE_H_

#include <string>

class TruthSample {

public:

    explicit TruthSample(const std::string& name);
    ~TruthSample() = default;

    TruthSample(const TruthSample& s) = default;
    TruthSample& operator=(const TruthSample& s) = default;
    TruthSample(TruthSample&& s) = default;
    TruthSample& operator=(TruthSample&& s) = default;

    inline const std::string& GetName() const {return fName;};
    inline void SetFillColor(const int c) {fFillColor = c;}
    inline int GetFillColor() const {return fFillColor;}
    inline void SetLineColor(const int c) {fLineColor = c;}
    inline int GetLineColor() const {return fLineColor;}
    inline void SetTruthDistributionPath(const std::string& s) {fTruthDistributionPath = s;}
    inline const std::string& GetTruthDistributionPath() const {return fTruthDistributionPath;}
    inline void SetTruthDistributionFile(const std::string& s) {fTruthDistributionFile = s;}
    inline const std::string& GetTruthDistributionFile() const {return fTruthDistributionFile;}
    inline void SetTruthDistributionName(const std::string& s) {fTruthDistributionName = s;}
    inline const std::string& GetTruthDistributionName() const {return fTruthDistributionName;}

private:

    std::string fName;
    int fFillColor;
    int fLineColor;
    std::string fTruthDistributionPath;
    std::string fTruthDistributionFile;
    std::string fTruthDistributionName;
};

#endif
