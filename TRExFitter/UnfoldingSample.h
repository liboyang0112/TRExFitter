#ifndef UNFOLDINGSAMPLE_H_
#define UNFOLDINGSAMPLE_H_

#include <string>
#include <vector>

class Region;
class Sample;

class UnfoldingSample {

public:
    explicit UnfoldingSample();
    ~UnfoldingSample() = default;

    UnfoldingSample(const UnfoldingSample& s) = default;
    UnfoldingSample(UnfoldingSample&& s) = default;
    UnfoldingSample& operator=(const UnfoldingSample& s) = default;
    UnfoldingSample& operator=(UnfoldingSample&& s) = default;

    inline void SetName(const std::string& s) {fName = s;}
    inline const std::string& GetName() const {return fName;}
    inline void SetTitle(const std::string& s) {fTitle = s;}
    inline const std::string& GetTitle() const {return fTitle;}
    inline void SetFillColor(const int c) {fFillColor = c;}
    inline int GetFillColor() const {return fFillColor;}
    inline void SetLineColor(const int c) {fLineColor = c;}
    inline int GetLineColor() const {return fLineColor;}

    Sample* ConvertToSample(const Region* reg) const;
    
    std::vector<std::string> fResponseMatrixFiles;
    std::vector<std::string> fResponseMatrixNames;
    std::vector<std::string> fResponseMatrixPaths;
    std::vector<std::string> fResponseMatrixFileSuffs;
    std::vector<std::string> fResponseMatrixNameSuffs;
    std::vector<std::string> fResponseMatrixPathSuffs;
    std::vector<std::string> fRegions;

private:
    std::string fName;
    std::string fTitle;
    int fFillColor;
    int fLineColor;
};

#endif
