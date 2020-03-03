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
    inline void SetHasResponse(const bool r) {fHasResponse = r;}
    inline bool GetHasResponse() const {return fHasResponse;}
    inline void SetHasAcceptance(const bool a) {fHasAcceptance = a;}
    inline bool GetHasAcceptance() const {return fHasAcceptance;}
    

    std::vector<Sample*> ConvertToSample(const Region* reg,
                                         const int bins,
                                         const std::string& name) const;
    
    std::vector<std::string> fResponseMatrixFiles;
    std::vector<std::string> fResponseMatrixNames;
    std::vector<std::string> fResponseMatrixPaths;
    std::vector<std::string> fResponseMatrixFileSuffs;
    std::vector<std::string> fResponseMatrixNameSuffs;
    std::vector<std::string> fResponseMatrixPathSuffs;
    
    std::vector<std::string> fAcceptanceFiles;
    std::vector<std::string> fAcceptanceNames;
    std::vector<std::string> fAcceptancePaths;
    std::vector<std::string> fAcceptanceFileSuffs;
    std::vector<std::string> fAcceptanceNameSuffs;
    std::vector<std::string> fAcceptancePathSuffs;
    
    std::vector<std::string> fSelectionEffFiles;
    std::vector<std::string> fSelectionEffNames;
    std::vector<std::string> fSelectionEffPaths;
    std::vector<std::string> fSelectionEffFileSuffs;
    std::vector<std::string> fSelectionEffNameSuffs;
    std::vector<std::string> fSelectionEffPathSuffs;
    
    std::vector<std::string> fMigrationFiles;
    std::vector<std::string> fMigrationNames;
    std::vector<std::string> fMigrationPaths;
    std::vector<std::string> fMigrationFileSuffs;
    std::vector<std::string> fMigrationNameSuffs;
    std::vector<std::string> fMigrationPathSuffs;
    
    std::vector<std::string> fRegions;

private:
    std::string fName;
    std::string fTitle;
    int fFillColor;
    int fLineColor;
    bool fHasResponse;
    bool fHasAcceptance;
};

#endif
