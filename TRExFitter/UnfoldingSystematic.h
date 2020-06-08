#ifndef UNFOLDINGSYSTEMATIC_H_
#define UNFOLDINGSYSTEMATIC_H_

#include "TRExFitter/HistoTools.h"

#include <memory>
#include <string>
#include <vector>

class Region;
class Sample;
class Systematic;

class UnfoldingSystematic {
public:

    explicit UnfoldingSystematic();
    ~UnfoldingSystematic() = default;

    UnfoldingSystematic(const UnfoldingSystematic& s) = default;
    UnfoldingSystematic& operator=(const UnfoldingSystematic& s) = default;
    UnfoldingSystematic(UnfoldingSystematic&& s) = default;
    UnfoldingSystematic& operator=(UnfoldingSystematic&& s) = default;

    inline void SetName(const std::string& name) {fName = name;}
    inline const std::string& GetName() const {return fName;}
    inline void SetTitle(const std::string& title) {fTitle = title;}
    inline const std::string& GetTitle() const {return fTitle;}
    inline void SetType(const int type) {fType = type;}
    inline void SetSymmetrisationType(const HistoTools::SymmetrizationType type) {fSymmetrisationType = type;}
    inline void SetSmoothOption(const HistoTools::SmoothOption opt) {fSampleSmoothingOption = opt;}
    inline void SetHasResponse(const bool r) {fHasResponse = r;}
    inline bool GetHasResponse() const {return fHasResponse;}
    inline void SetHasAcceptance(const bool a) {fHasAcceptance = a;}
    inline bool GetHasAcceptance() const {return fHasAcceptance;}

    std::vector<std::shared_ptr<Systematic> > ConvertToSystematic(const Region* reg,
                                                                  const int bins,
                                                                  const std::string& name,
                                                                  const std::string& unfoldingSampleName,
                                                                  std::vector<std::shared_ptr<Sample> >& samples) const;

    std::vector<std::string> fRegions;
    std::vector<std::string> fSamples;
    std::string fCategory;
    std::string fSubCategory;
    std::vector<std::string> fResponseMatrixPathsUp;
    std::vector<std::string> fResponseMatrixPathsDown;
    std::vector<std::string> fResponseMatrixNamesUp;
    std::vector<std::string> fResponseMatrixNamesDown;
    std::vector<std::string> fResponseMatrixFilesUp;
    std::vector<std::string> fResponseMatrixFilesDown;
    std::vector<std::string> fResponseMatrixPathSuffsUp;
    std::vector<std::string> fResponseMatrixPathSuffsDown;
    std::vector<std::string> fResponseMatrixNameSuffsUp;
    std::vector<std::string> fResponseMatrixNameSuffsDown;
    std::vector<std::string> fResponseMatrixFileSuffsUp;
    std::vector<std::string> fResponseMatrixFileSuffsDown;
    std::vector<std::string> fAcceptancePathsUp;
    std::vector<std::string> fAcceptancePathsDown;
    std::vector<std::string> fAcceptanceNamesUp;
    std::vector<std::string> fAcceptanceNamesDown;
    std::vector<std::string> fAcceptanceFilesUp;
    std::vector<std::string> fAcceptanceFilesDown;
    std::vector<std::string> fAcceptancePathSuffsUp;
    std::vector<std::string> fAcceptancePathSuffsDown;
    std::vector<std::string> fAcceptanceNameSuffsUp;
    std::vector<std::string> fAcceptanceNameSuffsDown;
    std::vector<std::string> fAcceptanceFileSuffsUp;
    std::vector<std::string> fAcceptanceFileSuffsDown;
    std::vector<std::string> fSelectionEffPathsUp;
    std::vector<std::string> fSelectionEffPathsDown;
    std::vector<std::string> fSelectionEffNamesUp;
    std::vector<std::string> fSelectionEffNamesDown;
    std::vector<std::string> fSelectionEffFilesUp;
    std::vector<std::string> fSelectionEffFilesDown;
    std::vector<std::string> fSelectionEffPathSuffsUp;
    std::vector<std::string> fSelectionEffPathSuffsDown;
    std::vector<std::string> fSelectionEffNameSuffsUp;
    std::vector<std::string> fSelectionEffNameSuffsDown;
    std::vector<std::string> fSelectionEffFileSuffsUp;
    std::vector<std::string> fSelectionEffFileSuffsDown;
    std::vector<std::string> fMigrationPathsUp;
    std::vector<std::string> fMigrationPathsDown;
    std::vector<std::string> fMigrationNamesUp;
    std::vector<std::string> fMigrationNamesDown;
    std::vector<std::string> fMigrationFilesUp;
    std::vector<std::string> fMigrationFilesDown;
    std::vector<std::string> fMigrationPathSuffsUp;
    std::vector<std::string> fMigrationPathSuffsDown;
    std::vector<std::string> fMigrationNameSuffsUp;
    std::vector<std::string> fMigrationNameSuffsDown;
    std::vector<std::string> fMigrationFileSuffsUp;
    std::vector<std::string> fMigrationFileSuffsDown;

    bool fHasUpVariation;
    bool fHasDownVariation;
    bool fSampleSmoothing;
    std::string fNuisanceParameter;

private:
    std::string fName;
    std::string fTitle;
    int fType;
    HistoTools::SymmetrizationType fSymmetrisationType;
    HistoTools::SmoothOption fSampleSmoothingOption;
    bool fHasResponse;
    bool fHasAcceptance;
};

#endif
