#ifndef UNFOLDINGSYSTEMATIC_H_
#define UNFOLDINGSYSTEMATIC_H_

#include "TRExFitter/HistoTools.h"

#include <string>
#include <vector>

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
    inline int GetType() const {return fType;}
    inline void SetSymmetrisationType(const HistoTools::SymmetrizationType type) {fSymmetrisationType = type;}
    inline int GetSymmetrisationType() const {return fSymmetrisationType;}
    inline void SetSmoothOption(const HistoTools::SmoothOption opt) {fSampleSmoothingOption = opt;}
    inline int GetSmoothOption() const {return fSampleSmoothingOption;}

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

};

#endif
