//
// Fit.h
// ==========
// Allows to implement all necessary informations and functions to perform likelihood fits.
//

#ifndef FITTINGTOOL_H
#define FITTINGTOOL_H

/// RooStats include
#include "RooStats/ModelConfig.h"

// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forward declaration
class RooFitResult;
class TString;
class RooAbsPdf;
class RooAbsData;

class FittingTool {

public:

    //
    // Standard C++ functions
    //
    explicit FittingTool();
    ~FittingTool();
    FittingTool(const FittingTool& rhs) = default;
    FittingTool& operator=(const FittingTool& rhs) = default;
    FittingTool(FittingTool&& rhs) = default;
    FittingTool& operator=(FittingTool&& rhs) = default;

    //
    // Gettters and setters
    //
    inline void SetDebug (const int debug){ m_debug = debug; }

    inline void MinimType (const TString& type){ m_minimType = type; }

    inline void ValPOI(const double value) { m_valPOI = value; }

    inline void ConstPOI(const bool constant) { m_constPOI = constant; }

    inline void NoGammas()      { m_noGammas=true;      }
    inline void NoSystematics() { m_noSystematics=true; }

    inline void SetRandomNP(const double rndNP,
                             const bool rndize,
                             const long int rndSeed=-999) { m_randomNP = rndNP; m_randomize = rndize; m_randSeed = rndSeed; }

    void SetSubCategories();
    void SetSystMap(const std::map<std::string, std::string>& subCategoryMap) { m_subCategoryMap = subCategoryMap; SetSubCategories(); } // fills both m_subCategoryMap and m_subCategories

    inline void ResetFixedNP() { m_constNP.clear(); m_constNPvalue.clear(); };
    inline void FixNP(const std::string& np, const double value) { m_constNP.emplace_back(np); m_constNPvalue.emplace_back(value); }
    inline void FixNPs(const std::vector<std::string>& nps, const std::vector<double>& values) { m_constNP = nps; m_constNPvalue = values; }
    inline void SetNPs(const std::vector<std::string>& nps, const std::vector<double>& values) { m_initialNP = nps; m_initialNPvalue = values; }

    inline void UseMinos(const std::vector<std::string>& minosvar){ m_useMinos = true; m_varMinos = minosvar; }

    inline void SetExternalConstraints(const RooArgSet* externalConstraints = 0){ m_externalConstraints = externalConstraints; }

    inline void SetStrategy(const int strategy){m_strategy = strategy;}
    
    inline void SetUseHesse(const bool flag){m_useHesse = flag;}

    //
    // Specific functions
    //
    double FitPDF(RooStats::ModelConfig* model,
                  RooAbsPdf* fitpdf,
                  RooAbsData* fitdata,
                  bool fastFit = false,
                  bool noFit = false,
                  bool saturatedModel = false );

    void SaveFitResult(const std::string& fileName);

    void ExportFitResultInTextFile(const std::string& finalName,
                                   const std::vector<std::string>& blinded);

    std::map < std::string, double > ExportFitResultInMap();

    int GetGroupedImpact(RooStats::ModelConfig* model,
                         RooAbsPdf* fitpdf,
                         RooAbsData* fitdata,
                         RooWorkspace* ws,
                         const std::string& categoryOfInterest,
                         const std::string& outFileName ) const;

    void FitExcludingGroup(bool excludeGammas,
                           bool statOnly,
                           RooAbsData*& fitdata,
                           RooAbsPdf*& fitpdf,
                           RooArgSet*& constrainedParams,
                           RooStats::ModelConfig* mc,
                           RooWorkspace* ws,
                           const std::string& category,
                           const std::vector<std::string>& affectedParams) const;

    void CheckUnderconstraint(const RooRealVar* const var) const;

private:
    TString m_minimType;
    int m_minuitStatus;
    int m_hessStatus;
    double m_edm;
    double m_valPOI;
    bool m_useMinos;
    std::vector<std::string> m_varMinos;
    bool m_constPOI;
    std::unique_ptr<RooFitResult> m_fitResult;
    int m_debug;
    bool m_noGammas;
    bool m_noSystematics;
    bool m_noNormFactors;
    bool m_noShapeFactors;
    double m_RangePOI_up;
    double m_RangePOI_down;
    bool m_randomize;
    double m_randomNP;
    long int m_randSeed;

    std::vector<std::string> m_constNP;
    std::vector<double> m_constNPvalue;
    std::vector<std::string> m_initialNP;
    std::vector<double> m_initialNPvalue;

    std::map<std::string, std::string> m_subCategoryMap;
    std::set<std::string> m_subCategories;

    const RooArgSet* m_externalConstraints;
    int m_strategy;
    bool m_useHesse;
};


#endif //FitTools
