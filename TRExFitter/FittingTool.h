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
#include <string>
#include <map>
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
    FittingTool();
    ~FittingTool();
    FittingTool( const FittingTool & );

    //
    // Gettters and setters
    //
    inline void SetDebug ( const int debug ){ m_debug = debug; }

    inline void MinimType ( const TString &type ){ m_minimType = type; }
    inline TString GetMinimType(){ return m_minimType; }

    inline int GetMinuitStatus() { return m_minuitStatus; }
    inline int GetHessStatus() { return m_hessStatus; }
    inline double GetEDM() { return m_edm; }

    inline void ValPOI( const double value ) { m_valPOI = value; }
    inline double GetValPOI() { return m_valPOI; }

    //inline void UseMinos( const bool use ) { m_useMinos = use; }

    inline void ConstPOI( const bool constant ) { m_constPOI = constant; }
    inline double GetConstPOI() { return m_constPOI; }

    inline void NoGammas()      { m_noGammas=true;      }
    inline void NoSystematics() { m_noSystematics=true; }
    inline void NoNormFators()  { m_noNormFactors=true; }
    inline void NoShapeFators() { m_noShapeFactors=true; }

    inline void SetRandomNP( const double rndNP, const bool rndize, const long int rndSeed=-999 ) { m_randomNP = rndNP; m_randomize = rndize; m_randSeed = rndSeed; }

    void SetSubCategories();
    inline void SetSystMap( std::map<std::string, std::string> subCategoryMap ) { m_subCategoryMap = subCategoryMap; SetSubCategories(); } // fills both m_subCategoryMap and m_subCategories

//     inline void FixNP( const TString &np, const double value ) { m_constNP.push_back(np); m_constNPvalue.push_back(value); }
    inline void ResetFixedNP() { m_constNP.clear(); m_constNPvalue.clear(); };
    inline void FixNP( std::string np, const double value ) { m_constNP.push_back(np); m_constNPvalue.push_back(value); }
    inline void FixNPs( std::vector<std::string> nps, std::vector<double> values ) { m_constNP = nps; m_constNPvalue = values; }
    inline void SetNPs( std::vector<std::string> nps, std::vector<double> values ) { m_initialNP = nps; m_initialNPvalue = values; }

    inline RooFitResult* GetFitResult() { return m_fitResult; }

    inline void RangePOI_up( const double value){m_RangePOI_up = value;}
    inline void RangePOI_down( const double value){m_RangePOI_down = value;}

    inline void UseMinos( const std::vector<std::string> minosvar){ m_useMinos = true; m_varMinos = minosvar; }

    //
    // Specific functions
    //
    double FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit = false, bool noFit = false );
    void ExportFitResultInTextFile( const std::string &finaName );
    std::map < std::string, double > ExportFitResultInMap();

    int GetGroupedImpact( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, RooWorkspace* ws, std::string categoryOfInterest, std::string outFileName );
    void FitExcludingGroup(bool excludeGammas, bool statOnly, RooAbsData*& fitdata, RooAbsPdf*& fitpdf, RooArgSet*& constrainedParams,
                           RooStats::ModelConfig* mc, RooWorkspace* ws, std::string category, std::vector<std::string> affectedParams);

private:
    TString m_minimType;
    int m_minuitStatus;
    int m_hessStatus;
    double m_edm;
    double m_valPOI;
    bool m_useMinos;
    std::vector<std::string> m_varMinos;
    bool m_constPOI;
    RooFitResult* m_fitResult;
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

//     TString m_constNP;
//     double m_constNPvalue;
    std::vector<std::string> m_constNP;
    std::vector<double> m_constNPvalue;
    std::vector<std::string> m_initialNP;
    std::vector<double> m_initialNPvalue;

    std::map<std::string, std::string> m_subCategoryMap;
    std::set<std::string> m_subCategories;
};


#endif //FitTools
