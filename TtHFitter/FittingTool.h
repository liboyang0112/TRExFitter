//
// Fit.h
// ==========
// Allows to implement all necessary informations and functions to perform likelihood fits.
//

#ifndef _FittingTool_
#define _FittingTool_

#include "TtHFitter/Common.h"

class RooFitResult;
class TString;
class RooAbsPdf;
class RooAbsData;

#include <string>
#include <map>
#include <vector>
#include "RooStats/ModelConfig.h"
#include "TVirtualFitter.h"

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
    inline void SetDebug ( const int debug ){ m_debug = debug>0; }

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
    float FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit = false, bool noFit = false );
    void ExportFitResultInTextFile( const std::string &finaName );
    std::map < std::string, double > ExportFitResultInMap();

private:
    TString m_minimType;
    int m_minuitStatus, m_hessStatus;
    double m_edm,m_valPOI,m_randomNP;
    long int m_randSeed;
    double m_RangePOI_up,m_RangePOI_down;
    bool m_useMinos,m_constPOI;
    bool m_noGammas,m_noSystematics,m_noNormFactors;
    bool m_noShapeFactors;
    std::vector<std::string> m_varMinos;
    RooFitResult* m_fitResult;
    bool m_debug,m_randomize;

//     TString m_constNP;
//     double m_constNPvalue;
    std::vector<std::string> m_constNP;
    std::vector<double> m_constNPvalue;
    std::vector<std::string> m_initialNP;
    std::vector<double> m_initialNPvalue;
};


#endif //FitTools
