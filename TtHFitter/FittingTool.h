//
// Fit.h
// ==========
// Allows to implement all necessary informations and functions to perform likelihood fits.
//

#ifndef _FittingTool_
#define _FittingTool_

class RooFitResult;
class TString;
class RooAbsPdf;
class RooAbsData;

#include <string>
#include <map>
#include "RooStats/ModelConfig.h"

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
    inline void MinimType ( const TString &type ){ m_minimType = type; }
    inline TString GetMinimType(){ return m_minimType; }
    
    inline int GetMinuitStatus() { return m_minuitStatus; }
    inline int GetHessStatus() { return m_hessStatus; }
    inline double GetEDM() { return m_edm; }
    
    inline void ValPOI( const double value ) { m_valPOI = value; }
    inline double GetValPOI() { return m_valPOI; }
    
    inline void UseMinos( const bool use ) { m_useMinos = use; }
    
    inline void ConstPOI( const bool constant ) { m_constPOI = constant; }
    inline double GetConstPOI() { return m_constPOI; }
    
    inline RooFitResult* GetFitResult() { return m_fitResult; }
    
    //
    // Specific functions
    //
    void FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata );
    void ExportFitResultInTextFile( const std::string &finaName );
    std::map < std::string, double > ExportFitResultInMap();
    
private:
    TString m_minimType;
    int m_minuitStatus, m_hessStatus;
    double m_edm,m_valPOI;
    bool m_useMinos,m_constPOI;
    RooFitResult* m_fitResult;
    bool m_debug;
};


#endif //FitTools
