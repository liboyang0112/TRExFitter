//
// Tool extracted from FitCrossCheckForLimits.C tool provided by the ATLAS Statistical Forum
// https://svnweb.cern.ch/cern/wsvn/atlasphys/Physics/StatForum/NuisanceCheck/trunk/FitCrossCheckForLimits.C
//

//Standard headers
#include <iostream>
#include <fstream>

//Root headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"

//Roostats headers
#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"

//Roofit headers
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooArgSet.h"

//TtHFitter includes
#include "TtHFitter/FittingTool.h"

using namespace std;

const bool debug = false;


//________________________________________________________________________
//
FittingTool::FittingTool():
m_minimType("Minuit2"),
m_minuitStatus(-1),
m_hessStatus(-1),
m_edm(-1.),
m_valPOI(0.),
m_useMinos(false),
m_constPOI(false),
m_fitResult(0),
m_debug(false),
m_constNP(""),
m_constNPvalue(0.),
m_RangePOI_up(100.),
m_RangePOI_down(-10.)
{}

//________________________________________________________________________
//
FittingTool::FittingTool( const FittingTool &q ){
    m_minimType     = q.m_minimType;
    m_minuitStatus  = q.m_minuitStatus;
    m_hessStatus    = q.m_hessStatus;
    m_edm           = q.m_edm;
    m_valPOI        = q.m_valPOI;
    m_useMinos      = q.m_useMinos;
    m_constPOI      = q.m_constPOI;
    m_fitResult     = q.m_fitResult;
    m_debug         = q.m_debug;
    m_RangePOI_up   = q.m_RangePOI_up;
    m_RangePOI_down = q.m_RangePOI_down;
}

//________________________________________________________________________
//
FittingTool::~FittingTool()
{}

//________________________________________________________________________
//
void FittingTool::FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata ) {
    
    if(m_debug) std::cout << "-> Entering in FitPDF function" << std::endl;
    
    //
    // Printing the whole model for information
    //
    if(m_debug) model->Print();
    
    //
    // Getting the list of model that can be constrained (parameters of the MC)
    //
    RooArgSet* constrainedParams = fitpdf->getParameters(*fitdata);
    RooStats::RemoveConstantParameters(constrainedParams);
    RooFit::Constrain(*constrainedParams);
    
    //
    // Get the global observables (nominal values)
    //
    const RooArgSet* glbObs = model->GetGlobalObservables();
    
    //
    // Getting the POI
    //
    RooRealVar * poi = (RooRealVar*) model->GetParametersOfInterest()->first();
    if(!poi){
        std::cout << "<!> In FittingTool::FitPDF(): Cannot find the parameter of interest !" << std::endl;
        return;
    }
    
    poi -> setVal(m_valPOI);
    poi -> setConstant(m_constPOI);
    //poi -> setRange(m_RangePOI_down,m_RangePOI_up); // Commented by Loic to avoid overwriting user's setting in config file
    
    if(m_debug){
        std::cout << "   -> Constant POI : " << poi->isConstant() << std::endl;
        std::cout << "   -> Value of POI : " << poi->getVal()     << std::endl;
    }
    
    //
    // Needed for Ranking plot
    //
    RooRealVar* var = NULL;
    RooArgSet* nuis = (RooArgSet*) model->GetNuisanceParameters();
    if(nuis){
        TIterator* it2 = nuis->createIterator();
        while( (var = (RooRealVar*) it2->Next()) ){
            string np = var->GetName();
            if( np == ("alpha_"+m_constNP) ){
                var->setVal(m_constNPvalue);
                var->setConstant(1);
            } else if( np.find("alpha_")!=string::npos ){
                var->setVal(0);
                var->setConstant(0);
            }
        }
    }
    
    //
    // Create the likelihood based on fitpdf, fitData and the parameters
    //
    RooAbsReal * nll = fitpdf->createNLL(*fitdata, RooFit::Constrain(*constrainedParams), RooFit::GlobalObservables(*glbObs), RooFit::Offset(1) );
    
    const double nllval = nll->getVal();
    if(m_debug){
        std::cout << "   -> Initial value of the NLL = " << nllval << std::endl;
        constrainedParams->Print("v");
    }
    
    //
    //
    // Safe fit loop
    //
    //
    static int nrItr = 0;
    const int maxRetries = 3;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(m_minimType);
    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    int save_strat = strat;
    RooMinimizer minim(*nll);
    minim.setStrategy(strat);
    minim.setPrintLevel(1);
    minim.setEps(1);
    
    TStopwatch sw; sw.Start();
    
    int status=-99;
    m_hessStatus=-99;
    m_edm = -99;
    RooFitResult * r;
    
    while (nrItr<maxRetries && status!=0 && status!=1){
        
        cout << endl;
        cout << endl;
        cout << endl;
        cout << "Fit try n°" << nrItr+1 << endl;
        cout << "======================" << endl;
        cout << endl;
        
        ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        m_hessStatus= minim.hesse();
        r = minim.save();
        m_edm = r->edm();
        
        //up the strategy
        bool FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if (FitIsNotGood && strat<2){
            cout << endl;
            cout << "   *******************************" << endl;
            cout << "   * Increasing Minuit strategy (was " << strat << ")" << endl;
            strat++;
            cout << "   * Fit failed with : " << endl;
            cout << "      - minuit status " << status << endl;
            cout << "      - hess status " << m_hessStatus << endl;
            cout << "      - Edm = " << m_edm << endl;
            cout << "   * Retrying with strategy " << strat << endl;
            cout << "   ********************************" << endl;
            cout << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            m_hessStatus= minim.hesse();
            r = minim.save();
            m_edm = r->edm();
        }
        
        FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if (FitIsNotGood && strat < 2){
            cout << endl;
            cout << "   ********************************" << endl;
            cout << "   * Increasing Minuit strategy (was " << strat << ")" << endl;
            strat++;
            cout << "   * Fit failed with : " << endl;
            cout << "      - minuit status " << status << endl;
            cout << "      - hess status " << m_hessStatus << endl;
            cout << "      - Edm = " << m_edm << endl;
            cout << "   * Retrying with strategy " << strat << endl;
            cout << "   ********************************" << endl;
            cout << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            r = minim.save();
            m_edm = r->edm();
        }
        
        FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if (FitIsNotGood && strat < 2){
            cout << endl;
            cout << "   *******************************" << endl;
            cout << "   * Increasing Minuit strategy (was " << strat << ")" << endl;
            strat++;
            cout << "   * Fit failed with : " << endl;
            cout << "      - minuit status " << status << endl;
            cout << "      - hess status " << m_hessStatus << endl;
            cout << "      - Edm = " << m_edm << endl;
            cout << "   * Retrying with strategy " << strat << endl;
            cout << "   ********************************" << endl;
            cout << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            m_hessStatus= minim.hesse();
            r = minim.save();
            m_edm = r->edm();
        }
        
        if(m_useMinos) { minim.minos(); }
        
        FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if ( FitIsNotGood ) nrItr++;
        if (nrItr == maxRetries) {
            cout << endl;
            cout << endl;
            cout << endl;
            cout << "***********************************************************" << endl;
            cout << "WARNING::Fit failure unresolved with status " << status << endl;
            cout << "   Please investigate your workspace" << endl;
            cout << "   Find a wall : you will need it to crash your head on it" << endl;
            cout << "***********************************************************" << endl;
            cout << endl;
            cout << endl;
            cout << endl;
            m_minuitStatus = status;
            m_fitResult = 0;
            return;
        }
        
    }
    
    r = minim.save();
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "***********************************************************" << endl;
    cout << "         FIT FINALIZED SUCCESSFULLY : " << endl;
    cout << "            - minuit status " << status << endl;
    cout << "            - hess status " << m_hessStatus << endl;
    cout << "            - Edm = " << m_edm << endl;
    cout << " -- " ; sw.Print();
    cout << "***********************************************************" << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    
    m_minuitStatus = status;
    m_fitResult = r;

}

//____________________________________________________________________________________
//
void FittingTool::ExportFitResultInTextFile( const std::string &fileName )
{
    if(!m_fitResult){
        std::cerr << "<!> ERROR in FittingTool::ExportFitResultInTextFile(): the FitResultObject seems not to be defined." << std::endl;
    }
    
    //
    // Printing the nuisance parameters post-fit values
    //
    ofstream nuisParAndCorr(fileName);
    nuisParAndCorr << "NUISANCE_PARAMETERS" << std::endl;
    
    RooRealVar* var(nullptr);
    TIterator* param = m_fitResult -> floatParsFinal().createIterator();
    while( (var = (RooRealVar*) param->Next()) ){
        
        // Not consider nuisance parameter being not associated to syst (yet)
        string varname = (string) var->GetName();
        //if ((varname.find("gamma_stat")!=string::npos)) continue;
        TString vname=var->GetName();
        vname.ReplaceAll("alpha_","");
        
        double pull  = var->getVal() / 1.0 ; // GetValue() return value in unit of sigma
        double errorHi = var->getErrorHi() / 1.0;
        double errorLo = var->getErrorLo() / 1.0;
        
        nuisParAndCorr << vname << "  " << pull << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "" << endl;
    }
    
    //
    // Correlation matrix
    //
    TH2* h2Dcorrelation = m_fitResult -> correlationHist();
    nuisParAndCorr << endl << endl << "CORRELATION_MATRIX" << endl;
    nuisParAndCorr << h2Dcorrelation->GetNbinsX() << "   " << h2Dcorrelation->GetNbinsY() << endl;
    for(int kk=1; kk < h2Dcorrelation->GetNbinsX()+1; kk++) {
        for(int ll=1; ll < h2Dcorrelation->GetNbinsY()+1; ll++) {
            nuisParAndCorr << h2Dcorrelation->GetBinContent(kk,ll) << "   ";
        }
        nuisParAndCorr << endl;
    }
    
    //
    // Closing the output file
    //
    nuisParAndCorr << endl;
    nuisParAndCorr.close();
}

//____________________________________________________________________________________
//
std::map < std::string, double > FittingTool::ExportFitResultInMap(){
    
    if(!m_fitResult){
        std::cerr << "<!> ERROR in FittingTool::ExportFitResultInMap(): the FitResultObject seems not to be defined." << std::endl;
    }
    std::map < std::string, double > result;
    RooRealVar* var(nullptr);
    TIterator* param = m_fitResult -> floatParsFinal().createIterator();
    while( (var = (RooRealVar*) param->Next()) ){
        // Not consider nuisance parameter being not associated to syst
        string varname = (string) var->GetName();
        double pull  = var->getVal() / 1.0 ;
        result.insert( std::pair < std::string, double >(varname, pull) );
    }
    return result;
}
