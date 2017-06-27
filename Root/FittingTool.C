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
#include "TRandom3.h"

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

//________________________________________________________________________
//
FittingTool::FittingTool():
m_minimType("Minuit2"),
m_minuitStatus(-1),
m_hessStatus(-1),
m_edm(-1.),
m_valPOI(0.),
m_useMinos(false),
m_varMinos(0),
m_constPOI(false),
m_fitResult(0),
m_debug(false),
m_noGammas(false),
m_noSystematics(false),
m_noNormFactors(false),
m_noShapeFactors(false),
// m_constNP(""),
// m_constNPvalue(0.),
m_RangePOI_up(100.),
m_RangePOI_down(-10.),
m_randomize(false),
m_randomNP(0.1)
{
    m_constNP.clear();
    m_constNPvalue.clear();
}

//________________________________________________________________________
//
FittingTool::FittingTool( const FittingTool &q ){
    m_minimType      = q.m_minimType;
    m_minuitStatus   = q.m_minuitStatus;
    m_hessStatus     = q.m_hessStatus;
    m_edm            = q.m_edm;
    m_valPOI         = q.m_valPOI;
    m_useMinos       = q.m_useMinos;
    m_varMinos       = q.m_varMinos;
    m_constPOI       = q.m_constPOI;
    m_fitResult      = q.m_fitResult;
    m_debug          = q.m_debug;
    m_RangePOI_up    = q.m_RangePOI_up;
    m_RangePOI_down  = q.m_RangePOI_down;
    m_noGammas       = q.m_noGammas;
    m_noSystematics  = q.m_noSystematics;
    m_noNormFactors  = q.m_noNormFactors;
    m_noShapeFactors = q.m_noShapeFactors;
}

//________________________________________________________________________
//
FittingTool::~FittingTool()
{}

//________________________________________________________________________
//
float FittingTool::FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit, bool noFit ) {
    
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
    // Create the likelihood based on fitpdf, fitData and the parameters
    //
//     RooAbsReal * nll = fitpdf->createNLL(*fitdata, RooFit::Constrain(*constrainedParams), RooFit::GlobalObservables(*glbObs), RooFit::Offset(1) );
    RooAbsReal * nll = fitpdf->createNLL(*fitdata, RooFit::Constrain(*constrainedParams), RooFit::GlobalObservables(*glbObs), RooFit::Offset(1), 
//     RooAbsReal * nll = fitpdf->createNLL(*fitdata, RooFit::Constrain(*constrainedParams), RooFit::GlobalObservables(*glbObs), RooFit::Offset(0), 
                                         RooFit::NumCPU(TtHFitter::NCPU,RooFit::Hybrid)
//                                          ,RooFit::Optimize(2)
//                                          ,RooFit::Extended(true)   // experimental
                                        );
    

    //
    // Needed for Ranking plot, but also to set random initial values for the NPs
    //
    if(m_randSeed == -999){
      gRandom->SetSeed(time(NULL));
    }
    else{
      gRandom->SetSeed(m_randSeed);
    }

    //
    // Getting the POI
    //
    RooRealVar * poi = (RooRealVar*) model->GetParametersOfInterest()->first();
    if(!poi){
        std::cout << "<!> In FittingTool::FitPDF(): Cannot find the parameter of interest !" << std::endl;
        return 0;
    }
    
    poi -> setConstant(m_constPOI);
    //poi -> setRange(m_RangePOI_down,m_RangePOI_up); // Commented by Loic to avoid overwriting user's setting in config file
    
    poi -> setVal(m_valPOI);
    std::cout << "Setting starting mu = " << m_valPOI << std::endl;
    // randomize the POI
    if(!m_constPOI && m_randomize){
        poi->setVal( m_valPOI + m_randomNP*(gRandom->Uniform(2)-1.) );
    }
    
    if(m_debug){
        std::cout << "   -> Constant POI : " << poi->isConstant() << std::endl;
        std::cout << "   -> Value of POI : " << poi->getVal()     << std::endl;
    }
    
    RooRealVar* var = NULL;
    RooArgSet* nuis = (RooArgSet*) model->GetNuisanceParameters();
    if(nuis){
        TIterator* it2 = nuis->createIterator();
        while( (var = (RooRealVar*) it2->Next()) ){
            string np = var->GetName();
//             if( np == ("alpha_"+m_constNP) || np == m_constNP ){
//                 var->setVal(m_constNPvalue);
//                 var->setConstant(1);
            bool found = false;
            //
            // first check if all systs, norm and gammas should be set to constant
            if(np.find("gamma_stat")!=string::npos && m_noGammas){
                if(m_debug) cout << "setting to constant : " << np <<" at value " << var->getVal() << endl;
                var->setConstant(1);
                found = true;
            }
            else if(np.find("alpha_")!=string::npos && m_noSystematics){
                if(m_debug) cout << "setting to constant : " << np <<" at value " << var->getVal() << endl;
                var->setConstant( 1 );
//                 var->setVal( 0 );
                found = true;
            }
            else if(m_noNormFactors){
                if(m_debug) cout << "setting to constant : " << np <<" at value " << var->getVal() << endl;
                var->setConstant( 1 );
//                 var->setVal( 1 );
                found = true;
            }
	    // FIXME SF
//             else if(m_noShapeFactors){
//                 if(m_debug) cout << "setting to constant : " << np <<" at value " << var->getVal() << endl;
//                 var->setConstant( 1 );
// //                 var->setVal( 1 );
//                 found = true;
//             }
            if(found) continue;
            //
            // loop on the NP specified to be constant
            for( unsigned int i_np = 0; i_np<m_constNP.size(); i_np++ ){
                if( np == ("alpha_"+m_constNP[i_np]) || np == m_constNP[i_np]
                    || np == ("gamma_"+m_constNP[i_np])
                ){
                    if(m_debug) cout << "setting to constant : " << np <<" at value " << m_constNPvalue[i_np] << endl;
                    var->setVal(m_constNPvalue[i_np]);
                    var->setConstant(1);
                    found = true;
                    break;
                }
            }
            //
            // loop on the NP specified to have custom starting value
            for( unsigned int i_np = 0; i_np<m_initialNP.size(); i_np++ ){
                if( np == ("alpha_"+m_initialNP[i_np]) || np == m_initialNP[i_np] ){
                    var->setVal(m_initialNPvalue[i_np]);
//                     var->setVal(m_initialNPvalue[i_np]/2.); // why was it like this?
                    std::cout << " ---> Setting " << m_initialNP[i_np] << " to " << m_initialNPvalue[i_np] << std::endl;
                    found = true;
                    break;
                }
            }
            //
            if(!found){
                if( np.find("alpha_")!=string::npos ){   // for syst NP
                    if(m_randomize) var->setVal( m_randomNP*(gRandom->Uniform(2)-1.) );
                    else            var->setVal(0);
                    var->setConstant(0);
                }
                else {  // for norm factors & gammas
                    if(m_randomize) var->setVal( 1 + m_randomNP*(gRandom->Uniform(2)-1.) );
                    else            var->setVal( 1 );
                }
            }
        }
        if(it2) delete it2;
    }
    
    double nllval = nll->getVal();
//     double nLLatMLE = 0.;//m_fitResult->minNll();
//     double nlloffset = nll->getVal() - nLLatMLE;
    
    if(m_debug){
        std::cout << "   -> Initial value of the NLL = " << nllval << std::endl;
//         std::cout << "   -> Initial value of offset  = " << nlloffset << std::endl;
//         std::cout << "   -> Initial NLL - offset     = " << nllval-nlloffset << std::endl;
        constrainedParams->Print("v");
    }
    
    //
    // return here if specified not to perform the fit
    if(noFit) return nllval;
    
    // 
    // Desperate attempts
//     TVirtualFitter::SetPrecision(1e-10); // default 1e-6
//     TVirtualFitter::SetMaxIterations(100000);  // default 5000
    
    //
    //
    // Safe fit loop
    //
    //
    static int nrItr = 0;
    const int maxRetries = 3;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(m_minimType);
    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    if(TtHFitter::OPTION["FitStrategy"]!=0){
        strat = TtHFitter::OPTION["FitStrategy"];
        if(TtHFitter::OPTION["FitStrategy"]<0) strat = 0;
    }
    int save_strat = strat;
    RooMinimizer minim(*nll);
    minim.setStrategy(strat);
    minim.setPrintLevel(1);
    minim.setEps(1);

    // experimental - playing around fit minimisation precision
//     minim.setEps(100);
//     minim.setMaxIterations(500*200*10);
//     minim.setMaxFunctionCalls(500*200*10);
//     minim.setMaxIterations(500*200*10);
//     minim.setMaxFunctionCalls(500*200*10);
//     minim.setOffsetting(true);
    
    //
    // fast fit - e.g. for ranking
    if(fastFit){
        minim.setStrategy(0);  // to be the same as ttH comb
//         minim.setEps(0.1);  // to balance and not to have crazy results...
        minim.setPrintLevel(0);
//         minim.setMaxIterations(100*10);
//         minim.setMaxFunctionCalls(100*10);
//         minim.setMaxIterations(100*minim.getNPar());
//         minim.setMaxFunctionCalls(100*minim.getNPar());
    }
    
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

        if(m_useMinos){
            TIterator* it3 = model->GetNuisanceParameters()->createIterator();
            TIterator* it4 = model->GetParametersOfInterest()->createIterator();
            RooArgSet* SliceNPs = new RooArgSet( *(model->GetNuisanceParameters()) );
            SliceNPs->add(*(model->GetParametersOfInterest()));
            RooRealVar* var = NULL;
            RooRealVar* var2 = NULL;
            std::cout << "Size of variables for MINOS: " << m_varMinos.size() << std::endl;
            
            if (m_varMinos.at(0)!="all"){
                while( (var = (RooRealVar*) it3->Next()) ){
                    TString vname=var->GetName();
                    bool isthere=false;
                    for (unsigned int m=0;m<m_varMinos.size();++m){
                        //std::cout << "MINOS var: " << m_varMinos.at(m) << std::endl;
                        if(vname.Contains(m_varMinos.at(m))) {isthere=true; break;}
                        //cout << " --> NP: " << vname << endl;
                    }
                    if (!isthere) SliceNPs->remove(*var, true, true);
                }
                while( (var2 = (RooRealVar*) it4->Next()) ){
                    TString vname=var2->GetName();
                    bool isthere=false;
                    for (unsigned int m=0;m<m_varMinos.size();++m){
                        //std::cout << "MINOS var: " << m_varMinos.at(m) << std::endl;
                        if(vname.Contains(m_varMinos.at(m))) {isthere=true; break;}
                        //cout << " --> POI: " << vname << endl;
                    }
                    if (!isthere) SliceNPs->remove(*var2, true, true);
                }
                minim.minos(*SliceNPs);
            }
            else
                minim.minos();
            
            if(SliceNPs) delete SliceNPs;
            if(it3) delete it3;
            if(it4) delete it4;
        }//end useMinos
        
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
            return 0;
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
//     m_fitResult = r;
    m_fitResult = (RooFitResult*)r->Clone();
    delete r;

    //
    // clean stuff
//     if(constrainedParams) delete constrainedParams;
    
    nllval = 0;
//     nLLatMLE = 0;
//     nlloffset = 0;
    if(nll) nllval = nll->getVal();
//     if(m_fitResult) nLLatMLE = m_fitResult->minNll();
//     if(nll) nlloffset = nll->getVal() - nLLatMLE;
    
//     RooArgList poiList; 
//     poiList.addClone(fNullParams); // make a clone list 
//     Double_t deltaNLL = std::max( nLLatCondMLE-nLLatMLE, 0.);
//     RemoveConstantParameters(poiList);
//     int ndf = poiList.getSize();
//     Double_t pvalue = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndf);
    
    if(m_debug){
//         RemoveConstantParameters(poiList);
//         int ndof = poiList.getSize();
        int ndof = 1;
        std::cout << "   -> Final value of the NLL = " << nllval << std::endl;
//         std::cout << "   -> Final value of offset  = " << nlloffset << std::endl;
//         std::cout << "   -> Final NLL - offset     = " << nllval-nlloffset << std::endl;
//         std::cout << "   -> Goodness of fit (NLL/ndof) = " << nllval/ndof << std::endl;
    }
    if(nll) delete nll;
//     delete poi;  // creates a crash
//     poi->~RooRealVar();  // creates a crash
//     delete var;
//     delete nuis;
//     nuis->~RooArgSet();
//     if(glbObs) delete glbObs;
//     glbObs->~RooArgSet();
    return nllval;
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
    if(param) delete param;
    
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
    if(param) delete param;
    return result;
}
