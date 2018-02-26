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
#include "TtHFitter/StatusLogbook.h"

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
m_randomNP(0.1),
m_randSeed(-999)
{
    m_constNP.clear();
    m_constNPvalue.clear();
    m_systSubCategoryMap.clear();
    m_subCategories.clear();
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
void FittingTool::SetSubCategories() {
    WriteDebugStatus("FittingTool::SetSubCategories", "finding unique SubCategories");
    // loop over m_systSubCategoryMap to find all unique SubCategories, save in m_subCategories set
    for(std::map<std::string, std::string>::iterator it = m_systSubCategoryMap.begin(); it != m_systSubCategoryMap.end(); ++it) {
        //std::cout << it->first << " " << it->second << endl;
        m_subCategories.insert(it->second);
    }
}

//________________________________________________________________________
//
float FittingTool::FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit, bool noFit ) {

    if (TtHFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    WriteDebugStatus("FittingTool::FitPDF", "-> Entering in FitPDF function");

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
        if (TtHFitter::DEBUGLEVEL < 2) std::cout.clear();
        WriteErrorStatus("FittingTool::FitPDF", "Cannot find the parameter of interest !");
        return 0;
    }

    poi -> setConstant(m_constPOI);
    //poi -> setRange(m_RangePOI_down,m_RangePOI_up); // Commented by Loic to avoid overwriting user's setting in config file

    poi -> setVal(m_valPOI);
    if (TtHFitter::DEBUGLEVEL < 2) std::cout.clear();
    WriteInfoStatus("FittingTool::FitPDF", "Setting starting mu = " + std::to_string(m_valPOI));
    if (TtHFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    // randomize the POI
    if(!m_constPOI && m_randomize){
        poi->setVal( m_valPOI + m_randomNP*(gRandom->Uniform(2)-1.) );
    }

    WriteDebugStatus("FittingTool::FitPDF", "   -> Constant POI : " + std::to_string(poi->isConstant()));
    WriteDebugStatus("FittingTool::FitPDF", "   -> Value of POI : " + std::to_string(poi->getVal()));

    RooRealVar* var = NULL;
    RooArgSet* nuis = (RooArgSet*) model->GetNuisanceParameters();
    if (TtHFitter::DEBUGLEVEL < 2) std::cout.clear();
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
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant(1);
                found = true;
            }
            else if(np.find("alpha_")!=string::npos && m_noSystematics){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
//                 var->setVal( 0 );
                found = true;
            }
            else if(m_noNormFactors){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
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
                    WriteInfoStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(m_constNPvalue[i_np]));
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
                    WriteInfoStatus("FittingTool::FitPDF", " ---> Setting " + m_initialNP[i_np] + " to "  +std::to_string(m_initialNPvalue[i_np]));
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
    double nLLatMLE = 0.;//m_fitResult->minNll();
    double nlloffset = nll->getVal() - nLLatMLE;

    WriteDebugStatus("FittingTool::FitPDF","   -> Initial value of the NLL = " +std::to_string(nllval));
    if(m_debug){
//         std::cout << "   -> Initial value of offset  = " << nlloffset << std::endl;
//         std::cout << "   -> Initial NLL - offset     = " << nllval-nlloffset << std::endl;
        if(TtHFitter::DEBUGLEVEL > 2) constrainedParams->Print("v");
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

    if (TtHFitter::DEBUGLEVEL < 2) std::cout.clear();
    while (nrItr<maxRetries && status!=0 && status!=1){

        WriteInfoStatus("FittingTool::FitPDF", "");
        WriteInfoStatus("FittingTool::FitPDF", "");
        WriteInfoStatus("FittingTool::FitPDF", "");
        WriteInfoStatus("FittingTool::FitPDF", "Fit try no." + std::to_string(nrItr+1));
        WriteInfoStatus("FittingTool::FitPDF", "======================");
        WriteInfoStatus("FittingTool::FitPDF", "");

        ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        m_hessStatus= minim.hesse();
        r = minim.save();
        m_edm = r->edm();

        //up the strategy
        bool FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if (FitIsNotGood && strat < 2){
            WriteInfoStatus("FittingTool::FitPDF", "");
            WriteInfoStatus("FittingTool::FitPDF", "   *******************************");
            WriteInfoStatus("FittingTool::FitPDF", "   * Increasing Minuit strategy (was " + std::to_string(strat) + ")");
            strat++;
            WriteInfoStatus("FittingTool::FitPDF", "   * Fit failed with : ");
            WriteInfoStatus("FittingTool::FitPDF", "      - minuit status " + std::to_string(status));
            WriteInfoStatus("FittingTool::FitPDF", "      - hess status " + std::to_string(m_hessStatus));
            WriteInfoStatus("FittingTool::FitPDF", "      - Edm = " + std::to_string(m_edm));
            WriteInfoStatus("FittingTool::FitPDF", "   * Retrying with strategy " + std::to_string(strat));
            WriteInfoStatus("FittingTool::FitPDF", "   ********************************");
            WriteInfoStatus("FittingTool::FitPDF", "");
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            m_hessStatus= minim.hesse();
            r = minim.save();
            m_edm = r->edm();
        }

        FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if (FitIsNotGood && strat < 2){
            WriteInfoStatus("FittingTool::FitPDF", "");
            WriteInfoStatus("FittingTool::FitPDF", "   ********************************");
            WriteInfoStatus("FittingTool::FitPDF", "   * Increasing Minuit strategy (was " + std::to_string(strat) + ")");
            strat++;
            WriteInfoStatus("FittingTool::FitPDF", "   * Fit failed with : ");
            WriteInfoStatus("FittingTool::FitPDF", "      - minuit status " + std::to_string(status));
            WriteInfoStatus("FittingTool::FitPDF", "      - hess status " + std::to_string(m_hessStatus));
            WriteInfoStatus("FittingTool::FitPDF", "      - Edm = " + std::to_string(m_edm));
            WriteInfoStatus("FittingTool::FitPDF", "   * Retrying with strategy " + std::to_string(strat));
            WriteInfoStatus("FittingTool::FitPDF", "   ********************************");
            WriteInfoStatus("FittingTool::FitPDF", "");
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            r = minim.save();
            m_edm = r->edm();
        }

        FitIsNotGood = ((status!=0 && status!=1) || (m_hessStatus!=0 && m_hessStatus!=1) || m_edm>1.0);
        if (FitIsNotGood && strat < 2){
            WriteInfoStatus("FittingTool::FitPDF", "");
            WriteInfoStatus("FittingTool::FitPDF", "   ********************************");
            WriteInfoStatus("FittingTool::FitPDF", "   * Increasing Minuit strategy (was " + std::to_string(strat) + ")");
            strat++;
            WriteInfoStatus("FittingTool::FitPDF", "   * Fit failed with : ");
            WriteInfoStatus("FittingTool::FitPDF", "      - minuit status " + std::to_string(status));
            WriteInfoStatus("FittingTool::FitPDF", "      - hess status " + std::to_string(m_hessStatus));
            WriteInfoStatus("FittingTool::FitPDF", "      - Edm = " + std::to_string(m_edm));
            WriteInfoStatus("FittingTool::FitPDF", "   * Retrying with strategy " + std::to_string(strat));
            WriteInfoStatus("FittingTool::FitPDF", "   ********************************");
            WriteInfoStatus("FittingTool::FitPDF", "");
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
            WriteDebugStatus("FittingTool::FitPDF", "Size of variables for MINOS: " + std::to_string(m_varMinos.size()));

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
            WriteWarningStatus("FittingTool::FitPDF", "");
            WriteWarningStatus("FittingTool::FitPDF", "");
            WriteWarningStatus("FittingTool::FitPDF", "");
            WriteWarningStatus("FittingTool::FitPDF", "");
            WriteWarningStatus("FittingTool::FitPDF", "***********************************************************");
            WriteWarningStatus("FittingTool::FitPDF", "Fit failure unresolved with status " + std::to_string(status));
            WriteWarningStatus("FittingTool::FitPDF", "   Please investigate your workspace");
            WriteWarningStatus("FittingTool::FitPDF", "   Find a wall : you will need it to crash your head on it");
            WriteWarningStatus("FittingTool::FitPDF", "***********************************************************");
            WriteWarningStatus("FittingTool::FitPDF", "");
            WriteWarningStatus("FittingTool::FitPDF", "");
            WriteWarningStatus("FittingTool::FitPDF", "");
            m_minuitStatus = status;
            m_fitResult = 0;
            return 0;
        }

    }

    r = minim.save();
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
    WriteInfoStatus("FittingTool::FitPDF", "         FIT FINALIZED SUCCESSFULLY : ");
    WriteInfoStatus("FittingTool::FitPDF", "            - minuit status " + std::to_string(status));
    WriteInfoStatus("FittingTool::FitPDF", "            - hess status " + std::to_string(m_hessStatus));
    WriteInfoStatus("FittingTool::FitPDF", "            - Edm = " + std::to_string(m_edm));
    WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
    if (TtHFitter::DEBUGLEVEL > 2) sw.Print();
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");

    m_minuitStatus = status;
//     m_fitResult = r;
    m_fitResult = (RooFitResult*)r->Clone();
    delete r;

    //
    // clean stuff
//     if(constrainedParams) delete constrainedParams;

    nllval = 0;
    nLLatMLE = 0;
    nlloffset = 0;
    if(nll) nllval = nll->getVal();
    if(m_fitResult) nLLatMLE = m_fitResult->minNll();
    if(nll) nlloffset = nll->getVal() - nLLatMLE;

//     RooArgList poiList;
//     poiList.addClone(fNullParams); // make a clone list
//     Double_t deltaNLL = std::max( nLLatCondMLE-nLLatMLE, 0.);
//     RemoveConstantParameters(poiList);
//     int ndf = poiList.getSize();
//     Double_t pvalue = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndf);

    if(m_debug){
//         RemoveConstantParameters(poiList);
        double redNLL = nllval - 1000000.0;
        std::stringstream redNLL_ss;
        redNLL_ss << std::fixed << std::setprecision(20) << redNLL;

        std::cout << std::fixed << std::setprecision(20);

        WriteInfoStatus("FittingTool::FitPDF", "   -> Reduced Final value of the NLL = " + redNLL_ss.str());
        WriteInfoStatus("FittingTool::FitPDF", "   -> Final value of the NLL = " + std::to_string(nllval));
        WriteInfoStatus("FittingTool::FitPDF", "   -> Final value of offset = " + std::to_string(nlloffset));
        WriteInfoStatus("FittingTool::FitPDF", "   -> Final NLL - offset = " + std::to_string(nllval-nlloffset));
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
        WriteErrorStatus("FittingTool::ExportFitResultInTextFile", "The FitResultObject seems not to be defined.");
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
        WriteErrorStatus("FittingTool::ExportFitResultInMap", "The FitResultObject seems not to be defined.");
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

//____________________________________________________________________________________
//
int FittingTool::GetGroupedImpact( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, RooWorkspace* ws ){

    // obtain constrainedParams and poi just as in FitPDF()
    RooArgSet* constrainedParams = fitpdf->getParameters(*fitdata);
    RooStats::RemoveConstantParameters(constrainedParams);
    RooFit::Constrain(*constrainedParams);

    RooRealVar * poi = (RooRealVar*) model->GetParametersOfInterest()->first();
    if(!poi){
        if (TtHFitter::DEBUGLEVEL < 2) std::cout.clear();
        WriteErrorStatus("FittingTool::FitPDF", "Cannot find the parameter of interest !");
        return -1;
    }
    poi -> setConstant(m_constPOI);
    poi -> setVal(m_valPOI);
    if(!m_constPOI && m_randomize){
        poi->setVal( m_valPOI + m_randomNP*(gRandom->Uniform(2)-1.) );
    }

    if (!poi->isConstant()) {
        // save snapshot of original workspace
        ws->saveSnapshot("snapshot_AfterFit_POI", *model->GetParametersOfInterest() );
        ws->saveSnapshot("snapshot_AfterFit_NP" , *(model->GetNuisanceParameters()) );
        ws->saveSnapshot("snapshot_AfterFit_GO" , *(model->GetGlobalObservables())  );
    }

    std::vector<std::string> associatedSysts;

    ScanSingleParamReversed(true, true , fitdata, fitpdf, constrainedParams, model, ws, "Stat_Norm", associatedSysts);
    ScanSingleParamReversed(true, false, fitdata, fitpdf, constrainedParams, model, ws, "Stat_Gamma", associatedSysts);
    ScanSingleParamReversed(true, true,  fitdata, fitpdf, constrainedParams, model, ws, "Stat_Stat", associatedSysts);

    // loop over unique SubCategories
    for (std::set<std::string>::iterator itCategories = m_subCategories.begin(); itCategories != m_subCategories.end(); ++itCategories){
        WriteInfoStatus("FittingTool::GetGroupedImpact","performing grouped systematics impact evaluation for: " + *itCategories);

        // find all associated systematics
        associatedSysts.clear();
        for(std::map<std::string, std::string>::iterator itSysts = m_systSubCategoryMap.begin(); itSysts != m_systSubCategoryMap.end(); ++itSysts) {
            if (itSysts->second == *itCategories) {
                WriteDebugStatus("FittingTool::GetGroupedImpact","  associated syst: " + itSysts->first);
                associatedSysts.push_back("alpha_" + itSysts->first); // note the required "alpha_" prefix
            }
        }

        if (!poi->isConstant()) {
            ScanSingleParamReversed(true, true , fitdata, fitpdf, constrainedParams, model, ws, *itCategories, associatedSysts );
        }
    }

    if (!poi->isConstant()) {
        ws->loadSnapshot("snapshot_AfterFit_GO");
        ws->loadSnapshot("snapshot_AfterFit_POI");
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }

    WriteInfoStatus("FittingTool::FitPDF","-----------------------------------------------------");
    WriteInfoStatus("FittingTool::GetGroupedImpact", "final results:");

    cout << setprecision(2) << fixed << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_Full");
    cout << "POI          is: " << poi->getVal() << " +/- " << poi->getError() << "  :  high: " << poi->getErrorHi() << "  low: " << poi->getErrorLo() << endl;
    float StatUp2=(poi->getErrorHi()*poi->getErrorHi());
    float StatLo2=(poi->getErrorLo()*poi->getErrorLo());
    float Stat2  =(poi->getError()*poi->getError());
    ws->loadSnapshot("snapshot_AfterFit_POI_Stat");
    cout << "POI StatOnly is: " << poi->getVal() << " +/- " << poi->getError() << "  :  high: " << poi->getErrorHi() << "  low: " << poi->getErrorLo() << endl;
    // float StatUp2=(poi->getErrorHi()*poi->getErrorHi());
    // float StatLo2=(poi->getErrorLo()*poi->getErrorLo());
    // float Stat2  =(poi->getError()*poi->getError());
    ws->loadSnapshot("snapshot_AfterFit_POI_Gamma");
    cout << "POI Stat+Gam is: " << poi->getVal() << " +/- " << poi->getError() << "  :  high: " << poi->getErrorHi() << "  low: " << poi->getErrorLo() << endl;
    //cout << "   --> MC stat.  : +/- " << sqrt( - pow(poi->getError(),2) + Stat2 ) << endl;
    cout << "   --> MC stat.:\t +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    // ws->loadSnapshot("snapshot_AfterFit_POI_CRstat");
    // cout << "   --> CRstat.            : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    //
    std::cout<<std::endl;

    for (std::set<std::string>::iterator itCategories = m_subCategories.begin(); itCategories != m_subCategories.end(); ++itCategories){
        ws->loadSnapshot(("snapshot_AfterFit_POI_" + *itCategories).c_str());
        std::cout << "POI          is: " << poi->getVal() << " +/- " << poi->getError() << "  :  high: " << poi->getErrorHi() << "  low: " << poi->getErrorLo() << std::endl;
        std::cout << "   --> " + *itCategories + " impact:\t +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << std::endl;
        std::cout<<std::endl;
    }

    /*
    ws->loadSnapshot("snapshot_AfterFit_POI_ttbNorm");
    cout << "   --> ttb normalization  : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_ttcNorm");
    cout << "   --> ttc normalization  : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_Theo");
    cout << "   --> ttH modelling      : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_FTAG");
    cout << "   --> Flavour tagging    : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_JE");
    cout << "   --> JET and JER        : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_ttb");
    cout << "   --> ttb modelling      : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_ttbGen");
    cout << "   --> ttbGen modelling      : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_ttc");
    cout << "   --> ttc modelling      : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_ttlight");
    cout << "   --> ttlight modelling  : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_oth");
    cout << "   --> Other backgrounds  : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_PRWjvt");
    cout << "   --> PRW and JVT        : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_lumi");
    cout << "   --> Luminosity         : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI_lepton");
    cout << "   --> Leptons            : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    // ws->loadSnapshot("snapshot_AfterFit_POI_MET");
    // cout << "   --> MET                : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    */

    // ws->loadSnapshot("snapshot_AfterFit_POI_Rest");
    // cout << "   --> Rest               : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    ws->loadSnapshot("snapshot_AfterFit_POI");
//    cout << "   --> Full sys. : +" << sqrt( - pow(poi->getErrorHi(),2) + StatUp2 ) << " ,  -" << sqrt( - pow(poi->getErrorLo(),2) + StatLo2 ) << endl;
    //cout << " MC stat. abs: " << sqrt( - pow(GammaUp,2) - - pow(poi->getErrorHi(),2) ) << " ,  - " << sqrt( - pow(GammaDo,2) - - pow(poi->getErrorLo(),2) ) << endl;
    //cout << " MC stat. rel: " << sqrt( - pow(GammaUp,2) - - pow(poi->getErrorHi(),2) )/poi->getVal()*100 << " ,  - " << sqrt( - pow(GammaDo,2) - - pow(poi->getErrorLo(),2) )/poi->getVal()*100 << endl;
    //ws->loadSnapshot("snapshot_AfterFit_POI");
//    cout << endl << endl;
    // ws->loadSnapshot("snapshot_AfterFit_POI_CRstat");
    // cout << "POI wCRstat is: " << poi->getVal() << " +/- " << poi->getError() << "  :  high: " << poi->getErrorHi() << "  low: " << poi->getErrorLo() << endl;
    // StatUp2=(poi->getErrorHi()*poi->getErrorHi());
    // StatLo2=(poi->getErrorLo()*poi->getErrorLo());
    // Stat2  =(poi->getError()*poi->getError());
    ws->loadSnapshot("snapshot_AfterFit_POI_Stat");
    cout << "POI StatOnly is: " << poi->getVal() << " +/- " << poi->getError() << "  :  high: " << poi->getErrorHi() << "  low: " << poi->getErrorLo() << endl;
    StatUp2=(poi->getErrorHi()*poi->getErrorHi());
    StatLo2=(poi->getErrorLo()*poi->getErrorLo());
    Stat2  =(poi->getError()*poi->getError());

    ws->loadSnapshot("snapshot_AfterFit_POI");
    cout << "   --> Full sys. : +" << sqrt(  pow(poi->getErrorHi(),2) - StatUp2 ) << " ,  -" << sqrt(  pow(poi->getErrorLo(),2) - StatLo2 ) << endl;


    float step=0.25;
    int npoints=(int)((4.5+0.5)/step);
    for (int i=0; i<npoints+1; i++) {
    //GetLH( -0.5+((float)i)*step, LHmin, fitdata, fitpdf, constrainedParams);
    }
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    WriteInfoStatus("FittingTool::FitPDF","finishing grouped ranking section");
    WriteInfoStatus("FittingTool::FitPDF","-----------------------------------------------------");

    return 0;
}


//____________________________________________________________________________________
//
float FittingTool::ScanSingleParamReversed(bool doStat, bool excludeGammas, RooAbsData*& fitdata, RooAbsPdf*& fitpdf, RooArgSet*& constrainedParams, RooStats::ModelConfig* mc,
                                           RooWorkspace* ws, std::string category, std::vector<std::string> affectedSysts, bool calib ) {

    string nameV = "";

    // (VD): use this to fix nuisance parameter before the fit
    const RooArgSet* glbObs = mc->GetGlobalObservables();
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    WriteInfoStatus("FittingTool::ScanSingleParamReversed", "-----------------------------------------------------");
    WriteInfoStatus("FittingTool::ScanSingleParamReversed", "           breakdown for " + category);
    WriteInfoStatus("FittingTool::ScanSingleParamReversed", "-----------------------------------------------------");

    cout << "Scan : " << category << endl;

    TIterator* it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar* var2 = NULL;
    if (!calib) {
        while( (var2 = (RooRealVar*) it->Next()) ){
            string varname = (string) var2->GetName();
            if (doStat) {
                // remove everything but gammas
                if (!excludeGammas) {
                    if (varname.find("gamma")!=string::npos) var2->setConstant(1);
                    continue;
                }

                // DISABLE EVERYTHING!!!
                var2->setConstant(0);


                // new behavior
                if (std::find(affectedSysts.begin(), affectedSysts.end(), varname) != affectedSysts.end()) {
                  WriteDebugStatus("FittingTool::ScanSingleParamReversed", "setting " + varname + " constant");
                  var2->setConstant(1);
                }

                if (category=="Stat_Stat") {
                    // enable the CR Stat.
                    var2->setConstant(1);
                    // ws->var("ttb_norm")->setVal(1);
                    // ws->var("ttc_norm")->setVal(1);
                }

                // EXAMPLE for NFs:
                //    ws->var("ttb_norm")->setConstant(1);

                // EXAMPLE "default" behavior
                //if (nameV=="Stat_ttlight") {
                //    // enable the CR Stat.
                //    if ( varname.find("alpha_ttlight_")!=string::npos || varname.find("alpha_tt_XS")!=string::npos ) var2->setConstant(1);
                //}

            }
            else {
                if (varname.find(nameV)!=string::npos) {
                    var2->setConstant(1);
                }
            }
        }
    }


    constrainedParams->Print("v");
    // repeat the fit here ....
    RooAbsReal* nll = fitpdf->createNLL(*fitdata, RooFit::Constrain(*constrainedParams), RooFit::GlobalObservables(*glbObs), RooFit::Offset(1), NumCPU(4, RooFit::Hybrid) );
    RooMinimizer minim2(*nll);
    minim2.setStrategy(1);
    // minim2.setPrintLevel(-1);
    minim2.setPrintLevel(1);
    minim2.setEps(1);
    int status = minim2.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    ////////////////RooFitResult * r=minim2.save();
    RooRealVar * thePOI = dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());

    bool HessStatus= minim2.hesse();
    RooArgSet minosSet(*thePOI);
    ///////////////////////////////////////if ( ws->var("normNF_ttb")!=0 ) minosSet.add( *(ws->var("normNF_ttb") ) );
    ///////////////////////////////////////if ( ws->var("normNF_ttc")!=0 ) minosSet.add( *(ws->var("normNF_ttc") ) );
    //if ( ws->var("mu_XS_ttH_tthbb")!=0 )  minosSet.add( *(ws->var("mu_XS_ttH_tthbb") ) );
    ///////////////////////////////////////if ( ws->var("mu_XS_ttH_tthlep")!=0 ) minosSet.add( *(ws->var("mu_XS_ttH_tthlep") ) );
    minim2.minos(minosSet);

    if (status!=0) cout << "I WAS UNABLE TO PERFORM THE FIT CORRECTLY ... HessStatus: " << HessStatus << endl;

    //cout << "VALERIO SAYS: after fit for configuration: " << nameV << " , " << doStat << " , " << excludeGammas << "   --> with MU= " << thePOI->getVal() <<  " +/- " << thePOI->getError() << endl;
    float newPOIerr =thePOI->getError();
    float newPOIerrU=thePOI->getErrorHi();
    float newPOIerrD=thePOI->getErrorLo();
    float newPOIVal=thePOI->getVal();
    if (category=="Stat_Norm")  ws->saveSnapshot("snapshot_AfterFit_POI_Full"  , *mc->GetParametersOfInterest() );
    if (category=="Stat_Gamma") ws->saveSnapshot("snapshot_AfterFit_POI_Gamma" , *mc->GetParametersOfInterest() );

    if (category=="Stat_Stat")   ws->saveSnapshot("snapshot_AfterFit_POI_Stat"  , *mc->GetParametersOfInterest() );

    if (nameV=="Stat_ttbNorm")   ws->saveSnapshot("snapshot_AfterFit_POI_ttbNorm"  , *mc->GetParametersOfInterest() );
    if (nameV=="Stat_ttlight")   ws->saveSnapshot("snapshot_AfterFit_POI_ttlight"  , *mc->GetParametersOfInterest() );


    // generic implementation
    if (nameV=="") ws->saveSnapshot(("snapshot_AfterFit_POI_" + category).c_str(), *mc->GetParametersOfInterest() );

    if (!doStat) ws->loadSnapshot("snapshot_AfterFit_POI_CALIB");
    else         ws->loadSnapshot("snapshot_AfterFit_POI");
    float oldPOIerr =thePOI->getError();
    float oldPOIerrU=thePOI->getErrorHi();
    float oldPOIerrD=thePOI->getErrorLo();
    float oldPOIVal=thePOI->getVal();
    //cout << "newPOIval= " << newPOIVal << "   ...oldPOIval: " << oldPOIVal << endl;
    // recalibrate ... if necessary
    /*
      newPOIerr += (oldPOIVal-newPOIVal);
      newPOIerrU+= (oldPOIVal-newPOIVal);
      newPOIerrD+= (oldPOIVal-newPOIVal);
    */

    /*
      cout << " How things changes: " << oldPOIerr << " --> " << newPOIerr
      << " ||||  UP: " << oldPOIerrU << " --> " << newPOIerrU
      << " ||||  DO: " << oldPOIerrD << " --> " << newPOIerrD << endl;
    */
    if (doStat && excludeGammas) {
        //cout << " Stat. rel Err: " << newPOIerrU/thePOI->getVal()*100 << "   ,   " << newPOIerrD/thePOI->getVal()*100 << endl;
    }
    else if (doStat && !excludeGammas) {
        //cout << " Stat+Gam Err : " << newPOIerrU/thePOI->getVal()*100 << "   ,   " << newPOIerrD/thePOI->getVal()*100 << endl;
    }
    else {
        if ( (fabs(newPOIerrU)>fabs(oldPOIerrU)) || (fabs(newPOIerrD)>fabs(oldPOIerrD)) ) {
            cout << " PROBLEM for sys: " << category << " .... please check" << endl;
            cout << "      Error: " << oldPOIerr << " --> " << newPOIerr
                 << " ||||  UP: " << oldPOIerrU << " --> " << newPOIerrU
                 << " ||||  DO: " << oldPOIerrD << " --> " << newPOIerrD << endl;
        }
        //cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << nameV << "    :   "
        //<< sqrt(oldPOIerrU*oldPOIerrU-newPOIerrU*newPOIerrU)/thePOI->getVal()*100 << "   ,   "
        //  << sqrt(oldPOIerrD*oldPOIerrD-newPOIerrD*newPOIerrD)/thePOI->getVal()*100 << endl << endl;
        //// then please do the average
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //if (!calib) { //this might be tricky when we have hidden NP that needs to be kept fixed ... but the snapshot should protect this
      //it = constrainedParams->createIterator();
      //var2 = NULL;
      //while( (var2 = (RooRealVar*) it->Next()) ){
      //string varname = (string) var2->GetName();
      //var2->setConstant(0);
      //
    //}
    //cout << "VALERIO: " << nameV << " after everything POI error: " << thePOI->getError() << endl;

    return newPOIerr;
}
