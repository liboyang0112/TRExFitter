// Class include
#include "TRExFitter/FittingTool.h"

//Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"

//ROOR includes
#include "TCanvas.h"
#include "TH2.h"
#include "TRandom3.h"

//Roostats includes
#include "RooMinimizer.h"
#include "Math/MinimizerOptions.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"

//Roofit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooArgSet.h"

//c++ includes
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

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
    m_constPOI(false),
    m_fitResult(nullptr),
    m_debug(1),
    m_noGammas(false),
    m_noSystematics(false),
    m_noNormFactors(false),
    m_noShapeFactors(false),
    m_RangePOI_up(100.),
    m_RangePOI_down(-10.),
    m_randomize(false),
    m_randomNP(0.1),
    m_randSeed(-999),
    m_externalConstraints(nullptr) {
}

//________________________________________________________________________
//
FittingTool::~FittingTool() {
}

//________________________________________________________________________
//
void FittingTool::SetSubCategories() {
    WriteDebugStatus("FittingTool::SetSubCategories", "finding unique SubCategories");
    // loop over m_subCategoryMap to find all unique SubCategories, save in m_subCategories set
    for(std::map<std::string, std::string>::iterator it = m_subCategoryMap.begin(); it != m_subCategoryMap.end(); ++it) {
        m_subCategories.insert(it->second);
    }
}

//________________________________________________________________________
//
double FittingTool::FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit, bool noFit, bool saturatedModel ) {

    if (m_debug < 1) std::cout.setstate(std::ios_base::failbit);
    WriteDebugStatus("FittingTool::FitPDF", "-> Entering in FitPDF function");

    //
    // Printing the whole model for information
    //
    if(m_debug >= 2) model->Print();

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
    std::unique_ptr<RooAbsReal> nll(fitpdf->createNLL(*fitdata,
                                    RooFit::Constrain(*constrainedParams),
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(1),
                                    RooFit::NumCPU(TRExFitter::NCPU,RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE),
                                    RooFit::ExternalConstraints(*m_externalConstraints)
                                    ));

    //
    // Needed for Ranking plot, but also to set random initial values for the NPs
    //
    if(m_randSeed == -999){
        gRandom->SetSeed(time(nullptr));
    }
    else{
        gRandom->SetSeed(m_randSeed);
    }

    //
    // Getting the POI
    //
    RooRealVar * poi = static_cast<RooRealVar*>(model->GetParametersOfInterest()->first());
    if(!poi){
        if (m_debug < 1) std::cout.clear();
        WriteErrorStatus("FittingTool::FitPDF", "Cannot find the parameter of interest !");
        return 0;
    }
    poi -> setConstant(m_constPOI || saturatedModel);
    poi -> setVal(m_valPOI);
    if (m_debug < 1) std::cout.clear();
    WriteDebugStatus("FittingTool::FitPDF", "Setting starting mu = " + std::to_string(m_valPOI));
    if (m_debug < 1) std::cout.setstate(std::ios_base::failbit);

    // randomize the POI
    if(!m_constPOI && m_randomize){
        poi->setVal( m_valPOI + m_randomNP*(gRandom->Uniform(2)-1.) );
    }

    WriteDebugStatus("FittingTool::FitPDF", "   -> Constant POI : " + std::to_string(poi->isConstant()));
    WriteDebugStatus("FittingTool::FitPDF", "   -> Value of POI : " + std::to_string(poi->getVal()));

    RooRealVar* var(nullptr);
    const RooArgSet* nuis = static_cast<const RooArgSet*>(model->GetNuisanceParameters());
    if(nuis){
        std::unique_ptr<TIterator> it2(nuis->createIterator());
        while( (var = static_cast<RooRealVar*>(it2->Next())) ){
            const std::string np = var->GetName();
            bool found = false;
            //
            // first check if all systs, norm and gammas should be set to constant
            if((np.find("gamma_stat")!=string::npos || np.find("gamma_shape_stat")!=string::npos) && m_noGammas){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                var->setVal( 1 );
                found = true;
            }
            else if((np.find("alpha_")!=string::npos || (np.find("gamma_shape")!=string::npos && np.find("gamma_shape_stat")==string::npos)) && m_noSystematics){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                var->setVal( 0 );
                found = true;
            }
            else if(np.find("alpha_")==string::npos && np.find("gamma_")==string::npos && (m_noNormFactors || saturatedModel)){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                found = true;
            }
            if(found) continue;
            //
            // set to constant the saturatedModel shape factor parameters if saturatedModel is left to false
            if(np.find("saturated_model_sf_")!=std::string::npos && !saturatedModel){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                found = true;
            }
            if(found) continue;
            //
            // loop on the NP specified to have custom starting value
            for( unsigned int i_np = 0; i_np<m_initialNP.size(); i_np++ ){
                if( np == ("alpha_"+m_initialNP[i_np]) || np == m_initialNP[i_np] ){
                    var->setVal(m_initialNPvalue[i_np]);
                    WriteInfoStatus("FittingTool::FitPDF", " ---> Setting " + m_initialNP[i_np] + " to "  +std::to_string(m_initialNPvalue[i_np]));
                    found = true;
                    break;
                }
            }
            //
            // loop on the NP specified to be constant - This should be after setting to initial value
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
    }

    double nllval = nll->getVal();
    double nLLatMLE = 0.;//m_fitResult->minNll();

    WriteDebugStatus("FittingTool::FitPDF","   -> Initial value of the NLL = " +std::to_string(nllval));
    if(m_debug >= 2) constrainedParams->Print("v");

    //
    // return here if specified not to perform the fit
    if(noFit) return nllval;

    //
    // Safe fit loop
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(m_minimType);
    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    if(TRExFitter::OPTION["FitStrategy"]!=0){
        strat = TRExFitter::OPTION["FitStrategy"];
        if(TRExFitter::OPTION["FitStrategy"]<0) strat = 0;
        if(TRExFitter::OPTION["FitStrategy"]>3){
            WriteWarningStatus("FittingTool::FitPDF","Set strategy > 3. Setting to default (1)");
            strat = 1;
        }
    }
    RooMinimizer minim(*nll);
    minim.setStrategy(strat);
    minim.setPrintLevel(1);
    minim.setEps(1);
    // it doesn't make sense to try more than 2 additional strategies
    const int maxRetries = 3 - strat;

    // experimental - playing around fit minimisation precision
//     minim.setEps(100);
//     minim.setMaxIterations(500*200*10);
//     minim.setMaxFunctionCalls(500*200*10);
//     minim.setMaxIterations(500*200*10);
//     minim.setMaxFunctionCalls(500*200*10);
//     minim.setOffsetting(true);

    // fast fit - e.g. for ranking
    if(fastFit){
        minim.setStrategy(0);  // to be the same as ttH comb
        minim.setPrintLevel(0);
//         minim.setEps(0.1);  // to balance and not to have crazy results...
//         minim.setMaxIterations(100*10);
//         minim.setMaxFunctionCalls(100*10);
//         minim.setMaxIterations(100*minim.getNPar());
//         minim.setMaxFunctionCalls(100*minim.getNPar());
    }

    TStopwatch sw;
    sw.Start();

    int status=-99;
    m_hessStatus=-99;
    m_edm = -99;

    // always run one fit
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "Fit try no." + std::to_string(1));
    WriteInfoStatus("FittingTool::FitPDF", "======================");
    WriteInfoStatus("FittingTool::FitPDF", "");

    ROOT::Math::MinimizerOptions::SetDefaultStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    m_hessStatus= minim.hesse();
    std::unique_ptr<RooFitResult> r(minim.save());
    m_edm = r->edm();

    // check if the fit converged
    bool FitIsNotGood = !(status==0 && m_hessStatus==0 && m_edm<0.001);

    // if not, loop maximum 2 times with increased strategy
    int nrItr = 0;
    while (nrItr<maxRetries && FitIsNotGood){
        WriteWarningStatus("FittingTool::FitPDF", "");
        WriteWarningStatus("FittingTool::FitPDF", "   *******************************");
        WriteWarningStatus("FittingTool::FitPDF", "   * Increasing Minuit strategy (was " + std::to_string(strat) + ")");
        strat++;
        WriteWarningStatus("FittingTool::FitPDF", "   * Fit failed with : ");
        WriteWarningStatus("FittingTool::FitPDF", "      - minuit status " + std::to_string(status));
        WriteWarningStatus("FittingTool::FitPDF", "      - hess status " + std::to_string(m_hessStatus));
        WriteWarningStatus("FittingTool::FitPDF", "      - Edm = " + std::to_string(m_edm));
        WriteWarningStatus("FittingTool::FitPDF", "   * Retrying with strategy " + std::to_string(strat));
        WriteWarningStatus("FittingTool::FitPDF", "   ********************************");
        WriteWarningStatus("FittingTool::FitPDF", "");
        minim.setStrategy(strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        m_hessStatus= minim.hesse();
        r.reset(minim.save());
        m_edm = r->edm();

        FitIsNotGood = !(status==0 && m_hessStatus==0 && m_edm<0.001);
        nrItr++;
    }

    // if the fit is not good even after retries print an error message
    if (FitIsNotGood) {
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "***********************************************************");
        WriteErrorStatus("FittingTool::FitPDF", "Fit failure unresolved with status " + std::to_string(status));
        WriteErrorStatus("FittingTool::FitPDF", "   Please investigate your workspace");
        WriteErrorStatus("FittingTool::FitPDF", "   Find a wall : you will need it to crash your head on it");
        WriteErrorStatus("FittingTool::FitPDF", "***********************************************************");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        m_minuitStatus = status;
        m_fitResult = nullptr;

        return 0;
    }

    if(m_useMinos){
        std::unique_ptr<TIterator> it3(model->GetNuisanceParameters()->createIterator());
        std::unique_ptr<TIterator> it4(model->GetParametersOfInterest()->createIterator());
        std::unique_ptr<RooArgSet> SliceNPs(new RooArgSet( *(model->GetNuisanceParameters()) ));
        SliceNPs->add(*(model->GetParametersOfInterest()));
        RooRealVar* var_tmp = nullptr;
        RooRealVar* var2 = nullptr;
        WriteDebugStatus("FittingTool::FitPDF", "Size of variables for MINOS: " + std::to_string(m_varMinos.size()));

        if (m_varMinos.at(0)!="all"){
            while( (var_tmp = static_cast<RooRealVar*>(it3->Next())) ){
                TString vname=var_tmp->GetName();
                bool isthere=false;
                for (unsigned int m=0;m<m_varMinos.size();++m){
                    if(vname.Contains(m_varMinos.at(m))) {isthere=true; break;}
                }
                if (!isthere) SliceNPs->remove(*var_tmp, true, true);
            }
            while( (var2 = static_cast<RooRealVar*>(it4->Next())) ){
                TString vname=var2->GetName();
                bool isthere=false;
                for (unsigned int m=0;m<m_varMinos.size();++m){
                    if(vname.Contains(m_varMinos.at(m))) {isthere=true; break;}
                }
                if (!isthere) SliceNPs->remove(*var2, true, true);
            }
            minim.minos(*SliceNPs);
        }
        else {
            minim.minos();
        }

    }//end useMinos

    r.reset(minim.save());
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
    WriteInfoStatus("FittingTool::FitPDF", "         FIT FINALIZED SUCCESSFULLY : ");
    WriteInfoStatus("FittingTool::FitPDF", "            - minuit status " + std::to_string(status));
    WriteInfoStatus("FittingTool::FitPDF", "            - hess status " + std::to_string(m_hessStatus));
    WriteInfoStatus("FittingTool::FitPDF", "            - Edm = " + std::to_string(m_edm));
    WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
    if (m_debug >= 2) sw.Print();
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");

    m_minuitStatus = status;
    if(r!=nullptr) m_fitResult = std::unique_ptr<RooFitResult>(static_cast<RooFitResult*>(r->Clone()));

    //
    // clean stuff

    nllval = 0;
    nLLatMLE = 0;
    double nlloffset = 0;
    nllval = nll->getVal();
    if(m_fitResult) nLLatMLE = m_fitResult->minNll();
    nlloffset = nll->getVal() - nLLatMLE;

    if(m_debug >= 1) {
        const double redNLL = nllval - 1000000.0;
        std::stringstream redNLL_ss;
        redNLL_ss << std::fixed << std::setprecision(20) << redNLL;

        std::streamsize ss = std::cout.precision();
        std::cout << std::fixed << std::setprecision(20);

        WriteInfoStatus("FittingTool::FitPDF", "   -> Reduced Final value of the NLL = " + redNLL_ss.str());
        WriteInfoStatus("FittingTool::FitPDF", "   -> Final value of the NLL = " + std::to_string(nllval));
        WriteInfoStatus("FittingTool::FitPDF", "   -> Final value of offset = " + std::to_string(nlloffset));
        WriteInfoStatus("FittingTool::FitPDF", "   -> Final NLL - offset = " + std::to_string(nllval-nlloffset));

        std::cout << resetiosflags( ios::fixed | ios::showpoint );
        std::cout << std::setprecision(ss);
    }
    if (m_debug < 1) std::cout.clear();
    return nllval;
}

//____________________________________________________________________________________
//
void FittingTool::ExportFitResultInTextFile( const std::string &fileName, const std::vector<std::string>& blinded )
{
    if(!m_fitResult){
        WriteErrorStatus("FittingTool::ExportFitResultInTextFile", "The FitResultObject seems not to be defined.");
        return;
    }

    //
    // Printing the nuisance parameters post-fit values
    //
    ofstream nuisParAndCorr(fileName);
    nuisParAndCorr << "NUISANCE_PARAMETERS\n";

    RooRealVar* var(nullptr);
    std::unique_ptr<TIterator> param(m_fitResult -> floatParsFinal().createIterator());
    while( (var = static_cast<RooRealVar*>(param->Next())) ){

        // Not consider nuisance parameter being not associated to syst (yet)
        std::string varname = var->GetName();
        TString vname=var->GetName();
        vname.ReplaceAll("alpha_","");

        const double pull  = var->getVal(); // GetValue() return value in unit of sigma
        const double errorHi = var->getErrorHi();
        const double errorLo = var->getErrorLo();

        if (blinded.size() == 0){
            FittingTool::CheckUnderconstraint(var);
            nuisParAndCorr << vname << "  " << pull << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "\n";
        } else {
            std::string vname_s = vname.Data();
            if (std::find(blinded.begin(), blinded.end(), vname_s) == blinded.end()){
                FittingTool::CheckUnderconstraint(var);
                nuisParAndCorr << vname << "  " << pull << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "\n";
            } else {
                const std::string& hex = Common::DoubleToPseudoHex(pull);
                nuisParAndCorr << vname << "  " << hex << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "\n";
            }
        }
    }

    //
    // Correlation matrix
    //
    TH2* h2Dcorrelation = m_fitResult -> correlationHist();
    nuisParAndCorr << "\n\nCORRELATION_MATRIX\n";
    nuisParAndCorr << h2Dcorrelation->GetNbinsX() << "   " << h2Dcorrelation->GetNbinsY() << "\n";
    for(int kk=1; kk < h2Dcorrelation->GetNbinsX()+1; kk++) {
        for(int ll=1; ll < h2Dcorrelation->GetNbinsY()+1; ll++) {
            nuisParAndCorr << h2Dcorrelation->GetBinContent(kk,ll) << "   ";
        }
        nuisParAndCorr << "\n";
    }

    //
    // NLL value
    //
    nuisParAndCorr << "\n\nNLL\n";
    nuisParAndCorr << m_fitResult -> minNll() << "\n";

    //
    // Closing the output file
    //
    nuisParAndCorr << "\n";
    nuisParAndCorr.close();
}

//____________________________________________________________________________________
//
std::map < std::string, double > FittingTool::ExportFitResultInMap(){
    std::map < std::string, double > result;
    if(!m_fitResult){
        WriteErrorStatus("FittingTool::ExportFitResultInMap", "The FitResultObject seems not to be defined.");
        return result;
    }
    RooRealVar* var(nullptr);
    std::unique_ptr<TIterator> param(m_fitResult -> floatParsFinal().createIterator());
    while( (var = static_cast<RooRealVar*>(param->Next())) ){
        // Not consider nuisance parameter being not associated to syst
        std::string varname = var->GetName();
        const double pull  = var->getVal();
        result.insert( std::pair < std::string, double >(varname, pull) );
    }
    return result;
}

//____________________________________________________________________________________
//
int FittingTool::GetGroupedImpact( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, RooWorkspace* ws, const std::string& categoryOfInterest, const std::string& outFileName ) const{
    // file to save results to
    std::ofstream outFile;
    outFile.open(outFileName.c_str());

    // obtain constrainedParams and poi like in FitPDF(), but make sure POI is not constant
    RooArgSet* constrainedParams = fitpdf->getParameters(*fitdata);
    RooStats::RemoveConstantParameters(constrainedParams);
    RooFit::Constrain(*constrainedParams);

    RooRealVar * poi = static_cast<RooRealVar*>(model->GetParametersOfInterest()->first());
    if(!poi){
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        WriteErrorStatus("FittingTool::GetGroupedImpact", "Cannot find the parameter of interest !");
        return -1;
    }
    poi -> setConstant(false);
    poi -> setVal(m_valPOI);
    if(!m_constPOI && m_randomize){
        poi->setVal( m_valPOI + m_randomNP*(gRandom->Uniform(2)-1.) );
    }

    // save snapshot of original workspace
    ws->saveSnapshot("snapshot_AfterFit_POI", *(model->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_AfterFit_NP" , *(model->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(model->GetGlobalObservables())    );

    std::vector<std::string> associatedParams; // parameters associated to a SubCategory

    // repeat the nominal fit - done so that the initial randomization is the exact same as for the following fit(s)
    // this should help avoid issues with fits ending up in different local minima for groups with very small impact on the POI
    FitExcludingGroup(false, false, fitdata, fitpdf, constrainedParams, model, ws, "Nominal", associatedParams);  // nothing held constant -> "snapshot_AfterFit_POI_Nominal"

    //
    // eventually do it once more
    if(TRExFitter::OPTION["GroupedImpactMoreFit"]>0){
        ws->saveSnapshot("snapshot_AfterFit_POI", *(model->GetParametersOfInterest()) );
        ws->saveSnapshot("snapshot_AfterFit_NP" , *(model->GetNuisanceParameters())   );
        ws->saveSnapshot("snapshot_AfterFit_GO" , *(model->GetGlobalObservables())    );
        FitExcludingGroup(false, false, fitdata, fitpdf, constrainedParams, model, ws, "Nominal", associatedParams);  // nothing held constant -> "snapshot_AfterFit_POI_Nominal"
    }

    // loop over unique SubCategories
    for (std::set<std::string>::iterator itCategories = m_subCategories.begin(); itCategories != m_subCategories.end(); ++itCategories){
        if(categoryOfInterest!="all" && *itCategories!=categoryOfInterest)
            continue; // if a category was specified via command line, only process that one

        WriteInfoStatus("FittingTool::GetGroupedImpact","performing grouped systematics impact evaluation for: " + *itCategories);

        // find all associated parameters per SubCategory
        associatedParams.clear();
        for(auto itSysts = m_subCategoryMap.begin(); itSysts != m_subCategoryMap.end(); ++itSysts) {
            if (itSysts->second == *itCategories) {
                //WriteDebugStatus("FittingTool::GetGroupedImpact","  associated parameter: " + itSysts->first);
                associatedParams.push_back(itSysts->first);
            }
        }

        // special case for gammas
        if(*itCategories=="Gammas") {
            FitExcludingGroup(true,  false, fitdata, fitpdf, constrainedParams, model, ws, *itCategories, associatedParams);
        }

        // special case for stat-only fit
        else if(*itCategories=="FullSyst") {
            FitExcludingGroup(true,  true,  fitdata, fitpdf, constrainedParams, model, ws, *itCategories, associatedParams);
        }

        // default: perform a fit where parameters in SubCategory are held constant
        else {
            FitExcludingGroup(false, false, fitdata, fitpdf, constrainedParams, model, ws, *itCategories, associatedParams);
        }
    }

    // load original workspace again
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    WriteInfoStatus("FittingTool::GetGroupedImpact","-----------------------------------------------------");

    // report replication of nominal fit
    ws->loadSnapshot("snapshot_AfterFit_POI_Nominal");
    WriteInfoStatus("FittingTool::GetGroupedImpact", "replicated nominal fit");
    WriteInfoStatus("FittingTool::GetGroupedImpact", "POI is:   " + std::to_string(poi->getVal()) + " +/- " + std::to_string(poi->getError()) +
                                                     "    ( +" + std::to_string(poi->getErrorHi()) + ", " + std::to_string(poi->getErrorLo()) + " )");

    const double NomUp2=(poi->getErrorHi()*poi->getErrorHi());
    const double NomLo2=(poi->getErrorLo()*poi->getErrorLo());
    const double Nom2  =(poi->getError()*poi->getError());

    // report impact calculations, impact is obtained by quadrature subtraction from replicated nominal fit
    for (std::set<std::string>::iterator itCategories = m_subCategories.begin(); itCategories != m_subCategories.end(); ++itCategories){
        if(categoryOfInterest!="all" && *itCategories!=categoryOfInterest)
            continue; // if a category was specified via command line, only process that one

        ws->loadSnapshot(("snapshot_AfterFit_POI_" + *itCategories).c_str());
        WriteInfoStatus("FittingTool::GetGroupedImpact","-----------------------------------------------------");
        WriteInfoStatus("FittingTool::GetGroupedImpact", "category: " + *itCategories + " (fixed to best-fit values for fit)");
        WriteInfoStatus("FittingTool::GetGroupedImpact", "POI is:   " + std::to_string(poi->getVal()) + " +/- " + std::to_string(poi->getError()) +
                                                         "    ( +" + std::to_string(poi->getErrorHi()) + ", " + std::to_string(poi->getErrorLo()) + " )");
        if(*itCategories=="FullSyst") WriteDebugStatus("FittingTool::GetGroupedImpact", "  (corresponds to a stat-only fit)");
        WriteInfoStatus("FittingTool::GetGroupedImpact", "           --> impact: " + std::to_string(std::sqrt(-(poi->getError()*poi->getError()) + Nom2)) +
                                                         "    ( +" + std::to_string(std::sqrt(- (poi->getErrorHi()*poi->getErrorHi()) + NomUp2)) + ", -" + std::to_string(std::sqrt(-(poi->getErrorLo()*poi->getErrorLo()) + NomLo2)) + " )" );

        // write results to file
        outFile << *itCategories << "    " << std::sqrt( -(poi->getError()*poi->getError()) + Nom2 ) << "  ( +" << std::sqrt( - (poi->getErrorHi()*poi->getErrorHi()) + NomUp2 ) << ", -" << std::sqrt( - (poi->getErrorLo()*poi->getErrorLo()) + NomLo2 ) << " )\n";
    }

    WriteInfoStatus("FittingTool::GetGroupedImpact", "-----------------------------------------------------");

    outFile.close();

    // load original workspace again
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    return 0;
}


//____________________________________________________________________________________
//
// perform a fit where all parameters in "affectedParams" (usually coming from SubGroup "category") are set to constant, optionally also gammas or all parameters
void FittingTool::FitExcludingGroup(bool excludeGammas, bool statOnly, RooAbsData*& fitdata, RooAbsPdf*& fitpdf, RooArgSet*& constrainedParams,
                                    RooStats::ModelConfig* mc, RooWorkspace* ws, const std::string& category, const std::vector<std::string>& affectedParams) const {

    // (VD): use this to fix nuisance parameter before the fit
    const RooArgSet* glbObs = mc->GetGlobalObservables();
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    WriteInfoStatus("FittingTool::FitExcludingGroup", "-----------------------------------------------------");
    WriteInfoStatus("FittingTool::FitExcludingGroup", "           breakdown for " + category);
    WriteInfoStatus("FittingTool::FitExcludingGroup", "-----------------------------------------------------");

    std::unique_ptr<TIterator> it(mc->GetNuisanceParameters()->createIterator());
    RooRealVar* var2(nullptr);

    while( (var2 = static_cast<RooRealVar*>(it->Next())) ){
        std::string varname = var2->GetName();

        // default: set everything non-constant (but the saturated-model norm-factors!)
        if (varname.find("saturated_model_sf_")==std::string::npos) var2->setConstant(0);
        else var2->setConstant(1);

        // if excludeGammas==true, set gammas to constant
        if (excludeGammas) {
            if (varname.find("gamma_stat")!=string::npos) var2->setConstant(1);
        }

        // set all affectedParams constant
        if (std::find(affectedParams.begin(), affectedParams.end(), varname) != affectedParams.end()) {
            var2->setConstant(1);
        }

        // for stat-only fits, set everything constant
        if (statOnly) {
            var2->setConstant(1);
        }
    }

    //constrainedParams->Print("v");
    // repeat the fit here ....
    std::unique_ptr<RooAbsReal> nll(fitpdf->createNLL(*fitdata,
                                    RooFit::Constrain(*constrainedParams),
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(1),
                                    NumCPU(TRExFitter::NCPU,RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE),
                                    RooFit::ExternalConstraints(*m_externalConstraints)
                                   ));
    RooMinimizer minim2(*nll);
    minim2.setStrategy(1);
    minim2.setPrintLevel(1); // set to -1 to reduce output
    minim2.setEps(1);
    const int status = minim2.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    RooRealVar* thePOI = static_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());

    const bool HessStatus = minim2.hesse();

    RooArgSet minosSet(*thePOI);
    if(m_useMinos && (m_varMinos.at(0)=="all" || find(m_varMinos.begin(),m_varMinos.end(),thePOI->GetName())<m_varMinos.end())){
        minim2.minos(minosSet);
    }

    if (status!=0) WriteErrorStatus("FittingTool::FitExcludingGroup", "unable to perform fit correctly! HessStatus: " + std::to_string(HessStatus));

    const double newPOIerr =thePOI->getError();
    const double newPOIerrU=thePOI->getErrorHi();
    const double newPOIerrD=thePOI->getErrorLo();

    const std::string snapshotName = "snapshot_AfterFit_POI_" + category;
    ws->saveSnapshot(snapshotName.c_str(), *mc->GetParametersOfInterest() );

    ws->loadSnapshot("snapshot_AfterFit_POI_Nominal");
    const double oldPOIerr =thePOI->getError();
    const double oldPOIerrU=thePOI->getErrorHi();
    const double oldPOIerrD=thePOI->getErrorLo();

    // check if uncertainties have increased compared to nominal fit, with 0.5% tolerance
    if ( (std::fabs(newPOIerrU)>std::fabs(oldPOIerrU)*1.005) || (std::fabs(newPOIerrD)>std::fabs(oldPOIerrD)*1.005) ) {
        WriteErrorStatus("FittingTool::FitExcludingGroup", "uncertainty has increased for " + category + "! please check the fit");
        WriteErrorStatus("FittingTool::FitExcludingGroup", "old: " + std::to_string(oldPOIerr) + " (+" + std::to_string(oldPOIerrU) + ", " + std::to_string(oldPOIerrD) + ")");
        WriteErrorStatus("FittingTool::FitExcludingGroup", "new: " + std::to_string(newPOIerr) + " (+" + std::to_string(newPOIerrU) + ", " + std::to_string(newPOIerrD) + ")");
    }
}

//____________________________________________________________________________________
//
// Check for underconstraints
void FittingTool::CheckUnderconstraint(const RooRealVar* const var) const {
    const std::string name = var->GetName();
    const double errorHi = var->getErrorHi();
    const double errorLo = var->getErrorLo();

    // dont check gamma parameters
    if (name.find("alpha_") == std::string::npos) return;

    if (errorHi > 1.001 || errorLo < -1.001){
        WriteWarningStatus("FittingTool::CheckUnderconstraint","NuisanceParameter: " + name + " is underconstrained! This may indicate fit convergence problems!");
    }
}
