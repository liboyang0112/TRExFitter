#include "TRExFitter/LikelihoodScanManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"

#include "TRandom3.h"

#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"

#include <algorithm>

LikelihoodScanManager::LikelihoodScanManager() :
    fScanMinX(-3),
    fScanMinY(-3),
    fStepsX(30),
    fScanMaxX(3),
    fScanMaxY(3),
    fStepsY(30),
    fUseOffset(true),
    fCPU(1),
    fParal2D(false),
    fParal2Dstep(-1)
{
    // shut-up RooFit!
    if(TRExFitter::DEBUGLEVEL<=1){
        if(TRExFitter::DEBUGLEVEL<=0) gErrorIgnoreLevel = kError;
        else gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Generation);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::LinkStateMgmt);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Optimization);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Tracing);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Contents);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
        RooMsgService::instance().setStreamStatus(1,false);
    }

}

//__________________________________________________________________________________
//
LikelihoodScanManager::~LikelihoodScanManager() {
}
    
//__________________________________________________________________________________
//
LikelihoodScanManager::scanResult1D LikelihoodScanManager::Run1DScan(const RooWorkspace* ws,
                                                                     const std::string& varName,
                                                                     RooDataSet* data) const {

    if (!ws || !data) {
        WriteErrorStatus("LikelihoodScanManager::Run1DScan", "Passed nullptr for ws");
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("LikelihoodScanManager::Run1DScan", "Passed nullptr for mc");
        exit(EXIT_FAILURE);
    }
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    double min(-3);
    double max(3);

    RooRealVar* var(nullptr);
    bool found(false);
    {
        std::unique_ptr<TIterator> it(mc->GetParametersOfInterest()->createIterator());
        while( (var = static_cast<RooRealVar*>(it->Next()))){
            const std::string vname = var->GetName();
            if (vname == varName) {
                WriteInfoStatus("LikelihoodScanManager::Run1DScan", "GetLikelihoodScan for NP = " + vname);
                found=true;
                min = var->getMin();
                max = var->getMax();
                break;
            }
        }
    }

    if (!found) {
        std::unique_ptr<TIterator> it(mc->GetNuisanceParameters()->createIterator());
        while( (var = static_cast<RooRealVar*>(it->Next()))){
            const std::string vname = var->GetName();
            if (vname == varName || vname == "alpha_"+varName) {
                WriteInfoStatus("LikelihoodScanManager::Run1DScan", "GetLikelihoodScan for NP = " + vname);
                found=true;
                break;
            }
        }
    }

    if (fScanMinX < 99999) { // is actually set
        min = fScanMinX;
    }

    if (fScanMaxX > -99999) { // is actually set
        max = fScanMaxX;
    }

    LikelihoodScanManager::scanResult1D result;
    if (!found) {
        WriteWarningStatus("LikelihoodScanManager::Run1DScan", "Cannot find NP: " + varName);
        return result;
    }

    result.first.resize(fStepsX);
    result.second.resize(fStepsX);

    std::unique_ptr<RooAbsReal> nll(simPdf->createNLL(*data,
                                                      RooFit::Constrain(*mc->GetNuisanceParameters()),
                                                      RooFit::Offset(fUseOffset),
                                                      NumCPU(fCPU, RooFit::Hybrid),
                                                      RooFit::Optimize(kTRUE)));
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this
    
    RooMinimizer m(*nll); // get MINUIT interface of fit
    m.setPrintLevel(-1);
    m.setStrategy(1);
    m.optimizeConst(2);
    m.setEps(tol);
    var->setConstant(kTRUE); // make POI constant in the fit
    
    double mnll = 9999999;
    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {
        WriteInfoStatus("LikelihoodScanManager::Run1DScan","Running LHscan for point " + std::to_string(ipoint+1) + " out of " + std::to_string(fStepsX) + " points");
        result.first[ipoint] = min+ipoint*(max-min)/(fStepsX - 1);
        *var = result.first[ipoint]; // set POI
        m.migrad(); // minimize again with new posSigXsecOverSM value
        std::unique_ptr<RooFitResult> r(m.save()); // save fit result
        result.second[ipoint] = r->minNll();
        if (result.second[ipoint] < mnll) mnll = result.second[ipoint];
    }
    var->setConstant(kFALSE);

    for (auto & iY : result.second) {
        iY = iY - mnll;
    }

    TRandom3 rand(1234567);
    if (std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varName) != fBlindedParameters.end()) {
        const double rndNumber = rand.Uniform(5);
        for (auto& ix : result.first) {
            ix += rndNumber;
        }
    }

    return result;
}

//__________________________________________________________________________________
//
LikelihoodScanManager::Result2D LikelihoodScanManager::Run2DScan(const RooWorkspace* ws,
                                                                 const std::pair<std::string, std::string>& varNames,
                                                                 RooDataSet* data) const {

    if (!ws || !data) {
        WriteErrorStatus("LikelihoodScanManager::Run2DScan", "Passed nullptr for ws");
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("LikelihoodScanManager::Run2DScan", "Passed nullptr for mc");
        exit(EXIT_FAILURE);
    }
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    int count = 0;
    RooRealVar* varX = nullptr;
    RooRealVar* varY = nullptr;
    {
        std::unique_ptr<TIterator> it(mc->GetNuisanceParameters()->createIterator());
        RooRealVar* var_tmp(nullptr);
        while ( (var_tmp = static_cast<RooRealVar*>(it->Next())) ){
            const std::string vname = var_tmp->GetName();
            if (vname == varNames.first || vname == "alpha_"+varNames.first){
                varX = var_tmp;
                ++count;
            }
            if (vname == varNames.second || vname == "alpha_"+varNames.second){
                varY = var_tmp;
                ++count;
            }
            if (count == 2) break;
        }
    }

    if (count != 2) {
        std::unique_ptr<TIterator> it(mc->GetParametersOfInterest()->createIterator());
        RooRealVar* var_tmp(nullptr);
        while ( (var_tmp = static_cast<RooRealVar*>(it->Next())) ){
            const std::string vname = var_tmp->GetName();
            if (vname == varNames.first || vname == "alpha_"+varNames.first){
                varX = var_tmp;
                ++count;
            }
            if (vname == varNames.second || vname == "alpha_"+varNames.second){
                varY = var_tmp;
                ++count;
            }
            if (count == 2) break;
        }
    }

    LikelihoodScanManager::Result2D result;

    if (count != 2) {
        WriteWarningStatus("LikelihoodScanManager::Run2DScan","Did not find the two parameters you want to use in the 2D likelihood scan");
        return result;
    }

    double minValX = varX->getMin();
    double maxValX = varX->getMax();
    double minValY = varY->getMin();
    double maxValY = varY->getMax();
    
    if (fScanMinX < 99999) { // is actually set
        minValX = fScanMinX;
    }
    if (fScanMinY < 99999) { // is actually set
        minValY = fScanMinY;
    }
    if (fScanMaxX > -99999) { // is actually set
        maxValX = fScanMaxX;
    }
    if (fScanMaxY > -99999) { // is actually set
        maxValY = fScanMaxY;
    }
    
    std::unique_ptr<RooAbsReal> nll(simPdf->createNLL(*data,
                                                      RooFit::Constrain(*mc->GetNuisanceParameters()),
                                                      RooFit::Offset(fUseOffset),
                                                      NumCPU(fCPU, RooFit::Hybrid),
                                                      RooFit::Optimize(kTRUE)));
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance();

    RooMinimizer m(*nll); // get MINUIT interface of fit
    m.optimizeConst(2);
    m.setErrorLevel(-1);
    m.setPrintLevel(-1);
    m.setStrategy(1); // set precision to high
    m.setEps(tol);
    //Set both POIs to constant
    varX->setConstant(kTRUE); // make POI constant in the fit
    varY->setConstant(kTRUE); // make POI constant in the fit

    result.x.resize(fStepsX);
    result.y.resize(fStepsY);
    result.z.resize(fStepsX);
    for (auto& iz : result.z) {
        iz.resize(fStepsY);
    }

    //values for parameter1, parameter2 and the NLL value
    double zmin = 9999999;
    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {
        if (fParal2D && ipoint != fParal2Dstep) continue;
        WriteInfoStatus("LikelihoodScanManager::Run2DScan","Running LHscan for point " + std::to_string(ipoint+1) + " out of " + std::to_string(fStepsX) + " points");
        result.x[ipoint] = minValX + ipoint * (maxValX - minValX) / (fStepsX - 1);
        *varX = result.x[ipoint]; // set POI
        for (int jpoint = 0; jpoint < fStepsY; ++jpoint) {
            WriteInfoStatus("LikelihoodScanManager::Run2DScan","Running LHscan for subpoint " + std::to_string(jpoint+1) + " out of " + std::to_string(fStepsY) + " points");
            result.y[jpoint] = minValY + jpoint * (maxValY - minValY) / (fStepsY - 1);
            *varY = result.y[jpoint]; // set POI
            m.migrad(); // minimize again with new posSigXsecOverSM value
            std::unique_ptr<RooFitResult> r(m.save()); // save fit result
            const double z_tmp = r->minNll();
            result.z[ipoint][jpoint] = z_tmp;

            // save the best values
            if (z_tmp < zmin) {
                zmin = z_tmp;
            }
        }
    }
    varX->setConstant(kFALSE);
    varY->setConstant(kFALSE);

    TRandom3 rand(1234567);
    const bool blindX = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.first) != fBlindedParameters.end();
    const bool blindY = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.second) != fBlindedParameters.end();
    const double rndNumber = rand.Uniform(5);
    if (blindX) {
        for (auto& ix : result.x) {
            ix += rndNumber;
        }
    }

    if (blindY) {
        for (auto& iy : result.y) {
            iy += rndNumber;
        }
    }

    return result;
}
