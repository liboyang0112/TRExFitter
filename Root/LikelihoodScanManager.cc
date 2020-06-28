#include "TRExFitter/LikelihoodScanManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/StatusLogbook.h"

#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"

LikelihoodScanManager::LikelihoodScanManager() :
    fScanMinX(-3),
    fScanMinY(-3),
    fStepsX(30),
    fScanMaxX(3),
    fScanMaxY(3),
    fStepsY(30),
    fUseOffset(true),
    fCPU(1)
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

LikelihoodScanManager::~LikelihoodScanManager() {
}
    
LikelihoodScanManager::scanResult1D LikelihoodScanManager::Run1DScan(const RooWorkspace* ws,
                                                                     const std::string& varName,
                                                                     RooDataSet* data,
                                                                     const std::vector<std::shared_ptr<NormFactor> >& nfs) const {

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

    for (const auto& inf : nfs) {
        if (inf->fName == varName) {
            min = inf->fMin;
            max = inf->fMax;
            break;
        }
    }

    if (fScanMinX < 99999) { // is actually set
        min = fScanMinX;
    }

    if (fScanMaxX > -99999) { // is actually set
        max = fScanMaxX;
    }

    RooRealVar* var(nullptr);
    bool found(false);
    {
        std::unique_ptr<TIterator> it(mc->GetParametersOfInterest()->createIterator());
        while( (var = static_cast<RooRealVar*>(it->Next()))){
            const std::string vname = var->GetName();
            if (vname == varName) {
                WriteInfoStatus("LikelihoodScanManager::Run1DScan", "GetLikelihoodScan for NP = " + vname);
                found=true;
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

    return result;
}
