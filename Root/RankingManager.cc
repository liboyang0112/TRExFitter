#include "TRExFitter/RankingManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/StatusLogbook.h"

#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooSimultaneous.h"

#include <algorithm>
#include <fstream>

RankingManager::RankingManager() :
    fOutputPath(""),
    fInjectGlobalObservables(false),
    fFitStrategy(1),
    fPOIName("POI"),
    fCPU(1),
    fRndRange(0.1),
    fUseRnd(false),
    fRndSeed(-999),
    fStatOnly(false)
{
}

RankingManager::~RankingManager() {
}

void RankingManager::AddNuisPar(const std::string& name, const bool isNF) {

    // check if the NP already exists by checking the names
    auto it = std::find_if(fNuisPars.begin(), fNuisPars.end(),
                [&name](const std::pair<std::string, bool>& element){return element.first == name;});

    if (it != fNuisPars.end()) {
        WriteWarningStatus("RankingManager::AddNuisPar", "NP " + name + " already exists in the list, not adding it");
        return;
    }
    if(name.find("_bin_") != std::string::npos) {
        fNuisPars.emplace_back("gamma_" + name, isNF);
    } else {
        fNuisPars.emplace_back(name, isNF);
    }
}

void RankingManager::RunRanking(FitResults* fitResults,
                                RooWorkspace* ws,
                                RooDataSet* data,
                                const std::vector<std::shared_ptr<NormFactor> >& nfs) const {

    if (fOutputPath == "") {
        WriteErrorStatus("RankingManager::RunRanking", "OutputPath not set, plese set it via SetOutputPath()");
        exit(EXIT_FAILURE);
    }

    if (!ws) {
        WriteErrorStatus("RankingManager::RunRanking", "Workspace is nullptr");
        exit(EXIT_FAILURE);
    }

    std::ofstream outFile(fOutputPath.c_str());

    if (!outFile.good() || !outFile.is_open()) {
        WriteErrorStatus("RankingManager::RunRanking", "Cannot open file at " + fOutputPath);
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());
    if (!mc || !simPdf || !data){
        WriteErrorStatus("TRExFit::ProduceNPRanking","At least one of the objects that is needed to run ranking is not present");
        exit(EXIT_FAILURE);
    }
    
    if (fInjectGlobalObservables && !fFitValues.empty()) {
        FitUtils::InjectGlobalObservables(ws, fFitValues);
    }

    ws->saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));
   
    double poiInitial(0.);
    for(const auto& inf : nfs) {
        if (inf->fName == fPOIName) {
            poiInitial = inf->fNominal;
            break;
        }
    }
 
    FittingTool fitTool{};
    fitTool.SetUseHesse(false);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetDebug(TRExFitter::DEBUGLEVEL);
    fitTool.ValPOI(poiInitial);
    fitTool.SetNCPU(fCPU);
    fitTool.ConstPOI(false);
    if(fStatOnly){
        fitTool.NoGammas();
        fitTool.NoSystematics();
    }
    fitTool.SetRandomNP(fRndRange, fUseRnd, fRndSeed);

    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : nfs) {
        if (inf->fName == fPOIName) continue;
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->fNominal);
    }
    fitTool.SetNPs(npNames,npValues);
    
    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, nfs);
    
    const double muhat = fitResults -> GetNuisParValue(fPOIName);
    for(const auto& iNP : fNuisPars){

        std::string npName = iNP.first;
        if (npName.find("_bin_") != std::string::npos) {
            npName = Common::ReplaceString(npName, "gamma_", "");
        }
        //
        // Getting the postfit values of the nuisance parameter
        const double central = fitResults -> GetNuisParValue(  npName);
        const double up      = fitResults -> GetNuisParErrUp(  npName);
        const double down    = fitResults -> GetNuisParErrDown(npName);
        
        outFile << iNP.first << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";
        //
        // Experimental: reduce the range of ranking

        RankingManager::RankingValues values;
        values.central = central;
        values.up = up;
        values.down = down;

        //
        double dMuUp   = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, true, false, values, muhat);
        double dMuDown = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, false, false, values, muhat);

        outFile << dMuUp << "   " << dMuDown << "  ";

        dMuUp   = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, true, true, values, muhat);
        dMuDown = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, false, true, values, muhat);
        
        outFile << dMuUp << "   " << dMuDown << " " << std::endl;

    }
 
    ws->loadSnapshot("tmp_snapshot");
    outFile.close();
}

double RankingManager::RunSingleFit(FittingTool* fitTool,
                                    RooWorkspace* ws,       
                                    RooStats::ModelConfig *mc,
                                    RooSimultaneous *simPdf,
                                    RooDataSet* data,
                                    const std::pair<std::string, bool>& np,
                                    const bool isUp,
                                    const bool isPrefit,
                                    const RankingManager::RankingValues& values,
                                    const double muhat) const {

    if (isPrefit && np.second) {
        return 0;
    }    

    ws->loadSnapshot("tmp_snapshot");
    fitTool->ResetFixedNP();
    if(fFitFixedNPs.size()>0){
        for(const auto& nuisParToFix : fFitFixedNPs){
            fitTool->FixNP(nuisParToFix.first,nuisParToFix.second);
        }
    }

    double shift = 1.0;
    if (isPrefit) {
       shift = isUp ? 1. : -1.;
    } else {
       shift = isUp ? std::fabs(values.up) : -std::fabs(values.down);
    }

    if(TRExFitter::OPTION["ReduceRanking"]!=0){
        shift *= TRExFitter::OPTION["ReduceRanking"];
    }

    fitTool->FixNP( np.first, values.central + shift);
    fitTool->FitPDF( mc, simPdf, data );

    double result = (fitTool->ExportFitResultInMap())[fPOIName] - muhat;

    if(TRExFitter::OPTION["ReduceRanking"]!=0){
        result /= TRExFitter::OPTION["ReduceRanking"];
    }

    return result;
}
