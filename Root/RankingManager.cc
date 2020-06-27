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
    
    std::map< std::string,double > muVarUp;
    std::map< std::string,double > muVarDown;
    std::map< std::string,double > muVarNomUp;
    std::map< std::string,double > muVarNomDown;

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
        double central = fitResults -> GetNuisParValue(  npName);
        double up      = fitResults -> GetNuisParErrUp(  npName);
        double down    = fitResults -> GetNuisParErrDown(npName);
        
        outFile << iNP.first << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";
        //
        // Experimental: reduce the range of ranking
        if(TRExFitter::OPTION["ReduceRanking"]!=0){
            up   *= TRExFitter::OPTION["ReduceRanking"];
            down *= TRExFitter::OPTION["ReduceRanking"];
        }
        //
        // Set the NP to its post-fit *up* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool.ResetFixedNP();
        // fix NPs that are fixed in the config
        // this has nothing to do with fixing NPs for the ranking
        // this is just needed to be compatible with the normal fit
        if(fFitFixedNPs.size()>0){
            for(const auto& nuisParToFix : fFitFixedNPs){
                fitTool.FixNP(nuisParToFix.first,nuisParToFix.second);
            }
        }

        fitTool.FixNP( iNP.first, central + std::abs(up));
        fitTool.FitPDF( mc, simPdf, data );
        muVarUp[iNP.first]   = (fitTool.ExportFitResultInMap())[fPOIName];

        // Set the NP to its post-fit *down* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool.ResetFixedNP();
        fitTool.FixNP(iNP.first, central - std::abs(down));
        if(fFitFixedNPs.size()>0){
            for(const auto& nuisParToFix : fFitFixedNPs){
                fitTool.FixNP(nuisParToFix.first,nuisParToFix.second);
            }
        }
        fitTool.FitPDF( mc, simPdf, data );
        muVarDown[iNP.first] = (fitTool.ExportFitResultInMap())[fPOIName];

        double dMuUp   = muVarUp[iNP.first]-muhat;
        double dMuDown = muVarDown[iNP.first]-muhat;

        // Experimental: reduce the range of ranking
        if(TRExFitter::OPTION["ReduceRanking"]!=0){
            dMuUp   /= TRExFitter::OPTION["ReduceRanking"];
            dMuDown /= TRExFitter::OPTION["ReduceRanking"];
        }

        outFile << dMuUp << "   " << dMuDown << "  ";

        if(iNP.second){
            muVarNomUp[  iNP.first] = muhat;
            muVarNomDown[iNP.first] = muhat;
        }
        else{
            up   = 1.;
            down = 1.;

            // Experimental: reduce the range of ranking
            if(TRExFitter::OPTION["ReduceRanking"]!=0){
                up   *= TRExFitter::OPTION["ReduceRanking"];
                down *= TRExFitter::OPTION["ReduceRanking"];
            }

            // Set the NP to its pre-fit *up* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool.ResetFixedNP();
            fitTool.FixNP( iNP.first, central + std::abs(up));
            fitTool.FitPDF( mc, simPdf, data );
            if(fFitFixedNPs.size()>0){
                for(const auto& nuisParToFix : fFitFixedNPs){
                    fitTool.FixNP(nuisParToFix.first,nuisParToFix.second);
                }
            }
            muVarNomUp[iNP.first]   = (fitTool.ExportFitResultInMap())[fPOIName];
            //
            // Set the NP to its pre-fit *down* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool.ResetFixedNP();
            fitTool.FixNP(iNP.first, central - std::abs(down));
            if(fFitFixedNPs.size()>0){
                for(const auto& nuisParToFix : fFitFixedNPs){
                    fitTool.FixNP(nuisParToFix.first,nuisParToFix.second);
                }
            }
            fitTool.FitPDF( mc, simPdf, data );
            //
            muVarNomDown[iNP.first] = (fitTool.ExportFitResultInMap())[fPOIName];
        }
        dMuUp   = muVarNomUp[iNP.first]-muhat;
        dMuDown = muVarNomDown[iNP.first]-muhat;
        //
        // Experimental: reduce the range of ranking
        if(TRExFitter::OPTION["ReduceRanking"]!=0){
            dMuUp   /= TRExFitter::OPTION["ReduceRanking"];
            dMuDown /= TRExFitter::OPTION["ReduceRanking"];
        }
        //
       outFile << dMuUp << "   " << dMuDown << " " << std::endl;

    }
 
    ws->loadSnapshot("tmp_snapshot");
    outFile.close();
}
