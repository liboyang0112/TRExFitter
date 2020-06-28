#include "TRExFitter/FitUtils.h"

#include "TRExFitter/FittingTool.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/StatusLogbook.h"

#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMultiVarGaussian.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "TVectorD.h"
#include "TMatrixDSym.h"

#include <sstream>

//__________________________________________________________________________________
//
void FitUtils::ApplyExternalConstraints(RooWorkspace* ws,
                                        FittingTool* fitTool,
                                        RooSimultaneous* simPdf,
                                        const std::vector<std::shared_ptr<NormFactor> >& normFactors) {

    // Tikhonov regularization (for unfolding)
    RooArgList l;
    std::vector<double> nomVec;
    std::vector<double> tauVec;
    std::vector<std::string> names; 
    for(const auto& nf : normFactors) {
        if(nf->fTau == 0) continue;
        // only add the unique NFs
        if (std::find(names.begin(), names.end(), nf->fName) != names.end()) continue;
        names.emplace_back(nf->fName);
        
        if(ws->var(nf->fName.c_str())) {
            l.add(*ws->var(nf->fName.c_str()));
        } else if (ws->function(nf->fName.c_str())) {
            l.add(*ws->function(nf->fName.c_str()));
        } else {
            WriteWarningStatus("FitUtils::ApplyExternalConstraints","Cannot apply tau to norm factor " + nf->fName);
            continue;
        }
        
        nomVec.push_back( nf->fNominal );
        tauVec.push_back( nf->fTau );
    }

    if(tauVec.empty()) return;

    TVectorD nominal(nomVec.size());
    TMatrixDSym cov(tauVec.size());
    for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
        nominal(i_tau) = nomVec[i_tau];
        cov(i_tau,i_tau) = (1./tauVec[i_tau]) * (1./tauVec[i_tau]);
    }
    RooMultiVarGaussian r("regularization","regularization",l,nominal,cov);
    ws->import(r);
    ws->defineSet("myConstraints","regularization");
    simPdf->setStringAttribute("externalConstraints","myConstraints");

    if(simPdf->getStringAttribute("externalConstraints")){
        WriteInfoStatus("FitUtils::ApplyExternalConstraints",Form("Building NLL with external constraints %s",simPdf->getStringAttribute("externalConstraints")));
        const RooArgSet* externalConstraints = ws->set(simPdf->getStringAttribute("externalConstraints"));
        fitTool->SetExternalConstraints( externalConstraints );
    }
}

//__________________________________________________________________________________
//
void FitUtils::SetBinnedLikelihoodOptimisation(RooWorkspace* ws) {
    RooFIter rfiter = ws->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = rfiter.next())) {
        if (arg->IsA() == RooRealSumPdf::Class()) {
            arg->setAttribute("BinnedLikelihood");
            const std::string temp_string = arg->GetName();
            WriteDebugStatus("FitUtils::SetBinnedLikelihoodOptimisation", "Activating binned likelihood attribute for " + temp_string);
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::InjectGlobalObservables(RooWorkspace* ws,
                                       const std::map< std::string, double >& npValues) {
    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("FitUtils::InjectGlobalObservables", "Cannot read ModelCOnfig");
        exit(EXIT_FAILURE);
    }
    RooArgSet mc_globs = *mc->GetGlobalObservables();

    WriteInfoStatus("FitUtils::InjectGlobalObservables", "Injecting the following NP values to global observables");
    for(const auto& np_value : npValues) {
        const std::string this_name = np_value.first;
        double this_value = np_value.second;
    
        std::ostringstream tmp;
        tmp << this_name << ": " << this_value;
        WriteInfoStatus("FitUtils::InjectGlobalObservables", tmp.str());
    
        // find the corresponding glob
        const std::string glob_name = "nom_" + this_name;
        const std::string glob_name_alpha = "nom_alpha_" + this_name;
        const std::string glob_name_gamma = "nom_gamma_" + this_name;
        std::unique_ptr<TIterator> gIter(mc_globs.createIterator());
        RooRealVar* glob = nullptr;
        RooRealVar* this_glob = nullptr;
        while ((glob = static_cast<RooRealVar*>(gIter->Next()))) {
            if(glob->GetName() == glob_name || glob->GetName() == glob_name_alpha || glob->GetName() == glob_name_gamma) {
                this_glob = glob;
                break;
            }
        }

        if(!this_glob) {
            WriteErrorStatus("FitUtils::InjectGlobalObservables", "Could not find global observable "+glob_name);
            continue;
        }

        // set gamma values to gamma*nom_gamma
        // cppcheck-suppress stlIfStrFind
        if(glob_name.find("nom_gamma_") == 0) {
            this_value = this_value * this_glob->getVal();
        }

        this_glob->setVal(this_value);
    }
}
