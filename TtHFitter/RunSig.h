#ifndef RUNSIG_H
#define RUNSIG_H

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooMinimizerFcn.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"

#include <string>

RooDataSet* makeAsimovData(RooStats::ModelConfig* mc, bool doConditional, RooWorkspace* w, RooNLLVar* conditioning_nll, double mu_val, std::string* mu_str, std::string* mu_prof_str, double mu_val_profile, bool doFit);

int minimize(RooNLLVar* nll, RooWorkspace* combWS = NULL);

void RunSig(const char* inFileName,
    const char* wsName = "combined",
    const char* modelConfigName = "ModelConfig",
    const char* dataName = "obsData",
    const char* asimov1DataName = "asimovData_1",
    const char* conditional1Snapshot = "conditionalGlobs_1",
    const char* nominalSnapshot = "nominalGlobs",
    std::string smass = "130",
    std::string folder = "test");

RooDataSet* makeData(RooDataSet* orig, RooSimultaneous* simPdf, const RooArgSet* observables, RooRealVar* firstPOI, double mass, double& mu_min);

void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
#endif
