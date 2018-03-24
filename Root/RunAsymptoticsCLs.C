/*
Author: Aaron Armbruster
Date:       2012-05-25
Email:  armbrusa@umich.edu
Description: Script to run asymptotic CLs.

--------
00-01-00
-First version with updated bands

--------
00-01-01
-Fixed problem in asimov data creation that affected +1,2sigma bands

--------
00-01-02
-(Re)added support for non-sim pdfs (still need to be extended)
-Fixed default doFit arg of makeAsimovData
-Added better output for unresolved fit failures
-Improved retry loop for fit failures


/////////////////////
//////PREAMBLE///////
/////////////////////

The script uses an iterative process to find the crossing of qmu with the qmu95(mu/sigma) curve,
where qmu95(mu/sigma) is found assuming asymptotic formula for the distribution of the
test statistic f(qmu|mu') (arxiv 1007.1727) and of the test statistic qmu (or tilde)

The sequence is

mu_i+1 = mu_i - gamma_i*(mu_i - mu'_i)

where gamma_i is a dynamic damping factor used for convergence (nominal gamma_i = 1), and mu'_i is
determined by extrapolating the test statistic to the qmu95 curve assuming qmu is parabolic:

qmu'_i = (mu'_i - muhat)^2 / sigma_i^2 = qmu95(mu'_i / sigma_i)

where sigma_i is determined by computing qmu_i (not '):

sigma_i = (mu_i - muhat) / sqrt(qmu_i)

At the crossing qmu_N = qmu95 the assumption that qmu is a parabola goes away,
so we're not ultimately dependent on this assumption beyond its use in the asymptotic formula.

The sequence ends when the relative correction factor gamma*(mu_i - mu'_i) / mu_i is less than some
specified precision (0.005 by default)




///////////////////////////
//////AFTER RUNNING////////
///////////////////////////


The results will be printed as well as stored in a root file in the folder 'root-files/<folder>', where <folder>
is specified by you (default 'test')

The root file has a 7-bin TH1D named 'limit', where each bin is filled with the upper limit values in this order:

1: Observed
2: Median
3: +2 sigma
4: +1 sigma
5: -1 sigma
6: -2 sigma
7: mu=0 fit status (only meaningful if asimov data is generated within the macro)

It will also store the result of the old bands procedure in a TH1D named 'limit_old'.




//////////////////////////


This version is functionally fully consistent with the previous tag.

NOTE: The script runs significantly faster when compiled
*/
#include "TtHFitter/RunAsymptoticsCLs.h"
#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/Common.h"

#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"

#include "RooWorkspace.h"
#include "RooNLLVar.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "RooMinimizerFcn.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooProduct.h"
#include "RooRealSumPdf.h"

#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace RooFit;
using namespace RooStats;

//band configuration
bool betterBands             = 1; //1 // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
bool betterNegativeBands     = 0; // (recommendation = 0) improve also the negative bands
bool profileNegativeAtZero   = 0; // (recommendation = 0) profile asimov for negative bands at zero


//other configuration
std::string defaultMinimizer     = "Minuit2";         // or "Minuit"
int defaultPrintLevel            = -1;                      // Minuit print level
int defaultStrategy              = 1; //1                           // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
bool killBelowFatal              = 1;                           // In case you want to suppress RooFit warnings further, set to 1
bool doBlind                     = 0;                           // in case your analysis is blinded
bool conditionalExpected         = 1 && !doBlind; // Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs
bool doTilde                     = 1;                           // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
bool doExp                       = 1;                           // compute expected limit
bool doObs                       = 1 && !doBlind; // compute observed limit
double precision                 = 0.005; // 0.005               // % precision in mu that defines iterative cutoff
bool usePredictiveFit            = 0; // 0                      // experimental, extrapolate best fit nuisance parameters based on previous fit results
bool extrapolateSigma            = 0; // 0                      // experimantal, extrapolate sigma based on previous fits
int maxRetries                   = 3; //3                           // number of minimize(fcn) retries before giving up
int numCPU                       = 8;                           // added by Michele

//don't touch!
std::map<RooNLLVar*, double> map_nll_muhat;
std::map<RooNLLVar*, double> map_muhat;
std::map<RooDataSet*, RooNLLVar*> map_data_nll;
std::map<RooNLLVar*, std::string> map_snapshots;
std::map<RooNLLVar*, std::map<double, double> > map_nll_mu_sigma;
RooWorkspace* w = NULL;
RooStats::ModelConfig* mc = NULL;
RooDataSet* data = NULL;
RooRealVar* firstPOI = NULL;
RooNLLVar* asimov_0_nll = NULL;
RooNLLVar* obs_nll = NULL;
int nrMinimize=0;
int direction=1;
int global_status=0;
double target_CLs=0.05;
// range of firstPOI from ModelConfig mc
double firstPOIMax = 0;
double firstPOIMin = 0;

void LimitsCLs::RunAsymptoticsCLs(const char* infile,
                                  const char* workspaceName,
                                  const char* modelConfigName,
                                  const char* dataName,
                                  string option,
                                  const char* asimovDataName,
                                  const char* conditionalSnapshot,
                                  const char* nominalSnapshot,
                                  string folder,
                                  string mass,
                                  double CL){
    stringstream smass;
    smass << mass;

    if(option.find("blind")!=string::npos) doBlind = 1;

    conditionalSnapshot = ""; // warningless compile
    nominalSnapshot = "";           // warningless compile

    RunAsymptoticsCLs(infile, workspaceName, modelConfigName, dataName, asimovDataName, folder, smass.str(), CL);
}



void LimitsCLs::RunAsymptoticsCLs(const char* infile,
                     const char* workspaceName,
                     const char* modelConfigName,
                     const char* dataName,
                     const char* asimovDataName,
                     string folder,
                     string mass,
                     double CL){
    TStopwatch timer;
    timer.Start();

    if (TtHFitter::DEBUGLEVEL < 2){
        gErrorIgnoreLevel = kError;
        RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }

    if (killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizer.c_str());
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(defaultPrintLevel);
    //RooNLLVar::SetIgnoreZeroEntries(1);

//check inputs
    TFile f(infile);
    w = (RooWorkspace*)f.Get(workspaceName);
    if (!w){
        std::string s = workspaceName;
        WriteErrorStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Workspace: " + s + " doesn't exist!");
        return;
    }
    RooFIter rfiter = w->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = rfiter.next())) {
        if (arg->IsA() == RooRealSumPdf::Class()) {
            arg->setAttribute("BinnedLikelihood");
        }
    }

    mc = (ModelConfig*)w->obj(modelConfigName);
    if (!mc){
        std::string s = modelConfigName;
        WriteErrorStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "ModelConfig: " + s + " doesn't exist!");
        return;
    }
    firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
    firstPOIMax = firstPOI->getMax();
    firstPOIMin = firstPOI->getMin();
    WriteDebugStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "get min and max of firstPOI");
    WriteDebugStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "firstPOIMin " + std::to_string(firstPOIMin));
    WriteDebugStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "firstPOIMax " + std::to_string(firstPOIMax));

    data = (RooDataSet*)w->data(dataName);
    if (!data){
        WriteErrorStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Dataset: doesn't exist!");
        return;
    }

    //RooAbsPdf* pdf = mc->GetPdf();
    obs_nll = createNLL(data);//(RooNLLVar*)pdf->createNLL(*data);
    map_snapshots[obs_nll] = "nominalGlobs";
    map_data_nll[data] = obs_nll;
    w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
    w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());

    global_status=0;
    RooDataSet* asimovData_0 = (RooDataSet*)w->data(asimovDataName);
    if (!asimovData_0)
    {
        asimovData_0 = makeAsimovData(conditionalExpected, obs_nll, 0);

        //asimovData_0 = makeAsimovData2((conditionalExpected ? obs_nll : (RooNLLVar*)NULL), 0., 0.);
    }
    int asimov0_status=global_status;

    asimov_0_nll = createNLL(asimovData_0);//(RooNLLVar*)pdf->createNLL(*asimovData_0);
    map_snapshots[asimov_0_nll] = "conditionalGlobs_0";
    map_data_nll[asimovData_0] = asimov_0_nll;
    setMu(0);
    map_muhat[asimov_0_nll] = 0;
    saveSnapshot(asimov_0_nll, 0);
    w->loadSnapshot("conditionalNuis_0");
    w->loadSnapshot("conditionalGlobs_0");
    map_nll_muhat[asimov_0_nll] = asimov_0_nll->getVal();


    target_CLs=1-CL;


    double med_limit = doExp ? getLimit(asimov_0_nll, 1.0) : 1.0;
    int med_status=global_status;

    double sigma = med_limit/sqrt(3.84); // pretty close
    double mu_up_p2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 2), 1) + 2);
    double mu_up_p1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 1), 1) + 1);
    double mu_up_n1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-1), 1) - 1);
    double mu_up_n2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-2), 1) - 2);

    double mu_up_p2 = mu_up_p2_approx;
    double mu_up_p1 = mu_up_p1_approx;
    double mu_up_n1 = mu_up_n1_approx;
    double mu_up_n2 = mu_up_n2_approx;

    firstPOI->setRange(-5*sigma, 5*sigma);
    map<int, int> N_status;
    if (betterBands && doExp){ // no better time than now to do this
        //find quantiles, starting with +2, since we should be at +1.96 right now

        double init_targetCLs = target_CLs;
        firstPOI->setRange(-5*sigma, 5*sigma);
        for (int N=2;N>=-2;N--){
            if (N < 0 && !betterNegativeBands) continue;
            if (N == 0) continue;
            target_CLs=2*(1-ROOT::Math::gaussian_cdf(fabs(N))); // change this so findCrossing looks for sqrt(qmu95)=2
            if (N < 0) direction = -1;

            //get the acual value
            double NtimesSigma = getLimit(asimov_0_nll, N*med_limit/sqrt(3.84)); // use N * sigma(0) as an initial guess
            N_status[N] += global_status;
            sigma = NtimesSigma/N;
            WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "");
            WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Found N * sigma = " + std::to_string(N) + " * " + std::to_string(sigma));

            string muStr,muStrPr;
            w->loadSnapshot("conditionalGlobs_0");
            double pr_val = NtimesSigma;
            if (N < 0 && profileNegativeAtZero) pr_val = 0;
            RooDataSet* asimovData_N = makeAsimovData(1, asimov_0_nll, NtimesSigma, &muStr, &muStrPr, pr_val, 0);
            //RooDataSet* asimovData_N = makeAsimovData2(asimov_0_nll, NtimesSigma, pr_val, &muStr, &muStrPr);


            RooNLLVar* asimov_N_nll = createNLL(asimovData_N);//(RooNLLVar*)pdf->createNLL(*asimovData_N);
            map_data_nll[asimovData_N] = asimov_N_nll;
            map_snapshots[asimov_N_nll] = "conditionalGlobs"+muStrPr;
            w->loadSnapshot(map_snapshots[asimov_N_nll].c_str());
            w->loadSnapshot(("conditionalNuis"+muStrPr).c_str());
            setMu(NtimesSigma);

            double nll_val = asimov_N_nll->getVal();
            saveSnapshot(asimov_N_nll, NtimesSigma);
            map_muhat[asimov_N_nll] = NtimesSigma;
            if (N < 0 && doTilde){
                setMu(0);
                firstPOI->setConstant(1);
                nll_val = getNLL(asimov_N_nll);
            }
            map_nll_muhat[asimov_N_nll] = nll_val;

            target_CLs = init_targetCLs;
            direction=1;
            double initial_guess = findCrossing(NtimesSigma/N, NtimesSigma/N, NtimesSigma);
            double limit = getLimit(asimov_N_nll, initial_guess);
            N_status[N] += global_status;

            if (N == 2) mu_up_p2 = limit;
            else if (N == 1) mu_up_p1 = limit;
            else if (N ==-1) mu_up_n1 = limit;
            else if (N ==-2) mu_up_n2 = limit;
            //return;
        }
        direction = 1;
        target_CLs = init_targetCLs;

    }

    w->loadSnapshot("conditionalNuis_0");
    firstPOI->setRange(firstPOIMin, firstPOIMax);
    double obs_limit = doObs ? getLimit(obs_nll, med_limit) : 0;
    int obs_status=global_status;

    bool hasFailures = false;
    if (obs_status != 0 || med_status != 0 || asimov0_status != 0) hasFailures = true;
    for (map<int, int>::iterator itr=N_status.begin();itr!=N_status.end();itr++){
        if (itr->second != 0) hasFailures = true;
    }
    if (hasFailures){
        WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "--------------------------------");
        WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Unresolved fit failures detected");
        WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Asimov0:  " + std::to_string(asimov0_status));
        for (map<int, int>::iterator itr=N_status.begin();itr!=N_status.end();itr++){
            WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "+" + std::to_string(itr->first) + "sigma:    " +  std::to_string(itr->first));
        }
        WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Median:    " + std::to_string(med_status));
        WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Observed:  " + std::to_string(obs_status));
        WriteWarningStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "--------------------------------");
    }

    if (betterBands) WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Guess for bands");
    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "+2sigma:  " + std::to_string(mu_up_p2_approx));
    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "+1sigma:  " + std::to_string(mu_up_p1_approx));
    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "-1sigma:  " + std::to_string(mu_up_n1_approx));
    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "-2sigma:  " + std::to_string(mu_up_n2_approx));
    if (betterBands){
        WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "");
        WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Correct bands");
        WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "+2sigma:  " + std::to_string(mu_up_p2));
        WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "+1sigma:  " + std::to_string(mu_up_p1));
        WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "-1sigma:  " + std::to_string(mu_up_n1));
        WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "-2sigma:  " + std::to_string(mu_up_n2));

    }

    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Median:    " + std::to_string(med_limit));
    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Observed:  " + std::to_string(obs_limit));
    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "");


//   system(("mkdir -vp root-files/" + folder).c_str());
    system(("mkdir -vp " + folder).c_str());

    stringstream fileName;
//   fileName << "root-files/" << folder << "/" << mass << ".root";
    fileName << folder << "/" << mass << ".root";
    TFile fout(fileName.str().c_str(),"recreate");

    TH1D* h_lim = new TH1D("limit","limit",7,0,7);
    h_lim->SetBinContent(1, obs_limit);
    h_lim->SetBinContent(2, med_limit);
    h_lim->SetBinContent(3, mu_up_p2);
    h_lim->SetBinContent(4, mu_up_p1);
    h_lim->SetBinContent(5, mu_up_n1);
    h_lim->SetBinContent(6, mu_up_n2);
    h_lim->SetBinContent(7, global_status);

    h_lim->GetXaxis()->SetBinLabel(1, "Observed");
    h_lim->GetXaxis()->SetBinLabel(2, "Expected");
    h_lim->GetXaxis()->SetBinLabel(3, "+2sigma");
    h_lim->GetXaxis()->SetBinLabel(4, "+1sigma");
    h_lim->GetXaxis()->SetBinLabel(5, "-1sigma");
    h_lim->GetXaxis()->SetBinLabel(6, "-2sigma");
    h_lim->GetXaxis()->SetBinLabel(7, "Global status"); // do something with this later

    TH1D* h_lim_old = new TH1D("limit_old","limit_old",7,0,7); // include also old approximation of bands
    h_lim_old->SetBinContent(1, obs_limit);
    h_lim_old->SetBinContent(2, med_limit);
    h_lim_old->SetBinContent(3, mu_up_p2_approx);
    h_lim_old->SetBinContent(4, mu_up_p1_approx);
    h_lim_old->SetBinContent(5, mu_up_n1_approx);
    h_lim_old->SetBinContent(6, mu_up_n2_approx);
    h_lim_old->SetBinContent(7, global_status);

    h_lim_old->GetXaxis()->SetBinLabel(1, "Observed");
    h_lim_old->GetXaxis()->SetBinLabel(2, "Expected");
    h_lim_old->GetXaxis()->SetBinLabel(3, "+2sigma");
    h_lim_old->GetXaxis()->SetBinLabel(4, "+1sigma");
    h_lim_old->GetXaxis()->SetBinLabel(5, "-1sigma");
    h_lim_old->GetXaxis()->SetBinLabel(6, "-2sigma");
    h_lim_old->GetXaxis()->SetBinLabel(7, "Global status");

    fout.Write();
    fout.Close();

    WriteInfoStatus("RunAsymptoticsCLs::RunAsymptoticsCLs", "Finished with " + std::to_string(nrMinimize) + " calls to minimize(nll)");
    timer.Print();
}

double LimitsCLs::getLimit(RooNLLVar* nll, double initial_guess){
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "------------------------");
    std::string temp_nll_s = nll->GetName();
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "Getting limit for nll: " + temp_nll_s);
    //get initial guess based on muhat and sigma(muhat)
    firstPOI->setConstant(0);
    global_status=0;

    if (nll == asimov_0_nll) {
        setMu(0);
        firstPOI->setConstant(1);
    }

    double muhat;
    if (map_nll_muhat.find(nll) == map_nll_muhat.end()){
        double nll_val = getNLL(nll);
        muhat = firstPOI->getVal();
        saveSnapshot(nll, muhat);
        map_muhat[nll] = muhat;
        if (muhat < 0 && doTilde){
            setMu(0);
            firstPOI->setConstant(1);
            nll_val = getNLL(nll);
        }

        map_nll_muhat[nll] = nll_val;
    }
    else{
        muhat = map_muhat[nll];
    }

    if (muhat < 0.1 || initial_guess != 0) setMu(initial_guess);
    double qmu,qmuA;
    double sigma_guess = getSigma(asimov_0_nll, firstPOI->getVal(), 0, qmu);
    double sigma_b = sigma_guess;
    double mu_guess = findCrossing(sigma_guess, sigma_b, muhat);
    double pmu = calcPmu(qmu, sigma_b, mu_guess);
    double pb = calcPb(qmu, sigma_b, mu_guess);
    double CLs = calcCLs(qmu, sigma_b, mu_guess);
    double qmu95 = getQmu95(sigma_b, mu_guess);
    setMu(mu_guess);

    WriteInfoStatus("RunAsymptoticCLs::getLimit", "Initial guess:  " + std::to_string(mu_guess));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "Sigma(obs):     " + std::to_string(sigma_guess));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "Sigma(mu,0):    " + std::to_string(sigma_b));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "muhat:          " + std::to_string(muhat));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "pmu:            " + std::to_string(pmu));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "qmu95:          " + std::to_string(qmu95));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "qmu:            " + std::to_string(qmu));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "1-pb:           " + std::to_string(pb));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "CLs:            " + std::to_string(CLs));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "");

    int nrDamping = 1;
    map<double, double> guess_to_corr;
    double damping_factor = 1.0;
    //double damping_factor_pre = damping_factor;
    int nrItr = 0;
    double mu_pre = muhat;//mu_guess-10*precision*mu_guess;
    double mu_pre2 = muhat;
    while (fabs(mu_pre-mu_guess) > precision*mu_guess*direction){
        std::string tmp_s = nll->GetName();
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "----------------------");
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "Starting iteration " + std::to_string(nrItr) + " of " + tmp_s);
        // do this to avoid comparing multiple minima in the conditional and unconditional fits
        if (nrItr == 0) loadSnapshot(nll, muhat);
        else if (usePredictiveFit) doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
        else loadSnapshot(asimov_0_nll, mu_pre);

        sigma_guess=getSigma(nll, mu_guess, muhat, qmu);
        saveSnapshot(nll, mu_guess);


        if (nll != asimov_0_nll){
            if (nrItr == 0) loadSnapshot(asimov_0_nll, map_nll_muhat[asimov_0_nll]);
            else if (usePredictiveFit) {
                if (nrItr == 1) doPredictiveFit(nll, map_nll_muhat[asimov_0_nll], mu_pre, mu_guess);
                else doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
            }
            else loadSnapshot(asimov_0_nll, mu_pre);

            sigma_b=getSigma(asimov_0_nll, mu_guess, 0, qmuA);
            saveSnapshot(asimov_0_nll, mu_guess);
        }
        else{
            sigma_b=sigma_guess;
            qmuA=qmu;
        }

        double corr = damping_factor*(mu_guess - findCrossing(sigma_guess, sigma_b, muhat));
        for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++){
            if (fabs(itr->first - (mu_guess-corr)) < direction*mu_guess*0.02 && fabs(corr) > direction*mu_guess*precision) {
                damping_factor *= 0.8;
                WriteInfoStatus("RunAsymptoticCLs::getLimit", "Changing damping factor to " +std::to_string( damping_factor) + ", nrDamping = " + std::to_string(nrDamping));
                if (nrDamping++ > 10){
                    nrDamping = 1;
                    damping_factor = 1.0;
                }
                corr *= damping_factor;
                break;
            }
        }

        //subtract off the difference in the new and damped correction
        guess_to_corr[mu_guess] = corr;
        mu_pre2 = mu_pre;
        mu_pre = mu_guess;
        mu_guess -= corr;


        pmu = calcPmu(qmu, sigma_b, mu_pre);
        pb = calcPb(qmu, sigma_b, mu_pre);
        CLs = calcCLs(qmu, sigma_b, mu_pre);
        qmu95 = getQmu95(sigma_b, mu_pre);


        tmp_s = nll->GetName();
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "NLL             " + tmp_s);
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "Previous guess  " + std::to_string(mu_pre));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "Sigma(obs):     " + std::to_string(sigma_guess));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "Sigma(mu,0):    " + std::to_string(sigma_b));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "muhat:          " + std::to_string(muhat));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "pmu:            " + std::to_string(pmu));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "1-pb:           " + std::to_string(pb));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "CLs:            " + std::to_string(CLs));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "qmu95:          " + std::to_string(qmu95));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "qmu:            " + std::to_string(qmu));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "qmuA0:          " + std::to_string(qmuA));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "Precision:      " + std::to_string(direction*mu_guess*precision));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "Correction:     " + std::to_string(-corr));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "New guess:      " + std::to_string(mu_guess));
        WriteInfoStatus("RunAsymptoticCLs::getLimit", "");

        nrItr++;
        if (nrItr > 25){
            WriteErrorStatus("RunAsymptoticCLs::getLimit", "Infinite loop detected in getLimit(). Please intervene.");
            break;
        }
    }

    std::string tmp_s = nll->GetName();
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "Found limit for nll " + tmp_s + ": " + std::to_string(mu_guess));
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "Finished in " + std::to_string(nrItr) + " iterations. ");
    WriteInfoStatus("RunAsymptoticCLs::getLimit", "");
    return mu_guess;
}


double LimitsCLs::getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu){
    qmu = getQmu(nll, mu);
    WriteDebugStatus("RunAsymptoticCLs::getLimit", "qmu= " + std::to_string(qmu));
    if (mu*direction < muhat) return fabs(mu-muhat)/sqrt(qmu);
    else if (muhat < 0 && doTilde) return sqrt(mu*mu-2*mu*muhat*direction)/sqrt(qmu);
    else return (mu-muhat)*direction/sqrt(qmu);
}

double LimitsCLs::getQmu(RooNLLVar* nll, double mu){
    double nll_muhat = map_nll_muhat[nll];
    bool isConst = firstPOI->isConstant();
    firstPOI->setConstant(1);
    setMu(mu);
    double nll_val = getNLL(nll);
    firstPOI->setConstant(isConst);
    return 2*(nll_val-nll_muhat);
}

void LimitsCLs::saveSnapshot(RooNLLVar* nll, double mu){
    stringstream snapshotName;
    snapshotName << nll->GetName() << "_" << mu;
    w->saveSnapshot(snapshotName.str().c_str(), *mc->GetNuisanceParameters());
}

void LimitsCLs::loadSnapshot(RooNLLVar* nll, double mu){
    stringstream snapshotName;
    snapshotName << nll->GetName() << "_" << mu;
    w->loadSnapshot(snapshotName.str().c_str());
}

void LimitsCLs::doPredictiveFit(RooNLLVar* nll, double mu1, double mu2, double mu){
    if (fabs(mu2-mu) < direction*mu*precision*4){
        loadSnapshot(nll, mu2);
        return;
    }

//extrapolate to mu using mu1 and mu2 assuming nuis scale linear in mu
    const RooArgSet* nuis = mc->GetNuisanceParameters();
    int nrNuis = nuis->getSize();
    double* theta_mu1 = new double[nrNuis];
    double* theta_mu2 = new double[nrNuis];

    TIterator* itr = nuis->createIterator();
    RooRealVar* var = 0;
    int counter = 0;
    loadSnapshot(nll, mu1);
    while ((var = (RooRealVar*)itr->Next())){
        theta_mu1[counter++] = var->getVal();
    }

    itr->Reset();
    counter = 0;
    loadSnapshot(nll, mu2);
    while ((var = (RooRealVar*)itr->Next())){
        theta_mu2[counter++] = var->getVal();
    }

    itr->Reset();
    counter = 0;
    while ((var = (RooRealVar*)itr->Next()))
    {
        double m = (theta_mu2[counter] - theta_mu1[counter])/(mu2-mu1);
        double b = theta_mu2[counter] - m*mu2;
        double theta_extrap = m*mu+b;

        var->setVal(theta_extrap);
        counter++;
    }

    delete itr;
    delete[] theta_mu1;
    delete[] theta_mu2;
}

RooNLLVar* LimitsCLs::createNLL(RooDataSet* _data){
    RooArgSet nuis = *mc->GetNuisanceParameters();
//   RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*_data, Constrain(nuis));
//   RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*_data, Constrain(nuis),RooFit::NumCPU(numCPU,RooFit::Hybrid) );
    RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*_data, Constrain(nuis),RooFit::Offset(1),RooFit::NumCPU(numCPU,RooFit::Hybrid) );
    return nll;
}

double LimitsCLs::getNLL(RooNLLVar* nll){
    string snapshotName = map_snapshots[nll];
    if (snapshotName != "") w->loadSnapshot(snapshotName.c_str());
    minimize(nll);
    double val = nll->getVal();
    w->loadSnapshot("nominalGlobs");
    return val;
}


double LimitsCLs::findCrossing(double sigma_obs, double sigma, double muhat){
    double mu_guess = muhat + ROOT::Math::gaussian_quantile(1-target_CLs,1)*sigma_obs*direction;
    int nrItr = 0;
    int nrDamping = 1;

    map<double, double> guess_to_corr;
    double damping_factor = 1.0;
    double mu_pre = mu_guess - 10*mu_guess*precision;
    while (fabs(mu_guess-mu_pre) > direction*mu_guess*precision){
        double sigma_obs_extrap = sigma_obs;
        double eps = 0;
        mu_pre = mu_guess;

        double qmu95 = getQmu95(sigma, mu_guess);
        double qmu;
        qmu = 1./sigma_obs_extrap/sigma_obs_extrap*(mu_guess-muhat)*(mu_guess-muhat);
        if (muhat < 0 && doTilde) qmu = 1./sigma_obs_extrap/sigma_obs_extrap*(mu_guess*mu_guess-2*mu_guess*muhat);

        double dqmu_dmu = 2*(mu_guess-muhat)/sigma_obs_extrap/sigma_obs_extrap - 2*qmu*eps;

        double corr = damping_factor*(qmu-qmu95)/dqmu_dmu;
        for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++){
            if (fabs(itr->first - mu_guess) < direction*mu_guess*precision) {
                damping_factor *= 0.8;
                WriteDebugStatus("RunAsymptoticCLs::findCrossing", "Changing damping factor to " + std::to_string(damping_factor) + ", nrDamping = " + std::to_string(nrDamping));
                if (nrDamping++ > 10){
                    nrDamping = 1;
                    damping_factor = 1.0;
                }
                corr *= damping_factor;
                break;
            }
        }
        guess_to_corr[mu_guess] = corr;

        mu_guess = mu_guess - corr;
        nrItr++;
        if (nrItr > 100){
            WriteErrorStatus("RunAsymptoticCLs::findCrossing", "Infinite loop detected in findCrossing. Please intervene.");
            exit(1);
        }
        WriteDebugStatus("RunAsymptoticCLs::findCrossing", "mu_guess = " + std::to_string(mu_guess) + ", mu_pre = " + std::to_string(mu_pre) + ", qmu = " + std::to_string(qmu) + ", qmu95 = " + std::to_string(qmu95) + ", sigma_obs_extrap = " + std::to_string(sigma_obs_extrap) + ", sigma = " + std::to_string(sigma) + ", direction*mu*prec = " + std::to_string(direction*mu_guess*precision));
    }

    return mu_guess;
}

void LimitsCLs::setMu(double mu){
    if (mu != mu){
        WriteErrorStatus("RunAsymptoticCLs::setMu", "POI gave nan. Please intervene.");
        exit(1);
    }
    if (mu > 0 && firstPOI->getMax() < mu) firstPOI->setMax(2*mu);
    if (mu < 0 && firstPOI->getMin() > mu) firstPOI->setMin(2*mu);
    firstPOI->setVal(mu);
}


double LimitsCLs::getQmu95_brute(double sigma, double mu){
    double step_size = 0.001;
    double start = step_size;
    if (mu/sigma > 0.2) start = 0;
    for (double qmu=start;qmu<20;qmu+=step_size){
        double CLs = calcCLs(qmu, sigma, mu);

        if (CLs < target_CLs) return qmu;
    }

    return 20;
}

double LimitsCLs::getQmu95(double sigma, double mu){
    double qmu95 = 0;
    //no sane man would venture this far down into |mu/sigma|
    double target_N = ROOT::Math::gaussian_cdf(1-target_CLs,1);
    if (fabs(mu/sigma) < 0.25*target_N){
        qmu95 = 5.83/target_N;
    }
    else{
        map<double, double> guess_to_corr;
        double qmu95_guess = pow(ROOT::Math::gaussian_quantile(1-target_CLs,1),2);
        int nrItr = 0;
        int nrDamping = 1;
        double damping_factor = 1.0;
        double qmu95_pre = qmu95_guess - 10*2*qmu95_guess*precision;
        while (fabs(qmu95_guess-qmu95_pre) > 2*qmu95_guess*precision){
            qmu95_pre = qmu95_guess;
            if (TtHFitter::DEBUGLEVEL > 2){
                WriteVerboseStatus("RunAsymptoticCLs::getQmu95", "qmu95_guess = " + std::to_string(qmu95_guess));
                WriteVerboseStatus("RunAsymptoticCLs::getQmu95", "CLs =         " + std::to_string(calcCLs(qmu95_guess, sigma, mu)));
                WriteVerboseStatus("RunAsymptoticCLs::getQmu95", "Derivative =  " + std::to_string(calcDerCLs(qmu95_guess, sigma, mu)));
            }

            double corr = damping_factor*(calcCLs(qmu95_guess, sigma, mu)-target_CLs)/calcDerCLs(qmu95_guess, sigma, mu);
            for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++){
                if (fabs(itr->first - qmu95_guess) < 2*qmu95_guess*precision) {
                    damping_factor *= 0.8;
                    if (TtHFitter::DEBUGLEVEL > 2) WriteVerboseStatus("RunAsymptoticCLs::getQmu95", "Changing damping factor to " + std::to_string(damping_factor) + ", nrDamping = " + std::to_string(nrDamping));
                    if (nrDamping++ > 10){
                        nrDamping = 1;
                        damping_factor = 1.0;
                    }
                    corr *= damping_factor;
                }
            }

            guess_to_corr[qmu95_guess] = corr;
            qmu95_guess = qmu95_guess - corr;

            if (TtHFitter::DEBUGLEVEL > 2) {
                WriteVerboseStatus("RunAsymptoticCLs::getQmu95", "next guess = " + std::to_string(qmu95_guess));
                WriteVerboseStatus("RunAsymptoticCLs::getQmu95", "precision =  " + std::to_string(2*qmu95_guess*precision));
            }
            nrItr++;
            if (nrItr > 200){
                WriteErrorStatus("RunAsymptoticCLs::getQmu95", "Infinite loop detected in getQmu95. Please intervene.");
                exit(1);
            }
        }
        qmu95 = qmu95_guess;
    }

    if (qmu95 != qmu95) {
        qmu95 = getQmu95_brute(sigma, mu);
    }
    WriteDebugStatus("RunAsymptoticCLs::getQmu95", "Returning qmu95 = " + std::to_string(qmu95));

    return qmu95;
}

double LimitsCLs::calcCLs(double qmu_tilde, double sigma, double mu){
    double pmu = calcPmu(qmu_tilde, sigma, mu);
    double pb = calcPb(qmu_tilde, sigma, mu);
    if (TtHFitter::DEBUGLEVEL > 2){ // need to explicitelly do this otherwise the code is super slow
        WriteVerboseStatus("RunAsymptoticCLs::calcCLs", "pmu = " + std::to_string(pmu));
        WriteVerboseStatus("RunAsymptoticCLs::calcCLs", "pb =  " + std::to_string(pb));
    }
    if (pb == 1) return 0.5;
    return pmu/(1-pb);
}

double LimitsCLs::calcPmu(double qmu, double sigma, double mu){
    double pmu;
    if (qmu < mu*mu/(sigma*sigma) || !doTilde){
        pmu = 1-ROOT::Math::gaussian_cdf(sqrt(qmu));
    }
    else{
        pmu = 1-ROOT::Math::gaussian_cdf((qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma)));
    }
    if (TtHFitter::DEBUGLEVEL > 2) WriteVerboseStatus("RunAsymptoticCLs::calcPmu", "for pmu, qmu = " + std::to_string(qmu) + ", sigma = " + std::to_string(sigma) + ", mu = " + std::to_string(mu) + ", pmu = " + std::to_string(pmu));
    return pmu;
}

double LimitsCLs::calcPb(double qmu, double sigma, double mu){
    if (qmu < mu*mu/(sigma*sigma) || !doTilde){
        return 1-ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qmu));
    }
    else{
        return 1-ROOT::Math::gaussian_cdf((mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma)));
    }
}

double LimitsCLs::calcDerCLs(double qmu, double sigma, double mu){
    double dpmu_dq = 0;
    double d1mpb_dq = 0;

    if (qmu < mu*mu/(sigma*sigma)){
        double zmu = sqrt(qmu);
        dpmu_dq = -1./(2*sqrt(qmu*2*TMath::Pi()))*exp(-zmu*zmu/2);
    }
    else {
        double zmu = (qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma));
        dpmu_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zmu*zmu/2);
    }

    if (qmu < mu*mu/(sigma*sigma)){
        double zb = fabs(mu/sigma)-sqrt(qmu);
        d1mpb_dq = -1./sqrt(qmu*2*TMath::Pi())*exp(-zb*zb/2);
    }
    else{
        double zb = (mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma));
        d1mpb_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zb*zb/2);
    }

    double pb = calcPb(qmu, sigma, mu);
    return dpmu_dq/(1-pb)-calcCLs(qmu, sigma, mu)/(1-pb)*d1mpb_dq;
}

int LimitsCLs::minimize(RooNLLVar* nll){
    nrMinimize++;
    RooAbsReal* fcn = (RooAbsReal*)nll;
    return minimize(fcn);
}

int LimitsCLs::minimize(RooAbsReal* fcn){
    static int nrItr = 0;
     // mc->GetGlobalObservables()->Print("v");


    int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    int save_strat = strat;
    RooMinimizer minim(*fcn);
    minim.setStrategy(strat);
    minim.setPrintLevel(printLevel);


    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


//up the strategy
    if (status != 0 && status != 1 && strat < 2){
        strat++;
        WriteWarningStatus("RunAsymptoticCLs::minimize", "Fit failed with status " + std::to_string(status) + ". Retrying with strategy " + std::to_string(strat));
        minim.setStrategy(strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1 && strat < 2){
        strat++;
        WriteWarningStatus("RunAsymptoticCLs::minimize", "Fit failed with status " + std::to_string(status) + ". Retrying with strategy " + std::to_string(strat));
        minim.setStrategy(strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }


// //switch minuit version and try again
    if (status != 0 && status != 1){
        string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
        string newMinType;
        if (minType == "Minuit2") newMinType = "Minuit";
        else newMinType = "Minuit2";

        WriteDebugStatus("RunAsymptoticCLs::minimize", "Switching minuit type from " + minType + " to " + newMinType);

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
        strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
        minim.setStrategy(strat);

        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


        if (status != 0 && status != 1 && strat < 2){
            strat++;
            WriteWarningStatus("RunAsymptoticCLs::minimize", "Fit failed with status " + std::to_string(status) + ". Retrying with strategy " + std::to_string(strat));
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        }

        if (status != 0 && status != 1 && strat < 2){
            strat++;
            WriteWarningStatus("RunAsymptoticCLs::minimize", "Fit failed with status " + std::to_string(status) + ". Retrying with strategy " + std::to_string(strat));
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        }

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
    }

    if (status != 0 && status != 1){
        nrItr++;
        if (nrItr > maxRetries){
            nrItr = 0;
            global_status++;
            WriteWarningStatus("RunAsymptoticCLs::minimize", "Fit failure unresolved with status " + std::to_string(status));
            return status;
        }
        else{
            if (nrItr == 0){ // retry with mu=0 snapshot
                w->loadSnapshot("conditionalNuis_0");
                return minimize(fcn);
            }
            else if (nrItr == 1){ // retry with nominal snapshot

                w->loadSnapshot("nominalNuis");
                return minimize(fcn);
            }
        }
    }

    if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);


    if (nrItr != 0) WriteInfoStatus("RunAsymptoticCLs::minimize", "Successful fit");
    nrItr=0;
    return status;
}


/*
RooDataSet* makeAsimovData2(RooDataSet* conditioningData, double mu_true, double mu_prof, string* mu_str, string* mu_prof_str)
{
    RooNLLVar* conditioningNLL = NULL;
    if (conditioningData)
    {
        conditioningNLL = (RooNLLVar*)mc->GetPdf()->createNLL(*conditioningData);
    }
    return makeAsimovData2(conditioningNLL, mu_true, mu_prof, mu_str, mu_prof_str);
}


RooDataSet* makeAsimovData2(RooNLLVar* conditioningNLL, double mu_true, double mu_prof, string* mu_str, string* mu_prof_str)
{
    if (mu_prof == -999) mu_prof = mu_true;
    bool doTest = 0;

    int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();

    int test = 0;

    int _printLevel = 0;

    stringstream muStr;
    muStr << setprecision(5);
    muStr << "_" << mu_true;
    if (mu_str) *mu_str = muStr.str();

    stringstream muStrProf;
    muStrProf << setprecision(5);
    muStrProf << "_" << mu_prof;
    if (mu_prof_str) *mu_prof_str = muStrProf.str();


    const RooArgSet* globs = mc->GetGlobalObservables();
    const RooArgSet* nuis = mc->GetNuisanceParameters();
    const RooArgSet* obs = mc->GetObservables();
    const RooArgSet* pois = mc->GetParametersOfInterest();

    TIterator* gItr = globs->createIterator();
    TIterator* nItr = nuis->createIterator();
    RooRealVar* var;

    RooArgSet emptySet;
    RooArgSet params(*nuis);
    params.add(*globs);
    params.add(*pois);
    w->saveSnapshot("initial_params", params);


//condition the MLEs
    if (conditioningNLL)
    {
        //get the conditional MLEs
        firstPOI->setVal(mu_prof);
        firstPOI->setConstant(1);
        minimize(conditioningNLL);
    }

    w->saveSnapshot(("conditionalNuis" +muStrProf.str()).c_str(),*nuis);


//to find the conditional globs, do a fit to the constraint only pdf with the globs floating and the MLEs constant
    RooArgSet obsCopy = *obs;
    RooArgSet nuisCopy = *nuis;

    RooArgSet constraints(*mc->GetPdf()->getAllConstraints(obsCopy, nuisCopy));
    RooRealVar minusOne("minusOne","minusOne",-1);
    constraints.add(minusOne);
    RooProduct constrFunc("constrFunc","constrFunc",constraints);

    while ((var = (RooRealVar*)gItr->Next()))
    {
        var->setConstant(false);
    }
    gItr->Reset();

    while ((var = (RooRealVar*)nItr->Next()))
    {
        var->setConstant(true);
    }
    nItr->Reset();

    minimize(&constrFunc);

    while ((var = (RooRealVar*)gItr->Next()))
    {
        var->setConstant(true);
    }
    gItr->Reset();

    while ((var = (RooRealVar*)nItr->Next()))
    {
        var->setConstant(false);
    }
    nItr->Reset();

    w->saveSnapshot(("conditionalGlobs"+muStrProf.str()).c_str(),*globs);





//make the asimov data
    const char* weightName="weightVar";
    RooArgSet obsAndWeight;
    obsAndWeight.add(*mc->GetObservables());

    RooRealVar* weightVar = NULL;
    if (!(weightVar = w->var(weightName)))
    {
        w->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
        weightVar = w->var(weightName);
    }
    obsAndWeight.add(*w->var(weightName));



    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
    map<string, RooDataSet*> asimovDataMap;

    //try fix for sim pdf
    RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;
    int nrIndices = 0;
    int iFrame=0;
    while((tt=(RooCatType*) iter->Next())) {
        nrIndices++;
    }
    for (int i=0;i<nrIndices;i++){
        channelCat->setIndex(i);
        iFrame++;
        // Get pdf associated with state from simpdf
        RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;

        // Generate observables defined by the pdf associated with this state
        RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

        if (_printLevel >= 1)
        {
            obstmp->Print();
        }

        RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
        RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
        double expectedEvents = pdftmp->expectedEvents(*obstmp);
        double thisNorm = 0;
        for(int jj=0; jj<thisObs->numBins(); ++jj){
            thisObs->setBin(jj);

            thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
            if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
        }

        if (_printLevel >= 1)
        {
            obsDataUnbinned->Print();
        }
        if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
            exit(1);
        }


        asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;

        if (_printLevel >= 1)
        {
            obsDataUnbinned->Print();
        }
    }

    RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
    if (w->data(("asimovData"+muStr.str()).c_str()))
    {
        w->import(*asimovData, true);
    }
    else
    {
        w->import(*asimovData);
    }




    w->loadSnapshot("initial_params");
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

    return asimovData;

}
*/

void LimitsCLs::unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter){
    if (counter > 50){
        WriteErrorStatus("RunAsymptoticCLs::unfoldConstraints", "Couldn't unfold constraints!");
        WriteErrorStatus("RunAsymptoticCLs::unfoldConstraints", "Initial");
        initial.Print("v");
        WriteErrorStatus("RunAsymptoticCLs::unfoldConstraints", "");
        WriteErrorStatus("RunAsymptoticCLs::unfoldConstraints", "Final: ");
        final.Print("v");
        exit(1);
    }
    TIterator* itr = initial.createIterator();
    RooAbsPdf* pdf;
    while ((pdf = (RooAbsPdf*)itr->Next())){
        RooArgSet nuis_tmp = nuis;
        RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
        //if (constraint_set.getSize() > 1)
        //{
        string className(pdf->ClassName());
        if (className != "RooGaussian" && className != "RooLognormal" && className != "RooGamma" && className != "RooPoisson" && className != "RooBifurGauss"){
            counter++;
            unfoldConstraints(constraint_set, final, obs, nuis, counter);
        }
        else{
            final.add(*pdf);
        }
    }
    delete itr;
}



RooDataSet* LimitsCLs::makeAsimovData(bool doConditional, RooNLLVar* conditioning_nll, double mu_val, string* mu_str, string* mu_prof_str, double mu_val_profile, bool doFit){
    if (mu_val_profile == -999) mu_val_profile = mu_val;


    WriteInfoStatus("RunAsymptoticCLs::makeAsimovData", "Creating asimov data at mu = " + std::to_string(mu_val) + ", profiling at mu = " +std::to_string( mu_val_profile));

    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    //int strat = ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
    //int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    //RooMinuit::SetMaxIterations(10000);
    //RooMinimizer::SetMaxFunctionCalls(10000);

////////////////////
//make asimov data//
////////////////////
    RooAbsPdf* combPdf = mc->GetPdf();

    stringstream muStr;
    muStr << setprecision(5);
    muStr << "_" << mu_val;
    if (mu_str) *mu_str = muStr.str();

    stringstream muStrProf;
    muStrProf << setprecision(5);
    muStrProf << "_" << mu_val_profile;
    if (mu_prof_str) *mu_prof_str = muStrProf.str();

    RooRealVar* mu = (RooRealVar*)mc->GetParametersOfInterest()->first();//w->var("mu");
    mu->setVal(mu_val);

    RooArgSet mc_obs = *mc->GetObservables();
    RooArgSet mc_globs = *mc->GetGlobalObservables();
    RooArgSet mc_nuis = *mc->GetNuisanceParameters();

//pair the nuisance parameter to the global observable
    RooArgSet mc_nuis_tmp = mc_nuis;
    RooArgList nui_list("ordered_nuis");
    RooArgList glob_list("ordered_globs");
    RooArgSet constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
    RooArgSet constraint_set;
    int counter_tmp = 0;
    unfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);

    TIterator* cIter = constraint_set.createIterator();
    RooAbsArg* arg;
    while ((arg = (RooAbsArg*)cIter->Next())){
        RooAbsPdf* pdf = (RooAbsPdf*)arg;
        if (!pdf) continue;
//       pdf->Print();
        TIterator* nIter = mc_nuis.createIterator();
        RooRealVar* thisNui = NULL;
        RooAbsArg* nui_arg;
        while ((nui_arg = (RooAbsArg*)nIter->Next())){
            if (pdf->dependsOn(*nui_arg)){
                thisNui = (RooRealVar*)nui_arg;
                break;
            }
        }
        delete nIter;

        //RooRealVar* thisNui = (RooRealVar*)pdf->getObservables();


//need this incase the observable isn't fundamental.
//in this case, see which variable is dependent on the nuisance parameter and use that.
        RooArgSet* components = pdf->getComponents();
//       components->Print();
        components->remove(*pdf);
        if (components->getSize()){
            TIterator* itr1 = components->createIterator();
            RooAbsArg* arg1;
            while ((arg1 = (RooAbsArg*)itr1->Next())){
                TIterator* itr2 = components->createIterator();
                RooAbsArg* arg2;
                while ((arg2 = (RooAbsArg*)itr2->Next())){
                    if (arg1 == arg2) continue;
                    if (arg2->dependsOn(*arg1)){
                        components->remove(*arg1);
                    }
                }
                delete itr2;
            }
            delete itr1;
        }
        if (components->getSize() > 1){
            WriteErrorStatus("RunAsymptoticCLs::makeAsimovData", "Couldn't isolate proper nuisance parameter");
            return NULL;
        }
        else if (components->getSize() == 1){
            thisNui = (RooRealVar*)components->first();
        }


        TIterator* gIter = mc_globs.createIterator();
        RooRealVar* thisGlob = NULL;
        RooAbsArg* glob_arg;
        while ((glob_arg = (RooAbsArg*)gIter->Next())){
            if (pdf->dependsOn(*glob_arg)){
                thisGlob = (RooRealVar*)glob_arg;
                break;
            }
        }
        delete gIter;

        if (!thisNui || !thisGlob){
            std::string s = pdf->GetName();
            WriteWarningStatus("RunAsymptoticCLs::makeAsimovData", "Couldn't find nui or glob for constraint: " + s);
            //return;
            continue;
        }

        std::string s1 = thisNui->GetName();
        std::string s2 = thisGlob->GetName();
        std::string s3 = pdf->GetName();;
        WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "Pairing nui: " + s1 + ", with glob: " + s2 + ", from constraint: " + s3);

        nui_list.add(*thisNui);
        glob_list.add(*thisGlob);

//       thisNui->Print();
//       thisGlob->Print();
    }
    delete cIter;


//save the snapshots of nominal parameters, but only if they're not already saved
    w->saveSnapshot("tmpGlobs",*mc->GetGlobalObservables());
    w->saveSnapshot("tmpNuis",*mc->GetNuisanceParameters());
    if (!w->loadSnapshot("nominalGlobs")){
        WriteInfoStatus("RunAsymptoticCLs::makeAsimovData", "nominalGlobs doesn't exist. Saving snapshot.");
        w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
    }
    else w->loadSnapshot("tmpGlobs");
    if (!w->loadSnapshot("nominalNuis")){
        WriteInfoStatus("RunAsymptoticCLs::makeAsimovData", "nominalNuis doesn't exist. Saving snapshot.");
        w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());
    }
    else w->loadSnapshot("tmpNuis");

    RooArgSet nuiSet_tmp(nui_list);

    mu->setVal(mu_val_profile);
    mu->setConstant(1);
    //int status = 0;
    if (doConditional && doFit){
        minimize(conditioning_nll);
        // mc->GetGlobalObservables()->Print("v");
        // RooAbsReal* nll;
        // if (!(nll = map_data_nll[combData])) nll = combPdf->createNLL(*combData, RooFit::Constrain(nuiSet_tmp));
        // RooMinimizer minim(*nll);
        // minim.setStrategy(0);
        // minim.setPrintLevel(1);
        // status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        // if (status != 0)
        // {
        // }

        //combPdf->fitTo(*combData,Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp));
    }
    mu->setConstant(0);
    mu->setVal(mu_val);


//loop over the nui/glob list, grab the corresponding variable from the tmp ws, and set the glob to the value of the nui
    int nrNuis = nui_list.getSize();
    if (nrNuis != glob_list.getSize()){
       WriteErrorStatus("RunAsymptoticCLs::makeAsimovData", "nui_list.getSize() != glob_list.getSize()!");
        return NULL;
    }

    for (int i=0;i<nrNuis;i++){
        RooRealVar* nui = (RooRealVar*)nui_list.at(i);
        RooRealVar* glob = (RooRealVar*)glob_list.at(i);


        glob->setVal(nui->getVal());
    }

//save the snapshots of conditional parameters
    w->saveSnapshot(("conditionalGlobs"+muStrProf.str()).c_str(),*mc->GetGlobalObservables());
    w->saveSnapshot(("conditionalNuis" +muStrProf.str()).c_str(),*mc->GetNuisanceParameters());

    if (!doConditional){
        w->loadSnapshot("nominalGlobs");
        w->loadSnapshot("nominalNuis");
    }

    WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "Making asimov");
//make the asimov data (snipped from Kyle)
    mu->setVal(mu_val);

    int iFrame=0;

    const char* weightName="weightVar";
    RooArgSet obsAndWeight;
    obsAndWeight.add(*mc->GetObservables());

    RooRealVar* weightVar = NULL;
    if (!(weightVar = w->var(weightName))){
        w->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
        weightVar = w->var(weightName);
    }
    obsAndWeight.add(*w->var(weightName));

    w->defineSet("obsAndWeight",obsAndWeight);


    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    // MAKE ASIMOV DATA FOR OBSERVABLES

    // dummy var can just have one bin since it's a dummy
    //if(w->var("ATLAS_dummyX"))    w->var("ATLAS_dummyX")->setBins(1);

    //RooDataSet* simData=NULL;
    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());

    RooDataSet* asimovData;
    if (!simPdf){
        // Get pdf associated with state from simpdf
        RooAbsPdf* pdftmp = mc->GetPdf();//simPdf->getPdf(channelCat->getLabel()) ;

        // Generate observables defined by the pdf associated with this state
        RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

        if (TtHFitter::DEBUGLEVEL > 1){
            obstmp->Print();
        }

        asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight),WeightVar(*weightVar));

        RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
        double expectedEvents = pdftmp->expectedEvents(*obstmp);
        double thisNorm = 0;
        for(int jj=0; jj<thisObs->numBins(); ++jj){
            thisObs->setBin(jj);

            thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
            if (thisNorm*expectedEvents <= 0){
                std::string s = thisObs->GetName();
                WriteWarningStatus("RunAsymptoticCLs::makeAsimovData", "Detected bin with zero expected events (" + std::to_string(thisNorm*expectedEvents) + ") ! Please check your inputs. Obs = " + s + ", bin = " + std::to_string(jj));
            }
            if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) asimovData->add(*mc->GetObservables(), thisNorm*expectedEvents);
        }

        if (TtHFitter::DEBUGLEVEL > 1){
            asimovData->Print();
            WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "sum entries " + std::to_string(asimovData->sumEntries()));
        }
        if(asimovData->sumEntries()!=asimovData->sumEntries()){
            WriteErrorStatus("RunAsymptoticCLs::makeAsimovData", "sum entries is nan");
            exit(1);
        }

        //((RooRealVar*)obstmp->first())->Print();

        w->import(*asimovData);

        if (TtHFitter::DEBUGLEVEL > 1){
            asimovData->Print();
            WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "");
        }
    }
    else{
        map<string, RooDataSet*> asimovDataMap;


        //try fix for sim pdf
        RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();//(RooCategory*)w->cat("master_channel");//(RooCategory*) (&simPdf->indexCat());
        //      TIterator* iter = simPdf->indexCat().typeIterator() ;
        TIterator* iter = channelCat->typeIterator() ;
        RooCatType* tt = NULL;
        int nrIndices = 0;
        while((tt=(RooCatType*) iter->Next())) {
            nrIndices++;
        }
        for (int i=0;i<nrIndices;i++){
            channelCat->setIndex(i);
            iFrame++;
            // Get pdf associated with state from simpdf
            RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;

            // Generate observables defined by the pdf associated with this state
            RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

            if (TtHFitter::DEBUGLEVEL > 1){
                obstmp->Print();
                std::string s = channelCat->getLabel();
                WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "on type " + s + " ");
            }

            RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
            RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
            double expectedEvents = pdftmp->expectedEvents(*obstmp);
            double thisNorm = 0;
            for(int jj=0; jj<thisObs->numBins(); ++jj){
                thisObs->setBin(jj);

                thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
                if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
            }

            if (TtHFitter::DEBUGLEVEL > 1){
                obsDataUnbinned->Print();
                WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "sum entries " + std::to_string(obsDataUnbinned->sumEntries()));
            }
            if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
                WriteErrorStatus("RunAsymptoticCLs::makeAsimovData", "sum entries is nan");
                exit(1);
            }

            // ((RooRealVar*)obstmp->first())->Print();

            asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;

            if (TtHFitter::DEBUGLEVEL > 1){
                std::string s = channelCat->getLabel();
                WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "channel: " + s + ", data: ");
                obsDataUnbinned->Print();
                WriteDebugStatus("RunAsymptoticCLs::makeAsimovData", "");
            }
        }

        asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
        w->import(*asimovData);
    }

//bring us back to nominal for exporting
    //w->loadSnapshot("nominalNuis");
    w->loadSnapshot("nominalGlobs");

    //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

    return asimovData;
}
