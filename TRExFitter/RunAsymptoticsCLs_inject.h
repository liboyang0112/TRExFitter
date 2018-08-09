#ifndef RUNASYMPTOTICSCLS_INJECT_H
#define RUNASYMPTOTICSCLS_INJECT_H

#include <string>

class RooNLLVar;
class RooDataSet;
class RooWorkspace;
class RooRealVar;
class RooAbsReal;
class RooArgSet;

namespace LimitsCLs_inject{
    //main
    void RunAsymptoticsCLs_inject(const char* infile,
    		       const char* workspaceName,
    		       const char* modelConfigName,
    		       const char* dataName,
    		       const char* asimovDataName,
    		       std::string folder,
    		       std::string mass,
    		       double CL);

     //for backwards compatibility
    void RunAsymptoticsCLs_inject(const char* infile,
    		       const char* workspaceName = "combWS",
    		       const char* modelConfigName = "ModelConfig",
    		       const char* dataName = "combData",
    		       std::string option = "",
    		       const char* asimovDataName = "asimovData_0",
    		       const char* conditionalSnapshot = "conditionalGlobs_0",
    		       const char* nominalSnapshot = "nominalGlobs",
    		       std::string folder = "test",
    		       std::string mass = "125",
    		       double CL = 0.95);

    double getLimit(RooNLLVar* nll, double initial_guess = 0);
    double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu);
    double getQmu(RooNLLVar* nll, double mu);
    void saveSnapshot(RooNLLVar* nll, double mu);
    void loadSnapshot(RooNLLVar* nll, double mu);
    void doPredictiveFit(RooNLLVar* nll, double mu1, double m2, double mu);
    RooNLLVar* createNLL(RooDataSet* _data);
    double getNLL(RooNLLVar* nll);
    double findCrossing(double sigma_obs, double sigma, double muhat);
    void setMu(double mu);
    double getQmu95_brute(double sigma, double mu);
    double getQmu95(double sigma, double mu);
    double calcCLs(double qmu_tilde, double sigma, double mu);
    double calcPmu(double qmu_tilde, double sigma, double mu);
    double calcPb(double qmu_tilde, double sigma, double mu);
    double calcDerCLs(double qmu, double sigma, double mu);
    int minimize(RooNLLVar* nll);
    int minimize(RooAbsReal* nll);

    void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
    RooDataSet* makeAsimovData(bool doConditional, RooNLLVar* conditioning_nll, double mu_val, std::string* mu_str = NULL, std::string* mu_prof_str = NULL, double mu_val_profile = -999, bool doFit = true);

}

#endif
