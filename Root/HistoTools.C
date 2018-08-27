/*

 HistoTools
 =========

 Contains all functions needed to handle properly the histograms:
    -> Variable Binning handling for fit/limits
    -> Symmetrisation of systematics
    -> Smoothing of systematics

 Call of the functions:
    -> #include "TRExFitter/HistoTools.C"
    -> Call of the function with HistoTools:: (Helps readability)

 Contact: Loic Valery <lvalery@cern.ch>

 */

// Class include
#include "TRExFitter/HistoTools.h"

// framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/SystematicHist.h"

// CommonStatTools includes
#include "CommonSystSmoothingTool/SmoothSystematics/SmoothHist.h"

// ROOT includes
#include "TH1.h"

// c++ includes
#include <algorithm>
#include <iostream>
#include <memory>

using namespace std;

//_________________________________________________________________________
//
TH1D* HistoTools::TranformHistogramBinning(TH1* originalHist){

    //
    // In RooStats, input histogram variable binning is not supported => convert to a constant binning
    // by creating an histogram with the same number of bins but with constant binning between 0 and 1
    //
    // - now in case some bins are < 0 (due to bin drop functionality), they are ignored for the regBin histos
    const unsigned int nBins = originalHist -> GetNbinsX();
    unsigned int nBinsNew = 0;
    for(unsigned int iBin = 1; iBin <= nBins; ++iBin){
        if(originalHist->GetBinContent(iBin)>=0) nBinsNew++;
    }
    TH1D *hFinal = new TH1D(originalHist->GetName()+(TString)"_regBin",originalHist->GetTitle(),nBinsNew,0,1);
    hFinal -> SetDirectory(0);
    unsigned int iBinNew = 1;
    for(unsigned int iBin = 1; iBin <= nBins; ++iBin){
        if(originalHist->GetBinContent(iBin)<0) continue;
        hFinal -> SetBinContent(iBinNew,originalHist->GetBinContent(iBin));
        hFinal -> SetBinError(iBinNew,originalHist->GetBinError(iBin));
        iBinNew++;
    }
    //
    return hFinal;
}

//_________________________________________________________________________
//
void HistoTools::ManageHistograms( int histOps,  TH1* hNom, TH1* originUp, TH1* originDown,
                                    TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown, const SmoothOption &smoothOpt, bool TtresSmoothing) {
    //
    // Only function called directly to handle operations on the histograms (symmetrisation and smoothing)
    //

    //Sanity checks
    if( histOps % 10 > 2){
        WriteWarningStatus("HistoTools::ManageHistograms", "histOpt %10 > 2 - the operations to perform are not allowed ");
        WriteWarningStatus("HistoTools::ManageHistograms", "   two different symmetrisations ! Please check.");
        return;
    }
    if( histOps/10 > 99 ){
        WriteWarningStatus("HistoTools::ManageHistograms", "histOpt/10 > 99 - the operations to perform are not allowed ");
        WriteWarningStatus("HistoTools::ManageHistograms", "   two different symmetrisations ! Please check.");
        return;
    }

    // if one-sided & symmetrization asked, do smoothing first and symmetrization after
    if( histOps % 10 == SYMMETRIZEONESIDED ){
        SmoothHistograms(    histOps,hNom,modifiedUp,modifiedDown, smoothOpt, TtresSmoothing);
        SymmetrizeHistograms(histOps,hNom,modifiedUp,modifiedDown,modifiedUp,modifiedDown,scaleUp,scaleDown);
    }
    // otherwise, first symmetrization and then smoothing
    else{
        SymmetrizeHistograms(histOps,hNom,originUp,originDown,modifiedUp,modifiedDown,scaleUp,scaleDown);
        SmoothHistograms(    histOps,hNom,modifiedUp,modifiedDown, smoothOpt, TtresSmoothing);
    }
}

//_________________________________________________________________________
//
void HistoTools::SymmetrizeHistograms( int histOps,  TH1* hNom, TH1* originUp, TH1* originDown,
                                    TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown){
    //##################################################
    //
    // FIRST STEP: SYMMETRISATION
    //
    //##################################################
    if( histOps % 10 == SYMMETRIZEONESIDED ) {
        bool isUp = true; //is the provided uncertainty the up or down variation (based on yield)
        if     (originUp==nullptr && originDown!=nullptr) isUp = false;
        else if(originUp!=nullptr && originDown==nullptr) isUp = true;
        else if(originUp==nullptr && originDown==nullptr){
            WriteWarningStatus("HistoTools::SymmetrizeHistograms", "Both up and down variations are empty.");
            WriteWarningStatus("HistoTools::SymmetrizeHistograms", "Will not symmetrize.");
            return;
        }
        // if both are non-empty, check the differences with the nominal
        else{
            double separationUp = Separation(hNom,originUp);
            double separationDown = Separation(hNom,originDown);
            if( separationUp > separationDown ) isUp = true;
            if( separationUp < separationDown ) isUp = false;
        }


        TH1D* temp;
        if(isUp){
            temp = SymmetrizeOneSided(hNom, originUp, isUp);
            modifiedUp = (TH1*)originUp -> Clone();
            modifiedDown = (TH1*)temp -> Clone();
        } else {
            temp = SymmetrizeOneSided(hNom, originDown, isUp);
            modifiedUp = (TH1*)temp -> Clone();
            modifiedDown = (TH1*)originDown -> Clone();
        }
        delete temp;
    } else if ( histOps % 10 == SYMMETRIZETWOSIDED ) {
        modifiedUp = SymmetrizeTwoSided(originUp, originDown, hNom);
        modifiedDown = InvertShift(modifiedUp,hNom);
    } else {
        modifiedUp = originUp;
        modifiedDown = originDown;
    }
    Scale(modifiedDown, hNom, scaleDown);
    Scale(modifiedUp, hNom, scaleUp);
    modifiedDown -> SetName(originDown->GetName());
    modifiedUp   -> SetName(originUp->GetName());
}

//_________________________________________________________________________
//
void HistoTools::SmoothHistograms( int histOps,  TH1* hNom,
                                    TH1* &modifiedUp, TH1* &modifiedDown, const SmoothOption &smoothOpt, bool TtresSmoothing){
    //##################################################
    //
    // SECOND STEP: SMOOTHING
    //
    //##################################################

    // Initialize common smoothing tool
    SmoothHist smoothTool;

    if(hNom->GetNbinsX()==1){
        std::string temp =  hNom->GetName();
        WriteDebugStatus("HistoTools::SmoothHistograms", "In HistoTools::ManageHistograms(): skipping smoothing for systematics on \"" + temp + "\" since just 1 bin.");
        return;
    }
    if (TtresSmoothing || smoothOpt == TTBARRESONANCE) {
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothTtresDependent");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothTtresDependent");
        } else if( ( histOps - ( histOps % 10 ) ) >= SMOOTH_INDEPENDENT && (histOps - ( histOps % 10 ) ) < UNKNOWN ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothTtresIndependent");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothTtresIndependent");
        }
    } else if (smoothOpt == MAXVARIATION){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            smoothTool.setTRExTolerance(0.08); // This was also default before
            const int smoothingLevel = (histOps - ( histOps % 10 ) ) / 10;
            smoothTool.setTRExNbins(smoothingLevel);
            modifiedUp   = smoothTool.Smooth(hNom, modifiedUp,   "smoothTRExDefault");
            modifiedDown = smoothTool.Smooth(hNom, modifiedDown, "smoothTRExDefault");
        }
    } else if (smoothOpt == COMMONTOOLSMOOTHMONOTONIC){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothRebinMonotonic");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothRebinMonotonic");
        }
    } else if (smoothOpt == COMMONTOOLSMOOTHPARABOLIC){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothRebinParabolic");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothRebinParabolic");
        }
    } else if (smoothOpt == KERNELRATIOUNIFORM){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothRatioUniformKernel");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothRatioUniformKernel");
        }
    } else if (smoothOpt == KERNELDELTAGAUSS){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothDeltaGaussKernel");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothDeltaGaussKernel");
        }
    } else if (smoothOpt == KERNELRATIOGAUSS){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.Smooth(hNom, modifiedUp,   "smoothRatioGaussKernel");
            modifiedDown    = smoothTool.Smooth(hNom, modifiedDown, "smoothRatioGaussKernel");
        }
    } else {
        WriteWarningStatus("HistoTools::SmoothHistograms", "Unknown smoothing option. Please check the config file.");
        WriteWarningStatus("HistoTools::SmoothHistograms", "No smoothing will be applied.");
        return;
    }
}

//_________________________________________________________________________
//
TH1D* HistoTools::SymmetrizeOneSided( TH1* h_nominal, TH1* h_syst, bool &isUp ){

    float yield_nominal     = h_nominal->Integral();
    float yield_syst        = h_syst->Integral();

    //Convention: one sided systematic leading to larger yield: "up" variation
    if(yield_syst>yield_nominal) isUp = true;

    if(Separation(h_nominal,h_syst)<1e-05){
        std::string temp1 =  h_nominal -> GetName();
        std::string temp2 =  h_syst -> GetName();
        WriteWarningStatus("HistoTools::SymmetrizeOneSided", "The two histograms for one-sided symmetrisation are the same (difference is < 1e-5)");
        WriteWarningStatus("HistoTools::SymmetrizeOneSided", "Will not symmetrize");
        WriteWarningStatus("HistoTools::SymmetrizeOneSided", "      --> Nominal : " + temp1);
        WriteWarningStatus("HistoTools::SymmetrizeOneSided", "      --> Syst    : " + temp2);
    }

    return InvertShift(h_syst,h_nominal);
}

//_________________________________________________________________________
//
TH1D* HistoTools::InvertShift(TH1* h_syst, TH1* h_nominal){

    //Sanity check
    if(!h_syst || !h_nominal) return 0;

    //Just to be sure, set sumw2() on histograms
    if(!h_nominal->GetSumw2())h_nominal -> Sumw2();
    if(!h_syst->GetSumw2())h_syst -> Sumw2();

    //Compute the symmetric histogram
    TH1D* result = (TH1D*)h_nominal -> Clone();
    result -> SetDirectory(0);
    result -> Add(h_syst,-1);
    result -> Add(h_nominal,1);
    for( int iBin = 1; iBin <= result -> GetNbinsX(); ++iBin ){
        result -> SetBinError(iBin, h_syst->GetBinError(iBin));
    }

    //Another sanity check: search for negative bins
    for( int iBin = 1; iBin <= result -> GetNbinsX(); ++iBin ){
        double content = result -> GetBinContent(iBin);
        if(content < 0){
            result -> SetBinContent(iBin, 0.);
        }
    }
    return result;
}

//_________________________________________________________________________
//
float HistoTools::Separation(TH1* h1,TH1* h2){
    float sep = 0;
    for(int i_bin=1;i_bin<=h1->GetNbinsX();i_bin++){
        sep += TMath::Abs( h1->GetBinContent(i_bin) - h2->GetBinContent(i_bin) );
    }
    return sep;
}

//_________________________________________________________________________
//
TH1D* HistoTools::SymmetrizeTwoSided(TH1* var1, TH1* var2, TH1* hnom) {

    bool isLarge = false;
    if (std::fabs(var1->Integral()/hnom->Integral()-1) > 0.005) isLarge = true;
    if (std::fabs(var2->Integral()/hnom->Integral()-1) > 0.005) isLarge = true;
    //
    // Symmetrize a variation that is already two sided to smooth out possible fluctuations
    //
    //Nominal
    TH1D* nom = (TH1D*)hnom->Clone();

    //Up variation
    TH1D* tmp1 = (TH1D* )var1->Clone();
    tmp1->Divide(nom);
    if(!tmp1->GetSumw2())tmp1->Sumw2();

    //Down variation
    TH1D* tmp2 = (TH1D* )var2->Clone();
    tmp2->Divide(nom);
    if(!tmp2->GetSumw2())tmp2->Sumw2();

    std::unique_ptr<TH1D> nom_div(static_cast<TH1D*>(hnom->Clone()));
    nom_div->Divide(hnom);

    if (isLarge){
        std::vector<int> binSameShift{};
        // check if the variations are not in the same direction
        for (int ibin = 1; ibin <= hnom->GetNbinsX(); ibin++){
            const double& nom_error = nom_div->GetBinError(ibin);
            const double& up_content = tmp1->GetBinContent(ibin);
            const double& up_error = tmp1->GetBinError(ibin);
            const double& down_content = tmp2->GetBinContent(ibin);
            const double& down_error = tmp2->GetBinError(ibin);

            const double total_uncert_up = std::sqrt(up_error*up_error + nom_error*nom_error);
            const double total_uncert_down = std::sqrt(down_error*down_error + nom_error*nom_error);

            // check if the variations are larger than stat uncertainty
            if ((std::fabs(up_content - 1) > total_uncert_up) && (std::fabs(down_content - 1) > total_uncert_down)){
                // check if the variations have the same shift
                // have same sign = same shift
                if ((up_content - 1) * (down_content - 1) > 0){
                    binSameShift.emplace_back(ibin);
                }
            }
        }
        if (binSameShift.size() > 0){
            const std::string& name_up = var1->GetName();
            const std::string& name_down = var2->GetName();
            std::string tmp = "";
            for (const auto& ibin : binSameShift){
                tmp+= std::to_string(ibin) + " ";
            }
            WriteWarningStatus("HistoTools::SymmetrizeTwoSided", "Histograms " + name_up + " and " + name_down + " in bins " + tmp + " have the same shift for up and down variation");
            WriteWarningStatus("HistoTools::SymmetrizeTwoSided", "You should check this");
        }
    }

    //Flat (content = 0) histogram to substract
    TH1D* unit = (TH1D* )nom->Clone();
    if(!unit->GetSumw2())unit->Sumw2();
    for (int bin=1; bin<= unit->GetNbinsX(); bin++){
        unit->SetBinContent(bin,1);
        unit->SetBinError(bin,0.0);
    }

    //Computes Var/Nom - 1
    tmp1->Add(unit,-1);
    tmp2->Add(unit,-1);

    //Computes Corrected = (DeltaUp-DeltaDown)/2 + 1
    tmp1->Add(tmp2,-1);
    tmp1->Scale(0.5);
    tmp1->Add(unit);

    //Computed the new histogram Corrected*Nominal
    tmp1->Multiply(nom);

    //Protection against negative bin
    for (int bin=1; bin<= unit->GetNbinsX(); bin++){
        double content = tmp1->GetBinContent(bin);
        if(content<0){
            tmp1->SetBinContent(bin,0.);
        }
        //tmp1->SetBinError(bin,nom->GetBinError(bin));
    }

    delete tmp2;
    delete unit;
    delete nom;

    return tmp1;
}

//_________________________________________________________________________
//
void HistoTools::Scale(TH1* h_syst, TH1* h_nominal, float factor){

    //Sanity check
    if(!h_syst || !h_nominal) return;

    //Just to be sure, set sumw2() on histograms
    if(!h_nominal->GetSumw2())h_nominal -> Sumw2();
    if(!h_syst->GetSumw2())h_syst -> Sumw2();

    // clone nominal histogram
    std::unique_ptr<TH1> h_nominal_tmp(static_cast<TH1*>(h_nominal->Clone()));
    // set errors to zero
    for(int i_bin=0;i_bin<=h_nominal_tmp->GetNbinsX()+1;i_bin++){
        h_nominal_tmp->SetBinError(i_bin,0.);
    }

    //scale difference to nominal
    h_syst -> Add(h_nominal_tmp.get(),-1);
    h_syst -> Scale(factor);
    h_syst -> Add(h_nominal_tmp.get(),1);

    //Another sanity check: search for negative bins
    for( int iBin = 1; iBin <= h_syst -> GetNbinsX(); ++iBin ){
        double content = h_syst -> GetBinContent(iBin);
        if(content < 0){
            h_syst -> SetBinContent(iBin, 0.);
        }
    }
}

double HistoTools::avgError(std::vector<Bin> &hist, bool independentVar) {
  int Nbins = hist.size();
  double avg = 0;
  std::vector<double> errs;
  for (int k = 0; k < Nbins; ++k) {
    double dM = 0;
    if (independentVar) dM = sqrt(hist[k].dS2 + hist[k].dN2);
    else dM = max(sqrt(hist[k].dN2), sqrt(hist[k].dS2));
    errs.push_back(dM/hist[k].N);
    avg += dM/hist[k].N;
  }
  std::sort(errs.begin(), errs.end());
  return avg/((double) Nbins);
}

bool HistoTools::systSmallerThanStat(std::vector<Bin> &hist, bool independentVar, double avgErr) {
    int Nbins = hist.size();
    for (int k = 0; k < Nbins; ++k) {
        double dM = 0;
        if (independentVar) dM = sqrt(hist[k].dS2 + hist[k].dN2);
        else dM = max(sqrt(hist[k].dN2), sqrt(hist[k].dS2));
        if (independentVar) {
            if (fabs(hist[k].S - hist[k].N) < dM && fabs(hist[k].S - hist[k].N)/hist[k].N > avgErr)
                //if (fabs(hist[k].S - hist[k].N) < dM)
                return true;
        } else {
            if (fabs(hist[k].S - hist[k].N) < dM)
                return true;
        }
    }
    return false;
}

//_________________________________________________________________________
//
bool HistoTools::HasShape(TH1* hnom, SystematicHist* sh, float threshold){

    //Save time
    if(hnom->GetNbinsX()==1) return false;

    //Integral -> Protects against nan's for Shape histograms
    double integralUp = sh->fHistShapeUp->Integral();
    if(integralUp!=integralUp) return false;

    double integralDown = sh->fHistShapeDown->Integral();
    if(integralDown!=integralDown) return false;

    //If at least one bin is the shape histogram is larger than the threshold, keep the uncertainty
    bool hasShape = false;
    for (int iBin = 1; iBin <= hnom->GetNbinsX(); ++iBin) {
        double nom = hnom->GetBinContent(iBin);
        double up = sh->fHistShapeUp->GetBinContent(iBin);
        double down = sh->fHistShapeDown->GetBinContent(iBin);

        double up_err = TMath::Abs((up-nom)/nom);
        double down_err = TMath::Abs((down-nom)/nom);
        if(up_err>=threshold || down_err>=threshold){
            hasShape = true;
            break;
        }
    }
    return hasShape;
}


//_________________________________________________________________________
//
bool HistoTools::CheckHistograms(TH1* nom, SystematicHist* sh, bool checkNullContent, bool causeCrash){

    //======================================================
    // This code performs sanity checks on the provided
    // histograms.
    //======================================================

    //
    // Two kinds of detection
    //    -> ERRORS: fundamental problems (non sense of fit, missing inputs, ...)
    //    -> WARNINGS: potential issues (possibly affecting the fits or the plots)
    //

    bool isGood = true;

    //
    // 1) Checks that the histograms exist
    //
    if(nom == nullptr){
        if (causeCrash){
            WriteErrorStatus("HistoTools::CheckHistograms", "The nominal histogram doesn't seem to exist !");
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("HistoTools::CheckHistograms", "The nominal histogram doesn't seem to exist !");
        }
        return false;
    }
    if(sh != nullptr){
        if(sh->fHistUp == nullptr){
            if (causeCrash){
                WriteErrorStatus("HistoTools::CheckHistograms", "The up variation histogram doesn't seem to exist !");
                exit(EXIT_FAILURE);
            } else {
                WriteWarningStatus("HistoTools::CheckHistograms", "The up variation histogram doesn't seem to exist !");
            }
            return false;
        }
        if(sh->fHistDown == nullptr){
            if (causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "The down variation histogram doesn't seem to exist !");
                exit(EXIT_FAILURE);
            } else {
                WriteWarningStatus("HistoTools::CheckHistograms", "The down variation histogram doesn't seem to exist !");
            }
            return false;
        }
    }
    else return false;

    //
    // 2) Checks the binning is the same for the three histograms
    //
    const int NbinsNom  = nom->GetNbinsX();
    const int NbinsUp   = sh->fHistUp->GetNbinsX();
    const int NbinsDown = sh->fHistDown->GetNbinsX();
    if( (NbinsNom != NbinsUp) || (NbinsNom != NbinsDown) || (NbinsUp != NbinsDown) ){
        if (causeCrash) {
            WriteErrorStatus("HistoTools::CheckHistograms", "The number of bins is found inconsistent ! Please check!");
            exit(EXIT_FAILURE);
        } else {
            WriteWarningStatus("HistoTools::CheckHistograms", "The number of bins is found inconsistent ! Please check!");
        }
        return false;
    }

    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double lowEdgeNom   = nom->GetBinLowEdge(iBin);
        double lowEdgeUp    = sh->fHistUp->GetBinLowEdge(iBin);
        double lowEdgeDown  = sh->fHistDown->GetBinLowEdge(iBin);

        if( abs(lowEdgeNom-lowEdgeUp)>1e-05 || abs(lowEdgeNom-lowEdgeDown)>1e-05 || abs(lowEdgeDown-lowEdgeUp)>1e-05 ){
            if (causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "The bin low edges are not consistent ! Please check !");
                exit(EXIT_FAILURE);
            } else {
                WriteWarningStatus("HistoTools::CheckHistograms", "The bin low edges are not consistent ! Please check !");
            }
            return false;
        }

        double binWidthNom   = nom->GetBinWidth(iBin);
        double binWidthUp    = sh->fHistUp->GetBinWidth(iBin);
        double binWidthDown  = sh->fHistDown->GetBinWidth(iBin);

        if( abs(binWidthNom-binWidthUp)>1e-05 || abs(binWidthNom-binWidthDown)>1e-05 || abs(binWidthDown-binWidthUp)>1e-05 ){
            if (causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "The bin widths are not consistent ! Please check !");
                exit(EXIT_FAILURE);
            } else {
                WriteWarningStatus("HistoTools::CheckHistograms", "The bin widths are not consistent ! Please check !");
            }
            return false;
        }
    }

    //
    // 3) Checks the absence on bins with 0 content (for nominal)
    //
    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double content = nom->GetBinContent(iBin);
        if( ( checkNullContent && content<=0 ) || ( !checkNullContent && content<0 ) ){
            std::string temp = nom->GetName();
            isGood = false;
            if(causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp + "\", bin " + std::to_string(iBin) + " has 0 content ! Please check");
                WriteErrorStatus("HistoTools::CheckHistograms", "Nominal: " + std::to_string(content));
              exit(EXIT_FAILURE);
            } else {
              WriteWarningStatus("HistoTools::CheckHistograms", "In histo \"" + temp + "\", bin " + std::to_string(iBin) + " has 0 content ! Please check");
              WriteWarningStatus("HistoTools::CheckHistograms", "Nominal: " + std::to_string(content));
              //Corrects the nominal
              WriteWarningStatus("HistoTools::CheckHistograms", "I set the bin content to 1e-05 pm 1e-06 ! Please check !");
              nom -> SetBinContent(iBin,1e-05);
              nom -> SetBinError(iBin, 1e-06);
            }
        }

        //Now, for those histograms, checks if a systematics also has 0 content to
        //avoid 100% down systematics
        if( TMath::Abs( nom->GetBinContent(iBin)-1e-05) < 1e-10  ){
            //This bin has most likely been changed to the default non-zero value
            sh->fHistUp->SetBinContent(iBin,1e-05);
            sh->fHistDown->SetBinContent(iBin,1e-05);
        }

    }//loop over the bins


    //
    // 4) Check the presence of abnormal bins (content-based)
    //
    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double contentNom   = nom->GetBinContent(iBin);
        double contentUp    = sh->fHistUp->GetBinContent(iBin);
        double contentDown  = sh->fHistDown->GetBinContent(iBin);

        //
        // 4.a) Checks the presence of negative bins for systematic
        //
        if(contentUp<0){
            std::string temp_string = sh->fHistUp->GetName();
            isGood = false;
            if(causeCrash){
                WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has negative content ! Please check");
                WriteErrorStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
                exit(EXIT_FAILURE);
            } else {
                WriteWarningStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has negative content ! Please check");
                WriteWarningStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
                WriteWarningStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
                WriteWarningStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
                WriteWarningStatus("HistoTools::CheckHistograms", "  => Setting Up to 1e-06");
                sh->fHistUp->SetBinContent(iBin,1e-06);
            }
        }
        if(contentDown<0){
            std::string temp_string = sh->fHistDown->GetName();
            isGood = false;
            if(causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has negative content ! Please check");
                WriteErrorStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
                exit(EXIT_FAILURE);
            }
            else{
                WriteWarningStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has negative content ! Please check");
                WriteWarningStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
                WriteWarningStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
                WriteWarningStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
                WriteWarningStatus("HistoTools::CheckHistograms", "  => Setting Up to 1e-06");
                sh->fHistDown->SetBinContent(iBin,1e-06);
            }
        }

        //
        // 4.b) Checks that the systematics are not crazy (too large ratio, nan returned, ...)
        //
        contentNom   = nom->GetBinContent(iBin);
        contentUp    = sh->fHistUp->GetBinContent(iBin);
        contentDown  = sh->fHistDown->GetBinContent(iBin);
        //
        double ratioUp   = 0.;
        double ratioDown = 0.;
        if(contentNom != 0 ){
            ratioUp   = contentUp  /contentNom;
            ratioDown = contentDown/contentNom;
        } else if( TMath::Abs(contentUp)>0 || TMath::Abs(contentDown)>0 ) {
            std::string temp_string = sh->fHistUp->GetName();
            WriteWarningStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string);
            WriteWarningStatus("HistoTools::CheckHistograms", "\", bin " + std::to_string(iBin) + " has null nominal content but systematics are >0! Hope you know that");
            continue;
        }

        // * first check is for nan
        // * second check is for negative values (normally impossible after step 3 and 4.a)
        // * third check is for abnormal systematic/nominal values (ratio > 100)
        if( (ratioUp!=ratioUp) || (ratioUp < 0) || (abs(ratioUp-1.) >= 100) ){
            std::string temp_string = sh->fHistUp->GetName();
            isGood = false;
            if(causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has weird content ! Please check");
                WriteErrorStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
                exit(EXIT_FAILURE);
            }
            WriteWarningStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has weird content ! Please check");
            WriteWarningStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
            WriteWarningStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
            WriteWarningStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
            // try to fix it, if not aborting
            if(ratioUp!=ratioUp) sh->fHistUp->SetBinContent(iBin,contentNom);
            else return isGood;
        }
        if( (ratioDown!=ratioDown) || (ratioDown < 0) || (abs(ratioDown-1.) >= 100) ){
            std::string temp_string = sh->fHistDown->GetName();
            isGood = false;
            if(causeCrash) {
                WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has weird content ! Please check");
                WriteErrorStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
                WriteErrorStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
                exit(EXIT_FAILURE);
            }
            WriteWarningStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has weird content ! Please check");
            WriteWarningStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
            WriteWarningStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
            WriteWarningStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
            // try to fix it, if not aborting
            if(ratioDown!=ratioDown) sh->fHistDown->SetBinContent(iBin,contentNom);
            else return isGood;
        }
    }

    //
    // 5) Warning level: check the presence of underflow/overflow (ignored in the code)
    //
    double underflowNom     = nom -> GetBinContent(0);
    double underflowUp      = sh->fHistUp -> GetBinContent(0);
    double underflowDown    = sh->fHistDown -> GetBinContent(0);
    if( abs(underflowNom)>0 || abs(underflowUp)>0 || abs(underflowDown)>0 ){
        std::string temp1 = nom -> GetName();
        std::string temp2 = sh->fHistUp -> GetName();
        std::string temp3 = sh->fHistDown -> GetName();
        WriteWarningStatus("HistoTools::CheckHistograms", "Underflow detected ! This will not be taken into account.");
        WriteWarningStatus("HistoTools::CheckHistograms", "Nominal (" + temp1 + "): " + std::to_string(underflowNom));
        WriteWarningStatus("HistoTools::CheckHistograms", "Nominal (" + temp2 + "): " + std::to_string(underflowUp));
        WriteWarningStatus("HistoTools::CheckHistograms", "Nominal (" + temp3 + "): " + std::to_string(underflowDown));
    }

    double overflowNom      = nom -> GetBinContent(NbinsNom+1);
    double overflowUp       = sh->fHistUp -> GetBinContent(NbinsUp+1);
    double overflowDown     = sh->fHistDown -> GetBinContent(NbinsDown+1);
    if( abs(overflowNom)>0 || abs(overflowUp)>0 || abs(overflowDown)>0 ){
        std::string temp1 = nom -> GetName();
        std::string temp2 = sh->fHistUp -> GetName();
        std::string temp3 = sh->fHistDown -> GetName();
        WriteWarningStatus("HistoTools::CheckHistograms", "Overflowflow detected ! This will not be taken into account.");
        WriteWarningStatus("HistoTools::CheckHistograms", "Nominal (" + temp1 + "): " + std::to_string(overflowNom));
        WriteWarningStatus("HistoTools::CheckHistograms", "Nominal (" + temp2 + "): " + std::to_string(overflowUp));
        WriteWarningStatus("HistoTools::CheckHistograms", "Nominal (" + temp3 + "): " + std::to_string(overflowDown));
    }
    return isGood;
}
