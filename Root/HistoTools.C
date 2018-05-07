/*

 HistoTools
 =========

 Contains all functions needed to handle properly the histograms:
    -> Variable Binning handling for fit/limits
    -> Symmetrisation of systematics
    -> Smoothing of systematics

 Call of the functions:
    -> #include "TtHFitter/HistoTools.C"
    -> Call of the function with HistoTools:: (Helps readability)

 Contact: Loic Valery <lvalery@cern.ch>

 */

#include <iostream>
#include "TH1.h"

#include "TtHFitter/HistoTools.h"
#include "TtHFitter/Common.h"
#include "TtHFitter/SystematicHist.h"
#include "TtHFitter/StatusLogbook.h"

#include "CommonSystSmoothingTool/SmoothSystematics/SmoothHist.h"

using namespace std;

//_________________________________________________________________________
//
TH1F* HistoTools::TranformHistogramBinning(TH1* originalHist){

    //
    // In RooStats, input histogram variable binning is not supported => convert to a constant binning
    // by creating an histogram with the same number of bins but with constant binning between 0 and 1
    //
//     const unsigned int nBins = originalHist -> GetNbinsX();
//     TH1F *hFinal = new TH1F(originalHist->GetName()+(TString)"_regBin",originalHist->GetTitle(),nBins,0,1);
//     hFinal -> SetDirectory(0);
//     for(unsigned int iBin = 1; iBin <= nBins; ++iBin){
//         hFinal -> SetBinContent(iBin,originalHist->GetBinContent(iBin));
//         hFinal -> SetBinError(iBin,originalHist->GetBinError(iBin));
//     }
    // Updated:
    // - now in case some bins are < 0 (due to bin drop functionality), they are ignored for the regBin histos
    const unsigned int nBins = originalHist -> GetNbinsX();
    unsigned int nBinsNew = 0;
    for(unsigned int iBin = 1; iBin <= nBins; ++iBin){
        if(originalHist->GetBinContent(iBin)>=0) nBinsNew++;
    }
    TH1F *hFinal = new TH1F(originalHist->GetName()+(TString)"_regBin",originalHist->GetTitle(),nBinsNew,0,1);
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
                                    TH1* &modifiedUp, TH1* &modifiedDown, float scaleUp, float scaleDown, const SmoothOption &smoothOpt, bool TtresSmoothing, std::string kernelOpt, std::string kernelSmoothType) {
    //
    // Only function called directly to handle operations on the histograms (symmetrisation and smoothing)
    //

    //Sanity checks
    if( histOps % 10 > 2){
        WriteErrorStatus("HistoTools::ManageHistograms", "histOpt %10 > 2 - the operations to perform are not allowed ");
        WriteErrorStatus("HistoTools::ManageHistograms", "   two different symmetrisations ! Please check.");
        return;
    }
    if( histOps/10 > 99 ){
        WriteErrorStatus("HistoTools::ManageHistograms", "histOpt/10 > 99 - the operations to perform are not allowed ");
        WriteErrorStatus("HistoTools::ManageHistograms", "   two different symmetrisations ! Please check.");
        return;
    }

    // if one-sided & symmetrization asked, do smoothing first and symmetrization after
    if( histOps % 10 == SYMMETRIZEONESIDED ){
        SmoothHistograms(    histOps,hNom,originUp,originDown,modifiedUp,modifiedDown, smoothOpt, TtresSmoothing, kernelOpt, kernelSmoothType);
        SymmetrizeHistograms(histOps,hNom,modifiedUp,modifiedDown,modifiedUp,modifiedDown,scaleUp,scaleDown);
    }
    // otherwise, first symmetrization and then smoothing
    else{
        SymmetrizeHistograms(histOps,hNom,originUp,originDown,modifiedUp,modifiedDown,scaleUp,scaleDown);
        SmoothHistograms(    histOps,hNom,originUp,originDown,modifiedUp,modifiedDown, smoothOpt, TtresSmoothing, kernelOpt, kernelSmoothType);
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
        if     (originUp==0x0 && originDown!=0x0) isUp = false;
        else if(originUp!=0x0 && originDown==0x0) isUp = true;
        else if(originUp==0x0 && originDown==0x0){
            WriteErrorStatus("HistoTools::SymmetrizeHistograms", "both up and down variations are empty.");
            return;
        }
        // if both are non-empty, check the differences with the nominal
        else{
            double separationUp = Separation(hNom,originUp);
            double separationDown = Separation(hNom,originDown);
            if( separationUp > separationDown ) isUp = true;
            if( separationUp < separationDown ) isUp = false;
        }


        TH1F* temp;
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
void HistoTools::SmoothHistograms( int histOps,  TH1* hNom, TH1* originUp, TH1* originDown,
                                    TH1* &modifiedUp, TH1* &modifiedDown, const SmoothOption &smoothOpt, bool TtresSmoothing, std::string kernelOpt, std::string kernelSmoothType){
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
            modifiedUp      = smoothTool.Smooth_Ttres(hNom, originUp,   false);
            modifiedDown    = smoothTool.Smooth_Ttres(hNom, originDown, false);
        } else if( ( histOps - ( histOps % 10 ) ) >= SMOOTH_INDEPENDENT && (histOps - ( histOps % 10 ) ) < UNKNOWN ){
            modifiedUp      = smoothTool.Smooth_Ttres(hNom, originUp,   true);
            modifiedDown    = smoothTool.Smooth_Ttres(hNom, originDown, true);
        }
    } else if (smoothOpt == MAXVARIATION){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            const int smoothingLevel = (histOps - ( histOps % 10 ) ) / 10;
            Smooth_maxVariations( modifiedUp,    hNom,   smoothingLevel );
            Smooth_maxVariations( modifiedDown,  hNom,   smoothingLevel );
        }
    } else if (smoothOpt == COMMONTOOLSMOOTH){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){
            modifiedUp      = smoothTool.smoothHistogram(hNom, originUp, true);
            modifiedDown    = smoothTool.smoothHistogram(hNom, originDown, true);
        }
    } else if (smoothOpt == KERNELFUNCTION){
        if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < SMOOTH_INDEPENDENT ){

            // choose kernel option
            smoothTool.setKernelOption(kernelOpt);
    
            // choose smoothing - ratio or difference
            smoothTool.setSmoothType(kernelSmoothType); // or delta

            // Set list of kernel radii
            //smoothTool.setSpans({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
            //                     1.0, 1.2, 1.4, 1.6, 1.8, 2.0
            //                    });
            
            /*In the default span list, maximum span is 2.
             * If you bin width is larger then 2, you should call this fuction.
             * */
            if (GetMaxBinWidth(hNom) > 2 ){
                smoothTool.useRalativeSpans(true);
            } else {
                smoothTool.useRalativeSpans(false);
            }

            modifiedUp      = smoothTool.smoothWithKernel(hNom, originUp);
            modifiedDown    = smoothTool.smoothWithKernel(hNom, originDown);
        }
    } else {
        WriteWarningStatus("HistoTools::SmoothHistograms", "Unknown smoothing option. Please check the config file.");
        WriteWarningStatus("HistoTools::SmoothHistograms", "No smoothing will be applied.");
        return; 
    }
}

//_________________________________________________________________________
//
TH1F* HistoTools::SymmetrizeOneSided( TH1* h_nominal, TH1* h_syst, bool &isUp ){

    float yield_nominal     = h_nominal->Integral();
    float yield_syst        = h_syst->Integral();

    //Convention: one sided systematic leading to larger yield: "up" variation
    if(yield_syst>yield_nominal) isUp = true;

    if(Separation(h_nominal,h_syst)<1e-05){
        std::string temp1 =  h_nominal -> GetName();
        std::string temp2 =  h_syst -> GetName();
        WriteErrorStatus("HistoTools::SymmetrizeOneSided", "The two histograms for one-sided symmetrisation are the same ...");
        WriteErrorStatus("HistoTools::SymmetrizeOneSided", "Please check.");
        WriteErrorStatus("HistoTools::SymmetrizeOneSided", "      --> Nominal : " + temp1);
        WriteErrorStatus("HistoTools::SymmetrizeOneSided", "      --> Syst    : " + temp2);
    }

    return InvertShift(h_syst,h_nominal);
}

//_________________________________________________________________________
//
TH1F* HistoTools::InvertShift(TH1* h_syst, TH1* h_nominal){

    //Sanity check
    if(!h_syst || !h_nominal) return 0;

    //Just to be sure, set sumw2() on histograms
    if(!h_nominal->GetSumw2())h_nominal -> Sumw2();
    if(!h_syst->GetSumw2())h_syst -> Sumw2();

    //Compute the symmetric histogram
    TH1F* result = (TH1F*)h_nominal -> Clone();
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
TH1F* HistoTools::SymmetrizeTwoSided(TH1* var1, TH1* var2, TH1* hnom) {

    //
    // Symmetrize a variation that is already two sided to smooth out possible fluctuations
    //
    //Nominal
    TH1F* nom = (TH1F*)hnom->Clone();

    //Up variation
    TH1F* tmp1 = (TH1F* )var1->Clone();
    tmp1->Divide(nom);
    if(!tmp1->GetSumw2())tmp1->Sumw2();

    //Down variation
    TH1F* tmp2 = (TH1F* )var2->Clone();
    tmp2->Divide(nom);
    if(!tmp2->GetSumw2())tmp2->Sumw2();

    //Flat (content = 0) histogram to substract
    TH1F* unit = (TH1F* )nom->Clone();
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

    //scale difference to nominal
    h_syst -> Add(h_nominal,-1);
    h_syst -> Scale(factor);
    h_syst -> Add(h_nominal,1);

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
  //int s = errs.size();
  //if (s % 2 == 0) return (errs[s/2-1] + errs[s/2+1])/2.0;
  //return errs[s/2];
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
void HistoTools::Smooth_maxVariations(TH1* hsyst, TH1* hnom, int nbins){

    //
    // General idea: avoid having more than "nbins" slope variations in the systematic histogram
    //
    float systIntegral = hsyst->Integral();

    double tolerance = 0.08;

    int nVar = rebin_getMaxVar(hnom,hsyst,tolerance);

    WriteVerboseStatus("HistoTools::Smooth_maxVariations", "---: " + std::to_string(tolerance) + " " + std::to_string(nVar));

    //
    // Iterates the smoothing of the systematic histogram until the number a slope changes is lower than "nbins"
    //
    while (nVar > nbins){
        tolerance = tolerance/2.;
        nVar = rebin_getMaxVar(hnom,hsyst,tolerance);
        WriteVerboseStatus("HistoTools::Smooth_maxVariations", "---: " + std::to_string(tolerance) + " " + std::to_string(nVar) );
        if(tolerance==0){
            std::string temp1 = hnom->GetName();
            std::string temp2 = hnom->GetTitle();
            WriteErrorStatus("HistoTools::Smooth_maxVariations", "Buuuuuuuuuuug: infinite while");
            WriteErrorStatus("HistoTools::Smooth_maxVariations", temp1 + " " + temp2 + " nbins: " + std::to_string(nbins));
            break;
        }
    }
    WriteVerboseStatus("HistoTools::Smooth_maxVariations", "Final: " + std::to_string(tolerance) + " " + std::to_string(nVar));

    TH1F *ratio = (TH1F *) hsyst->Clone();
    ratio->Divide(hnom);
    hsyst->Add(hnom,-1);
    hsyst->Divide(hnom);

    //First non-empty bin
    int minbin = 0;
    for(int i=1; i < ratio->GetNbinsX()+1;i++){
        if(ratio->GetBinContent(i)!=0){
            minbin = i;
            break;
        }
    }

    //Last non-empty bin
    int maxbin = ratio->GetNbinsX();
    for(int i=maxbin; i>= 1;i--){
        if(ratio->GetBinContent(i)!=0){
            maxbin = i;
            break;
        }
    }

    // Smooth only works well for positive entries: shifts all entries by an offset of 100
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        hsyst->SetBinContent( i, hsyst->GetBinContent(i) + 100 );
        //         hsyst->SetBinContent( i, hsyst->GetBinContent(i) + 1000 );
    }

    // Due to the rebinning, some bins can have the same content. Call the ROOT smooth function to avoid this.
    int binwidth = getBinWidth(ratio);
    if(binwidth>4) binwidth=4;
    hsyst->GetXaxis()->SetRange(minbin,maxbin);
    if(binwidth*2<maxbin-minbin){
        hsyst->Smooth(binwidth*2,"R");
    }

    // Removes the 100 offset
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        hsyst->SetBinContent( i, hsyst->GetBinContent(i) - 100 );
        //         hsyst->SetBinContent( i, hsyst->GetBinContent(i) - 1000 );
    }
    hsyst->Multiply(hnom);
    hsyst->Add(hnom);
    if(hsyst->Integral()!=0){
        hsyst->Scale(systIntegral/hsyst->Integral());
    }
    // Checks if any bin with < 0 content exists
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        double content = hsyst->GetBinContent( i );
        if(content < 0){
            hsyst -> SetBinContent(i, 0.);
        }
    }
    delete ratio;
}

//_________________________________________________________________________
//
int HistoTools::getBinWidth(TH1 *ratio){

    //
    // Returns the minimal number of consecutive bins with the same content
    //
    float prev=0;
    int count=1, mincount=99;

    for(int iBin = 1;iBin <= ratio->GetNbinsX(); ++iBin){

        if(ratio->GetBinContent(iBin)==0) continue;

        if(prev!=0){
            if( TMath::Abs(prev-ratio->GetBinContent(iBin))< 1e-05 ){
                count++;
            } else {
                if(count < mincount) mincount=count;
                count=1;
            }
        }
        prev = ratio->GetBinContent(iBin);
    }
    if(count < mincount) mincount=count;
    return mincount;
}

//_________________________________________________________________________
//
int HistoTools::rebin_getMaxVar(TH1* hnom,TH1* hsyst, double tolerance){

  //
  // Recompute the systematic histogram based on a new binning (built based on the relative stat uncertainty to suppress
  // statistical fluctuations)
  //

  for(int i=1;i<=hsyst->GetNbinsX();i++){
    WriteVerboseStatus("HistoTools::rebin_getMaxVar", "In: " + std::to_string(hnom->GetBinContent(i)) + " " + std::to_string(hsyst->GetBinContent(i)));
  }

  std::vector<double> binLimit;
  binLimit.push_back( hnom->GetXaxis()->GetXmin() );
  double relErr=20000;
  double cumulIntSyst=0;
  double cumulInt=0;
  double cumulErr=0;
  int thisBin=0;

  do {

    do { //while (relErr > tolerance && thisBin!=hnom->GetNbinsX() );

          //
          // Compute the relative statistical uncertainty of a group of bins. Performs this operation until
          // the relative statistical uncertainty is lower than the tolerance (or the number of bins)
          //

          thisBin++;

          cumulInt+=fabs(hnom->GetBinContent(thisBin));
          cumulErr+=hnom->GetBinError(thisBin)*hnom->GetBinError(thisBin);
          cumulIntSyst += hsyst->GetBinContent(thisBin);
          if ( cumulInt!=0 && cumulIntSyst!=0 ){
              relErr= sqrt(cumulErr)/cumulInt;
          }
          if (relErr==0) relErr=20000;

          WriteVerboseStatus("HistoTools::rebin_getMaxVar", std::to_string(thisBin) + " " + std::to_string(hnom->GetBinContent(thisBin)) + " " + std::to_string(hnom->GetBinError(thisBin)));
          WriteVerboseStatus("HistoTools::rebin_getMaxVar", std::to_string(std::sqrt(cumulErr)) + " " + std::to_string(cumulInt) + " " + std::to_string(relErr) + " " + std::to_string(tolerance));

        } while (relErr > tolerance && thisBin!=hnom->GetNbinsX() );

        if((relErr < tolerance) || binLimit.size() == 1){//a group of bins with a sufficiently low stat error has been found, let's add it
            binLimit.push_back(hnom->GetBinCenter(thisBin)+hnom->GetBinWidth(thisBin)/2);
        } else {//no such group of bins has been found: merge with the last found bin
            binLimit.back() = hnom->GetBinCenter(thisBin)+hnom->GetBinWidth(thisBin)/2;
        }

        WriteVerboseStatus("HistoTools::rebin_getMaxVar", "Push back bin: " + std::to_string(thisBin));

        cumulInt=0;
        cumulErr=0;
        relErr=20000;

    } while ( thisBin!=hnom->GetNbinsX() );

    //
    // Define the array for new bin ranges and a template histogram with this binning
    //
    double* Bins;
    Bins=new double[binLimit.size()];
    for ( unsigned int i=0; i< binLimit.size(); i++) {
        Bins[i]=binLimit[i];
    }

    TH1F* rebinTemplate=new TH1F( "binRef", "binRef", binLimit.size()-1, Bins);
    rebinTemplate->Sumw2();
    for (int V=1; V<=rebinTemplate->GetNbinsX(); V++) {
        rebinTemplate->SetBinContent(V,V);
    }

    //
    // Performs a rebin "by hand" of the nominal and systematic histograms
    //
    //nominal
    TH1F* hnomBinned=new TH1F(*rebinTemplate);
    hnomBinned->Reset();
    for (int bin=0; bin <=hnom->GetNbinsX(); bin++) {
        int bigBin=hnomBinned->FindBin( hnom->GetBinCenter(bin));
        hnomBinned->SetBinContent( bigBin, hnomBinned->GetBinContent(bigBin)+ hnom->GetBinContent(bin) );
        hnomBinned->SetBinError( bigBin, sqrt( pow(hnomBinned->GetBinError(bigBin),2) + pow(hnom->GetBinError(bin),2) ) );
    }
    //systematics
    TH1F* hsystBinned=new TH1F(*rebinTemplate);
    hsystBinned->Reset();
    for (int bin=0; bin <=hsyst->GetNbinsX(); bin++) {
        int bigBin=hsystBinned->FindBin( hsyst->GetBinCenter(bin));
        hsystBinned->SetBinContent( bigBin, hsystBinned->GetBinContent(bigBin)+ hsyst->GetBinContent(bin) );
        hsystBinned->SetBinError( bigBin, sqrt( pow(hsystBinned->GetBinError(bigBin),2) + pow(hsyst->GetBinError(bin),2) ) );
    }
    hsystBinned->Divide(hnomBinned);//now, hsystBinned is the relative systematic uncertainty

    //
    // Modify the systematic uncertainty histogram based on the "stat reliable" ratio to avoid statistical fluctuations
    //
    for(int i=1;i<=hsyst->GetNbinsX();i++){
        // systematic_new = nominal(with normal binning) * relative_systematic_uncertainty(with new binning)
        hsyst->SetBinContent(i,hnom->GetBinContent(i)*hsystBinned->GetBinContent(hsystBinned->FindBin( hsyst->GetBinCenter(i))));
    }

    //
    // Computes the number of slope variations in the new systematic histogram
    //
    int nVar = get_nVar(hsystBinned);

    delete rebinTemplate;
    delete hnomBinned;
    delete hsystBinned;
    delete [] Bins;

    return nVar;
}

//_________________________________________________________________________
//
int HistoTools::get_nVar(TH1* hratio){

    //
    // Counts the number of slope changes
    //

    int nVar=0;
    bool thisUp = true;
    bool goingUp = true;
    double prevBin = 0, thisBin = 0;
    int usedBins = 0;
    bool hasChange = false;

    for (int bin=1; bin <=hratio->GetNbinsX(); bin++) {

        WriteVerboseStatus("HistoTools::get_nVar", "Bin/content: " + std::to_string(bin) + " " + std::to_string(hratio->GetBinContent(bin)));

        if(hratio->GetBinContent(bin)==0) continue;

        prevBin = thisBin;
        thisBin = hratio->GetBinContent(bin);


        WriteVerboseStatus("HistoTools::get_nVar", "Prev/this: " + std::to_string(prevBin) + " " + std::to_string(thisBin));

        if(fabs(thisBin-prevBin)<1e-7) continue;//skip bins with same content

        usedBins++;

        thisUp = (thisBin > prevBin);
        if (usedBins > 2 && ((thisUp && !goingUp) || (!thisUp && goingUp))){
            nVar++;
            if(hasChange)
                nVar++;
            hasChange = true;
        }
        else
            hasChange = false;

        WriteVerboseStatus("HistoTools::get_nVar", "Var, usedBins, up/wasup: " + std::to_string(nVar) + " " + std::to_string(usedBins) + " " + std::to_string(thisUp) + std::to_string(goingUp));

        goingUp = thisUp;
    }

    WriteVerboseStatus("HistoTools::get_nVar", "Out get_nVar:" + std::to_string(nVar));

    return nVar;

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
    if(!nom){
        WriteErrorStatus("HistoTools::CheckHistograms", "The nominal histogram doesn't seem to exist !");
        if(causeCrash) abort();
        isGood = false;
        return isGood;
    }
    if(sh){
        if(!sh->fHistUp){
            WriteErrorStatus("HistoTools::CheckHistograms", "The up variation histogram doesn't seem to exist !");
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }
        if(!sh->fHistDown){
            WriteErrorStatus("HistoTools::CheckHistograms", "The down variation histogram doesn't seem to exist !");
            if(causeCrash) abort();
            isGood = false;
            return isGood;
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
        WriteErrorStatus("HistoTools::CheckHistograms", "The number of bins is found inconsistent ! Please check!");
        if(causeCrash) abort();
        isGood = false;
        return isGood;
    }

    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double lowEdgeNom   = nom->GetBinLowEdge(iBin);
        double lowEdgeUp    = sh->fHistUp->GetBinLowEdge(iBin);
        double lowEdgeDown  = sh->fHistDown->GetBinLowEdge(iBin);

        if( abs(lowEdgeNom-lowEdgeUp)>1e-05 || abs(lowEdgeNom-lowEdgeDown)>1e-05 || abs(lowEdgeDown-lowEdgeUp)>1e-05 ){
            WriteErrorStatus("HistoTools::CheckHistograms", "The bin low edges are not consistent ! Please check !");
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }

        double binWidthNom   = nom->GetBinWidth(iBin);
        double binWidthUp    = sh->fHistUp->GetBinWidth(iBin);
        double binWidthDown  = sh->fHistDown->GetBinWidth(iBin);

        if( abs(binWidthNom-binWidthUp)>1e-05 || abs(binWidthNom-binWidthDown)>1e-05 || abs(binWidthDown-binWidthUp)>1e-05 ){
            WriteErrorStatus("HistoTools::CheckHistograms", "The bin widths are not consistent ! Please check !");
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }
    }

    //
    // 3) Checks the absence on bins with 0 content (for nominal)
    //
    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double content = nom->GetBinContent(iBin);
        if( ( checkNullContent && content<=0 ) || ( !checkNullContent && content<0 ) ){
            std::string temp = nom->GetName();
            WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp + "\", bin " + std::to_string(iBin) + " has 0 content ! Please check");
            WriteInfoStatus("HistoTools::CheckHistograms", "Nominal: " + std::to_string(content));
            isGood = false;
            if(causeCrash){
              abort();
            } else {
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
            WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has negative content ! Please check");
            WriteDebugStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
            isGood = false;
            if(causeCrash) abort();
            else{
                WriteDebugStatus("HistoTools::CheckHistograms", "  => Setting Up to 1e-06");
                sh->fHistUp->SetBinContent(iBin,1e-06);
            }
        }
        if(contentDown<0){
            std::string temp_string = sh->fHistDown->GetName();
            WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has negative content ! Please check");
            WriteDebugStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
            isGood = false;
            if(causeCrash) abort();
            else{
                WriteDebugStatus("HistoTools::CheckHistograms", "  => Setting Up to 1e-06");
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
            WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has weird content ! Please check");
            WriteDebugStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
            isGood = false;
            if(causeCrash) abort();
            // try to fix it, if not aborting
            if(ratioUp!=ratioUp) sh->fHistUp->SetBinContent(iBin,contentNom);
            else return isGood;
        }
        if( (ratioDown!=ratioDown) || (ratioDown < 0) || (abs(ratioDown-1.) >= 100) ){
            std::string temp_string = sh->fHistDown->GetName();
            WriteErrorStatus("HistoTools::CheckHistograms", "In histo \"" + temp_string + "\", bin " + std::to_string(iBin) + " has weird content ! Please check");
            WriteDebugStatus("HistoTools::CheckHistograms", "  Nominal: " + std::to_string(contentNom));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Up: " + std::to_string(contentUp));
            WriteDebugStatus("HistoTools::CheckHistograms", "  Down: " + std::to_string(contentDown));
            isGood = false;
            if(causeCrash) abort();
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
        WriteDebugStatus("HistoTools::CheckHistograms", "Nominal (" + temp1 + "): " + std::to_string(underflowNom));
        WriteDebugStatus("HistoTools::CheckHistograms", "Nominal (" + temp2 + "): " + std::to_string(underflowUp));
        WriteDebugStatus("HistoTools::CheckHistograms", "Nominal (" + temp3 + "): " + std::to_string(underflowDown));
    }

    double overflowNom      = nom -> GetBinContent(NbinsNom+1);
    double overflowUp       = sh->fHistUp -> GetBinContent(NbinsUp+1);
    double overflowDown     = sh->fHistDown -> GetBinContent(NbinsDown+1);
    if( abs(overflowNom)>0 || abs(overflowUp)>0 || abs(overflowDown)>0 ){
        std::string temp1 = nom -> GetName();
        std::string temp2 = sh->fHistUp -> GetName();
        std::string temp3 = sh->fHistDown -> GetName();
        WriteWarningStatus("HistoTools::CheckHistograms", "Overflowflow detected ! This will not be taken into account.");
        WriteDebugStatus("HistoTools::CheckHistograms", "Nominal (" + temp1 + "): " + std::to_string(overflowNom));
        WriteDebugStatus("HistoTools::CheckHistograms", "Nominal (" + temp2 + "): " + std::to_string(overflowUp));
        WriteDebugStatus("HistoTools::CheckHistograms", "Nominal (" + temp3 + "): " + std::to_string(overflowDown));
    }
    return isGood;
}

//_________________________________________________________________________
//
int HistoTools::GetMaxBinWidth(TH1* hist){
    float max = -999;

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ibin++){
        float width = hist->GetBinWidth(ibin);
        if (max < width) max = width;
    }

    return max;
}
