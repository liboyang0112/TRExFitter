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

//_________________________________________________________________________
//
TH1F* HistoTools::TranformHistogramBinning(TH1* originalHist){
    
    //
    // In RooStats, input histogram variable binning is not supported => convert to a constant binning
    // by creating an histogram with the same number of bins but with constant binning between 0 and 1
    //
    const unsigned int nBins = originalHist -> GetNbinsX();
    TH1F *hFinal = new TH1F(originalHist->GetName()+(TString)"_regBin",originalHist->GetTitle(),nBins,0,1);
    hFinal -> SetDirectory(0);
    for(unsigned int iBin = 1; iBin <= nBins; ++iBin){
        hFinal -> SetBinContent(iBin,originalHist->GetBinContent(iBin));
        hFinal -> SetBinError(iBin,originalHist->GetBinError(iBin));
    }
    return hFinal;
}

//_________________________________________________________________________
//
void HistoTools::ManageHistograms( int histOps,  TH1* hNom, TH1* originUp, TH1* originDown,
                                    TH1* &modifiedUp, TH1* &modifiedDown){
    
    //
    // Only function called directly to handle operations on the histograms (symmetrisation and smoothing)
    //
    
    //Sanity checks
    if( histOps % 10 > 2){
        std:: cerr << "\033[1;31m<!> ERROR in HistoTools::ManageHistograms() the operations to perform are not allowed ";
        std:: cerr << "(two different symmetrisations) ! Please check. \033[0m" << std::endl;
        return;
    }
    if( histOps/10 > 9 ){
        std:: cerr << "\033[1;31m<!> ERROR in HistoTools::ManageHistograms() the operations to perform are not allowed ";
        std:: cerr << "(non recognised smoothing) ! Please check. \033[0m" << std::endl;
        return;
    }
    
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
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::ManageHistograms() both up and down variations are empty.  \033[0m" << std::endl;;
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
    modifiedDown -> SetName(originDown->GetName());
    modifiedUp -> SetName(originUp->GetName());
    
    
    //##################################################
    //
    // SECOND STEP: SMOOTHING
    //
    //##################################################
    if(hNom->GetNbinsX()==1){
        if(TtHFitter::DEBUGLEVEL>3){
            std::cout << "In HistoTools::ManageHistograms(): skipping smoothing for systematics on \"" << hNom->GetName() << "\" since just 1 bin." << std::endl;
        }
        return;
    }
    if( ( histOps - ( histOps % 10 ) ) >= SMOOTH && (histOps - ( histOps % 10 ) ) < UNKNOWN ){
        const int smoothingLevel = (histOps - ( histOps % 10 ) ) / 10;
        Smooth_maxVariations( modifiedUp,    hNom,   smoothingLevel );
        Smooth_maxVariations( modifiedDown,  hNom,   smoothingLevel );
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
        std::cerr << "\033[1;31m<!> ERROR in HistoTools::SymmetrizeOneSided() the two histograms for one-sided symmetrisation are the same ... ";
        std::cerr << "Crash is coming ! Please check. \033[0m" << std::endl;
        std::cerr << "      --> Nominal : " << h_nominal -> GetName() << std::endl;
        std::cerr << "      --> Syst    : " << h_syst -> GetName() << std::endl;
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
    
    //Another sanity check: search for negative bins
    for( unsigned int iBin = 1; iBin <= result -> GetNbinsX(); ++iBin ){
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
        tmp1->SetBinError(bin,nom->GetBinError(bin));
    }
    
    delete tmp2;
    delete unit;
    delete nom;

    return tmp1;
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
    
    if(TtHFitter::DEBUGLEVEL>3)
        std::cout << "---: " << tolerance << " " << nVar << std::endl;
    
    //
    // Iterates the smoothing of the systematic histogram until the number a slope changes is lower than "nbins"
    //
    while (nVar > nbins){
        tolerance = tolerance/2.;
        nVar = rebin_getMaxVar(hnom,hsyst,tolerance);
        if(TtHFitter::DEBUGLEVEL>3)
            std::cout << "---: " << tolerance << " " << nVar << std::endl;
        if(tolerance==0){
            std::cout << "Buuuuuuuuuuug: infinite while" << std::endl;
            std::cout << hnom->GetName() << " " << hnom->GetTitle() << " nbins: " << nbins << std::endl;
            break;
        }
    }
    if(TtHFitter::DEBUGLEVEL>3)
        std::cout << "Final: " << tolerance << " " << nVar << std::endl;
    
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
    hsyst->Scale(systIntegral/hsyst->Integral());
    
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
    
    if(TtHFitter::DEBUGLEVEL>3){
        for(int i=1;i<=hsyst->GetNbinsX();i++){
            std::cout << "In: " << hnom->GetBinContent(i) << " " << hsyst->GetBinContent(i) << std::endl;
        }
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
            
            if(TtHFitter::DEBUGLEVEL>3){
                std::cout << thisBin << " " << hnom->GetBinContent(thisBin) << " "<< hnom->GetBinError(thisBin)<<std::endl;
                std::cout << sqrt(cumulErr) << " " << cumulInt << " " << relErr << " " << tolerance << std::endl;
            }
            
        } while (relErr > tolerance && thisBin!=hnom->GetNbinsX() );
        
        if(relErr < tolerance || binLimit.size() == 1){//a group of bins with a sufficiently low stat error has been found, let's add it
            binLimit.push_back(hnom->GetBinCenter(thisBin)+hnom->GetBinWidth(thisBin)/2);
        } else {//no such group of bins has been found: merge with the last found bin
            binLimit.back() = hnom->GetBinCenter(thisBin)+hnom->GetBinWidth(thisBin)/2;
        }
        
        if(TtHFitter::DEBUGLEVEL>3)
            std::cout << "push back bin: " << thisBin <<std::endl;
        
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
    
    if(TtHFitter::DEBUGLEVEL>3)
        std::cout << "get_nVar" << std::endl;
    
    for (int bin=1; bin <=hratio->GetNbinsX(); bin++) {
        
        if(TtHFitter::DEBUGLEVEL>3)
            std::cout << "Bin/content: " << bin << " " << hratio->GetBinContent(bin) << std::endl;
        
        if(hratio->GetBinContent(bin)==0) continue;
        
        prevBin = thisBin;
        thisBin = hratio->GetBinContent(bin);
        
        if(TtHFitter::DEBUGLEVEL>3)
            std::cout << "Prev/this: "<<prevBin << " " << thisBin << std::endl;
        
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
        
        if(TtHFitter::DEBUGLEVEL>3)
            std::cout << "Var, usedBins, up/wasup: " <<nVar << " "<< usedBins << " " << thisUp<<goingUp << std::endl;
        
        goingUp = thisUp;
    }
    
    if(TtHFitter::DEBUGLEVEL>3)
        std::cout << "out get_nVar:" << nVar << std::endl;
    
    return nVar;

}

//_________________________________________________________________________
//
bool HistoTools::HasShape(TH1* hnom, SystematicHist* sh, float threshold){
    
    //Save time
    if(hnom->GetNbinsX()==1) return false;
    
    //If at least one bin is the shape histogram is larger than the threshold, keep the uncertainty
    bool hasShape = false;
    for (unsigned int iBin = 1; iBin < hnom->GetNbinsX(); ++iBin) {
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
        std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms() the nominal histogram doesn't seem to exist !\033[0m" << std::endl;
        if(causeCrash) abort();
        isGood = false;
        return isGood;
    }
    if(sh){
        if(!sh->fHistUp){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms() the up variation histogram doesn't seem to exist !\033[0m" << std::endl;
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }
        if(!sh->fHistDown){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms() the up variation histogram doesn't seem to exist ! \033[0m" << std::endl;
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }
    }
    
    //
    // 2) Checks the binning is the same for the three histograms
    //
    const int NbinsNom  = nom->GetNbinsX();
    const int NbinsUp   = sh->fHistUp->GetNbinsX();
    const int NbinsDown = sh->fHistDown->GetNbinsX();
    if( (NbinsNom != NbinsUp) || (NbinsNom != NbinsDown) || (NbinsUp != NbinsDown) ){
        std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): The number of bins is found inconsistent ! Please check !\033[0m" << std::endl;
        if(causeCrash) abort();
        isGood = false;
        return isGood;
    }
    
    for( unsigned int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double lowEdgeNom   = nom->GetBinLowEdge(iBin);
        double lowEdgeUp    = sh->fHistUp->GetBinLowEdge(iBin);
        double lowEdgeDown  = sh->fHistDown->GetBinLowEdge(iBin);
        
        if( abs(lowEdgeNom-lowEdgeUp)>1e-05 || abs(lowEdgeNom-lowEdgeDown)>1e-05 || abs(lowEdgeDown-lowEdgeUp)>1e-05 ){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): The bin low edges are not consistent ! Please check !\033[0m" << std::endl;
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }

        double binWidthNom   = nom->GetBinWidth(iBin);
        double binWidthUp    = sh->fHistUp->GetBinWidth(iBin);
        double binWidthDown  = sh->fHistDown->GetBinWidth(iBin);
        
        if( abs(binWidthNom-binWidthUp)>1e-05 || abs(binWidthNom-binWidthDown)>1e-05 || abs(binWidthDown-binWidthUp)>1e-05 ){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): The bin widths are not consistent ! Please check !\033[0m" << std::endl;
            if(causeCrash) abort();
            isGood = false;
            return isGood;
        }
    }
    
    //
    // 3) Checks the absence on bins with 0 content (for nominal)
    //
    for( unsigned int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double content = nom->GetBinContent(iBin);
        if( ( checkNullContent && content<=0 ) || ( !checkNullContent && content<0 ) ){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): In histo \""<< nom->GetName() << "\", bin " << iBin << " has 0 content ! Please check !\033[0m" << std::endl;
            std::cout << "Nominal: " << content << std::endl;
            isGood = false;
            if(causeCrash){
                abort();
            } else {
                std:: cerr << "\033[1;33m<!> WARNING in HistoTools::CheckHistograms(): I set the bin content to 1e-05 ! Please check !\033[0m" << std::endl;
                nom -> SetBinContent(iBin,1e-05);
            }
        }
    }
    
    //
    // 4) Check the presence of abnormal bins (content-based)
    //
    for( unsigned int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double contentNom   = nom->GetBinContent(iBin);
        double contentUp    = sh->fHistUp->GetBinContent(iBin);
        double contentDown  = sh->fHistDown->GetBinContent(iBin);
        
        //
        // 4.a) Checks the presence of negative bins for systematic
        //
        if(contentUp<0){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): In histo \""<< sh->fHistUp->GetName() << "\", bin " << iBin << " has negative content ! Please check !\033[0m" << std::endl;
            std::cout << "  Nominal: " << contentNom << std::endl;
            std::cout << "  Up: " << contentUp << std::endl;
            std::cout << "  Down: " << contentDown << std::endl;
            isGood = false;
            if(causeCrash) abort();
            else{
                std::cout << "  => Setting Up to 1e-06" << std::endl;
                sh->fHistUp->SetBinContent(iBin,1e-06);
            }
        }
        if(contentDown<0){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): In histo \"" << sh->fHistDown->GetName() << "\", bin " << iBin << " has negative content ! Please check !\033[0m" << std::endl;
            std::cout << "  Nominal: " << contentNom << std::endl;
            std::cout << "  Up: " << contentUp << std::endl;
            std::cout << "  Down: " << contentDown << std::endl;
            isGood = false;
            if(causeCrash) abort();
            else{
                std::cout << "  => Setting Down to 1e-06" << std::endl;
                sh->fHistDown->SetBinContent(iBin,1e-06);
            }
        }

        //
        // 4.b) Checks that the systematics are not crazy (too large ratio, nan returned, ...)
        //
        double ratioUp   = 0.;
        double ratioDown = 0.;
        if(contentNom != 0 ){
            ratioUp   = contentUp  /contentNom;
            ratioDown = contentDown/contentNom;
        } else if( TMath::Abs(contentUp)>0 || TMath::Abs(contentDown)>0 ) {
            std:: cerr << "\033[1;33m<!> WARNING in HistoTools::CheckHistograms(): In histo \""<< sh->fHistUp->GetName() << "\", bin " << iBin << " has null nominal content but systematics are >0! Hope you know that !\033[0m" << std::endl;
            continue;
        }
        
        // * first check is for nan
        // * second check is for negative values (normally impossible after step 3 and 4.a)
        // * third check is for abnormal systematic/nominal values (ratio > 100)
        if( (ratioUp!=ratioUp) || (ratioUp < 0) || (abs(ratioUp-1.) >= 100) ){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): In histo \""<< sh->fHistUp->GetName() << "\", bin " << iBin << " has weird content ! Please check !\033[0m" << std::endl;
            std::cout << "Nominal: " << contentNom << std::endl;
            std::cout << "Up: " << contentUp << std::endl;
            std::cout << "Down: " << contentDown << std::endl;
            isGood = false;
            if(causeCrash) abort();
            // try to fix it, if not aborting
            if(ratioUp!=ratioUp) sh->fHistUp->SetBinContent(iBin,contentNom);
            else return isGood;
//             return false;
        }
        if( (ratioDown!=ratioDown) || (ratioDown < 0) || (abs(ratioDown-1.) >= 100) ){
            std:: cerr << "\033[1;31m<!> ERROR in HistoTools::CheckHistograms(): In histo \""<< sh->fHistDown->GetName() << "\", bin " << iBin << " has weird content ! Please check !\033[0m" << std::endl;
            std::cout << "Nominal: " << contentNom << std::endl;
            std::cout << "Up: " << contentUp << std::endl;
            std::cout << "Down: " << contentDown << std::endl;
            isGood = false;
            if(causeCrash) abort();
            // try to fix it, if not aborting
            if(ratioDown!=ratioDown) sh->fHistDown->SetBinContent(iBin,contentNom);
            else return isGood;
//             return false;
        }
    }
    
    //
    // 5) Warning level: check the presence of underflow/overflow (ignored in the code)
    //
    double underflowNom     = nom -> GetBinContent(0);
    double underflowUp      = sh->fHistUp -> GetBinContent(0);
    double underflowDown    = sh->fHistDown -> GetBinContent(0);
    if( abs(underflowNom)>0 || abs(underflowUp)>0 || abs(underflowDown)>0 ){
        std:: cerr << "\033[1;33m<!> WARNING in HistoTools::CheckHistograms(): Underflow detected ! This will not be taken into account.\033[0m" << std::endl;
        std::cout << "Nominal (" << nom -> GetName() << "): " << underflowNom << std::endl;
        std::cout << "Up (" << sh->fHistUp -> GetName() << "): " << underflowUp << std::endl;
        std::cout << "Down (" << sh->fHistDown -> GetName() << "): " << underflowDown << std::endl;
    }
    
    double overflowNom      = nom -> GetBinContent(NbinsNom+1);
    double overflowUp       = sh->fHistUp -> GetBinContent(NbinsUp+1);
    double overflowDown     = sh->fHistDown -> GetBinContent(NbinsDown+1);
    if( abs(overflowNom)>0 || abs(overflowUp)>0 || abs(overflowDown)>0 ){
        std:: cerr << "\033[1;33m<!> WARNING in HistoTools::CheckHistograms(): Overflow detected ! This will not be taken into account.\033[0m" << std::endl;
        std::cout << "Nominal (" << nom -> GetName() << "): " << overflowNom << std::endl;
        std::cout << "Up (" << sh->fHistUp -> GetName() << "): " << overflowUp << std::endl;
        std::cout << "Down (" << sh->fHistDown -> GetName() << "): " << overflowDown << std::endl;
    }
    
//     return true;
    return isGood;
}

