#include "TRExFitter/NtupleReader.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/TRExFit.h"

#include "TChain.h"

NtupleReader::NtupleReader(TRExFit* fitter) :
    fFitter(fitter)
{
}

NtupleReader::~NtupleReader() {
}

void NtupleReader::ReadNtuples(){
    WriteInfoStatus("NtupleReader::ReadNtuples", "-------------------------------------------");
    WriteInfoStatus("NtupleReader::ReadNtuples", "Reading ntuples...");
    TH1D* h = nullptr;
    TH1D* hUp = nullptr;
    TH1D* hDown = nullptr;
    std::string variable;
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;
    SampleHist *sh;
    //
    // Import custom functions from .C files
    //
    for(const auto& file : fFitter->fCustomFunctions){
        WriteInfoStatus("NtupleReader::ReadNtuples", "  Loading function from " + file + " ...");
        gROOT->ProcessLineSync((".L "+file+"+").c_str());
    }
    //
    // Loop on regions
    //
    for(int i_ch=0;i_ch<fFitter->fNRegions;i_ch++){
        WriteInfoStatus("NtupleReader::ReadNtuples", "  Region region " + fFitter->fRegions[i_ch]->fName + " ...");
        //
        if(TRExFitter::SPLITHISTOFILES) fFitter->fFiles[i_ch]->cd();
        //
        if(fFitter->fRegions[i_ch]->fBinTransfo != "") fFitter->ComputeBinning(i_ch);
        if(fFitter->fRegions[i_ch]->fCorrVar1 != ""){
            if(fFitter->fRegions[i_ch]->fCorrVar2 == ""){
                WriteWarningStatus("NtupleReader::ReadNtuples", "Only first correlation variable defined, do not read region : " + fFitter->fRegions[i_ch]->fName);
                continue;
            }
            WriteDebugStatus("NtupleReader::ReadNtuples", "calling the function 'DefineVariable(i_ch)'");
            DefineVariable(i_ch);
        }
        else if(fFitter->fRegions[i_ch]->fCorrVar2 != ""){
            WriteWarningStatus("NtupleReader::ReadNtuples", "Only second correlation variable defined, do not read region : " + fFitter->fRegions[i_ch]->fName);
            continue;
        }
        //
        // Loop on samples
        //
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            WriteInfoStatus("NtupleReader::ReadNtuples","    Reading sample " + fFitter->fSamples[i_smp]->fName);
            //
            // eventually skip sample / region combination
            //
            if( Common::FindInStringVector(fFitter->fSamples[i_smp]->fRegions,fFitter->fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            // set variables, selection, weight and paths
            variable      = fFitter->Variable(       fFitter->fRegions[i_ch],fFitter->fSamples[i_smp]);
            fullSelection = fFitter->FullSelection(  fFitter->fRegions[i_ch],fFitter->fSamples[i_smp]);
            fullMCweight  = fFitter->FullWeight(     fFitter->fRegions[i_ch],fFitter->fSamples[i_smp]);
            fullPaths     = fFitter->FullNtuplePaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp]);
            //
            h = nullptr;
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                TH1D* htmp = nullptr;
                if(fFitter->fRegions[i_ch]->fHistoBins.size() > 0){
                    htmp = Common::HistFromNtupleBinArr( fullPaths[i_path],
                                                 variable,
                                                 fFitter->fRegions[i_ch]->fHistoNBinsRebin,
                                                 &(fFitter->fRegions[i_ch]->fHistoBins[0]),
                                                 fullSelection,
                                                 fullMCweight,
                                                 fFitter->fDebugNev);
                }
                else{
                    htmp = Common::HistFromNtuple( fullPaths[i_path],
                                           variable,
                                           fFitter->fRegions[i_ch]->fNbins,
                                           fFitter->fRegions[i_ch]->fXmin,
                                           fFitter->fRegions[i_ch]->fXmax,
                                           fullSelection,
                                           fullMCweight,
                                           fFitter->fDebugNev);
                    //Pre-processing of histograms (rebinning, lumi scaling)
                    if(fFitter->fRegions[i_ch]->fHistoNBinsRebin != -1){
                        htmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin);
                    }
                }
                //
                if(fFitter->fSamples[i_smp]->fNormalizedByTheory && fFitter->fSamples[i_smp]->fType!=Sample::DATA) htmp -> Scale(fFitter->fLumi);
                //
                if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path)  htmp -> Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
                else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
                //
                if(i_path==0) h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fFitter->fRegions[i_ch]->fName.c_str(),fFitter->fSamples[i_smp]->fName.c_str())));
                else h->Add(htmp);

                delete htmp;
            }
            //
            // Save the original histogram
            TH1* h_orig = static_cast<TH1*>(h->Clone( Form("%s_orig",h->GetName()) ));
            //
            // Importing the histogram in TRExFitter
            sh = fFitter->fRegions[i_ch]->SetSampleHist( fFitter->fSamples[i_smp], h );
            sh->fHist_orig.reset(h_orig);
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            //
            //  -----------------------------------
            //
            // read norm factors
            for(int i_norm=0;i_norm<fFitter->fSamples[i_smp]->fNNorm;i_norm++){
                NormFactor *nf = fFitter->fSamples[i_smp]->fNormFactors[i_norm].get();
                //
                // eventually skip norm factor / region combination
                if( nf->fRegions.size()>0 && Common::FindInStringVector(nf->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( nf->fExclude.size()>0 && Common::FindInStringVector(nf->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("NtupleReader::ReadNtuples", "Adding norm " + nf->fName);
                //
                sh->AddNormFactor( nf );
            }

            //
            //  -----------------------------------
            //
            // read shape factors
            for(int i_shape=0;i_shape<fFitter->fSamples[i_smp]->fNShape;i_shape++){
                ShapeFactor *sf = fFitter->fSamples[i_smp]->fShapeFactors[i_shape].get();
                //
                // eventually skip shape factor / region combination
                if( sf->fRegions.size()>0 && Common::FindInStringVector(sf->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( sf->fExclude.size()>0 && Common::FindInStringVector(sf->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("NtupleReader::ReadNtuples", "Adding shape " + sf->fName);
                //
                sh->AddShapeFactor( sf );
            }

            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fFitter->fSamples[i_smp]->fNSyst;i_syst++){
                Systematic * syst = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,
                                                                                       fFitter->fRegions[i_ch]->fName,
                                                                                       fFitter->fSamples[i_smp]->fName)>=0 ) continue;
                //
                WriteDebugStatus("NtupleReader::ReadNtuples", "Adding syst " + syst->fName);
                //
                Region *reg = fFitter->fRegions[i_ch];
                Sample *smp = fFitter->fSamples[i_smp];
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL){
                    SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddOverallSyst(syst->fName,syst->fStoredName,syst->fOverallUp,syst->fOverallDown);
                    syh->fSystematic = syst;
                    syh->fScaleUp = syst->fScaleUp;
                    if(syst->fScaleUpRegions.size()!=0)
                        if(syst->fScaleUpRegions[reg->fName]!=0)
                            syh->fScaleUp *= syst->fScaleUpRegions[reg->fName];
                    syh->fScaleDown = syst->fScaleDown;
                    if(syst->fScaleDownRegions.size()!=0)
                        if(syst->fScaleDownRegions[reg->fName]!=0)
                            syh->fScaleDown *= syst->fScaleDownRegions[reg->fName];
                    continue;
                }
                // if Stat uncertainty on MC sample
                if(syst->fType == Systematic::STAT){
                    SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddStatSyst(syst->fName,syst->fStoredName,syst->fBins[0]);
                    syh->fSystematic = syst;
                    continue;
                }
                // else ...
                if(Common::FindInStringVector(syst->fDummyForSamples,smp->fName)>=0){
                    WriteInfoStatus("NtupleReader::ReadNtuples", "Systematic " + syst->fName + " set as dummy for sample " + smp->fName + " (region " + reg->fName + ")");
                    hUp   = (TH1D*)sh->fHist->Clone(Form("h_%s_%s_%sUp",  reg->fName.c_str(),smp->fName.c_str(),syst->fStoredName.c_str()));
                    hDown = (TH1D*)sh->fHist->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),smp->fName.c_str(),syst->fStoredName.c_str()));
                    SystematicHist *syh = sh->AddHistoSyst(syst->fName,syst->fStoredName,hUp,hDown);
                    syh->fSystematic = syst;
                    syh->fScaleUp = syst->fScaleUp;
                    if(syst->fScaleUpRegions.size()!=0)
                        if(syst->fScaleUpRegions[reg->fName]!=0)
                            syh->fScaleUp *= syst->fScaleUpRegions[reg->fName];
                    syh->fScaleDown = syst->fScaleDown;
                    if(syst->fScaleDownRegions.size()!=0)
                        if(syst->fScaleDownRegions[reg->fName]!=0)
                            syh->fScaleDown *= syst->fScaleDownRegions[reg->fName];
                    continue;
                }
                //
                if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                    if(fFitter->GetSample(syst->fReferenceSample)!=nullptr) smp = fFitter->GetSample(syst->fReferenceSample);
                    else{
//                         WriteWarningStatus("TRExFit::ReadNtuples", "Reference Sample "+syst->fReferenceSample+" for systematic "+syst->fName+", sample "+fSamples[i_smp]->fName+", region "+fRegions[i_ch]->fName+" not found. Ignoring.");
                        WriteErrorStatus("NtupleReader::ReadNtuples", "Reference sample: " + syst->fReferenceSample + " does not exist for region: " + reg->fName + ". Please check this!");
                        WriteErrorStatus("NtupleReader::ReadNtuples", "This probably means that you run over a specific sample, you need to run over the reference sample as well.");
                        WriteErrorStatus("NtupleReader::ReadNtuples", "Ignoring SeparateSample setting.");
                    }
                }
                //
                // Up
                //
                hUp = nullptr;
                if(syst->fHasUpVariation){
                    // set variables, selection, weight and paths
                    fullMCweight  = fFitter->FullWeight(     fFitter->fRegions[i_ch],fFitter->fSamples[i_smp],syst,true);
                    fullPaths     = fFitter->FullNtuplePaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp],syst,true);
                    WriteDebugStatus("NtupleReader::ReadNtuples", "  Syst Up full weight: " + fullMCweight);
                    for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                        TH1D* htmp = nullptr;
                        if(reg->fHistoBins.size() > 0){
                            htmp = Common::HistFromNtupleBinArr(fullPaths[i_path],
                                                        variable,
                                                        reg->fHistoNBinsRebin,
                                                        &reg->fHistoBins[0],
                                                        fullSelection,
                                                        fullMCweight,
                                                        fFitter->fDebugNev);
                        }
                        else{
                            htmp = Common::HistFromNtuple( fullPaths[i_path],
                                                   variable,
                                                   reg->fNbins,
                                                   reg->fXmin,
                                                   reg->fXmax,
                                                   fullSelection,
                                                   fullMCweight,
                                                   fFitter->fDebugNev);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1){
                                htmp->Rebin(reg->fHistoNBinsRebin);
                            }
                        }
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fFitter->fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);

                        //
                        // Importing histogram in TRExFitter
                        if(i_path==0){
                            hUp = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fFitter->fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str())));
                        }
                        else hUp->Add(htmp);
                        delete htmp;
                    } // end loop over files


                    // BW
                    // pulled this out of the file loop to apply it only to the fully constructed histogram insead of file by file

                    // obtain relative variation and apply it to proper sample
                    // & try to keep also the same total relative variation
                    if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr){
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist.get();
                        TH1* hnom = reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();

                        // Protection added: fix empty bins before starting to divide and multiply

                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin< hUp->GetNbinsX()+2;i_bin++) if(hUp ->GetBinContent(i_bin)<=1e-6) hUp ->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) hUp ->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                        //
                        double relVar   = hUp->Integral(0,hUp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);

                        // get copies with no error
                        auto hrefTmp = Common::GetHistCopyNoError(href);
                        auto hnomTmp = Common::GetHistCopyNoError(hnom);
                        hUp->Divide(   hrefTmp.get() );
                        hUp->Multiply( hnomTmp.get() );
                        double newVar   = hUp->Integral(0,hUp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                        if( syst->fKeepReferenceOverallVar && std::abs(relVar-1) > 0.0001 && std::abs(newVar) > 0.0001) hUp->Scale( relVar / newVar );
                    }
                    // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                    else if (syst->fReferenceSample!="" && syst->fReferenceSample != fFitter->fSamples[i_smp]->fName && syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr) {
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist.get();
                        TH1* href_up = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistUp.get();
                        TH1* hnom = reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();

                        // Protection added: fix empty bins before starting to divide and multiply

                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin< hUp->GetNbinsX()+2;i_bin++) if( hUp->GetBinContent(i_bin)<=1e-6) hUp->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) hUp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6

                        // Formula: UpHisto = [1+(up-nom)/nom-(DataUp-Data)/Data]*nom = up+nom+DataUp/Data*nom
                        TH1* href_up_Tmp = static_cast<TH1*>(href_up->Clone(Form("%s_Tmp", href_up->GetName())));
                        // get copies with no error
                        auto hrefTmp = Common::GetHistCopyNoError(href);
                        auto hnomTmp = Common::GetHistCopyNoError(hnom);
                        href_up_Tmp->Divide(hrefTmp.get());
                        href_up_Tmp->Multiply(hnomTmp.get());
                        hUp->Add(hnomTmp.get());
                        auto href_up_TmpNoError = Common::GetHistCopyNoError(href_up_Tmp);
                        hUp->Add(href_up_TmpNoError.get(),-1);

                        delete href_up_Tmp;// it's a clone, and it's the purpose of clones to die
                    }


                //--------------------------------------

                }  // end Up variation
                //
                // Down
                //
                hDown = nullptr;
                if(syst->fHasDownVariation){
                    fullMCweight  = fFitter->FullWeight(     fFitter->fRegions[i_ch],fFitter->fSamples[i_smp],syst,false);
                    fullPaths     = fFitter->FullNtuplePaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp],syst,false);
                    for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                        TH1D* htmp = nullptr;
                        if(reg->fHistoBins.size() > 0){
                            htmp = Common::HistFromNtupleBinArr(fullPaths[i_path],
                                                        variable,
                                                        reg->fHistoNBinsRebin,
                                                        &reg->fHistoBins[0],
                                                        fullSelection,
                                                        fullMCweight,
                                                        fFitter->fDebugNev);
                        }
                        else{
                            htmp = Common::HistFromNtuple( fullPaths[i_path],
                                                   variable,
                                                   reg->fNbins,
                                                   reg->fXmin,
                                                   reg->fXmax,
                                                   fullSelection,
                                                   fullMCweight,
                                                   fFitter->fDebugNev);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1){
                                htmp->Rebin(reg->fHistoNBinsRebin);
                            }
                        }
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fFitter->fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);

                        //
                        // Importing histogram in TRExFitter
                        if(i_path==0){
                            hDown = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fFitter->fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str())));
                        }
                        else hDown->Add(htmp);
                        delete htmp;
                    }  // end loop over files

                    // BW
                    // pulled this out of the file loop to apply it only to the fully constructed histogram insead of file by file
                    //
                    // obtain relative variation and apply it to proper sample
                    // & try to keep also the same total relative variation
                    if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr){
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist.get();
                        TH1* hnom = reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();

                        // Protection added: fix empty bins before starting to divide and multiply

                        for(int i_bin=0;i_bin< href->GetNbinsX()+2;i_bin++) if(href ->GetBinContent(i_bin)<=1e-6) href ->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<hDown->GetNbinsX()+2;i_bin++) if(hDown->GetBinContent(i_bin)<=1e-6) hDown->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin< href->GetNbinsX()+2;i_bin++) if(href ->GetBinContent(i_bin)<=1e-6) hDown->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                        //
                        double relVar   = hDown->Integral(0,hDown->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                        hDown->Divide(   href );
                        hDown->Multiply( hnom );
                        double newVar   = hDown->Integral(0,hDown->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                        if( syst->fKeepReferenceOverallVar && std::abs(relVar-1) > 0.0001 && std::abs(newVar-1) > 0.0001) hDown->Scale( relVar / newVar );
                    }
                    // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                    else if (syst->fReferenceSample!="" && syst->fReferenceSample != fFitter->fSamples[i_smp]->fName && syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr) {
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist.get();
                        TH1* href_down = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistDown.get();
                        TH1* hnom = reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();

                        // Protection added: fix empty bins before starting to divide and multiply

                        for(int i_bin=0;i_bin<href ->GetNbinsX()+2;i_bin++) if(href ->GetBinContent(i_bin)<=1e-6) href ->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<hDown->GetNbinsX()+2;i_bin++) if(hDown->GetBinContent(i_bin)<=1e-6) hDown->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<href ->GetNbinsX()+2;i_bin++) if(href ->GetBinContent(i_bin)<=1e-6) hDown->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6

                        // Formula: UpHisto = [1+(down-nom)/nom-(DataDown-Data)/Data]*nom = down+nom+DataDown/Data*nom
                        TH1* href_down_Tmp = (TH1*) href_down->Clone(Form("%s_Tmp", href_down->GetName()));
                        href_down_Tmp->Divide(href);
                        href_down_Tmp->Multiply(hnom);
                        hDown->Add(hnom);
                        hDown->Add(href_down_Tmp,-1);

                        delete href_down_Tmp;// it's a clone, and it's the purpose of clones to die
                    }

                }  // end Down variation
                //
                if(hUp==nullptr)   hUp   = static_cast<TH1D*>(reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get());
                if(hDown==nullptr) hDown = static_cast<TH1D*>(reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get());
                //
                SystematicHist *syh = sh->AddHistoSyst(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fName,
                                                       fFitter->fSamples[i_smp]->fSystematics[i_syst]->fStoredName,hUp,hDown);
                syh->fSystematic = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                syh->fScaleUp = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0)
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName]!=0)
                        syh->fScaleUp *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName];
                syh->fScaleDown = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0)
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName]!=0)
                        syh->fScaleDown *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName];
            }
        }
    }
}

void NtupleReader::DefineVariable(int regIter){
    TH1::StatOverflows(true);  //////  What is the defaut in root for this ???
    WriteDebugStatus("NtupleReader::DefineVariable", "//////// --------");
    WriteDebugStatus("NtupleReader::DefineVariable", "// DEBUG CORR VAR");
    TH1* h1 = new TH1D("h1","h1",1,-2000.,1000.);
    TH1* h2 = new TH1D("h2","h2",1,-2000.,1000.);
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;

    // copy of NtupleReading function.
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        WriteDebugStatus("NtupleReader::DefineVariable", "Processing sample : " + fFitter->fSamples[i_smp]->fName);
        if(fFitter->fSamples[i_smp]->fType==Sample::DATA) continue;
        if(Common::FindInStringVector(fFitter->fSamples[i_smp]->fRegions,fFitter->fRegions[regIter]->fName)<0 ) continue;
        WriteDebugStatus("NtupleReader::DefineVariable", " -> is used in the considered region");
        //
        // set selection, weight and paths (no variables)
        fullSelection = fFitter->FullSelection(  fFitter->fRegions[regIter],fFitter->fSamples[i_smp]);
        fullMCweight  = fFitter->FullWeight(     fFitter->fRegions[regIter],fFitter->fSamples[i_smp]);
        fullPaths     = fFitter->FullNtuplePaths(fFitter->fRegions[regIter],fFitter->fSamples[i_smp]);
        //
        for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
            WriteDebugStatus("TRExFit::DefineVariable", " -> Retrieving : " + fFitter->fRegions[regIter]->fCorrVar1 +
                                                        " w/ weight " + fullMCweight + "*" + fullSelection  +
                                                        " from " +  fullPaths[i_path]);
            TH1* htmp1 = new TH1D("htmp1","htmp1",1,-2000.,1000.);
            TH1* htmp2 = new TH1D("htmp2","htmp2",1,-2000.,1000.);
            TChain *t = new TChain();
            t->Add(fullPaths[i_path].c_str());
            t->Draw(Form("%s>>htmp1",fFitter->fRegions[regIter]->fCorrVar1.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            t->Draw( Form("%s>>htmp2",fFitter->fRegions[regIter]->fCorrVar2.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            delete t;
            //
            if(fFitter->fSamples[i_smp]->fType!=Sample::DATA && fFitter->fSamples[i_smp]->fNormalizedByTheory) htmp1 -> Scale(fFitter->fLumi);
            if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path)  htmp1 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
            else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1) htmp1 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
            //
            if(fFitter->fSamples[i_smp]->fType!=Sample::DATA && fFitter->fSamples[i_smp]->fNormalizedByTheory) htmp2 -> Scale(fFitter->fLumi);
            if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path)  htmp2 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
            else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1) htmp2 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
            //
            h1->Add(htmp1);
            h2->Add(htmp2);
            delete htmp1;
            delete htmp2;
        }
    }
    double mean1 = h1->GetMean();
    double rms1 = h1->GetRMS();
    double mean2 = h2->GetMean();
    double rms2 = h2->GetRMS();
    WriteDebugStatus("NtupleReader::DefineVariable", "the new variable : ( ( (" + fFitter->fRegions[regIter]->fCorrVar1 + ") - " + std::to_string(mean1) + " )*( (" +fFitter->fRegions[regIter]->fCorrVar2 + ")-" + std::to_string(mean2) + " ) )/( " + std::to_string(rms1) + " * " + std::to_string(rms2) + ")");
    fFitter->fRegions[regIter]->fVariable = Form("( ( (%s)-%f )*( (%s)-%f ) )/( %f * %f )",fFitter->fRegions[regIter]->fCorrVar1.c_str(),mean1,fFitter->fRegions[regIter]->fCorrVar2.c_str(),mean2,rms1,rms2);
    TH1::StatOverflows(false);  //////  What is the defaut in root for this ???
    delete h1;
    delete h2;
}
