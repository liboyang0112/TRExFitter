#include "TRExFitter/HistoReader.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExFit.h"

#include "TH1.h"

HistoReader::HistoReader(TRExFit* fitter) :
    fFitter(fitter)
{
}

HistoReader::~HistoReader() {
}

void HistoReader::ReadHistograms(){
    TH1* hUp = nullptr;
    TH1* hDown = nullptr;

    //
    // Loop on regions and samples
    //
    for(int i_ch=0;i_ch<fFitter->fNRegions;i_ch++){
        WriteInfoStatus("HistoReader::ReadHistograms", "  Region " + fFitter->fRegions[i_ch]->fName + " ...");
        //
        if(TRExFitter::SPLITHISTOFILES) fFitter->fFiles[i_ch]->cd();
        //
        if(fFitter->fRegions[i_ch]->fBinTransfo != "") fFitter->ComputeBinning(i_ch);
        // first we must read the DATA samples
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            if(fFitter->fSamples[i_smp]->fType!=Sample::DATA) continue;
            WriteDebugStatus("HistoReader::ReadHistograms", "  Reading DATA sample " + fFitter->fSamples[i_smp]->fName);
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fFitter->fSamples[i_smp]->fRegions,fFitter->fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            std::vector<std::string> fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp]);

            TH1D* h = ReadSingleHistogram(fullPaths, nullptr, i_ch, i_smp, true, false); // is nominal and not MC
            //
            // Save the original histogram
            TH1* h_orig = static_cast<TH1*>(h->Clone( Form("%s_orig",h->GetName()) ));
            //
            // Importing the histogram in TRExFitter
            SampleHist *sh = fFitter->fRegions[i_ch]->SetSampleHist( fFitter->fSamples[i_smp], h );
            sh->fHist_orig.reset(h_orig);
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            // in fact DATA can be used for systs that have SubtractRefSampleVar: TRUE
            for(int i_syst=0;i_syst<fFitter->fSamples[i_smp]->fNSyst;i_syst++){
                Systematic *syst = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                // only relevant for systs that have this sample as reference
                if (!syst->fSubtractRefSampleVar || syst->fReferenceSample != fFitter->fSamples[i_smp]->fName) continue;

                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fFitter->fRegions[i_ch]->fName, fFitter->fSamples[i_smp]->fName)>=0 ) continue;
                //
                WriteDebugStatus("HistoReader::ReadHistograms", "Adding syst " + syst->fName);
                //

                //
                // Up
                //
                hUp = nullptr;
                if(syst->fHasUpVariation){
                    fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp],syst,true);
                    hUp = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, true, false); // is up variation and not MC
                }
                //
                // Down
                //
                hDown = nullptr;
                if(syst->fHasDownVariation){
                    fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp],syst,false);
                    hDown = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, false, false); // is down variation and not MC
                }
                //
                if(hUp==nullptr){
                    hUp   = fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();
                }
                if(hDown==nullptr){
                    hDown = fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();
                }
                //
                SystematicHist *syh = sh->AddHistoSyst(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fName,
                                                       fFitter->fSamples[i_smp]->fSystematics[i_syst]->fStoredName,hUp,hDown);
                syh->fSystematic = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                syh->fScaleUp = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0) {
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[fFitter->fRegions[i_ch]->fName]!=0){
                        syh->fScaleUp *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[fFitter->fRegions[i_ch]->fName];
                    }
                }
                syh->fScaleDown = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0) {
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[fFitter->fRegions[i_ch]->fName]!=0){
                        syh->fScaleDown *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[fFitter->fRegions[i_ch]->fName];
                    }
                }
            }
        }

        // then we can read the other samples
        std::set < std::string > files_names;
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            if(fFitter->fSamples[i_smp]->fType==Sample::DATA) continue;
            WriteDebugStatus("HistoReader::ReadHistograms", "  Reading " + fFitter->fSamples[i_smp]->fName);
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fFitter->fSamples[i_smp]->fRegions,fFitter->fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            std::vector<std::string>fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch],fFitter->fSamples[i_smp]);
            for (const auto& ipath : fullPaths){
                files_names.insert(ipath);
            }
            TH1D* h = ReadSingleHistogram(fullPaths, nullptr, i_ch, i_smp, true, true); // is MC
            //
            // Save the original histogram
            TH1* h_orig = static_cast<TH1*>(h->Clone( Form("%s_orig",h->GetName())));
            //
            // Importing the histogram in TRExFitter
            SampleHist *sh = fFitter->fRegions[i_ch]->SetSampleHist(fFitter->fSamples[i_smp], h);
            sh->fHist_orig.reset(h_orig);
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            // end here no systematics allowed (e.g. generally for GHOST samples)
            if (!(fFitter->fSamples[i_smp]->fUseSystematics)) continue;

            //
            //  -----------------------------------
            //
            // read norm factors
            for(int i_norm=0;i_norm<fFitter->fSamples[i_smp]->fNNorm;i_norm++){
                NormFactor *nf = fFitter->fSamples[i_smp]->fNormFactors[i_norm].get();
                //
                // eventually skip systematic / region combination
                if( nf->fRegions.size()>0 && FindInStringVector(nf->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( nf->fExclude.size()>0 && FindInStringVector(nf->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("HistoReader::ReadHistograms", "Adding norm " + nf->fName);
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
                // eventually skip systematic / region combination
                if( sf->fRegions.size()>0 && FindInStringVector(sf->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( sf->fExclude.size()>0 && FindInStringVector(sf->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("HistoReader::ReadHistograms", "Adding shape " + sf->fName);
                //
                sh->AddShapeFactor(sf);
            }

            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fFitter->fSamples[i_smp]->fNSyst;i_syst++){
                Systematic *syst = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,
                                                                                       fFitter->fRegions[i_ch]->fName,
                                                                                       fFitter->fSamples[i_smp]->fName)>=0) continue;
                //
                WriteDebugStatus("HistoReader::ReadHistograms", "Adding syst " + syst->fName);
                //
                Region *reg = fFitter->fRegions[i_ch];
                Sample *smp = fFitter->fSamples[i_smp];
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL) {
                    SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddOverallSyst(syst->fName,
                                                                                         syst->fStoredName,
                                                                                         syst->fOverallUp,
                                                                                         syst->fOverallDown);
                    syh->fSystematic = syst;
                    syh->fScaleUp = syst->fScaleUp;
                    if(syst->fScaleUpRegions.size()!=0) {
                        if(syst->fScaleUpRegions[reg->fName]!=0) {
                            syh->fScaleUp *= syst->fScaleUpRegions[reg->fName];
                        }
                    }
                    syh->fScaleDown = syst->fScaleDown;
                    if(syst->fScaleDownRegions.size()!=0) {
                        if(syst->fScaleDownRegions[reg->fName]!=0) {
                            syh->fScaleDown *= syst->fScaleDownRegions[reg->fName];
                        }
                    }
                    continue;
                }
                // else ...
                //
                if(syst->fReferenceSample!="") smp = fFitter->GetSample(syst->fReferenceSample);
                //
                // Up
                //
                hUp = nullptr;
                if(syst->fHasUpVariation){
                    fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch],
                                                            fFitter->fSamples[i_smp],
                                                            syst,
                                                            true);
                    for (const auto& ipath : fullPaths){
                        files_names.insert(ipath);
                    }
                    hUp = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, true, true); // isUp and isMC
                }
                //
                // Down
                //
                hDown = nullptr;
                if(syst->fHasDownVariation){
                    fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch],
                                                            fFitter->fSamples[i_smp],
                                                            syst,
                                                            false);
                    for (const auto& ipath : fullPaths){
                        files_names.insert(ipath);
                    }
                    hDown = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, false, true); // isUp and isMC
                }
                //
                if(hUp==nullptr)   hUp   = static_cast<TH1D*>(reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get());
                if(hDown==nullptr) hDown = static_cast<TH1D*>(reg->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get());
                //
                SystematicHist *syh = sh->AddHistoSyst(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fName,
                                                       fFitter->fSamples[i_smp]->fSystematics[i_syst]->fStoredName,
                                                       hUp,
                                                       hDown);
                syh->fSystematic = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                syh->fScaleUp = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0) {
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName]!=0) {
                        syh->fScaleUp *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName];
                    }
                }
                syh->fScaleDown = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0) {
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName]!=0) {
                        syh->fScaleDown *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName];
                    }
                }
            }
            //closing the files for this sample
            CloseFiles(files_names);
            files_names.clear();
        }
    }
}

TH1D* HistoReader::ReadSingleHistogram(const std::vector<std::string>& fullPaths,
                                       Systematic* syst,
                                       int i_ch,
                                       int i_smp,
                                       bool isUp,
                                       bool isMC) {
    TH1D* h = nullptr;
    for(unsigned int i_path = 0; i_path < fullPaths.size(); ++i_path){
        std::unique_ptr<TH1> htmp = HistFromFile( fullPaths.at(i_path) );
        if (!htmp) {
            WriteErrorStatus("HistoReader::ReadSingleHistogram", "Histo pointer is nullptr, cannot continue running the code");
            exit(EXIT_FAILURE);
        }
        //Pre-processing of histograms (rebinning, lumi scaling)
        if(fFitter->fRegions[i_ch]->fHistoBins.size() > 0){
            const char *hname = htmp->GetName();
            std::unique_ptr<TH1> tmp_copy(static_cast<TH1D*>(htmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin, "tmp_copy", &(fFitter->fRegions[i_ch]->fHistoBins[0]))));
            htmp.reset(tmp_copy.release());
            htmp->SetName(hname);
            if(TRExFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(htmp.get());
        }
        else if(fFitter->fRegions[i_ch]->fHistoNBinsRebin != -1) {
            htmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin);
        }

        if (isMC){
            if(fFitter->fSamples[i_smp]->fNormalizedByTheory){
                htmp->Scale(fFitter->fLumi);
            }
        }

        if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path){
             htmp->Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
        }
        else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1){
            htmp->Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
        }

        if (isMC && syst != nullptr){
            // obtain relative variation and apply it to proper sample
            // & try to keep also the same total relative variation
            if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                // check if the reference sample exists
                if (fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample) == nullptr){
                    WriteErrorStatus("HistoReader::ReadSingleHistogram", "Reference sample: " + syst->fReferenceSample + " does not exist for region: " + fFitter->fRegions[i_ch]->fName + ". Please check this!");
                    exit(EXIT_FAILURE);
                }
                TH1* href = fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->fHist.get();
                TH1* hnom = fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();
                // Protection added: fix empty bins before starting to divide and multiply
                for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++){
                    if(href->GetBinContent(i_bin)<=1e-6){
                        href->SetBinContent(i_bin,1e-6);
                    }
                }
                for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++){
                    if(htmp->GetBinContent(i_bin)<=1e-6){
                        htmp->SetBinContent(i_bin,1e-6);
                    }
                }
                for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++){
                    if(href->GetBinContent(i_bin)<=1e-6){
                        htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                    }
                }
                //
                const double relVar = htmp->Integral(0,htmp->GetNbinsX()+1) /
                     href->Integral(0,href->GetNbinsX()+1);

                // get copies with no error
                auto hrefTmp = GetHistCopyNoError(href);
                auto hnomTmp = GetHistCopyNoError(hnom);
                htmp->Divide(   hrefTmp.get() );
                htmp->Multiply( hnomTmp.get() );
                const double newVar = htmp->Integral(0,htmp->GetNbinsX()+1) /
                    hnom->Integral(0,hnom->GetNbinsX()+1);
                if( syst->fKeepReferenceOverallVar && (std::abs(relVar-1) > 0.0001) &&
                    (std::abs(newVar-1) > 0.0001)){
                    htmp->Scale( relVar / newVar );
                }
            }
            // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
            else if (syst->fReferenceSample!="" && syst->fSubtractRefSampleVar) {
                // check if the reference sample exists
                if (fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample) == nullptr){
                    WriteErrorStatus("HistoReader::ReadSingleHistogram", "Reference sample: " + syst->fReferenceSample + " does not exist for region: " + fFitter->fRegions[i_ch]->fName + ". Please check this!");
                    exit(EXIT_FAILURE);
                }
                TH1* href = fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->fHist.get();
                TH1* href_upDown = nullptr;
                if (isUp){
                    href_upDown = fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->
                        GetSystematic(syst->fName)->fHistUp.get();
                } else {
                    href_upDown = fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->
                        GetSystematic(syst->fName)->fHistDown.get();
                }
                TH1* hnom = fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();
                // Protection added: fix empty bins before starting to divide and multiply
                for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++){
                    if(href->GetBinContent(i_bin)<=1e-6){
                        href->SetBinContent(i_bin,1e-6);
                    }
                }
                for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++){
                    if(htmp->GetBinContent(i_bin)<=1e-6){
                        htmp->SetBinContent(i_bin,1e-6);
                    }
                }
                for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++){
                    if(href->GetBinContent(i_bin)<=1e-6){
                        htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                    }
                }
                // Formula: UpHisto = [1+(up-nom)/nom-(DataUp-Data)/Data]*nom = up+nom+DataUp/Data*nom
                TH1* href_upDown_Tmp = static_cast<TH1*>(href_upDown->Clone(
                    Form("%s_Tmp", href_upDown->GetName())));
                // get copies with no error
                auto hrefTmp = GetHistCopyNoError(href);
                auto hnomTmp = GetHistCopyNoError(hnom);
                href_upDown_Tmp->Divide(hrefTmp.get());
                href_upDown_Tmp->Multiply(hnomTmp.get());
                htmp->Add(hnomTmp.get());
                auto href_upDown_TmpNoErr = GetHistCopyNoError(href_upDown_Tmp);
                htmp->Add(href_upDown_TmpNoErr.get(),-1);

                delete href_upDown_Tmp;// it's a clone, and it's the purpose of clones to die
            }
        }

        if(i_path == 0){
            if (syst == nullptr){ // is nominal
                h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fFitter->fRegions[i_ch]->fName.c_str(),
                    fFitter->fSamples[i_smp]->fName.c_str())));
            } else { // is syst
                if (isUp){ // up variation
                    h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s_%sUp",fFitter->fRegions[i_ch]->fName.c_str(),
                        fFitter->fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str())));
                } else { // down variation
                    h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s_%sDown",fFitter->fRegions[i_ch]->fName.c_str(),
                        fFitter->fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str())));
                }
            }
        }
        else{
            h->Add(htmp.get());
        }
    }
    return h;
}


