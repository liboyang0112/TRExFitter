#include "TRExFitter/HistoReader.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExFit.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

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

void HistoReader::ReadTRExProducedHistograms() {
    std::string fileName = "";
    std::string fileNameBootstrap = "";
    SampleHist *sh = nullptr;
    SystematicHist *syh = nullptr;
    std::string regionName;
    std::string sampleName;
    std::string normName;
    std::string shapeName;
    //
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    if(singleOutputFile){
        if(fFitter->fInputFolder!="") fileName = fFitter->fInputFolder           + fFitter->fInputName + "_histos.root";
        else                          fileName = fFitter->fName + "/Histograms/" + fFitter->fInputName + "_histos.root";
        // Bootstrap
        if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0){
            fileName = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fFitter->fBootstrapIdx));
        }
        WriteInfoStatus("HistoReader::ReadTRExProducedHistograms", "-----------------------------");
        WriteInfoStatus("HistoReader::ReadTRExProducedHistograms", "Reading histograms from file " + fileName + " ...");
    }
    //
    std::vector< TH2* > histPrun;
    std::unique_ptr<TFile> filePrun = nullptr;
    if( fFitter->fKeepPruning ){
        filePrun = std::unique_ptr<TFile>(new TFile( (fFitter->fName+"/Pruning.root").c_str() ));
        if(!filePrun) fFitter->fKeepPruning = false;
    }
    //
    // when we multply/divide by or subtract/add other samples, need to add systematics on the other samples
    Systematic * tmpsyst;
    for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
        if(!fFitter->fSamples[i_smp]->fUseSystematics) continue;
        if(fFitter->fSamples[i_smp]->fDivideBy!=""){
            Sample* smp = fFitter->GetSample(fFitter->fSamples[i_smp]->fDivideBy);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fFitter->fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", " The sample " + fFitter->fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst].get())[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fFitter->fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fFitter->fNSyst;j_syst++){
                       if(fFitter->fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if( FindInStringVector(fFitter->fSystematics[j_syst]->fSamples,
                                                 fFitter->fSamples[i_smp]->fName)<0 ) fFitter->fSystematics[j_syst]->fSamples.push_back(fFitter->fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
        if(fFitter->fSamples[i_smp]->fMultiplyBy!=""){
            Sample* smp = fFitter->GetSample(fFitter->fSamples[i_smp]->fMultiplyBy);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fFitter->fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", " The sample " + fFitter->fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst].get())[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fFitter->fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fFitter->fNSyst;j_syst++){
                       if(fFitter->fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if( FindInStringVector(fFitter->fSystematics[j_syst]->fSamples,
                                                 fFitter->fSamples[i_smp]->fName)<0 ) fFitter->fSystematics[j_syst]->fSamples.push_back(fFitter->fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
        for(const auto& sample : fFitter->fSamples[i_smp]->fSubtractSamples){
            Sample* smp = fFitter->GetSample(sample);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fFitter->fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", " The sample " + fFitter->fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst].get())[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fFitter->fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fFitter->fNSyst;j_syst++){
                       if(fFitter->fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if(FindInStringVector(fFitter->fSystematics[j_syst]->fSamples,
                                                fFitter->fSamples[i_smp]->fName)<0 ) fFitter->fSystematics[j_syst]->fSamples.push_back(fFitter->fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
        for(const auto& sample : fFitter->fSamples[i_smp]->fAddSamples){
            Sample* smp = fFitter->GetSample(sample);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fFitter->fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", " The sample " + fFitter->fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst].get())[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fFitter->fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fFitter->fNSyst;j_syst++){
                       if(fFitter->fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if(FindInStringVector(fFitter->fSystematics[j_syst]->fSamples,
                                                fFitter->fSamples[i_smp]->fName)<0 ) fFitter->fSystematics[j_syst]->fSamples.push_back(fFitter->fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
    }
    //
    // Syst for morphing samples inherited from nominal sample
    if (fFitter->fPropagateSystsForMorphing){
        for(const auto& par : fFitter->fMorphParams){
            double nominalValue = 0.;
            for(const auto& norm : fFitter->fNormFactors){
                if(norm->fName == par) nominalValue = norm->fNominal;
            }
            Sample *smpNominal = nullptr;
            for(const auto& smp : fFitter->fSamples){
                if(!smp->fIsMorph[par]) continue;
                if(smp->fMorphValue[par] == nominalValue){ // FIXME: eventually add something to flag a sample as nominal for morphing
                    smpNominal = smp;
                    break;
                }
            }
            for(const auto& smp : fFitter->fSamples){
                if(!smp->fIsMorph[par]) continue;
                if(smp == smpNominal) continue;
                for(const auto& syst : smpNominal->fSystematics){
                    smp->AddSystematic(syst.get());
                }
            }
        }
    }
    //
    // fSystFromSample
    for(const auto& smp : fFitter->fSamples){
        if(smp->fSystFromSample!=""){
            Sample *smpReference = nullptr;
            for(const auto& smp2 : fFitter->fSamples){
                if(smp2->fName == smp->fName) continue;
                if(smp2->fName == smp->fSystFromSample){
                    smpReference = smp2;
                    break;
                }
            }
            if(smpReference!=nullptr){
                for(const auto& syst : smpReference->fSystematics){
                    smp->AddSystematic(syst.get());
                }
            }
        }
    }
    //
    for(int i_ch=0;i_ch<fFitter->fNRegions;i_ch++){
        if(fFitter->fKeepPruning){
            histPrun.push_back( static_cast<TH2*>(filePrun->Get( Form("h_prun_%s_toSave", fFitter->fRegions[i_ch]->fName.c_str()))));
        }
        regionName = fFitter->fRegions[i_ch]->fName;
        WriteDebugStatus("HistoReader::ReadTRExProducedHistograms","  Reading region " + regionName);
        //
        if(!singleOutputFile){
            if(fFitter->fInputFolder!="") fileName = fFitter->fInputFolder           + fFitter->fInputName + "_" + regionName + "_histos.root";
            else                          fileName = fFitter->fName + "/Histograms/" + fFitter->fInputName + "_" + regionName + "_histos.root";
            // Bootstrap
            if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0){
                if(fFitter->fBootstrapSyst == "") {
                    fileName = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fFitter->fBootstrapIdx));
                } else {
                    fileNameBootstrap = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fFitter->fBootstrapIdx));
                }
            }
            WriteInfoStatus("HistoReader::ReadTRExProducedHistograms", "-----------------------------");
            WriteInfoStatus("HistoReader::ReadTRExProducedHistograms", "Reading histograms from file " + fileName + " ...");
        }
        //
        for(int i_smp=0;i_smp<fFitter->fNSamples;i_smp++){
            //
            // eventually skip sample / region combination
            //
            if(FindInStringVector(fFitter->fSamples[i_smp]->fRegions,regionName)<0 && fFitter->fSamples[i_smp]->fName.find("customAsimov_")==std::string::npos ) continue;
            //
            sampleName = fFitter->fSamples[i_smp]->fName;
            WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "    Reading sample " + sampleName);
            fFitter->fRegions[i_ch]->SetSampleHist(fFitter->fSamples[i_smp],regionName+"_"+sampleName,fileName);
            sh = fFitter->fRegions[i_ch]->GetSampleHist(sampleName);
            if(sh==nullptr) continue;
            //
            // separate gammas -> Add systematic
            if(fFitter->fSamples[i_smp]->fSeparateGammas){
                std::string systName = "stat_"+fFitter->fSamples[i_smp]->fName;
                std::string systStoredName = systName;
                WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "adding separate gammas as SHAPE systematic " + systName);
                SystematicHist *syh_tmp = sh->AddHistoSyst(systName,
                                                           systStoredName,
                                                           Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                           Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                           0
                                                        );
                if(syh_tmp==nullptr){
                    WriteWarningStatus("HistoReader::ReadTRExProducedHistograms", "No histogram found for separate gamma, but may be you will create it right now.");
                }
                else{
                    Systematic *gamma = nullptr;
                    if(FindInStringVector(fFitter->fSystematicNames,systName)>=0) gamma = fFitter->fSystematics[FindInStringVector(fFitter->fSystematicNames,systName)];  //GetSystematic(systName);
                    if(gamma==nullptr) gamma = fFitter->NewSystematic(systName);
                    WriteDebugStatus("TRExFit::ReadHistos", "adding separate gammas as SHAPE systematic " + systName);
                    gamma->fType = Systematic::SHAPE;
                    gamma->fRegions.clear();
                    gamma->fRegions.push_back(fFitter->fRegions[i_ch]->fName);
                    syh_tmp->fSystematic = gamma;
                    gamma->fNuisanceParameter = gamma->fName;
                    TRExFitter::NPMAP[gamma->fName] = gamma->fNuisanceParameter;
                }
            }
            //
            // norm factors
            for(int i_norm=0;i_norm<fFitter->fSamples[i_smp]->fNNorm;i_norm++){
                //
                // eventually skip norm factor / region combination
                if(fFitter->fSamples[i_smp]->fNormFactors[i_norm]->fRegions.size()>0 && FindInStringVector(fFitter->fSamples[i_smp]->fNormFactors[i_norm]->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if(fFitter->fSamples[i_smp]->fNormFactors[i_norm]->fExclude.size()>0 && FindInStringVector(fFitter->fSamples[i_smp]->fNormFactors[i_norm]->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                normName = fFitter->fSamples[i_smp]->fNormFactors[i_norm]->fName;
                WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "      Reading norm " + normName);
                // norm only
                sh->AddNormFactor(fFitter->fSamples[i_smp]->fNormFactors[i_norm].get());
            }
            //
            // shape factors
            for(int i_shape=0;i_shape<fFitter->fSamples[i_smp]->fNShape;i_shape++){
                //
                // eventually skip shape factor / region combination
                if( fFitter->fSamples[i_smp]->fShapeFactors[i_shape]->fRegions.size()>0 && FindInStringVector(fFitter->fSamples[i_smp]->fShapeFactors[i_shape]->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( fFitter->fSamples[i_smp]->fShapeFactors[i_shape]->fExclude.size()>0 && FindInStringVector(fFitter->fSamples[i_smp]->fShapeFactors[i_shape]->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                shapeName = fFitter->fSamples[i_smp]->fShapeFactors[i_shape]->fName;
                WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "      Reading shape " + shapeName);
                // shape only
                sh->AddShapeFactor(fFitter->fSamples[i_smp]->fShapeFactors[i_shape].get());
            }
            //
            // systematics
            for(int i_syst=0;i_syst<fFitter->fSamples[i_smp]->fNSyst;i_syst++){
                //
                // eventually skip systematic / region combination
                if( fFitter->fSamples[i_smp]->fSystematics[i_syst]->fRegions.size()>0 && FindInStringVector(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if( fFitter->fSamples[i_smp]->fSystematics[i_syst]->fExclude.size()>0 && FindInStringVector(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                if( fFitter->fSamples[i_smp]->fSystematics[i_syst]->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fExcludeRegionSample,fFitter->fRegions[i_ch]->fName, fFitter->fSamples[i_smp]->fName)>=0 ) continue;
                //
                std::string systName       = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fName;
                std::string systStoredName = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fStoredName; // if no StoredName specified in the config, this should be == fName
                //
                // eventually skip systematics if pruned
                int binContent(0);
                if( fFitter->fKeepPruning && histPrun[i_ch]!=nullptr ){
                    const int xbin = histPrun[i_ch]->GetXaxis()->FindBin( sampleName.c_str() ); // sample
                    const int ybin = histPrun[i_ch]->GetYaxis()->FindBin( systName.c_str() ); // syst
                    const int bin = histPrun[i_ch]->GetBin(xbin,ybin);
                    binContent = histPrun[i_ch]->GetBinContent(bin);
                    if( binContent <= -4 || binContent == -1 || binContent >= 3 ){
                        WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "SKIPPING systematic " + systName);
                        continue;
                    }
                }
                WriteDebugStatus("HistoReader::ReadTRExProducedHistograms", "      Reading syst " + systName);
                // norm only
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    if( fFitter->fKeepPruning ){
                        if( binContent == -2 || binContent == 2 ) continue;
                    }
                    syh = sh->AddOverallSyst(systName,
                                             systStoredName,
                                             fFitter->fSamples[i_smp]->fSystematics[i_syst]->fOverallUp,
                                             fFitter->fSamples[i_smp]->fSystematics[i_syst]->fOverallDown);
                }
                // histo syst
                else{
                    int pruned = 0;
                    if(fFitter->fKeepPruning){
                        if(binContent==1 || binContent==-2) pruned = 1;
                        if(binContent==2 || binContent==-3) pruned = 2;
                    }
                    if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0 && fFitter->fBootstrapSyst == systName ){
                        syh = sh->AddHistoSyst(systName,
                                               systStoredName,
                                               Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileNameBootstrap,
                                               Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileNameBootstrap,
                                               pruned
                                              );
                    }
                    else{
                        syh = sh->AddHistoSyst(systName,
                                               systStoredName,
                                               Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                               Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                               pruned
                                              );
                    }
                    if(syh==nullptr){
                        if (!pruned) WriteWarningStatus("HistoReader::ReadTRExProducedHistograms", "No syst histo found for syst " + systName + ", sample " + sampleName + ", region " + regionName);
                        continue;
                    }
                }
                // for both
                syh->fSystematic = fFitter->fSamples[i_smp]->fSystematics[i_syst].get();
                syh->fHistoNameShapeUp   = Form("%s_%s_%s_Shape_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                syh->fHistoNameShapeDown = Form("%s_%s_%s_Shape_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0 && fFitter->fBootstrapSyst == systName ){
                    syh->fFileNameShapeUp    = fileNameBootstrap;
                    syh->fFileNameShapeDown  = fileNameBootstrap;
                }
                else{
                    syh->fFileNameShapeUp    = fileName;
                    syh->fFileNameShapeDown  = fileName;
                }
                syh->fScaleUp = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0){
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName]!=0){
                        syh->fScaleUp *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName];
                    }
                }
                syh->fScaleDown = fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0){
                    if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName]!=0){
                        syh->fScaleDown *= fFitter->fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName];
                    }
                }
                //
                if(fFitter->fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    syh->fNormUp   *= syh->fScaleUp;
                    syh->fNormDown *= syh->fScaleDown;
                }
            }
        }
    }

    if (filePrun != nullptr){
        filePrun->Close();
    }
}
