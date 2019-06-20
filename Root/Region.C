// Class include
#include "TRExFitter/Region.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/TRExPlot.h"

// ROOT includes
#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TPaveText.h"
#include "TSystem.h"

// c++ includes
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

// -------------------------------------------------------------------------------------------------
// class Region

//__________________________________________________________________________________
//
Region::Region(string name){
    fName = name;
    fLabel = name;
    fTexLabel = "";
    fShortLabel = name;
    fRegionType = CONTROL;
    fNBkg = 0;
    fNSig = 0;
    fNSyst = 0;
    fNNorm = 0;
    fNShape = 0;
    fHasData = false;
    fHasSig = false;
    fNSamples = 0;
    fUseStatErr = false;
    fHistoBins = 0;
    fHistoNBinsRebin = -1;
    fHistoBinsPost = 0;
    fHistoNBinsRebinPost = -1;
    fYTitle = "";
    fYmaxScale = 0;
    fYmax = 0;
    fYmin = 0;
    fRatioYmax = 2.;
    fRatioYmin = 0;
    fRatioYmaxPostFit = 1.5;
    fRatioYminPostFit = 0.5;
    fRatioYtitle = "";
    fRatioType = "DATA/MC";

    int canvasWidth = 600;
    int canvasHeight = 700;
    string cName = "c_"+fName;
    if(TRExFitter::OPTION["CanvasWidth"]!=0)  canvasWidth  = TRExFitter::OPTION["CanvasWidth"];
    if(TRExFitter::OPTION["CanvasHeight"]!=0) canvasHeight = TRExFitter::OPTION["CanvasHeight"];
    fPlotPreFit = new TRExPlot(cName,canvasWidth,canvasHeight,TRExFitter::NORATIO);
    fPlotPreFit->fShowYields = TRExFitter::SHOWYIELDS;
    cName = "c_"+fName+"_postFit";
    fPlotPostFit = new TRExPlot(cName,canvasWidth,canvasHeight,TRExFitter::NORATIO);
    fPlotPostFit->fShowYields = TRExFitter::SHOWYIELDS;

    fSampleHists.clear();
    fSamples.clear();
    fSystematics.clear();
    fNormFactors.clear();
    fShapeFactors.clear();

    fAlternativeVariables.clear();
    fAlternativeSelections.clear();

    fIntCode_overall = 4;
    fIntCode_shape = 0;

    fCorrVar1 = "";
    fCorrVar2 = "";

    fFitName = "";
    fFitType = TRExFit::SPLUSB;
    fPOI = "";
    fFitLabel = "";

    fSelection = "1";
    fMCweight = "1";

    fLogScale = false;
    fSkipSmoothing = false;

    fBinWidth = 0;

    fLumiScale = 1.;

    fBlindingThreshold = -1;

    fRegionDataType = REALDATA;

    fBinTransfo = "";
    fTransfoDzBkg = 0.;
    fTransfoDzSig = 0.;
    fTransfoFzBkg = 0.;
    fTransfoFzSig = 0.;
    fTransfoJpar1 = 0.;
    fTransfoJpar2 = 0.;
    fTransfoJpar3 = 0.;
    fAutoBinBkgsInSig.clear();

    fATLASlabel = "Internal";

    fSuffix = "";
    fGroup = "";

    fBlindedBins = nullptr;
    fKeepPrefitBlindedBins = false;
    fGetChi2 = 0;

    fDropBins.clear();

    fBinLabels.clear();

    fChi2val = -1;
    fNDF = -1;
    fChi2prob = -1;

    fUseGammaPulls = false;

    fTot = nullptr;

    fLabelX = -1;
    fLabelY = -1;
    fLegendX1 = -1;
    fLegendX2 = -1;
    fLegendY = -1;

    fLegendNColumns = 2;
}

//__________________________________________________________________________________
//
Region::~Region(){
    if(fPlotPreFit) delete fPlotPreFit;
    if(fPlotPostFit) delete fPlotPostFit;

    for(unsigned int i = 0; i < fSampleHists.size(); ++i){
        if(fSampleHists[i]){
            delete fSampleHists[i];
        }
    }
    fSampleHists.clear();
}

//__________________________________________________________________________________
//
SampleHist* Region::SetSampleHist(Sample *sample, string histoName, string fileName){
    fSampleHists.push_back(new SampleHist( sample, histoName, fileName ));
    if(sample->fType==Sample::DATA){
        fHasData = true;
        fData = fSampleHists[fNSamples];
    }
    else if(sample->fType==Sample::SIGNAL){
        fHasSig = true;
        fSig[fNSig] = fSampleHists[fNSamples];
        fNSig ++;
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg[fNBkg] = fSampleHists[fNSamples];
        fNBkg ++;
    }
    else if(sample->fType==Sample::GHOST){
        WriteDebugStatus("Region::SetSampleHist", "Adding GHOST sample.");
    }
    else{
        WriteErrorStatus("Region::SetSampleHist", "SampleType not supported.");
    }
    fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists[fNSamples]->fRegionName = fName;
    fSampleHists[fNSamples]->fRegionLabel = fLabel;
    fSampleHists[fNSamples]->fFitName = fFitName;
    fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
    fNSamples++;
    return fSampleHists[fNSamples-1];
}

//__________________________________________________________________________________
//
SampleHist* Region::SetSampleHist(Sample *sample, TH1* hist ){
    fSampleHists.push_back(new SampleHist( sample, hist ));
    if(sample->fType==Sample::DATA){
        fHasData = true;
        fData = fSampleHists[fNSamples];
    }
    else if(sample->fType==Sample::SIGNAL){
        fHasSig = true;
        fSig[fNSig] = fSampleHists[fNSamples];
        fNSig ++;
    }
    else if(sample->fType==Sample::BACKGROUND){
        fBkg[fNBkg] = fSampleHists[fNSamples];
        fNBkg ++;
    }
    else if(sample->fType==Sample::GHOST){
        WriteDebugStatus("Region::SetSampleHist", "Adding GHOST sample.");
    }
    else{
        WriteErrorStatus("Region::SetSampleHist", "SampleType not supported.");
    }
    fSampleHists[fNSamples]->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists[fNSamples]->fRegionName = fName;
    fSampleHists[fNSamples]->fRegionLabel = fLabel;
    fSampleHists[fNSamples]->fFitName = fFitName;
    fSampleHists[fNSamples]->fVariableTitle = fVariableTitle;
    fNSamples++;
    return fSampleHists[fNSamples-1];
}

//__________________________________________________________________________________
//
void Region::AddSample(Sample* sample){
    fSamples.push_back(sample);
    fNSamples++;
}

//__________________________________________________________________________________
//
void Region::AddSystematic(Systematic *syst){
    fSystematics.push_back(syst);
    fNSyst++;
}

//__________________________________________________________________________________
//
void Region::SetBinning(int N, double *bins){
    fNbins = fHistoNBinsRebin = N;
    fHistoBins = new double[N+1];
    for(int i=0; i<=N; ++i) fHistoBins[i] = bins[i];
}

//__________________________________________________________________________________
//
void Region::Rebin(int N){
    fNbins = fHistoNBinsRebin = N;
}

//__________________________________________________________________________________
//
void Region::SetRebinning(int N, double *bins){
    fHistoNBinsRebinPost = N;
    fHistoBinsPost = new double[N+1];
    for(int i=0; i<=N; ++i) fHistoBinsPost[i] = bins[i];
}

//__________________________________________________________________________________
//
void Region::SetRegionType( RegionType type ){
    fRegionType = type;
}

//__________________________________________________________________________________
//
void Region::SetRegionDataType( DataType type ){
    fRegionDataType = type;
}

//__________________________________________________________________________________
//
SampleHist* Region::GetSampleHist(const std::string &sampleName) const{
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSampleHists[i_smp]->fName == sampleName) return fSampleHists[i_smp];
    }
    return nullptr;
}

//__________________________________________________________________________________
//
void Region::BuildPreFitErrorHist(){
    WriteInfoStatus("Region::BuildPreFitErrorHist", "Building pre-fit plot for region " + fName + " ...");
    //
    double yieldNominal(0.), yieldUp(0.), yieldDown(0.);
    double diffUp(0.), diffDown(0.);
    //
    fSystNames.clear();
    fNpNames.clear();
    std::map<std::string,std::string> origShapeSystName;
    std::map<std::string,bool> systIsThere;
    systIsThere.clear();
    TH1* hUp = nullptr;
    TH1* hDown = nullptr;
    string systName = "";
    string systNuisPar = "";
    SystematicHist *sh = nullptr;

    //
    // Collect all the systematics on all the samples
    // (in parallel also colect all the nuisance parameters)
    //
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType == Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType == Sample::GHOST) continue;

        //
        // non-SHAPE Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
            systName = fSampleHists[i]->fSyst[i_syst]->fName;
            if(fSampleHists[i]->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
            if(fSampleHists[i]->fSyst[i_syst]->fSystematic!=nullptr)
                systNuisPar = fSampleHists[i]->fSyst[i_syst]->fSystematic->fNuisanceParameter;
            if(FindInStringVector(fNpNames,systNuisPar)<0){
                fNpNames.push_back(systNuisPar);
            }
        }

        //
        // SHAPE Systematics
        // (but remove the MC-stat ones (for samples with fSeparateGammas)!!
        //
        for(int i_syst=0;i_syst<fSampleHists[i]->fNSyst;i_syst++){
            systName = fSampleHists[i]->fSyst[i_syst]->fName;
            if(fSampleHists[i]->fSyst[i_syst]->fSystematic->fType!=Systematic::SHAPE) continue;
            if(systName.find("stat_")!=std::string::npos) continue;
            for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
                std::string gammaName = Form("shape_%s_%s_bin_%d",systName.c_str(),fName.c_str(),i_bin-1);
                if(!systIsThere[gammaName]){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                    origShapeSystName[gammaName] = systName;
                    TRExFitter::NPMAP[gammaName] = gammaName;
                }
                if(fSampleHists[i]->fSyst[i_syst]->fSystematic!=nullptr)
                    systNuisPar = gammaName;
                if(FindInStringVector(fNpNames,systNuisPar)<0){
                    fNpNames.push_back(systNuisPar);
                }
            }
        }
    }

    //
    // Build pre-fit error hists (for each sample and each systematic)
    //

    // - loop on samples
    //
    for(int i=0;i<fNSamples;i++){

        // skip data
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;

        WriteDebugStatus("Region::BuildPreFitErrorHist", "  Sample: " + fSampleHists[i]->fName);

        // - loop on systematics
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            WriteDebugStatus("Region::BuildPreFitErrorHist", "    Systematic: " + fSystNames[i_syst]);
            systName = fSystNames[i_syst];

            // get SystematicHist
            sh = fSampleHists[i]->GetSystematic(systName);

            // hack: add a systematic hist if not there... FIXME
            if(sh==nullptr){
                fSampleHists[i]->AddHistoSyst(systName,fSampleHists[i]->fHist,fSampleHists[i]->fHist);
                sh = fSampleHists[i]->GetSystematic(systName);
                // initialize the up and down variation histograms
                // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
                sh->fHistUp   = (TH1*)fSampleHists[i]->fHist->Clone(Form("%s_%s_Up",  fSampleHists[i]->fHist->GetName(),systName.c_str()));
                sh->fHistDown = (TH1*)fSampleHists[i]->fHist->Clone(Form("%s_%s_Down",fSampleHists[i]->fHist->GetName(),systName.c_str()));
            }
            //
            // - loop on bins
            for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
                diffUp = 0.;
                diffDown = 0.;
                hUp = nullptr;
                hDown = nullptr;
                yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
                // if it's a systematic (NB: skip Norm-Factors!!)
                if(fSampleHists[i]->HasSyst(fSystNames[i_syst])){
                    hUp   = sh->fHistUp;
                    hDown = sh->fHistDown;
                    if(hUp!=nullptr)    yieldUp     = hUp  ->GetBinContent(i_bin);
                    else            yieldUp     = yieldNominal;
                    if(hDown!=nullptr)  yieldDown   = hDown->GetBinContent(i_bin);
                    else            yieldDown   = yieldNominal;
                    diffUp   += yieldUp   - yieldNominal;
                    diffDown += yieldDown - yieldNominal;
                }
                // for shape systematics
                if(origShapeSystName[systName]!="" && systName.find(Form("%s_bin_%d",fName.c_str(),i_bin-1))!=std::string::npos){
                    SystematicHist *shOrig = fSampleHists[i]->GetSystematic(origShapeSystName[systName]);
                    if(!shOrig) continue;
                    yieldUp   = shOrig->fHistUp->GetBinContent(i_bin);
                    yieldDown = 2*yieldNominal - yieldUp;
                    sh->fHistUp  ->SetBinContent(i_bin,yieldUp);
                    sh->fHistDown->SetBinContent(i_bin,yieldDown);
                    diffUp   += yieldUp   - yieldNominal;
                    diffDown += yieldDown - yieldNominal;
                }
                WriteDebugStatus("Region::BuildPreFitErrorHist", "        Bin " + std::to_string(i_bin) + ":  " + " \t +" + std::to_string(100*diffUp/yieldNominal)
                 + "%\t " + std::to_string(100*diffDown/yieldNominal) + "%");
            }
        }
    }

    // at this point all the sample-by-sample pre-fit variation histograms should be filled
    //
    // Now build the total prediction variations, for each systematic
    // - loop on systematics
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        systName = fSystNames[i_syst];

        // initialize the tot variation hists
        fTotUp[i_syst]   = (TH1*)fTot->Clone(Form("h_%s_tot_%s_Up",  fName.c_str(), systName.c_str()));
        fTotDown[i_syst] = (TH1*)fTot->Clone(Form("h_%s_tot_%s_Down",fName.c_str(), systName.c_str()));
        // - loop on bins
        for(int i_bin=1;i_bin<fTot->GetNbinsX()+1;i_bin++){
            diffUp = 0.;
            diffDown = 0.;
            // - loop on samples
            for(int i=0;i<fNSamples;i++){
                //
                // scale according to NormFactors
                float scale = 1.;
                for(unsigned int i_nf=0;i_nf<fSampleHists[i]->fSample->fNormFactors.size();i_nf++){
                    NormFactor *nf = fSampleHists[i]->fSample->fNormFactors[i_nf];
                    // if this norm factor is a morphing one
                    if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                        std::string formula = TRExFitter::SYSTMAP[nf->fName];
                        std::string name = TRExFitter::NPMAP[nf->fName];
                        WriteDebugStatus("Region::BuildPreFitErrorHist", "formula: " +formula);
                        WriteDebugStatus("Region::BuildPreFitErrorHist", "name: " +name);
                        std::vector < std::pair < std::string,std::vector<double> > > nameS;
                        if(nf->fName.find("morph_")!=std::string::npos){
                            nameS.push_back(std::make_pair(name,std::vector<double>{double(nf->fNominal),double(nf->fMin),double(nf->fMax)}));
                        }
                        else{
                            nameS = processString(name);
                        }
                        std::vector <double> nfNominalvec;
                        for (unsigned int j = 0; j<nameS.size(); j++){
                            formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                            nfNominalvec.push_back(nameS[j].second[0]);
                        }
                        double *nfNominal = nfNominalvec.data();
                        WriteDebugStatus("Region::BuildPreFitErrorHist", "formula: " +formula);
                        for(unsigned int j = 0; j<nameS.size(); j++){
                            WriteDebugStatus("Region::BuildPreFitErrorHist", "nfNominal["+std::to_string(j)+"]: "+std::to_string(nfNominal[j]));
                        }
                        TFormula f_morph ("f_morph",formula.c_str());
                        scale *= f_morph.EvalPar(nfNominal,nullptr);
                    }
                    else {
                        scale *= fSampleHists[i]->fSample->fNormFactors[i_nf]->fNominal;
                    }
                }
                //
                // skip data
                if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
                if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
                if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;
                // get SystematicHist
                sh = fSampleHists[i]->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);
                diffUp   += (sh->fHistUp  ->GetBinContent(i_bin) - yieldNominal)*scale;
                diffDown += (sh->fHistDown->GetBinContent(i_bin) - yieldNominal)*scale;
            }
            // add the proper bin content to the variation hists
            fTotUp[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }

    WriteDebugStatus("Region::BuildPreFitErrorHist", "----");

    //
    // build the vectors of variations (sum histograms for systematics with the same NP)
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    for(size_t i_syst=0;i_syst< fSystNames.size(); ++i_syst){
        // look for all the already stored systematics, to find if one had the same NP
        bool found = false;
        for(size_t j_syst=0;j_syst<i_syst;++j_syst){
            if(TRExFitter::NPMAP[fSystNames[i_syst]]==TRExFitter::NPMAP[fSystNames[j_syst]]){
                found = true;
                const int whichsyst = FindInStringVector(fNpNames,TRExFitter::NPMAP[fSystNames[i_syst]]);
                auto h_diff_up   = std::unique_ptr<TH1> (static_cast<TH1*> (h_up[whichsyst]  ->Clone(Form("%s_%s","clone_",h_up[whichsyst]  ->GetName()))));
                auto h_diff_down = std::unique_ptr<TH1> (static_cast<TH1*> (h_down[whichsyst]->Clone(Form("%s_%s","clone_",h_down[whichsyst]->GetName()))));
                h_diff_up  ->Add(fTotUp[  i_syst],fTot,1,-1);
                h_diff_down->Add(fTotDown[i_syst],fTot,1,-1);
                h_up[   FindInStringVector(fNpNames,TRExFitter::NPMAP[fSystNames[i_syst]])]->Add(h_diff_up.get());
                h_down[ FindInStringVector(fNpNames,TRExFitter::NPMAP[fSystNames[i_syst]])]->Add(h_diff_down.get());
                break;
            }
        }
        if(found) continue;
        //
        // if shape-only systematics to show, normalize each syst variation
        if(TRExFitter::OPTION["ShapeOnlySystBand"]>0){
            fTotUp[i_syst]  ->Scale(fTot->Integral()/fTotUp[i_syst]  ->Integral());
            fTotDown[i_syst]->Scale(fTot->Integral()/fTotDown[i_syst]->Integral());
        }
        h_up.  push_back( fTotUp[i_syst]   );
        h_down.push_back( fTotDown[i_syst] );
    }
    fErr = BuildTotError( fTot, h_up, h_down, fNpNames );
    fErr->SetName("g_totErr");
    // at this point fTot and fErr should be ready

    //
    // Goodness of pre-fit
    if(fGetChi2){
        // remove blinded bins
        if (!fHasData){
            WriteWarningStatus("Region::BuildPreFitErrorHist", "Data histogram is nullptr, cannot calculate Chi2 agreement.");
            WriteWarningStatus("Region::BuildPreFitErrorHist", "Maybe you do not have data sample defined?");
            return;
        }
        auto h_data = std::unique_ptr<TH1>(static_cast<TH1*>(fData->fHist->Clone()));
        if(fBlindedBins!=nullptr){
            for(int i_bin=1;i_bin<=h_data->GetNbinsX();i_bin++){
                if(fBlindedBins->GetBinContent(i_bin)>0) h_data->SetBinContent(i_bin,-1);
                if(find(fDropBins.begin(),fDropBins.end(),i_bin-1)!=fDropBins.end()) h_data->SetBinContent(i_bin,-1);
            }
        }
        if(fGetChi2==1) fNpNames.clear();
        std::pair<double,int> res = GetChi2Test( h_data.get(), fTot, h_up, fNpNames );
        fChi2val = res.first;
        fNDF = res.second;
        fChi2prob = ROOT::Math::chisquared_cdf_c( res.first, res.second);
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- PRE-FIT AGREEMENT EVALUATION -----------------------");
        if(fGetChi2==1)
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- -------- STAT-ONLY --------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "--- REGION " + fName + ":");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  chi2        = " + std::to_string(fChi2val));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  ndof        = " + std::to_string(fNDF));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "  probability = " + std::to_string(fChi2prob));
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPreFitErrorHist", "----------------------- ---------------------------- -----------------------");
    }
}

//__________________________________________________________________________________
//
TRExPlot* Region::DrawPreFit(const std::vector<int>& canvasSize, string opt){

    TRExPlot *p{};
    if (canvasSize.size() == 0){
        p = fPlotPreFit;
    } else {
        p = new TRExPlot(("c_"+fName).c_str(), canvasSize.at(0), canvasSize.at(1),TRExFitter::NORATIO);
        p->fShowYields = TRExFitter::SHOWYIELDS;
    }
    p->SetXaxisRange(fXaxisRange);
    if(fYmaxScale==0) p->SetYmaxScale(1.8);
    else              p->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) p->fYmax = fYmax;
    if(fYmin!=0) p->fYmin = fYmin;
    p->fRatioYmax = fRatioYmax;
    p->fRatioYmin = fRatioYmin;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") p->SetYaxis(fYTitle);
    //
    // For 4-top-style plots
    if(TRExFitter::OPTION["FourTopStyle"]>0){
        if(fRegionType==CONTROL) p->AddLabel("#font[62]{Control Region}");
        else if(fRegionType==SIGNAL) p->AddLabel("#font[62]{Signal Region}");
        else if(fRegionType==VALIDATION) p->AddLabel("#font[62]{Validation Region}");
        else p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        p->AddLabel("#font[52]{Pre-fit}");
    }
    //
    // old-style plots
    else{
        p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        if(TRExFitter::OPTION["NoPrePostFit"]==0) p->AddLabel("Pre-Fit");
    }
    //
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    p->fLegendNColumns = fLegendNColumns;
    if(fBlindingThreshold>=0) p->SetBinBlinding(true,fBlindingThreshold,fBlindingType);

    if(fBinLabels.size() && ((int)fBinLabels.size()==fNbins)) {
      for(int i_bin=0; i_bin<fNbins; i_bin++) {
        p->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
      }
    }

    //
    // build h_tot
    //
    fTot = nullptr;
    string title;
    TH1* h = nullptr;
    if(fHasData && opt.find("blind")==string::npos) p->SetData(fData->fHist,fData->fSample->fTitle);
    for(int i=0;i<fNSig;i++){
        title = fSig[i]->fSample->fTitle;
        if(fSig[i]->fSample->fGroup != "") title = fSig[i]->fSample->fGroup;
        h = (TH1*)fSig[i]->fHist->Clone();
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if(!fSig[i]->fSample->fUseMCStat && !fSig[i]->fSample->fSeparateGammas){
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++) h->SetBinError(i_bin,0.);
        }
        // else still check if the value is reasonable
        else{
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++){
                if(h->GetBinError(i_bin)>10*h->GetBinContent(i_bin)) h->SetBinError(i_bin,h->GetBinContent(i_bin));
            }
        }
        // scale it according to NormFactors
        for(unsigned int i_nf=0;i_nf<fSig[i]->fSample->fNormFactors.size();i_nf++){
            NormFactor *nf = fSig[i]->fSample->fNormFactors[i_nf];
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                std::string formula = TRExFitter::SYSTMAP[nf->fName];
                std::string name = TRExFitter::NPMAP[nf->fName];
                WriteDebugStatus("Region::DrawPreFit", "formula: " +formula);
                WriteDebugStatus("Region::DrawPreFit", "name: " +name);
                std::vector < std::pair < std::string,std::vector<double> > > nameS;
                if(nf->fName.find("morph_")!=std::string::npos){
                    nameS.push_back(std::make_pair(name,std::vector<double>{double(nf->fNominal),double(nf->fMin),double(nf->fMax)}));
                }
                else{
                    nameS = processString(name);
                }
                std::vector <double> nfNominalvec;
                for (unsigned int j = 0; j<nameS.size(); j++){
                    formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                    nfNominalvec.push_back(nameS[j].second[0]);
                }
                double *nfNominal = nfNominalvec.data();
                WriteDebugStatus("Region::DrawPreFit", "formula: " +formula);
                for(unsigned int j = 0; j<nameS.size(); j++){
                    WriteDebugStatus("Region::DrawPreFit", "nfNominal["+std::to_string(j)+"]: "+std::to_string(nfNominal[j]));
                }
                TFormula f_morph("f_morph",formula.c_str());
                float scale = 1;
                scale = f_morph.EvalPar(nfNominal,nullptr);
                h->Scale(scale);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fSig[i]->fSample->fName + " by " + std::to_string(scale));
            }
            else{
                h->Scale(nf->fNominal);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fSig[i]->fSample->fName + " by " + std::to_string(fSig[i]->fSample->fNormFactors[i_nf]->fNominal));
            }
        }
        if(TRExFitter::SHOWSTACKSIG)   p->AddSignal(    h,title);
        if(TRExFitter::SHOWNORMSIG){
            if( (TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
             || !TRExFitter::OPTION["NormSigSRonly"] )
                p->AddNormSignal(h,title);
        }
        else{
            if(TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(h,title);
        }
        if(TRExFitter::SHOWOVERLAYSIG) p->AddOverSignal(h,title);
        if(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG){
            if(fTot==nullptr) fTot = (TH1*)h->Clone("h_tot");
            else              fTot->Add(h);
        }
    }
    for(int i=0;i<fNBkg;i++){
        title = fBkg[i]->fSample->fTitle;
        if(fBkg[i]->fSample->fGroup != "") title = fBkg[i]->fSample->fGroup;
        h = (TH1*)fBkg[i]->fHist->Clone();
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if(!fBkg[i]->fSample->fUseMCStat && !fBkg[i]->fSample->fSeparateGammas){
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++) h->SetBinError(i_bin,0.);
        }
        // else still check if the value is reasonable
        else{
            for(int i_bin=0;i_bin<h->GetNbinsX()+2;i_bin++){
                if(h->GetBinError(i_bin)>10*h->GetBinContent(i_bin)) h->SetBinError(i_bin,h->GetBinContent(i_bin));
            }
        }
        // scale it according to NormFactors
        for(unsigned int i_nf=0;i_nf<fBkg[i]->fSample->fNormFactors.size();i_nf++){
            NormFactor *nf = fBkg[i]->fSample->fNormFactors[i_nf];
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                std::string formula = TRExFitter::SYSTMAP[nf->fName];
                std::string name = TRExFitter::NPMAP[nf->fName];
                WriteDebugStatus("Region::DrawPreFit", "formula: " +formula);
                WriteDebugStatus("Region::DrawPreFit", "name: " +name);
                std::vector < std::pair < std::string,std::vector<double> > > nameS;
                if(nf->fName.find("morph_")!=std::string::npos){
                    nameS.push_back(std::make_pair(name,std::vector<double>{double(nf->fNominal),double(nf->fMin),double(nf->fMax)}));
                }
                else{
                    nameS = processString(name);
                }
                std::vector <double> nfNominalvec;
                for (unsigned int j = 0; j<nameS.size(); j++){
                    formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                    nfNominalvec.push_back(nameS[j].second[0]);
                }
                double *nfNominal = nfNominalvec.data();
                WriteDebugStatus("Region::DrawPreFit", "formula: " +formula);
                for(unsigned int j = 0; j<nameS.size(); j++){
                    WriteDebugStatus("Region::DrawPreFit", "nfNominal["+std::to_string(j)+"]: "+std::to_string(nfNominal[j]));
                }
                TFormula f_morph("f_morph",formula.c_str());
                float scale = 1;
                scale = f_morph.EvalPar(nfNominal,nullptr);
                h->Scale(scale);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fBkg[i]->fSample->fName + " by " + std::to_string(scale));
            }
            else{
                h->Scale(nf->fNominal);
                WriteDebugStatus("Region::DrawPreFit", nf->fName + " => Scaling " + fBkg[i]->fSample->fName + " by " + std::to_string(fBkg[i]->fSample->fNormFactors[i_nf]->fNominal));
            }
        }
        p->AddBackground(h,title);
        if(fTot==nullptr) fTot = (TH1*)h->Clone("h_tot");
        else          fTot->Add(h);
    }

    //
    // set error to 0 if no MCstat
    //
    if(!fUseStatErr){
        for(int i_bin=1;i_bin<=fTot->GetNbinsX();i_bin++){
            fTot->SetBinError(i_bin,0);
        }
    }
    else{
        for(int i_bin=0;i_bin<fTot->GetNbinsX()+2;i_bin++){
            if(fTot->GetBinError(i_bin)>10*fTot->GetBinContent(i_bin)) fTot->SetBinError(i_bin,fTot->GetBinContent(i_bin));
        }
    }

    p->SetTotBkg((TH1*)fTot);
    p->BlindData();
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);
    fBlindedBins = p->h_blinding;

    //
    // Computes the uncertainty bands arround the h_tot histogram
    //
    BuildPreFitErrorHist();

    //
    // Print chi2 info
    //
    if(fGetChi2 && TRExFitter::SHOWCHI2)  p->SetChi2KS(fChi2prob,-1,fChi2val,fNDF);

    //
    // Sets the last ingredients in the TRExPlot object
    //
    p->SetTotBkgAsym(fErr);
    p->fATLASlabel = fATLASlabel;
    p->fRatioYtitle = fRatioYtitle;
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType=="DATA/MC"){
        p->fRatioType = "DATA/BKG";
    }
    p->fLabelX = fLabelX;
    p->fLabelY = fLabelY;
    p->fLegendX1 = fLegendX1;
    p->fLegendX2 = fLegendX2;
    p->fLegendY = fLegendY;
    if(fLogScale) opt += " log";
    p->Draw(opt);
    return p;
}

//__________________________________________________________________________________
//
double Region::GetMultFactors( FitResults *fitRes, std::ofstream& pullTex,
                                const int i /*sample*/, const int i_bin /*bin number*/,
                                const double binContent0,
                                const std::string &var_syst_name,
                                const bool isUp ) const{
    double multNorm = 1.;
    double multShape = 0.;
    float systValue = 0.;
    SampleHist *sh = fSampleHists[i];
    for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
        SystematicHist *syh = sh->fSyst[i_syst];
        std::string systName = syh->fName;
        TString systNameNew(systName); // used in pull tables
        Systematic *syst = syh->fSystematic;
        if(syst!=nullptr){
            systName = syst->fNuisanceParameter;
            if(syst->fType==Systematic::SHAPE){
                continue;
            }
            else {
                if(var_syst_name!="" && systName==var_syst_name){
                    if(isUp) systValue = fitRes->GetNuisParValue(systName) + fitRes->GetNuisParErrUp(systName);
                    else     systValue = fitRes->GetNuisParValue(systName) + fitRes->GetNuisParErrDown(systName);
                }
                else {
                    systValue = fitRes->GetNuisParValue(systName);
                }
            }
        }
        else{
            systValue = fitRes->GetNuisParValue(systName);
        }
        
        //
        // Normalisation component: use the exponential interpolation and the multiplicative combination
        //
        if(syh->fIsOverall){
            double binContentUp   = (syh->fNormUp+1) * binContent0;
            double binContentDown = (syh->fNormDown+1) * binContent0;
            double factor = GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_overall);
            multNorm *= factor;
            if (fSampleHists[i]->fSample->fBuildPullTable>0){
                if (((factor > 1.01) || (factor < 0.99)) && (i_bin==1)) {
                    WriteDebugStatus("Region::DrawPostFit", "Syst " + systName +" in bin " + std::to_string(i_bin) + " has norm effect " + std::to_string( factor*100 ));
                    pullTex << setprecision(2) << "norm " << systNameNew.ReplaceAll("_","-") << "&" << (factor-1)*100 << " \\% \\\\\n" << std::endl;
                }
            }
        }

        //
        // Shape component: use the linear interpolation and the additive combination
        //
        if(syh->fIsShape){
            double binContentUp   = syh->fHistShapeUp->GetBinContent(i_bin);
            double binContentDown = syh->fHistShapeDown->GetBinContent(i_bin);
            double factor = GetDeltaN(systValue, binContent0, binContentUp, binContentDown, fIntCode_shape);
            multShape += factor - 1;
            if (fSampleHists[i]->fSample->fBuildPullTable==2){
                if (((factor-1) > 0.03) || ((factor-1) < - 0.03)) {
                    WriteDebugStatus("Region::DrawPostFit", "Syst " + systName +" in bin " + std::to_string(i_bin) + " has shape effect " + std::to_string( (factor-1)*100 ));
                    pullTex << setprecision(2) << "shape " << systNameNew.ReplaceAll("_","-") << " bin " << i_bin << "&" << (factor-1)*100  << " \\% \\\\\n" << std::endl;
                }
            }
        }
    }
    return multNorm*(multShape+1.);
}

//__________________________________________________________________________________
//
void Region::BuildPostFitErrorHist(FitResults *fitRes, const std::vector<std::string>& morph_names){

    WriteInfoStatus("Region::BuildPostFitErrorHist", "Building post-fit plot for region " + fName + " ...");
    double diffUp(0.), diffDown(0.);

    //
    // 0) Collect all the systematics on all the samples
    //
    fSystNames.clear();
    std::map<string,bool> systIsThere;
    systIsThere.clear();
    double systValue(0.);
    double systErrUp(0.);
    double systErrDown(0.);
    string systName = "";
    string systNuisPar = "";
    SystematicHist *sh = nullptr;
    string systNameSF = "";
    int iBinSF = 0;

    for(int i_sample=0;i_sample<fNSamples;i_sample++){
        if(fSampleHists[i_sample]->fSample->fType == Sample::DATA) continue;
        if(fSampleHists[i_sample]->fSample->fType == Sample::GHOST) continue;

        //
        // Norm factors
        //
        for(int i_norm=0;i_norm<fSampleHists[i_sample]->fNNorm;i_norm++){
            NormFactor *nf = fSampleHists[i_sample]->fNormFactors[i_norm];
            systName = nf->fName;
            // if this norm factor is a morphing one => save the nuis.par
            // skip POI if B-only fit FIXME
            if(fFitType==TRExFit::BONLY && systName==fPOI) continue;
            if(nf->fConst) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }

        //
        // Shape factors
        //
        // extract number of bins
        // loop over shape factors
        for(int i_shape=0;i_shape<fSampleHists[i_sample]->fNShape;i_shape++){
            systName = fSampleHists[i_sample]->fShapeFactors[i_shape]->fName;
            // add syst name for each bin
            for(int i_bin = 0; i_bin < fSampleHists[i_sample]->fHist->GetNbinsX(); i_bin++){
                systNameSF = systName + "_bin_" + std::to_string(i_bin);
                // the shape factor naming used i_bin - 1 for the first bin
                // add it as one syst per bin
                if(!systIsThere[systNameSF]){
                    fSystNames.push_back(systNameSF);
                    systIsThere[systNameSF] = true;
                }
            }
        }

        //
        // Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i_sample]->fNSyst;i_syst++){
            if(!fSampleHists[i_sample]->fSyst[i_syst]->fSystematic) continue;
            systName = fSampleHists[i_sample]->fSyst[i_syst]->fName;
            if(fSampleHists[i_sample]->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE) continue;
            if(!systIsThere[systName]){
                fSystNames.push_back(systName);
                systIsThere[systName] = true;
            }
        }

        //
        // SHAPE Systematics
        //
        for(int i_syst=0;i_syst<fSampleHists[i_sample]->fNSyst;i_syst++){
            if(!fSampleHists[i_sample]->fSyst[i_syst]->fSystematic) continue;
            systName = fSampleHists[i_sample]->fSyst[i_syst]->fName;
            if(systName.find("stat_")!=std::string::npos) continue; // fSeparateGammas already added later
            if(fSampleHists[i_sample]->fSyst[i_syst]->fSystematic->fType!=Systematic::SHAPE) continue;
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                std::string gammaName = Form("shape_%s_%s_bin_%d",systName.c_str(),fName.c_str(),i_bin-1);
                if(!systIsThere[gammaName]){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                }
            }
        }

        //
        // Gammas
        //
        if(fUseGammaPulls && (fSampleHists[i_sample]->fSample->fUseMCStat || fSampleHists[i_sample]->fSample->fSeparateGammas)){
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                std::string gammaName = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                if(fSampleHists[i_sample]->fSample->fSeparateGammas)
                    gammaName = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i_sample]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                if(!systIsThere[gammaName]){
                    fSystNames.push_back(gammaName);
                    systIsThere[gammaName] = true;
                }
            }
        }
    }
    //
    // 1) Build post-fit error hists (for each sample and each systematic):
    //
    std::vector<double> morph_scale;
    std::vector<double> morph_scale_nominal;
    PrepareMorphScales(fitRes, &morph_scale, &morph_scale_nominal);

    // - loop on systematics
    for(size_t i_syst=0;i_syst<fSystNames.size();++i_syst){
        WriteVerboseStatus("Region::BuildPostFitErrorHist", "    Systematic: " + fSystNames[i_syst]);

        int i_morph_sample = 0;
        //
        // Get fit result
        //
        systName    = fSystNames[i_syst];
        if(TRExFitter::NPMAP[systName]=="") TRExFitter::NPMAP[systName] = systName;
        systValue   = fitRes->GetNuisParValue(TRExFitter::NPMAP[systName]);
        systErrUp   = fitRes->GetNuisParErrUp(TRExFitter::NPMAP[systName]);
        systErrDown = fitRes->GetNuisParErrDown(TRExFitter::NPMAP[systName]);
        
        WriteVerboseStatus("Region::BuildPostFitErrorHist", "      alpha = " + std::to_string(systValue) + " +" + std::to_string(systErrUp) + " " + std::to_string(systErrDown));

        // needed for morph samples
        const int nbins = fTot_postFit->GetNbinsX();

        std::vector<double> morph_nominal(nbins, 0.0);
        std::vector<double> morph_nominal_postfit(nbins, 0.0);
        std::vector<double> morph_up_postfit(nbins, 0.0);
        std::vector<double> morph_down_postfit(nbins, 0.0);

        std::vector<double> morph_syst_up(nbins, 0.0);
        std::vector<double> morph_syst_down(nbins, 0.0);

        bool isMorph = false;
        int morph_index = -1;

        // - loop on samples
        for(int i=0;i<fNSamples;i++){
            // skip data
            if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
            if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
            WriteVerboseStatus("Region::BuildPostFitErrorHist", "  Sample: " + fSampleHists[i]->fName);

            //
            // Get SystematicHist
            //
            sh = fSampleHists[i]->GetSystematic(systName);

            // hack: add a systematic hist if not there
            if(sh==nullptr){
                fSampleHists[i]->AddHistoSyst(systName,fSampleHists[i]->fHist,fSampleHists[i]->fHist);
                sh = fSampleHists[i]->GetSystematic(systName);
            }
            
            //
            // initialize the up and down variation histograms
            // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
            //
            sh->fHistUp_postFit   = (TH1*)fSampleHists[i]->fHist_postFit->Clone(Form("%s_%s_Up_postFit",  fSampleHists[i]->fHist->GetName(),systName.c_str()));
            sh->fHistDown_postFit = (TH1*)fSampleHists[i]->fHist_postFit->Clone(Form("%s_%s_Down_postFit",fSampleHists[i]->fHist->GetName(),systName.c_str()));

            // check if the sample is morph sample
            for (const auto& i_morph : morph_names){
                if (fSampleHists[i]->fIsMorph[i_morph]){
                    isMorph = true;
                    morph_index = i;
                    break;
                }
            }

            // - loop on bins
            for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                diffUp = 0.;
                diffDown = 0.;
                double yieldNominal = fSampleHists[i]->fHist->GetBinContent(i_bin);  // store nominal yield for this bin
                double yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);  // store nominal yield for this bin, but do it post fit

                if (isMorph){
                    float scaleNom  = morph_scale_nominal.at(i_morph_sample);
                    float scaleNom_postfit  = morph_scale.at(i_morph_sample);
                    morph_nominal_postfit.at(i_bin-1)+= yieldNominal*scaleNom_postfit;
                    morph_nominal.at(i_bin-1)+= yieldNominal*scaleNom;
                }

                size_t posTmp = systName.find("_bin_");
                std::string gammaName      = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                std::string gammaNameShape = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                //
                // if it's a gamma
                if(gammaName==fSystNames[i_syst] && fSampleHists[i]->fSample->fUseMCStat && !fSampleHists[i]->fSample->fSeparateGammas){
                    diffUp   += yieldNominal_postFit*systErrUp;
                    diffDown += yieldNominal_postFit*systErrDown;
                    if (isMorph){
                        morph_up_postfit.at(i_bin-1)  += yieldNominal_postFit*(systErrUp+1);
                        morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                    }
                }
                //
                // if it's a specific-sample gamma
                else if(gammaNameShape==fSystNames[i_syst] && fSampleHists[i]->fSample->fSeparateGammas){
                    diffUp   += yieldNominal_postFit*systErrUp;
                    diffDown += yieldNominal_postFit*systErrDown;
                    if (isMorph){
                        morph_up_postfit.at(i_bin-1)  += yieldNominal_postFit*(systErrUp+1);
                        morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                    }
                }
                //
                // if it's a shape-systematic gamma
                else if(fSystNames[i_syst].find("shape_")!=std::string::npos && fSystNames[i_syst].find(Form("%s_bin_%d",fName.c_str(),i_bin-1))!=std::string::npos){
                    for(auto syh : fSampleHists[i]->fSyst){
                        Systematic *syst = syh->fSystematic;
                        if(!syst) continue;
                        if(syst->fType==Systematic::SHAPE){
                            std::string gammaNameShapeSyst = Form("shape_%s_%s_bin_%d",syst->fName.c_str(),fName.c_str(),i_bin-1);
                            if(gammaNameShapeSyst==fSystNames[i_syst]){
                                diffUp   += yieldNominal_postFit*systErrUp;
                                diffDown += yieldNominal_postFit*systErrDown;
                                if (isMorph){
                                    morph_up_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrUp+1);
                                    morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                                }
                            }
                        }
                    }
                }
                //
                // if it's a norm factor
                else if(fSampleHists[i]->HasNorm(fSystNames[i_syst])){
                    // if this norm factor is a morphing one
                    if(fSystNames[i_syst].find("morph_")!=string::npos || fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fExpression.first!=""){
                        std::string formula = TRExFitter::SYSTMAP[fSystNames[i_syst]];
                        std::string name = TRExFitter::NPMAP[fSystNames[i_syst]];
                        WriteDebugStatus("Region::BuildPostFitErrorHist", "formula: " +formula);
                        WriteDebugStatus("Region::BuildPostFitErrorHist", "name: " +name);
                        std::vector < std::pair < std::string,std::vector<double> > > nameS;
                        if(fSystNames[i_syst].find("morph_")!=std::string::npos){
                            nameS.push_back(std::make_pair(name,std::vector<double>{double(fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fNominal),
                            double(fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fMin),double(fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fMax)}));
                        }
                        else {
                            nameS = processString(name);
                        }
                        std::vector <double> nfValuevec, nfUpvec, nfDownvec;
                        for (unsigned int j = 0; j<nameS.size(); j++){
                            formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                            nfValuevec.push_back(fitRes->GetNuisParValue(nameS[j].first));
                            nfUpvec.push_back(fitRes->GetNuisParValue(nameS[j].first) + fitRes->GetNuisParErrUp(nameS[j].first));
                            nfDownvec.push_back(fitRes->GetNuisParValue(nameS[j].first) + fitRes->GetNuisParErrDown(nameS[j].first));
                        }
                        double *nfValue = nfValuevec.data();
                        WriteDebugStatus("Region::BuildPostFitErrorHist", "formula: " +formula);
                        for(unsigned int j = 0; j<nameS.size(); j++){
                            WriteDebugStatus("Region::BuildPostFitErrorHist", "nfValue["+std::to_string(j)+"]: "+std::to_string(nfValue[j]));
                        }
                        TFormula f_morph ("f_morph",formula.c_str());
                        float scaleUp = f_morph.EvalPar(nfValue,nullptr); // nominal value
                        float scaleDown = f_morph.EvalPar(nfValue,nullptr); // nominal value
                        if(fSystNames[i_syst].find("morph_")!=std::string::npos){
                            double *nfUpValue = nfUpvec.data();
                            double *nfDownValue = nfDownvec.data();
                            scaleUp = f_morph.EvalPar(nfUpValue,nullptr); // 1-parameter up variation
                            scaleDown = f_morph.EvalPar(nfDownValue,nullptr); // 1-parameter down variation
                        }
                        // multi-parameter dependence => find the combinatio with largest/smallest "express" (definition of exprUp exprDown)
                        else {
                            for (int ii = 0; ii < (1 << nameS.size()); ii++) {
                                std::vector <double> exprvec;
                                for(unsigned int j=0;j<nameS.size();j++){
                                    if(ii & (1<<j)) exprvec.push_back(nfUpvec[j]);
                                    else            exprvec.push_back(nfDownvec[j]);
                                }
                                double *expr = exprvec.data();
                                scaleUp = (f_morph.EvalPar(expr,nullptr) > scaleUp) ? f_morph.EvalPar(expr,nullptr) : scaleUp;
                                scaleDown = (f_morph.EvalPar(expr,nullptr) < scaleDown) ? f_morph.EvalPar(expr,nullptr) : scaleDown;
                            }
                        }
                        morph_syst_up.at(i_bin-1)   += yieldNominal*scaleUp;
                        morph_syst_down.at(i_bin-1) += yieldNominal*scaleDown;
                    }
                    else{
                        diffUp   += yieldNominal_postFit*systErrUp/systValue;
                        diffDown += yieldNominal_postFit*systErrDown/systValue;
                        if (isMorph){
                            morph_up_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrUp+1);
                            morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                        }
                    }
                }
                //
                // ShapeFactor have to get NP per bin
                else if(posTmp != std::string::npos){
                    // get the shape factor name without bin index
                    systNameSF = systName.substr(0, posTmp);
                    // get the shape factor bin as integer
                    iBinSF = std::atoi(systName.substr(posTmp + 5).c_str()) + 1;
                    // FIXME could still be a problem with pruning?
                    if(fSampleHists[i]->HasShapeFactor(systNameSF) && iBinSF == i_bin){
                        diffUp   += yieldNominal_postFit*systErrUp;
                        diffDown += yieldNominal_postFit*systErrDown;
                        if (isMorph){
                            morph_up_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrUp+1);
                            morph_down_postfit.at(i_bin-1)+= yieldNominal_postFit*(systErrDown+1);
                        }
                    }
                }
                //
                // Systematics treatment
                //
                else if(fSampleHists[i]->HasSyst(fSystNames[i_syst])){
                    std::ofstream dummy;
                    double multNom  = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal );
                    double multUp   = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal, TRExFitter::NPMAP[systName], true);
                    double multDown = GetMultFactors( fitRes, dummy, i, i_bin, yieldNominal, TRExFitter::NPMAP[systName], false);
                    if (isMorph){
                        morph_up_postfit.at(i_bin-1)+=(multUp/multNom)*yieldNominal_postFit;
                        morph_down_postfit.at(i_bin-1)+= (multDown/multNom)*yieldNominal_postFit;
                    } else {
                        diffUp   += (multUp/multNom   - 1.)*yieldNominal_postFit;
                        diffDown += (multDown/multNom - 1.)*yieldNominal_postFit;
                    }
                }

                WriteVerboseStatus("Region::BuildPostFitErrorHist", "        Bin " + std::to_string(i_bin) + ":   " + "\t +" + std::to_string(100*diffUp/yieldNominal)
                    + "%\t " + std::to_string(100*diffDown/yieldNominal) + "%");

                //
                // Add the proper bin content to the variation hists (coming from post-fit total histogram)
                //
                if (!isMorph){
                    sh->fHistUp_postFit  ->AddBinContent( i_bin, diffUp   );
                    sh->fHistDown_postFit->AddBinContent( i_bin, diffDown );
                }

            } // loop over bins
        } // loop over samples
        if(isMorph) i_morph_sample++;

        if (isMorph){
            // now apply the corrections from morph
            // Use it only ONCE for all combines morph templates
            sh = fSampleHists[morph_index]->GetSystematic(systName);

            // hack: add a systematic hist if not there
            if(sh==nullptr){
                fSampleHists[morph_index]->AddHistoSyst(systName,fSampleHists[morph_index]->fHist,fSampleHists[morph_index]->fHist);
                sh = fSampleHists[morph_index]->GetSystematic(systName);
            }

            //
            // initialize the up and down variation histograms
            // (note: do it even if the syst is not there; in this case the variation hist will be = to the nominal)
            //
            sh->fHistUp_postFit   = (TH1*)fSampleHists[morph_index]->fHist_postFit->Clone(Form("%s_%s_Up_postFit",  fSampleHists[morph_index]->fHist->GetName(),systName.c_str()));
            sh->fHistDown_postFit = (TH1*)fSampleHists[morph_index]->fHist_postFit->Clone(Form("%s_%s_Down_postFit",fSampleHists[morph_index]->fHist->GetName(),systName.c_str()));

            // loop over bins
            for (int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                sh->fHistUp_postFit    ->AddBinContent( i_bin, (morph_up_postfit.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
                sh->fHistDown_postFit  ->AddBinContent( i_bin, (morph_down_postfit.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
            }

            // uncertainty from the morph uncertainty
            if (fSystNames[i_syst].find("morph_")!=string::npos){
                for (int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
                    sh->fHistUp_postFit    ->AddBinContent( i_bin, (morph_syst_up.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
                    sh->fHistDown_postFit  ->AddBinContent( i_bin, (morph_syst_down.at(i_bin-1)-morph_nominal_postfit.at(i_bin-1)) );
                }
            }
        }
    }// loop over systs

    // at this point all the sample-by-sample post-fit variation histograms should be filled
    //
    // Now build the total prediction variations, for each systematic
    //

    // - loop on systematics
    for(size_t i_syst=0;i_syst<fSystNames.size();++i_syst){
        systName = fSystNames[i_syst];
        //
        // Initialize the tot variation hists
        //
        fTotUp_postFit[i_syst]   = (TH1*)fTot_postFit->Clone(Form("h_tot_%s_Up_postFit",  systName.c_str()));
        fTotUp_postFit[i_syst] -> Scale(0); //initialising the content to 0
        fTotDown_postFit[i_syst] = (TH1*)fTot_postFit->Clone(Form("h_tot_%s_Down_postFit",systName.c_str()));
        fTotDown_postFit[i_syst] -> Scale(0); //initialising the content to 0

        // - loop on bins
        for(int i_bin=1;i_bin<fTot_postFit->GetNbinsX()+1;i_bin++){
            diffUp = 0.;
            diffDown = 0.;
            // - loop on samples
            for(int i=0;i<fNSamples;i++){
                // skip data
                if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
                if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
                if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;
                // skip signal if Bkg only
                if(fFitType==TRExFit::BONLY && fSampleHists[i]->fSample->fType==Sample::SIGNAL) continue;
                // get SystematicHist
                sh = fSampleHists[i]->GetSystematic(systName);
                // increase diffUp/Down according to the previously stored histograms
                // yieldNominal_postFit = fSampleHists[i]->fHist_postFit->GetBinContent(i_bin);
                if(!sh) continue;
                diffUp   += sh->fHistUp_postFit  ->GetBinContent(i_bin);
                diffDown += sh->fHistDown_postFit->GetBinContent(i_bin);
            }
            // add the proper bin content to the variation hists
            fTotUp_postFit[i_syst]  ->AddBinContent( i_bin, diffUp   );
            fTotDown_postFit[i_syst]->AddBinContent( i_bin, diffDown );
        }
    }

    // at this point all the total expectation post-fit variation histograms should be filled
    //
    // Build the vectors of variations (necessary to call BuildTotError)
    //
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    std::vector<std::string> systNuisPars;
    for(size_t i_syst=0;i_syst<fSystNames.size();++i_syst){
        h_up.  push_back( fTotUp_postFit[i_syst]   );
        h_down.push_back( fTotDown_postFit[i_syst] );
        systNuisPars.push_back(TRExFitter::NPMAP[fSystNames[i_syst]]);
    }
    fErr_postFit = BuildTotError( fTot_postFit, h_up, h_down, systNuisPars, fitRes->fCorrMatrix );
    fErr_postFit->SetName("g_totErr_postFit");
    // at this point fTot and fErr _postFit should be ready

    //
    // Goodness of post-fit
    if(fGetChi2){
        if (!fHasData){
            WriteWarningStatus("Region::BuildPostFitErrorHist", "Data histogram is nullptr, cannot calculate Chi2 agreement.");
            WriteWarningStatus("Region::BuildPostFitErrorHist", "Maybe you do not have data sample defined?");
            return;
        }
        // remove blinded bins
        auto h_data = std::unique_ptr<TH1>(static_cast<TH1*>(fData->fHist->Clone()));
        if(fBlindedBins!=nullptr){
            for(int i_bin=1;i_bin<=h_data->GetNbinsX();++i_bin){
                if(fBlindedBins->GetBinContent(i_bin)>0) h_data->SetBinContent(i_bin,-1);
                if(find(fDropBins.begin(),fDropBins.end(),i_bin-1)!=fDropBins.end()) h_data->SetBinContent(i_bin,-1);
            }
        }
        if(fGetChi2==1) fSystNames.clear();
        std::pair<double,int> res = GetChi2Test( h_data.get(), fTot_postFit, h_up, fSystNames, fitRes->fCorrMatrix );
        fChi2val = res.first;
        fNDF = res.second;
        fChi2prob = ROOT::Math::chisquared_cdf_c( res.first, res.second);
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- POST-FIT AGREEMENT EVALUATION -----------------------");
        if(fGetChi2==1)
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- -------- STAT-ONLY --------- -----------------------");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "--- REGION " + fName + ":");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "  chi2        = " + std::to_string(fChi2val));
        WriteInfoStatus("Region::BuildPostFitErrorHist", "  ndof        = " + std::to_string(fNDF));
        WriteInfoStatus("Region::BuildPostFitErrorHist", "  probability = " + std::to_string(fChi2prob));
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- ---------------------------- -----------------------");
        WriteInfoStatus("Region::BuildPostFitErrorHist", "----------------------- ---------------------------- -----------------------");
    }

}

//__________________________________________________________________________________
//
TRExPlot* Region::DrawPostFit(FitResults *fitRes,ofstream& pullTex, const std::vector<std::string> &morph_names, const std::vector<int>& canvasSize, string opt){

    if(TRExFitter::PREFITONPOSTFIT){
        fPlotPostFit->h_tot_bkg_prefit = (TH1*)fPlotPreFit->GetTotBkg()->Clone("h_tot_bkg_prefit");
    }

    TRExPlot *p{};
    if (canvasSize.size() == 0){
        p = fPlotPostFit;
        p->fShowYields = TRExFitter::SHOWYIELDS;
    } else {
        p = new TRExPlot(("c_"+fName).c_str(), canvasSize.at(0), canvasSize.at(1),TRExFitter::NORATIO);
    }

    p->SetXaxisRange(fXaxisRange);
    if(fYmaxScale==0) p->SetYmaxScale(1.8);
    else              p->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) p->fYmax = fYmax;
    if(fYmin!=0) p->fYmin = fYmin;
    p->fRatioYmax = fRatioYmaxPostFit;
    p->fRatioYmin = fRatioYminPostFit;
    p->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") p->SetYaxis(fYTitle);
    //
    // For 4-top-style plots
    if(TRExFitter::OPTION["FourTopStyle"]>0){
        if(fRegionType==CONTROL) p->AddLabel("#font[62]{Control Region}");
        else if(fRegionType==SIGNAL) p->AddLabel("#font[62]{Signal Region}");
        else if(fRegionType==VALIDATION) p->AddLabel("#font[62]{Validation Region}");
        else p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        p->AddLabel("#font[52]{Post-fit}");
    }
    //
    // old-style plots
    else{
        p->AddLabel(fFitLabel);
        p->AddLabel(fLabel);
        if(TRExFitter::OPTION["NoPrePostFit"]==0) p->AddLabel("Post-Fit");
    }
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    p->fLegendNColumns = fLegendNColumns;

    if(fBinLabels.size() && ((int)fBinLabels.size()==fNbins)) {
        for(int i_bin=0; i_bin<fNbins; i_bin++) {
            p->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
        }
    }

    //
    // 0) Create a new hist for each sample
    //
    TH1* hSmpNew[MAXsamples];
    for(int i=0;i<fNSamples;i++){
        hSmpNew[i] = (TH1*)fSampleHists[i]->fHist->Clone();
        // set to 0 uncertainty in each bin if MCstat set to FALSE
        if((!fSampleHists[i]->fSample->fUseMCStat && !fSampleHists[i]->fSample->fSeparateGammas) || fUseGammaPulls){
            for(int i_bin=0;i_bin<hSmpNew[i]->GetNbinsX()+2;i_bin++) hSmpNew[i]->SetBinError(i_bin,0.);
        }
    }

    //
    // 1) Propagates the post-fit NP values to the central value (pulls)
    //
    string systName;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        //
        if (fSampleHists[i]->fSample->fBuildPullTable>0){
            WriteDebugStatus("Region::DrawPostFit", "Propagating post-fit to Sample " + fSampleHists[i]->fSample->fTitle);
            TString sampleTex= fSampleHists[i]->fSample->fTexTitle;
            pullTex << "\\hline\n" << endl;
            pullTex << "{\\color{blue}{$\\rightarrow \\,$ "<< sampleTex << "}} & \\\\\n"<< endl;
        }
        TH1* hNew = nullptr;
        hNew = (TH1*)hSmpNew[i]->Clone();
        for(int i_bin=1;i_bin<=hNew->GetNbinsX();i_bin++){
            double binContent0 = hSmpNew[i]->GetBinContent(i_bin);
            double mult_factor = GetMultFactors(fitRes, pullTex, i, i_bin, binContent0);

            //
            // Final computation
            //
            double binContentNew = binContent0*mult_factor;

            //
            // stat gammas
            if(fUseGammaPulls && (fSampleHists[i]->fSample->fUseMCStat || fSampleHists[i]->fSample->fSeparateGammas)){
                // find the gamma for this bin of this distribution in the fit results
                std::string gammaName = Form("stat_%s_bin_%d",fName.c_str(),i_bin-1);
                if(fSampleHists[i]->fSample->fSeparateGammas)
                    gammaName = Form("shape_stat_%s_%s_bin_%d",fSampleHists[i]->fSample->fName.c_str(),fName.c_str(),i_bin-1);
                WriteDebugStatus("Region::DrawPostFit", "Looking for gamma " + gammaName);
                float gammaValue = fitRes->GetNuisParValue(gammaName);
                WriteDebugStatus("Region::DrawPostFit", "  -->  pull = " + std::to_string(gammaValue));
                // linear effect
                if(gammaValue>0) binContentNew *= gammaValue;
            }
            // gammas from SHAPE systematics
            for(auto syh : fSampleHists[i]->fSyst){
                Systematic *syst = syh->fSystematic;
                if(!syst) continue;
                if(syst->fType==Systematic::SHAPE){
                    std::string gammaName = Form("shape_%s_%s_bin_%d",syst->fName.c_str(),fName.c_str(),i_bin-1);
                    WriteDebugStatus("Region::DrawPostFit", "Looking for gamma " + gammaName);
                    float gammaValue = fitRes->GetNuisParValue(gammaName);
                    WriteDebugStatus("Region::DrawPostFit", "  -->  pull = " + std::to_string(gammaValue));
                    // linear effect
                    if(gammaValue>0) binContentNew *= gammaValue;
                }
            }

            //
            // Setting to new values
            //
            hNew->SetBinContent(i_bin,binContentNew);
        }
        hSmpNew[i] = (TH1*)hNew->Clone();
        fSampleHists[i]->fHist_postFit = hSmpNew[i];
        delete hNew;
    }
    //
    // 2) Scale all samples by norm factors
    //    Done after the propagation of the NP (avoids nans due to "0" value of some NormFactors)
    //    Seems consistent with application in Roostats
    //
    string nfName;
    float nfValue;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        for(int i_norm=0;i_norm<fSampleHists[i]->fNNorm;i_norm++){
            NormFactor *nf = fSampleHists[i]->fNormFactors[i_norm];
            nfName = nf->fName;
            if(nf->fConst) nfValue = nf->fNominal;
            else           nfValue = fitRes->GetNuisParValue(TRExFitter::NPMAP[nfName]);
            //
            // if this norm factor is a morphing one
            if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                std::string formula = TRExFitter::SYSTMAP[nfName];
                std::string name = TRExFitter::NPMAP[nfName];
                WriteDebugStatus("Region::DrawPostFit", "formula: " +formula);
                WriteDebugStatus("Region::DrawPostFit", "name: " +name);
                std::vector < std::pair < std::string,std::vector<double> > > nameS;
                if(nf->fName.find("morph_")!=std::string::npos){
                    nameS.push_back(std::make_pair(name,std::vector<double>{double(nf->fNominal),double(nf->fMin),double(nf->fMax)}));
                }
                else {
                    nameS = processString(name);
                }
                std::vector <double> nfValuevec;
                for (unsigned int j = 0; j<nameS.size(); j++){
                    formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                    nfValuevec.push_back(fitRes->GetNuisParValue(nameS[j].first));
                }
                double *nfValuemorph = nfValuevec.data();
                WriteDebugStatus("Region::DrawPostFit", "formula: " +formula);
                for(unsigned int j = 0; j<nameS.size(); j++){
                    WriteDebugStatus("Region::DrawPostFit", "nfValuemorph["+std::to_string(j)+"]: "+std::to_string(nfValuemorph[j]));
                }
                TFormula f_morph ("f_morph",formula.c_str());
                float scale = 1;
                scale = f_morph.EvalPar(nfValuemorph,nullptr);
                hSmpNew[i]->Scale(scale);
            }
            else{
                hSmpNew[i]->Scale(nfValue);
            }
        }
    }


    //
    // 2)b) Scale all samples by shape factors factors
    //

    string sfName;
    string sfNameBin;
    float sfValue;
    double binContentSFNew = 0;
    int iBinSF = 0;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        for(int i_shape=0;i_shape<fSampleHists[i]->fNShape;i_shape++){
            ShapeFactor *sf = fSampleHists[i]->fShapeFactors[i_shape];
            sfName = sf->fName;
            if(sfName.find("saturated_model_sf")!=std::string::npos) continue;
            // loop over bins
            // there should be a NP per bin in the fit file
            // already checked by GetNuisParValue()
            // the shape factor naming used i_bin - 1 for the first bin
            for(int i_bin = 1; i_bin <= hSmpNew[i]->GetNbinsX(); i_bin++){
                iBinSF = i_bin - 1;
                sfNameBin = sfName + "_bin_" + std::to_string(iBinSF);
                if(sf->fConst) sfValue = sf->fNominal;
                else           sfValue = fitRes->GetNuisParValue(sfNameBin);
                // scale bin content by shape factor
                binContentSFNew = hSmpNew[i]->GetBinContent(i_bin) * sfValue;
                // Setting to new value
                hSmpNew[i]->SetBinContent(i_bin, binContentSFNew);
            }
        }
    }

    // Michele:
    // Scale samples acording to fScaleSamplesToData:
    if(fScaleSamplesToData.size()>0){
        TH1* hTot = nullptr;
        for(int i=0;i<fNSamples;i++){
            if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
            if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
            if(hTot==nullptr) hTot = (TH1*)hSmpNew[i]->Clone("hTotPostFit");
            else              hTot->Add(hSmpNew[i]);
        }
        float totPred = hTot->Integral();
        float totData = fData->fHist->Integral();
        float totToScale = 0;
        std::vector<int> shIdxToScale;
        for(int i=0;i<fNSamples;i++){
            SampleHist *sh = fSampleHists[i];
            if(sh->fHist==nullptr) continue;
            if(sh->fSample->fType==Sample::GHOST){
                WriteWarningStatus("Region::DrawPostFit","Requested to scale to data a GHOST sample, " + sh->fSample->fName + ". Skipping this sample.");
                continue;
            }
            if(FindInStringVector(fScaleSamplesToData,sh->fSample->fName)>=0){
                shIdxToScale.emplace_back(i);
                totToScale += hSmpNew[i]->Integral();
            }
        }
        if(totToScale>0 && shIdxToScale.size()>0){
            float scale = (totData-(totPred-totToScale))/totToScale;
            for(auto idx : shIdxToScale){
                WriteInfoStatus("Region::DrawPostFit","Scaling sample " + fSampleHists[idx]->fSample->fName + " by " + std::to_string(scale) + " in region " + fName);
                hSmpNew[idx]->Scale(scale);
            }
        }
    }
    
    //
    // 3) Add the new Sig and Bkg to plot
    //
    string title;
    TH1* hBkgNew[MAXsamples];
    TH1* hSigNew[MAXsamples];
    for(int i=0, i_bkg=0, i_sig=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::BACKGROUND){
            hBkgNew[i_bkg] = hSmpNew[i];
            i_bkg++;
        }
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL){
            hSigNew[i_sig] = hSmpNew[i];
            i_sig++;
        }
    }
    if(fHasData && opt.find("blind")==string::npos) p->SetData(fData->fHist,fData->fSample->fTitle);
    for(int i=0;i<fNSig;i++){
        title = fSig[i]->fSample->fTitle;
        if(fSig[i]->fSample->fGroup != "") title = fSig[i]->fSample->fGroup;
        if(TRExFitter::SHOWSTACKSIG)    p->AddSignal(    hSigNew[i],title);
        if(TRExFitter::SHOWNORMSIG){
            if( (TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
             || !TRExFitter::OPTION["NormSigSRonly"] )
                p->AddNormSignal(hSigNew[i],title);
        }
        else{
            if(TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) p->AddNormSignal(hSigNew[i],title);
        }
        if(TRExFitter::SHOWOVERLAYSIG){
            ScaleNominal(fSig[i],hSigNew[i]);
            p->AddOverSignal(hSigNew[i],title);
        }
    }
    for(int i=0;i<fNBkg;i++){
        title = fBkg[i]->fSample->fTitle;
        if(fBkg[i]->fSample->fGroup != "") title = fBkg[i]->fSample->fGroup;
        p->AddBackground(hBkgNew[i],title);
    }

    //
    // 4) Build post-fit error band
    //    Build total histogram (fTot_postFit)
    //
    int j = 0;
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        if(fSampleHists[i]->fSample->fType==Sample::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;
        if(j==0) fTot_postFit = (TH1*)hSmpNew[i]->Clone("h_tot_postFit");
        else fTot_postFit->Add(hSmpNew[i]);
        j++;
    }

    //
    // MCstat set to 0 if disabled
    //
    if(!fUseStatErr || fUseGammaPulls){
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
            fTot_postFit->SetBinError(i_bin,0);
        }
    }

    p->SetTotBkg(fTot_postFit);
    if(fBinWidth>0) p->SetBinWidth(fBinWidth);

    //
    // blinding bins
    //
    if(fBlindingThreshold>=0){
        p->SetBinBlinding(true,fBlindingThreshold,fBlindingType);
        if(fKeepPrefitBlindedBins && fBlindedBins!=nullptr) p->SetBinBlinding(true,fBlindedBins,fBlindingType);
    }
    p->BlindData();

    //
    // Build error band
    //
    BuildPostFitErrorHist(fitRes, morph_names);

    //
    // Print chi2 info
    //
    if(fGetChi2 && TRExFitter::SHOWCHI2)  p->SetChi2KS(fChi2prob,-1,fChi2val,fNDF);

    //
    // 5) Finishes configuration of TRExPlot objects
    //
    p->SetTotBkgAsym(fErr_postFit);
    p->fATLASlabel = fATLASlabel;
    p->fRatioYtitle = fRatioYtitle;
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType=="DATA/MC"){
        p->fRatioType = "DATA/BKG";
    }
    p->fLabelX = fLabelX;
    p->fLabelY = fLabelY;
    p->fLegendX1 = fLegendX1;
    p->fLegendX2 = fLegendX2;
    p->fLegendY = fLegendY;
    if(fLogScale) opt += " log";

    p->Draw(opt);

    //
    // Print bin content and errors
    //
    WriteDebugStatus("Region::DrawPostFit", "--------------------");
    WriteDebugStatus("Region::DrawPostFit", "Final bin contents");
    WriteDebugStatus("Region::DrawPostFit", "--------------------");
        for(int i_bin=1;i_bin<=fTot_postFit->GetNbinsX();i_bin++){
        WriteDebugStatus("Region::DrawPostFit", std::to_string(i_bin) + ":\t" + std::to_string(fTot_postFit->GetBinContent(i_bin)) + " +" +
        std::to_string( fErr_postFit->GetErrorYhigh(i_bin-1)) + " -" + std::to_string(fErr_postFit->GetErrorYlow(i_bin-1)));
    }

    //
    // Save in a root file...
    //
    gSystem->mkdir((fFitName+"/Histograms").c_str());
    WriteInfoStatus("Region::DrawPostFit", "Writing file " + fFitName+"/Histograms/"+fName+fSuffix+"_postFit.root");
    TFile *f = new TFile((fFitName+"/Histograms/"+fName+fSuffix+"_postFit.root").c_str(),"RECREATE");
    f->cd();
    fErr_postFit->Write("",TObject::kOverwrite);
    fTot_postFit->Write("",TObject::kOverwrite);
    if( fPlotPostFit->h_tot_bkg_prefit)
      fPlotPostFit->h_tot_bkg_prefit->Write("",TObject::kOverwrite);
    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        if(fTotUp_postFit[i_syst])   fTotUp_postFit[i_syst]  ->Write("",TObject::kOverwrite);
        if(fTotDown_postFit[i_syst]) fTotDown_postFit[i_syst]->Write("",TObject::kOverwrite);
    }
    for(int i=0;i<fNSamples;i++){
        if(fSampleHists[i]->fSample->fType == Sample::DATA){
            fSampleHists[i]->fHist->Write(Form("h_%s",fSampleHists[i]->fName.c_str()),TObject::kOverwrite);
            continue;
        }
        if(fSampleHists[i]->fHist_postFit){
            fSampleHists[i]->fHist_postFit->Write(Form("h_%s_postFit",fSampleHists[i]->fName.c_str()),TObject::kOverwrite);
            for(unsigned int i_syst=0;i_syst<fSampleHists[i]->fSyst.size();i_syst++){
                if(fSampleHists[i]->fSyst[i_syst])
                    if(fSampleHists[i]->fSyst[i_syst]->fHistUp_postFit)
                        fSampleHists[i]->fSyst[i_syst]->fHistUp_postFit  ->Write(
                          Form("h_%s_%s_Up_postFit",fSampleHists[i]->fName.c_str(),fSampleHists[i]->fSyst[i_syst]->fName.c_str()), TObject::kOverwrite);
                if(fSampleHists[i]->fSyst[i_syst])
                    if(fSampleHists[i]->fSyst[i_syst]->fHistDown_postFit)
                        fSampleHists[i]->fSyst[i_syst]->fHistDown_postFit->Write(
                          Form("h_%s_%s_Down_postFit",fSampleHists[i]->fName.c_str(),fSampleHists[i]->fSyst[i_syst]->fName.c_str()),TObject::kOverwrite);
            }
        }
    }
    f->Close();
    delete f;
    //
    return p;
}

//__________________________________________________________________________________
//
void Region::AddSelection(const std::string& selection){
    if(selection=="") return;
    if(fSelection=="1" || fSelection=="") fSelection = selection;
    else fSelection += " && "+selection;
}

//__________________________________________________________________________________
//
void Region::AddMCweight(const std::string& weight){
    if(weight=="") return;
    if(fMCweight=="1" || fMCweight=="") fMCweight = weight;
    else fMCweight += " * "+weight;
}

//__________________________________________________________________________________
//
void Region::SetVariable(const std::string& variable,int nbin,float xmin,float xmax,string corrVar1,string corrVar2){
    fVariable = variable;
    fCorrVar1 = corrVar1;
    fCorrVar2 = corrVar2;
    fNbins = nbin;
    fXmin = xmin;
    fXmax = xmax;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeVariable(const std::string& variable, const std::string& sample){
    fAlternativeVariables[sample] = variable;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeSelection(const std::string& selection, const std::string& sample){
    fAlternativeSelections[sample] = selection;
}

//__________________________________________________________________________________
//
bool Region::UseAlternativeVariable(const std::string& sample){
    std::vector<std::string> tmpVec;
    for(auto tmp : fAlternativeVariables){
        tmpVec.push_back(tmp.first);
    }
    if (FindInStringVector(tmpVec,sample)<0){
        return false;
    }
    else {
        return true;
    }
}

//__________________________________________________________________________________
//
bool Region::UseAlternativeSelection(const std::string& sample){
    std::vector<std::string> tmpVec;
    for(auto tmp : fAlternativeSelections){
        tmpVec.push_back(tmp.first);
    }
    if (FindInStringVector(tmpVec,sample)<0){
        return false;
    }
    else {
        return true;
    }
}

//__________________________________________________________________________________
//
std::string Region::GetAlternativeVariable(const std::string& sample) const{
    std::vector<std::string> tmpVec;
    std::vector<std::string> tmpVec2;
    for(auto tmp : fAlternativeVariables){
        tmpVec.push_back(tmp.first);
        tmpVec2.push_back(tmp.second);
    }
    int idx = FindInStringVector(tmpVec,sample);
    if(idx<0){
        return "";
    }
    else {
        return tmpVec2[idx];
    }
}

//__________________________________________________________________________________
//
std::string Region::GetAlternativeSelection(const std::string& sample) const{
    std::vector<std::string> tmpVec;
    std::vector<std::string> tmpVec2;
    for(auto tmp : fAlternativeSelections){
        tmpVec.push_back(tmp.first);
        tmpVec2.push_back(tmp.second);
    }
    int idx = FindInStringVector(tmpVec,sample);
    if(idx<0){
        return "";
    }
    else {
        return tmpVec2[idx];
    }
}

//__________________________________________________________________________________
//
void Region::SetVariableTitle(const std::string& name){
    fVariableTitle = name;
}

//__________________________________________________________________________________
//
void Region::SetLabel(const std::string& label, std::string shortLabel){
    fLabel = label;
    if(shortLabel=="") fShortLabel = label;
    else fShortLabel = shortLabel;
}

//__________________________________________________________________________________
//
void Region::Print() const{
    WriteInfoStatus("Region::Print", "    Region: " + fName);
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        fSampleHists[i_smp]->Print();
    }
}

//__________________________________________________________________________________
//
void Region::PrintSystTable(FitResults *fitRes, string opt) const{
    bool isPostFit  = false; if(opt.find("post")!=string::npos)     isPostFit  = true;
    bool doClean    = false; if(opt.find("clean")!=string::npos)    doClean    = true;
    bool doCategory = false; if(opt.find("category")!=string::npos) doCategory = true;
    bool standalone = false; if(opt.find("standalone")!=string::npos)standalone= true;
    bool landscape  = false; if(opt.find("landscape")!=string::npos) landscape = true;
    bool footnotesize=false; if(opt.find("footnotesize")!=string::npos)footnotesize=true;;
    //
    ofstream out;
    ofstream texout;
    ofstream out_cat;
    ofstream texout_cat;
    gSystem->mkdir(fFitName.c_str());
    gSystem->mkdir((fFitName+"/Tables").c_str());
    if(isPostFit){
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.tex").c_str());
        }
    }
    else{
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.tex").c_str());
        }
    }
    Sample *s = nullptr;
    SampleHist *sh = nullptr;
    SystematicHist *syh = nullptr;


    out << " | ";
    if (standalone) {
        texout << "\\documentclass[10pt]{article}" << endl;
        texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
        texout << "\\begin{document}" << endl;
    }
    if (landscape) texout << "\\begin{landscape}" << endl;
    texout << "\\begin{table}[htbp]" << endl;
    texout << "\\begin{center}" << endl;
    if (footnotesize) texout << "\\footnotesize" << endl;
    texout << "\\begin{tabular}{|c" ;

    if(doCategory){
        out_cat << " | ";
        if (standalone) {
            texout_cat << "\\documentclass[10pt]{article}" << endl;
            texout_cat << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
            texout_cat << "\\begin{document}" << endl;
        }
        texout_cat << "\\begin{table}[htbp]" << endl;
        texout_cat << "\\begin{center}" << endl;
        texout_cat << "\\begin{tabular}{|c" ;
    }


    float Ncol = 2.;
    for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
        sh = fSampleHists[i_smp];
        s = sh->fSample;
        if(s->fType==Sample::DATA) continue;
        if(s->fType==Sample::GHOST) continue;
        texout << "|c";
        if(doCategory) texout_cat << "|c";
        Ncol+=1;
    }
    texout << "|}" << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        texout_cat << "|}" << endl;
        texout_cat << "\\hline " << endl;
    }
    //
    float i_col = 1;
    for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
        sh = fSampleHists[i_smp];
        s = sh->fSample;
        if(s->fType==Sample::DATA) continue;
        if(s->fType==Sample::GHOST) continue;
        std::string title = s->fTitle;
        if(s->fTexTitle!="") title = s->fTexTitle;
        out << "      | " << s->fTitle;
        texout << "      & " << title;
        if(doCategory){
            out_cat << "      | " << s->fTitle;
            texout_cat << "      & " << title;
        }
        i_col+=1;
    }
    out << " |" << endl;
    texout << " \\\\ " << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        out_cat << " |" << endl;
        texout_cat << " \\\\ " << endl;
        texout_cat << "\\hline " << endl;
    }

    for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
        if(TRExFitter::SYSTMAP[fSystNames[i_syst]]!="") out << " | " << TRExFitter::SYSTMAP[fSystNames[i_syst]];
        else                                           out << " | " << fSystNames[i_syst];
        if(TRExFitter::SYSTTEX[fSystNames[i_syst]]!="")      texout << "  " << TRExFitter::SYSTTEX[fSystNames[i_syst]];
        else if(TRExFitter::SYSTMAP[fSystNames[i_syst]]!=""){
            string fixedTitle = TRExFitter::SYSTMAP[fSystNames[i_syst]];
            fixedTitle = ReplaceString(fixedTitle,"#geq","$\\geq$");
            texout << "  " << fixedTitle;
        }
        else                                                texout << " " << fSystNames[i_syst];
        for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
            sh = fSampleHists[i_smp];
            s = sh->fSample;
            if(s->fType==Sample::DATA) continue;
            if(s->fType==Sample::GHOST) continue;
            syh = sh->GetSystematic(fSystNames[i_syst]);
            if(syh==nullptr){
                out << " |    nan   ";
                texout << " &    nan   ";
            }
            else{
              double normUp;
              double normDown;
              if(isPostFit){
                normUp   = (syh->fHistUp_postFit->Integral()   - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
                normDown = (syh->fHistDown_postFit->Integral() - sh->fHist_postFit->Integral()) / sh->fHist_postFit->Integral();
              }
              else{
                normUp = syh->fNormUp;
                normDown = syh->fNormDown;
              }

              out << " | " << normUp;
              texout << setprecision(3) << " & " << normUp;
              out << " / " << normDown;
              texout << setprecision(3) << " / " << normDown;
              //                 pt[i_smp]->AddText(Form("%.2f / %.2f",syh->fNormUp,syh->fNormDown) );
            }
        }

        out << " |" << endl;
        texout << " \\\\ " << endl;
    }


    if(doCategory){
        //--- Systematic tables per category:
        //--- Get systematic names per category and sample
        std::set<std::string> category_names;
        std::map<std::string, std::map<std::string, std::vector<std::string> > > category_syst_names;
        for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
            sh = fSampleHists[i_smp];
            s = sh->fSample;
            for(int i_samplesyst=0; i_samplesyst<(int)s->fSystematics.size();i_samplesyst++){
                if(s->fType==Sample::DATA) continue;
                if(s->fType==Sample::GHOST) continue;
                std::string category = s->fSystematics[i_samplesyst]->fCategory;
                if (category!=""){
                    category_names.insert(category);
                    std::vector<std::string> systePerSample;
                    // Check for systematics existing not in all regions...
                    if (std::find(fSystNames.begin(), fSystNames.end(), s->fSystematics.at(i_samplesyst)->fName)!=fSystNames.end()){
                        std::vector<std::string> sample_syste = category_syst_names[s->fName][category];
                        if (std::find(sample_syste.begin(), sample_syste.end(), s->fSystematics.at(i_samplesyst)->fName) == sample_syste.end()){
                            category_syst_names[s->fName][category].push_back(s->fSystematics.at(i_samplesyst)->fName);
                        }
                    }
                }
            }
        }
        //--- Loop over categories and samples, to get systematic effects using BuildTotError function
        for (auto category : category_names){
            out_cat << " | " << category;
            texout_cat << " " << category;
            for(int i_smp=0;i_smp<(int)fSampleHists.size();i_smp++){
                sh = fSampleHists[i_smp];
                s = sh->fSample;
                if(s->fType==Sample::DATA) continue;
                if(s->fType==Sample::GHOST) continue;

                std::vector<TH1*> category_histo_up;
                std::vector<TH1*> category_histo_down;
                std::vector<std::string> sample_syste = category_syst_names[s->fName][category];
                for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
                    if(!sh->HasSyst(fSystNames[i_syst]))
                      continue;

                    if (std::find(sample_syste.begin(), sample_syste.end(), fSystNames.at(i_syst)) == sample_syste.end())
                      continue;

                    syh = sh->GetSystematic(fSystNames[i_syst]);

                    if(isPostFit){
                        TH1 *h_up = (TH1*) syh->fHistUp_postFit;
                        h_up->SetDirectory(0);
                        TH1 *h_down = (TH1*) syh->fHistDown_postFit;
                        h_down->SetDirectory(0);
                        category_histo_up.push_back(h_up);
                        category_histo_down.push_back(h_down);
                    }
                    else{
                        TH1 *h_up = (TH1*) syh->fHistUp;
                        h_up->SetDirectory(0);
                        TH1 *h_down = (TH1*) syh->fHistDown;
                        h_down->SetDirectory(0);
                        category_histo_up.push_back(h_up);
                        category_histo_down.push_back(h_down);
                    }
                }

                TGraphAsymmErrors *g_err;
                double err = 0.;
                if (isPostFit){
                    g_err = BuildTotError(sh->fHist_postFit, category_histo_up, category_histo_down, category_syst_names[s->fName][category], fitRes->fCorrMatrix);
                    if (category_histo_up.size()>0 && sh->fHist_postFit->Integral()>0.){
                        for (int ibin=1; ibin<sh->fHist_postFit->GetNbinsX()+1; ibin++){
                            if (pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist_postFit->GetBinError(ibin), 2)>=0.){ // dummy check
                                //Need to substract the statistical unc.
                                err += (sqrt(pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist_postFit->GetBinError(ibin), 2)))/sh->fHist_postFit->Integral();
                            }
                        }
                    }
                }
                else{
                    g_err = BuildTotError(sh->fHist, category_histo_up, category_histo_down, category_syst_names[s->fName][category]);
                    if (category_histo_up.size()>0 && sh->fHist->Integral()>0.){
                        for (int ibin=1; ibin<sh->fHist->GetNbinsX()+1; ibin++){
                            if (pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist->GetBinError(ibin), 2)>=0.){
                                //Need to substract the statistical unc.
                                err += (sqrt(pow(g_err->GetErrorYhigh(ibin-1), 2) - pow(sh->fHist->GetBinError(ibin), 2)))/sh->fHist->Integral();
                            }
                        }
                    }
                }
                out_cat << setprecision(3) << " | " << err;
                texout_cat << setprecision(3) << " & " << err;
            }
            out_cat << " |" << endl;
            texout_cat << " \\\\ " << endl;
        }
    }

    texout << "\\hline " << endl;
    texout << "\\end{tabular} " << endl;
    texout << "\\caption{Relative effect of each systematic on the yields.} " << endl;
    texout << "\\end{center} " << endl;
    texout << "\\end{table} " << endl;

    if(doCategory){
        texout_cat << "\\hline " << endl;
        texout_cat << "\\end{tabular} " << endl;
        texout_cat << "\\caption{Realtive effect of each group of systematics on the yields.} " << endl;
        texout_cat << "\\end{center} " << endl;
        texout_cat << "\\end{table} " << endl;
    }
    if (landscape) texout << "\\end{landscape}" << endl;
    if (standalone) {
        texout << "\\end{document}" << endl;
    }

    if(doClean){
        std::string shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
        if(isPostFit) shellcommand += "_postFit";
        shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
        shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst_clean.tex";
        gSystem->Exec(shellcommand.c_str());
        if(doCategory){
            shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
            shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += "_clean.tex";
            gSystem->Exec(shellcommand.c_str());
        }
    }
}

// --------------- Functions --------------- //

double GetDeltaN(double alpha, double Iz, double Ip, double Imi, int intCode){
    // protection against negative values
    if(Ip<=0)  Ip  = 0.00000001*Iz;
    if(Imi<=0) Imi = 0.00000001*Iz;

    double deltaN = 0.;
    if(alpha>0)      deltaN = Ip;
    else if(alpha<0) deltaN = Imi;
    else             return 1.;

    if(intCode==4){

        //////////////////////////////////////////////////////////////
        // NORMALISATION
        // =============
        // Exponential-polynomial
        // Equation solved with Mathematica
        //////////////////////////////////////////////////////////////

        if(TMath::Abs(alpha)>1){
            deltaN /= Iz; // divide h_tmp by the nominal
            deltaN = pow( deltaN, TMath::Abs(alpha) );  // d -> d^(|a|)
        } else {
            double logImiIz = TMath::Log(Imi/Iz);
            double logImiIzSqr = logImiIz*logImiIz;
            double logIpIz = TMath::Log(Ip/Iz);
            double logIpIzSqr = logIpIz*logIpIz;
            // polinomial: equations solved with Mathematica
            double a1 = -(15*Imi - 15*Ip - 7*Imi*logImiIz + Imi*logImiIzSqr + 7*Ip*logIpIz - Ip*logIpIzSqr)/(16.*Iz);
            double a2 = -3 + (3*Imi)/(2.*Iz) + (3*Ip)/(2.*Iz) - (9*Imi*logImiIz)/(16.*Iz) + (Imi*logImiIzSqr)/(16.*Iz) -
            (9*Ip*logIpIz)/(16.*Iz) + (Ip*logIpIzSqr)/(16.*Iz);
            double a3 = (5*Imi)/(8.*Iz) - (5*Ip)/(8.*Iz) - (5*Imi*logImiIz)/(8.*Iz) + (Imi*logImiIzSqr)/(8.*Iz) + (5*Ip*logIpIz)/(8.*Iz) -
            (Ip*logIpIzSqr)/(8.*Iz);
            double a4 = 3 - (3*Imi)/(2.*Iz) - (3*Ip)/(2.*Iz) + (7*Imi*logImiIz)/(8.*Iz) -
            (Imi*logImiIzSqr)/(8.*Iz) + (7*Ip*logIpIz)/(8.*Iz) - (Ip*logIpIzSqr)/(8.*Iz);
            double a5 = (-3*Imi)/(16.*Iz) + (3*Ip)/(16.*Iz) + (3*Imi*logImiIz)/(16.*Iz) - (Imi*logImiIzSqr)/(16.*Iz) -
            (3*Ip*logIpIz)/(16.*Iz) + (Ip*logIpIzSqr)/(16.*Iz);
            double a6 = -1 + Imi/(2.*Iz) + Ip/(2.*Iz) - (5*Imi*logImiIz)/(16.*Iz) + (Imi*logImiIzSqr)/(16.*Iz) - (5*Ip*logIpIz)/(16.*Iz) +
            (Ip*logIpIzSqr)/(16.*Iz);
            double a = alpha;
            deltaN = 1 + a1*a + a2*a*a + a3*a*a*a + a4*a*a*a*a + a5*a*a*a*a*a + a6*a*a*a*a*a*a;
        }

    }
    else if (intCode==0) {

        //////////////////////////////////////////////////////////////
        // SHAPES
        // ======
        // Linear-polynomial
        // Extracted and adapted from RooFit
        //////////////////////////////////////////////////////////////

        if(TMath::Abs(alpha)>1){
            deltaN = 1. + TMath::Abs(alpha)*(deltaN - Iz)/Iz;
        } else {
            double eps_plus = Ip - Iz;
            double eps_minus = Iz - Imi;
            double S = 0.5 * (eps_plus + eps_minus);
            double A = 0.0625 * (eps_plus - eps_minus);
            double val = Iz + alpha * (S + alpha * A * ( 15 + alpha * alpha * (-10 + alpha * alpha * 3  ) ) );
            if (val < 0) val = 0.;
            if(Iz == Iz && Iz>0) deltaN = val/Iz;
            else deltaN = 0.;
        }
    }

    if(deltaN!=deltaN) deltaN = 1;  // to avoid nan
    if(deltaN<=0) deltaN = 0; //protection against negative values (can happen for linear extrapolation)

    return deltaN;
}

//___________________________________________________________
//
std::map < int , double > GetDeltaNForUncertainties(double alpha, double alpha_errUp, double alpha_errDown, double Iz, double Ip, double Imi, int intCode){
    double nominal = GetDeltaN(alpha, Iz, Ip, Imi, intCode);
    double up = GetDeltaN(alpha+alpha_errUp, Iz, Ip, Imi, intCode);
    double down = GetDeltaN(alpha+alpha_errDown, Iz, Ip, Imi, intCode);
    return {{-1,down},{0,nominal},{1,up}};
}


//___________________________________________________________
// function to get pre/post-fit agreement
std::pair<double,int> GetChi2Test( TH1* h_data, TH1* h_nominal, std::vector< TH1* > h_up, std::vector< string > fSystNames, CorrelationMatrix *matrix ){
    const unsigned int nbins = h_nominal->GetNbinsX();
    int ndf = 0;
    for(unsigned int i=0;i<nbins;++i){
        const double ydata_i = h_data->GetBinContent(i+1);
        if(ydata_i<0) continue; // skip dropped / blinded bins
        ndf ++;
    }
    //
    //Speed Up: remove irrelevant systematics (which would give in any case 0 correlation)
    std::vector< string > EffectiveSystNames;
    std::vector< unsigned int > EffectiveSystIndex;
    for(unsigned int n=0;n<fSystNames.size();++n){
      if(matrix!=nullptr){
        if (matrix->fNuisParIsThere[fSystNames[n]]) {
           EffectiveSystNames.push_back(fSystNames[n]);
           EffectiveSystIndex.push_back(n);
        }
        else WriteDebugStatus("GetChi2Test"," will skip syst. "+fSystNames[n]);
      }
      else {
        EffectiveSystNames.push_back(fSystNames[n]);
        EffectiveSystIndex.push_back(n);
      }
    }
    const unsigned int nsyst=EffectiveSystNames.size();
    //
    TMatrixD C(ndf,ndf);
    //
    int ibin = 0;
    int jbin = 0;
    for(unsigned int i=0;i<nbins;++i){
        double ydata_i = h_data->GetBinContent(i+1);
        if(ydata_i<0) continue; // skip dropped / blinded bins
        const double ynom_i = h_nominal->GetBinContent(i+1);
        jbin = 0;
        for(unsigned int j=0;j<nbins;++j){
            double sum = 0.;
            const double ydata_j = h_data->GetBinContent(j+1);
            if(ydata_j<0) continue; // skip dropped / blinded bins
            const double ynom_j = h_nominal->GetBinContent(j+1);
            for(unsigned int n=0;n<nsyst;++n){ //n!=m, run only across correlated systs
                if (EffectiveSystNames[n].find("saturated_model") != std::string::npos) continue;
                const double ysyst_i_n = h_up[EffectiveSystIndex[n]]->GetBinContent(i+1);
                for(unsigned int m=0;m<nsyst;++m){
                    if (n==m) continue;
                    const double ysyst_j_m = h_up[EffectiveSystIndex[m]]->GetBinContent(j+1);
                    // more than Bill's suggestion: add correlation between systematics!!
                    double corr = 0.;
                    if(n==m) corr = 1.;
                    else if(matrix!=nullptr) corr = matrix->GetCorrelation(EffectiveSystNames[n],EffectiveSystNames[m]);
                    else continue;
                    sum += (ysyst_i_n-ynom_i) * corr * (ysyst_j_m-ynom_j);
                }
            }
            for(unsigned int n=0;n<fSystNames.size();++n){ // n==m, all systs, corr=1
                if (fSystNames[n].find("saturated_model") != std::string::npos) continue;
                const double ysyst_i_n = h_up[n]->GetBinContent(i+1), ysyst_j_n = h_up[n]->GetBinContent(j+1);
                sum += (ysyst_i_n-ynom_i) /* * 1.0 */ * (ysyst_j_n-ynom_j);
            }
            if(i==j && ynom_i>0) sum += ynom_i;  // add stat uncertainty to diagonal
            if(i==j && ynom_i>0) sum += h_nominal->GetBinError(i+1) * h_nominal->GetBinError(i+1); // add MC stat as well
            C[ibin][jbin] = sum;
            ++jbin;
        }
        ++ibin;
    }
    //
    if(TRExFitter::DEBUGLEVEL > 1) C.Print();
    //
    // Invert the matrix
    C.Invert();
    if(TRExFitter::DEBUGLEVEL > 1) C.Print();
    //
    double chi2 = 0.;
    ibin = 0;
    jbin = 0;
    for(unsigned int i=0;i<nbins;++i){
        const double ydata_i = h_data->GetBinContent(i+1);
        if(ydata_i<0) continue; // skip dropped / blinded bins
        const double ynom_i = h_nominal->GetBinContent(i+1);
        jbin = 0;
        for(unsigned int j=0;j<nbins;j++){
            const double ydata_j = h_data->GetBinContent(j+1);
            if(ydata_j<0) continue; // skip dropped / blinded bins
            const double ynom_j = h_nominal->GetBinContent(j+1);
            chi2 += (ydata_i - ynom_i)*C[ibin][jbin]*(ydata_j-ynom_j);
            ++jbin;
        }
        ++ibin;
    }
    //
    return std::make_pair(chi2,ndf);
}

//___________________________________________________________
// function to build total error band from:
// - a nominal histo (tot exepcted)
// - syst variation histos (eventually already scaled by post-fit pulls)
// - correlation matrix
// Note: if matrix = nullptr => no correlation, i.e. matrix = 1 (used for pre-fit, or to neglect correlation)
TGraphAsymmErrors* BuildTotError( TH1* h_nominal, std::vector< TH1* > h_up, std::vector< TH1* > h_down, std::vector< string > fSystNames, CorrelationMatrix *matrix ){
    if(!h_nominal){
        WriteErrorStatus("BuildTotError","h_nominal not defined.");
        exit(EXIT_FAILURE);
    }
    if(h_up.size()!=h_down.size()){
        WriteErrorStatus("BuildTotError","h_up and h_down have different size.");
        exit(EXIT_FAILURE);
    }
    if(h_up.size()!=fSystNames.size()){
        WriteErrorStatus("BuildTotError","h_up and fSystNames have different size.");
        exit(EXIT_FAILURE);
    }
    // FIXME FIXME FIXME...
    //
    //Speed Up: remove irrelevant systematics (which would give in any case 0 correlation)
    std::vector< string > EffectiveSystNames;
    std::vector< unsigned int > EffectiveSystIndex;
    for(unsigned int n=0;n<fSystNames.size();++n){
      if(matrix!=nullptr){
        if (matrix->fNuisParIsThere[fSystNames[n]]) {
           EffectiveSystNames.push_back(fSystNames[n]);
           EffectiveSystIndex.push_back(n);
        }
      }
      else {
          EffectiveSystNames.push_back(fSystNames[n]);
          EffectiveSystIndex.push_back(n);
      }
    }
    //
    TGraphAsymmErrors *g_totErr = new TGraphAsymmErrors( h_nominal );
    float finalErrPlus(0.);
    float finalErrMinus(0.);
    float corr(0.);
    float errUp_i(0.),   errUp_j(0.);
    float errDown_i(0.), errDown_j(0.);
    //
    // - loop on bins
    for(int i_bin=1;i_bin<h_nominal->GetNbinsX()+1;i_bin++){
        finalErrPlus = 0;
        finalErrMinus = 0;
        corr = 0;
        // yieldNominal = h_nominal->GetBinContent(i_bin);
        // - loop on the syst, two by two, to include the correlations
        for(unsigned int i_syst=0;i_syst<EffectiveSystNames.size();i_syst++){
            for(unsigned int j_syst=0;j_syst<EffectiveSystNames.size();j_syst++){
                if (i_syst==j_syst) continue;
                if(matrix!=nullptr){
                    corr = matrix->GetCorrelation(EffectiveSystNames[i_syst],EffectiveSystNames[j_syst]);
                }
                else{
                    if(EffectiveSystNames[i_syst]==EffectiveSystNames[j_syst]) corr = 1.;
                    else                                                       corr = 0.;
                }
                errUp_i   = h_up[EffectiveSystIndex[i_syst]]  ->GetBinContent(i_bin);// - yieldNominal;
                errDown_i = h_down[EffectiveSystIndex[i_syst]]->GetBinContent(i_bin);// - yieldNominal;
                errUp_j   = h_up[EffectiveSystIndex[j_syst]]  ->GetBinContent(i_bin);// - yieldNominal;
                errDown_j = h_down[EffectiveSystIndex[j_syst]]->GetBinContent(i_bin);// - yieldNominal;

                //
                // Symmetrize (seems to be done in Roostats ??)
                //
                double err_i = (errUp_i - errDown_i)/2.;
                double err_j = (errUp_j - errDown_j)/2.;

                //
                // Compute the + and - variations
                //
                finalErrPlus  += err_i * err_j * corr;
                finalErrMinus += err_i * err_j * corr;
            }
        }
        // now all diagonal el. of all systematics, corr = 1;
        for(unsigned int i_syst=0;i_syst<fSystNames.size();i_syst++){
            errUp_i   = h_up[i_syst]  ->GetBinContent(i_bin);// - yieldNominal;
            errDown_i = h_down[i_syst]->GetBinContent(i_bin);// - yieldNominal;

            //
            // Symmetrize (seems to be done in Roostats ??)
            //
            double err_i = (errUp_i - errDown_i)/2.;

            //
            // Compute the + and - variations
            //
            finalErrPlus  += err_i * err_i ;
            finalErrMinus += err_i * err_i ;
        }
        // add stat uncertainty, which should have been stored as orignal bin errors in the h_nominal (if fUseStatErr is true)
        finalErrPlus  += pow( h_nominal->GetBinError(i_bin), 2 );
        finalErrMinus += pow( h_nominal->GetBinError(i_bin), 2 );

        g_totErr->SetPointEYhigh(i_bin-1,sqrt(TMath::Abs(finalErrPlus )));
        g_totErr->SetPointEYlow( i_bin-1,sqrt(TMath::Abs(finalErrMinus)));
    }

    return g_totErr;
}

//___________________________________________________________
//
void Region::PrepareMorphScales(FitResults *fitRes, std::vector<double> *morph_scale, std::vector<double> *morph_scale_nominal) const{
    for(int i=0;i<fNSamples;i++){
        // skip data
        if(fSampleHists[i]->fSample->fType==Sample::DATA) continue;
        if(fSampleHists[i]->fSample->fType==Sample::GHOST) continue;
        for(int i_syst=0;i_syst<(int)fSystNames.size();i_syst++){
            std::string systName    = fSystNames[i_syst];
            if(TRExFitter::NPMAP[systName]=="") TRExFitter::NPMAP[systName] = systName;

            if(fSampleHists[i]->HasNorm(fSystNames[i_syst])){
                // if this norm factor is a morphing one
                if(fSystNames[i_syst].find("morph_")!=string::npos || fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fExpression.first!=""){
                    std::string formula = TRExFitter::SYSTMAP[fSystNames[i_syst]];
                    std::string name = TRExFitter::NPMAP[fSystNames[i_syst]];
                    WriteDebugStatus("Region::PrepareMorphScales", "formula: " +formula);
                    WriteDebugStatus("Region::PrepareMorphScales", "name: " +name);
                    std::vector < std::pair < std::string,std::vector<double> > > nameS;
                    if(fSystNames[i_syst].find("morph_")!=std::string::npos){
                        nameS.push_back(std::make_pair(name,std::vector<double>{double(fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fNominal),
                            double(fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fMin),double(fSampleHists[i]->GetNormFactor(fSystNames[i_syst])->fMax)}));
                    }
                    else{
                        nameS = processString(name);
                    }
                    std::vector <double> nfValuevec;
                    for (unsigned int j = 0; j<nameS.size(); j++){
                        formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                        nfValuevec.push_back(fitRes->GetNuisParValue(nameS[j].first));
                    }
                    double *nfValue = nfValuevec.data();
                    WriteDebugStatus("Region::PrepareMorphScales", "formula: " +formula);
                    for(unsigned int j = 0; j<nameS.size(); j++){
                        WriteDebugStatus("Region::PrepareMorphScales", "nfValue["+std::to_string(j)+"]: "+std::to_string(nfValue[j]));
                    }
                    TFormula f_morph ("f_morph",formula.c_str());
                    float scaleNom = f_morph.EvalPar(nfValue,nullptr);
                    morph_scale->emplace_back(scaleNom);
                }
            }
            for(unsigned int i_nf=0;i_nf<fSampleHists[i]->fSample->fNormFactors.size();i_nf++){
                NormFactor *nf = fSampleHists[i]->fSample->fNormFactors[i_nf];
                // if this norm factor is a morphing one
                if(nf->fName.find("morph_")!=string::npos || nf->fExpression.first!=""){
                    std::string formula = TRExFitter::SYSTMAP[nf->fName];
                    std::string name = TRExFitter::NPMAP[nf->fName];
                    WriteDebugStatus("Region::PrepareMorphScales", "formula: " +formula);
                    WriteDebugStatus("Region::PrepareMorphScales", "name: " +name);
                    std::vector < std::pair < std::string,std::vector<double> > > nameS;
                    if(nf->fName.find("morph_")!=std::string::npos){
                        nameS.push_back(std::make_pair(name,std::vector<double>{double(nf->fNominal),double(nf->fMin),double(nf->fMax)}));
                    }
                    else{
                        nameS = processString(name);
                    }
                    std::vector <double> nfNominalvec;
                    for (unsigned int j = 0; j<nameS.size(); j++){
                        formula = ReplaceString(formula,nameS[j].first,"x["+std::to_string(j)+"]");
                        nfNominalvec.push_back(nameS[j].second[0]);
                    }
                    double *nfNominal = nfNominalvec.data();
                    WriteDebugStatus("Region::PrepareMorphScales", "formula: " +formula);
                    for(unsigned int j = 0; j<nameS.size(); j++){
                        WriteDebugStatus("Region::PrepareMorphScales", "nfNominal["+std::to_string(j)+"]: "+std::to_string(nfNominal[j]));
                    }
                    TFormula f_morph("f_morph",formula.c_str());
                    float scaleNom = f_morph.EvalPar(nfNominal,nullptr);
                    morph_scale_nominal->emplace_back(scaleNom);
                }
            }
        }
    }
}

//___________________________________________________________
//
void Region::SystPruning(PruningUtil *pu){
    TH1* hTot = nullptr;
    if(pu->fStrategy==1){
        hTot = GetTotHist(false); // don't include signal
    }
    else if(pu->fStrategy==2){
        hTot = GetTotHist(true); // include signal
    }
    for(auto sh : fSampleHists){
        sh->SystPruning(pu,hTot);
        //
        // flag overall systematics as no shape also for pruning purposes
        for(auto syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            if(syh->fSystematic->fType==Systematic::OVERALL){
                syh->fShapePruned = true;
            }
        }
        //
        // add by-hand dropping of shape or norm
        for(auto syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            if( FindInStringVector(syh->fSystematic->fDropShapeIn,fName)>=0 ){
                syh->fShapePruned = true;
            }
            if( FindInStringVector(syh->fSystematic->fDropNormIn,fName)>=0 ){
                syh->fNormPruned = true;
            }
        }
    }
    //
    // reference pruning
    for(auto sh : fSampleHists){
        for(auto syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            Systematic *syst = syh->fSystematic;
            if(syst->fReferencePruning == "") continue;
            SampleHist *refSmpH = GetSampleHist(syst->fReferencePruning);
            if(refSmpH==nullptr){
                WriteWarningStatus("Region::SystPruning", "Cannot find reference pruning sample: " + syst->fReferencePruning + " in region: " + fName);
                continue;
            }
            SystematicHist *refSysH = refSmpH->GetSystematic(syst->fName);
            if(refSysH==nullptr){
                WriteWarningStatus("Region::SystPruning", "Cannot find systematic " + syst->fName + " for reference pruning sample " + syst->fReferencePruning);
                continue;
            }
            syh->fNormPruned = refSysH->fNormPruned;
            syh->fShapePruned = refSysH->fShapePruned;
            syh->fBadShape = refSysH->fBadShape;
            syh->fBadNorm = refSysH->fBadNorm;
        }
    }
}

//___________________________________________________________
//
TH1* Region::GetTotHist(bool includeSignal){
    TH1* hTot = nullptr;
    for(auto sh : fSampleHists){
        if(!sh->fSample) continue;
        if(sh->fSample->fType==Sample::GHOST) continue;
        if(sh->fSample->fType==Sample::DATA) continue;
        if(!includeSignal && sh->fSample->fType==Sample::SIGNAL) continue;
        if(!sh->fHist) continue;
        TH1* hTmp = (TH1*)sh->fHist->Clone(("hTot_"+fName).c_str());
        // scale accoring to nominal SF (considering morphing as well)
        const float& scale = GetNominalMorphScale(sh);
        hTmp->Scale(scale);
        if(!hTot) hTot = hTmp;
        else{
            hTot->Add(hTmp);
            delete hTmp;
        }
    }
    return hTot;
}
