// Class include
#include "TRExFitter/TRExFit.h"

// Framework includes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/LikelihoodScanManager.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExPlot.h"
#include "TRExFitter/RankingManager.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/PruningUtil.h"
#include "TRExFitter/TruthSample.h"
#include "TRExFitter/UnfoldingSample.h"
#include "TRExFitter/UnfoldingSystematic.h"
#include "TRExFitter/YamlConverter.h"
#include "TRExFitter/CorrelationMatrix.h"

// UnfoldingCode includes
#include "UnfoldingCode/UnfoldingCode/UnfoldingTools.h"
#include "UnfoldingCode/UnfoldingCode/UnfoldingResult.h"

// CommonStatTiils includes
#include "CommonStatTools/runSig.h"
#include "CommonStatTools/runAsymptoticsCLs.h"

//Roofit headers
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"

//HistFactory headers
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

// ATLAS stuff
#include "AtlasUtils/AtlasStyle.h"
#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

// ROOT includes
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFormula.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPie.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"

// c++ includes
#include <algorithm>
#include <cctype>
#include <iomanip>

using namespace RooFit;

// -------------------------------------------------------------------------------------------------
// class TRExFit

//__________________________________________________________________________________
//
TRExFit::TRExFit(std::string name) :
    fName(name),
    fDir(""),
    fLabel(""),
    fInputFolder(""),
    fInputName(name),
    fFitResultsFile(""),
    fNRegions(0),
    fNSamples(0),
    fPOI(""),
    fPOIunit(""),
    fUseStatErr(false),
    fStatErrThres(0.05),
    fUseGammaPulls(false),
    fLumi(1.),
    fLumiScale(1.),
    fLumiErr(0.000001),
    fThresholdSystPruning_Normalisation(-1),
    fThresholdSystPruning_Shape(-1),
    fThresholdSystLarge(-1),
    fMCweight("1"),
    fSelection("1"),
    fFitResults(nullptr),
    fWithPullTables(false),
    fIntCode_overall(4),
    fIntCode_shape(0),
    fInputType(HIST),
    fSystDataPlot_upFrame(false),
    fStatOnly(false),
    fGammasInStatOnly(false),
    fStatOnlyFit(false),
    fFixNPforStatOnlyFit(false),
    fYmin(0),
    fYmax(0),
    fRatioYmin(0.5),
    fRatioYmax(1.5),
    fRatioYminPostFit(0.5),
    fRatioYmaxPostFit(1.5),
    fRatioYtitle(""),
    fRatioType("DATA/MC"),
    fLumiLabel("XX fb^{-1}"),
    fCmeLabel("13 TeV"),
    fSuffix(""),
    fSaveSuffix(""),
    fUpdate(false),
    fKeepPruning(false),
    fBlindingThreshold(-1),
    fBlindingType(Common::SOVERB),
    fRankingMaxNP(10),
    fRankingOnly("all"),
    fRankingPlot("Merge"),
    fImageFormat("png"),
    fAtlasLabel("Internal"),
    fDoSummaryPlot(true),
    fDoMergedPlot(false),
    fDoTables(true),
    fDoSignalRegionsPlot(true),
    fDoPieChartPlot(true),
    fGroupedImpactCategory("all"),
    fSummaryPrefix(""),
    fFitType(UNDEFINED),
    fFitRegion(CRSR),
    fFitNPValuesFromFitResults(""),
    fInjectGlobalObservables(false),
    fFitPOIAsimov(0),
    fFitIsBlind(false),
    fUseRnd(false),
    fRndRange(0.1),
    fRndSeed(-999),
    fLHscanMin(999999),
    fLHscanMax(-999999),
    fLHscanSteps(30),
    fLHscanMinY(999999),
    fLHscanMaxY(-999999),
    fLHscanStepsY(30),
    fParal2D(false),
    fParal2Dstep(-1),
    fWorkspaceFileName(""),
    fDoGroupedSystImpactTable(false),
    fLimitType(ASYMPTOTIC),
    fLimitIsBlind(false),
    fSignalInjection(false),
    fSignalInjectionValue(0),
    fLimitParamName("parameter"),
    fLimitParamValue(0),
    fLimitOutputPrefixName("myLimit"),
    fLimitsConfidence(0.95),
    fSignificanceIsBlind(false),
    fSignificancePOIAsimov(0),
    fSignificanceParamName("parameter"),
    fSignificanceParamValue(0),
    fSignificanceOutputPrefixName("mySignificance"),
    fCleanTables(false),
    fSystCategoryTables(false),
    fKeepPrefitBlindedBins(false),
    fBlindedBins(nullptr),
    fCustomAsimov(""),
    fTableOptions("STANDALONE"),
    fGetGoodnessOfFit(false),
    fGetChi2(0), // 0: no, 1: stat-only, 2: with syst
    fSmoothOption(HistoTools::SmoothOption::MAXVARIATION),
    fSuppressNegativeBinWarnings(false),
    fTemplateInterpolationOption(TRExFit::LINEAR),
    fBootstrap(""),
    fBootstrapSyst(""),
    fBootstrapSample(""),
    fBootstrapIdx(-1),
    fDecorrSuff("_decor"),
    fDoNonProfileFit(false),
    fNonProfileFitSystThreshold(0),
    fFitToys(0),
    fToysHistoMin(9999),
    fToysHistoMax(-9999),
    fToysHistoNbins(50),
    fToysPseudodataNP(""),
    fToysPseudodataNPShift(1.),
    fSmoothMorphingTemplates(""),
    fPOIPrecision(2),
    fRankingPOIName("#mu"),
    fUseATLASRounding(false),
    fUseATLASRoundingTxt(false),
    fUseATLASRoundingTex(false),
    fuseGammasForCorr(false),
    fPropagateSystsForMorphing(false),
    fPruningType(SEPARATESAMPLE),
    fLabelX(-1),
    fLabelY(-1),
    fLegendX1(-1),
    fLegendX2(-1),
    fLegendY(-1),
    fLabelXSummary(-1),
    fLabelYSummary(-1),
    fLegendX1Summary(-1),
    fLegendX2Summary(-1),
    fLegendYSummary(-1),
    fLabelXMerge(-1),
    fLabelYMerge(-1),
    fLegendX1Merge(-1),
    fLegendX2Merge(-1),
    fLegendYMerge(-1),
    fLegendNColumns(2),
    fLegendNColumnsSummary(3),
    fLegendNColumnsMerge(3),
    fShowRatioPad(true),
    fShowRatioPadSummary(true),
    fShowRatioPadMerge(true),
    fExcludeFromMorphing(""),
    fSaturatedModel(false),
    fDoSystNormalizationPlots(true),
    fDebugNev(-1),
    fCPU(1),
    fMatrixOrientation(FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS),
    fTruthDistributionPath(""),
    fTruthDistributionFile(""),
    fTruthDistributionName(""),
    fNumberUnfoldingTruthBins(0),
    fNumberUnfoldingRecoBins(0),
    fUnfoldingResultMin(0),
    fUnfoldingResultMax(2),
    fHasAcceptance(false),
    fUnfoldingTitleX("X axis"),
    fUnfoldingTitleY("Y axis"),
    fUnfoldingRatioYmax(1.5),
    fUnfoldingRatioYmin(0.5),
    fUnfoldingScaleRangeY(-1.),
    fUnfoldingLogX(false),
    fUnfoldingLogY(false),
    fUnfoldingTitleOffsetX(3.5),
    fUnfoldingTitleOffsetY(1.6),
    fNominalTruthSample("SetMe"),
    fAlternativeAsimovTruthSample(""),
    fMigrationTitleX("X"),
    fMigrationTitleY("Y"),
    fMigrationLogX(false),
    fMigrationLogY(false),
    fMigrationTitleOffsetX(1),
    fMigrationTitleOffsetY(1.5),
    fPlotSystematicMigrations(false),
    fMigrationZmin(0),
    fMigrationZmax(1),
    fResponseZmin(0),
    fResponseZmax(1),
    fMigrationText(false),
    fUnfoldingDivideByBinWidth(false),
    fPruningShapeOption(PruningUtil::SHAPEOPTION::MAXBIN),
    fSummaryLogY(true),
    fUseInFit(true),
    fUseInComparison(true),
    fReorderNPs(false),
    fBlindSRs(false),
    fHEPDataFormat(false),
    fAlternativeShapeHistFactory(false),
    fFitStrategy(-1),
    fBinnedLikelihood(false),
    fRemoveLargeSyst(true),
    fRemoveSystOnEmptySample(false),
    fValidationPruning(false),
    fUnfoldNormXSec(false),
    fUnfoldNormXSecBinN(-1)
{
    TRExFitter::IMAGEFORMAT.emplace_back("png");
}

//__________________________________________________________________________________
//
TRExFit::~TRExFit() {
    delete fFitResults;

    for(auto ireg : fRegions) {
        delete ireg;
    }
}

//__________________________________________________________________________________
//
void TRExFit::SetPOI(std::string name){
    fPOI = name;
}

//__________________________________________________________________________________
//
void TRExFit::SetStatErrorConfig(bool useIt, double thres, std::string cons){
    fUseStatErr = useIt;
    fStatErrThres = thres;
    fStatErrCons = cons;
}

//__________________________________________________________________________________
//
void TRExFit::SetLumiErr(double err){
    fLumiErr = err;
}

//__________________________________________________________________________________
//
void TRExFit::SetLumi(const double lumi){
    fLumi = lumi;
}

//__________________________________________________________________________________
//
void TRExFit::SetFitType(FitType type){
    fFitType = type;
}

//__________________________________________________________________________________
//
void TRExFit::SetLimitType(LimitType type){
    fLimitType = type;
}

//__________________________________________________________________________________
//
void TRExFit::SetFitRegion(FitRegion region){
    fFitRegion = region;
}

//__________________________________________________________________________________
//
std::shared_ptr<Sample> TRExFit::NewSample(const std::string& name,int type){
    fSamples.emplace_back(new Sample(name,type));
    //
    fNSamples++;
    return fSamples.back();
}

//__________________________________________________________________________________
//
std::shared_ptr<Systematic> TRExFit::NewSystematic(const std::string& name){
    fSystematics.emplace_back(new Systematic(name));
    return fSystematics.back();
}

//__________________________________________________________________________________
//
Region* TRExFit::NewRegion(const std::string& name){
    fRegions.push_back(new Region(name));
    //
    fRegions[fNRegions]->fFitName = fName;
    fRegions[fNRegions]->fSuffix = fSuffix;
    fRegions[fNRegions]->fFitLabel = fLabel;
    fRegions[fNRegions]->fFitType = fFitType;
    fRegions[fNRegions]->fPOI = fPOI;
    fRegions[fNRegions]->fIntCode_overall = fIntCode_overall;
    fRegions[fNRegions]->fIntCode_shape   = fIntCode_shape;
    fRegions[fNRegions]->fLumiScale = fLumiScale;
    fRegions[fNRegions]->fBlindingThreshold = fBlindingThreshold;
    fRegions[fNRegions]->fBlindingType = fBlindingType;
    fRegions[fNRegions]->fKeepPrefitBlindedBins = fKeepPrefitBlindedBins;
    fRegions[fNRegions]->fRatioYmax = fRatioYmax;
    fRegions[fNRegions]->fRatioYmin = fRatioYmin;
    fRegions[fNRegions]->fRatioYmaxPostFit = fRatioYmaxPostFit;
    fRegions[fNRegions]->fRatioYminPostFit = fRatioYminPostFit;
    fRegions[fNRegions]->fRatioYtitle = fRatioYtitle;
    fRegions[fNRegions]->fRatioType = fRatioType;
    fRegions[fNRegions]->fLabelX = fLabelX;
    fRegions[fNRegions]->fLabelY = fLabelY;
    fRegions[fNRegions]->fLegendX1 = fLegendX1;
    fRegions[fNRegions]->fLegendX2 = fLegendX2;
    fRegions[fNRegions]->fLegendY = fLegendY;
    fRegions[fNRegions]->fLegendNColumns = fLegendNColumns;
    //
    fNRegions ++;
    return fRegions[fNRegions-1];
}

//__________________________________________________________________________________
//
void TRExFit::AddNtuplePath(const std::string& path){
    fNtuplePaths.push_back(path);
}

//__________________________________________________________________________________
//
void TRExFit::SetMCweight(const std::string &weight){
    fMCweight = weight;
}

//__________________________________________________________________________________
//
void TRExFit::SetSelection(const std::string& selection){
    fSelection = selection;
}

//__________________________________________________________________________________
//
void TRExFit::SetNtupleName(const std::string& name){
    fNtupleNames.clear();
    fNtupleNames.push_back(name);
}

//__________________________________________________________________________________
//
void TRExFit::SetNtupleFile(const std::string& name){
    fNtupleFiles.clear();
    fNtupleFiles.push_back(name);
}

//__________________________________________________________________________________
//
void TRExFit::AddHistoPath(const std::string& path){
    fHistoPaths.push_back(path);
}

//__________________________________________________________________________________
// apply smoothing to systematics
void TRExFit::SmoothSystematics(std::string syst){
    WriteInfoStatus("TRExFit::SmoothSystematics", "-------------------------------------------");
    WriteInfoStatus("TRExFit::SmoothSystematics", "Smoothing and/or Symmetrising Systematic Variations ...");

    for(int i_ch=0; i_ch<fNRegions; ++i_ch){
        //
        // Scale systematics according to smoothing of another region (inter-region smoothing)
        // Loop on previous regions
        for(int j_ch=0; j_ch<i_ch; ++j_ch){
            if(fRegions[i_ch]->fIsBinOfRegion[fRegions[j_ch]->fName]>0){ // NB: these bins have to start with 1, not 0 !! 0 means not filled (due to map implementation)
                int binIdx = fRegions[i_ch]->fIsBinOfRegion[fRegions[j_ch]->fName];
                for(auto& sh : fRegions[i_ch]->fSampleHists){
                    for(auto& syh : sh->fSyst){
                        float scaleUp = 1.;
                        float scaleDown = 1.;
                        // get scale factors to apply according to reference bin of reference region
                        if(fRegions[j_ch]->GetSampleHist(sh->fSample->fName)!=nullptr){
                            std::shared_ptr<SampleHist> sh_ref = fRegions[j_ch]->GetSampleHist(sh->fSample->fName);
                            std::shared_ptr<SystematicHist> syh_ref = sh_ref->GetSystematic(syh->fSystematic->fName);
                            if(syh_ref!=nullptr){
                                float systVarUp_orig = syh_ref->fHistUp_orig->GetBinContent(binIdx);
                                float systVarDown_orig = syh_ref->fHistDown_orig->GetBinContent(binIdx);
                                float systVarUp = syh_ref->fHistUp->GetBinContent(binIdx);
                                float systVarDown = syh_ref->fHistDown->GetBinContent(binIdx);
                                if(systVarUp_orig!=0) scaleUp = systVarUp/systVarUp_orig;
                                else WriteWarningStatus("TRExFit::SmoothSystematics","In inter-region smoothing attempting to divide by zero. Skipping scaling.");
                                if(syh->fSystematic->fSymmetrisationType==HistoTools::SYMMETRIZEONESIDED) scaleDown = scaleUp;
                                else scaleDown = systVarDown/systVarDown_orig;
                            }
                        }
                        else{
                            WriteWarningStatus("TRExFit::SmoothSystematics","Sample not found in region indicated as reference for inter-region smoothing.");
                        }
                        // scale
                        syh->fHistUp->Scale(scaleUp);
                        syh->fHistDown->Scale(scaleDown);
                    }
                }
                continue;
            }
        }

        // collect information which systematics contain reference smoothing samples
        std::vector<std::string> referenceSmoothSysts{};
        for (const auto isyst : fSystematics){
            if (std::find(isyst->fRegions.begin(), isyst->fRegions.end(), fRegions[i_ch]->fName) == isyst->fRegions.end()) continue;
            if (isyst->fReferenceSmoothing != ""){
                referenceSmoothSysts.emplace_back(isyst->fName);
            }
        }

        // if there are no reference smoothing samples, proceed as usual
        if (referenceSmoothSysts.size() == 0){
            for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
                fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(fSmoothOption, fAlternativeShapeHistFactory, syst, false);
            }
        } else {
            std::vector<std::size_t> usedSysts{};
            for (int i_smp=0; i_smp<fRegions[i_ch]->fNSamples; ++i_smp){
                for (std::size_t i_syst = 0; i_syst < fSystematics.size(); ++i_syst){
                    if (fSystematics.at(i_syst) == nullptr) continue;
                    // check only systematics for the samples that are specified
                    if (std::find(fSystematics.at(i_syst)->fSamples.begin(), fSystematics.at(i_syst)->fSamples.end(), fRegions[i_ch]->fSampleHists[i_smp]->GetSample()->fName) == fSystematics.at(i_syst)->fSamples.end()) continue;
                    // take only systematics that belong to this region
                    if (std::find(fSystematics.at(i_syst)->fRegions.begin(), fSystematics.at(i_syst)->fRegions.end(), fRegions[i_ch]->fName) == fSystematics.at(i_syst)->fRegions.end()) continue;
                    if (fSystematics.at(i_syst)->fReferenceSmoothing == "") {
                        // the systemtic is not using special smoothing
                        fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(fSmoothOption, fAlternativeShapeHistFactory, fSystematics.at(i_syst)->fName, true);
                    } else {
                        // check if the syst has been smoothed already
                        if (std::find(usedSysts.begin(), usedSysts.end(), i_syst) != usedSysts.end()) continue;
                        // Need to apply special smoothing
                        // smooth the reference sample
                        std::shared_ptr<SampleHist> sh = GetSampleHistFromName(fRegions[i_ch], fSystematics.at(i_syst)->fReferenceSmoothing);
                            if (sh == nullptr){
                            WriteErrorStatus("TRExFit::SmoothSystematics","Cannot find ReferenceSmoothing in the list of samples!");
                            exit(EXIT_FAILURE);
                        }

                        std::unique_ptr<TH1> nominal_cpy = nullptr;
                        std::unique_ptr<TH1> up_cpy = nullptr;
                        std::unique_ptr<TH1> down_cpy = nullptr;

                        int systIndex = -1;
                        // smooth on the sample that is specified in ReferenceSmoothing
                        for (int i_sample=0; i_sample<fRegions[i_ch]->fNSamples; ++i_sample){
                            if (fRegions[i_ch]->fSampleHists[i_sample]->GetSample()->fName == fSystematics.at(i_syst)->fReferenceSmoothing){
                                sh->SmoothSyst(fSmoothOption, fAlternativeShapeHistFactory, fSystematics.at(i_syst)->fName, true);

                                // save the smoothed histograms
                                nominal_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(fRegions[i_ch]->fSampleHists[i_sample]->fHist->Clone()));
                                systIndex = GetSystIndex(fRegions[i_ch]->fSampleHists[i_sample].get(), fSystematics.at(i_syst)->fName);
                                if (systIndex < 0){
                                    WriteWarningStatus("TRExFit::SmoothSystematics", "Cannot find systematic in the list wont smooth!");
                                    return;
                                }
                                up_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistUp->Clone()));
                                down_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistDown->Clone()));
                                break;
                            }
                        }

                        // finally, apply the same smoothing to all other samples, bin-by-bin
                        for (int i_sample=0; i_sample<fRegions[i_ch]->fNSamples; ++i_sample){
                            // skip samples that do not belong to this systematics
                            if (std::find(fSystematics.at(i_syst)->fSamples.begin(), fSystematics.at(i_syst)->fSamples.end(), fRegions[i_ch]->fSampleHists[i_sample]->GetSample()->fName) ==
                                fSystematics.at(i_syst)->fSamples.end()) continue;
                            // skip the one that has already been smoothed, the ReferenceSmoothing
                            if (fRegions[i_ch]->fSampleHists[i_sample]->GetSample()->fName == fSystematics.at(i_syst)->fReferenceSmoothing) continue;

                            if (systIndex < 0){
                                WriteWarningStatus("TRExFit::SmoothSystematics", "Cannot find systematic in the list wont smooth!");
                                return;
                            }
                            fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistUp.reset(CopySmoothedHisto(fRegions[i_ch]->fSampleHists[i_sample].get(),nominal_cpy.get(),up_cpy.get(),down_cpy.get(),true));
                            fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistDown.reset(CopySmoothedHisto(fRegions[i_ch]->fSampleHists[i_sample].get(),nominal_cpy.get(),up_cpy.get(),down_cpy.get(),false));
                        }

                        usedSysts.emplace_back(i_syst);
                    }
                } // loop over systs
            } // loop over samples
        }
    }
}

//
// Try to split root file creation and histogram wiriting
//__________________________________________________________________________________
// create new root file(s)
void TRExFit::CreateRootFiles(){
    bool recreate = !fUpdate;
    gSystem->mkdir( fName.c_str());
    gSystem->mkdir( (fName + "/Histograms/").c_str() );
    std::string fileName("");
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    //
    if(singleOutputFile){
        if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos" + fSaveSuffix + ".root";
        else                 fileName = fName + "/Histograms/" + fInputName + "_histos" + fSaveSuffix + ".root";
        // Bootstrap
        if(fBootstrap!="" && fBootstrapIdx>=0){
            fileName = Common::ReplaceString(fileName,"_histos.root",Form("_histos__%s%d.root",fBootstrapSample.c_str(),fBootstrapIdx));
        }
        WriteInfoStatus("TRExFit::CreateRootFiles","-------------------------------------------");
        WriteInfoStatus("TRExFit::CreateRootFiles","Creating/updating file " + fileName + " ...");
        if(recreate) fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"RECREATE")));
        else         fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"UPDATE")));
        TRExFitter::TFILEMAP.insert(std::make_pair(fileName,fFiles.back()));
    }
    else{
        for(int i_ch=0;i_ch<fNRegions;i_ch++){
            if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
            else                 fileName = fName + "/Histograms/" + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
            // Bootstrap
            if(fBootstrap!="" && fBootstrapIdx>=0){
                fileName = Common::ReplaceString(fileName,"_histos.root",Form("_histos__%s%d.root",fBootstrapSample.c_str(),fBootstrapIdx));
            }
            WriteInfoStatus("TRExFit::CreateRootFiles","-------------------------------------------");
            WriteInfoStatus("TRExFit::CreateRootFiles","Creating/updating file " + fileName + " ...");
            if(recreate) fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"RECREATE")));
            else         fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"UPDATE")));
            TRExFitter::TFILEMAP.insert(std::make_pair(fileName,fFiles.back()));
        }
    }
}

//__________________________________________________________________________________
// fill files with all the histograms
void TRExFit::WriteHistos(bool reWriteOrig) const{
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    std::string fileName("");
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        //
        if(singleOutputFile) fileName = fFiles[0]   ->GetName();
        else                 fileName = fFiles[i_ch]->GetName();
        WriteInfoStatus("TRExFit::WriteHistos","-------------------------------------------");
        WriteInfoStatus("TRExFit::WriteHistos","Writing histograms to file " + fileName + " ...");
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            std::shared_ptr<SampleHist> sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if(!sh){
                WriteDebugStatus("TRExFit::WriteHistos", "SampleHist[" + std::to_string(i_smp) + "] for sample " + fSamples[i_smp]->fName + " not there.");
                continue;
            }
            // set file and histo names for nominal
            sh->fHistoName = sh->fHist->GetName();
            sh->fFileName = fileName;
            // set file and histo names for systematics
            for(std::size_t i_syst=0;i_syst<sh->fSyst.size();i_syst++){
                if(sh->fSyst[i_syst]->fHistUp  ==nullptr) continue;
                if(sh->fSyst[i_syst]->fHistDown==nullptr) continue;
                sh->fSyst[i_syst]->fFileNameUp    = fileName;
                sh->fSyst[i_syst]->fHistoNameUp   = sh->fSyst[i_syst]->fHistUp->GetName();
                sh->fSyst[i_syst]->fFileNameDown  = fileName;
                sh->fSyst[i_syst]->fHistoNameDown = sh->fSyst[i_syst]->fHistDown->GetName();
                if(sh->fSyst[i_syst]->fIsShape){
                    sh->fSyst[i_syst]->fFileNameShapeUp    = fileName;
                    sh->fSyst[i_syst]->fHistoNameShapeUp   = sh->fSyst[i_syst]->fHistShapeUp->GetName();
                    sh->fSyst[i_syst]->fFileNameShapeDown  = fileName;
                    sh->fSyst[i_syst]->fHistoNameShapeDown = sh->fSyst[i_syst]->fHistShapeDown->GetName();
                }
            }
            if(singleOutputFile) sh->WriteToFile(fFiles[0]   ,reWriteOrig);
            else                 sh->WriteToFile(fFiles[i_ch],reWriteOrig);
        }
    }
    WriteInfoStatus("TRExFit::WriteHistos","-------------------------------------------");
}

//__________________________________________________________________________________
// Draw morphing plots
void TRExFit::DrawMorphingPlots(const std::string& name) const{
    for(auto reg : fRegions){
        TCanvas c("c","c",600,600);
        TPad p0("p0","p0",0,0.35,1,1);
        TPad p1("p1","p1",0,0,1,0.35);
        p0.SetBottomMargin(0);
        p1.SetTopMargin(0);
        p1.SetBottomMargin(0.3);
        p0.Draw();
        p1.Draw();
        p0.cd();
        int nTemp = 0;
        std::vector<std::unique_ptr<TH1> > hVec;
        std::vector<std::unique_ptr<TH1> > hVecRatio;
        for(auto& sh : reg->fSampleHists){
            Sample* smp = sh->fSample;
            // if the sample has morphing
            if(smp->fIsMorph[name]){
                std::unique_ptr<TH1> h(static_cast<TH1*>(sh->fHist->Clone(("h_temp_"+smp->fName).c_str())));
                if(h->GetFillColor()!=0) h->SetLineColor(h->GetFillColor());
                h->SetFillStyle(0);
                h->SetLineWidth(2);
                h->Scale(1./h->Integral());
                if(nTemp==0) h->Draw("HIST");
                else         h->Draw("HIST same");
                hVec.push_back(std::move(h));
                nTemp++;
            }
        }
        if(hVec.size()>0){
            hVec[0]->SetMaximum(1.25*hVec[0]->GetMaximum());
            hVec[0]->GetYaxis()->SetTitle("Fraction of events");
            hVec[0]->GetYaxis()->SetTitleOffset(1.75);
            // ratio
            p1.cd();
            for(const auto& hh : hVec){
                hVecRatio.push_back(std::move(std::unique_ptr<TH1>(static_cast<TH1*>(hh->Clone()))));
                hVecRatio[hVecRatio.size()-1]->Divide(hVec[0].get());
                if(hVecRatio.size()-1==0) hVecRatio[hVecRatio.size()-1]->Draw("HIST");
                else                      hVecRatio[hVecRatio.size()-1]->Draw("HIST same");
            }
            hVecRatio[0]->SetMinimum(0.91);
            hVecRatio[0]->SetMaximum(1.09);
            hVecRatio[0]->GetXaxis()->SetTitle(reg->fVariableTitle.c_str());
            hVecRatio[0]->GetYaxis()->SetTitle("Ratio");
            hVecRatio[0]->GetYaxis()->SetTitleOffset(1.75);
            hVecRatio[0]->GetXaxis()->SetTitleOffset(3);
            for(const auto& format : TRExFitter::IMAGEFORMAT) {
                c.SaveAs((fName+"/Morphing/Templates_"+name+"_"+reg->fName+"."+format).c_str());
            }
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots
void TRExFit::DrawSystPlots() const{
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
            std::shared_ptr<SampleHist> sh = fRegions[i_ch]->fSampleHists[i_smp];
            sh->DrawSystPlot("all");
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots for sombined samples
void TRExFit::DrawSystPlotsSumSamples() const{
    WriteInfoStatus("TRExFit::DrawSystPlotsSumSamples", "-------------------------------------------");
    WriteInfoStatus("TRExFit::DrawSystPlotsSumSamples", "Drawing combined plots of syst effects on data...");
    std::unique_ptr<TH1> h_dataCopy(nullptr);
    for(const auto& reg : fRegions){
        SampleHist hist{};
        bool empty = true;
        std::set<std::string> systNames;
        for(int i_regSmp=0; i_regSmp<reg->fNSamples; i_regSmp++){
            for(std::size_t i_smSyst=0; i_smSyst<reg->fSampleHists[i_regSmp]->fSyst.size(); i_smSyst++){
                systNames.insert(reg->fSampleHists[i_regSmp]->fSyst[i_smSyst]->fName);
            }
        }
        for(int i_smp=0;i_smp<reg->fNSamples;i_smp++){
            if(reg->fSampleHists[i_smp]->fSample->fType==Sample::DATA) h_dataCopy=std::unique_ptr<TH1>(static_cast<TH1*>(reg->fSampleHists[i_smp]->fHist->Clone()));
            else if(reg->fSampleHists[i_smp]->fSample->fType==Sample::GHOST) continue;
            else {
                double scale = Common::GetNominalMorphScale(reg->fSampleHists[i_smp].get());
                if(empty){
                    hist.CloneSampleHist(reg->fSampleHists[i_smp].get(),systNames, scale);
                    hist.fName = reg->fName + "_Combined";
                    empty=false;
                } else {
                    hist.SampleHistAdd(reg->fSampleHists[i_smp].get(), scale);
                }
            }
        }
        hist.DrawSystPlot("all", h_dataCopy.get(), true, fSystDataPlot_upFrame);
    }
}

//__________________________________________________________________________________
// this method takes care of rebinning, smoothing, fixing
void TRExFit::CorrectHistograms(){
    //
    // loop on regions, and then perform a set of operations for each of them
    for(auto reg : fRegions){
        //
        // 1. Reset histograms to the ones save as "_orig" (both for nominal and systematics
        for(auto smp : fSamples){
            //
            // eventually skip sample / region combination
            if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            TH1* h_orig = sh->fHist_orig.get();
            TH1* h = nullptr;
            if(h_orig!=nullptr) h = (TH1*)h_orig->Clone(sh->fHist->GetName());
            sh->fHist = std::unique_ptr<TH1>(h);
            if(sh->fHist==nullptr) continue;
            sh->fHist->SetLineColor(linecolor);
            sh->fHist->SetFillColor(fillcolor);

            // Scale MC stat
            Common::ScaleMCstatInHist(sh->fHist.get(), smp->fMCstatScale);

            // loop on systematics
            for(auto& syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                // skip also separate gamma systs
                if(syst->fName.find("stat_")!=std::string::npos) continue;
                //
                // get the original syst histograms & reset the syst histograms
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic( syst->fName );
                //
                if(syh==nullptr) continue;
                TH1* hUp_orig   = syh->fHistUp_orig.get();
                TH1* hDown_orig = syh->fHistDown_orig.get();
                //
                // if Overall only => fill SystematicHist
                if(syst->fType==Systematic::OVERALL){
                    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
                        if(hUp_orig!=nullptr)   hUp_orig  ->SetBinContent(i_bin,h_orig->GetBinContent(i_bin)*(1.+syst->fOverallUp));
                        if(hDown_orig!=nullptr) hDown_orig->SetBinContent(i_bin,h_orig->GetBinContent(i_bin)*(1.+syst->fOverallDown));
                    }
                }
                //
                TH1* hUp   = nullptr;
                TH1* hDown = nullptr;
                if(hUp_orig!=nullptr)   hUp   = static_cast<TH1*>(hUp_orig->Clone( syh->fHistUp->GetName()));
                if(hDown_orig!=nullptr) hDown = static_cast<TH1*>(hDown_orig->Clone(syh->fHistDown->GetName()));
                syh->fHistUp.reset(hUp);
                syh->fHistDown.reset(hDown);
            }
        }
        //
        // 2. Rebin
        for(auto smp : fSamples){
            //
            // eventually skip sample / region combination
            if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            //
            // Rebinning (FIXME: better to introduce a method Region::Rebin() ?)
            if(reg->fHistoNBinsRebinPost>0){
                WriteDebugStatus("TRExFit::CorrectHistograms", "rebinning " + smp->fName + " to " + std::to_string(reg->fHistoNBinsRebinPost) + " bins.");
                sh->fHist = std::unique_ptr<TH1>(sh->fHist->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]));
                for(auto& syh : sh->fSyst){
                    WriteDebugStatus("TRExFit::CorrectHistograms", "  systematic " + syh->fName + " to " + std::to_string(reg->fHistoNBinsRebinPost) + " bins.");
                    if(syh==nullptr) continue;
                    if(syh->fSystematic->fSampleUp==""   && syh->fSystematic->fHasUpVariation   && syh->fHistUp!=nullptr)   syh->fHistUp.reset(syh->fHistUp  ->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]));
                    else                                                                                                    syh->fHistUp.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistUp->GetName())));
                    if(syh->fSystematic->fSampleDown=="" && syh->fSystematic->fHasDownVariation && syh->fHistDown!=nullptr) syh->fHistDown.reset(syh->fHistDown->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]));
                    else                                                                                                    syh->fHistDown.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistDown->GetName())));
                }
                //
                // rebin also separate-gamma hists!
                if(smp->fSeparateGammas){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic( "stat_"+smp->fName );
                    if(syh==nullptr) continue;
                    if(syh->fHistUp!=nullptr)   syh->fHistUp  ->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]);
                    if(syh->fHistDown!=nullptr) syh->fHistDown->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]);
                }
            }
        }

        // Randomize MC (before add/multiply/scale)
        if(TRExFitter::OPTION["RandomizeMC"]!=0){
            gRandom->SetSeed(TRExFitter::OPTION["RandomizeMC"]);
            for(auto smp : fSamples){
                //
                // eventually skip sample / region combination
                if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                //
                std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                if(sh->fHist==nullptr) continue;
                //
                if(smp->fUseMCStat){
                    TH1* hTmp = sh->fHist.get();
                    for(int i_bin=1;i_bin<=hTmp->GetNbinsX();i_bin++){
                        hTmp->SetBinContent(i_bin,gRandom->Poisson( hTmp->GetBinContent(i_bin) ));
                    }
                    for(auto& syh : sh->fSyst){
                        for(int i_ud=0;i_ud<2;i_ud++){
                            if(i_ud==0) hTmp = syh->fHistUp.get();
                            else        hTmp = syh->fHistDown.get();
                            for(int i_bin=1;i_bin<=hTmp->GetNbinsX();i_bin++){
                                hTmp->SetBinContent(i_bin,gRandom->Poisson( hTmp->GetBinContent(i_bin) ));
                            }
                        }
                    }
                }
            }
        }

        // 3. Add/Multiply/Scale
        for(auto smp : fSamples){
            //
            // eventually skip sample / region combination
            if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            //
            // Subtraction / Addition of sample
            for(auto sample : smp->fSubtractSamples){
                WriteDebugStatus("TRExFit::CorrectHistograms"," subtracting sample " + sample + " from sample " + smp->fName);
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(sample);
                if(smph0!=nullptr) sh->Add(smph0.get(),-1);
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+sample+" not found ...");
            }
            for(auto sample : smp->fAddSamples){
                WriteDebugStatus("TRExFit::CorrectHistograms", "adding sample " + sample + " to sample " + smp->fName);
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(sample);
                if(smph0!=nullptr) sh->Add(smph0.get());
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+sample+" not found ...");
            }
            // Division & Multiplication by other samples
            if(smp->fMultiplyBy!=""){
                WriteDebugStatus("TRExFit::CorrectHistograms", "multiplying " + smp->fName  + " by sample " + smp->fMultiplyBy);
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(smp->fMultiplyBy);
                if(smph0!=nullptr) sh->Multiply(smph0.get());
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+smp->fMultiplyBy+" not found ...");
            }
            if(smp->fDivideBy!=""){
                WriteDebugStatus("TRExFit::CorrectHistograms", "dividing " + smp->fName  + " by sample " + smp->fDivideBy + " from sample " + smp->fName);
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(smp->fDivideBy);
                if(smph0!=nullptr) sh->Divide(smph0.get());
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+smp->fDivideBy+" not found ...");
            }
            // Norm to sample
            if(smp->fNormToSample!=""){
                WriteDebugStatus("TRExFit::CorrectHistograms", "normalizing " + smp->fName  + " to sample " + smp->fNormToSample);
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(smp->fNormToSample);
                if(smph0!=nullptr) sh->Scale(smph0->fHist->Integral()/sh->fHist->Integral());
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+smp->fNormToSample+" not found ...");
            }

            //
            // For SampleUp / SampleDown
            for(auto& syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                // get the original syst histograms & reset the syst histograms
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic( syst->fName );
                //
                // if syst defined with SampleUp / SampleDown
                if( syst->fSampleUp != "" || syst->fSampleDown != "" ){
                    WriteDebugStatus("TRExFit::CorrectHistograms", "SampleUp/SampleDown set for systematic " + syst->fName + ".");
                    bool isDummy = ( syst->fDummyForSamples.size()>0 && Common::FindInStringVector(syst->fDummyForSamples,smp->fName)>=0 );
                    std::unique_ptr<TH1> h_up = nullptr;
                    if(syst->fSampleUp   !="" && !isDummy){
                        if(reg->GetSampleHist(syst->fSampleUp  )){
                            h_up.reset(static_cast<TH1*>(reg->GetSampleHist(syst->fSampleUp  )->fHist->Clone("h_tmp_up")));
                        }
                    }
                    else{
                        h_up.reset(static_cast<TH1*>(sh->fHist->Clone("h_tmp_up")));
                    }
                    std::unique_ptr<TH1> h_down = nullptr;
                    if(syst->fSampleDown !="" && !isDummy){
                        if(reg->GetSampleHist(syst->fSampleDown)){
                            h_down.reset(static_cast<TH1*>(reg->GetSampleHist(syst->fSampleDown)->fHist.get()->Clone("h_tmp_down")));
                        }
                    }
                    else{
                        h_down.reset(static_cast<TH1*>(sh->fHist.get()->Clone("h_tmp_down")));
                    }
                    //
                    // if systematic also uses ReferenceSample, produce syst variations according to the refefence sample instead of nominal
                    if(syst->fReferenceSample!=""){
                        WriteDebugStatus("TRExFit::CorrectHistograms", "ReferenceSample set for a systematic with SampleUp/SampleDown. Building proper systematic variation.");
                        std::shared_ptr<SampleHist> refSh = reg->GetSampleHist(syst->fReferenceSample);
                        if(refSh!=nullptr){
                            if(syst->fSampleUp != ""){
                                h_up->Divide(refSh->fHist.get());
                                h_up->Multiply(sh->fHist.get());
                            }
                            if(syst->fSampleDown != ""){
                                h_down->Divide(refSh->fHist.get());
                                h_down->Multiply(sh->fHist.get());
                            }
                        }
                        else{
                            WriteWarningStatus("TRExFit::CorrectHistograms","ReferenceSample for systematc " + syst->fName + " set but no corresponding sample found. Ignoring.");
                        }
                    }
                    syh = sh->AddHistoSyst(syst->fName,syst->fStoredName,h_up.get(),h_down.get());
                    syh->fSystematic = syst;
                }
            }

            //
            // Save to _preSmooth histograms (to be shown in syst plots) at this point
            sh->fHist_preSmooth.reset(static_cast<TH1*>(sh->fHist->Clone(Form("%s_preSmooth",sh->fHist->GetName()))));
            sh->fHist_preSmooth->SetDirectory(nullptr);
            for(auto& syh : sh->fSyst){
                if(syh!=nullptr){
                    if(syh->fHistUp!=nullptr)   syh->fHistUp_preSmooth.reset(static_cast<TH1*>(syh->fHistUp->Clone(Form("%s_preSmooth",syh->fHistUp->GetName()))));
                    else                        syh->fHistUp_preSmooth.reset(static_cast<TH1*>(sh->fHist_preSmooth->Clone()));
                    syh->fHistUp_preSmooth->SetDirectory(nullptr);
                    if(syh->fHistDown!=nullptr) syh->fHistDown_preSmooth.reset(static_cast<TH1*>(syh->fHistDown->Clone(Form("%s_preSmooth",syh->fHistDown->GetName()))));
                    else                        syh->fHistDown_preSmooth.reset(static_cast<TH1*>(sh->fHist_preSmooth->Clone()));
                    syh->fHistDown_preSmooth->SetDirectory(nullptr);
                }
            }

            //
            // Fix empty bins
            if(smp->fType!=Sample::DATA && smp->fType!=Sample::SIGNAL){
                sh->FixEmptyBins(fSuppressNegativeBinWarnings);
            }

            //
            // Eventually smooth nominal histogram  (use with caution...)
            TH1* h_correction = nullptr;
            bool isFlat = false;
            if(smp->fSmooth && !reg->fSkipSmoothing){
                h_correction = (TH1*)sh->fHist->Clone( Form("%s_corr",sh->fHist->GetName()) );
                TH1* h0 = (TH1*)sh->fHist->Clone( Form("%s_orig0",sh->fHist->GetName()) );
                if (fSmoothOption == HistoTools::TTBARRESONANCE) {
                    isFlat = false;
                    Common::SmoothHistogramTtres( sh->fHist.get() );
                } else {
                    isFlat = Common::SmoothHistogram( sh->fHist.get() );
                }
                h_correction->Divide( h0 );
            }

            //
            // Systematics
            for(auto& syst : smp->fSystematics){
                if(syst==nullptr) continue;
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic( syst->fName );
                if(syh==nullptr) continue;
                TH1* hUp   = syh->fHistUp.get();
                TH1* hDown = syh->fHistDown.get();
                //
                // if Overall only, re-create it if smoothing was applied
                if(syst->fType==Systematic::OVERALL){
                    if(h_correction!=nullptr && smp->fSmooth){
                        for(int i_bin=1;i_bin<=sh->fHist->GetNbinsX();i_bin++){
                            hUp  ->SetBinContent(i_bin,sh->fHist->GetBinContent(i_bin)*(1.+syst->fOverallUp));
                            hDown->SetBinContent(i_bin,sh->fHist->GetBinContent(i_bin)*(1.+syst->fOverallDown));
                        }
                    }
                    continue;
                }
                //
                // correct according to the sample nominal smoothing
                if(h_correction!=nullptr && smp->fSmooth){
                    if(hUp!=nullptr  ) Common::SmoothHistogram( hUp  , isFlat );
                    if(hDown!=nullptr) Common::SmoothHistogram( hDown, isFlat );
                }

                //
                // Histogram smoothing, Symmetrisation, Massaging...
                if(!reg->fSkipSmoothing) syh -> fSmoothType = syst -> fSmoothType;
                else                                syh -> fSmoothType = 0;
                syh -> fSymmetrisationType = syst -> fSymmetrisationType;

            }  // end syst loop
            //
            // Histograms checking
            for(auto& syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                if( sh->GetSystematic( syst->fName )==nullptr ) continue;
                //
                HistoTools::CheckHistograms( reg->GetSampleHist(smp->fName)->fHist.get() /*nominal*/,
                                            sh->GetSystematic( syst->fName ).get() /*systematic*/,
                                            smp->fType!=Sample::SIGNAL/*check bins with content=0*/,
                                            TRExFitter::HISTOCHECKCRASH /*cause crash if problem*/);
            }

            // set the fill color
            sh->fHist->SetFillColor(fillcolor);
            sh->fHist->SetLineColor(linecolor);
        } // end sample loop
        //
    } // end region loop

    //
    // Morph smoothing
    if(fSmoothMorphingTemplates!=""){
        for(auto par : fMorphParams){
            WriteInfoStatus("TRExFit::CorrectHistograms","Smoothing morphing templates for parameter "+par);
            if(fSmoothMorphingTemplates=="TRUE") SmoothMorphTemplates(par);
            else SmoothMorphTemplates(par,fSmoothMorphingTemplates);
            // to add: possibility to set initial values of parameters
        }
    }

    //
    // Plot Morphing templates
    if (fMorphParams.size() > 0){
        gSystem->mkdir((fName+"/Morphing").c_str());
        for(const auto& par : fMorphParams){
            DrawMorphingPlots(par);
        }
    }

    // drop normalisation part of systematic according to fDropNormIn
    for(auto reg : fRegions){
        for(auto& sh : reg->fSampleHists){
            if(sh->fHist==nullptr) continue;
            for(auto syst : fSystematics){
                if(  Common::FindInStringVector(syst->fDropNormIn, reg->fName)>=0
                  || Common::FindInStringVector(syst->fDropNormIn, sh->fSample->fName)>=0
                  || Common::FindInStringVector(syst->fDropNormIn, "all")>=0
                  ){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    if(sh->fHist->Integral()!=0){
                        WriteDebugStatus("TRExFit::CorrectHistograms", "  Normalising syst " + syst->fName + " for sample " + sh->fSample->fName);
                        Common::DropNorm(syh->fHistUp.get(), syh->fHistDown.get(), sh->fHist.get());
                    }
                }
            }
        }
    }

    // drop shape part of systematic according to fDropShapeIn
    for(auto reg : fRegions){
        for(auto& sh : reg->fSampleHists){
            if(sh->fHist==nullptr) continue;
            for(auto syst : fSystematics){
                if(  Common::FindInStringVector(syst->fDropShapeIn, reg->fName)>=0
                  || Common::FindInStringVector(syst->fDropShapeIn, sh->fSample->fName)>=0
                  || Common::FindInStringVector(syst->fDropShapeIn, "all")>=0
                  ){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    WriteDebugStatus("TRExFit::CorrectHistograms", "  Removing shape component of syst " + syst->fName + " for sample " + sh->fSample->fName);
                    Common::DropShape(syh->fHistUp.get(), syh->fHistDown.get(), sh->fHist.get());
                }
            }
        }
    }

    //
    // Smooth systematics
    SmoothSystematics("all");

    //
    // Artifificially set all systematics not to affect overall normalisation for sample or set of samples
    // (the form should be KeepNormForSamples: ttlight+ttc+ttb,wjets
    //
    for(auto reg : fRegions){
        for(auto syst : fSystematics){
            if(syst->fKeepNormForSamples.size()==0) continue;
            for(unsigned int ii=0;ii<syst->fKeepNormForSamples.size();ii++){
                std::vector<std::string> subSamples = Vectorize(syst->fKeepNormForSamples[ii],'+');
                // get nominal yield and syst yields for this sum of samples
                double yieldNominal = 0.;
                double yieldUp = 0.;
                double yieldDown = 0.;
                for(auto smp : fSamples){
                    if(Common::FindInStringVector(subSamples,smp->fName)<0) continue;
                    std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                    if(sh==nullptr) continue;
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    yieldNominal += sh ->fHist    ->Integral();
                    yieldUp      += syh->fHistUp  ->Integral();
                    yieldDown    += syh->fHistDown->Integral();
                }
                // scale each syst variation
                for(auto smp : fSamples){
                    if(Common::FindInStringVector(subSamples,smp->fName)<0) continue;
                    std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                    if(sh==nullptr) continue;
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    WriteDebugStatus("TRExFit::CorrectHistograms", "  Normalising syst " + syst->fName + " for sample " + smp->fName);
                    WriteDebugStatus("TRExFit::CorrectHistograms", "scaling by " + std::to_string(yieldNominal/yieldUp) + " (up), " + std::to_string(yieldNominal/yieldDown) + " (down)");
                    syh->fHistUp  ->Scale(yieldNominal/yieldUp);
                    syh->fHistDown->Scale(yieldNominal/yieldDown);
                }
            }
        }
    }

    // Systematics for morphing samples inherited from nominal sample
    if (fPropagateSystsForMorphing){
        for(auto par : fMorphParams){
            for(auto reg : fRegions){
                // find nominal morphing sample Hist
                double nominalValue = 0.;
                for(auto norm : fNormFactors){
                    if(norm->fName==par) nominalValue = norm->fNominal;
                }
                std::shared_ptr<SampleHist> shNominal = nullptr;
                for(auto& sh : reg->fSampleHists){
                    if(!sh->fSample->fIsMorph[par]) continue;
                    if(sh->fSample->fMorphValue[par]==nominalValue){ // FIXME: eventually add something to flag a sample as nominal for morphing
                        shNominal = sh;
                        break;
                    }
                }
                // loop on all other samples
                for(auto& sh : reg->fSampleHists){
                    if(!sh->fSample->fIsMorph[par]) continue;
                    if(sh != shNominal){
                        for(auto& syh : shNominal->fSyst){
                            std::shared_ptr<Systematic> syst = syh->fSystematic;
                            if(syst->fIsNormOnly){
                                std::shared_ptr<SystematicHist> syhNew = sh->AddOverallSyst(syst->fName,syst->fStoredName,syst->fOverallUp,syst->fOverallDown);
                                syhNew->fSystematic = syst;
                            }
                            else{
                                TH1* hUpNew   = (TH1*)syh->fHistUp->Clone();
                                TH1* hDownNew = (TH1*)syh->fHistDown->Clone();
                                hUpNew->Divide(shNominal->fHist.get());
                                hUpNew->Multiply(sh->fHist.get());
                                hDownNew->Divide(shNominal->fHist.get());
                                hDownNew->Multiply(sh->fHist.get());
                                std::shared_ptr<SystematicHist> syhNew = sh->AddHistoSyst(syst->fName,syst->fStoredName,hUpNew,hDownNew);
                                syhNew->fSystematic = syst;
                                sh->fSample->fUseSystematics = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Propagate all systematics from another sample
    for(auto reg : fRegions){
        for(auto smp : fSamples){
            if(smp->fSystFromSample != ""){
                // eventually skip sample / region combination
                if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                sh->fSample->fUseSystematics = true;
                //
                std::shared_ptr<SampleHist> shReference = reg->GetSampleHist(smp->fSystFromSample);
                for(auto& syh : shReference->fSyst){
                    std::shared_ptr<Systematic> syst = syh->fSystematic;
                    if(syst->fIsNormOnly){
                        std::shared_ptr<SystematicHist> syhNew = sh->AddOverallSyst(syst->fName,syst->fStoredName,syst->fOverallUp,syst->fOverallDown);
                        syhNew->fSystematic = syst;
                    }
                    else{
                        TH1* hUpNew   = (TH1*)syh->fHistUp->Clone();
                        TH1* hDownNew = (TH1*)syh->fHistDown->Clone();
                        hUpNew->Divide(shReference->fHist.get());
                        hUpNew->Multiply(sh->fHist.get());
                        hDownNew->Divide(shReference->fHist.get());
                        hDownNew->Multiply(sh->fHist.get());
                        std::shared_ptr<SystematicHist> syhNew = sh->AddHistoSyst(syst->fName,syst->fStoredName,hUpNew,hDownNew);
                        syhNew->fSystematic = syst;
                    }
                }
            }
        }
    }

    //
    // set the hasData flag
    bool hasData = false;
    for(const auto& smp : fSamples){
        if(smp->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    //
    // Poissonize data
    if(hasData && TRExFitter::OPTION["PoissonizeData"]!=0){
        for(auto reg : fRegions){
            if(reg->fData!=nullptr){
                if(reg->fData->fHist!=nullptr){
                    TH1 *hdata = reg->fData->fHist.get();
                    if(TRExFitter::OPTION["PoissonizeData"]>0) gRandom->SetSeed(TRExFitter::OPTION["PoissonizeData"]);
                    for(int i_bin=1;i_bin<=hdata->GetNbinsX();i_bin++){
                        hdata->SetBinContent(i_bin,gRandom->Poisson( hdata->GetBinContent(i_bin) ));
                        hdata->SetBinError(i_bin,std::sqrt(hdata->GetBinContent(i_bin)));
                    }
                }
            }
        }
    }

    // Manually change the shape of systematics
    RunForceShape();

    DropBins();
}

//__________________________________________________________________________________
//
void TRExFit::CloseInputFiles(){
    //
    // Close all input files
    for(auto& it : TRExFitter::TFILEMAP){
        TDirectory *dir = gDirectory;
        std::shared_ptr<TFile> f = it.second;
        if(f) {
            dir->cd();
            f->Close();
        }
    }
    TRExFitter::TFILEMAP.clear();
}

//__________________________________________________________________________________
//
void TRExFit::DrawAndSaveAll(std::string opt){
    bool isPostFit = opt.find("post")!=std::string::npos;
    //
    // Scale sample(s) to data (only pre-fit)
    if(!isPostFit && fScaleSamplesToData.size()>0){
        for(auto& reg : fRegions){
            std::unique_ptr<TH1> hTot = reg->GetTotHist(true);
            const double totPred = hTot->Integral();
            const double totData = reg->fData->fHist->Integral();
            double totToScale = 0;
            std::vector<std::shared_ptr<SampleHist> > shToScale;
            for(auto& sh : reg->fSampleHists){
                if(sh->fHist==nullptr) continue;
                if(sh->fSample->fType==Sample::GHOST){
                    WriteWarningStatus("TRExFit::CorrectHistograms","Requested to scale to data a GHOST sample, " + sh->fSample->fName + ". Skipping this sample.");
                    continue;
                }
                if(Common::FindInStringVector(fScaleSamplesToData,sh->fSample->fName)>=0){
                    shToScale.emplace_back(sh);
                    const double morph_scale = Common::GetNominalMorphScale(sh.get());
                    totToScale += morph_scale*sh->fHist->Integral();
                }
            }
            if(totToScale<=0 || shToScale.size()==0) continue;
            const double scale = (totData-(totPred-totToScale))/totToScale;
            for(auto& sh : shToScale){
                WriteInfoStatus("TRExFit::CorrectHistograms","Scaling sample " + sh->fSample->fName + " by " + std::to_string(scale) + " in region " + reg->fName);
                sh->Scale(scale);
            }

        }
    }
    else if(fScaleSamplesToData.size()>0){
        for(auto& reg : fRegions){
            reg->fScaleSamplesToData = fScaleSamplesToData;
        }
    }
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    if(TRExFitter::POISSONIZE) opt += " poissonize";
    if(isPostFit){
        if(fFitResultsFile!=""){
            ReadFitResults(fFitResultsFile);
        }
        else {
            ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
        }
    }
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        std::shared_ptr<TRExPlot> p(nullptr);
        fRegions[i_ch]->fUseStatErr = fUseStatErr;
        fRegions[i_ch]->fATLASlabel = fAtlasLabel;
        //
        if(fCustomAsimov!=""){
            std::string name = "customAsimov_"+fCustomAsimov;
            std::shared_ptr<SampleHist> cash = fRegions[i_ch]->GetSampleHist(name);
            if(cash==nullptr){
                WriteWarningStatus("TRExFit::DrawAndSaveAll", "No Custom Asimov " + fCustomAsimov + " available. Taking regular Asimov.");
            }
            else{
                std::string s = cash->fHist->GetName();
                WriteDebugStatus("TRExFit::DrawAndSaveAll", "  Adding Custom-Asimov Data: " + s);
                fRegions[i_ch]->fData = cash;
            }
        }
        for (auto& ireg : fRegions) {
            ireg->fFolder = fName;
            ireg->fHEPDataFormat = fHEPDataFormat;
        }
        //
        if(isPostFit){
            std::ofstream pullTex;
            if(fWithPullTables){
                gSystem->mkdir((fName+"/Tables").c_str()); // need to create directory, as it may not exist yet
                pullTex.open((fName+"/Tables/Pulls_"+fSuffix+fRegions[i_ch]->fName+".tex").c_str());
                pullTex << "\\documentclass[10pt]{article}" << std::endl;
                pullTex << "\\usepackage{siunitx}" << std::endl;
                pullTex << "\\usepackage{xcolor}" << std::endl;
                pullTex << "\\usepackage[margin=0.1in,landscape,papersize={210mm,100mm}]{geometry}" << std::endl;
                pullTex << "\\begin{document}" << std::endl;

                pullTex << "\\begin{tabular}{|lr|}\n" << std::endl;
                pullTex << "\\hline\\hline\n" << std::endl;
                TString region(fRegions[i_ch]->fName);
                pullTex << "\\multicolumn{2}{|c|}{" << fRegions[i_ch]->fTexLabel << "} \\\\\n"<< std::endl;
            }

            gSystem->mkdir( (fName + "/Histograms/").c_str() );
            if(fRegions[i_ch]->fRegionDataType==Region::ASIMOVDATA) p = fRegions[i_ch]->DrawPostFit(fFitResults,pullTex,fMorphParams,fPrePostFitCanvasSize,opt+" blind");
            else                                                    p = fRegions[i_ch]->DrawPostFit(fFitResults,pullTex,fMorphParams,fPrePostFitCanvasSize,opt);
            for(const auto& format : TRExFitter::IMAGEFORMAT) {
                p->SaveAs((fName+"/Plots/"+fRegions[i_ch]->fName+"_postFit"+fSuffix+"."+format).c_str());
            }

            if(fWithPullTables){
                pullTex << "\\hline\\hline\n" << std::endl;
                pullTex << "\\end{tabular}"  << std::endl;
                pullTex << "\\end{document}" << std::endl;
                pullTex.close();
            }
        }
        else{
            if(fRegions[i_ch]->fRegionDataType==Region::ASIMOVDATA) p = fRegions[i_ch]->DrawPreFit(fPrePostFitCanvasSize, opt+" blind");
            else                                                    p = fRegions[i_ch]->DrawPreFit(fPrePostFitCanvasSize, opt);
            // this line to fix the y-axis maximum getting doubled in some cases (FIXME)
            if((fRegions[i_ch]->fYmin==0) && (fRegions[i_ch]->fYmax==0) && (fRegions[i_ch]->fYmaxScale==0)){
                if(!fRegions[i_ch]->fLogScale) p->h_dummy->GetYaxis()->SetRangeUser(p->h_dummy->GetYaxis()->GetXmin(),p->h_dummy->GetMaximum());
                else                           p->h_dummy->GetYaxis()->SetRangeUser(1                                ,p->h_dummy->GetMaximum());
            }
            for(const auto& format : TRExFitter::IMAGEFORMAT) {
                p->SaveAs((fName+"/Plots/"+fRegions[i_ch]->fName+fSuffix+"."+format).c_str());
            }
        }
    }
}

//__________________________________________________________________________________
//
std::shared_ptr<TRExPlot> TRExFit::DrawSummary(std::string opt, std::shared_ptr<TRExPlot> prefit_plot) {
    WriteInfoStatus("TRExFit::DrawSummary", "-------------------------------------------");
    WriteInfoStatus("TRExFit::DrawSummary", "Building Summary Plot...");
    gSystem->mkdir(fName.c_str(),true);
    const bool isPostFit = opt.find("post")!=std::string::npos;
    const bool checkVR = opt.find("valid")!=std::string::npos;
    if(TRExFitter::POISSONIZE) opt += " poissonize";
    // build one bin per region
    std::unique_ptr<TH1D> h_data = nullptr;
    std::vector<std::unique_ptr<TH1D> > h_sig;
    std::vector<std::unique_ptr<TH1D> > h_bkg;
    std::unique_ptr<TGraphAsymmErrors> g_err(nullptr);
    int Nsig = 0;
    int Nbkg = 0;
    //
    std::string name;
    std::string title;
    int lineColor;
    int fillColor;
    int lineWidth;
    double integral;
    double intErr; // to store the integral error
    //
    // Building region - bin correspondence
    //
    std::vector<int> regionVec;
    std::vector<int> divisionVec;
    //
    if(checkVR){
        if(fSummaryPlotValidationRegions.size()==0){
            for(int i_ch=0;i_ch<fNRegions;i_ch++){
                if(fRegions[i_ch]->fRegionType==Region::VALIDATION){
                    regionVec.push_back(i_ch);
                }
            }
        }
        else{
            for(int i_reg=0;i_reg<(int)fSummaryPlotValidationRegions.size();i_reg++){
                if(fSummaryPlotValidationRegions[i_reg]=="|"){
                    divisionVec.push_back(regionVec.size());
                    continue;
                }
                for(int i_ch=0;i_ch<fNRegions;i_ch++){
                    if(fSummaryPlotValidationRegions[i_reg]==fRegions[i_ch]->fName){
                        regionVec.push_back(i_ch);
                        break;
                    }
                }
            }
        }
    }
    else{
        if(fSummaryPlotRegions.size()==0){
            for(int i_ch=0;i_ch<fNRegions;i_ch++){
                if(fRegions[i_ch]->fRegionType!=Region::VALIDATION){
                    regionVec.push_back(i_ch);
                }
            }
        }
        else{
            for(int i_reg=0;i_reg<(int)fSummaryPlotRegions.size();i_reg++){
                if(fSummaryPlotRegions[i_reg]=="|"){
                    divisionVec.push_back(regionVec.size());
                    continue;
                }
                for(int i_ch=0;i_ch<fNRegions;i_ch++){
                    if(fSummaryPlotRegions[i_reg]==fRegions[i_ch]->fName){
                        regionVec.push_back(i_ch);
                        break;
                    }
                }
            }
        }
    }
    //
    if(regionVec.size()==0) return 0;
    //
    int Nbin = (int)regionVec.size();
    if(Nbin<=0) return nullptr;
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        std::shared_ptr<SampleHist> sh = nullptr;
        name = (fSamples[i_smp]->fName).c_str();
        title = fSamples[i_smp]->fTitle.c_str();
        if(fSamples[i_smp]->fGroup != "") title = fSamples[i_smp]->fGroup.c_str();
        // look for the first SampleHist defined for this sample
        for(int i_ch=0;i_ch<(int)regionVec.size();i_ch++){
            sh = fRegions[regionVec[i_ch]]->GetSampleHist( name );
            if(sh!=nullptr) break;
        }
        // skip sample if no SampleHist found
        if(sh==nullptr) continue;
        if(sh->fHist==nullptr) continue;
        //
        lineColor = sh->fHist->GetLineColor();
        fillColor = sh->fHist->GetFillColor();
        lineWidth = sh->fHist->GetLineWidth();
        //
        if(fSamples[i_smp]->fType==Sample::SIGNAL){
            h_sig.emplace_back(new TH1D(name.c_str(),title.c_str(), Nbin,0,Nbin));
            std::string temp_string = h_sig[Nsig]->GetTitle();
            WriteDebugStatus("TRExFit::DrawSummary", "Adding Signal: " + temp_string);
            h_sig.back()->SetLineColor(lineColor);
            h_sig.back()->SetFillColor(fillColor);
            h_sig.back()->SetLineWidth(lineWidth);
            for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
                std::unique_ptr<TH1D> h;
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=nullptr){
                    if(isPostFit)  h.reset(static_cast<TH1D*>(sh->fHist_postFit->Clone())); // Michele
                    else           h.reset(static_cast<TH1D*>(sh->fHist->Clone())); // Michele
                    //
                    if(!isPostFit){
                        // FIXME SF
                        // scale it according to NormFactors
                        for(unsigned int i_nf=0;i_nf<sh->fSample->fNormFactors.size();i_nf++){
                            h->Scale(sh->fSample->fNormFactors[i_nf]->fNominal);
                            WriteDebugStatus("TRExFit::DrawSummary", "Scaling " + sh->fSample->fName + " by " + std::to_string(sh->fSample->fNormFactors[i_nf]->fNominal));
                        }
                    }
                    //
                    integral = h->IntegralAndError(1,h->GetNbinsX(),intErr);
                    // this becuase MC stat is taken into account by the gammas
                    if( (isPostFit && fUseGammaPulls) || !fUseStatErr || (!sh->fSample->fUseMCStat && !sh->fSample->fSeparateGammas))
                        intErr = 0.;
                }
                else{
                    integral = 0.;
                    intErr   = 0.;
                }
                h_sig.back()->SetBinContent( i_bin,integral );
                h_sig.back()->SetBinError( i_bin,intErr );
            }
            Nsig++;
        }
        else if(fSamples[i_smp]->fType==Sample::BACKGROUND){
            h_bkg.emplace_back(new TH1D(name.c_str(),title.c_str(), Nbin,0,Nbin));
            std::string temp_string = h_bkg[Nbkg]->GetTitle();
            WriteDebugStatus("TRExFit::DrawSummary", "Adding Bkg:    " + temp_string);
            h_bkg.back()->SetLineColor(lineColor);
            h_bkg.back()->SetFillColor(fillColor);
            h_bkg.back()->SetLineWidth(lineWidth);
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                std::unique_ptr<TH1D> h;
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=nullptr){
                    if(isPostFit)  h.reset(static_cast<TH1D*>(sh->fHist_postFit->Clone())); // Michele
                    else           h.reset(static_cast<TH1D*>(sh->fHist->Clone())); // Michele
                    //
                    if(!isPostFit){
                        // scale it according to NormFactors but only for the correct regions
                        for(unsigned int i_nf=0; i_nf<sh->fSample->fNormFactors.size(); ++i_nf){
                            const std::vector<std::string>& regions = sh->fSample->fNormFactors[i_nf]->fRegions;
                            if (std::find(regions.begin(), regions.end(), fRegions[regionVec[i_bin-1]]->fName) != regions.end() ) {
                                h->Scale(sh->fSample->fNormFactors[i_nf]->fNominal);
                                WriteDebugStatus("TRExFit::DrawSummary", "Scaling " + sh->fSample->fName + " by " + std::to_string(sh->fSample->fNormFactors[i_nf]->fNominal));
                            }
                        }
                    }
                    //
                    integral = h->IntegralAndError(1,h->GetNbinsX(),intErr);
                    //
                    // this becuase MC stat is taken into account by the gammas
                    if( (isPostFit && fUseGammaPulls) || !fUseStatErr || (!sh->fSample->fUseMCStat && !sh->fSample->fSeparateGammas))
                        intErr = 0.;
                }
                else{
                    integral = 0.;
                    intErr = 0.;
                }
                h_bkg.back()->SetBinContent( i_bin,integral );
                h_bkg.back()->SetBinError( i_bin,intErr );
            }
            Nbkg++;
        }
        else if(fSamples[i_smp]->fType==Sample::DATA){
            h_data = std::make_unique<TH1D>(name.c_str(),title.c_str(), Nbin,0,Nbin);
            const std::string temp_string = h_data->GetTitle();
            WriteDebugStatus("TRExFit::DrawSummary", "Adding Data:   " + temp_string);
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                if(fRegions[regionVec[i_bin-1]]->fRegionDataType==Region::ASIMOVDATA)
                    h_data->SetBinContent( i_bin,0 );
                else
                    h_data->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fData->fHist->Integral() );
            }
        }
    }
    //
    std::shared_ptr<TRExPlot> p;
    //
    if (fSummaryCanvasSize.size() == 0){
        p = std::make_shared<TRExPlot>(fInputName+"_summary",900,700,TRExFitter::NORATIO);
    } else {
        p = std::make_shared<TRExPlot>(fInputName+"_summary",fSummaryCanvasSize.at(0),fSummaryCanvasSize.at(1),TRExFitter::NORATIO);
    }
    if(fYmin!=0) p->fYmin = fYmin;
    else         p->fYmin = 1;
    if(fYmax!=0) p->fYmax = fYmax;
    else         p->SetYmaxScale(2);
    p->SetXaxis("",false);
    p->AddLabel(fLabel);
    if(TRExFitter::OPTION["NoPrePostFit"]==0){
        if(isPostFit) p->AddLabel("Post-Fit");
        else          p->AddLabel("Pre-Fit");
    }
    //
    if(isPostFit) p->fRatioYmax = fRatioYmaxPostFit;
    else          p->fRatioYmax = fRatioYmax;
    if(isPostFit) p->fRatioYmin = fRatioYminPostFit;
    else          p->fRatioYmin = fRatioYmin;
    //
    // propagate settings from Job to plot
    p->fRatioYtitle = fRatioYtitle;
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType=="DATA/MC"){
        p->fRatioType = "DATA/BKG";
    }
    p->fATLASlabel = fAtlasLabel;
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    p->fLabelX = fLabelXSummary;
    p->fLabelY = fLabelYSummary;
    p->fLegendX1 = fLegendX1Summary;
    p->fLegendX2 = fLegendX2Summary;
    p->fLegendY = fLegendYSummary;
    p->fLegendNColumns = fLegendNColumnsSummary;
    if(fBlindingThreshold >= 0) {
        std::unique_ptr<TH1> signal(nullptr);
        if (Nsig > 0) {
            signal.reset(static_cast<TH1*>(h_sig[0]->Clone()));
        }
        std::unique_ptr<TH1> bkg(nullptr);
        for (int i = 0; i < Nbkg; ++i) {
            if (!h_bkg[i]) continue;
            if (!bkg) bkg.reset(static_cast<TH1*>(h_bkg[i]->Clone()));
            else      bkg->Add(h_bkg[i].get());
        }
        const std::vector<int>& blindedBins = Common::ComputeBlindedBins(signal.get(),
                                                                         bkg.get(),
                                                                         fBlindingType,
                                                                         fBlindingThreshold);
        p->SetBinBlinding(blindedBins);
        if(isPostFit && fKeepPrefitBlindedBins && fBlindedBins) p->SetBlindingHisto(fBlindedBins.get());
    }
    //
    if(h_data) p->SetData(h_data.get(), h_data->GetTitle());
    for(int i=0;i<Nsig;i++){
        if(TRExFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    h_sig[i].get(),h_sig[i]->GetTitle());
        if(TRExFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(h_sig[i].get(),h_sig[i]->GetTitle());
        if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(h_sig[i].get(),h_sig[i]->GetTitle());
    }
    for(int i=0;i<Nbkg;i++){
        p->AddBackground(h_bkg[i].get(),h_bkg[i]->GetTitle());
    }

    if( TRExFitter::PREFITONPOSTFIT && isPostFit) {
      p->h_tot_bkg_prefit = (TH1*)prefit_plot->GetTotBkg()->Clone("h_tot_bkg_prefit");
    }

    //
    // Build tot
    //
    std::unique_ptr<TH1D> h_tot = std::make_unique<TH1D>("h_Tot_summary","h_Tot_summary", Nbin,0,Nbin);

    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        double mc_stat_err;
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot_postFit->IntegralAndError(1,fRegions[regionVec[i_bin-1]]->fTot_postFit->GetNbinsX(),mc_stat_err) );
        else          h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot->IntegralAndError(1,fRegions[regionVec[i_bin-1]]->fTot->GetNbinsX(),mc_stat_err) );
        if(!fUseStatErr || (isPostFit && fUseGammaPulls)) h_tot->SetBinError( i_bin,0. );
        else                                              h_tot->SetBinError( i_bin,mc_stat_err );
    }
    //
    //   Build error band
    // build the vectors of variations
    std::vector< std::shared_ptr<TH1> > h_up;
    std::vector< std::shared_ptr<TH1> > h_down;
    std::unique_ptr<TH1> h_tmp_Up;
    std::unique_ptr<TH1> h_tmp_Down;
    std::vector<std::string> systNames;
    std::vector<std::string> npNames;
    int i_np = -1;
    // actual systematics
    for(const auto& isyst : fSystematics) {
        if(isPostFit && isyst->fType == Systematic::SHAPE) continue;
        std::string systName = isyst->fName;
        std::string systNuisPar = systName;
        systNames.push_back( systName );
        if(isyst!=nullptr)
            systNuisPar = isyst->fNuisanceParameter;
        if(Common::FindInStringVector(npNames,systNuisPar)<0){
            npNames.push_back(systNuisPar);
            i_np++;
        }
        else
            continue;
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
            // find the systematic in the region
            int syst_idx = -1;
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(systNuisPar==TRExFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
                    syst_idx = j_syst;
                }
            }
            //
            if(isPostFit){
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp_postFit[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown_postFit[syst_idx]->Clone()));
                }
            }
            else{
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown[syst_idx]->Clone()));
                }
            }
            if(i_bin==1){
                h_up.  emplace_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,systName.c_str()), Form("h_Tot_%s_Up_TMP",  systName.c_str()), Nbin,0,Nbin) );
                h_down.emplace_back( new TH1D(Form("h_Tot_%s_Down_TMP",systName.c_str()), Form("h_Tot_%s_Down_TMP",systName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
            //
            // look for other syst with the same np
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(j_syst==syst_idx) continue;
                if(systNuisPar==TRExFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
                    std::unique_ptr<TH1> h_tmp = nullptr;
                    if(isPostFit){
                        h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp_postFit[j_syst]->Clone()));
                        h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown_postFit[j_syst]->Clone()));
                        h_tmp.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                    }
                    else{
                        h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp[j_syst]->Clone()));
                        h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown[j_syst]->Clone()));
                        h_tmp.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                    }
                    h_up[i_np]  ->AddBinContent( i_bin,h_tmp_Up  ->Integral()-h_tmp->Integral() );
                    h_down[i_np]->AddBinContent( i_bin,h_tmp_Down->Integral()-h_tmp->Integral() );
                }
            }
        }
    }
    // add the gammas (only if post-fit)
    if(isPostFit && fUseGammaPulls){
        // loop on regions
        for(int i_ch=1;i_ch<=Nbin;i_ch++){
            Region *region = fRegions[regionVec[i_ch-1]];
            if(region==nullptr) continue;
            if(region->fTot_postFit==nullptr) continue;
            // loop on bins
            for(int i_bin=1;i_bin<=region->fTot_postFit->GetNbinsX();i_bin++){
                // set gamma name
                std::string gammaName = Form("stat_%s_bin_%d",region->fName.c_str(),i_bin-1);
                npNames.push_back(gammaName);
                i_np++;
                systNames.push_back( gammaName );
                // find the systematic in the region
                int syst_idx = -1;
                for(int j_syst=0;j_syst<(int)region->fSystNames.size();j_syst++){
                    if(gammaName==region->fSystNames[j_syst]){
                        syst_idx = j_syst;
                    }
                }
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(region->fTot_postFit->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(region->fTot_postFit->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(region->fTotUp_postFit[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(region->fTotDown_postFit[syst_idx]->Clone()));
                }
                h_up.  emplace_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,gammaName.c_str()), Form("h_Tot_%s_Up_TMP",  gammaName.c_str()), Nbin,0,Nbin) );
                h_down.emplace_back( new TH1D(Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Nbin,0,Nbin) );
                h_up[i_np]  ->SetBinContent( i_ch,h_tmp_Up  ->Integral() );
                h_down[i_np]->SetBinContent( i_ch,h_tmp_Down->Integral() );
            }
        }

        // loop on regions
        for(int i_ch=1;i_ch<=Nbin;i_ch++){
            Region *region = fRegions[regionVec[i_ch-1]];
            if(region==nullptr) continue;
            if(region->fTot_postFit==nullptr) continue;
            // loop on bins
            for(int i_bin=1;i_bin<=region->fTot_postFit->GetNbinsX();i_bin++){
                for(auto sample : fSamples){
                    // set gamma name
                    if(!sample->fSeparateGammas) continue;
                    std::string gammaName = Form("shape_stat_%s_%s_bin_%d",sample->fName.c_str(),region->fName.c_str(),i_bin-1);
                    npNames.push_back(gammaName);
                    i_np++;
                    systNames.push_back( gammaName );
                    // find the systematic in the region
                    int syst_idx = -1;
                    for(int j_syst=0;j_syst<(int)region->fSystNames.size();j_syst++){
                        if(gammaName==region->fSystNames[j_syst]){
                            syst_idx = j_syst;
                        }
                    }
                    if(syst_idx<0){
                        h_tmp_Up.reset(static_cast<TH1*>(region->fTot_postFit->Clone()));
                        h_tmp_Down.reset(static_cast<TH1*>(region->fTot_postFit->Clone()));
                    }
                    else{
                        h_tmp_Up.reset(static_cast<TH1*>(region->fTotUp_postFit[syst_idx]->Clone()));
                        h_tmp_Down.reset(static_cast<TH1*>(region->fTotDown_postFit[syst_idx]->Clone()));
                    }
                    h_up.  emplace_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,gammaName.c_str()), Form("h_Tot_%s_Up_TMP",  gammaName.c_str()), Nbin,0,Nbin) );
                    h_down.emplace_back( new TH1D(Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Nbin,0,Nbin) );
                    h_up[i_np]  ->SetBinContent( i_ch,h_tmp_Up  ->Integral() );
                    h_down[i_np]->SetBinContent( i_ch,h_tmp_Down->Integral() );
                }
            }
        }
    }
    // add the norm factors
    for(const auto& inorm : fNormFactors) {
        const std::string normName = inorm->fName;
        if(Common::FindInStringVector(npNames,normName)<0){
            npNames.push_back(normName);
            i_np++;
        }
        else {
            continue;
        }
        systNames.push_back( normName );
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
            // find the systematic in the region
            int syst_idx = -1;
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(normName==fRegions[regionVec[i_bin-1]]->fSystNames[j_syst]){
                    syst_idx = j_syst;
                }
            }
            //
            if(isPostFit){
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp_postFit[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown_postFit[syst_idx]->Clone()));
                }
            }
            else{
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown[syst_idx]->Clone()));
                }
            }
            if(i_bin==1){
                h_up.  emplace_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,normName.c_str()), Form("h_Tot_%s_Up_TMP",  normName.c_str()), Nbin,0,Nbin) );
                h_down.emplace_back( new TH1D(Form("h_Tot_%s_Down_TMP",normName.c_str()), Form("h_Tot_%s_Down_TMP",normName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
        }
    }
    //
    if(isPostFit)  g_err = BuildTotError( h_tot.get(), h_up, h_down, npNames, fFitResults->fCorrMatrix.get() );
    else           g_err = BuildTotError( h_tot.get(), h_up, h_down, npNames );
    //
    p->SetTotBkg(h_tot.get());
    p->BlindData();
    p->SetTotBkgAsym(g_err.get());
    //
    p->ResizeBinLabel(Nbin+1);
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        p->SetBinLabel(i_bin,fRegions[regionVec[i_bin-1]]->fShortLabel.c_str());
    }
    p->Draw(opt);
    if(!isPostFit && p->GetBlindingHisto()) {
        fBlindedBins.reset(static_cast<TH1D*>(p->GetBlindingHisto()->Clone("blinding_trexfit")));
        fBlindedBins->SetDirectory(nullptr);
    }
    //
    if(divisionVec.size()>0){
        p->pad0->cd();
        TLine line(0,1,0,1);
        line.SetNDC(0);
        line.SetLineStyle(7);
        line.SetLineColor(kBlack);
        line.SetLineWidth(2);
        line.DrawLine(divisionVec[0],((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMinimum(),
                      divisionVec[0],std::pow(((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(),0.73) );
        p->pad1->cd();
        line.DrawLine(divisionVec[0],((TH1D*)p->pad1->GetPrimitive("h_dummy2"))->GetMinimum(),
                      divisionVec[0],((TH1D*)p->pad1->GetPrimitive("h_dummy2"))->GetMaximum() );
    }
    //
    p->pad0->cd();
    if(!checkVR){
        if(fSummaryPlotLabels.size()>0){
            TLatex tex;
            tex.SetNDC(0);
            tex.SetTextAlign(20);
            //
            for(unsigned int ii=0;ii<=divisionVec.size();ii++){
                if(fSummaryPlotLabels.size()<ii+1) break;
                double xmax = Nbin;
                double xmin = 0.;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                double xpos = xmin + 0.5*(xmax - xmin);
                double ypos = std::pow(((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(), 0.61 );
                tex.DrawLatex(xpos,ypos,fSummaryPlotLabels[ii].c_str());
            }
        }
    }
    else{
        if(fSummaryPlotValidationLabels.size()>0){
            TLatex tex;
            tex.SetNDC(0);
            tex.SetTextAlign(20);
            //
            for(unsigned int ii=0;ii<=divisionVec.size();ii++){
                if(fSummaryPlotValidationLabels.size()<ii+1) break;
                double xmax = Nbin;
                double xmin = 0.;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                double xpos = xmin + 0.5*(xmax - xmin);
                double ypos = std::pow(((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(), 0.61 );
                tex.DrawLatex(xpos,ypos,fSummaryPlotValidationLabels[ii].c_str());
            }
        }
    }
    //
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        WriteDebugStatus("TRExFit::DrawSummary", std::to_string(i_bin) + ":\t" + std::to_string(h_tot->GetBinContent(i_bin)) + "\t+" +
             std::to_string(g_err->GetErrorYhigh(i_bin-1)) + "\t-" + std::to_string(g_err->GetErrorYlow(i_bin-1)));
    }
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        if(fSummaryPrefix!=""){
            if(isPostFit)  p->SaveAs((fName+"/Plots/"+fSummaryPrefix+"_Summary_postFit"+(checkVR?"_VR":"")+fSuffix+"."+format).c_str());
            else           p->SaveAs((fName+"/Plots/"+fSummaryPrefix+"_Summary"        +(checkVR?"_VR":"")+fSuffix+"."+format).c_str());
        }
        else{
            if(isPostFit)  p->SaveAs((fName+"/Plots/Summary_postFit"+(checkVR?"_VR":"")+fSuffix+"."+format).c_str());
            else           p->SaveAs((fName+"/Plots/Summary"        +(checkVR?"_VR":"")+fSuffix+"."+format).c_str());
        }
    }
    return p;
}

//__________________________________________________________________________________
//
void TRExFit::DrawMergedPlot(std::string opt,std::string group) const{
    std::vector<Region*> regions;
    if(group=="") regions = fRegions;
    else{
        for(auto region : fRegions){
            if(region->fGroup == group) regions.push_back(region);
        }
    }
    bool isPostFit = false;
    if(opt.find("post")!=std::string::npos) isPostFit = true;
    if(TRExFitter::POISSONIZE) opt += " poissonize";
    // start with total prediction, which should be always there
    // build a vector of histograms
    int i_ch = 0;
    std::vector<TH1*> hTotVec;
    std::vector<double> edges;
    std::vector<TGaxis*> xaxis;
    std::vector<TGaxis*> yaxis;
    //
    double ymax0 = -1.; // ymax0 is the max y of the first region
    double ymax  = -1.;
    for(auto region : regions){
        TH1* h_tmp  = nullptr;
        if(isPostFit) h_tmp = (TH1*)region->fTot_postFit->Clone();
        else          h_tmp = (TH1*)region->fTot->Clone();
        TH1* h_data = nullptr;
        if(region->fData!=nullptr) h_data = (TH1*)region->fData->fHist->Clone();
        //
        for(int i_bin=1;i_bin<=h_tmp->GetNbinsX();i_bin++){
            if(isPostFit) h_tmp->SetBinError( i_bin,region->fErr_postFit->GetErrorY(i_bin-1) );
            else          h_tmp->SetBinError( i_bin,region->fErr->GetErrorY(i_bin-1) );
        }
        // find max y (don't rely on GetMaximum, since it could have been modified
        double ymaxTmp = 0;
        for(int i_bin=1;i_bin<=h_tmp->GetNbinsX();i_bin++){
            if(h_tmp->GetBinContent(i_bin)>ymaxTmp) ymaxTmp = h_tmp->GetBinContent(i_bin);
            if(h_data!=nullptr)
                if(h_data->GetBinContent(i_bin)+h_data->GetBinError(i_bin)>ymaxTmp) ymaxTmp = h_data->GetBinContent(i_bin)+h_data->GetBinError(i_bin);
        }
        h_tmp->SetMaximum(ymaxTmp);
        // set max for first hist
        if(ymax0<0){
            ymax0 = ymaxTmp;
            if(TRExFitter::OPTION["MergeYfrac"]==0) TRExFitter::OPTION["MergeYfrac"] = 0.66;
            ymax  = TRExFitter::OPTION["MergeYfrac"]*ymax0;
        }
        hTotVec.push_back(h_tmp);
        if(i_ch>0) edges.push_back(h_tmp->GetXaxis()->GetBinUpEdge(h_tmp->GetNbinsX())-h_tmp->GetXaxis()->GetBinLowEdge(1) + edges[i_ch-1]);
        else       edges.push_back(h_tmp->GetXaxis()->GetBinUpEdge(h_tmp->GetNbinsX()));
        // get xaxes
        if(i_ch>0) xaxis.push_back(new TGaxis(edges[i_ch-1],                      0,edges[i_ch],0, h_tmp->GetXaxis()->GetBinLowEdge(1),h_tmp->GetXaxis()->GetBinUpEdge(h_tmp->GetNbinsX()) ,510,"+"));
        else       xaxis.push_back(new TGaxis(h_tmp->GetXaxis()->GetBinLowEdge(1),0,edges[i_ch],0, h_tmp->GetXaxis()->GetBinLowEdge(1),h_tmp->GetXaxis()->GetBinUpEdge(h_tmp->GetNbinsX()) ,510,"+"));
        // get yaxes
        if(i_ch>0) yaxis.push_back(new TGaxis(edges[i_ch-1],                      0,edges[i_ch-1],                      ymax, 0,ymaxTmp, 510,"-"));
        else       yaxis.push_back(new TGaxis(h_tmp->GetXaxis()->GetBinLowEdge(1),0,h_tmp->GetXaxis()->GetBinLowEdge(1),ymax, 0,ymaxTmp, 510,"-"));
        i_ch ++;
    }
    // then proceed with data, singnal and bkg
    std::vector<TH1*> hDataVec;
    std::vector<std::vector<TH1*>> hSignalVec;
    std::vector<std::vector<TH1*>> hBackgroundVec;
    for(auto sample : fSamples){
        if(sample->fType==Sample::GHOST) continue;
        std::vector<TH1*> tmpVec;
        i_ch = 0;
        for(auto region : regions){
            TH1* h_tmp = nullptr;
            for(auto& sampleHist : region->fSampleHists){
                if(sampleHist->fSample->fName == sample->fName){
                    if(isPostFit && sample->fType!=Sample::DATA){
                        if(sampleHist->fHist_postFit!=nullptr){
                            h_tmp = (TH1*)sampleHist->fHist_postFit->Clone();
                        }
                    }
                    else{
                        if(sampleHist->fHist!=nullptr){
                            h_tmp = (TH1*)sampleHist->fHist->Clone();
                            if(!sampleHist->fSample->fUseMCStat && !sampleHist->fSample->fSeparateGammas){
                                for(int i_bin=0;i_bin<h_tmp->GetNbinsX()+2;i_bin++) h_tmp->SetBinError(i_bin,0.);
                            }
                            // scale it according to NormFactors
                            Common::ScaleNominal(sampleHist.get(), h_tmp);
                        }
                    }
                    break;
                }
            }
            // if the sample was not in the region...
            if(h_tmp==nullptr){
                h_tmp = (TH1*)hTotVec[i_ch]->Clone();
                h_tmp->Scale(0);
            }
            if(sample->fGroup!="") h_tmp->SetTitle(sample->fGroup.c_str());
            else                   h_tmp->SetTitle(sample->fTitle.c_str());
            tmpVec.push_back(h_tmp);
            //
            i_ch ++;
        }
        if(sample->fType==Sample::DATA)            hDataVec = tmpVec;
        else if(sample->fType==Sample::SIGNAL)     hSignalVec.push_back(tmpVec);
        else if(sample->fType==Sample::BACKGROUND) hBackgroundVec.push_back(tmpVec);
    }
    //
    // scale them (but the first region)
    for(unsigned int i_channel=1;i_channel<regions.size();i_channel++){
        double scale = ymax/hTotVec[i_channel]->GetMaximum();
        hTotVec[i_channel]->Scale( scale );
        for(int i_bin=1;i_bin<=hDataVec[i_channel]->GetNbinsX();i_bin++){
            hDataVec[i_channel]->SetBinError(i_bin,std::sqrt(hDataVec[i_channel]->GetBinContent(i_bin)));
        }
        hDataVec[i_channel]->Scale( scale );
        for(auto hVec : hSignalVec)     hVec[i_channel]->Scale( scale );
        for(auto hVec : hBackgroundVec) hVec[i_channel]->Scale( scale );
    }
    //
    // merge and plot them
    std::shared_ptr<TRExPlot> p;
    int cWidthMerge = 1200;
    int cHeightMerge = 600;
    if(fMergeCanvasSize.size()>1){
        cWidthMerge = fMergeCanvasSize.at(0);
        cHeightMerge = fMergeCanvasSize.at(1);
    }
    p = std::make_shared<TRExPlot>(fInputName+"_merge",cWidthMerge,cHeightMerge,TRExFitter::NORATIO);
    //
    p->SetData(Common::MergeHistograms(hDataVec),"");
    for(unsigned int i_sig=0;i_sig<hSignalVec.size();i_sig++){
        if(TRExFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    Common::MergeHistograms(hSignalVec[i_sig]),"");
        if(TRExFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(Common::MergeHistograms(hSignalVec[i_sig]),"");
        if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(Common::MergeHistograms(hSignalVec[i_sig]),"");

    }
    for(unsigned int i_bkg=0;i_bkg<hBackgroundVec.size();i_bkg++) p->AddBackground(Common::MergeHistograms(hBackgroundVec[i_bkg]),"");
    p->SetTotBkg(Common::MergeHistograms(hTotVec));
    //
    p->SetCME(fCmeLabel);
    p->SetLumi(fLumiLabel);
    p->fATLASlabel = fAtlasLabel;
    p->fLabelX = fLabelXMerge;
    p->fLabelY = fLabelYMerge;
    p->fLegendX1 = fLegendX1Merge;
    p->fLegendX2 = fLegendX2Merge;
    p->fLegendY = fLegendYMerge;
    p->fLegendNColumns = fLegendNColumnsMerge;
    p->SetXaxis(regions[0]->fVariableTitle);
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType=="DATA/MC"){
        p->fRatioType = "DATA/BKG";
    }
    if(fBlindingThreshold >= 0) {
        std::unique_ptr<TH1> signal(nullptr);
        if (hSignalVec.size() > 0) {
            signal.reset(static_cast<TH1*>(Common::MergeHistograms(hSignalVec[0])->Clone()));
        }
        std::unique_ptr<TH1> bkg(nullptr);
        for (std::size_t i = 0; i < hBackgroundVec.size(); ++i) {
            if (Common::MergeHistograms(hBackgroundVec[i])) continue;
            if (!bkg) bkg.reset(static_cast<TH1*>(Common::MergeHistograms(hBackgroundVec[i])->Clone()));
            else      bkg->Add(Common::MergeHistograms(hBackgroundVec[i]));
        }
        const std::vector<int>& blindedBins = Common::ComputeBlindedBins(signal.get(),
                                                                         bkg.get(),
                                                                         fBlindingType,
                                                                         fBlindingThreshold);
        p->SetBinBlinding(blindedBins);
    }

    if(TRExFitter::OPTION["MergeYmaxScale"]==0) TRExFitter::OPTION["MergeYmaxScale"] = 1.25;
    p->fYmax = TRExFitter::OPTION["MergeYmaxScale"]*ymax0;
    if(fRatioYmax>0) p->fRatioYmax = fRatioYmax;
    if(fRatioYmin>0) p->fRatioYmin = fRatioYmin;
    if(isPostFit && fRatioYmaxPostFit>0) p->fRatioYmax = fRatioYmaxPostFit;
    if(isPostFit && fRatioYminPostFit>0) p->fRatioYmin = fRatioYminPostFit;
    p->Draw(opt);
    //
    // manipulate canvas / pad
    //
    // dahsed line in ratio
    p->pad1->cd();
    std::vector<TLine*> l;
    for(auto edge : edges){
        TLine *l_tmp = new TLine(edge,((TH1*)p->pad1->GetPrimitive("h_dummy2"))->GetMinimum(),
                                 edge,((TH1*)p->pad1->GetPrimitive("h_dummy2"))->GetMaximum());
        l_tmp->SetLineStyle(kDashed);
        l_tmp->Draw("same");
        l.push_back(l_tmp);
    }
    //
    // (dahsed) line in main pad
    p->pad0->cd();
    for(auto edge : edges){
        TLine *l_tmp = new TLine(edge,0,edge,1.25*ymax);
        if(TRExFitter::OPTION["MergeScaleY"]!=0) l_tmp->SetLineStyle(kDashed);
        l_tmp->Draw("same");
    }
    //
    // y-axis
    p->pad0->cd();
    int i_yaxis = 0;
    for(auto a : yaxis){
        if(i_yaxis>0){
            a->SetLabelFont(gStyle->GetTextFont());
            a->SetLabelSize(gStyle->GetTextSize());
            a->SetNdivisions(805);
            // the following lines require newer version of ROOT
            a->ChangeLabel(1,-1,-1,-1,-1,-1," ");
            a->Draw();
        }
        i_yaxis ++;
    }
    //
    // x-axis
    p->pad0->cd();
    p->pad0->SetTickx(0);
    p->pad0->SetTicky(0);
    TH1* h_dummy = (TH1*)p->pad0->GetPrimitive("h_dummy");
    h_dummy->GetXaxis()->SetLabelSize(0);
    h_dummy->GetXaxis()->SetTickLength(0);
    p->pad0->RedrawAxis();
    for(auto a : xaxis){
        a->SetLabelFont(gStyle->GetTextFont());
        a->SetLabelSize(gStyle->GetTextSize());
        a->SetNdivisions(805);
        ((TGaxis*)(a->DrawClone()))->SetLabelSize(0);
    }
    // ratio
    p->pad1->cd();
    p->pad1->SetTickx(0);
    TH1* h_dummy2 = (TH1*)p->pad1->GetPrimitive("h_dummy2");
    h_dummy2->GetXaxis()->SetLabelSize(0);
    h_dummy2->GetXaxis()->SetTickLength(0);
    p->pad1->RedrawAxis();
    unsigned int i_reg = 0;
    for(auto a : xaxis){
        a->SetLabelFont(gStyle->GetTextFont());
        a->SetLabelSize(gStyle->GetTextSize());
        a->SetNdivisions(805);
        TGaxis *ga = (TGaxis*)a->DrawClone();
        if(fRatioYmin!=0) { ga->SetY1(fRatioYmin); ga->SetY2(fRatioYmin); }
        if(isPostFit && fRatioYminPostFit!=0) { ga->SetY1(fRatioYminPostFit); ga->SetY2(fRatioYminPostFit); }
        ga->SetTitle(regions[i_reg]->fVariableTitle.c_str());
        ga->SetTitleOffset(h_dummy2->GetXaxis()->GetTitleOffset()*0.45*(1200./600.)*(TRExFitter::OPTION["CanvasHeight"]/TRExFitter::OPTION["CanvasWidthMerge"]));
        ga->SetTitleSize(gStyle->GetTextSize());
        ga->SetTitleFont(gStyle->GetTextFont());
        // the following lines require newer version of ROOT
        if(i_reg<regions.size()-1){
            ga->ChangeLabel(-1,-1,-1,-1,-1,-1," "); // shut up the last label ;)
        }
        i_reg ++;
    }
    //
    // additional labels
    p->pad0->cd();
    TLatex *tex = new TLatex();
    tex->SetTextSize(gStyle->GetTextSize());
    tex->SetTextFont(gStyle->GetTextFont());
    for(unsigned int i_channel=0;i_channel<regions.size();i_channel++){
        tex->SetNDC(0);
        tex->DrawLatex(edges[i_channel],1.15*ymax,("#kern[-1]{"+regions[i_channel]->fLabel+" }").c_str());
    }
    //
    tex->SetNDC(1);
    double textHeight = 0.05*(672./p->pad0->GetWh());
    double labelY = 1-0.08*(700./p->c->GetWh());
    if(p->fLabelY>=0) labelY = p->fLabelY;
    labelY -= textHeight - 0.015;
    tex->DrawLatex(0.33,labelY,fLabel.c_str());
    if(isPostFit) tex->DrawLatex(0.33,labelY-textHeight,"Post-fit");
    else          tex->DrawLatex(0.33,labelY-textHeight,"Pre-fit");
    //
    // save image
    std::string saveName = fName+"/Plots/";
    if(fSummaryPrefix!="") saveName += fSummaryPrefix+"_";
    saveName += "Merge";
    if(group!="") saveName += "_"+group;
    if(isPostFit) saveName += "_postFit";
    saveName += fSuffix;
    for(const auto& format : TRExFitter::IMAGEFORMAT){
        p->SaveAs((saveName+"."+format).c_str());
    }
}

//__________________________________________________________________________________
//
void TRExFit::BuildYieldTable(std::string opt, std::string group) const{
    WriteInfoStatus("TRExFit::BuildYieldTable", "-------------------------------------------");
    WriteInfoStatus("TRExFit::BuildYieldTable", "Building Yields Table...");
    if(!TRExFitter::SHOWSTACKSIG || !TRExFitter::ADDSTACKSIG) WriteWarningStatus("TRExFit::BuildYieldTable", "Signal samples not added to \"Tot\" because of \"PlotOptions\" in config file.");
    bool isPostFit = opt.find("post")!=std::string::npos;
    std::ofstream out;
    std::ofstream texout;
    gSystem->mkdir(fName.c_str(),true);
    gSystem->mkdir((fName+"/Tables").c_str());
    std::string suffix = "";
    if(group!="") suffix += "_"+group;
    suffix += fSuffix;
    if(!isPostFit){
        out.open(   (fName+"/Tables/Yields"+suffix+".txt").c_str());
        texout.open((fName+"/Tables/Yields"+suffix+".tex").c_str());
    }
    else{
        out.open(   (fName+"/Tables/Yields_postFit"+suffix+".txt").c_str());
        texout.open((fName+"/Tables/Yields_postFit"+suffix+".tex").c_str());
    }
    // build one bin per region
    std::vector<std::unique_ptr<TH1D> > h_smp(fNSamples);
    std::vector<std::unique_ptr<TGraphAsymmErrors> > g_err(fNSamples);
    std::unique_ptr<TGraphAsymmErrors> g_err_tot(nullptr);
    //
    std::string name;
    std::string title;
    //
    double intErr; // to store the integral error
    TH1* h0; // to store varius histograms temporary

    YamlConverter::TableContainer container;

    //
    // Building region - bin correspondence
    //
    std::vector<int> regionVec;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
            if(group!="" && fRegions[i_ch]->fGroup!=group) continue;
            regionVec.push_back(i_ch);
    }
    if(regionVec.size()==0) return;
    int Nbin = regionVec.size();
    //
    out << " |       | ";
    for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
        out << fRegions[regionVec[i_bin-1]]->fLabel << " | ";
        container.regionNames.emplace_back(fRegions[regionVec[i_bin-1]]->fLabel);
    }
    out << std::endl;
    if(fTableOptions.find("STANDALONE")!=std::string::npos){
        texout << "\\documentclass[10pt]{article}" << std::endl;
        texout << "\\usepackage{siunitx}" << std::endl;
        texout << "\\sisetup{separate-uncertainty,table-format=6.3(6)}  % hint: modify table-format to best fit your tables" << std::endl;
        texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << std::endl;
        texout << "\\begin{document}" << std::endl;
    }
    // if not STANDALONE, add a comment in the tex saying that one needs to include siunitx
    else{
        texout << "% NB: add to main document: " << std::endl;
        texout << "% \\usepackage{siunitx} " << std::endl;
        texout << "% \\sisetup{separate-uncertainty,table-format=6.3(6)}  % hint: modify table-format to best fit your tables" << std::endl;
    }
    if(fTableOptions.find("LANDSCAPE")!=std::string::npos){
        texout << "\\begin{landscape}" << std::endl;
    }
    texout << "\\begin{table}[htbp]" << std::endl;
    texout << "\\begin{center}" << std::endl;
    if(fTableOptions.find("FOOTNOTESIZE")!=std::string::npos){
        texout << "\\footnotesize" << std::endl;
    }
    texout << "\\begin{tabular}{|l" ;
    for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
        texout << "|S";
    }
    texout << "|}" << std::endl;
    texout << "\\hline " << std::endl;
    for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
        if(fRegions[regionVec[i_bin-1]]->fTexLabel!="") texout << " & {" << fRegions[regionVec[i_bin-1]]->fTexLabel << "}";
        else                                            texout << " & {" << fRegions[regionVec[i_bin-1]]->fLabel    << "}";
    }
    texout << "\\\\" << std::endl;
    texout << "\\hline " << std::endl;
    //
    std::vector< std::string > titleVec;
    std::vector< int > idxVec;
    std::shared_ptr<SampleHist> sh = nullptr;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        name = fSamples[i_smp]->fName;
        title = fSamples[i_smp]->fTitle;
        //
        int idx = Common::FindInStringVector(titleVec,title);
        if(idx>=0){
            idxVec.push_back(idx);
        }
        else{
            idxVec.push_back(i_smp);
            h_smp[idxVec[i_smp]].reset(new TH1D(("h_"+name).c_str(),title.c_str(), Nbin,0,Nbin));
        }
        for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
            sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
            if(sh!=nullptr){
                if(isPostFit && fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::GHOST)
                    h0 = sh->fHist_postFit.get();
                else
                    h0 = sh->fHist.get();
                double tmpErr = h_smp[idxVec[i_smp]]->GetBinError(i_bin); // Michele -> get the error before adding content to bin, to avoid ROOT automatically increasing it!
                double scale = 1.;
                if (!isPostFit){
                    scale = Common::GetNominalMorphScale(sh.get());
                }
                h_smp[idxVec[i_smp]]->AddBinContent( i_bin,scale*h0->IntegralAndError(1,h0->GetNbinsX(),intErr) );
                intErr*=scale;
                if( (isPostFit && fUseGammaPulls) || !fUseStatErr || (!sh->fSample->fUseMCStat && !sh->fSample->fSeparateGammas))
                    h_smp[idxVec[i_smp]]->SetBinError(i_bin,0.);
                else
                    h_smp[idxVec[i_smp]]->SetBinError(i_bin, std::hypot(tmpErr, intErr));
            }
        }
        titleVec.push_back(title);
    }
    //
    // build a global list of systematics (including gammas), from the lists already created from each region
    // - for pre-fit the list contains only regular systematics and SHAPE systematics (no norm-factors, no stat gammas)
    // - for post-fit also norm-factors and stat gammas (only in the case of UseGammaPulls)
    std::vector<std::string> globalSystNames;
    std::vector<std::string> globalNpNames;
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        Region *reg = fRegions[regionVec[i_bin-1]];
        for(int i_syst=0;i_syst<(int)reg->fSystNames.size();i_syst++){
            std::string systName = reg->fSystNames[i_syst];
            std::string systNuisPar = systName;
            if(TRExFitter::NPMAP[systName]!="") systNuisPar = TRExFitter::NPMAP[systName];
            if (std::find(globalNpNames.begin(), globalNpNames.end(), systNuisPar) != globalNpNames.end()) continue;
            globalSystNames.push_back( systName );
            globalNpNames.push_back(systNuisPar);
        }
    }
    //
    // add tot uncertainty on each sample
    int i_np = -1;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        if(idxVec[i_smp]!=i_smp) continue;
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        name = fSamples[i_smp]->fName;
        // build the vectors of variations
        std::vector< std::shared_ptr<TH1> > h_up;
        std::vector< std::shared_ptr<TH1> > h_down;
        std::unique_ptr<TH1> h_tmp_Up;
        std::unique_ptr<TH1> h_tmp_Down;
        std::vector<std::string> npNames;
        i_np = -1;
        //
        // loop on the global list of systematics
        for(int i_syst=0;i_syst<(int)globalSystNames.size();i_syst++){
            std::string systName    = globalSystNames.at(i_syst);
            std::string systNuisPar = globalNpNames.at(i_syst);
            // if post-fit but no UseGammaPulls, skip stat gammas
            if(isPostFit && systName.find("stat_")!=std::string::npos && !fUseGammaPulls){
                continue;
            }
            npNames.push_back( systNuisPar );
            i_np++;
            for(int i_bin=1;i_bin<=Nbin;i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                //
                // find the systematic in the region
                int syst_idx = -1;
                for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                    if(systName==fRegions[regionVec[i_bin-1]]->fSystNames[j_syst]){
                        syst_idx = j_syst;
                    }
                }
                //
                if(sh!=nullptr){
                    if(isPostFit){
                        if(syst_idx<0 || sh->GetSystematic(systName)==nullptr){
                            h_tmp_Up.reset(static_cast<TH1*>(sh->fHist_postFit->Clone()));
                            h_tmp_Down.reset(static_cast<TH1*>(sh->fHist_postFit->Clone()));
                        }
                        else{
                            h_tmp_Up.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistUp_postFit->Clone()));
                            h_tmp_Down.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistDown_postFit->Clone()));
                        }
                    }
                    else {
                        if(syst_idx<0 || sh->GetSystematic(systName)==nullptr){
                            h_tmp_Up.reset(static_cast<TH1*>(sh->fHist->Clone()));
                            h_tmp_Down.reset(static_cast<TH1*>(sh->fHist->Clone()));
                        }
                        else{
                            h_tmp_Up.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistUp->Clone()));
                            h_tmp_Down.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistDown->Clone()));
                        }
                    }
                }
                else {
                    h_tmp_Up.reset(new TH1D(Form("h_DUMMY_%s_up_%i",  systName.c_str(),i_bin-1),"h_dummy",1,0,1));
                    h_tmp_Down.reset(new TH1D(Form("h_DUMMY_%s_down_%i",systName.c_str(),i_bin-1),"h_dummy",1,0,1));
                }
                if(i_bin==1){
                    h_up.  emplace_back( new TH1D(Form("h_%s_%s_Up_TMP",  name.c_str(),systName.c_str()),Form("h_%s_%s_Up_TMP",  name.c_str(),systName.c_str()), Nbin,0,Nbin) );
                    h_down.emplace_back( new TH1D(Form("h_%s_%s_Down_TMP",name.c_str(),systName.c_str()),Form("h_%s_%s_Down_TMP",name.c_str(),systName.c_str()), Nbin,0,Nbin) );
                }
                double scale = 1.;
                if (!isPostFit){
                    scale = Common::GetNominalMorphScale(sh.get());
                }
                h_up[i_np]  ->SetBinContent( i_bin,(h_tmp_Up  ->Integral(1,h_tmp_Up  ->GetNbinsX()))*scale );
                h_down[i_np]->SetBinContent( i_bin,(h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()))*scale );
                //
                // eventually add any other samples with the same title
                for(int j_smp=0;j_smp<fNSamples;j_smp++){
                    sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( fSamples[j_smp]->fName );
                    if(sh==nullptr) continue;
                    if(idxVec[j_smp]==i_smp && i_smp!=j_smp){
                        if(isPostFit){
                            if(syst_idx<0 || sh->GetSystematic(systName)==nullptr){
                                h_tmp_Up.reset(static_cast<TH1*>(sh->fHist_postFit->Clone()));
                                h_tmp_Down.reset(static_cast<TH1*>(sh->fHist_postFit->Clone()));
                            }
                            else{
                                h_tmp_Up.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistUp_postFit->Clone()));
                                h_tmp_Down.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistDown_postFit->Clone()));
                            }
                        }
                        else{
                            if(syst_idx<0 || sh->GetSystematic(systName)==nullptr){
                                h_tmp_Up.reset(static_cast<TH1*>(sh->fHist->Clone()));
                                h_tmp_Down.reset(static_cast<TH1*>(sh->fHist->Clone()));
                            }
                            else{
                                h_tmp_Up.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistUp->Clone()));
                                h_tmp_Down.reset(static_cast<TH1*>(sh->GetSystematic(systName)->fHistDown->Clone()));
                            }
                        }
                        double morph_scale = 1.;
                        if (!isPostFit) {
                            morph_scale = Common::GetNominalMorphScale(sh.get());
                        }
                        h_up[i_np]  ->AddBinContent( i_bin,(h_tmp_Up  ->Integral(1,h_tmp_Up->GetNbinsX()))*morph_scale );
                        h_down[i_np]->AddBinContent( i_bin,(h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()))*morph_scale );
                    }
                }
            }
        }
        //
        if(isPostFit)  g_err[i_smp] = BuildTotError( h_smp[i_smp].get(), h_up, h_down, npNames, fFitResults->fCorrMatrix.get() );
        else           g_err[i_smp] = BuildTotError( h_smp[i_smp].get(), h_up, h_down, npNames );
    }
    //
    // Print samples except ghosts, data for blind fits, signal for B-only...
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if( fSamples[i_smp]->fType==Sample::GHOST ) continue;
        if( fSamples[i_smp]->fType==Sample::DATA  ) continue;
        if( fSamples[i_smp]->fType==Sample::SIGNAL && (fFitType==FitType::BONLY && isPostFit) ) continue;
        if(idxVec[i_smp]!=i_smp) continue;
        //
        // print values
        out << " | " << fSamples[i_smp]->fTitle << " | ";
        container.sampleNames.emplace_back(fSamples[i_smp]->fTitle);
        std::vector<double> yamlTmpYields;
        std::vector<double> yamlTmpErrors;
        if(fSamples[i_smp]->fType==Sample::DATA) texout << "\\hline " << std::endl;
        if(fSamples[i_smp]->fTexTitle!="") texout << "  " << fSamples[i_smp]->fTexTitle << "  ";
        else                               texout << "  " << fSamples[i_smp]->fTitle << "  ";
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
            double mean = h_smp[i_smp]->GetBinContent(i_bin);
            double uncertainty = ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
            double mean_rounded = mean;
            double uncertainty_rounded = uncertainty;
            yamlTmpYields.emplace_back(mean);
            yamlTmpErrors.emplace_back(uncertainty);
            int n = -1; // this will contain the number of decimal places
            if (fUseATLASRoundingTxt || fUseATLASRoundingTex){
                n = Common::ApplyATLASrounding(mean_rounded, uncertainty_rounded);
            }
            if(fUseATLASRoundingTxt){
                out << mean_rounded << " pm " << uncertainty_rounded << " | ";
            }
            else{
                out << mean << " pm " << uncertainty << " | ";
            }
            if(fUseATLASRoundingTex){
                texout << " & ";
                if(n<0) texout << mean_rounded;
                else    texout << Form(("%."+std::to_string(n)+"f").c_str(),mean_rounded);
                if(uncertainty==0){ // to fix Latex siunitx issue
                    if(n<0) texout << " (" << uncertainty_rounded << ")";
                    else    texout << " (" << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded) << ")";
                }
                else{
                    if(n<0) texout << " \\pm " << uncertainty_rounded;
                    else    texout << " \\pm " << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded);
                }
            }
            else{
                texout << " & " << mean << " \\pm " << uncertainty;
            }
        }
        out << std::endl;
        texout << " \\\\ ";
        texout << std::endl;
        container.mcYields.emplace_back(yamlTmpYields);
        container.mcErrors.emplace_back(yamlTmpErrors);
    }

    //
    // Build tot
    //
    std::unique_ptr<TH1D> h_tot = std::make_unique<TH1D>("h_Tot_","h_Tot", Nbin,0,Nbin);
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot_postFit->IntegralAndError(1,fRegions[regionVec[i_bin-1]]->fTot_postFit->GetNbinsX(),intErr) );
        else          h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot->IntegralAndError(        1,fRegions[regionVec[i_bin-1]]->fTot->GetNbinsX(),        intErr) );
        h_tot->SetBinError( i_bin, intErr );
    }
    //
    //   Build error band
    // build the vectors of variations
    std::vector< std::shared_ptr<TH1> > h_up;
    std::vector< std::shared_ptr<TH1> > h_down;
    std::unique_ptr<TH1> h_tmp_Up;
    std::unique_ptr<TH1> h_tmp_Down;
    std::vector<std::string> npNames;
    i_np = -1;
    //
    // loop on the global list of systematics
    for(int i_syst=0;i_syst<(int)globalSystNames.size();i_syst++){
        std::string systName    = globalSystNames.at(i_syst);
        std::string systNuisPar = globalNpNames.at(i_syst);
        // if post-fit but no UseGammaPulls, skip stat gammas
        if(isPostFit && systName.find("stat_")!=std::string::npos && !fUseGammaPulls){
            continue;
        }
        npNames.push_back( systNuisPar );
        i_np++;
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
            // find the systematic in the region
            int syst_idx = -1;
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(systName==fRegions[regionVec[i_bin-1]]->fSystNames[j_syst]){
                    syst_idx = j_syst;
                }
            }
            //
            if(isPostFit){
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp_postFit[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown_postFit[syst_idx]->Clone()));
                }
            }
            else{
                if(syst_idx<0){
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                }
                else{
                    h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp[syst_idx]->Clone()));
                    h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown[syst_idx]->Clone()));
                }
            }
            if(i_bin==1){
                h_up.  emplace_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,systName.c_str()), Form("h_Tot_%s_Up_TMP",  systName.c_str()), Nbin,0,Nbin) );
                h_down.emplace_back( new TH1D(Form("h_Tot_%s_Down_TMP",systName.c_str()), Form("h_Tot_%s_Down_TMP",systName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
            //
            // look for other syst with the same np
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(j_syst==syst_idx) continue;
                if(systNuisPar==TRExFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
                    std::unique_ptr<TH1> h_tmp = nullptr;
                    if(isPostFit){
                        h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp_postFit[j_syst]->Clone()));
                        h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown_postFit[j_syst]->Clone()));
                        h_tmp.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot_postFit->Clone()));
                    }
                    else{
                        h_tmp_Up.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotUp[j_syst]->Clone()));
                        h_tmp_Down.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTotDown[j_syst]->Clone()));
                        h_tmp.reset(static_cast<TH1*>(fRegions[regionVec[i_bin-1]]->fTot->Clone()));
                    }
                    h_up[i_np]  ->AddBinContent( i_bin,h_tmp_Up  ->Integral()-h_tmp->Integral() );
                    h_down[i_np]->AddBinContent( i_bin,h_tmp_Down->Integral()-h_tmp->Integral() );
                }
            }
        }
    }
    //
    if(isPostFit)  g_err_tot = BuildTotError( h_tot.get(), h_up, h_down, npNames, fFitResults->fCorrMatrix.get() );
    else           g_err_tot = BuildTotError( h_tot.get(), h_up, h_down, npNames );
    //
    if(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) out << " | Total | ";
    else                                                    out << " | Tot.Bkg. | ";
    texout << "\\hline " << std::endl;
    if(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) texout << "  Total ";
    else                                                    texout << "  Total background ";
    if (TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) container.sampleNames.emplace_back("Total");
    else                                                     container.sampleNames.emplace_back("Total background");
    std::vector<double> yamlTmpYields;
    std::vector<double> yamlTmpErrors;
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        double mean = h_tot->GetBinContent(i_bin);
        double uncertainty = ( g_err_tot->GetErrorYhigh(i_bin-1) + g_err_tot->GetErrorYlow(i_bin-1) )/2.;
        double mean_rounded = mean;
        double uncertainty_rounded = uncertainty;
        yamlTmpYields.emplace_back(mean);
        yamlTmpErrors.emplace_back(uncertainty);
        int n = -1; // this will contain the number of decimal places
        if (fUseATLASRoundingTxt || fUseATLASRoundingTex){
            n = Common::ApplyATLASrounding(mean_rounded, uncertainty_rounded);
        }
        if(fUseATLASRoundingTxt){
            out << mean_rounded << " pm " << uncertainty_rounded << " | ";
        }
        else{
            out << mean << " pm " << uncertainty << " | ";
        }
        if(fUseATLASRoundingTex){
            texout << " & ";
            if(n<0) texout << mean_rounded;
            else    texout << Form(("%."+std::to_string(n)+"f").c_str(),mean_rounded);
            if(uncertainty==0){ // to fix Latex siunitx issue
                if(n<0) texout << " (" << uncertainty_rounded << ")";
                else    texout << " (" << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded) << ")";
            }
            else{
                if(n<0) texout << " \\pm " << uncertainty_rounded;
                else    texout << " \\pm " << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded);
            }
        }
        else{
            texout << " & " << mean << " \\pm " << uncertainty;
        }
    }
    container.mcYields.emplace_back(yamlTmpYields);
    container.mcErrors.emplace_back(yamlTmpErrors);
    out << std::endl;
    texout << " \\\\ ";
    texout << std::endl;

    //
    // Print data
    if( !fFitIsBlind ){
        texout << "\\hline " << std::endl;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            std::vector<double> yamlTempDataYields;
            if( fSamples[i_smp]->fType!=Sample::DATA  ) continue;
            if(idxVec[i_smp]!=i_smp) continue;
            //
            // print values
            out << " | " << fSamples[i_smp]->fTitle << " | ";
            if(fSamples[i_smp]->fTexTitle!="") texout << "  " << fSamples[i_smp]->fTexTitle << "  ";
            else                               texout << "  " << fSamples[i_smp]->fTitle << "  ";
            container.dataNames.emplace_back(fSamples[i_smp]->fTitle);
            for(int i_bin=1;i_bin<=Nbin;i_bin++){
                yamlTempDataYields.emplace_back(h_smp[i_smp]->GetBinContent(i_bin));
                texout << " & ";
                out << h_smp[i_smp]->GetBinContent(i_bin);
                texout << Form("%.0f",h_smp[i_smp]->GetBinContent(i_bin));
                out << " | ";
            }
            out << std::endl;
            texout << " \\\\ ";
            texout << std::endl;
            container.dataYields.emplace_back(yamlTempDataYields);
        }
    }
    //
    texout << "\\hline " << std::endl;
    texout << "\\end{tabular} " << std::endl;
    texout << "\\caption{Yields of the analysis} " << std::endl;
    texout << "\\end{center} " << std::endl;
    texout << "\\end{table} " << std::endl;
    if(fTableOptions.find("LANDSCAPE")!=std::string::npos){
        texout << "\\end{landscape}" << std::endl;
    }
    if(fTableOptions.find("STANDALONE")!=std::string::npos){
        texout << "\\end{document}" << std::endl;
    }
    if(fCleanTables){
        std::string shellcommand = "cat "+fName+"/Tables/Yields"+suffix+".tex|sed -e \"s/\\#/ /g\" > "+fName+"/Tables/Yields";
        if(isPostFit) shellcommand += "_postFit";
        shellcommand += suffix+"_clean.tex";
        gSystem->Exec(shellcommand.c_str());
    }

    // Write YAML
    YamlConverter converter{};
    converter.WriteTables(container, fName, isPostFit);
    if (fHEPDataFormat) {
        converter.WriteTablesHEPData(container, fName, isPostFit);
    }
}

//__________________________________________________________________________________
//
void TRExFit::DrawSignalRegionsPlot(int nCols,int nRows) const{
    std::vector< Region* > vRegions;
    if(fRegionsToPlot.size()>0){
        nCols = 1;
        nRows = 1;
        // first loop
        int nRegInRow = 0;
        for(unsigned int i=0;i<fRegionsToPlot.size();i++){
            WriteDebugStatus("TRExFit::DrawSignalRegionsPlot", "Regions to Plot: " + fRegionsToPlot[i]);
            if(fRegionsToPlot[i].find("ENDL")!=std::string::npos){
                nRows++;
                if(nRegInRow>nCols) nCols = nRegInRow;
                nRegInRow = 0;
            }
            else{
                vRegions.push_back( GetRegion(fRegionsToPlot[i]) );
                nRegInRow ++;
            }
        }
    }
    else{
        vRegions = fRegions;
    }
    DrawSignalRegionsPlot(nCols,nRows,vRegions);
}

//__________________________________________________________________________________
//
void TRExFit::DrawSignalRegionsPlot(int nCols,int nRows, std::vector < Region* > &regions) const{
    gSystem->mkdir(fName.c_str(), true);
    double Hp = 250.; // height of one mini-plot, in pixels
    double Wp = 200.; // width of one mini-plot, in pixels
    double H0 = 100.; // height of the top label pad
    if(TRExFitter::OPTION["FourTopStyle"]!=0) {
        H0 = 75; // height of the top label pad
        Hp = 200;
    }
    if(TRExFitter::OPTION["SignalRegionSize"]!=0){
        Hp = TRExFitter::OPTION["SignalRegionSize"];
        Wp = (200./250.)*TRExFitter::OPTION["SignalRegionSize"];
    }
    double H = H0 + nRows*Hp; // tot height of the canvas
    double W = nCols*Wp; // tot width of the canvas
    if(TRExFitter::OPTION["FourTopStyle"]!=0) W += 50.; // FIXME
    else W += 0.1; // to fix eps format (why is this needed?)

    TCanvas c("c","c",W,H);
    TPad pTop("c0","c0",0,1-H0/H,1,1);
    pTop.Draw();
    pTop.cd();
    ATLASLabel(0.1/(W/200.),1.-0.3*(100./H0),fAtlasLabel.c_str());
    myText(    0.1/(W/200.),1.-0.6*(100./H0),1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(fLabel!="-") myText(    0.1/(W/200.),1.-0.9*(100./H0),1,Form("%s",fLabel.c_str()));

    std::unique_ptr<TLegend> leg = nullptr;
    if(TRExFitter::OPTION["FourTopStyle"]!=0){
        leg = std::make_unique<TLegend>(0.35,1.-0.6*(100./H0),1,1.-0.3*(100./H0));
        leg->SetNColumns(3);
        leg->SetFillStyle(0.);
        leg->SetBorderSize(0.);
        leg->SetTextSize( gStyle->GetTextSize() );
        leg->SetTextFont( gStyle->GetTextFont() );
        leg->Draw();
    }

    c.cd();

    TPad pLeft("c1","c1",0,0,0+(W-nCols*Wp)/W,1-H0/H);
    pLeft.Draw();
    pLeft.cd();
    TLatex tex0{};
    tex0.SetNDC();
    tex0.SetTextAngle(90);
    tex0.SetTextAlign(23);
    tex0.DrawLatex(0.4,0.5,"S / #sqrt{ B }");

    c.cd();

    TPad pBottom("c1","c1",0+(W-nCols*Wp)/W,0,1,1-H0/H);
    pBottom.Draw();
    pBottom.cd();

    pBottom.Divide(nCols,nRows);
    unsigned int Nreg = static_cast<unsigned int> (nRows*nCols);
    if(Nreg>regions.size()) Nreg = regions.size();
    std::vector<TH1D> h;
    h.reserve(Nreg);
    std::vector<double> S(Nreg);
    std::vector<double> B(Nreg);
    std::vector<double> xbins = {0.0,0.1,0.9,1.0};
    TLatex tex{};
    tex.SetNDC();
    tex.SetTextSize(gStyle->GetTextSize());
    pBottom.cd(1);

    //
    // Get the values
    //
    for(unsigned int i=0;i<Nreg;i++){
        S[i] = 0.;
        B[i] = 0.;
        if(regions[i]==nullptr) continue;
        for(int i_sig=0;i_sig<regions[i]->fNSig;i_sig++) {
            if(regions[i]->fSig[i_sig]!=nullptr) {
                const double scale = Common::GetNominalMorphScale(regions[i]->fSig[i_sig].get());
                S[i] += scale * regions[i]->fSig[i_sig]->fHist->Integral();
            }
        }
        for(int i_bkg=0;i_bkg<regions[i]->fNBkg;i_bkg++){
            if(regions[i]->fBkg[i_bkg]!=nullptr) {
                const double scale = Common::GetNominalMorphScale(regions[i]->fBkg[i_bkg].get());
                B[i] += scale * regions[i]->fBkg[i_bkg]->fHist->Integral();
            }
        }
        // to avoid nan or inf...
        if(B[i]==0) B[i] = 1e-10;
        // scale up for projections
        if(fLumiScale!=1){
            S[i]*=fLumiScale;
            B[i]*=fLumiScale;
        }
    }
    //
    double yMax = 0;
    //
    bool hasSR = false;
    bool hasCR = false;
    bool hasVR = false;
    for(unsigned int i=0;i<Nreg;i++){
        if(regions[i]==nullptr) continue;
        pBottom.cd(i+1);
        if(TRExFitter::OPTION["LogSignalRegionPlot"]) gPad->SetLogy();
        if(TRExFitter::OPTION["FourTopStyle"]!=0){
            gPad->SetLeftMargin(0.1);
            gPad->SetRightMargin(0.);
        }
        std::string label = regions[i]->fShortLabel;
        h.emplace_back(Form("h[%d]",i),label.c_str(),3,&xbins[0]);
        h.back().SetBinContent(2,S[i]/std::sqrt(B[i]));
        if(TRExFitter::OPTION["FourTopStyle"]==0) h.back().GetYaxis()->SetTitle("S / #sqrt{B}");
        h.back().GetYaxis()->CenterTitle();
        h.back().GetYaxis()->SetLabelOffset(1.5*h.back().GetYaxis()->GetLabelOffset() / (Wp/200.));
        h.back().GetYaxis()->SetTitleOffset(9*nRows/4. );
        if(Wp<200) h.back().GetYaxis()->SetTitleOffset( h.back().GetYaxis()->GetTitleOffset()*0.90 );
        h.back().GetYaxis()->SetLabelSize( h.back().GetYaxis()->GetLabelSize() * (Wp/200.) );
        if(TRExFitter::OPTION["FourTopStyle"]!=0) h.back().GetYaxis()->SetLabelSize( h.back().GetYaxis()->GetLabelSize() * 1.1 );
        h.back().GetXaxis()->SetTickLength(0);
        if(TRExFitter::OPTION["LogSignalRegionPlot"]==0) h.back().GetYaxis()->SetNdivisions(3);
        else TGaxis::SetMaxDigits(5);
        yMax = std::max(yMax,h.back().GetMaximum());
        h.back().GetXaxis()->SetLabelSize(0);
        h.back().SetLineWidth(1);
        h.back().SetLineColor(kBlack);
        if(regions[i]->fRegionType==Region::SIGNAL)          h.back().SetFillColor(kRed+1);
        else if(regions[i]->fRegionType==Region::VALIDATION) h.back().SetFillColor(kGray);
        else                                                 h.back().SetFillColor(kAzure-4);
        if(leg!=nullptr){
            if(regions[i]->fRegionType==Region::CONTROL && !hasCR)    {
                leg->AddEntry(&h.back(),"Control Regions","f");
                hasCR = true;
            }
            if(regions[i]->fRegionType==Region::VALIDATION && !hasVR) {
                leg->AddEntry(&h.back(),"Validation Regions","f");
                hasVR = true;
            }
            if(regions[i]->fRegionType==Region::SIGNAL && !hasSR)     {
                leg->AddEntry(&h.back(),"Signal Regions","f");
                hasSR = true;
            }
        }
        h.back().Draw();
        gPad->SetLeftMargin( gPad->GetLeftMargin()*2.4 );
        gPad->SetRightMargin(gPad->GetRightMargin()*0.1);
        gPad->SetTicky(0);
        gPad->RedrawAxis();
        if(TRExFitter::OPTION["FourTopStyle"]==0) tex.DrawLatex(0.42,0.85,label.c_str());
        else                                      tex.DrawLatex(0.27,0.85,label.c_str());
        const double SoB = S[i]/B[i];
        std::string SB = Form("%.1f%%",(100.*SoB));
        if(TRExFitter::OPTION["FourTopStyle"]!=0){
            if( (100.*SoB)<0.1 ){
                SB = Form("%.0e%%",SoB);
                if(SB.find("0")!=std::string::npos) SB.replace(SB.find("0"), 1, "");
                if(SB.find("e")!=std::string::npos) SB.replace(SB.find("e"), 1, "#scale[0.75]{#times}10^{");
                if(SB.find("%")!=std::string::npos) SB.replace(SB.find("%"), 1, "}");
            }
        }
        SB = "#scale[0.75]{S/B} = "+SB;
        if(TRExFitter::OPTION["FourTopStyle"]==0) tex.DrawLatex(0.42,0.72,SB.c_str());
        else                                      tex.DrawLatex(0.27,0.72,SB.c_str());
    }
    //
    for(unsigned int i=0;i<Nreg;i++){
        if(regions[i]==nullptr) continue;
        if ((h.size() - 1)  <= i) break;
        if(TRExFitter::OPTION["LogSignalRegionPlot"]!=0){
            h[i].SetMaximum(yMax*200);
            h[i].SetMinimum(2e-4);
        }
        else{
            h[i].SetMaximum(yMax*1.5);
            h[i].SetMinimum(0.);
        }
    }
    //
    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fName+"/SignalRegions"+fSuffix+"."+format).c_str());
    }

}

//__________________________________________________________________________________
//
void TRExFit::DrawPieChartPlot(const std::string &opt, int nCols,int nRows) const{

    std::vector< Region* > vRegions;
    if(fRegionsToPlot.size()>0){
        nCols = 1;
        nRows = 1;
        // first loop
        int nRegInRow = 0;
        for(unsigned int i=0;i<fRegionsToPlot.size();i++){
            WriteDebugStatus("TRExFit::DrawPieChartPlot", "Regions to plot: " + fRegionsToPlot[i]);
            if(fRegionsToPlot[i].find("ENDL")!=std::string::npos){
                nRows++;
                if(nRegInRow>nCols) nCols = nRegInRow;
                nRegInRow = 0;
            }
            else{
                vRegions.push_back( GetRegion(fRegionsToPlot[i]) );
                nRegInRow ++;
            }
        }
    }
    else{
        vRegions = fRegions;
    }
    DrawPieChartPlot(opt, nCols,nRows,vRegions);

}


//__________________________________________________________________________________
//
void TRExFit::DrawPieChartPlot(const std::string &opt, int nCols,int nRows, std::vector < Region* > &regions ) const{

    double Hp = 250.; // height of one mini-plot, in pixels
    double Wp = 250.; // width of one mini-plot, in pixels
    double H0 = 100.; // height of the top label pad
    if(TRExFitter::OPTION["FourTopStyle"]>0) H0 = 75; // height of the top label pad

    if(TRExFitter::OPTION["PieChartSize"]!=0){
        Hp = TRExFitter::OPTION["PieChartSize"];
        Wp = TRExFitter::OPTION["PieChartSize"];
    }

    double H = H0 + nRows*Hp; // tot height of the canvas
    double W = nCols*Wp; // tot width of the canvas

    bool isPostFit = opt.find("post")!=std::string::npos;

    //
    // Create the canvas
    //
    if (fPieChartCanvasSize.size() != 0){
        W = fPieChartCanvasSize.at(0);
        H = fPieChartCanvasSize.at(1);
    }
    TCanvas c("c","c",W,H);
    TPad pTop("c0","c0",0,1-H0/H,1,1);
    pTop.Draw();
    pTop.cd();

    if(TRExFitter::OPTION["FourTopStyle"]>0){
        ATLASLabel(0.1/(W/200.),1.-0.3*(100./H0),fAtlasLabel.c_str());
        myText(    0.1/(W/200.),1.-0.6*(100./H0),1,Form("#sqrt{s} = %s",fCmeLabel.c_str()));
        if(fLabel!="-") myText(    0.1/(W/200.),1.-0.9*(100./H0),1,Form("%s",fLabel.c_str()));
    }
    else{
        ATLASLabel(0.05 / (W/200),0.7,fAtlasLabel.c_str());
        myText(    0.05 / (W/200),0.4,1,Form("#sqrt{s} = %s",fCmeLabel.c_str()));
        if(fLabel!="-") myText(    0.05 / (W/200),0.1,1,Form("%s",fLabel.c_str()));
    }

    c.cd();
    TPad pBottom("c1","c1",0,0,1,1-H0/H);
    pBottom.Draw();
    pBottom.cd();
    pBottom.Divide(nCols,nRows);
    int Nreg = nRows*nCols;
    if(Nreg>(int)regions.size()) Nreg = regions.size();
    TLatex tex{};
    tex.SetNDC();
    tex.SetTextSize(gStyle->GetTextSize());
    pBottom.cd(1);

    //
    // Create the map to store all the needed information
    //
    std::map < std::string, int > map_for_legend;
    std::vector < std::map < std::string, double > > results;
    std::vector < std::map < std::string, int > > results_color;

    //
    // Get the values
    //
    for(int i=0;i<Nreg;i++){
        std::map < std::string, double > temp_map_for_region;
        std::map < std::string, int > temp_map_for_region_color;

        if(regions[i]!=nullptr){
            for(int i_bkg=regions[i]->fNBkg-1;i_bkg>=0;i_bkg--){
                if(regions[i]->fBkg[i_bkg]!=nullptr){
                    std::string title = regions[i]->fBkg[i_bkg]->fSample->fTitle;
                    if(regions[i]->fBkg[i_bkg]->fSample->fGroup != "") title = regions[i]->fBkg[i_bkg]->fSample->fGroup.c_str();

                    double integral = 0;
                    if(!isPostFit) integral = regions[i]->fBkg[i_bkg]->fHist->Integral() * fLumiScale;
                    else integral = regions[i]->fBkg[i_bkg]->fHist_postFit->Integral() * fLumiScale;

                    if(temp_map_for_region.find(title)!=temp_map_for_region.end()){
                        temp_map_for_region[title] += integral;
                    } else {
                        temp_map_for_region.insert( std::pair < std::string, double > (title,integral) );
                        temp_map_for_region_color.insert( std::pair < std::string, int > (title, regions[i]->fBkg[i_bkg]->fSample->fFillColor) );
                    }
                    map_for_legend[title] = regions[i]->fBkg[i_bkg]->fSample->fFillColor;
                }
            }
        }
        results.push_back(temp_map_for_region);
        results_color.push_back(temp_map_for_region_color);
    }

    //
    // Finally writting the pie chart
    //
    std::vector<std::unique_ptr<TPie> > pie;
    for(int i=0;i<Nreg;i++){
        if(regions[i]==nullptr) continue;
        pBottom.cd(i+1);
        std::string label = regions[i]->fShortLabel;

        const unsigned int back_n = results[i].size();
        std::vector<double> values(back_n);
        std::vector<int> colors(back_n);
        for( unsigned int iTemp = 0; iTemp < back_n; ++iTemp ){
            values[iTemp] = 0.;
            colors[iTemp] = 0;
        }

        int count = 0;
        for ( std::pair < std::string, double > temp_pair : results[i] ){
            values[count] = temp_pair.second;
            colors[count] = results_color[i][temp_pair.first];
            count++;
        }

        pie.emplace_back(new TPie(("pie_"+label).c_str()," ",back_n, &values[0], &colors[0]));
        pie.back()->SetRadius( pie.back()->GetRadius() * 0.8 );
        for(int iEntry = 0; iEntry < pie.back()->GetEntries(); ++iEntry) {
            pie.back()->SetEntryLabel(iEntry,"");
        }
        pie.back()->Draw();
        tex.DrawLatex(0.1,0.85,label.c_str());
    }

    c.cd();

    //
    // Adding the legend in the top panel
    //
    pTop.cd();
    std::unique_ptr<TLegend> leg(nullptr);
    if(TRExFitter::OPTION["FourTopStyle"]>0 || TRExFitter::OPTION["TRExbbStyle"]>0){
        leg = std::make_unique<TLegend>(0.5,0.1,0.95,0.90);
        leg->SetNColumns(3);
    }
    else{
        leg = std::make_unique<TLegend>(0.7,0.1,0.95,0.90);
        if(map_for_legend.size()>4){
            leg->SetNColumns(2);
        }
    }

    leg->SetLineStyle(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetBorderSize(0);
    leg->SetTextFont( gStyle->GetTextFont() );
    leg->SetTextSize( gStyle->GetTextSize() );

    std::vector<std::string> legVec;
    for ( const std::pair < std::string, int > legend_entry : map_for_legend ) {
        legVec.push_back(legend_entry.first);
    }
    std::vector<std::unique_ptr<TH1D> > dummy;
    for(int i_leg=legVec.size()-1;i_leg>=0;i_leg--){
        dummy.emplace_back(new TH1D(("legend_entry_" + legVec[i_leg]).c_str(), "",1,0,1));
        dummy.back()->SetFillColor(map_for_legend[legVec[i_leg]]);
        dummy.back()->SetLineColor(kBlack);
        dummy.back()->SetLineWidth(1);
        leg->AddEntry(dummy.back().get(),legVec[i_leg].c_str(),"f");
    }
    leg->Draw();

    //
    // Stores the pie chart in the desired format
    //
    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fName+"/PieChart" + fSuffix + (isPostFit ? "_postFit" : "") + "."+format).c_str());
    }
}

//__________________________________________________________________________________
// called before w in case of CustomAsimov
void TRExFit::CreateCustomAsimov() const {
    WriteDebugStatus("TRExFit::CreateCustomAsimov", "Running CreateCustomAsimov");
    // get a list of all CustomAsimov to create
    std::vector<std::string> customAsimovList;
    for(int i_smp=0; i_smp<fNSamples; ++i_smp) {
        if(fSamples[i_smp]->fAsimovReplacementFor.first!="" && Common::FindInStringVector(customAsimovList,fSamples[i_smp]->fAsimovReplacementFor.first)<0) {
            customAsimovList.push_back(fSamples[i_smp]->fAsimovReplacementFor.first);
        }
    }
    //
    // fill a different CustomAsimov data-set for each element in the list
    for(const auto& customAsimov : customAsimovList){
        WriteDebugStatus("TRExFit::CreateCustomAsimov", "CustomAsimov: " + customAsimov);
        std::shared_ptr<Sample> ca = GetSample("customAsimov_"+customAsimov);
        // create a new data sample taking the nominal S and B
        for(int i_ch=0; i_ch<fNRegions; ++i_ch) {
            Region *reg = fRegions[i_ch];
            // Now we need to clone a histogram, but need to find one that is valid
            std::shared_ptr<SampleHist> sample_hist = reg->fData;
            if (!sample_hist) {
                // try to clone signal
                for (const auto& isig : reg->fSig) {
                    if (isig != nullptr) {
                        sample_hist = isig;
                        break;
                    }
                }
            }
            if (!sample_hist) {
                // try to clone background
                for (const auto& ibkg : reg->fBkg) {
                    if (ibkg != nullptr) {
                        sample_hist = ibkg;
                        break;
                    }
                }
            }

            if (!sample_hist) {
                WriteErrorStatus("TRExFit::CreateCustomAsimov","Cannot copy a valid sample hist!");
                exit(EXIT_FAILURE);
            }

            std::shared_ptr<SampleHist> cash = reg->SetSampleHist(ca.get(),static_cast<TH1*>(sample_hist->fHist->Clone()));
            cash->fHist_orig->SetName( Form("%s_orig",cash->fHist->GetName()) ); // fix the name
            cash->fHist->Scale(0.);
            //
            std::vector<std::string> smpToExclude;
            for(int i_smp=0; i_smp<fNSamples; ++i_smp){
                std::shared_ptr<SampleHist> h = reg->GetSampleHist(fSamples[i_smp]->fName);
                if(!h) continue;
                if(h->fSample->fType==Sample::DATA) continue;
                if(h->fSample->fType==Sample::GHOST) {
                    if(h->fSample->fAsimovReplacementFor.first!=customAsimov) continue;
                    if(h->fSample->fAsimovReplacementFor.second!="" ) smpToExclude.push_back(h->fSample->fAsimovReplacementFor.second);
                }
                if( Common::FindInStringVector( smpToExclude,fSamples[i_smp]->fName )>=0 ) continue;
                //
                // bug-fix: change normalisation factors to nominal value!
                double factor = 1.;
                for(const auto& norm : fSamples[i_smp]->fNormFactors){
                    if (std::find(norm->fRegions.begin(), norm->fRegions.end(), fRegions[i_ch]->fName) != norm->fRegions.end()) {
                        WriteDebugStatus("TRExFit::CreateCustomAsimov", "setting norm factor to " + std::to_string(norm->fNominal));
                        factor *= norm->fNominal;
                    }
                }
                //
                cash->fHist->Add(h->fHist.get(),factor);
            }
            cash->fHist->Sumw2(false);
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::UnfoldingAlternativeAsimov() {
    if (fFitType != TRExFit::FitType::UNFOLDING) return;
    if (fAlternativeAsimovTruthSample == "")     return;

    WriteInfoStatus("TRExFit::UnfoldingAlternativeAsimov", "Replacing data with alternative asimov");

    // loop over regions
    for (auto& ireg : fRegions) {
        std::shared_ptr<SampleHist> sh = ireg->fData;

        // Try getting signal
        if (!sh) {
            for (const auto& isig : ireg->fSig) {
                if (!isig) {
                    sh = isig;
                    break;
                }
            }
        }

        if (!sh) {
            WriteErrorStatus("TRExFit::UnfoldingAlternativeAsimov", "No data or signal found, this should not happen!");
            exit(EXIT_FAILURE);
        }

        // now replace fData
        TH1* hist = static_cast<TH1*>(sh->fHist->Clone());
        hist->Reset();

        // Add new signal
        std::shared_ptr<Sample> newAsimov = GetSample("AlternativeSignal_"+ireg->fName);
        if (!newAsimov) {
            WriteErrorStatus("TRExFit::UnfoldingAlternativeAsimov", "Cannot read the new asimov sample");
            exit(EXIT_FAILURE);
        }
        std::shared_ptr<SampleHist> newsh = ireg->GetSampleHist(newAsimov->fName);
        if (!newsh) {
            WriteErrorStatus("TRExFit::UnfoldingAlternativeAsimov", "Cannot read the new asimov SampleHist");
            exit(EXIT_FAILURE);
        }
        hist->Add(newsh->fHist.get());

        // add bkgs bkg
        for (const auto& isample : fSamples) {
            std::shared_ptr<SampleHist> sampleHist = ireg->GetSampleHist(isample->fName);
            if (!sampleHist) continue;
            if (sampleHist->fSample->fType != Sample::BACKGROUND) continue;
            if(Common::FindInStringVector(isample->fRegions,ireg->fName) <0) continue;
            hist->Add(sampleHist->fHist.get());
        }

        hist->Sumw2(false);
        ireg->fData->fHist.reset(static_cast<TH1*>(hist->Clone()));
    }
}

//__________________________________________________________________________________
// turn to RooStat::HistFactory
void TRExFit::ToRooStat(bool makeWorkspace, bool exportOnly) const {

    WriteInfoStatus("TRExFit::ToRooStat", "-------------------------------------------");
    WriteInfoStatus("TRExFit::ToRooStat", "Exporting to RooStats...");

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);

    RooStats::HistFactory::Measurement meas((fInputName+fSuffix).c_str(), (fInputName+fSuffix).c_str());
    if(fBootstrap!="" && fBootstrapIdx>=0) {
        meas.SetOutputFilePrefix((fName+"/RooStats/"+fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName).c_str());
    } else {
        meas.SetOutputFilePrefix((fName+"/RooStats/"+fInputName).c_str());
    }
    meas.SetExportOnly(exportOnly);
    meas.SetPOI(fPOI.c_str());
    meas.SetLumi(fLumiScale);
    if(fLumiErr==0){
        meas.AddConstantParam("Lumi");
        meas.SetLumiRelErr(0.1);
    } else {
        meas.SetLumiRelErr(fLumiErr);
    }

    for(int i_ch=0;i_ch<fNRegions;i_ch++) {

        if(fRegions[i_ch]->fRegionType==Region::VALIDATION) continue;

        WriteDebugStatus("TRExFit::ToRooStat", "Adding Channel: " + fRegions[i_ch]->fName);
        const RooStats::HistFactory::Channel chan = OneChannelToRooStats(&meas, i_ch);

        meas.AddChannel(chan);
    }
    // Experimental: turn off constraints for given systematics
    for(const auto& isyst : fSystematics) {
        if(isyst->fIsFreeParameter) meas.AddUniformSyst(isyst->fName.c_str());
    }

    // morphing
    for(const TRExFit::TemplateWeight& itemp : fTemplateWeightVec){
        const std::string normName = "morph_"+itemp.name+"_"+Common::ReplaceString(std::to_string(itemp.value),"-","m");
        WriteDebugStatus("TRExFit::ToRooStat", "Morphing: normName: " + normName);
        meas.AddPreprocessFunction(normName, itemp.function, itemp.range);
    }
    for(const auto& nf : fNormFactors){
        if(nf->fExpression.first!=""){
            meas.AddPreprocessFunction(nf->fName,nf->fExpression.first,nf->fExpression.second);
        }
    }
    //
    if(fBootstrap!="" && fBootstrapIdx>=0) {
        meas.PrintXML((fName+"/RooStats/"+fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx)+"/").c_str());
    } else {
        meas.PrintXML((fName+"/RooStats/").c_str());
    }
    meas.CollectHistograms();
    meas.PrintTree();

    if(makeWorkspace) RooStats::HistFactory::MakeModelAndMeasurementFast(meas);

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
}

//__________________________________________________________________________________
//
RooStats::HistFactory::Channel TRExFit::OneChannelToRooStats(RooStats::HistFactory::Measurement* meas, const int i_ch) const {
    RooStats::HistFactory::Channel chan(fRegions[i_ch]->fName.c_str());

    //Suffix used for the regular bin transformed histogram
    static const std::string suffix_regularBinning("_regBin");

    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0; i_smp<fNSamples; ++i_smp){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    if(fCustomAsimov != "") {
        const std::string name = "customAsimov_"+fCustomAsimov;
        std::shared_ptr<SampleHist> cash = fRegions[i_ch]->GetSampleHist(name);
        if(cash==nullptr){
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            WriteWarningStatus("TRExFit::OneChannelToRooStats", "No Custom Asimov " + fCustomAsimov + " available. Taking regular Asimov.");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        } else{
            const std::string temp_string = cash->fHist->GetName();
            WriteDebugStatus("TRExFit::OneChannelToRooStats", "  Adding Custom-Asimov Data: " + temp_string);
            chan.SetData(cash->fHistoName+suffix_regularBinning, cash->fFileName);
        }
    } else if(hasData) {
        const std::string temp_string = fRegions[i_ch]->fData->fHist->GetName();
        WriteDebugStatus("TRExFit::OneChannelToRooStats", "  Adding Data: " + temp_string);
        chan.SetData(fRegions[i_ch]->fData->fHistoName+suffix_regularBinning, fRegions[i_ch]->fData->fFileName);
    } else {
        chan.SetData("", "");
    }

    // fStatErrCons is upper case after config reading if the MCstatThreshold option is used, otherwise it defaults to "Poisson"
    // HistFactory expects the constraint not in all uppercase, but in form "Poisson"/"Gaussian" instead
    if(fStatErrCons=="Poisson" || fStatErrCons=="POISSON") chan.SetStatErrorConfig(fStatErrThres, "Poisson");
    else if(fStatErrCons=="GAUSSIAN")                      chan.SetStatErrorConfig(fStatErrThres, "Gaussian");

    for(int i_smp=0; i_smp<fNSamples; ++i_smp) {
        std::shared_ptr<SampleHist> h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
        if (!h) continue;
        if (h->fSample->fType == Sample::DATA) continue;
        if (h->fSample->fType == Sample::GHOST) continue;

        WriteDebugStatus("TRExFit::OneChannelToRooStats", "  Adding Sample: " + fSamples[i_smp]->fName);
        const RooStats::HistFactory::Sample sample = OneSampleToRooStats(meas, h.get(), i_ch, i_smp);;
        chan.AddSample(sample);
    }

    return chan;
}

//__________________________________________________________________________________
//
RooStats::HistFactory::Sample TRExFit::OneSampleToRooStats(RooStats::HistFactory::Measurement* meas,
                                                           const SampleHist* h,
                                                           const int i_ch,
                                                           const int i_smp) const {
    RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());

    static const std::string suffix_regularBinning("_regBin");

    if(fUseStatErr && fSamples[i_smp]->fUseMCStat) sample.ActivateStatError();
    sample.SetHistoName(h->fHistoName+suffix_regularBinning);
    sample.SetInputFile(h->fFileName);
    sample.SetNormalizeByTheory(fSamples[i_smp]->fNormalizedByTheory);
    // norm factors
    for(const auto& inorm : h->fSample->fNormFactors) {
        WriteDebugStatus("TRExFit::OneSampleToRooStats", "    Adding NormFactor: " + inorm->fName + ", " + std::to_string(inorm->fNominal));
        sample.AddNormFactor(inorm->fName,
                             inorm->fNominal,
                             inorm->fMin,
                             inorm->fMax);
        if (inorm->fConst) meas->AddConstantParam(inorm->fName);
        if (fStatOnly && fFixNPforStatOnlyFit && inorm->fName!=fPOI) {
            meas->AddConstantParam(inorm->fName);
        }
    }

    // shape factors
    for(const auto& ishape : h->fSample->fShapeFactors) {
        WriteDebugStatus("TRExFit::OneSampleToRooStats", "    Adding ShapeFactor: " + ishape->fName + ", " + std::to_string(ishape->fNominal));
        sample.AddShapeFactor(ishape->fName);
        if (ishape->fConst
            || (fStatOnly && fFixNPforStatOnlyFit && ishape->fName!=fPOI) ) {
            for(int i_bin=0; i_bin < ishape->fNbins; ++i_bin) {
                meas->AddConstantParam( "gamma_" + ishape->fName + "_bin_" + std::to_string(i_bin) );
            }
        }
    }
    // systematics
    if(fStatOnly) {
        sample.AddOverallSys( "Dummy",1,1 );
        return sample;
    }
    for(std::size_t i_syst = 0; i_syst < h->fSyst.size(); ++i_syst) {
        // add normalization part
        WriteDebugStatus("TRExFit::OneSampleToRooStats", "    Adding Systematic: " + h->fSyst[i_syst]->fName);
        if ( h->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
            std::string npName = "shape_" + h->fSyst[i_syst]->fSystematic->fNuisanceParameter+"_";
            std::string regionName = fRegions[i_ch]->fName;
            if(h->fSyst[i_syst]->fSystematic->fNuisanceParameter.find("stat_")!=std::string::npos) {
                // see if there are regions to correlate with others
                for(auto& set : h->fSample->fCorrelateGammasInRegions){
                    for(std::size_t i_reg = 0; i_reg < set.size(); ++i_reg){
                        if(i_reg != 0 && regionName == set[i_reg]) {
                            regionName = set[0];
                            break;
                        }
                    }
                }
                // eventually correlate MC stat with other samples
                if(h->fSample->fCorrelateGammasWithSample != ""){
                    npName = "shape_stat_"+h->fSample->fCorrelateGammasWithSample+"_";
                }
            }
            npName += regionName;
            sample.AddShapeSys( npName,
                                RooStats::HistFactory::Constraint::Poisson,
                                h->fSyst[i_syst]->fHistoNameUp+"_Var"+suffix_regularBinning,
                                h->fSyst[i_syst]->fFileNameUp, "");
        } else {
            if ( !h->fSyst[i_syst]->fSystematic->fIsShapeOnly &&
                !h->fSyst[i_syst]->fNormPruned &&
                !h->fSyst[i_syst]->fBadNorm
              ) {
                sample.AddOverallSys( h->fSyst[i_syst]->fSystematic->fNuisanceParameter,
                                      1+h->fSyst[i_syst]->fNormDown,
                                      1+h->fSyst[i_syst]->fNormUp);
            }
            // eventually add shape part
            if ( h->fSyst[i_syst]->fIsShape &&
                !h->fSyst[i_syst]->fSystematic->fIsNormOnly &&
                !h->fSyst[i_syst]->fShapePruned &&
                !h->fSyst[i_syst]->fBadShape
              ){
                sample.AddHistoSys( h->fSyst[i_syst]->fSystematic->fNuisanceParameter,
                                    h->fSyst[i_syst]->fHistoNameShapeDown+suffix_regularBinning,
                                    h->fSyst[i_syst]->fFileNameShapeDown, "",
                                    h->fSyst[i_syst]->fHistoNameShapeUp+suffix_regularBinning,
                                    h->fSyst[i_syst]->fFileNameShapeUp,   ""  );
            }
        }
    }

    return sample;
}

//__________________________________________________________________________________
//
void TRExFit::SystPruning() const {
    WriteInfoStatus("TRExFit::SystPruning", "------------------------------------------------------");
    WriteInfoStatus("TRExFit::SystPruning", "Apply Systematics Pruning ...");
    if (fPruningShapeOption == PruningUtil::SHAPEOPTION::KSTEST) {
        WriteInfoStatus("TRExFit::SystPruning", "Will run KS test to determine the shape pruning. This is slow compared to the default option (MAXBIN). Patience young padawan.");
    }
    if(fSystematics.size()==0 || fStatOnly){
        WriteInfoStatus("TRExFit::SystPruning", "No systematics => No Pruning applied.");
        return;
    }

    PruningUtil pu{};
    pu.SetShapeOption(fPruningShapeOption);
    pu.SetStrategy(static_cast<int>(fPruningType));
    pu.SetThresholdNorm(fThresholdSystPruning_Normalisation);
    pu.SetThresholdShape(fThresholdSystPruning_Shape);
    pu.SetThresholdIsLarge(fThresholdSystLarge);
    pu.SetRemoveLargeSyst(fRemoveLargeSyst);
    pu.SetRemoveSystOnEmptySample(fRemoveSystOnEmptySample);

    for(auto& reg : fRegions){
        // if want to skip validation regions from pruning, add a condition here
        reg->SystPruning(&pu);
    }

    // Draw plot with normalisation effect of each systematic
    if (fDoSystNormalizationPlots) {
        DrawSystematicNormalisationSummary();
    }

    // it also writes the txt file actually
    DrawPruningPlot();
}

//__________________________________________________________________________________
//
void TRExFit::DrawSystematicNormalisationSummary() const {

    // First get the dimensions right
    const std::vector<std::string> uniqueSysts = GetUniqueSystNamesWithoutGamma();
    const std::vector<Region*> regions = GetNonValidationRegions();
    const std::vector<std::shared_ptr<Sample> > samples = GetNonDataNonGhostSamples();

    const int xBins = regions.size()*samples.size();

    TH2D histo("", "", xBins, 0, xBins, 2*uniqueSysts.size(), 0, 2*uniqueSysts.size());
    histo.GetXaxis()->SetTickSize(0);
    histo.GetYaxis()->SetTickSize(0);
    for (std::size_t ireg = 0; ireg < regions.size(); ++ireg) {
        for (std::size_t ibin = 0; ibin < samples.size(); ++ibin) {
            histo.GetXaxis()->SetBinLabel(ireg*samples.size()+ibin+1,
                                          (regions.at(ireg)->fName+"_"+samples.at(ibin)->fName).c_str());
        }
    }
    for (std::size_t isyst = 0; isyst < uniqueSysts.size(); ++isyst) {
        histo.GetYaxis()->SetBinLabel(2*isyst+1, (uniqueSysts.at(isyst)+"_Up").c_str());
        histo.GetYaxis()->SetBinLabel(2*isyst+2, (uniqueSysts.at(isyst)+"_Dn").c_str());
    }
    histo.GetXaxis()->LabelsOption("v");

    // Fill the histograms
    for (std::size_t ireg = 0; ireg < regions.size(); ++ireg) {
        for (std::size_t isample = 0; isample < samples.size(); ++isample) {
            std::shared_ptr<SampleHist> sh = regions.at(ireg)->GetSampleHist(samples.at(isample)->fName);
            for (std::size_t isyst = 0; isyst < uniqueSysts.size(); ++isyst) {
                if (!sh) {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                    continue;
                }

                // actually calcualte the normalisation effect
                const TH1* nominal = sh->fHist.get();
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic(uniqueSysts.at(isyst));
                if (!syh) {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                    continue;
                }

                if (!nominal) {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                    continue;
                }

                const TH1* up   = syh->fHistUp.get();
                const TH1* down = syh->fHistDown.get();

                if (up) {
                    const double norm = 100*(Common::EffIntegral(up) - Common::EffIntegral(nominal))/Common::EffIntegral(nominal);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        norm);
                } else {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                }

                if (down) {
                    const double norm = 100*(Common::EffIntegral(down) - Common::EffIntegral(nominal))/Common::EffIntegral(nominal);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        norm);
                } else {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                }
            }
        }
    }

    static const int upSize = 50;
    static const int loSize = 150;
    static const int leftSize = 250;
    static const int separation = 10;
    const int regionSize = 20*samples.size();
    const int mainHeight = 2*uniqueSysts.size()*20;
    const int mainWidth = regions.size()*(regionSize+separation);

    TCanvas c("","", leftSize+mainWidth, upSize+mainHeight+loSize);
    c.SetTopMargin(2./(2.*uniqueSysts.size()));
    c.SetBottomMargin(200./(2600.+2.*uniqueSysts.size()));
    c.SetLeftMargin(0.3);
    c.SetRightMargin(0.0);
    gStyle->SetPalette(87);
    c.SetGrid();

    histo.SetMarkerSize(450);
    histo.GetXaxis()->SetLabelOffset(0.3*histo.GetXaxis()->GetLabelOffset());
    gStyle->SetPaintTextFormat(".1f");
    histo.Draw("col TEXT");
    c.RedrawAxis("g");

    for(const auto& iformat : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fName+"/NormalisationPlot"+fSuffix+"."+iformat).c_str());
    }
}

//__________________________________________________________________________________
//
void TRExFit::DrawPruningPlot() const{
    //
    std::ofstream out;
    out.open((fName+"/PruningText.txt").c_str());
    out << "-------///////                 ///////-------" << std::endl ;
    out << "-------/////// IN PRUNING PLOT ///////-------" << std::endl ;
    out << "-------///////                 ///////-------" << std::endl ;
    //
    std::vector< std::unique_ptr<TH2F> > histPrun;
    std::vector< std::unique_ptr<TH2F> > histPrun_toSave;
    int iReg = 0;
    int nSmp = 0;
    // make a list of non-data, non-ghost samples
    std::vector<std::shared_ptr<Sample> > samplesVec;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        samplesVec.emplace_back(fSamples[i_smp]);
        nSmp++;
    }
    // make a list of non-gamma systematics only
    std::vector<std::shared_ptr<Systematic> > nonGammaSystematics;
    const std::vector<std::string>& uniqueSyst = GetUniqueSystNamesWithoutGamma();
    for(const auto& isyst : fSystematics) {
        if(isyst->fType == Systematic::SHAPE) continue;

        nonGammaSystematics.push_back(isyst);
    }
    const size_t NnonGammaSyst = nonGammaSystematics.size();
    if(NnonGammaSyst==0){
        WriteInfoStatus("TRExFit::DrawPruningPlot", "No non-gamma systematics found => No Pruning plot generated.");
        return;
    }
    //
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        if(!fValidationPruning && fRegions[i_reg]->fRegionType==Region::VALIDATION) continue;

        out << "In Region : " << fRegions[i_reg]->fName << std::endl ;
        histPrun.emplace_back(new TH2F (Form("h_prun_%s", fRegions[i_reg]->fName.c_str()  ),fRegions[i_reg]->fShortLabel.c_str(),nSmp,0,nSmp, uniqueSyst.size(),0,uniqueSyst.size()));
        histPrun.back()->SetDirectory(0);

        for(int i_smp=0;i_smp<nSmp;i_smp++){
            out << " -> In Sample : " << samplesVec[i_smp]->fName << std::endl;

            for(std::size_t uniqueIndex = 0; uniqueIndex < uniqueSyst.size(); ++uniqueIndex){
               histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -1 );
            }

            std::shared_ptr<SampleHist> sh = fRegions[i_reg]->GetSampleHist(samplesVec[i_smp]->fName);
            if (sh == nullptr) continue;

            for(size_t i_syst=0;i_syst<NnonGammaSyst;i_syst++){
                // find the corresponding index of unique syst
                auto it = std::find(uniqueSyst.begin(), uniqueSyst.end(), nonGammaSystematics.at(i_syst)->fName);
                const std::size_t uniqueIndex = std::distance(uniqueSyst.begin(), it);
                out << " --->>  " << nonGammaSystematics[i_syst]->fName << "     " ;
                if( (Common::FindInStringVector(nonGammaSystematics[i_syst]->fSamples,samplesVec[i_smp]->fName)>=0 || nonGammaSystematics[i_syst]->fSamples[0] == "all")
                    && sh->HasSyst(nonGammaSystematics[i_syst]->fName)
                ){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(nonGammaSystematics[i_syst]->fName);
                    histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 0 );
                    const bool forgeDropShape = nonGammaSystematics[i_syst]->fIsNormOnly;
                    const bool forgeDropNorm  = nonGammaSystematics[i_syst]->fIsShapeOnly;
                    //
                    if(syh->fShapePruned && syh->fNormPruned) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 3 );
                    else if(syh->fShapePruned || forgeDropShape) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 1 );
                    else if(syh->fNormPruned || forgeDropNorm) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 2 );
                    //
                    if(syh->fBadShape && syh->fBadNorm) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -4 );
                    else if(syh->fBadShape) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -3 );
                    else if(syh->fBadNorm) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -2 );
                    //
                }
                if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -1 ) out << " is not present" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 0 ) out << " is kept" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 1 ) out << " is norm only" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 2 ) out << " is shape only" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 3 ) out << " is dropped" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -2 ) out << " has bad norm" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -3 ) out << " has bad shape" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -4 ) out << " is bad" << std::endl;
            }
        }
        //
        histPrun_toSave.emplace_back(std::move(std::unique_ptr<TH2F>(static_cast<TH2F*>(histPrun[iReg]->Clone(Form("%s_toSave",histPrun[iReg]->GetName()))))) );
        histPrun_toSave[iReg]->SetDirectory(0);
        //
        iReg++;
    }
    //
    // draw the histograms
    int upSize = 50;
    int loSize = 150;
    int mainHeight = uniqueSyst.size()*20;
    int leftSize = 250;
    int regionSize = 20*nSmp;
    int separation = 10;
    int mainWidth = iReg*(regionSize+separation);
    //
    TCanvas c("c_pruning","Canvas - Pruning",leftSize+mainWidth,upSize+mainHeight+loSize);
    std::vector<Int_t> colors = {kBlack,6,kBlue, kGray, 8, kYellow, kOrange-3, kRed}; // #colors >= #levels - 1
    gStyle->SetPalette(colors.size(), &colors[0]);
    TPad pUp("pUp","Pad High",0,(1.*loSize+mainHeight)/(upSize+mainHeight+loSize),1,1);
    pUp.Draw();
    c.cd();
    std::vector<std::unique_ptr<TPad> > pReg(100);
    for(std::size_t i_reg=0;i_reg<histPrun.size();i_reg++){
        c.cd();
        if(i_reg==0){
            pReg[i_reg] = std::move(std::unique_ptr<TPad> (new TPad(Form("pReg[%zu]",i_reg),"Pad Region",
                                  0,   0,
                                  (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize) )));
            pReg[i_reg]->SetLeftMargin( (1.*leftSize) / (1.*leftSize+regionSize) );
        }
        else{
            pReg[i_reg] = std::move(std::unique_ptr<TPad> (new TPad(Form("pReg[%zu]",i_reg),"Pad Region",
                                  (leftSize+1.*i_reg*(regionSize+separation))           /(leftSize+mainWidth),   0,
                                  (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize) )));
            pReg[i_reg]->SetLeftMargin(0);
        }
        pReg[i_reg]->SetBottomMargin( (1.*loSize) / (1.*loSize+mainHeight) );
        pReg[i_reg]->Draw();
        pReg[i_reg]->cd();
        gPad->SetGridy();
        for(int i_bin=1;i_bin<=histPrun[i_reg]->GetNbinsX();i_bin++){
            histPrun[i_reg]       ->GetXaxis()->SetBinLabel(i_bin,samplesVec[i_bin-1]->fTitle.c_str());
            histPrun_toSave[i_reg]->GetXaxis()->SetBinLabel(i_bin,samplesVec[i_bin-1]->fName.c_str());
        }
        for(int i_bin=1;i_bin<=histPrun[i_reg]->GetNbinsY();i_bin++){
            if(i_reg==0) {
                histPrun[i_reg]       ->GetYaxis()->SetBinLabel(i_bin,TRExFitter::SYSTMAP[uniqueSyst[i_bin-1]].c_str());
            }
            else {
                histPrun[i_reg]->GetYaxis()->SetBinLabel(i_bin,"");
            }
            histPrun_toSave[i_reg]->GetYaxis()->SetBinLabel(i_bin,uniqueSyst[i_bin-1].c_str());
        }
        histPrun[i_reg]->Draw("COL");
        histPrun[i_reg]->GetYaxis()->SetLabelOffset(0.03);
        gPad->SetTopMargin(0);
        gPad->SetRightMargin(0);
        histPrun[i_reg]->GetXaxis()->LabelsOption("v");
        histPrun[i_reg]->GetXaxis()->SetLabelSize( histPrun[i_reg]->GetXaxis()->GetLabelSize()*0.75 );
        histPrun[i_reg]->GetYaxis()->SetLabelSize( histPrun[i_reg]->GetYaxis()->GetLabelSize()*0.75 );
        gPad->SetTickx(0);
        gPad->SetTicky(0);
        histPrun[i_reg]->SetMinimum(-4);
        histPrun[i_reg]->SetMaximum( 3.1);
        histPrun[i_reg]->GetYaxis()->SetTickLength(0);
        histPrun[i_reg]->GetXaxis()->SetTickLength(0);
        gPad->SetGrid();
        //
        pUp.cd();
        myText((leftSize+1.*i_reg*(regionSize+separation))/(leftSize+mainWidth),0.1 ,1,histPrun[i_reg]->GetTitle());
    }
    c.cd();
    TPad pLo("pLo","Pad Low",0,0,(1.*leftSize)/(leftSize+mainWidth),(1.*loSize)/(upSize+mainHeight+loSize));
    pLo.Draw();
    //
    c.cd();
    pUp.cd();
    myText(0.01,0.5,1,fLabel.c_str());
    //
    pLo.cd();
    TLegend leg(0.005,0,0.95,0.95);
    TH1D hGray   ("hGray"  ,"hGray"  ,1,0,1);    hGray.SetFillColor(kGray);         hGray.SetLineWidth(0);
    TH1D hYellow ("hYellow","hYellow",1,0,1);    hYellow.SetFillColor(kYellow);     hYellow.SetLineWidth(0);
    TH1D hOrange ("hOrange","hOrange",1,0,1);    hOrange.SetFillColor(kOrange-3);   hOrange.SetLineWidth(0);
    TH1D hRed    ("hRed"   ,"hRed"   ,1,0,1);    hRed.SetFillColor(kRed);           hRed.SetLineWidth(0);
    TH1D hGreen  ("hGreen" ,"hGree"  ,1,0,1);    hGreen.SetFillColor(8);            hGreen.SetLineWidth(0);
    TH1D hBlue   ("hBlue"  ,"hBlue"  ,1,0,1);    hBlue.SetFillColor(kBlue);         hBlue.SetLineWidth(0);
    TH1D hPurple ("hPurple","hPurple",1,0,1);    hPurple.SetFillColor(6);           hPurple.SetLineWidth(0);
    TH1D hBlack  ("hBlack" ,"hBlack" ,1,0,1);    hBlack.SetFillColor(kBlack);       hBlack.SetLineWidth(0);
    std::string sysLarg="Dropped as >"+std::to_string((int)(fThresholdSystLarge*100))+"%";
    leg.SetBorderSize(0);
    leg.SetMargin(0.1);
    leg.SetFillStyle(0);
    leg.AddEntry(&hGray,"Not present","f");
    leg.AddEntry(&hGreen,"Kept","f");
    leg.AddEntry(&hYellow, "Shape dropped","f");
    leg.AddEntry(&hOrange, "Norm. dropped","f");
    leg.AddEntry(&hRed, "Dropped","f");
    if (fThresholdSystLarge > -1) {
        leg.AddEntry(&hBlue  , sysLarg.c_str() ,"f");
        leg.AddEntry(&hPurple, "Bad shape" ,"f");
        leg.AddEntry(&hBlack , "Bad shape & norm." ,"f");
    }
    leg.SetTextSize(0.85*gStyle->GetTextSize());
    leg.Draw();
    //
    for(const auto& iformat : TRExFitter::IMAGEFORMAT) {
        c.SaveAs( (fName+"/Pruning"+fSuffix+"."+iformat).c_str() );
    }

    //
    // Save prunign hist for future usage
    std::unique_ptr<TFile> filePrun = nullptr;
    // - checking if Pruning.root exists
    // if yes
    if(!gSystem->AccessPathName( (fName+"/Pruning.root").c_str() )){
        // ...
        filePrun.reset(TFile::Open( (fName+"/Pruning.root").c_str() ));
    }
    else{
        filePrun.reset(TFile::Open( (fName+"/Pruning.root").c_str(),"RECREATE" ));
        for(std::size_t i_reg=0;i_reg<histPrun.size();i_reg++){
            histPrun_toSave[i_reg]->Write("",TObject::kOverwrite);
        }
    }
    if (filePrun != nullptr){
        filePrun->Close();
    }
}

//__________________________________________________________________________________
//
void TRExFit::Fit(bool isLHscanOnly){

    std::unique_ptr<RooDataSet> data(nullptr);
    std::unique_ptr<RooWorkspace> ws(nullptr);
    
    //
    // Read NPvalues from fit-result file
    //
    if(fFitNPValuesFromFitResults!=""){
        WriteInfoStatus("TRExFit::Fit","Setting NP values for Asimov data-set creation from fit results stored in file " + fFitNPValuesFromFitResults + "...");
        fFitNPValues = NPValuesFromFitResults(fFitNPValuesFromFitResults);
    }

    //
    // If fDoNonProfileFit => set stat-only
    //
    if(fDoNonProfileFit){
        WriteInfoStatus("TRExFit::Fit","In non-profile mode => Setting to stat-only.");
        fStatOnly = true;
    }

    //
    // If there's a workspace specified, go on with simple fit, without looking for separate workspaces per region
    //
    if(fWorkspaceFileName!=""){
        WriteInfoStatus("TRExFit::Fit","");
        WriteInfoStatus("TRExFit::Fit","-------------------------------------------");
        WriteInfoStatus("TRExFit::Fit","Performing nominal fit on pre-specified workspace...");
        TFile *rootFile = TFile::Open(fWorkspaceFileName.c_str(),"read");
        ws = std::unique_ptr<RooWorkspace>(static_cast<RooWorkspace*>(rootFile->Get("combined")));
        if(!ws){
            WriteErrorStatus("TRExFit::Fit", "The workspace (\"combined\") cannot be found in file " + fWorkspaceFileName + ". Please check !");
            exit(EXIT_FAILURE);
        }
        if(!fFitIsBlind) data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("obsData")));
        else             data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("asimovData")));
    }
    //
    // Otherwise go on with normal fit
    //
    else{
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std:: string > regionsToFit = ListRegionsToFit();
        std::map < std::string, int > regionDataType = MapRegionDataTypes(regionsToFit,fFitIsBlind);
        std::map < std::string, double > npValues;
        //
        // flag if mixed Data / Asimov fit required
        // if mixed fit, perform a first fit on the regions with data only
        bool isMixedFit = DoingMixedFitting();
        if(isMixedFit && !fFitIsBlind){
            if (!isLHscanOnly){
                WriteInfoStatus("TRExFit::Fit","");
                WriteInfoStatus("TRExFit::Fit","-------------------------------------------");
                WriteInfoStatus("TRExFit::Fit","Performing fit on regions with DataType = DATA to get NPs to inject in Asimov...");
            }
            std::vector < std:: string > regionsForDataFit = ListRegionsToFit(Region::REALDATA);
            std::map < std::string, int > regionForDataFitDataType = MapRegionDataTypes(regionsForDataFit);
            //
            // Creates a combined workspace with the regions to be used *in the data fit*
            //
            WriteInfoStatus("TRExFit::Fit","Creating ws for regions with real data only...");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
            std::unique_ptr<RooWorkspace> ws_forFit = PerformWorkspaceCombination( regionsForDataFit );
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            if (!ws_forFit){
                WriteErrorStatus("TRExFit::Fit","Cannot retrieve the workspace, exiting!");
                exit(EXIT_FAILURE);
            }
            //
            // Calls the PerformFit() function to actually do the fit
            //
            WriteInfoStatus("TRExFit::Fit","Performing a B-only fit in regions with real data only...");
            npValues = PerformFit( ws_forFit.get(), data.get(), FitType::BONLY, false, TRExFitter::DEBUGLEVEL);
            WriteInfoStatus("TRExFit::Fit","Now will use the fit results to create the Asimov in the regions without real data!");
        }
        else{
            npValues = fFitNPValues;
        }
        //
        if (!isLHscanOnly){
            WriteInfoStatus("TRExFit::Fit","");
            WriteInfoStatus("TRExFit::Fit","-------------------------------------------");
            WriteInfoStatus("TRExFit::Fit","Performing nominal fit...");
        }
        //
        // Create the final asimov dataset for fit
        //
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        ws = PerformWorkspaceCombination( regionsToFit );
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        if (!ws){
            WriteErrorStatus("TRExFit::Fit","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npValues, npValues.find(fPOI)==npValues.end() ? fFitPOIAsimov : npValues[fPOI] ));
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    }
    
    //
    // Set the global observables to the specified NPValues if the corresponding flag is TRUE
    //
    if(fInjectGlobalObservables && !fFitNPValues.empty()) {
        FitUtils::InjectGlobalObservables(ws.get(), fFitNPValues);
    }
    
    //
    // Calls the PerformFit() function to actually do the fit
    //
    if (!isLHscanOnly) PerformFit( ws.get(), data.get(), fFitType, true, TRExFitter::DEBUGLEVEL);

    //
    // Toys
    //
    if(fFitToys>0 && !isLHscanOnly){
        RunToys();
    }

    //
    // Fit result on Asimov with shifted systematics
    //
    if(fDoNonProfileFit && !isLHscanOnly){
        //
        std::map<std::string,double> systGroups;
        std::vector<std::string> systGroupNames;
        std::map < std::string, double > npValues;
        //
        WriteInfoStatus("TRExFit::Fit","");
        WriteInfoStatus("TRExFit::Fit","-------------------------------------------");
        WriteInfoStatus("TRExFit::Fit","Scan of systematics for non-profile fit...");
        std::ofstream out;
        std::ofstream tex;
        std::ofstream out2;
        std::ofstream tex2;
        std::vector < std:: string > regionsToFit;
        if(fBootstrap!="" && fBootstrapIdx>=0){
            out.open((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts.txt").c_str());
            tex.open((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts.tex").c_str());
            tex2.open((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts_grouped.tex").c_str());
        }
        else{
            out.open((fName+"/Fits/"+fName+fSuffix+"_nonProfiledSysts.txt").c_str());
            tex.open((fName+"/Fits/"+fName+fSuffix+"_nonProfiledSysts.tex").c_str());
            tex2.open((fName+"/Fits/"+fName+fSuffix+"_nonProfiledSysts_grouped.tex").c_str());
        }
        std::map < std::string, int > regionDataType;
        for(auto reg : fRegions) regionDataType[reg->fName] = Region::ASIMOVDATA;
        for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
            if ( fFitRegion == CRONLY && fRegions[i_ch] -> fRegionType == Region::CONTROL )
                regionsToFit.push_back( fRegions[i_ch] -> fName );
            else if ( fFitRegion == CRSR && (fRegions[i_ch] -> fRegionType == Region::CONTROL || fRegions[i_ch] -> fRegionType == Region::SIGNAL) )
                regionsToFit.push_back( fRegions[i_ch] -> fName );
        }
        ws = PerformWorkspaceCombination( regionsToFit );
        if (!ws){
            WriteErrorStatus("TRExFit::Fit","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        // nominal fit on Asimov
        RooStats::ModelConfig *mc = (RooStats::ModelConfig*)ws -> obj("ModelConfig");
        ws->saveSnapshot("InitialStateModelGlob",   *mc->GetGlobalObservables());
        ws->saveSnapshot("InitialStateModelNuis",   *mc->GetNuisanceParameters());
        WriteInfoStatus("TRExFit::Fit","Fitting nominal Asimov...");
        data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, fFitNPValues, fFitPOIAsimov ));
        npValues = PerformFit( ws.get(), data.get(), fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
        double nominalPOIval = npValues[fPOI];
        RooRealVar* poiVar = (RooRealVar*) (& ws->allVars()[fPOI.c_str()]);
        double statUp = poiVar->getErrorHi();
        double statDo = poiVar->getErrorLo();
        // temporary switch off minos (not needed)
        std::vector<std::string> varMinosTmp = fVarNameMinos;
        fVarNameMinos.clear();
        // MC stat
        double MCstatUp = 0.;
        double MCstatDo = 0.;
        if(fUseStatErr){
            fGammasInStatOnly = true;
            WriteInfoStatus("TRExFit::Fit","MC stat...");
            std::map < std::string, double > npVal;
            for(auto reg : fRegions){
                TH1* hTot = nullptr;
                for(auto& sh : reg->fSampleHists){
                    if(sh->fSample->fSeparateGammas) continue;
                    if(sh->fSample->fType==Sample::GHOST) continue;
                    if(sh->fSample->fType==Sample::DATA) continue;
                    if(!sh->fSample->fUseMCStat) continue; // need to fix something for separate gammas
                    bool skip = false;
                    for(auto morphPar : fMorphParams){
                        // find nominal value of this morph parameter
                        double nfVal = 0.;
                        for(auto nf : fNormFactors){
                            if(nf->fName==morphPar){
                                nfVal = nf->fNominal;
                                break;
                            }
                        }
                        if(sh->fSample->fIsMorph[morphPar] && sh->fSample->fMorphValue[morphPar]!=nfVal) skip = true;
                    }
                    if(skip) continue;
                    WriteInfoStatus("TRExFit::Fit","  Including sample "+sh->fSample->fName);
                    if(hTot==nullptr && sh->fHist!=nullptr) hTot = static_cast<TH1*>(sh->fHist->Clone("h_tot"));
                    else if(            sh->fHist!=nullptr) hTot->Add(sh->fHist.get());
                }
                if(hTot==nullptr) continue;
                for(int i_bin=1;i_bin<=hTot->GetNbinsX();i_bin++){
                    double statErr = hTot->GetBinError(i_bin)/hTot->GetBinContent(i_bin);
                    std::string gammaName = "gamma_stat_"+reg->fName+"_bin_"+std::to_string(i_bin-1);
                    // up
                    npVal = fFitNPValues;
                    npVal[gammaName] = 1+statErr;
                    WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1+statErr));
                    data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npVal, fFitPOIAsimov ));
                    npVal[gammaName] = 1;
                    ws->loadSnapshot("InitialStateModelGlob");
                    ws->loadSnapshot("InitialStateModelNuis");
                    npValues = PerformFit( ws.get(), data.get(), fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                    MCstatUp = std::hypot(MCstatUp, (npValues[fPOI]-nominalPOIval));
                    // down
                    npVal = fFitNPValues;
                    npVal[gammaName] = 1-statErr;
                    WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1-statErr));
                    data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npVal, fFitPOIAsimov ));
                    npVal[gammaName] = 1;
                    ws->loadSnapshot("InitialStateModelGlob");
                    ws->loadSnapshot("InitialStateModelNuis");
                    npValues = PerformFit( ws.get(), data.get(), fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                    MCstatDo = -std::hypot(MCstatDo,(npValues[fPOI]-nominalPOIval));
                }
            }
            fGammasInStatOnly = false;
        }
        // MC-stat for specific samples (those using separate gammas)
        std::map<std::string,double> MCstatUpSample;
        std::map<std::string,double> MCstatDoSample;
        std::map<std::string,std::string> smpTexTitle;
        if(fUseStatErr){
            std::map < std::string, double > npVal;
            for(auto reg : fRegions){
                for(auto& sh : reg->fSampleHists){
                    TH1* hTot = nullptr;
                    if(sh->fSample->fType==Sample::GHOST) continue;
                    if(sh->fSample->fType==Sample::DATA) continue;
                    if(!sh->fSample->fSeparateGammas) continue;
                    bool skip = false;
                    for(auto morphPar : fMorphParams){
                        // find nominal value of this morph parameter
                        double nfVal = 0.;
                        for(auto nf : fNormFactors){
                            if(nf->fName==morphPar){
                                nfVal = nf->fNominal;
                                break;
                            }
                        }
                        if(sh->fSample->fIsMorph[morphPar] && sh->fSample->fMorphValue[morphPar]!=nfVal) skip = true;
                    }
                    if(skip) continue;
                    WriteInfoStatus("TRExFit::Fit","MC stat for sample "+sh->fSample->fName+"...");
                    hTot = (TH1*)sh->fHist->Clone("h_tot");
                    if(hTot==nullptr) continue;
                    for(int i_bin=1;i_bin<=hTot->GetNbinsX();i_bin++){
                        double statErr = hTot->GetBinError(i_bin)/hTot->GetBinContent(i_bin);
                        std::string gammaName = "gamma_shape_stat_"+sh->fSample->fName+"_"+reg->fName+"_bin_"+std::to_string(i_bin-1);
                        npVal = fFitNPValues;
                        npVal[gammaName] = 1+statErr;
                        WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1+statErr));
                        data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npVal, fFitPOIAsimov ));
                        npVal[gammaName] = 1;
                        ws->loadSnapshot("InitialStateModelGlob");
                        ws->loadSnapshot("InitialStateModelNuis");
                        npValues = PerformFit( ws.get(), data.get(), fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                        MCstatUpSample[sh->fSample->fName] = std::hypot(MCstatUpSample[sh->fSample->fName], (npValues[fPOI]-nominalPOIval));
                        npVal = fFitNPValues;
                        npVal[gammaName] = 1-statErr;
                        WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1-statErr));
                        data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npVal, fFitPOIAsimov ));
                        npVal[gammaName] = 1;
                        ws->loadSnapshot("InitialStateModelGlob");
                        ws->loadSnapshot("InitialStateModelNuis");
                        npValues = PerformFit( ws.get(), data.get(), fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                        MCstatDoSample[sh->fSample->fName] = -std::hypot(MCstatDoSample[sh->fSample->fName],(npValues[fPOI]-nominalPOIval));
                        smpTexTitle[sh->fSample->fName] = sh->fSample->fTexTitle;
                    }
                }
            }
        }
        //
        std::map < std::string, double > newPOIvalUp;
        std::map < std::string, double > newPOIvalDo;
        std::vector < std::string > npList;
        for(auto syst : fSystematics){
            if(Common::FindInStringVector(npList,syst->fNuisanceParameter)<0){
                npList.push_back(syst->fNuisanceParameter);
                for(int ud=0;ud<2;ud++){
                    //Be sure to take the initial values of the NP
                    ws->loadSnapshot("InitialStateModelGlob");
                    ws->loadSnapshot("InitialStateModelNuis");
                    // - create Asimov with that NP fixed to +/-1sigma
                    std::map < std::string, double > npVal;
                    npVal = fFitNPValues;
                    if(ud==0) npVal["alpha_"+syst->fNuisanceParameter] =  1;
                    if(ud==1) npVal["alpha_"+syst->fNuisanceParameter] = -1;
                    WriteInfoStatus("TRExFit::Fit","Systematic "+syst->fNuisanceParameter+"...");
                    data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npVal, fFitPOIAsimov ));
                    // - again a stat-only fit to that Asimov
                    fFitFixedNPs[syst->fNuisanceParameter] = 0;
                    fFitFixedNPs["alpha_"+syst->fNuisanceParameter] = 0;
                    npValues = PerformFit( ws.get(), data.get(), fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                    double newPOIval = npValues[fPOI];
                    if(ud==0) newPOIvalUp[syst->fNuisanceParameter] = newPOIval;
                    if(ud==1) newPOIvalDo[syst->fNuisanceParameter] = newPOIval;
                }
                std::string category = syst->fCategory;
                if(syst->fSubCategory!="") category = syst->fSubCategory;
                if(Common::FindInStringVector(systGroupNames,category)<0) systGroupNames.push_back(category);
                if(fabs(newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval)>fNonProfileFitSystThreshold
                || fabs(newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval)>fNonProfileFitSystThreshold)
                    systGroups[category] = std::hypot(systGroups[category], ((std::fabs(newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval)+fabs(newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval))/2));
            }
        }
        // print:
        std::cout << "Results of non-profile fit:" << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "StatisticalError\t" << statUp << "\t" << statDo << std::endl;
        out       << "StatisticalError\t" << statUp << "\t" << statDo << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        tex       << "\\begin{tabular}{lr}" << std::endl;
        tex       << "\\hline" << std::endl;
        tex       << "\\hline" << std::endl;
        tex       << "  Source & Shift up / down";
        if(fPOIunit!="") tex       << " [" << fPOIunit << "]";
        tex       << "\\\\" << std::endl;
        tex       << "\\hline" << std::endl;
        tex       << "  Statistical & $+" << Form("%.2f",statUp) << "$ / $" << Form("%.2f",statDo) << "$ \\\\" << std::endl;
        tex       << "\\hline" << std::endl;
        //
        double totUp = 0.;
        double totDo = 0.;
        npList.clear();
        // MC stat
        std::cout << "Stat.MC\t" << MCstatUp << "\t" << MCstatDo << std::endl;
        out       << "Stat.MC\t" << MCstatUp << "\t" << MCstatDo << std::endl;
        tex       << "  MC-stat & $+" << Form("%.2f",MCstatUp) << "$ / $" << Form("%.2f",MCstatDo) << "$ \\\\" << std::endl;
        tex       << "\\hline" << std::endl;
        totUp = std::hypot(totUp,MCstatUp);
        totDo = std::hypot(totDo,MCstatDo);
        // MC stat for separate gamma samples
        for(auto sepGammaPair : MCstatUpSample){
            std::string smpName = sepGammaPair.first;
            std::cout << "Stat." << smpName << "\t" << MCstatUpSample[smpName] << "\t" << MCstatDoSample[smpName] << std::endl;
            out       << "Stat." << smpName << "\t" << MCstatUpSample[smpName] << "\t" << MCstatDoSample[smpName] << std::endl;
            tex       << "  Stat (" << smpTexTitle[smpName] << ") & $+" << Form("%.2f",MCstatUpSample[smpName]) << "$ / $" << Form("%.2f",MCstatDoSample[smpName]) << "$ \\\\" << std::endl;
            tex       << "\\hline" << std::endl;
            totUp = std::hypot(totUp,MCstatUpSample[smpName]);
            totDo = std::hypot(totDo,MCstatDoSample[smpName]);
        }
        std::cout << "-----------------------------------" << std::endl;
        // systematics
        for(auto syst : fSystematics){
            if(syst->fName.find("stat_")!=std::string::npos && syst->fType==Systematic::SHAPE) continue;
            if(Common::FindInStringVector(npList,syst->fNuisanceParameter)<0){
                npList.push_back(syst->fNuisanceParameter);
                std::cout << syst->fNuisanceParameter;
                out       << syst->fNuisanceParameter;
                if(TRExFitter::SYSTTEX[syst->fNuisanceParameter]!="") tex << "  " << TRExFitter::SYSTTEX[syst->fNuisanceParameter];
                else                                                  tex << "  " << TRExFitter::SYSTMAP[syst->fNuisanceParameter];
                // - up and down
                for(int ud=0;ud<2;ud++){
                    double valUp = newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval;
                    double valDo = newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval;
                    if(ud==0) std::cout << "\t" << valUp;
                    if(ud==1) std::cout << "\t" << valDo;
                    if(ud==0) out       << "\t" << valUp;
                    if(ud==1) out       << "\t" << valDo;
                    if(ud==0) tex       << " & " << Form("$%s%.2f$",(valUp>=0 ? "+" : "-"),fabs(valUp));
                    if(ud==1) tex       << " / " << Form("$%s%.2f$",(valDo>=0 ? "+" : "-"),fabs(valDo));
                    if(fabs(valUp)>fNonProfileFitSystThreshold || fabs(valDo)>fNonProfileFitSystThreshold){
                        if(ud==0 && valUp>0) totUp = std::hypot(totUp,valUp);
                        if(ud==0 && valUp<0) totDo = std::hypot(totDo,valUp);
                        if(ud==1 && valDo>0) totUp = std::hypot(totUp,valDo);
                        if(ud==1 && valDo<0) totDo = std::hypot(totDo,valDo);
                    }
                }
                std::cout << std::endl;
                out       << std::endl;
                tex       << "\\\\" << std::endl;
            }
        }
        totUp = std::hypot(totUp,fNonProfileFitSystThreshold);
        totDo = std::hypot(totDo,fNonProfileFitSystThreshold);
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "TotalSystematic\t" << totUp << "\t-" << totDo << std::endl;
        out       << "TotalSystematic\t" << totUp << "\t-" << totDo << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "TotalStat+Syst\t" << std::hypot(totUp,statUp) << "\t-" << std::hypot(totDo,statDo) << std::endl;
        out       << "TotalStat+Syst\t" << std::hypot(totUp,statUp) << "\t-" << std::hypot(totDo,statDo) << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        out.close();
        tex << "\\hline" << std::endl;
        tex << "  Total systematics & $+" << Form("%.2f",totUp) << "$ / $-" << Form("%.2f",totDo) << "$ \\\\" << std::endl;
        tex << "\\hline" << std::endl;
        tex << "  Total stat+syst & $+"   << Form("%.2f",std::hypot(totUp,statUp))     << "$ / $-" << Form("%.2f",std::hypot(totDo,statDo))     << "$ \\\\" << std::endl;
        tex << "\\hline" << std::endl;
        tex << "\\hline" << std::endl;
        tex << "\\end{tabular}" << std::endl;
        tex.close();
        fVarNameMinos = varMinosTmp; // retore Minos settings
        //
        // Systematics merged according to syst groups
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "- Systematic impact per category  -" << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        tex2      << "\\begin{tabular}{lr}" << std::endl;
        tex2      << "\\hline" << std::endl;
        tex2      << "\\hline" << std::endl;
        std::cout << "Data statistics" << "\t" << (fabs(statUp)+fabs(statDo))/2. << std::endl;
        out2      << "Data statistics" << "\t" << (fabs(statUp)+fabs(statDo))/2. << std::endl;
        tex2      << "Data statistics" << " & " << Form("$%.2f$",(fabs(statUp)+fabs(statDo))/2.) << " \\\\" << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        tex2      << "\\hline" << std::endl;
        std::cout << "MC background stat." << "\t" << (fabs(MCstatUp)+fabs(MCstatDo))/2. << std::endl;
        out2      << "MC background stat." << "\t" << (fabs(MCstatUp)+fabs(MCstatDo))/2. << std::endl;
        tex2      << "MC background stat." << " & " << Form("$%.2f$",(fabs(MCstatUp)+fabs(MCstatDo))/2.) << " \\\\" << std::endl;
        for(auto sepGammaPair : MCstatUpSample){
            std::string smpName = sepGammaPair.first;
            std::cout << smpTexTitle[smpName] << " stat." << "\t" << (fabs(MCstatUpSample[smpName])+fabs(MCstatDoSample[smpName]))/2. << std::endl;
            out2      << smpTexTitle[smpName] << " stat." << "\t" << (fabs(MCstatUpSample[smpName])+fabs(MCstatDoSample[smpName]))/2. << std::endl;
            tex2      << smpTexTitle[smpName] << " stat." << " & " << Form("$%.2f$",(fabs(MCstatUpSample[smpName])+fabs(MCstatDoSample[smpName]))/2.) << " \\\\" << std::endl;
        }
        std::cout << "-----------------------------------" << std::endl;
        tex2      << "\\hline" << std::endl;
        for(auto systGroupName : systGroupNames){
            if(systGroupName=="") continue;
            std::cout << systGroupName << "\t" << systGroups[systGroupName] << std::endl;
            out2      << systGroupName << "\t" << systGroups[systGroupName] << std::endl;
            tex2      << systGroupName << " & " << Form("$%.2f$",systGroups[systGroupName]) << " \\\\" << std::endl;
        }
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "TotalSystematic\t" << (std::fabs(totUp)+std::fabs(totDo))/2. << std::endl;
        std::cout << "TotalStat+Syst\t" << (std::hypot(totUp,statUp)+std::hypot(totDo,statDo))/2. << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        tex2 << "\\hline" << std::endl;
        tex2 << "Total systematic uncertainty & " << Form("$%.2f$",(std::fabs(totUp)+std::fabs(totDo))/2.) << " \\\\" << std::endl;
        tex2 << "\\hline" << std::endl;
        tex2 << "Total & "   << Form("$%.2f$",(std::hypot(totUp,statUp)+std::hypot(totDo,statDo))/2.) << " \\\\" << std::endl;
        tex2 << "\\hline" << std::endl;
        tex2 << "\\hline" << std::endl;
        tex2 << "\\end{tabular}" << std::endl;
        tex2.close();
    }

    //
    // Calls the  function to create LH scan with respect to a parameter
    //
    if(fVarNameLH.size()>0 && !isLHscanOnly && !fParal2D){
        //
        // Don't do it if you did a non-profile fit (FIXME)
        if(fDoNonProfileFit){
            WriteWarningStatus("TRExFit::Fit","Better not to perform LH scan if you did non-profile fit with scan on systematics. Skipping LH scan.");
        }
        else{
            if (fVarNameLH[0]=="all"){
                for(std::map<std::string,std::string>::iterator it=TRExFitter::SYSTMAP.begin(); it!=TRExFitter::SYSTMAP.end(); ++it){
                    GetLikelihoodScan( ws.get(), it->first, data.get());
                }
            }
            else{
                for(unsigned int i=0; i<fVarNameLH.size(); ++i){
                    GetLikelihoodScan( ws.get(), fVarNameLH[i], data.get());
                }
            }
        }
    }
    if (isLHscanOnly && !fParal2D){
        if (fVarNameLH.size() == 0){
            WriteErrorStatus("TRExFit::Fit","Did not provide any LH scan parameter and running LH scan only. This is not correct.");
            exit(EXIT_FAILURE);
        }
        if (fVarNameLH[0]=="all"){
            WriteWarningStatus("TRExFit::Fit","You are running LHscan only option but running it for all parameters. Will not parallelize!.");
            for(std::map<std::string,std::string>::iterator it=TRExFitter::SYSTMAP.begin(); it!=TRExFitter::SYSTMAP.end(); ++it){
                GetLikelihoodScan( ws.get(), it->first, data.get());
            }
        } else {
            GetLikelihoodScan( ws.get(), fVarNameLH[0], data.get());
        }
    }

    // run 2D likelihood scan
    if(fVarName2DLH.size()>0){
        for (const auto & ipair : fVarName2DLH) {
            Get2DLikelihoodScan( ws.get(), ipair, data.get());
        }
    }
}

//__________________________________________________________________________________
//
RooDataSet* TRExFit::DumpData( RooWorkspace *ws,  std::map < std::string, int > &regionDataType, std::map < std::string, double > &npValues, const double poiValue  ){
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);

    //
    // This function dumps a RooDataSet object using the input informations provided by the user
    //    |-> Used when testing Fit response (inject one NP in data and check fit result)
    //    |-> Used when using fit results in some regions to generate Asimov data in blinded regions
    //
    WriteDebugStatus("TRExFit::DumpData", "Dumping data with the following parameters");
    WriteDebugStatus("TRExFit::DumpData", "    * Regions data type ");
    for( const std::pair < std::string, int >& dataType : regionDataType ){
        WriteDebugStatus("TRExFit::DumpData", "       - Region: " + dataType.first + "       DataType: " + std::to_string(dataType.second));
    }
    if(npValues.size()){
        WriteDebugStatus("TRExFit::DumpData", "    * Injected NP values ");
        for ( const std::pair < std::string, double >& npValue : npValues ){
            WriteDebugStatus("TRExFit::DumpData", "       - NP: " + npValue.first + "       Value: " + std::to_string(npValue.second));
        }
    }
    else {
        WriteDebugStatus("TRExFit::DumpData", "    * No NP values injected ");
    }
    WriteDebugStatus("TRExFit::DumpData", "    * POI value: " + std::to_string(poiValue) );

    RooStats::ModelConfig *mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    //Save the initial values of the NP
    ws->saveSnapshot("InitialStateModelGlob",   *mc->GetGlobalObservables());
    if (!(fStatOnly && fFitIsBlind)){
        ws->saveSnapshot("InitialStateModelNuis",   *mc->GetNuisanceParameters());
    }

    //Be sure to take the initial values of the NP
    ws->loadSnapshot("InitialStateModelGlob");
    if (!fStatOnly){
        ws->loadSnapshot("InitialStateModelNuis");
    }

    if (fBinnedLikelihood) {
        FitUtils::SetBinnedLikelihoodOptimisation(ws);
    }

    //Creating a set
    static const std::string weightName="weightVar";
    RooArgSet obsAndWeight;
    obsAndWeight.add(*mc->GetObservables());

    RooRealVar* weightVar = nullptr;
    if ( !(weightVar = ws->var(weightName.c_str())) ){
        ws->import(*(new RooRealVar(weightName.c_str(), weightName.c_str(), 1,0,10000000)));
        weightVar = ws->var(weightName.c_str());
    }
    obsAndWeight.add(*ws->var(weightName.c_str()));
    ws->defineSet("obsAndWeight",obsAndWeight);

    //
    // Getting observed data (in case some regions are unblinded)
    //
    RooDataSet* realData = static_cast<RooDataSet*>(ws->data("obsData"));

    //
    // Set some parameters for the Asimov production
    //     |-> Values of NPs
    //     |-> Values of POI
    //

    //-- POI
    RooRealVar * poi = static_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
    if (!poi){
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        WriteErrorStatus("TRExFit::DumpData", "Cannot find POI in workspace, exiting...");
        exit(EXIT_FAILURE);
    }
    poi -> setVal(poiValue);

    //-- Nuisance parameters
    RooRealVar* var(nullptr);
    std::unique_ptr<TIterator> npIterator(mc->GetNuisanceParameters()->createIterator());
    while( (var = static_cast<RooRealVar*>(npIterator->Next()))) {
        std::map < std::string, double >::const_iterator it_npValue = npValues.find( var -> GetName() );
        if( it_npValue != npValues.end() ){
            var -> setVal(it_npValue -> second);
        }
    }

    //Looping over regions
    std::map<std::string, RooDataSet*> asimovDataMap;
    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
    RooCategory* channelCat = (RooCategory*)(&simPdf->indexCat());
    std::unique_ptr<TIterator> iter(channelCat->typeIterator());
    RooCatType* tt = nullptr;
    int iFrame = 0;
    int i = 0;
    while( (tt = static_cast<RooCatType*>(iter->Next()))) {

        channelCat->setIndex(i);
        iFrame++;
        i++;

        //Check the type of data to store for this region !
        int dataType = Region::ASIMOVDATA;//default is AsimovData
        std::map < std::string, int >::const_iterator it_dataType = regionDataType.find( channelCat->getLabel() );
        if( it_dataType == regionDataType.end() ){
            const std::string temp_string = channelCat->getLabel();
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            WriteWarningStatus("TRExFit::DumpData", "The following region is not specified in the inputs to the function (" + temp_string + "): use Asimov");
            WriteWarningStatus("TRExFit::DumpData", "   This SHOULD NOT HAPPEN ! Please check if everything is fine !");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        }
        else {
            dataType = regionDataType[channelCat->getLabel()];
        }

        //A protection: if there is no real observed data, use only ASIMOV (but print a warning)
        if(dataType==Region::REALDATA && !realData){
            std::string temp_string = channelCat->getLabel();
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            WriteWarningStatus("TRExFit::DumpData", "You want real data for channel " + temp_string + " but none is available in the workspace. Using Asimov instead.");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
            dataType = Region::ASIMOVDATA;
        }

        if(dataType==Region::ASIMOVDATA){
            // Get pdf associated with state from simpdf
            RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;

            // Generate observables defined by the pdf associated with this state
            RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

            RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),RooFit::WeightVar(*weightVar));
            RooRealVar* thisObs = static_cast<RooRealVar*>(obstmp->first());
            const double expectedEvents = pdftmp->expectedEvents(*obstmp);

            for(int jj=0; jj<thisObs->numBins(); ++jj){
                thisObs->setBin(jj);
                const double thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
                if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < 1e18) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
            }
            obsDataUnbinned->Print();
            if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
                exit(1);
            }
            asimovDataMap[std::string(channelCat->getLabel())] = obsDataUnbinned;

        } else if(dataType==Region::REALDATA) {
            RooAbsData *datatmp = realData->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),tt->GetName()));
            asimovDataMap[std::string(channelCat->getLabel())] = static_cast<RooDataSet*>(datatmp);
        }
    }

    RooDataSet *asimovData = new RooDataSet("newasimovData",
                                            "newasimovData",
                                            RooArgSet(obsAndWeight,*channelCat),
                                            Index(*channelCat),
                                            Import(asimovDataMap),
                                            WeightVar(*weightVar));

    ws->loadSnapshot("InitialStateModelGlob");
    if (!fStatOnly){
        ws->loadSnapshot("InitialStateModelNuis");
    }

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    return asimovData;
}

//__________________________________________________________________________________
//
std::map < std::string, double > TRExFit::PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType, bool save, int debugLevel ){

    if (debugLevel < 1) std::cout.setstate(std::ios_base::failbit);
    std::map < std::string, double > result;

    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////

    // prepare vectors for starting point of normfactors
    std::vector<std::string> NPnames;
    std::vector<double> NPvalues;
    double poiInitial(0.);
    for(const auto& inf : fNormFactors) {
        if (inf->fName == fPOI) {
            poiInitial = inf->fNominal;
            continue;
        }
        NPnames. emplace_back(inf->fName);
        NPvalues.emplace_back(inf->fNominal);
    }
    //
    // Fit configuration (SPLUSB or BONLY)
    //
    FittingTool fitTool{};
    fitTool.SetUseHesse(true);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetDebug(debugLevel);
    fitTool.SetNCPU(fCPU);
    if(fitType==BONLY){
        fitTool.ValPOI(0.);
        fitTool.ConstPOI(true);
    } else if(fitType==SPLUSB){
        fitTool.ValPOI(poiInitial);
        fitTool.ConstPOI(false);
    }
    fitTool.SetNPs( NPnames,NPvalues );
    fitTool.SetRandomNP(fRndRange, fUseRnd, fRndSeed);
    if(fStatOnly){
        if(!fGammasInStatOnly) fitTool.NoGammas();
        fitTool.NoSystematics();
    }

    //
    // Fit starting from custom point
    if(fFitResultsFile!=""){
        ReadFitResults(fFitResultsFile);
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitResults->fNuisPar) {
            npNames.push_back( inp->fName);
            npValues.push_back(inp->fFitValue);
        }
        fitTool.SetNPs(npNames,npValues);
    }

    //
    // Set Minos
    if(fVarNameMinos.size()>0){
        WriteDebugStatus("TRExFit::PerformFit", "Setting the variables to use MINOS with");
        fitTool.UseMinos(fVarNameMinos);
    }

    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    if(inputData){
        data = inputData;
    } else {
        WriteWarningStatus("TRExFit::PerformFit", "You didn't provide inputData => will use the observed data !");
        data = static_cast<RooDataSet*>(ws->data("obsData"));
        if(data==nullptr){
            WriteWarningStatus("TRExFit::PerformFit", "No observedData found => will use the Asimov data !");
            data = static_cast<RooDataSet*>(ws->data("asimovData"));
        }
        inputData = data;
    }

    //
    // For stat-only fit on data:
    // - read fit resutls
    // - fix all NP to fitted ones before fitting
    if(fStatOnlyFit){
        WriteDebugStatus("TRExFit::PerformFit", "Fitting stat-only: reading fit results from full fit from file: ");
        WriteDebugStatus("TRExFit::PerformFit", "  " + fName+"/Fits/"+fInputName+fSuffix+".txt");
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp: fFitResults->fNuisPar) {
            if(!fFixNPforStatOnlyFit && Common::FindInStringVector(fNormFactorNames,inp->fName)>=0) continue;
            if(!fFixNPforStatOnlyFit && Common::FindInStringVector(fShapeFactorNames,inp->fName)>=0) continue;
            npNames.push_back( inp->fName);
            npValues.push_back(inp->fFitValue);
        }
        fitTool.FixNPs(npNames,npValues);
    }

    // FixNP
    if(fFitFixedNPs.size()>0){
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& nuisParToFix : fFitFixedNPs){
            npNames.push_back( nuisParToFix.first );
            npValues.push_back( nuisParToFix.second );
        }
        fitTool.FixNPs(npNames,npValues);
    }

    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, fNormFactors);

    // save snapshot before fit
    ws->saveSnapshot("snapshot_BeforeFit_POI", *(mc->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_BeforeFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_BeforeFit_GO" , *(mc->GetGlobalObservables())    );

    //
    // Get initial ikelihood value from Asimov
    double nll0 = 0.;
    if (fBlindedParameters.size() > 0) std::cout.setstate(std::ios_base::failbit);
    if(fGetGoodnessOfFit) nll0 = fitTool.FitPDF( mc, simPdf, static_cast<RooDataSet*>(ws->data("asimovData")), false, true );

    // save snapshot after fit
    ws->saveSnapshot("snapshot_AfterFit_POI", *(mc->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_AfterFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(mc->GetGlobalObservables())    );

    //
    // Get number of degrees of freedom
    // - number of bins
    int ndof = inputData->numEntries();
    // - minus number of free & non-constant parameters
    int nNF = 0;
    for(const auto& inf : fNormFactors) {
        if(inf->fConst) continue;
        if(fFitType==BONLY && fPOI==inf->fName) continue;
        // skip if it's a morphing parameter
        if(inf->fName.find("morph_")!=std::string::npos) continue;
        // skip if it has an "Expression"
        if(inf->fExpression.first!="") continue;
        // skip if not in the ws (e.g. because assigned to a sample or region not present in the fit)
        if(!ws->obj(inf->fName.c_str())) continue;
        nNF++;
    }
    ndof -= nNF;

    // Performs the fit
    const double nll = fitTool.FitPDF( mc, simPdf, data );
    if (debugLevel < 1 && fBlindedParameters.size() == 0) std::cout.clear();
    if(save){
        if(fBootstrap!="" && fBootstrapIdx>=0){
            gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
            if(fStatOnlyFit) fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+"_statOnly.txt", fBlindedParameters);
            else             fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".txt", fBlindedParameters);
        }
        else{
            gSystem -> mkdir((fName+"/Fits/").c_str(),true);
            if(fStatOnlyFit) fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_statOnly.txt", fBlindedParameters);
            else             fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+".txt", fBlindedParameters);
        }
    }
    result = fitTool.ExportFitResultInMap();
    if (fBlindedParameters.size() > 0) std::cout.clear();

    // If SaturatedModel used, superseed Asimov-based GOF
    if(fSaturatedModel && fGetGoodnessOfFit && !fDoGroupedSystImpactTable){
        ws->loadSnapshot("snapshot_BeforeFit_POI");
        ws->loadSnapshot("snapshot_BeforeFit_GO");
        ws->loadSnapshot("snapshot_BeforeFit_NP");
        //
        // perform fit to saturated model and store resulting nll as nll0
        nll0 = fitTool.FitPDF( mc, simPdf, data, false, false, true );
        // could be removed, but kept for debugging FIXME
        if(save){
            if(fStatOnlyFit) fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_saturatedModel_statOnly.txt", fBlindedParameters);
            else             fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_saturatedModel.txt", fBlindedParameters);
        }
    }

    //
    // Goodness of fit
    if(fGetGoodnessOfFit && !fDoGroupedSystImpactTable){
        const double deltaNLL = nll-nll0;
        const double prob = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndof);
        WriteInfoStatus("TRExFit::PerformFit", "----------------------- -------------------------- -----------------------");
        WriteInfoStatus("TRExFit::PerformFit", "----------------------- GOODNESS OF FIT EVALUATION -----------------------");
        WriteInfoStatus("TRExFit::PerformFit", "  NLL0        = " + std::to_string(nll0));
        WriteInfoStatus("TRExFit::PerformFit", "  NLL         = " + std::to_string(nll));
        WriteInfoStatus("TRExFit::PerformFit", "  ndof        = " + std::to_string(ndof));
        WriteInfoStatus("TRExFit::PerformFit", "  dNLL        = " + std::to_string(deltaNLL));
        WriteInfoStatus("TRExFit::PerformFit", "  2dNLL/nof   = " + std::to_string(2.*deltaNLL/ndof));
        WriteInfoStatus("TRExFit::PerformFit", "  probability = " + std::to_string(prob));
        WriteInfoStatus("TRExFit::PerformFit", "----------------------- -------------------------- -----------------------");
        WriteInfoStatus("TRExFit::PerformFit", "----------------------- -------------------------- -----------------------");
    }

    //
    // Load snapshots after nominal fit (needed in case GetGoodnessOfFit test with saturated model is evaluated)
    if(fSaturatedModel && fGetGoodnessOfFit && !fDoGroupedSystImpactTable){
        ws->loadSnapshot("snapshot_AfterFit_POI");
        ws->loadSnapshot("snapshot_AfterFit_GO");
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }

    //
    // grouped systematics impact
    if(fDoGroupedSystImpactTable){
        // name of file to write results to
        std::string outNameGroupedImpact = fName+"/Fits/GroupedImpact"+fSuffix;
        if(fBootstrap!="" && fBootstrapIdx>=0){
            gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
            outNameGroupedImpact = fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+"GroupedImpact"+fSuffix;
        }
        if(fGroupedImpactCategory!="all") outNameGroupedImpact += "_"+fGroupedImpactCategory;
        outNameGroupedImpact += ".txt";

        ProduceSystSubCategoryMap();                        // fill fSubCategoryImpactMap first
        fitTool.SetSystMap( fSubCategoryImpactMap );     // hand over the map to the FittingTool
        fitTool.GetGroupedImpact( mc, simPdf, data, ws, fGroupedImpactCategory, outNameGroupedImpact);
    }

    delete fFitResults;
    fFitResults = nullptr;
    return result;
}

//__________________________________________________________________________________
//
std::unique_ptr<RooWorkspace> TRExFit::PerformWorkspaceCombination( std::vector < std::string > &regionsToFit ) const{

    //
    // Definition of the fit regions
    //
    std::vector < RooWorkspace* > vec_ws;
    std::vector < std::string > vec_chName;
    std::unique_ptr<RooStats::HistFactory::Measurement> measurement(nullptr);
    //
    // Take the measurement from the combined workspace, to be sure to have all the systematics (even the ones which are not there in the first region)
    std::unique_ptr<TFile> rootFileCombined(nullptr);
    if(fBootstrap!="" && fBootstrapIdx>=0) {
        rootFileCombined.reset(TFile::Open( (fName+"/RooStats/"+fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read"));
    } else {
        rootFileCombined.reset(TFile::Open( (fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read"));
    }
    if(rootFileCombined) measurement = std::unique_ptr<RooStats::HistFactory::Measurement>(static_cast<RooStats::HistFactory::Measurement*>(rootFileCombined -> Get( (fInputName+fSuffix).c_str())));
    //
    std::vector<std::unique_ptr<TFile> > file_vec;
    for(int i_ch=0; i_ch < fNRegions; ++i_ch) {
        bool isToFit = false;
        for(unsigned int iRegion = 0; iRegion < regionsToFit.size(); ++iRegion){
            if(fRegions[i_ch] -> fName == regionsToFit[iRegion]){
                isToFit = true;
                break;
            }
        }
        if (!isToFit) continue;
        std::string fileName = fName+"/RooStats/"+fInputName+"_"+fRegions[i_ch]->fName+"_"+fInputName+fSuffix+"_model.root";
        if(fBootstrap!="" && fBootstrapIdx>=0) fileName = fName+"/RooStats/"+fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName+"_"+fRegions[i_ch]->fName+"_"+fInputName+fSuffix+"_model.root";
        file_vec.emplace_back(std::move(TFile::Open(fileName.c_str())));
        if (!file_vec.back()) {
            WriteErrorStatus("TRExFit::PerformWorkspaceCombination", "Input file with the workspace doesnt exist!");
            exit(EXIT_FAILURE);
        }
        RooWorkspace *tmp_ws = static_cast<RooWorkspace*>(file_vec.back()->Get((fRegions[i_ch]->fName).c_str()));
        if(!tmp_ws){
            WriteErrorStatus("TRExFit::PerformWorkspaceCombination", "The workspace (\"" + fRegions[i_ch] -> fName + "\") cannot be found in file " + fileName + ". Please check !");
            return nullptr;
        }
        vec_ws.emplace_back(std::move(tmp_ws));
        vec_chName.emplace_back(fRegions[i_ch] -> fName);
        // if failed to get the measurement from the combined ws, take it from the first region
        if(!measurement){
            measurement = std::unique_ptr<RooStats::HistFactory::Measurement>(static_cast<RooStats::HistFactory::Measurement*>(file_vec.back() -> Get( (fInputName+fSuffix).c_str())));
        }
    }


    //
    // Create the HistoToWorkspaceFactoryFast object to perform safely the combination
    //
    if(!measurement){
        WriteErrorStatus("TRExFit::PerformWorkspaceCombination", "The measurement object has not been retrieved ! Please check.");
        return nullptr;
    }
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    RooStats::HistFactory::HistoToWorkspaceFactoryFast factory(*measurement);

    // Creating the combined model
    std::unique_ptr<RooWorkspace> ws(factory.MakeCombinedModel( vec_chName, vec_ws ));

    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws.get(), *measurement );
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    for (auto& iws : vec_ws) {
        delete iws;
    }

    for (auto& ifile : file_vec) {
        ifile->Close();
    }

    if (rootFileCombined) rootFileCombined->Close();

    return ws;
}

//__________________________________________________________________________________
//
void TRExFit::PlotFittedNP(){
    if(fStatOnly || fStatOnlyFit){
        WriteInfoStatus("TRExFit::PlotFittedNP", "Stat only fit => No NP Pull plots generated.");
    }
    //
    // plot the NP fit pull plot
    //
    if(fStatOnlyFit){
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+"_statOnly.txt");
    }
    else{
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
    }
    if(fFitResults){
        fFitResults->fNuisParToHide = fVarNameHide;
        std::set < std::string > npCategories;
        for(const auto& syst : fSystematics) {
            npCategories.insert(syst->fCategory);
        }
        for(const auto& format: TRExFitter::IMAGEFORMAT) {
            if(!fStatOnly && !fStatOnlyFit){
                fFitResults->DrawNPPulls(fName+"/NuisPar"+fSuffix+"."+format,"all",fNormFactors, fBlindedParameters);
                fFitResults->DrawGammaPulls(fName+"/Gammas"+fSuffix+"."+format, fBlindedParameters);
            }
            if(fStatOnlyFit){
                fFitResults->DrawNormFactors(fName+"/NormFactors"+fSuffix+"_statOnly."+format,fNormFactors, fBlindedParameters);
            }
            else{
                fFitResults->DrawNormFactors(fName+"/NormFactors"+fSuffix+"."+format,fNormFactors, fBlindedParameters);
            }

        }
        if(npCategories.size()>1 && !fStatOnly && !fStatOnlyFit){
            for(const std::string& cat : npCategories){
                std::string cat_for_name = cat;
                std::replace(cat_for_name.begin(), cat_for_name.end(), ' ', '_');
                std::replace(cat_for_name.begin(), cat_for_name.end(), '#', '_');
                std::replace(cat_for_name.begin(), cat_for_name.end(), '{', '_');
                std::replace(cat_for_name.begin(), cat_for_name.end(), '}', '_');
                for(const auto& format : TRExFitter::IMAGEFORMAT) {
                  fFitResults->DrawNPPulls(fName+"/NuisPar_"+cat_for_name+fSuffix+"."+format,cat,fNormFactors, fBlindedParameters);
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::PlotCorrelationMatrix(){
    if(fStatOnly || fStatOnlyFit){
        WriteInfoStatus("TRExFit::PlotCorrelationMatrix", "Stat only fit => No Correlation Matrix generated.");
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TRExFitter::CORRELATIONTHRESHOLD)
    ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
    if(fFitResults){
        fFitResults->fNuisParToHide = fVarNameHide;
        fFitResults->fOutFolder = fName;
        std::vector<std::string> formats;
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            formats.emplace_back(fName+"/CorrMatrix"+fSuffix+"."+format);
        }
        fFitResults->DrawCorrelationMatrix(formats,
                                           fuseGammasForCorr,
                                           fHEPDataFormat,
                                           TRExFitter::CORRELATIONTHRESHOLD);
    }
}

//__________________________________________________________________________________
//
void TRExFit::PlotUnfoldedData() const {
    if (fFitType != TRExFit::FitType::UNFOLDING) return;

    WriteInfoStatus("TRExFit::PlotUnfoldedData", "Producing unfolded plots...");

    std::unique_ptr<TFile> input(TFile::Open((fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
    if (!input) {
        WriteErrorStatus("TRExFit::PlotUnfoldedData", "Cannot read file from " + fName + "/UnfoldingHistograms/FoldedHistograms.root");
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<TH1D> truth(dynamic_cast<TH1D*>(input->Get("truth_distribution")));
    if (!truth) {
        WriteErrorStatus("TRExFit::PlotUnfoldedData", "Cannot read the truth distribution");
        exit(EXIT_FAILURE);
    }
    truth->SetDirectory(nullptr);

    UnfoldingResult unfolded;
    unfolded.SetTruthDistribution(truth.get());
    unfolded.SetNormXSec(fUnfoldNormXSec);

    // pass the fit results to the tool
    for (int i = 0; i < fNumberUnfoldingTruthBins; ++i) {
        std::string name = "Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
        // if it's last bin and getting norm xsec:
        if(fUnfoldNormXSec && i==fUnfoldNormXSecBinN-1){
            // get value and error of non-independent parameter from workspace and fit result
            // - get workspace
            const std::string wsFileName = fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root";
            std::unique_ptr<TFile> wsFile(TFile::Open(wsFileName.c_str(),"read"));
//             std::unique_ptr<RooWorkspace> ws(static_cast<RooWorkspace*>(wsFile->Get("combined"))); // for some reasons induces a crash
            RooWorkspace* ws = static_cast<RooWorkspace*>(wsFile->Get("combined"));
            // - get fit-result
            const std::string frFileName = fName+"/Fits/"+fName+fSuffix+".root";
            std::unique_ptr<TFile> frFile(TFile::Open(frFileName.c_str(),"read"));
            std::unique_ptr<RooFitResult> fr(nullptr);
            // loop on keys in file
            for(auto *key : *gDirectory->GetListOfKeys()){
                std::string keyName = key->GetName();
                if(keyName.find("nll_simPdf")!=std::string::npos){
                    fr.reset(static_cast<RooFitResult*>(frFile->Get(keyName.c_str())));
                    break;
                }
            }
            // calculate propagated value and error in 2 ways:
            // 1. by hand:
            double num = 1.;
            double den = 1.;
            const double N = truth->Integral();
            const double Ni = truth->GetBinContent(i+1);
            for (int j = 0; j < fNumberUnfoldingTruthBins; ++j){
                if(j==fUnfoldNormXSecBinN-1) continue;
                const double Nj = truth->GetBinContent(j+1);
                num -= fFitResults->GetNuisParValue("Bin_" + Common::IntToFixLenStr(j+1) + "_mu") * Nj/N;
                den -= Nj/N;
            }
            const double mean0 = num/den;
            double error0 = 0.;
            for (int j = 0; j < fNumberUnfoldingTruthBins; ++j){
                if(j==fUnfoldNormXSecBinN-1) continue;
                const std::string par_j = "Bin_" + Common::IntToFixLenStr(j+1) + "_mu";
                const double Nj = truth->GetBinContent(j+1);
                for (int k = 0; k < fNumberUnfoldingTruthBins; ++k){
                    if(k==fUnfoldNormXSecBinN-1) continue;
                    const std::string par_k = "Bin_" + Common::IntToFixLenStr(k+1) + "_mu";
                    const double Nk = truth->GetBinContent(k+1);
                    double rho = fFitResults->fCorrMatrix->GetCorrelation(par_j,par_k);
                    double err_j = fFitResults->GetNuisParErrUp(par_j);
                    double err_k = fFitResults->GetNuisParErrUp(par_k);
                    error0 += rho*(Nj/Ni)*(Nk/Ni)*err_j*err_k;
                }
            }
            error0 = sqrt(error0);
            const double up0   = mean0 + error0;
            const double down0 = mean0 - error0;
            // 2. within RooFit:
            std::unique_ptr<RooFormulaVar> muN(nullptr);
            for(auto *arg : ws->allFunctions()){
                std::string argName = arg->GetName();
                if(argName=="Bin_" + Common::IntToFixLenStr(fUnfoldNormXSecBinN) + "_mu"){
                    muN.reset(static_cast<RooFormulaVar*>(arg));
                }
            }
            // propagate fit results to ws
            for(auto *par : fr->floatParsFinal()){
                ws->var( par->GetName() )->setVal( (static_cast<RooRealVar*>(par))->getVal() );
                ws->var( par->GetName() )->setError( (static_cast<RooRealVar*>(par))->getError() );
            }
            const double mean = muN!=nullptr ? muN->getVal() : 0.;
            const double error = muN!=nullptr ? muN->getPropagatedError(*fr) : 0.;
            const double up   = mean + error;
            const double down = mean - error;
            // Choose which one to take
            if(ws==nullptr || fr==nullptr || muN==nullptr){
                // if problems in getting roofit fit result, just use the by-hand one:
                WriteWarningStatus("TRExFit::PlotUnfoldedData","No suitable RooFitResult or norm factor in WS found: taking by-hand calculation.");
                unfolded.AddFitValue(mean0, up0, down0);
            }
            else{
                // otherwise, check compatibility:
                // - for mean (1% threshold)
                if(fabs(mean-mean0)/(0.5*(mean+mean0))>0.01) WriteWarningStatus("TRExFit::PlotUnfoldedData","Propagation of nominal unfolding result to non-independent bin giving incompatible results: " + std::to_string(mean0) + " vs. " + std::to_string(mean) + ". Note that the latter will be used.");
                // - for error (1% threshold)
                if(fabs(error-error0)/(0.5*(error+error0))>0.01) WriteWarningStatus("TRExFit::PlotUnfoldedData","Propagation of unfolding error to non-independent bin giving incompatible results: " + std::to_string(error0) + " vs. " + std::to_string(error) + ". Note that the latter will be used.");
                // in any case use the roofit one
                unfolded.AddFitValue(mean, up, down);
            }
            wsFile->Close();
        }
        else{
            const double mean = fFitResults->GetNuisParValue(name);
            const double up   = mean + fFitResults->GetNuisParErrUp(name);
            const double down = mean + fFitResults->GetNuisParErrDown(name);
            unfolded.AddFitValue(mean, up, down);
        }
    }

    std::unique_ptr<TH1D> data               = unfolded.GetUnfoldedResult();
    std::unique_ptr<TGraphAsymmErrors> error = unfolded.GetUnfoldedResultErrorBand();

    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    converter.WriteUnfolding(error.get(), fName);
    if (fHEPDataFormat) {
        converter.WriteUnfoldingHEPData(error.get(), fUnfoldingTitleX, fName);
    }

    PlotUnfold(data.get(), error.get());

    // Dump results into a text file
    auto text = std::make_unique<std::ofstream>();
    text->open(fName + "/Fits/UnfoldedResults.txt");
    if (!text->is_open() || !text->good()) {
        WriteErrorStatus("TRExFit::PlotUnfoldedData", "Cannot open text file in: " + fName + "/Fits/UnfoldedResults.txt");
        exit(EXIT_FAILURE);
    }

    unfolded.DumpResults(text.get());
    text->close();

    input->Close();
}

//__________________________________________________________________________________
//
void TRExFit::GetLimit(){
    std::unique_ptr<RooDataSet> data(nullptr);
    
    //
    // Read NPvalues from fit-result file
    //
    if(fFitNPValuesFromFitResults!=""){
        WriteInfoStatus("TRExFit::Fit","Setting NP values for Asimov data-set creation from fit results stored in file " + fFitNPValuesFromFitResults + "...");
        fFitNPValues = NPValuesFromFitResults(fFitNPValuesFromFitResults);
    }
    
    //
    // If a workspace file name is specified, do simple limit
    //
    int sigDebug = 3 - TRExFitter::DEBUGLEVEL;
    if (sigDebug < 0) sigDebug = 0;
    if(fWorkspaceFileName!=""){
        std::string dataName = "obsData";
        if(fLimitIsBlind) dataName = "asimovData";
        runAsymptoticsCLs(fWorkspaceFileName.c_str(), "combined", "ModelConfig", dataName.c_str(), fLimitParamName.c_str(), fLimitParamValue, (fLimitOutputPrefixName+fSuffix).c_str(), (fName+"/Limits/").c_str(), fLimitIsBlind, fLimitsConfidence, "asimovData_0", fSignalInjection, fSignalInjectionValue, sigDebug);
    }
    else{
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std::string > regionsForLimit = ListRegionsToFit();
        std::map < std::string, int > regionsForLimitDataType = MapRegionDataTypes(regionsForLimit,fLimitIsBlind);
        std::map < std::string, double > npValues;
        //
        // flag if mixed Data / Asimov limit required
        // if mixed limit, perform a first fit on the regions with data only
        bool isMixedFit = DoingMixedFitting();
        if(isMixedFit && !fLimitIsBlind){
            std::vector < std:: string > regionsForFit = ListRegionsToFit(Region::REALDATA);
            //
            // Creates a combined workspace with the regions to be used *in the fit*
            //
            WriteInfoStatus("TRExFit::GetLimit","Creating ws for regions with real data only...");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
            std::unique_ptr<RooWorkspace> ws_forFit = PerformWorkspaceCombination( regionsForFit );
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            if (!ws_forFit){
                WriteErrorStatus("TRExFit::GetLimit","Cannot retrieve the workspace, exiting!");
                exit(EXIT_FAILURE);
            }
            //
            // Calls the PerformFit() function to actually do the fit
            //
            WriteInfoStatus("TRExFit::GetLimit","Performing a fit in regions with real data only...");
            npValues = PerformFit( ws_forFit.get(), data.get(), FitType::BONLY, false, TRExFitter::DEBUGLEVEL);
            WriteInfoStatus("TRExFit::GetLimit","Now will use the fit results to create the Asimov in the regions without real data!");
        }
        else{
            npValues = fFitNPValues; // FIXME: maybe should add a set of configurable NPs for limit only??
        }

        //
        // Create the final asimov dataset for limit setting
        //
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        std::unique_ptr<RooWorkspace> ws_forLimit = PerformWorkspaceCombination( regionsForLimit );
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        if (!ws_forLimit){
            WriteErrorStatus("TRExFit::GetLimit","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        data = std::unique_ptr<RooDataSet>(DumpData( ws_forLimit.get(), regionsForLimitDataType, npValues, npValues.find(fPOI)==npValues.end() ? 0 : npValues[fPOI] ));
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        
        //
        // Set the global observables to the specified NPValues if the corresponding flag is TRUE
        //
        if(fInjectGlobalObservables && !fFitNPValues.empty()) {
            FitUtils::InjectGlobalObservables(ws_forLimit.get(), fFitNPValues);
        }

        //
        FitUtils::DisableSaturatedModel(ws_forLimit.get());

        //
        // Gets the measurement object in the original combined workspace (created with the "w" command)
        //
        const std::string originalCombinedFile = fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root";
        std::unique_ptr<TFile> f_origin(TFile::Open(originalCombinedFile.c_str(), "read"));
        RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get((fInputName+fSuffix).c_str());
        TString outputName = f_origin->GetName();
        f_origin -> Close();

        //
        // Creating the rootfile used as input for the limit setting :-)
        //
        outputName = outputName.ReplaceAll(".root","_forLimits.root");
        std::unique_ptr<TFile> f_clone(TFile::Open( outputName, "recreate"));
        ws_forLimit -> import(*data,Rename("ttHFitterData"));
        originalMeasurement -> Write();
        ws_forLimit -> Write();
        f_clone -> Close();

        //
        // Finally computing the limit
        //
        std::string outputName_s = static_cast<std::string> (outputName);
        runAsymptoticsCLs(outputName_s.c_str(), "combined", "ModelConfig", "ttHFitterData", fLimitParamName.c_str(), fLimitParamValue, (fLimitOutputPrefixName+fSuffix).c_str(), (fName+"/Limits/").c_str(), fLimitIsBlind, fLimitsConfidence, "asimovData_0", fSignalInjection, fSignalInjectionValue, sigDebug);
    }
}

//__________________________________________________________________________________
//
void TRExFit::GetSignificance(){
    std::unique_ptr<RooDataSet> data(nullptr);
    
    //
    // Read NPvalues from fit-result file
    //
    if(fFitNPValuesFromFitResults!=""){
        WriteInfoStatus("TRExFit::Fit","Setting NP values for Asimov data-set creation from fit results stored in file " + fFitNPValuesFromFitResults + "...");
        fFitNPValues = NPValuesFromFitResults(fFitNPValuesFromFitResults);
    }

    //
    // If a workspace file name is specified, do simple significance
    //
    int sigDebug = 3 - TRExFitter::DEBUGLEVEL;
    if (sigDebug < 0) sigDebug = 0;
    if(fWorkspaceFileName!=""){
        std::string dataName = "obsData";
        if(fSignificanceIsBlind) dataName = "asimovData";
        runSig(fWorkspaceFileName.c_str(), "combined", "ModelConfig", dataName.c_str(), fSignificanceParamName.c_str(), fSignificanceParamValue, fSignificanceOutputPrefixName.c_str(), (fName+"/Significance").c_str(), fSignificanceIsBlind, "asimovData_1", "conditionalGlobs_1", "nominalGlobs", false, fSignificancePOIAsimov, sigDebug);
    }
    else{
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std::string > regionsForSign = ListRegionsToFit();
        std::map < std::string, int > regionsForSignDataType = MapRegionDataTypes(regionsForSign,fSignificanceIsBlind);
        std::map < std::string, double > npValues;
        //
        // flag if mixed Data / Asimov limit required
        // if mixed limit, perform a first fit on the regions with data only
        bool isMixedFit = DoingMixedFitting();
        if(isMixedFit && !fSignificanceIsBlind){
            std::vector < std:: string > regionsForFit = ListRegionsToFit(Region::REALDATA);
            //
            // Creates a combined workspace with the regions to be used *in the fit*
            //
            WriteInfoStatus("TRExFit::GetSignificance","Creating ws for regions with real data only...");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
            std::unique_ptr<RooWorkspace> ws_forFit = PerformWorkspaceCombination( regionsForFit );
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            if (!ws_forFit){
                WriteErrorStatus("TRExFit::GetSignificance","Cannot retrieve the workspace, exiting!");
                exit(EXIT_FAILURE);
            }
            //
            // Calls the PerformFit() function to actually do the fit
            //
            WriteInfoStatus("TRExFit::GetSignificance","Performing a fit in regions with real data only...");
            npValues = PerformFit( ws_forFit.get(), data.get(), FitType::BONLY, false, TRExFitter::DEBUGLEVEL);
            WriteInfoStatus("TRExFit::GetSignificance","Now will use the fit results to create the Asimov in the regions without real data!");
        }
        else{
            npValues = fFitNPValues; // FIXME: maybe should add a set of configurable NPs for significance only??
        }

        //
        // Create the final asimov dataset for limit setting
        //
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        std::unique_ptr<RooWorkspace> ws_forSignificance = PerformWorkspaceCombination( regionsForSign );
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        if (!ws_forSignificance){
            WriteErrorStatus("TRExFit::GetSignificance","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        data = std::unique_ptr<RooDataSet>(DumpData( ws_forSignificance.get(), regionsForSignDataType, npValues, npValues.find(fPOI)==npValues.end() ? fSignificancePOIAsimov : npValues[fPOI] ));
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        
        //
        // Set the global observables to the specified NPValues if the corresponding flag is TRUE
        //
        if(fInjectGlobalObservables && !fFitNPValues.empty()) {
            FitUtils::InjectGlobalObservables(ws_forSignificance.get(), fFitNPValues);
        }

        //
        // Set all saturated model factors to constant
        FitUtils::DisableSaturatedModel(ws_forSignificance.get());

        //
        // Gets the measurement object in the original combined workspace (created with the "w" command)
        //
        const std::string originalCombinedFile = fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root";
        std::unique_ptr<TFile> f_origin(TFile::Open(originalCombinedFile.c_str(), "read"));
        RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get((fInputName + fSuffix).c_str());
        TString outputName = f_origin->GetName();
        f_origin -> Close();

        //
        // Creating the rootfile used as input for the significance computation
        //
        outputName = outputName.ReplaceAll(".root","_forSignificance.root");
        std::unique_ptr<TFile> f_clone(TFile::Open( outputName, "recreate" ));
        ws_forSignificance -> import(*data,Rename("ttHFitterData"));
        originalMeasurement -> Write();
        ws_forSignificance -> Write();
        f_clone -> Close();

        //
        // Finally computing the significance
        //
        std::string outputName_s = static_cast<std::string> (outputName);
        runSig(outputName_s.c_str(), "combined", "ModelConfig", "ttHFitterData", fSignificanceParamName.c_str(), fSignificanceParamValue, fSignificanceOutputPrefixName.c_str(), (fName+"/Significance").c_str(), fSignificanceIsBlind, "asimovData_1", "conditionalGlobs_1", "nominalGlobs", false, fSignificancePOIAsimov, sigDebug);
    }
}

//__________________________________________________________________________________
//
void TRExFit::ReadFitResults(const std::string& fileName){
    WriteInfoStatus("TRExFit::ReadFitResults", "------------------------------------------------------");
    WriteInfoStatus("TRExFit::ReadFitResults",  "Reading fit results from file ");
    if(fFitResults) delete fFitResults;
    fFitResults = new FitResults();
    fFitResults->SetPOIPrecision(fPOIPrecision);

    if(fileName.find(".txt")!=std::string::npos)
    {
        fFitResults->ReadFromTXT(fileName, fBlindedParameters);
    }
    // make a list of systematics from all samples...
    // ...
    // assign to each NP in the FitResults a title, and a category according to the syst in the fitter

    // note: some NPs are assigned to multiple systematics (those which are correlated)
    // we will just keep overwriting, so the title and catagory
    // will be from the last systematic with that NP name

    for( unsigned int i_np=0; i_np < fFitResults->fNuisPar.size(); ++i_np )
    {

        for( unsigned int j_sys=0; j_sys < fSystematics.size(); ++j_sys )
        {
            // the systematic fName doesn't necessarily equate to the NP fName
            // compare the fSystematics[j_sys]->fNuisanceParameter instead!

            if( fSystematics[j_sys]->fNuisanceParameter == fFitResults->fNuisPar[i_np]->fName )
            {
                fFitResults->fNuisPar[i_np]->fTitle = fSystematics[j_sys]->fTitle;
                fFitResults->fNuisPar[i_np]->fCategory = fSystematics[j_sys]->fCategory;
            }
        }
        for(unsigned int j=0;j<fNormFactors.size();j++){
            if(fNormFactors[j]->fName == fFitResults->fNuisPar[i_np]->fName){
                fFitResults->fNuisPar[i_np]->fTitle = fNormFactors[j]->fTitle;
                fFitResults->fNuisPar[i_np]->fCategory = fNormFactors[j]->fCategory;
            }
        }
        // FIXME SF probably there are several NPs associated to it
        for(unsigned int j=0;j<fShapeFactors.size();j++){
            if(fShapeFactors[j]->fName == fFitResults->fNuisPar[i_np]->fName){
                fFitResults->fNuisPar[i_np]->fTitle = fShapeFactors[j]->fTitle;
                fFitResults->fNuisPar[i_np]->fCategory = fShapeFactors[j]->fCategory;
            }
        }
    }
    
    // create an ordered list of NPs
    if(fReorderNPs){
        fFitResults->fNuisParList.clear();
        for(auto norm : fNormFactors){
            if(norm->fExpression.first==""){
                if(Common::FindInStringVector(fFitResults->fNuisParList,norm->fNuisanceParameter)<0) fFitResults->fNuisParList.emplace_back(norm->fNuisanceParameter);
            }
        }
        for(auto syst : fSystematics){
            // if non-SHAPE
            if(syst->fType!=Systematic::SHAPE){
                if(Common::FindInStringVector(fFitResults->fNuisParList,syst->fNuisanceParameter)<0) fFitResults->fNuisParList.emplace_back(syst->fNuisanceParameter);
            }
            // FIXME: something to add for actual SHAPE systematics (not "separate gammas")
        }
        // FIXME: something to add for shape factors
    }
}

//__________________________________________________________________________________
//
void TRExFit::PrintConfigSummary() const{
    WriteInfoStatus("TRExFit::PrintConfigSummary", "-------------------------------------------");
    WriteInfoStatus("TRExFit::PrintConfigSummary", "Job name: "+fName);
    WriteInfoStatus("TRExFit::PrintConfigSummary", "Reading the following regions:");
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->Print();
    }
    WriteInfoStatus("TRExFit::PrintConfigSummary", "Reading the following samples:");
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        WriteInfoStatus("Sample::Print:","     "+fSamples[i_smp]->fName);
    }
    WriteInfoStatus("TRExFit::PrintConfigSummary", "Reading the following systematics:");
    std::vector<std::string> tmp{};
    for(const auto& isyst : fSystematics) {
        if (std::find(tmp.begin(), tmp.end(), isyst->fName) == tmp.end()){
            WriteInfoStatus("TRExFit::PrintConfigSummary"," "+isyst->fName);
            tmp.emplace_back(isyst->fName);
        }
    }
    WriteInfoStatus("TRExFit::PrintConfigSummary", "-------------------------------------------");
}

//__________________________________________________________________________________
//
Region* TRExFit::GetRegion(const std::string& name) const{
    for(unsigned int i=0;i<fRegions.size();i++){
        if(fRegions[i]->fName == name) return fRegions[i];
    }
    return nullptr;
}

//__________________________________________________________________________________
//
std::shared_ptr<Sample> TRExFit::GetSample(const std::string& name) const{
    for(unsigned int i=0;i<fSamples.size();i++){
        if(fSamples[i]->fName == name) return fSamples[i];
    }
    return nullptr;
}

//__________________________________________________________________________________
//
std::size_t TRExFit::GetSampleIndex(const std::string& name) const{
    for(std::size_t i=0; i<fSamples.size(); ++i){
        if(fSamples[i]->fName == name) return i;
    }
    return 99999;
}

//__________________________________________________________________________________
//
void TRExFit::DrawAndSaveSeparationPlots() const{

    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    gSystem->mkdir((fName+"/Plots/Separation").c_str());


    // loop over regions
    for(unsigned int i_ch=0; i_ch < fRegions.size(); i_ch++){
        // begin plotting
        TCanvas dummy3 ("dummy3", "dummy3", 600,600);
        dummy3.cd();

        if(fRegions[i_ch]->fNSig==0){
            WriteErrorStatus("TRExFit::DrawAndSaveSeparationPlots", "No Signal Found");
            continue;
        }

        std::unique_ptr<TH1D> sig(static_cast<TH1D*>(fRegions[i_ch]->fSig[0]->fHist->Clone()));

        std::unique_ptr<TH1D> bkg (static_cast<TH1D*>(fRegions[i_ch]->fBkg[0]->fHist->Clone())); // clone the first bkg
        for(int i_bkg=1; i_bkg< fRegions[i_ch] -> fNBkg; i_bkg++){
            bkg->Add(fRegions[i_ch]->fBkg[i_bkg]->fHist.get()); // add the rest
        }

        sig->SetLineColor( 2 );
        sig->SetLineWidth( 3 );
        sig->SetFillStyle( 0 );
        sig->SetLineStyle( 2 );

        bkg->SetLineColor( kBlue );
        bkg->SetLineWidth( 3 );
        bkg->SetFillStyle( 0 );
        bkg->SetLineStyle( 1 );

        TLegend legend3(0.55,0.77,0.94,0.87);
        legend3.SetTextFont(gStyle->GetTextFont());
        legend3.SetTextSize(gStyle->GetTextSize());
        legend3.AddEntry(bkg.get(), "Total background" , "l");
        legend3.AddEntry(sig.get(), fRegions[i_ch]->fSig[0]->fSample->fTitle.c_str() , "l");
        legend3.SetFillStyle(0) ;
        legend3.SetBorderSize(0);

        std::string xaxis = fRegions[i_ch]->fVariableTitle;

        sig->GetYaxis()->SetTitle("Arbitrary units");
        sig->GetXaxis()->SetTitle(xaxis.c_str());

        sig->GetYaxis()->SetTitleOffset(1.6);

        bkg->GetYaxis()->SetTitle("Arbitrary units");
        bkg->GetXaxis()->SetTitle(xaxis.c_str());

        bkg->GetYaxis()->SetTitleOffset(1.6);

        sig->GetYaxis()->SetNdivisions(506);
        bkg->GetYaxis()->SetNdivisions(506);

        sig->Scale(1./sig->Integral());
        bkg->Scale(1./bkg->Integral());


        if(bkg->GetMaximum() > sig->GetMaximum()){
            bkg->GetYaxis()->SetRangeUser(0.,bkg->GetMaximum()*1.5);
            bkg->Draw("hist");
            sig->Draw("histsame");
        }
        else {
            sig->GetYaxis()->SetRangeUser(0.,sig->GetMaximum()*1.5);
            sig->Draw("hist");
            bkg->Draw("histsame");
            sig->Draw("histsame");
        }

        legend3.Draw("same");

        myText(0.20,0.78,1,fLabel.c_str());
        myText(0.20,0.73,1,fRegions[i_ch]->fLabel.c_str());

        std::string cme = fRegions[i_ch]->fCmeLabel;
        std::string lumi = fRegions[i_ch]->fLumiLabel;

        myText(0.20,0.83,1,Form("#sqrt{s} = %s, %s", cme.c_str(), lumi.c_str()));

        if(fAtlasLabel!="none") ATLASLabelNew(0.20,0.84+0.04,(fAtlasLabel+"  Simulation").c_str(), kBlack, gStyle->GetTextSize());

        std::ostringstream SEP;
        SEP.precision(3);
        SEP << "Separation: " << Common::GetSeparation(sig.get(),bkg.get())*100 << "%";
        myText(0.55,0.73,1,SEP.str().c_str());

        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            dummy3.SaveAs((fName+"/Plots/Separation/"+fRegions[i_ch]->fName+fSuffix+"."+format).c_str());
        }

    }// regions

   return;
}

//____________________________________________________________________________________
//
void TRExFit::ProduceNPRanking(const std::string& NPnames) {

    if(fFitType==BONLY){
        WriteErrorStatus("TRExFit::ProduceNPRanking", "For ranking plots, the SPLUSB FitType is needed.");
        exit(EXIT_FAILURE);
    }

    //
    // List of systematics to check
    //
    RankingManager manager{};
    std::vector<std::string> systNames_unique;
    for(const auto& isyst : fSystematics) {
        if(NPnames=="all" || NPnames==isyst->fNuisanceParameter ) {
            if(isyst->fType == Systematic::SHAPE) continue;
            if (std::find(systNames_unique.begin(), systNames_unique.end(), isyst->fNuisanceParameter) == systNames_unique.end())
                systNames_unique.push_back(isyst->fNuisanceParameter);
            else continue;
            manager.AddNuisPar(isyst->fNuisanceParameter, false);
        }
    }
    for(const auto& inorm : fNormFactors) {
        if(fPOI == inorm->fName) continue;
        if(NPnames=="all" || NPnames==inorm->fName){
            manager.AddNuisPar(inorm->fName, true);
        }
    }

    //
    // Text files containing information necessary for drawing of ranking plot
    //
    std::string outName = fName+"/Fits/NPRanking"+fSuffix;
    if(fBootstrap!="" && fBootstrapIdx>=0){
        gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
        outName = fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+"NPRanking"+fSuffix;
    }
    if(NPnames!="all") outName += "_"+NPnames;
    outName += ".txt";
    manager.SetOutputPath(outName);

    //
    // Creating the combined model
    //
    std::unique_ptr<RooDataSet> data(nullptr);
    std::unique_ptr<RooWorkspace> ws(nullptr);
    std::unique_ptr<TFile> customWSfile(nullptr);
    
    //
    // Read NPvalues from fit-result file
    //
    if(fFitNPValuesFromFitResults!=""){
        WriteInfoStatus("TRExFit::Fit","Setting NP values for Asimov data-set creation from fit results stored in file " + fFitNPValuesFromFitResults + "...");
        fFitNPValues = NPValuesFromFitResults(fFitNPValuesFromFitResults);
    }
    
    //
    // If there's a workspace specified, go on with simple fit, without looking for separate workspaces per region
    //
    if(fWorkspaceFileName!=""){
        customWSfile.reset(TFile::Open(fWorkspaceFileName.c_str(),"read"));
        ws = std::unique_ptr<RooWorkspace>(static_cast<RooWorkspace*>(customWSfile->Get("combined")));
        if(!fFitIsBlind) data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("obsData")));
        else             data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("asimovData")));
    }
    else{
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std:: string > regionsToFit = ListRegionsToFit();
        std::map < std::string, int > regionDataType = MapRegionDataTypes(regionsToFit,fFitIsBlind);
        std::map < std::string, double > npValues;
        //
        // flag if mixed Data / Asimov fit required
        // if mixed fit, perform a first fit on the regions with data only
        bool isMixedFit = DoingMixedFitting();
        if(isMixedFit && !fFitIsBlind){
            WriteInfoStatus("TRExFit::ProduceNPRanking","");
            WriteInfoStatus("TRExFit::ProduceNPRanking","-------------------------------------------");
            WriteInfoStatus("TRExFit::ProduceNPRanking","Performing fit on regions with DataType = DATA to get NPs to inject in Asimov...");
            std::vector < std:: string > regionsForDataFit = ListRegionsToFit(Region::REALDATA);
            //
            // Creates a combined workspace with the regions to be used *in the data fit*
            //
            WriteInfoStatus("TRExFit::ProduceNPRanking","Creating ws for regions with real data only...");
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
            std::unique_ptr<RooWorkspace> ws_forFit = PerformWorkspaceCombination( regionsForDataFit );
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            if (!ws_forFit){
                WriteErrorStatus("TRExFit::ProduceNPRanking","Cannot retrieve the workspace, exiting!");
                exit(EXIT_FAILURE);
            }
            //
            // Calls the PerformFit() function to actually do the fit
            //
            WriteInfoStatus("TRExFit::ProduceNPRanking","Performing a B-only fit in regions with real data only...");
            npValues = PerformFit( ws_forFit.get(), data.get(), FitType::BONLY, false, TRExFitter::DEBUGLEVEL);
            WriteInfoStatus("TRExFit::ProduceNPRanking","Now will use the fit results to create the Asimov in the regions without real data!");
        }
        else{
            npValues = fFitNPValues;
        }
        //
        // Create the final asimov dataset for fit
        //
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        ws = PerformWorkspaceCombination( regionsToFit );
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        if (!ws){
            WriteErrorStatus("TRExFit::ProduceNPRanking","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npValues, npValues.find(fPOI)==npValues.end() ? fFitPOIAsimov : npValues[fPOI] ));
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    }

    RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("TRExFit::ProduceNPRanking", "Model config is nullptr");
        exit(EXIT_FAILURE);
    }

    // Loop on NPs to find gammas and add to the list to be ranked
    if(NPnames=="all" || NPnames.find("gamma")!=std::string::npos){
        RooRealVar* var(nullptr);
        const RooArgSet* nuis = static_cast<const RooArgSet*>(mc->GetNuisanceParameters());
        if(nuis){
            std::unique_ptr<TIterator> it2(nuis->createIterator());
            while( (var = static_cast<RooRealVar*>(it2->Next())) ){
                const std::string& np = var->GetName();
                if (np.find("saturated_model") != std::string::npos) continue;
                if(np.find("gamma")!=std::string::npos){
                    // add the nuisance parameter to the list nuisPars if it's there in the ws
                    // remove "gamma"...
                    if(np==NPnames || NPnames=="all"){
                        manager.AddNuisPar(Common::ReplaceString(np,"gamma_",""), true);
                        if(NPnames!="all") break;
                    }
                }
            }
        }
    }

    ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");

    manager.SetInjectGlobalObservables(fInjectGlobalObservables);
    manager.SetNPValues(fFitNPValues);
    manager.SetFixedNPs(fFitFixedNPs);
    manager.SetFitStrategy(fFitStrategy);
    manager.SetPOIName(fPOI);
    manager.SetNCPU(fCPU);
    manager.SetRng(fRndRange, fUseRnd, fRndSeed);
    manager.SetStatOnly(fStatOnly);

    manager.RunRanking(fFitResults, ws.get(), data.get(), fNormFactors);

    if(customWSfile!=nullptr) customWSfile->Close();
}

//____________________________________________________________________________________
//
void TRExFit::PlotNPRankingManager() const{
  if(fRankingPlot=="Merge"  || fRankingPlot=="all") PlotNPRanking(true,true);
  if(fRankingPlot=="Systs"  || fRankingPlot=="all") PlotNPRanking(true,false);
  if(fRankingPlot=="Gammas" || fRankingPlot=="all") PlotNPRanking(false,true);
}

//____________________________________________________________________________________
//
void TRExFit::PlotNPRanking(const bool flagSysts, const bool flagGammas) const{
    //
    const std::string fileToRead = fName+"/Fits/NPRanking"+fSuffix+".txt";

    std::ifstream in(fileToRead.c_str());
    if (!in.good()) { // file doesnt exist
        const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fName+"/Fits/","NPRanking"+fSuffix+"_");
        Common::MergeTxTFiles(inPaths, fileToRead);
    }

    RankingManager manager{};
    manager.SetOutputPath(fileToRead);
    manager.SetAtlasLabel(fAtlasLabel);
    manager.SetLumiLabel(fLumiLabel);
    manager.SetCmeLabel(fCmeLabel);
    manager.SetUseHEPDataFormat(fHEPDataFormat);
    manager.SetName(fName);
    manager.SetSuffix(fSuffix);
    manager.SetMaxNPPlot(fRankingMaxNP);
    manager.SetRankingPOIName(fRankingPOIName);
    manager.SetRankingCanvasSize(fNPRankingCanvasSize);

    manager.PlotRanking(fRegions, flagSysts, flagGammas);
}

//____________________________________________________________________________________
//
void TRExFit::PrintSystTables(std::string opt) const{
    WriteInfoStatus("TRExFit::PrintSystTables", "Printing syst tables");
    if(fCleanTables) opt += "clean";
    if(fSystCategoryTables) opt += "category";
    if(fTableOptions.find("STANDALONE")!=std::string::npos) opt += "standalone";
    if(fTableOptions.find("LANDSCAPE")!=std::string::npos) opt +="landscape";
    if(fTableOptions.find("FOOTNOTESIZE")!=std::string::npos) opt +="footnotesize";
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        fRegions[i_reg]->PrintSystTable(fFitResults,opt);
    }
}


//____________________________________________________________________________________
// this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter
void TRExFit::MergeSystematics(){
    // loop on systematics, see if any of them has name != nuisance
    for(const auto& syst : fSystematics){
        if(syst->fName == syst->fNuisanceParameter) continue;
        // if so, loop on other systematics to find one with name = nuisance parameter
        for(const auto& syst1 : fSystematics){
            if(syst->fName==syst1->fName) continue;
            if(!(syst->fNuisanceParameter==syst1->fName && syst1->fNuisanceParameter==syst1->fName)) continue;
            // now merge all SystematicHist in all regions
            WriteDebugStatus("TRExFit::MergeSystematics", "Found NP(syst) " + syst->fNuisanceParameter + "(" + syst->fName + ") = to syst name " + syst1->fName );
            for(const auto& reg : fRegions){
                WriteDebugStatus("TRExFit::MergeSystematics", "Region: " + reg->fName);
                for(const auto& sh : reg->fSampleHists){
                    std::shared_ptr<SystematicHist> syh  = sh->GetSystematic(syst ->fName);
                    std::shared_ptr<SystematicHist> syh1 = sh->GetSystematic(syst1->fName);
                    if(syh==nullptr || syh1 ==nullptr) continue;
                    // FIXME...
                    // the issue here is that to combine uncertainties one has to act differently depending on the fact that the different sources come from a multiplication/division or not...
                    syh1 ->Add(syh.get());
                    syh1 ->Add(sh->fHist.get(),-1);
                    WriteDebugStatus("TRExFit::MergeSystematics", "Adding syst of " + syh->fName +  " to " + syh1->fName);
                    WriteDebugStatus("TRExFit::MergeSystematics", "Setting to 0 all Up/Down of " +  syh->fName);
                    //
                    // set to zero the other syst
                    syh->fHistUp.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistUp->GetName())));
                    syh->fHistDown.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistDown->GetName())));
                    syh->fHistShapeUp.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistShapeUp->GetName())));
                    syh->fHistShapeDown.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistShapeDown->GetName())));
                    syh->fNormUp   = 0.;
                    syh->fNormDown = 0.;
                    syh->fNormPruned  = true;
                    syh->fShapePruned = true;
                }
            }
        }
    }
}

//____________________________________________________________________________________
// this will combine special systematics into a single systematic (e.g. envelope)
void TRExFit::CombineSpecialSystematics() {
    WriteInfoStatus("TRExFit::CombineSpecialSystematics", "Combining special systematics");
    std::vector<std::string> combineNames;

    // First we need to know how many special systematics there are
    for (const auto& syst : fSystematics) {
        if (syst->fCombineName == "") continue;
        if (std::find(combineNames.begin(), combineNames.end(), syst->fCombineName) == combineNames.end()) {
            combineNames.emplace_back(syst->fCombineName);
        }
    }

    if (combineNames.size() == 0) return;

    for (const auto& icombine : combineNames) {
        Systematic::COMBINATIONTYPE type(Systematic::COMBINATIONTYPE::ENVELOPE);
        // loop over regions
        std::vector<std::string> dropNorm;
        for (std::size_t iReg = 0; iReg < fRegions.size(); ++iReg) {
            // loop over systematics and get the list of systematics to combine
            std::vector<std::string> names;
            for (const auto& isyst : fSystematics) {
                if (isyst->fCombineName != icombine) continue;
                if (std::find(names.begin(), names.end(), isyst->fName) != names.end()) continue;
                names.emplace_back(isyst->fName);
                type = isyst->fCombineType;

                if (dropNorm.size() == 0) {
                    if (isyst->fDropNormSpecialIn.size() != 0) {
                        dropNorm = isyst->fDropNormSpecialIn;
                    }
                }
            }

            const bool dropFromAll = std::find(dropNorm.begin(), dropNorm.end(), "all") != dropNorm.end() ? true : false;

            // now loop over the systematics to combine
            std::vector<std::vector<std::shared_ptr<SystematicHist> > > sh_vec;

            for (const auto& sh : fRegions[iReg]->fSampleHists) {
                std::vector<std::shared_ptr<SystematicHist> > tmp;
                for (std::size_t ispecial = 0; ispecial < names.size(); ++ispecial) {
                    auto sampleHist = sh->GetSystematic(names.at(ispecial));
                    if (!sampleHist) continue;
                    tmp.emplace_back(sh->GetSystematic(names.at(ispecial)));
                }
                sh_vec.emplace_back(tmp);
            }

            if (sh_vec.size() == 0) {
                WriteWarningStatus("TRExFit::CombineSpecialSystematics", "No systematics to combine for " + icombine + " as no samples are affected");
                return;
            }

            // Now we have a list of systematicHist to combine and we change the first one, and then remove the rest
            for (std::size_t ish = 0; ish < fRegions[iReg]->fSampleHists.size(); ++ish) {
                auto newSampleHist = fRegions[iReg]->fSampleHists.at(ish)->GetSystematic(names.at(0));
                newSampleHist = CombineSpecialHistos(newSampleHist, sh_vec.at(ish), type, fRegions[iReg]->fSampleHists.at(ish).get());

                // apply the drop normalisation
                if (!newSampleHist) continue;
                if (dropFromAll || std::find(dropNorm.begin(), dropNorm.end(), fRegions[iReg]->fName) != dropNorm.end()) {
                    const SampleHist* sh = fRegions[iReg]->fSampleHists.at(ish).get();
                    if (!sh) continue;
                    const TH1* nominal = sh->fHist.get();
                    if (!nominal) continue;

                    newSampleHist->fHistUp  ->Scale(nominal->Integral()/newSampleHist->fHistUp  ->Integral());
                    newSampleHist->fHistDown->Scale(nominal->Integral()/newSampleHist->fHistDown->Integral());
                    newSampleHist->fNormPruned = true;
                    newSampleHist->fShapePruned = false;
                    newSampleHist->fNormUp = 0.;
                    newSampleHist->fNormDown = 0.;
                }
            }

            // And set the remaining systematicHist to zero
            for (std::size_t ispecial = 1; ispecial < names.size(); ++ispecial) {
                for (const auto& sh : fRegions[iReg]->fSampleHists) {
                    auto sampleHist = sh->GetSystematic(names.at(ispecial));
                    if (!sampleHist) continue;
                    sampleHist->fHistUp.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistUp->GetName())));
                    sampleHist->fHistDown.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistDown->GetName())));
                    sampleHist->fHistShapeUp.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistShapeUp->GetName())));
                    sampleHist->fHistShapeDown.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistShapeDown->GetName())));
                    sampleHist->fNormUp   = 0.;
                    sampleHist->fNormDown = 0.;
                    sampleHist->fNormPruned  = true;
                    sampleHist->fShapePruned = true;
                }
            }
        }
    }
}

//____________________________________________________________________________________
//
void TRExFit::ComputeBinning(int regIter){
    //
    //Creating histograms to rebin
    TH1D* hsig = nullptr;
    TH1D* hbkg = nullptr;
    bool nDefSig=true;
    bool nDefBkg=true;
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;
    bool bkgReg=false;
    bool flatBkg=false;
    if(fRegions[regIter]->fRegionType==Region::CONTROL) bkgReg=true;
    if(bkgReg && fRegions[regIter]->fTransfoDzSig<1e-3) flatBkg=true;
    //
    WriteDebugStatus("TRExFit::ComputeBinning", "Will compute binning with the following options:");

    std::string tmp_string = std::to_string(fRegions[regIter]->fTransfoFzSig);
    if((fRegions[regIter]->fTransfoDzSig>1e-3 || fRegions[regIter]->fTransfoDzBkg>1e-3) )
        WriteDebugStatus("TRExFit::ComputeBinning", " TransfoD - zSig=" + tmp_string + " - zBkg=" + std::to_string(fRegions[regIter]->fTransfoDzBkg));
    if((fRegions[regIter]->fTransfoFzSig>1e-3 || fRegions[regIter]->fTransfoFzBkg>1e-3) )
        WriteDebugStatus("TRExFit::ComputeBinning", " TransfoF - zSig=" + tmp_string + " - zBkg=" + std::to_string(fRegions[regIter]->fTransfoFzBkg));
    if((fRegions[regIter]->fTransfoJpar1>1e-3 || fRegions[regIter]->fTransfoJpar2>1e-3 || fRegions[regIter]->fTransfoJpar3>1e-3) )
        WriteDebugStatus("TRExFit::ComputeBinning", " TransfoJ - z1=" + std::to_string(fRegions[regIter]->fTransfoJpar1)
             + " - z2=" + std::to_string(fRegions[regIter]->fTransfoJpar2) + " - z3=" + std::to_string(fRegions[regIter]->fTransfoJpar3));

    if(bkgReg) WriteDebugStatus("TRExFit::ComputeBinning", " - bkg reg");
    else WriteDebugStatus("TRExFit::ComputeBinning", " - sig reg");

    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        //
        // using NTuples
        if(fInputType==1){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            if(fSamples[i_smp]->fType==Sample::GHOST) continue;
            if( Common::FindInStringVector(fSamples[i_smp]->fRegions,fRegions[regIter]->fName)<0 ) continue;
            //
            fullSelection = FullSelection(  fRegions[regIter],fSamples[i_smp].get());
            fullMCweight  = FullWeight(     fRegions[regIter],fSamples[i_smp].get());
            fullPaths     = FullNtuplePaths(fRegions[regIter],fSamples[i_smp].get());
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                int tmp_debugLevel=TRExFitter::DEBUGLEVEL;
                TRExFitter::SetDebugLevel(0);
                TH1D* htmp = nullptr;
                htmp = Common::HistFromNtuple( fullPaths[i_path],
                                    fRegions[regIter]->fVariable, 10000, fRegions[regIter]->fXmin, fRegions[regIter]->fXmax,
                                    fullSelection, fullMCweight, fDebugNev);
                TRExFitter::SetDebugLevel(tmp_debugLevel);
                //
                // Pre-processing of histograms (rebinning, lumi scaling)
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                //
                if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                //
                // Importing the histogram in TRExFitter
                if(fSamples[i_smp]->fType==Sample::SIGNAL){
                    if(nDefSig){
                        hsig = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                        nDefSig=false;
                    }
                    else hsig->Add(htmp);
                }
                else{
                    if(bkgReg && !flatBkg){
                        bool usedInSig=false;
                        for(unsigned int i_bkgs=0; i_bkgs<fRegions[regIter]->fAutoBinBkgsInSig.size(); ++i_bkgs){
                            if(fSamples[i_smp]->fName==fRegions[regIter]->fAutoBinBkgsInSig[i_bkgs]){
                                usedInSig=true;
                                break;
                            }
                        }
                        if(usedInSig){
                            WriteDebugStatus("TRExFit::ComputeBinning", "Using " + fSamples[i_smp]->fName + " as signal");
                            if(nDefSig){
                                hsig = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                                nDefSig=false;
                            }
                            else hsig->Add(htmp);
                        }
                        else{
                            if(nDefBkg){
                                hbkg = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                                nDefBkg=false;
                            }
                            else hbkg->Add(htmp);
                        }
                    }
                    else{
                        if(nDefBkg){
                            hbkg = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp);
                    }
                }
                delete htmp;
            }
        }
        //
        // Input with hists
        else if(fInputType == 0){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            if(fSamples[i_smp]->fType==Sample::GHOST) continue;
            //
            fullPaths     = FullHistogramPaths(fRegions[regIter],fSamples[i_smp].get());
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                int tmp_debugLevel=TRExFitter::DEBUGLEVEL;
                TRExFitter::SetDebugLevel(0);
                std::unique_ptr<TH1> htmp = Common::HistFromFile( fullPaths[i_path] );
                if (!htmp) {
                    WriteErrorStatus("TRExFit::ReadHistograms", "Histo pointer is empty cannot continue running the code");
                    exit(EXIT_FAILURE);
                }
                TRExFitter::SetDebugLevel(tmp_debugLevel);
                //
                // Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[regIter]->fHistoBins.size() > 0){
                    const std::string hname = htmp->GetName();
                    std::unique_ptr<TH1> tmp_copy(static_cast<TH1*>(htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin, "tmp_copy", &fRegions[regIter]->fHistoBins[0])));
                    htmp.reset(tmp_copy.release());
                    htmp->SetName(hname.c_str());
                    if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(htmp.get());
                }
                else if(fRegions[regIter]->fHistoNBinsRebin != -1) {
                    htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin);
                }
                //
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                //
                if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                //
                // apply histogram to signal or background
                if(fSamples[i_smp]->fType==Sample::SIGNAL){
                    if(nDefSig){
                        hsig = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                        nDefSig=false;
                    }
                    else hsig->Add(htmp.get());
                }
                else{
                if(bkgReg && !flatBkg){
                    bool usedInSig=false;
                    for(unsigned int i_bkgs=0; i_bkgs<fRegions[regIter]->fAutoBinBkgsInSig.size(); ++i_bkgs){
                        if(fSamples[i_smp]->fName==fRegions[regIter]->fAutoBinBkgsInSig[i_bkgs]){
                            usedInSig=true;
                            break;
                        }
                    }
                    if(usedInSig){
                        WriteDebugStatus("TRExFit::ComputeBinning", "Using " + fSamples[i_smp]->fName + " as signal");
                        if(nDefSig){
                            hsig = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                            nDefSig=false;
                        }
                        else hsig->Add(htmp.get());
                    }
                    else{
                        if(nDefBkg){
                            hbkg = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp.get());
                    }
                }
                    else{
                        if(nDefBkg){
                            hbkg = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp.get());
                    }
                }
            }
        }
    }
    //
    //computing new bins
    //
    // get a vector of bins where to rebin to get an uncertainty <= 5% per bin.
    // starting from highest bin!
    // the numbers give the lowest bin included in the new bin
    // overflowbin+1 and underflow bins are returned as the first and last element in the vector, respectively.
    std::vector<int> bins_vec;
    //
    if (!hbkg || !hsig) {
        WriteErrorStatus("TRExFit::ComputeBinning", "Please provide signal and background histograms!");
        gSystem -> Exit(1);
    }
    int iBin2 = 0;
    int nBins_Vec = hbkg -> GetNbinsX();
    int iBin = nBins_Vec; // skip overflow bin
    bins_vec.push_back(nBins_Vec + 1);
    double nBkg = hbkg -> Integral(1, nBins_Vec );
    double nSig = hsig -> Integral(1, nBins_Vec );
    bool jBreak = false;
    double jTarget = 1e5;
    double jGoal = fRegions[regIter]->fTransfoJpar1;
    //
    while (iBin > 0) {
        double sumBkg = 0;
        double sumSig = 0;
        double sumSigL = 0;
        double err2Bkg = 0;
        bool pass = false;
        int binCount = 1;
        double dist = 1e10;
        double distPrev = 1e10;
        //
        while (!pass && iBin > 0) {
            double nBkgBin = hbkg -> GetBinContent(iBin);
            double nSigBin = hsig -> GetBinContent(iBin);
            sumBkg += nBkgBin;
            sumSig += nSigBin;
            if (nBkgBin > 0 && nSigBin > 0) {
                sumSigL += nSigBin * std::log1p(nSigBin / nBkgBin);
            }
            err2Bkg += (hbkg -> GetBinError(iBin))*(hbkg -> GetBinError(iBin));
            //
            double err2RelBkg = 1;
            if (sumBkg != 0) {
                err2RelBkg = err2Bkg / (sumBkg*sumBkg);
            }
            //
            double err2Rel = 1.;
            if(fRegions[regIter]->fBinTransfo == "TransfoD"){
                // "trafo D"
                if (sumBkg != 0 && sumSig != 0)
                  err2Rel = 1. / (sumBkg / (nBkg / fRegions[regIter]->fTransfoDzBkg) + sumSig / (nSig / fRegions[regIter]->fTransfoDzSig));
                else if (sumBkg != 0)
                  err2Rel = (nBkg / fRegions[regIter]->fTransfoDzBkg) / sumBkg;
                else if (sumSig != 0)
                  err2Rel = (nSig / fRegions[regIter]->fTransfoDzSig) / sumSig;
                pass = std::sqrt(err2Rel) < 1;
                // distance
                dist = std::fabs(err2Rel - 1);
            }
            else if(fRegions[regIter]->fBinTransfo == "TransfoF"){
                // "trafo F" with 5% bkg stat unc
                if (sumBkg != 0 && sumSigL != 0)
                  err2Rel = 1 / (std::sqrt(sumBkg / (nBkg / fRegions[regIter]->fTransfoFzBkg)) + std::sqrt(sumSigL / (1 / fRegions[regIter]->fTransfoFzSig)));
                else if (sumBkg != 0)
                  err2Rel = std::sqrt((nBkg / fRegions[regIter]->fTransfoFzBkg) / sumBkg);
                else if (sumSigL != 0)
                  err2Rel = std::sqrt((1 / fRegions[regIter]->fTransfoFzSig) / sumSigL);
                pass = std::sqrt(err2Rel) < 1 && std::sqrt(err2RelBkg) < 0.10;
            }
            else if(fRegions[regIter]->fBinTransfo == "TransfoJ"){
                if (!jBreak) pass = (sumBkg >  jGoal);
                else pass = (sumBkg > jTarget);
                if( pass && !jBreak ){
                    if( (sumSig/sumBkg) <  fRegions[regIter]->fTransfoJpar2*(nSig/nBkg) ){
                        jBreak = true;
                        jTarget = hbkg->Integral(0,iBin)/ fRegions[regIter]->fTransfoJpar3;
                    }
                    else{
                        jGoal = jGoal+1;
                    }
                }
            }
            else{
                WriteErrorStatus("TRExFit::ComputeBinning", "transformation method '" + fRegions[regIter]->fBinTransfo + "' unknown, try again!");
                exit(1);
            }
            if (!(pass && dist > distPrev)) {
                binCount++;
                iBin--;
            } // else use previous bin
            distPrev = dist;
        }
        iBin2++;
        // remove last bin
        if (iBin == 0 && bins_vec.size() > 1) {
            if (fRegions[regIter]->fBinTransfo == "TransfoF") {
                bins_vec.pop_back();
            }
            else if (fRegions[regIter]->fBinTransfo == "TransfoD" && bins_vec.size() > fRegions[regIter]->fTransfoDzSig + fRegions[regIter]->fTransfoDzBkg + 0.01) {
                // remove last bin if Nbin > Zsig + Zbkg
                // (1% threshold to capture rounding issues)
                bins_vec.pop_back();
            }
        }
        bins_vec.push_back(iBin + 1);
    }
    //
    //transform bin numbers in histo edges
    int nBins = bins_vec.size();
    double *bins = new double[nBins];
    bins[0] = hbkg->GetBinLowEdge(1);
    for(unsigned int i=1; i<bins_vec.size()-1; ++i){
        bins[i] = hbkg->GetBinLowEdge(bins_vec[nBins-i-1]);
    }
    bins[nBins-1]=hbkg->GetBinLowEdge( hbkg->GetNbinsX() + 1 );
    WriteInfoStatus("TRExFit::ComputeBinning", "Your final binning from automatic binning function is:");
    std::string temp_string = "";
    for(unsigned int i_bins=0; i_bins<bins_vec.size(); ++i_bins){
      temp_string+= std::to_string(bins[i_bins]) + " - ";
    }
    WriteInfoStatus("TRExFit::ComputeBinning", "  " + temp_string);
    //
    delete hsig;
    delete hbkg;
    fRegions[regIter]->SetBinning(nBins-1, bins);
    delete[] bins;
}

//__________________________________________________________________________________
//
void TRExFit::GetLikelihoodScan(RooWorkspace *ws,
                                const std::string& varName,
                                RooDataSet* data) const {
    WriteInfoStatus("TRExFit::GetLikelihoodScan", "Running likelihood scan for the parameter = " + varName);

    LikelihoodScanManager manager{};
    manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps);
    manager.SetNCPU(fCPU);
    manager.SetOffSet(true);
    manager.SetBlindedParameters(fBlindedParameters);
    
    auto scanResult = manager.Run1DScan(ws, varName, data);

    const std::vector<double> x = scanResult.first;
    const std::vector<double> y = scanResult.second;

    if (x.empty() || y.empty() ) {
        WriteWarningStatus("TRExFit::GetLikelihoodScan", "Skipping LHscan");
        return;
    }

    TCanvas can("NLLscan");

    TGraph graph(fLHscanSteps, &x[0], &y[0]);
    graph.Draw("ALP");
    graph.GetXaxis()->SetRangeUser(x[0],x.back());

    // y axis
    graph.GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());

    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(varName);

    can.SetTitle(cname);
    can.SetName(cname);
    can.cd();

    TLatex tex{};
    tex.SetTextColor(kGray+2);

    TLine l1s(x[0],0.5,x.back(),0.5);
    l1s.SetLineStyle(kDashed);
    l1s.SetLineColor(kGray);
    l1s.SetLineWidth(2);
    if(graph.GetMaximum()>2){
        l1s.Draw();
        tex.DrawLatex(x.back(),0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }

    if(graph.GetMaximum()>2){
        TLine l2s(x[0],2,x.back(),2);
        l2s.SetLineStyle(kDashed);
        l2s.SetLineColor(kGray);
        l2s.SetLineWidth(2);
        l2s.Draw();
        tex.DrawLatex(x.back(),2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
    }
    //
    if(graph.GetMaximum()>4.5){
        TLine l3s(x[0],4.5,x.back(),4.5);
        l3s.SetLineStyle(kDashed);
        l3s.SetLineColor(kGray);
        l3s.SetLineWidth(2);
        l3s.Draw();
        tex.DrawLatex(x.back(),4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
    }
    //
    TLine lv0(0,graph.GetMinimum(),0,graph.GetMaximum());
    lv0.Draw();
    //
    TLine lh0(x[0],0,x.back(),0);
    lh0.Draw();

    TString LHDir("LHoodPlots/");
    system(TString("mkdir -vp ")+fName+"/"+LHDir);

    can.RedrawAxis();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        can.SaveAs(fName+"/"+LHDir+"NLLscan_"+varName+fSuffix+"."+format);
    }

    // write it to a ROOT file as well
    std::unique_ptr<TFile> f(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varName+fSuffix+"_curve.root","UPDATE"));
    if (!f) {
        WriteWarningStatus("TRExFit::GetLikelihoodScan", "Cannot open ROOT file for likelihood scan!");
        return;
    }
    f->cd();
    graph.Write("LHscan",TObject::kOverwrite);
    f->Close();
}

//____________________________________________________________________________________
//
void TRExFit::Get2DLikelihoodScan( RooWorkspace *ws, const std::vector<std::string>& varNames, RooDataSet* data) const{
    if (varNames.size() != 2){
        WriteErrorStatus("TRExFit::Get2DLikelihoodScan", "Wrong number of parameters provided for 2D likelihood scan, returning");
        return;
    }
    
    LikelihoodScanManager manager{};
    manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps);
    manager.SetScanParamsY(fLHscanMinY, fLHscanMaxY, fLHscanStepsY);
    manager.SetNCPU(fCPU);
    manager.SetOffSet(!fParal2D);
    manager.SetBlindedParameters(fBlindedParameters);
    if (fParal2D) {
        manager.SetParalel2Dstep(fParal2Dstep);
    }

    const std::pair<std::string, std::string> names = std::make_pair(varNames.at(0), varNames.at(1));
    auto scanResult = manager.Run2DScan(ws, names, data);

    const std::vector<double> x = scanResult.x;
    const std::vector<double> y = scanResult.y;
    std::vector<std::vector<double> > z = scanResult.z;

    if (x.empty() || y.empty() || z.empty()) {
        WriteWarningStatus("TRExFit::Get2DLikelihoodScan", "Skipping LHscan");
        return;
    }

    // find minimium in z
    auto findMin = [](const std::vector<std::vector<double> >& vec) {
        double min(9999999);
        for (const auto& i : vec) {
            for (const auto& j : i) {
                if (j < min) min = j;
            }
        }
        return min;
    };

    const double zmin = findMin(z);

    // make plots
    TCanvas can("NLLscan_2D_");
    can.cd();

    TGraph2D graph(fLHscanSteps * fLHscanStepsY);

    TH2D h_nll("NLL", "NLL", fLHscanSteps, x[0], x.back(), fLHscanStepsY, y[0], y.back());
    unsigned int i=0;
    for (int ipoint = 0; ipoint < fLHscanSteps; ++ipoint) {
        if (fParal2D && ipoint!=fParal2Dstep) // if you are parallelizing, only run the point corresponding to the one passed from command line
            continue;
        for (int jpoint = 0; jpoint < fLHscanStepsY; ++jpoint) {
            if (!fParal2D) { // if you are paralellizing, no knowledge of the absolute minimum in each job
                // shift the likelihood values to zero
                z[ipoint][jpoint] -= zmin;
            }
            h_nll.SetBinContent(ipoint+1, jpoint+1, z[ipoint][jpoint]);
            i = ipoint * fLHscanStepsY + jpoint;
            graph.SetPoint(i,x[ipoint],y[jpoint],z[ipoint][jpoint]);
        }
    }

    TString LHDir("LHoodPlots/");
    system(TString("mkdir -vp ")+fName+"/"+LHDir);

    if (!fParal2D) { // Only draw and save graph when not running parallel
        gStyle->SetPalette(57); // Reset Palette to default (Pruning or Correlation matrinx changes this)
        graph.Draw("colz");
        graph.GetXaxis()->SetRangeUser(x[0],x.back());
        graph.GetYaxis()->SetRangeUser(y[0],y.back());

        // y axis
        graph.GetXaxis()->SetTitle(varNames.at(0).c_str());
        graph.GetYaxis()->SetTitle(varNames.at(1).c_str());

        // Print the canvas
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            can.SaveAs(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"."+format);
        }

        // write it to a ROOT file as well
        std::unique_ptr<TFile> f(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"_curve.root","UPDATE"));
        f->cd();
        graph.Write(("LHscan_2D_"+varNames.at(0)+"_"+varNames.at(1)).c_str(),TObject::kOverwrite);
        f->Close();
    }

    // Write histogram to Root file as well
    if (fParal2D) {
        std::ostringstream step_os;
        step_os << fParal2Dstep;
        std::string paral2Dstep_str=step_os.str();
        std::unique_ptr<TFile> f2(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"_step"+paral2Dstep_str+fSuffix+"_histo.root","UPDATE"));
        h_nll.Write("NLL",TObject::kOverwrite);
        f2->Close();
    } else {
        std::unique_ptr<TFile> f2(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"_histo.root","UPDATE"));
        h_nll.Write("NLL",TObject::kOverwrite);
        f2->Close();
    }
}

//__________________________________________________________________________________
//
void TRExFit::AddTemplateWeight(const std::string& name, double value){
    std::pair<double, std::string> temp = std::make_pair(value, name);
    fTemplatePair.push_back(temp);
}

//__________________________________________________________________________________
//
std::vector<TRExFit::TemplateWeight> TRExFit::GetTemplateWeightVec(const TRExFit::TemplateInterpolationOption& opt){
    std::vector<TRExFit::TemplateWeight> vec;
    for(auto name : fMorphParams){
        // create map only for values of the specified parameter
        std::vector<std::pair<double,std::string> > templatePair; templatePair.clear();
        for(auto tp : fTemplatePair){
            if(tp.second==name) templatePair.push_back(tp);
        }
        // first sort vector of inputs for templates
        if (templatePair.size() < 2){
            WriteErrorStatus("TRExFit::GetTemplateWeightVec", "You need to provide at least 2 templates for template fit to work, but you provided: " + std::to_string(fTemplatePair.size()));
        }
        std::sort(templatePair.begin(), templatePair.end());
        // find min and max for range
        double min = templatePair.at(0).first;
        double max = templatePair.at(templatePair.size() -1).first;
        for (unsigned int itemp = 0; itemp < (templatePair.size() ); itemp++){
            WriteDebugStatus("TRExFit::GetTemplateWeightVec", "Morphing: Template " + std::to_string(itemp));
            TRExFit::TemplateWeight tmp;
            tmp.name = templatePair.at(itemp).second;
            tmp.value = templatePair.at(itemp).first;
            WriteDebugStatus("TRExFit::GetTemplateWeightVec", "Morphing:   " + tmp.name + " = " + std::to_string(tmp.value));
            tmp.range = tmp.name+"["+std::to_string(min)+","+std::to_string(min)+","+std::to_string(max)+"]";
            // calculate the actual function
            tmp.function = TRExFit::GetWeightFunction(templatePair, itemp, opt);
            vec.push_back(tmp);
        }
    }
    return vec;
}

//__________________________________________________________________________________
//
std::string TRExFit::GetWeightFunction(std::vector<std::pair<double,std::string> > templatePair, unsigned int itemp, const TRExFit::TemplateInterpolationOption& opt) const{
    std::string fun = "";
    double x_i;
    std::string name;
    if (itemp < templatePair.size()){
        x_i = templatePair.at(itemp).first;
        name = templatePair.at(itemp).second;
    }
    else return fun;
    //
    if (opt == TRExFit::LINEAR){
        double deltaXp(-1.); // |x(i+1)-x(i)|
        double deltaXm(-1.); // |x(i-1)-x(i)|
        if ((itemp+1) < templatePair.size() ){
            deltaXp = std::fabs(templatePair.at(itemp+1).first - templatePair.at(itemp).first);
        }
        if (((int)itemp-1) >=0 ){
            deltaXm = std::fabs(templatePair.at(itemp-1).first - templatePair.at(itemp).first);
        }
        if(deltaXp<0 && deltaXm<0){
            WriteErrorStatus("TRExFit::GetWeightFunction", "Morphing: delta X = " + std::to_string(deltaXp) + ", " + std::to_string(deltaXm));
            return fun;
        }
        fun  = "(";
        if(deltaXm>0) fun += "((("+name+"-"+std::to_string(x_i)+")< 0)&&(fabs("+name+"-"+std::to_string(x_i)+")<"+std::to_string(deltaXm)+"))*(1.-(fabs("+name+"-"+std::to_string(x_i)+"))/"+std::to_string(deltaXm)+")";
        else fun += "0.";
        fun += "+";
        if(deltaXp>0) fun += "((("+name+"-"+std::to_string(x_i)+")>=0)&&(fabs("+name+"-"+std::to_string(x_i)+")<"+std::to_string(deltaXp)+"))*(1.-(fabs("+name+"-"+std::to_string(x_i)+"))/"+std::to_string(deltaXp)+")";
        else fun += "0.";
        fun += ")";
    } else if (opt == TRExFit::SMOOTHLINEAR) {
        // this will return a string that represents integral of hyperbolic tangent function that
        // approximates absolute value
        fun = GetSmoothLinearInterpolation(itemp);
    } else if (opt == TRExFit::SQUAREROOT) {
        fun = GetSquareRootLinearInterpolation(itemp);
    }
    // ...
    fun = Common::ReplaceString(fun,"--","+");
    WriteDebugStatus("TRExFit::GetWeightFunction", "Morphing:   weight function = " + fun);
    return fun;
}

//__________________________________________________________________________________
//
std::string TRExFit::GetSmoothLinearInterpolation(unsigned int itemp) const {
    // The idea is simple: use integral of hyperbolic tangent
    //
    // -2/width*(1/k)*corr*(log(e^(k*(x-x_mean)+e^(-(k*(x-x_mean))))) - ln(2)) + 1
    //
    // but the function is split into two parts to allow also non-equal steps in
    // the templates used for the template fit
    // "width" represents the x-axis size that is used for that particular function
    // "corr" represents correction to the function so that for x_min/x_max the functional value
    // is 0, without this correction it won't exactly be zero, because it is an anproximation

    if (itemp >= fTemplatePair.size()){
        return "";
    }

    // parameter that controls how close to a linear function we want to be
    static const double k_init(80.);

    double x_left = -99999.;
    double x_right = -99999.;
    double corr_left = 1.;
    double corr_right = 1.;

    if (itemp == 0) { // first template
        x_left = 2*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp+1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    } else if (itemp == (fTemplatePair.size()-1)) { // last template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = 2*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp-1).first;
    } else { // general template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    }

    double x_mean = fTemplatePair.at(itemp).first;
    double width_left = 2*std::fabs(x_mean - x_left);
    double width_right = 2*std::fabs(x_mean - x_right);

    // apply correction to the k parameter depending on the range of the x axis
    const double k_left = k_init/width_left;
    const double k_right = k_init/width_right;

    // calculate correction to the function to get y= 0 at x_min and x_max
    // use iterative process to find something which is close enough
    corr_left= 1+GetCorrection(k_left, width_left, x_mean, x_left);
    corr_left+= GetCorrection(k_left, width_left, x_mean, x_left, corr_left);
    corr_left+= GetCorrection(k_left, width_left, x_mean, x_left, corr_left);
    corr_right= 1+GetCorrection(k_right, width_right, x_mean, x_right);
    corr_right+= GetCorrection(k_right, width_right, x_mean, x_right, corr_right);
    corr_right+= GetCorrection(k_right, width_right, x_mean, x_right, corr_right);

    // prepare the actual string as "function" + "step function"
    std::string name = fTemplatePair.at(itemp).second;
    std::string step_left = "";
    std::string step_right = "";

    // the step function
    if (itemp == 0) {
        step_left = "(("+name+"-"+std::to_string(x_mean)+"<0)&&("+name+"-"+std::to_string(x_mean)+">0))";
        step_right = "(("+name+"-"+std::to_string(x_mean)+">=0) && ("+name+"<"+std::to_string(x_right)+"))";
    } else if (itemp == (fTemplatePair.size()-1)) {
        step_left = "(("+name+">="+std::to_string(x_left)+")&&("+name+"<"+std::to_string(x_mean)+"))";
        step_right = "(("+name+"-"+std::to_string(x_mean)+"<0)&&("+name+">"+std::to_string(x_right)+"))";
    } else {
        step_left = "((("+name+"-"+std::to_string(x_mean)+")<=0)&&("+name+">"+std::to_string(x_left)+"))";
        step_right = "((("+name+"-"+std::to_string(x_mean)+")>0)&&("+name+"<"+std::to_string(x_right)+"))";
    }

    // the functions
    std::string fun_left = "(-2/("+std::to_string(width_left*k_left)+")*"+std::to_string(corr_left)+"*(log(exp("+std::to_string(k_left)+"*("+name+"-"+std::to_string(x_mean)+")) + exp(-"+std::to_string(k_left)+"*("+name+"-"+std::to_string(x_mean)+")"+")) - log(2)) +1) * " + step_left;

    std::string fun_right = "(-2/("+std::to_string(width_right*k_right)+")*"+std::to_string(corr_right)+"*(log(exp("+std::to_string(k_right)+"*("+name+"-"+std::to_string(x_mean)+")) + exp(-"+std::to_string(k_right)+"*("+name+"-"+std::to_string(x_mean)+")"+")) - log(2)) +1) * " + step_right;

    return ("(("+fun_left+") + (" +fun_right+"))");
}

//__________________________________________________________________________________
//
double TRExFit::GetCorrection(double k, double width, double x_mean, double x_left, double init) const {
    double logterm = 0;
    double corr = 0;

    // since we are calculating logarithm of potentially wery large numbers (e^20)
    // we need to help the code in case we get overflow, when this happens we simply discard
    // the smaller contribution and set log(e^(x)) = x manually;

    //check if the number is inf
    if (exp(k*(x_left-x_mean)) > exp(k*(x_left-x_mean))) {
        logterm = k*(x_left-x_mean);
    } else if (exp(-k*(x_left-x_mean)) > exp(-k*(x_left-x_mean))) {
        logterm = -k*(x_left-x_mean);
    } else {
        logterm = log (exp(k*(x_left-x_mean)) + exp(-k*(x_left-x_mean)));
    }
    corr = ((-2/(k*width))*init*(logterm - log(2))+1);

    return corr;
}

//__________________________________________________________________________________
//
std::string TRExFit::GetSquareRootLinearInterpolation(unsigned int itemp) const {
    double epsilon = 0.0000001;

    double x_i = fTemplatePair.at(itemp).first;
    double x_left = -99999.;
    double x_right = -99999.;

    if (itemp == 0) { // first template
        x_left = 2.*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp+1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    } else if (itemp == (fTemplatePair.size()-1)) { // last template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = 2.*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp-1).first;
    } else { // general template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    }

    //apply correction
    double a_left = 0;
    double b_left = 0;
    double a_right = 0;
    double b_right = 0;

    GetSquareCorrection(&a_left, &b_left, x_i, x_left, epsilon);
    GetSquareCorrection(&a_right, &b_right, x_i, x_right, epsilon);

    // prepare the actual string as "function" + "step function"
    std::string name = fTemplatePair.at(itemp).second;
    std::string step_left = "";
    std::string step_right = "";

    // the step function
    if (itemp == 0) {
        step_left = "(("+name+"-"+std::to_string(x_i)+"<0)&&("+name+"-"+std::to_string(x_i)+">0))";
        step_right = "(("+name+"-"+std::to_string(x_i)+">=0) && ("+name+"<"+std::to_string(x_right)+"))";
    } else if (itemp == (fTemplatePair.size()-1)) {
        step_left = "(("+name+">="+std::to_string(x_left)+")&&("+name+"<"+std::to_string(x_i)+"))";
        step_right = "(("+name+"-"+std::to_string(x_i)+"<0)&&("+name+">"+std::to_string(x_right)+"))";
    } else {
        step_left = "((("+name+"-"+std::to_string(x_i)+")<=0)&&("+name+">"+std::to_string(x_left)+"))";
        step_right = "((("+name+"-"+std::to_string(x_i)+")>0)&&("+name+"<"+std::to_string(x_right)+"))";
    }

    std::string fun_left = "("+step_left+")*(-"+std::to_string(a_left)+"*sqrt(("+name+"-"+std::to_string(x_i)+")*("+name+"-"+std::to_string(x_i)+")+"+std::to_string(epsilon)+")+"+std::to_string(b_left)+")";
    std::string fun_right = "("+step_right+")*(-"+std::to_string(a_right)+"*sqrt(("+name+"-"+std::to_string(x_i)+")*("+name+"-"+std::to_string(x_i)+")+"+std::to_string(epsilon)+")+"+std::to_string(b_right)+")";

    return ("("+fun_left+"+"+fun_right+")");
}

//__________________________________________________________________________________
//
void TRExFit::GetSquareCorrection(double *a, double *b, double x_i, double x_left, double epsilon) const {
    if (x_left == 0) {
        x_left = 2*x_i;
    }

    // this can be analytically calculated
    double k = std::sqrt(((x_i - x_left)*(x_i - x_left)/epsilon) + 1);
    *b = k/(k-1);
    *a = (*b ) / std::sqrt((x_i-x_left)*(x_i-x_left) + epsilon);
}

//__________________________________________________________________________________
//
void TRExFit::SmoothMorphTemplates(const std::string& name,const std::string& formula,double *p) const{
    TCanvas c("c","c",600,600);
    // find NF associated to this morph param
    std::shared_ptr<NormFactor> nf = nullptr;
    for(const auto& norm : fNormFactors){
        if(norm->fName == name) nf = norm;
    }
    // get one histogram per bin (per region)
    for(auto reg : fRegions){
        std::map<double,TH1*> hMap; // map (paramater-value,histogram)
        TH1* h_tmp = nullptr;
        int nTemplates = 0;
        double min = -999.;
        double max = -999.;
        for(auto& sh : reg->fSampleHists){
            Sample* smp = sh->fSample;
            // if the sample has morphing
            if(smp->fIsMorph[name]){
                hMap[smp->fMorphValue[name]] = sh->fHist.get();
                if(smp->fMorphValue[name]<min) min = smp->fMorphValue[name];
                if(smp->fMorphValue[name]>max) max = smp->fMorphValue[name];
                nTemplates++;
                h_tmp = sh->fHist.get();
            }
        }
        if(h_tmp==nullptr) return;
        for(int i_bin=1;i_bin<=h_tmp->GetNbinsX();i_bin++){
            TGraphErrors g_bin(nTemplates);
            int i_pt = 0;
            for(auto vh : hMap){
                g_bin.SetPoint(i_pt,vh.first,vh.second->GetBinContent(i_bin));
                g_bin.SetPointError(i_pt,0,vh.second->GetBinError(i_bin));
                // if it's the nominal sample, set error to very small value => forced not to change nominal!
                if(nf!=nullptr) if(nf->fNominal==vh.first) g_bin.SetPointError(i_pt,0,vh.second->GetBinError(i_bin)*0.001);
                i_pt++;
            }
            c.cd();
            g_bin.Draw("epa");
            TF1 l("l",formula.c_str(),min,max);
            if(p!=0x0) l.SetParameters(p);
            g_bin.Fit("l","RQN");
            l.SetLineColor(kRed);
            l.Draw("same");
            gSystem->mkdir((fName+"/Morphing/").c_str());
            for(const auto& format : TRExFitter::IMAGEFORMAT) {
                c.SaveAs((fName+"/Morphing/g_"+name+"_"+reg->fName+"_bin"+std::to_string(i_bin)+"."+format).c_str());
            }
            for(auto vh : hMap){
                vh.second->SetBinContent(i_bin,l.Eval(vh.first));
            }
        }
    }
}

//____________________________________________________________________________________
//
bool TRExFit::MorphIsAlreadyPresent(const std::string& name, const double value) const {
    for (const std::pair<double, std::string> itemp : fTemplatePair){
        if ((itemp.second == name) && (itemp.first == value)){
            return true;
        }
    }
    return false;
}

//____________________________________________________________________________________
// create a map associating parameters to their SubCategory
void TRExFit::ProduceSystSubCategoryMap(){
    WriteDebugStatus("TRExFit::ProduceSystSubCategoryMap", "filling SubCategory map");

    // special treatment needed for two cases:
    // 1) stat-only fit where all parameters are fixed, see FittingTool::GetGroupedImpact()
    // 2) fit with all Gammas fixed, see FittingTool::GetGroupedImpact()
    fSubCategoryImpactMap.insert(std::make_pair("DUMMY_STATONLY", "FullSyst"));
    fSubCategoryImpactMap.insert(std::make_pair("DUMMY_GAMMAS", "Gammas"));

    // add all systematics, here an "alpha_" prefix is needed
    for(const auto& isyst : fSystematics) {
        if(isyst->fSubCategory=="Gammas" || isyst->fSubCategory=="FullSyst" || isyst->fSubCategory=="combine")
             WriteWarningStatus("TRExFit::ProduceSystSubCategoryMap"," use of \"Gammas\", \"FullSyst\" or \"combine\" as SubCategory names is not supported, you will likely run into issues");
        if(isyst->fType!=Systematic::SHAPE){
            fSubCategoryImpactMap.insert(std::make_pair(("alpha_" + isyst->fNuisanceParameter).c_str(), isyst->fSubCategory));
        }
        else{
            // treat SHAPE systematics separately, since they are not prefixed with "alpha_", but "gamma_shape_" instead
            // need one per bin per region
            for(const auto& reg : fRegions){
                if(reg->fNSamples < 1){
                    WriteErrorStatus("TRExFit::ProduceSystSubCategoryMap", "Can not determine binning (no samples assigned to region?), exiting");
                    exit(EXIT_FAILURE);
                }

                // determine the amount of bins in this region, requires that ReadHistos() was used
                int nRegionBins = reg->fSampleHists[0]->fHist->GetNbinsX();

                // sanity check, this should match for ntuples, for histogram configs fNbins is initialized as 0
                if((reg->fNbins>0) && (reg->fNbins != nRegionBins)){
                    WriteErrorStatus("TRExFit::ProduceSystSubCategoryMap","Inconsistent number of bins determined for region " + reg->fName);
                }

                for(int i_bin=1; i_bin < nRegionBins+1; i_bin++){
                    fSubCategoryImpactMap.insert(std::make_pair(Form("gamma_shape_%s_%s_bin_%d",(isyst->fNuisanceParameter).c_str(),(reg->fName).c_str(),i_bin-1), isyst->fSubCategory));
                }
            }
        }
    }

    // also add norm factors, no "alpha_" needed
    for(const auto& inorm : fNormFactors) {
        if(inorm->fSubCategory=="Gammas" || inorm->fSubCategory=="FullSyst" || inorm->fSubCategory=="combine")
             WriteWarningStatus("TRExFit::ProduceSystSubCategoryMap"," use of \"Gammas\", \"FullSyst\" or \"combine\" as SubCategory names is not supported, you will likely run into issues");
        if (inorm->fName != fPOI) {
            fSubCategoryImpactMap.insert(std::make_pair(inorm->fNuisanceParameter, inorm->fSubCategory));
        }
    }
}

//____________________________________________________________________________________
// combine individual results from grouped impact evaluation into one table
void TRExFit::BuildGroupedImpactTable() const{
    WriteInfoStatus("TRExFit::BuildGroupedImpactTable", "merging grouped impact evaluations");
    const std::string targetName = fName+"/Fits/GroupedImpact"+fSuffix+".txt";

    if(std::ifstream(targetName).good()){
        WriteWarningStatus("TRExFit::BuildGroupedImpactTable","file " + targetName + " already exists, will not overwrite");
    }
    else{
        const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fName+"/Fits/","GroupedImpact"+fSuffix+"_");
        Common::MergeTxTFiles(inPaths, targetName);
    }
}

//____________________________________________________________________________________
//
void TRExFit::RunToys(){
        gSystem->mkdir( (fName+"/Toys").c_str());
        // temporary switch off minos (not needed)
        std::vector<std::string> varMinosTmp = fVarNameMinos;
        fVarNameMinos.clear();
        WriteInfoStatus("TRExFit::RunToys","");
        WriteInfoStatus("TRExFit::RunToys","-------------------------------------------");
        WriteInfoStatus("TRExFit::RunToys","Generating and fitting toys...");
        WriteInfoStatus("TRExFit::RunToys","-------------------------------------------");
        // set all regions as ASIMOVDATA and create combined ws
        std::vector < std:: string > regionsToFit;
        std::map < std::string, int > regionDataType;
        for(const auto& reg : fRegions) {
            regionDataType[reg->fName] = Region::ASIMOVDATA;
        }
        for(const auto& ireg : fRegions) {
            if (fFitRegion == CRONLY && ireg->fRegionType == Region::CONTROL) {
                regionsToFit.emplace_back(ireg->fName);
            } else if (fFitRegion == CRSR && (ireg->fRegionType == Region::CONTROL || ireg->fRegionType == Region::SIGNAL)) {
                regionsToFit.emplace_back(ireg->fName);
            }
        }
        std::unique_ptr<RooWorkspace> ws = PerformWorkspaceCombination(regionsToFit);
        if (!ws){
            WriteErrorStatus("TRExFit::RunToys","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        if (fBinnedLikelihood) {
            FitUtils::SetBinnedLikelihoodOptimisation(ws.get());
        }
        // create map to store fit results
        // create histogram to store fitted POI values
        auto&& POInf = fNormFactors[Common::FindInStringVector(fNormFactorNames,fPOI)];
        double min = POInf->fMin;
        double max = POInf->fMax;
        if (fToysHistoMin < 9900 && fToysHistoMax > -9000){
            min = fToysHistoMin;
            max = fToysHistoMax;
        }
        // get RooStats stuff
        RooStats::ModelConfig mc = *(static_cast<RooStats::ModelConfig*>(ws -> obj("ModelConfig")));
        RooSimultaneous simPdf = *(static_cast<RooSimultaneous*>((mc.GetPdf())));
        RooAbsPdf *pdf = mc.GetPdf();
        RooArgSet obsSet = *(mc.GetObservables());
        std::vector<RooRealVar*> nfs;
        for (const auto& nf: fNormFactors) {
            const std::string& iname = nf->fName;
            if(nf->fExpression.first!="") continue;
            nfs.emplace_back(static_cast<RooRealVar*>((& ws->allVars()[iname.c_str()])));
        }

        //Get the desired NP
        RooRealVar* NPtoShift = nullptr;
        if (fToysPseudodataNP != "") {
            NPtoShift = static_cast<RooRealVar*>((& ws->allVars()[fToysPseudodataNP.c_str()]));
        }

        //For the loop over NPs
        std::string varname{};

        //Create NLL only once
        RooDataSet* dummy = pdf->generate(obsSet, RooFit::Extended());

        // apply external constraints
        RooArgList l;
        std::vector<double> nomVec;
        std::vector<double> tauVec;
        for(const auto& nf : fNormFactors){
            if(nf->fTau!=0){
                if(ws->var(nf->fName.c_str())!=nullptr) l.add(*ws->var(nf->fName.c_str()));
                else if(ws->function(nf->fName.c_str())!=nullptr) l.add(*ws->function(nf->fName.c_str()));
                else WriteWarningStatus("TRExFit::RunToys","Cannot apply tau to norm factor " + nf->fName);
                nomVec.push_back( nf->fNominal );
                tauVec.push_back( nf->fTau );
            }
        }

        const RooArgSet* externalConstraints(nullptr);

        if(!tauVec.empty()) {
            TVectorD nominal(nomVec.size());
            TMatrixDSym cov(tauVec.size());
            for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
                nominal(i_tau) = nomVec[i_tau];
                cov(i_tau,i_tau) = (1./tauVec[i_tau]) * (1./tauVec[i_tau]);
            }
            RooMultiVarGaussian r("regularization","regularization",l,nominal,cov);
            ws->import(r);
            ws->defineSet("myConstraints","regularization");
            simPdf.setStringAttribute("externalConstraints","myConstraints");

            if(simPdf.getStringAttribute("externalConstraints")){
                WriteInfoStatus("TRExFit::RunToys",Form("Building NLL with external constraints %s",simPdf.getStringAttribute("externalConstraints")));
                externalConstraints = ws->set(simPdf.getStringAttribute("externalConstraints"));
            }
        }


        const RooArgSet* glbObs = mc.GetGlobalObservables();
        std::unique_ptr<RooAbsReal> nll(simPdf.createNLL(*dummy,
                                                         RooFit::Constrain(*mc.GetNuisanceParameters()),
                                                         RooFit::GlobalObservables(*glbObs),
                                                         RooFit::Offset(1),
                                                         RooFit::NumCPU(fCPU, RooFit::Hybrid),
                                                         RooFit::Optimize(kTRUE),
                                                         RooFit::ExternalConstraints(*externalConstraints)));

        auto isUnfolding = [&](std::size_t index) {
            if (fNormFactors.at(index)->fName.find("Bin_") != std::string::npos) {
                return true;
            }
            return false;
        };

        std::vector<TH1D> h_toys;
        std::vector<TH1D> h_pulls;
        for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
            h_toys.emplace_back(("h_toys_nf_"+std::to_string(inf)).c_str(),("h_toys_nf_"+std::to_string(inf)).c_str(),fToysHistoNbins,min,max);
            if (fFitType == TRExFit::FitType::UNFOLDING && isUnfolding(inf)) {
                h_pulls.emplace_back(("h_pulls_nf_"+std::to_string(inf)).c_str(),("h_pulls_nf_"+std::to_string(inf)).c_str(),fToysHistoNbins,-3,3);
            }
        }
        for(int i_toy = 0; i_toy < fFitToys; ++i_toy) {

            if (fToysPseudodataNP != "") {
                std::unique_ptr<TIterator> it(mc.GetNuisanceParameters()->createIterator());
                RooRealVar* var = nullptr;
                while( (var = static_cast<RooRealVar*>(it->Next()))){
                    varname = var->GetName();
                    if (varname.find("alpha_")!=std::string::npos) {
                        var->setConstant(1);
                        var->setVal(0);
                    } else {
                        var->setConstant(1);
                        var->setVal(1);
                    }
                }
                delete var;
            }

            // setting POI to constant, not to allow it to fluctuate in toy creation
            for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
                nfs.at(inf)->setConstant(1);
                if (fNormFactorNames.at(inf) == fPOI) {
                    nfs.at(inf)->setVal(fFitPOIAsimov);
                } else {
                    auto&& nf = fNormFactors[inf];
                    nfs.at(inf)->setVal(nf->fNominal);
                }
            }
            if (fToysPseudodataNP != "") {
                NPtoShift->setConstant(1);
                NPtoShift->setVal(fToysPseudodataNPShift);
            }

            WriteInfoStatus("TRExFit::RunToys","Generating toy n. " + std::to_string(i_toy+1) + " out of " + std::to_string(fFitToys) + " toys");
            RooDataSet* toyData = pdf->generate(obsSet, RooFit::Extended());
            // re-set POI to free-floating, and to nominal value
            if (fToysPseudodataNP != "") {
                const RooArgSet* nuis = static_cast<const RooArgSet*>(mc.GetNuisanceParameters());
                if (nuis){
                    RooRealVar* vartmp = nullptr;
                    std::unique_ptr<TIterator> it2(nuis->createIterator());
                    while( (vartmp = static_cast<RooRealVar*>(it2->Next()))){
                        const std::string& np = vartmp->GetName();
                        if (np.find("alpha_")!=std::string::npos) {
                            vartmp->setConstant(0);
                            vartmp->setVal(0);
                        }
                        else if(np.find("gamma_")!=std::string::npos){
                            vartmp->setVal(1);
                            vartmp->setConstant(0);
                        }
                        else {  // for norm factors
                            vartmp->setVal( 1 );
                        }
                    }
                    delete vartmp;
                }
            }
            for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
                if (fFitType==FitType::BONLY && fNormFactorNames.at(inf) == fPOI) {
                    nfs.at(inf)->setVal(0);
                    nfs.at(inf)->setConstant(true);
                } else {
                    auto&& nf = fNormFactors[inf];
                    nfs.at(inf)->setConstant(false);
                    nfs.at(inf)->setVal(nf->fNominal);
                }
            }

            // NP is fixed constant for each fit, and to nominal value
            if (fToysPseudodataNP != "") {
                NPtoShift->setConstant(1);
                NPtoShift->setVal(0);
            }
            // extract POI from fit result and fill histogram
            WriteInfoStatus("TRExFit::RunToys","Fitting toy n. " + std::to_string(i_toy+1));
            //Set new dataset for NLL
            nll->setData(*toyData);
            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
            ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
            ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
            const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this
            RooMinimizer m(*nll); // get MINUIT interface of fit
            m.optimizeConst(2);
            m.setPrintLevel(-1);
            m.setStrategy(1);
            m.setEps(tol);
            m.migrad();
            RooFitResult* r = m.save(); // save fit result

            std::size_t unfIndex(0);
            for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
                h_toys.at(inf).Fill(nfs.at(inf)->getVal());
//                 WriteInfoStatus("TRExFit::RunToys","Toy n. " + std::to_string(i_toy+1) + ", fitted value of NF: " + fNormFactorNames.at(inf) + ": " + std::to_string(nfs.at(inf)->getVal()));
                WriteInfoStatus("TRExFit::RunToys","Toy n. " + std::to_string(i_toy+1) + ", fitted value of NF: " + nfs.at(inf)->GetName() + ": " + std::to_string(nfs.at(inf)->getVal()) + " +/- " + std::to_string(nfs.at(inf)->getError()));

                if (fFitType == TRExFit::FitType::UNFOLDING && isUnfolding(inf)) {
                    const double value = nfs.at(inf)->getError() > 1e-6 ? (nfs.at(inf)->getVal() - 1.) / nfs.at(inf)->getError() : -9999;
                    h_pulls.at(unfIndex).Fill(value);
                    ++unfIndex;
                }
            }

            delete r;
            delete toyData;
        }

        std::unique_ptr<TFile> out (TFile::Open((fName+"/Toys/Toys"+fSuffix+".root").c_str(), "RECREATE"));
        for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
            out->cd();
            TCanvas c("c","c",600,600);
            h_toys.at(inf).Draw("E");
            TF1 g("g","gaus",POInf->fMin,POInf->fMax);
            g.SetLineColor(kRed);
            h_toys.at(inf).Fit("g","RQ");
            g.Draw("same");
            h_toys.at(inf).GetXaxis()->SetTitle(nfs.at(inf)->GetName());
            h_toys.at(inf).GetYaxis()->SetTitle("Pseudo-experiements");
            myText(0.60,0.90,1,Form("Mean  = %.2f #pm %.2f",g.GetParameter(1), g.GetParError(1)));
            myText(0.60,0.85,1,Form("Sigma = %.2f #pm %.2f",g.GetParameter(2),g.GetParError(2)));
            myText(0.60,0.80,1,Form("#chi^{2}/ndf = %.2f / %d",g.GetChisquare(),g.GetNDF()));
            for(const auto& format : TRExFitter::IMAGEFORMAT) {
                c.SaveAs((fName+"/Toys/ToysPlot_"+nfs.at(inf)->GetName()+"."+format).c_str());
            }
            fVarNameMinos = varMinosTmp; // retore Minos settings

            out->cd();
            // Also create a ROOT file
            h_toys.at(inf).Write();
        }

        if (fFitType == TRExFit::FitType::UNFOLDING) {
            DrawToyPullPlot(h_pulls, out.get());
        }

        out->Close();
}

//__________________________________________________________________________________
// Computes the variable string to be used when reading ntuples, for a given region, sample combination
std::string TRExFit::Variable(Region *reg,Sample *smp){
    // protection against nullptr
    if(reg==nullptr){
        WriteErrorStatus("TRExFit::Variable","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        WriteErrorStatus("TRExFit::Variable","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    //
    std::string variable = "";
    // from Region
    if(reg->UseAlternativeVariable(smp->fName)){
        variable = reg->GetAlternativeVariable(smp->fName);
    }
    else{
        variable = reg->fVariable;
    }
    // check the final expression
    if(!Common::CheckExpression(variable)){
        WriteErrorStatus("TRExFit::Variable","Variable expression not valid. Please check: "+variable);
        exit(EXIT_FAILURE);
    }
    //
    return variable;
}

//__________________________________________________________________________________
// Computes the full selection string to be used when reading ntuples, for a given region, sample combination
std::string TRExFit::FullSelection(Region *reg,Sample *smp){
    // protection against nullptr
    if(reg==nullptr){
        WriteErrorStatus("TRExFit::FullSelection","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        WriteErrorStatus("TRExFit::FullSelection","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    //
    std::string selection = "";
    // from Job
    if(fSelection!="" && fSelection!="1"){
        selection = "("+fSelection+")";
    }
    // from Region
    if(reg->UseAlternativeSelection(smp->fName)){
        if(selection!="") selection += " && ";
        selection += "("+reg->GetAlternativeSelection(smp->fName)+")";
    }
    else if(reg->fSelection!="" && reg->fSelection!="1"){
        if(selection!="") selection += " && ";
        selection += "("+reg->fSelection+")";
    }
    // eventually apply IgnoreSelection Sample-option
    if(smp->fIgnoreSelection=="TRUE"){
        selection = "";
    }
    else if(smp->fIgnoreSelection!="FALSE" && smp->fIgnoreSelection!=""){
        selection = Common::ReplaceString(selection,smp->fIgnoreSelection,"1");
    }
    // from Sample
    if(smp->fSelection!="" && smp->fSelection!="1"){
        if(selection!="") selection += " && ";
        selection += "("+smp->fSelection+")";
    }
    // check the final expression
    if(!Common::CheckExpression(selection)){
        WriteErrorStatus("TRExFit::FullSelection","Full selection expression not valid. Please check: "+selection);
        exit(EXIT_FAILURE);
    }
    //
    if(selection=="") selection = "1";
    return selection;
}

//__________________________________________________________________________________
// Computes the full weight string to be used when reading ntuples, for a given region, sample and systematic combination
std::string TRExFit::FullWeight(Region *reg,Sample *smp,Systematic *syst,bool isUp){
    // protection against nullptr
    if(reg==nullptr){
        WriteErrorStatus("TRExFit::FullWeight","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        WriteErrorStatus("TRExFit::FullWeight","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    // if it's Data, just return "1"
    if(smp->fType==Sample::DATA) return "1";
    //
    std::string weight = "";
    // from Job (only for fNormalizedByTheory samples)
    if(fMCweight!="" && fMCweight!="1" && smp->fNormalizedByTheory){
        weight += "("+fMCweight+")";
    }
    // from Region
    if(reg->fMCweight!="" && reg->fMCweight!="1"){
        if(weight!="") weight += " * ";
        weight += "("+reg->fMCweight+")";
    }
    // eventually apply IgnoreWeight Sample-option
    if(smp->fIgnoreWeight=="TRUE"){
        weight = "";
    }
    else if(smp->fIgnoreWeight!="FALSE" && smp->fIgnoreWeight!=""){
        weight = Common::ReplaceString(weight,smp->fIgnoreWeight,"1");
    }
    // from Sample (nominal...
    std::string sampleWeight = "";
    if(smp->fMCweight!="" && smp->fMCweight!="1"){
        sampleWeight = "("+smp->fMCweight+")";
    }
    // ... and systematics)
    if(syst!=nullptr){
        if(syst->fIgnoreWeight!=""){
            weight = Common::ReplaceString(weight,syst->fIgnoreWeight,"1");
            sampleWeight = Common::ReplaceString(sampleWeight,syst->fIgnoreWeight,"1");
        }
        if(isUp){
            if(syst->fWeightUp!=""){
                sampleWeight = syst->fWeightUp;
            }
            else if(syst->fWeightSufUp!=""){
                if(sampleWeight!="") sampleWeight += " * ";
                sampleWeight += "("+syst->fWeightSufUp+")";
            }
        }
        else{
            if(syst->fWeightDown!=""){
                sampleWeight = syst->fWeightDown;
            }
            else if(syst->fWeightSufDown!=""){
                if(sampleWeight!="") sampleWeight += " * ";
                sampleWeight += "("+syst->fWeightSufDown+")";
            }
        }
    }
    if(sampleWeight!=""){
        if(weight!="") weight += " * ";
        weight += "("+sampleWeight+")";
    }
    // add Bootstrap weights
    // fBootstrapSyst means only systematic variation is bootstrapped, not nominal
    // WiP fBootstrapSample means bootstrap on sample and all the correlated systs (so no fntupleFiles -> not sure it captures everything)
    if(fBootstrap!="" && fBootstrapIdx>=0 
      && !(fBootstrapSyst!="" && syst==nullptr)
      && (fBootstrapSample=="" || ( smp->fName==fBootstrapSample && ( !syst || syst->fIsCorrelated ) ) ) ){
        if(weight!="") weight += " * ";
        weight += "("+Common::ReplaceString(fBootstrap,"BootstrapIdx",Form("%d",fBootstrapIdx))+")";
        gRandom->SetSeed(fBootstrapIdx);
    }
    // check the final expression
    WriteDebugStatus("TRExFit::FullWeight","Full weight expression : "+weight);
    if(!Common::CheckExpression(weight)){
        WriteErrorStatus("TRExFit::FullWeight","Full weight expression not valid. Please check: "+weight);
        exit(EXIT_FAILURE);
    }
    //
    if(weight=="") weight = "1";
    return weight;
}

//__________________________________________________________________________________
// Computes the full list of path + file-name + ntuple-name string to be used when reading ntuples, for a given region, sample and systematic combination
std::vector<std::string> TRExFit::FullNtuplePaths(Region *reg,Sample *smp,Systematic *syst,bool isUp){
    // protection against nullptr
    if(reg==nullptr){
        WriteErrorStatus("TRExFit::FullPaths","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        WriteErrorStatus("TRExFit::FullPaths","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    bool isData = (smp->fType==Sample::DATA);
    if(syst!=nullptr){
        if(isUp){
            if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
                if(syst->fNtuplePathsUp.size()  >0) paths = syst->fNtuplePathsUpRefSample;
                if(syst->fNtupleFilesUp.size()  >0) files = syst->fNtupleFilesUpRefSample;
                if(syst->fNtupleNamesUp.size()  >0) names = syst->fNtupleNamesUpRefSample;
            }
            else{
                if(syst->fNtuplePathsUp.size()  >0) paths = syst->fNtuplePathsUp;
                if(syst->fNtupleFilesUp.size()  >0) files = syst->fNtupleFilesUp;
                if(syst->fNtupleNamesUp.size()  >0) names = syst->fNtupleNamesUp;
            }
        }
        else{
            if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
                if(syst->fNtuplePathsDown.size()>0) paths = syst->fNtuplePathsDownRefSample;
                if(syst->fNtupleFilesDown.size()>0) files = syst->fNtupleFilesDownRefSample;
                if(syst->fNtupleNamesDown.size()>0) names = syst->fNtupleNamesDownRefSample;
            }
            else{
                if(syst->fNtuplePathsDown.size()>0) paths = syst->fNtuplePathsDown;
                if(syst->fNtupleFilesDown.size()>0) files = syst->fNtupleFilesDown;
                if(syst->fNtupleNamesDown.size()>0) names = syst->fNtupleNamesDown;
            }
        }
    }
    if(paths.size()==0 && smp->fNtuplePaths.size()>0) paths = smp->fNtuplePaths;
    if(files.size()==0 && smp->fNtupleFiles.size()>0) files = smp->fNtupleFiles;
    if(names.size()==0 && smp->fNtupleNames.size()>0) names = smp->fNtupleNames;
    //
    if(paths.size()==0 && reg->fNtuplePaths.size()>0) paths = reg->fNtuplePaths;
    if(files.size()==0 && reg->fNtupleFiles.size()>0) files = reg->fNtupleFiles;
    if(names.size()==0 && reg->fNtupleNames.size()>0) names = reg->fNtupleNames;
    //
    if(paths.size()==0 && fNtuplePaths.size()>0) paths = fNtuplePaths;
    if(files.size()==0 && fNtupleFiles.size()>0) files = fNtupleFiles;
    if(names.size()==0 && fNtupleNames.size()>0) names = fNtupleNames;
    //
    // now combining suffs, with all the combinations instead of giving priority
    // (same order as above used for suffix order - OK? FIXME)
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), Common::ToVec(syst->fNtuplePathSufUpRefSample) );
            else     pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), Common::ToVec(syst->fNtuplePathSufDownRefSample) );
        }
        else{
            if(isUp) pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), Common::ToVec(syst->fNtuplePathSufUp) );
            else     pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), Common::ToVec(syst->fNtuplePathSufDown) );
        }
    }
    else{
        pathSuffs = Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufUpRefSample) );
            else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufDownRefSample) );
        }
        else{
            if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufUp) );
            else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufDown) );
        }
    }
    else{
        fileSuffs = Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufUpRefSample) );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufDownRefSample) );
        }
        else{
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufUp) );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufDown) );
        }
    }
    else{
        nameSuffs = Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs );
    }
    //
    // And finally put everything together
    fullPaths = Common::CreatePathsList( paths,pathSuffs, files,fileSuffs, names,nameSuffs );
    return fullPaths;
}

//__________________________________________________________________________________
// Computes the full list of path + file-name + histogram-name string to be used when reading ntuples, for a given region, sample and systematic combination
std::vector<std::string> TRExFit::FullHistogramPaths(Region *reg,Sample *smp,Systematic *syst,bool isUp){
    // protection against nullptr
    if(reg==nullptr){
        WriteErrorStatus("TRExFit::FullHistogramPaths","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        WriteErrorStatus("TRExFit::FullHistogramPaths","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    bool isData = (smp->fType==Sample::DATA);
    if(syst!=nullptr){
        if(isUp){
            if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
                if(syst->fHistoPathsUp.size()  >0) paths = syst->fHistoPathsUpRefSample;
                if(syst->fHistoFilesUp.size()  >0) files = syst->fHistoFilesUpRefSample;
                if(syst->fHistoNamesUp.size()  >0) names = syst->fHistoNamesUpRefSample;
            }
            else{
                if(syst->fHistoPathsUp.size()  >0) paths = syst->fHistoPathsUp;
                if(syst->fHistoFilesUp.size()  >0) files = syst->fHistoFilesUp;
                if(syst->fHistoNamesUp.size()  >0) names = syst->fHistoNamesUp;
            }
        }
        else{
            if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
                if(syst->fHistoPathsDown.size()>0) paths = syst->fHistoPathsDownRefSample;
                if(syst->fHistoFilesDown.size()>0) files = syst->fHistoFilesDownRefSample;
                if(syst->fHistoNamesDown.size()>0) names = syst->fHistoNamesDownRefSample;
            }
            else{
                if(syst->fHistoPathsDown.size()>0) paths = syst->fHistoPathsDown;
                if(syst->fHistoFilesDown.size()>0) files = syst->fHistoFilesDown;
                if(syst->fHistoNamesDown.size()>0) names = syst->fHistoNamesDown;
            }
        }
    }
    if(paths.size()==0 && smp->fHistoPaths.size()>0) paths = smp->fHistoPaths;
    if(files.size()==0 && smp->fHistoFiles.size()>0) files = smp->fHistoFiles;
    if(names.size()==0 && smp->fHistoNames.size()>0) names = smp->fHistoNames;
    //
    if(paths.size()==0 && reg->fHistoPaths.size()>0) paths = reg->fHistoPaths;
    if(files.size()==0 && reg->fHistoFiles.size()>0) files = reg->fHistoFiles;
    if(names.size()==0 && reg->fHistoNames.size()>0) names = reg->fHistoNames;
    //
    if(paths.size()==0 && fHistoPaths.size()>0) paths = fHistoPaths;
    if(files.size()==0 && fHistoFiles.size()>0) files = fHistoFiles;
    if(names.size()==0 && fHistoNames.size()>0) names = fHistoNames;
    //
    // now combining suffs, with all the combinations instead of giving priority
    // (same order as above used for suffix order - OK? FIXME)
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), Common::ToVec(syst->fHistoPathSufUpRefSample) );
            else     pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), Common::ToVec(syst->fHistoPathSufDownRefSample) );
        }
        else{
            if(isUp) pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), Common::ToVec(syst->fHistoPathSufUp) );
            else     pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), Common::ToVec(syst->fHistoPathSufDown) );
        }
    }
    else{
        pathSuffs = Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), Common::ToVec(syst->fHistoFileSufUpRefSample) );
            else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), Common::ToVec(syst->fHistoFileSufDownRefSample) );
        }
        else{
            if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), Common::ToVec(syst->fHistoFileSufUp) );
            else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), Common::ToVec(syst->fHistoFileSufDown) );
        }
    }
    else{
        fileSuffs = Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), Common::ToVec(syst->fHistoNameSufUpRefSample) );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), Common::ToVec(syst->fHistoNameSufDownRefSample) );
        }
        else{
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), Common::ToVec(syst->fHistoNameSufUp) );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), Common::ToVec(syst->fHistoNameSufDown) );
        }
    }
    else{
      nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs), fHistoNamesNominal );
    }
    //
    // And finally put everything together
    fullPaths = Common::CreatePathsList( paths,pathSuffs, files,fileSuffs, names,nameSuffs );
    return fullPaths;
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> TRExFit::GetSampleHistFromName(const Region* const reg, const std::string& name) const{
    for (int i_smp = 0; i_smp < reg->fNSamples; ++i_smp){
       if (reg->fSampleHists[i_smp]->fName == name){
            return reg->fSampleHists[i_smp];
        }
    }
    return nullptr;
}

//__________________________________________________________________________________
//
TH1* TRExFit::CopySmoothedHisto(const SampleHist* const sh, const TH1* const nominal, const TH1* const up, const TH1* const down, const bool isUp) const{
    TH1* currentNominal = sh->fHist.get();

    if (currentNominal->GetNbinsX() != nominal->GetNbinsX()){
        WriteErrorStatus("TRExFit::CopySmoothedHisto", "Histograms to smooth have different binning (nominals)");
        exit(EXIT_FAILURE);
    }
    if (currentNominal->GetNbinsX() != up->GetNbinsX()){
        WriteErrorStatus("TRExFit::CopySmoothedHisto", "Histograms to smooth have different binning (nominal vs up)");
        exit(EXIT_FAILURE);
    }
    if (currentNominal->GetNbinsX() != down->GetNbinsX()){
        WriteErrorStatus("TRExFit::CopySmoothedHisto", "Histograms to smooth have different binning (nominal vs down)");
        exit(EXIT_FAILURE);
    }

    TH1* result = static_cast<TH1*>(currentNominal->Clone());

    for (int ibin = 1; ibin <= currentNominal->GetNbinsX(); ++ibin){
        double ratio = 1;
        if (nominal->GetBinContent(ibin) != 0){
            if (isUp){
                ratio = up->GetBinContent(ibin)/nominal->GetBinContent(ibin);
            } else {
                ratio = down->GetBinContent(ibin)/nominal->GetBinContent(ibin);
            }
        }
        ratio*= currentNominal->GetBinContent(ibin);

        result->SetBinContent(ibin, ratio);
    }

    return result;
}

//__________________________________________________________________________________
//
int TRExFit::GetSystIndex(const SampleHist* const sh, const std::string& name) const{
    for (std::size_t i = 0; i < sh->fSyst.size(); ++i){
        if (sh->fSyst[i]->fName == name){
            return i;
        }
    }

    return -1;
}

//__________________________________________________________________________________
//
std::shared_ptr<SystematicHist> TRExFit::CombineSpecialHistos(std::shared_ptr<SystematicHist> orig,
                                                              const std::vector<std::shared_ptr<SystematicHist> >& vec,
                                                              Systematic::COMBINATIONTYPE type,
                                                              const SampleHist* sh) const {
    if (vec.size() == 0) return nullptr;
    if (!orig) return nullptr;
    if (!orig->fHistUp) return nullptr;
    if (!orig->fHistDown) return nullptr;
    if (!sh) return nullptr;
    if (!sh->fHist) return nullptr;

    const int nbins = sh->fHist->GetNbinsX();
    if (type == Systematic::COMBINATIONTYPE::ENVELOPE) {
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            std::vector<double> hist_max;
            std::vector<double> hist_min;
            for (const auto& isyst : vec) {
                if (!isyst) continue;
                if (!isyst->fHistUp) continue;
                if (!isyst->fHistDown) continue;
                const double diff_up = isyst->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                const double diff_down = isyst->fHistDown->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                if (diff_up > 0) {
                    hist_max.emplace_back(diff_up);
                } else {
                    hist_min.emplace_back(diff_up);
                }
                if (diff_down > 0) {
                    hist_max.emplace_back(diff_down);
                } else {
                    hist_min.emplace_back(diff_down);
                }
            }

            double max(0);
            double min(0);

            if (hist_max.size() > 0) {
                max = *std::max_element(hist_max.begin(), hist_max.end());
            }
            if (hist_min.size() > 0) {
                min = *std::min_element(hist_min.begin(), hist_min.end());
            }

            if (max >= 0) {
                orig->fHistUp->SetBinContent(ibin, max + sh->fHist->GetBinContent(ibin));
            }
            if (min <= 0) {
                orig->fHistDown->SetBinContent(ibin, min + sh->fHist->GetBinContent(ibin));
            }
        }
    } else if (type == Systematic::COMBINATIONTYPE::STANDARDDEVIATION) {
        for (int ibin = 1; ibin <= nbins; ++ ibin) {
            std::vector<double> content;
            for (const auto& isyst : vec) {
                if (!isyst) continue;
                if (!isyst->fHistUp) continue;
                if (!isyst->fHistDown) continue;
                const double up = isyst->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                const double down = isyst->fHistDown->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                // take the larger (signed) variation
                const double max = (std::fabs(up) > std::fabs(down) ? up : down);
                content.emplace_back(max);
            }

            if (content.size() < 2) continue;

            // calculate the standard deviation in each bin
            double sum(0);
            double sigma(0);

            for (const auto& value : content) {
                sum+= value;
            }

            // get mean
            const double mean = static_cast<double>(sum / content.size());

            // get std
            for (const auto& value : content) {
                sigma+= (value-mean)*(value-mean);
            }

            sigma = static_cast<double>(sigma / (content.size() - 1));
            sigma = std::sqrt(sigma);

            // now set the bin content
            orig->fHistUp->SetBinContent(ibin,    sigma + sh->fHist->GetBinContent(ibin));
            orig->fHistDown->SetBinContent(ibin, -sigma + sh->fHist->GetBinContent(ibin));
        }
    } else {
        WriteErrorStatus("TRExFit::CombineSpecialHistos","Unknown combination type, this should not happen!");
        exit(EXIT_FAILURE);
    }

    return orig;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::GetUniqueSystNamesWithoutGamma() const {
    std::vector<std::string> result;
    for(const auto& isyst : fSystematics) {
        if (isyst->fType == Systematic::SHAPE) continue;
        if (std::find(result.begin(), result.end(), isyst->fName) == result.end()) {
            result.emplace_back(isyst->fName);
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<Region*> TRExFit::GetNonValidationRegions() const {
    std::vector<Region*> result;
    for (const auto& ireg : fRegions) {
        if (ireg->fRegionType == Region::VALIDATION) continue;
        result.emplace_back(ireg);
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<Sample> > TRExFit::GetNonDataNonGhostSamples() const {
    std::vector<std::shared_ptr<Sample> > result;

    for (const auto& isample : fSamples) {
        if (isample->fType == Sample::DATA) continue;
        if (isample->fType == Sample::GHOST) continue;

        result.emplace_back(isample);
    }

    return result;
}

//__________________________________________________________________________________
//
void TRExFit::DropBins() {

    for(const auto& reg : fRegions){

        const std::vector<int>& blindedBins = reg->GetAutomaticDropBins() ?
                Common::GetBlindedBins(reg,fBlindingType,fBlindingThreshold) : reg->fDropBins;
        if (blindedBins.size() == 0) continue;
        for(const auto& smp : fSamples){
            // eventually skip sample / region combination
            if(Common::FindInStringVector(smp->fRegions,reg->fName)<0) continue;
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(!sh) continue;
            Common::DropBins(sh->fHist.get(), blindedBins);
            for(auto& syst : smp->fSystematics) {
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic( syst->fName );
                if(!syh) continue;
                Common::DropBins(syh->fHistUp.get(),   blindedBins);
                Common::DropBins(syh->fHistDown.get(), blindedBins);
                if(syh->fHistShapeUp)   Common::DropBins(syh->fHistShapeUp.get(),   blindedBins);
                if(syh->fHistShapeDown) Common::DropBins(syh->fHistShapeDown.get(), blindedBins);
            }
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::PrepareUnfolding() {
    // Prepare the folder structure
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/UnfoldingHistograms").c_str());

    // Open the otuput ROOT file
    std::unique_ptr<TFile> outputFile(TFile::Open((fName+"/UnfoldingHistograms/FoldedHistograms.root").c_str(), "RECREATE"));
    if (!outputFile) {
        WriteErrorStatus("TRExFit::PrepareUnfolding", "Cannot open the output file at: " + fName + "/UnfoldingHistograms/FoldedHistograms.root");
        exit(EXIT_FAILURE);
    }

    FoldingManager manager{};
    manager.SetMatrixOrientation(fMatrixOrientation);

    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);

    std::unique_ptr<TH1> alternativeTruth(nullptr);

    {
        std::unique_ptr<TH1> truth(nullptr);
        for (const auto& itruth : fTruthSamples) {
            if (itruth->GetName() == fNominalTruthSample) {
                truth = itruth->GetHisto(this);
                break;
            }
        }
        if (!truth) {
            exit(EXIT_FAILURE);
        }
        if (truth->GetNbinsX() != fNumberUnfoldingTruthBins) {
            WriteErrorStatus("TRExFit::PrepareUnfolding", "The number of truth bins doesnt match the value from the config");
            exit(EXIT_FAILURE);
        }
        manager.SetTruthDistribution(truth.get());
        manager.WriteTruthToHisto(outputFile.get(), "", "truth_distribution");

        // read the alternative truth sample
        if (fAlternativeAsimovTruthSample != "") {
            for (const auto& itruth : fTruthSamples) {
                if (itruth->GetName() == fAlternativeAsimovTruthSample) {
                    alternativeTruth = itruth->GetHisto(this);
                    break;
                }
            }
        }
    }

    // loop over regions
    for (const auto& ireg : fRegions) {

        // only signal regions are processed at this step
        if (ireg->fRegionType != Region::RegionType::SIGNAL) continue;

        // loop over all samples
        for (const auto& isample : fUnfoldingSamples) {

            // skip samples not associated to the region
            if(isample->fRegions[0] != "all" &&
                Common::FindInStringVector(isample->fRegions, ireg->fName) < 0) continue;

            // first process nominal
            if (isample->GetHasResponse()) {
                const std::vector<std::string>& fullResponsePaths = FullResponseMatrixPaths(ireg, isample.get());

                std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(fullResponsePaths);
                if (!matrix) {
                    WriteErrorStatus("TRExFit::PrepareUnfolding", "Cannot read the response matrix!");
                    exit(EXIT_FAILURE);
                }
                const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
                const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
                if (nRecoBins != ireg->fNumberUnfoldingRecoBins) {
                    WriteErrorStatus("TRExFit::PrepareUnfolding", "Number of reco bins do not match the number of reco bins for the response matrix in region: " + ireg->fName);
                    exit(EXIT_FAILURE);
                }
                if (nTruthBins != fNumberUnfoldingTruthBins) {
                    WriteErrorStatus("TRExFit::PrepareUnfolding", "Number of truth bins do not match the number of truth bins for the response matrix in regoin: " + ireg->fName);
                    exit(EXIT_FAILURE);
                }

                PlotMigrationResponse(matrix.get(), false, ireg->fName, "");

                manager.SetResponseMatrix(matrix.get());
                outputFile->cd();
                matrix->Write((ireg->fName + "_" + isample->GetName() + "_response").c_str());
            } else {
                // need to add acceptance, selection and migration
                {
                    const std::vector<std::string>& fullMigrationMatrixPaths = FullMigrationMatrixPaths(ireg, isample.get());
                    std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(fullMigrationMatrixPaths);
                    if (!matrix) {
                        exit(EXIT_FAILURE);
                    }
                    const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
                    const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
                    if (nRecoBins != ireg->fNumberUnfoldingRecoBins) {
                        WriteErrorStatus("TRExFit::PrepareUnfolding", "Number of reco bins do not match the number of reco bins for the migration matrix in region: " + ireg->fName);
                        exit(EXIT_FAILURE);
                    }
                    if (nTruthBins != fNumberUnfoldingTruthBins) {
                        WriteErrorStatus("TRExFit::PrepareUnfolding", "Number of truth bins do not match the number of truth bins for the migration matrix in region: " + ireg->fName);
                        exit(EXIT_FAILURE);
                    }

                    UnfoldingTools::NormalizeMatrix(matrix.get(), !horizontal);

                    PlotMigrationResponse(matrix.get(), true, ireg->fName, "");

                    // pass the migration to the tool
                    manager.SetMigrationMatrix(matrix.get(), false);
                    outputFile->cd();
                    matrix->Write((ireg->fName + "_" + isample->GetName() + "_migration").c_str());
                }

                // add selection eff
                {
                    const std::vector<std::string>& fullSelectionEffPaths = FullSelectionEffPaths(ireg, isample.get());
                    std::unique_ptr<TH1> eff = Common::CombineHistosFromFullPaths(fullSelectionEffPaths);
                    if (!eff) {
                        exit(EXIT_FAILURE);
                    }
                    const int nbins = eff->GetNbinsX();
                    if (nbins != fNumberUnfoldingTruthBins) {
                        WriteErrorStatus("TRExFit::PrepareUnfolding", "Number of efficiency selection bins doesnt match the number of truth bins");
                        exit(EXIT_FAILURE);
                    }

                    manager.SetSelectionEfficiency(eff.get());
                }

                // add acceptance
                if (fHasAcceptance || isample->GetHasAcceptance() || ireg->fHasAcceptance) {
                    const std::vector<std::string>& fullAcceptancePaths = FullAcceptancePaths(ireg, isample.get());
                    std::unique_ptr<TH1> acc = Common::CombineHistosFromFullPaths(fullAcceptancePaths);
                    if (!acc) {
                        exit(EXIT_FAILURE);
                    }
                    const int nbins = acc->GetNbinsX();
                    if (nbins != ireg->fNumberUnfoldingRecoBins) {
                        WriteErrorStatus("TRExFit::PrepareUnfolding", "Number of acceptance bins doesnt match the number of reco bins in region " + ireg->fName);
                        exit(EXIT_FAILURE);
                    }

                    manager.SetAcceptance(acc.get());
                }

                manager.CalculateResponseMatrix(true);
                PlotMigrationResponse(manager.GetResponseMatrix(), false, ireg->fName, "");
                outputFile->cd();
                manager.GetResponseMatrix()->Write((ireg->fName + "_" + isample->GetName() + "_response").c_str());
            }

            manager.FoldTruth();


            // create the folder structure
            TDirectory* dir = dynamic_cast<TDirectory*>(outputFile->Get("nominal"));
            if (!dir) {
                outputFile->cd();
                outputFile->mkdir("nominal");
            }

            const std::string histoName = ireg->fName + "_" + isample->GetName();
            manager.WriteFoldedToHisto(outputFile.get(), "nominal", histoName);

            // fold and store the altenative truth sample
            if (fAlternativeAsimovTruthSample != "") {
                std::unique_ptr<TH1> tmp = manager.TotalFold(alternativeTruth.get());
                outputFile->cd();
                tmp->Write((ireg->fName + "_AlternativeAsimov").c_str());
            }

            // Process systematics
            for (const auto& isyst : fUnfoldingSystematics) {
                if (!isyst) continue;
                if (isyst->GetName() == "Dummy") continue; // wtf?

                if(isyst->fRegions.at(0) != "all" &&
                     Common::FindInStringVector(isyst->fRegions, ireg->fName) < 0) continue;
                if(isyst->fSamples.at(0) != "all" &&
                     Common::FindInStringVector(isyst->fSamples, isample->GetName()) < 0) continue;

                ProcessUnfoldingSystematics(&manager,
                                            outputFile.get(),
                                            ireg,
                                            isample.get(),
                                            isyst.get());
            }
        }
    }

    outputFile->Close();
}

//__________________________________________________________________________________
//
void TRExFit::ProcessUnfoldingSystematics(FoldingManager* manager,
                                          TFile* file,
                                          const Region* reg,
                                          const UnfoldingSample* sample,
                                          const UnfoldingSystematic* syst) const {

    const bool horizontal = (fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
    // lambda for processing one (up or down) variation
    auto ProcessOneVariation = [&](const bool isUp) {
        if (syst->GetHasResponse()) {
            const std::vector<std::string>& paths = FullResponseMatrixPaths(reg, sample, syst, true);
            std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(paths);
            if (!matrix) {
                exit(EXIT_FAILURE);
            }
            const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
            const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
            if (nRecoBins != reg->fNumberUnfoldingRecoBins) {
                WriteErrorStatus("TRExFit::ProcessUnfoldingSystematics", "Number of reco bins do not match the number of reco bins for the response matrix in region: " + reg->fName);
                exit(EXIT_FAILURE);
            }
            if (nTruthBins != fNumberUnfoldingTruthBins) {
                WriteErrorStatus("TRExFit::ProcessUnfoldingSystematics", "Number of truth bins do not match the number of truth bins for the response matrix: " + reg->fName);
                exit(EXIT_FAILURE);
            }
            if (fPlotSystematicMigrations) {
                PlotMigrationResponse(matrix.get(), false, reg->fName, syst->GetName());
            }

            manager->SetResponseMatrix(matrix.get());
            
            if (fPlotSystematicMigrations) {
                file->cd();
                matrix->Write((reg->fName + "_" + syst->GetName() + "_response").c_str());
            }
        } else {
            {
                /// migration first
                const std::vector<std::string>& paths = FullMigrationMatrixPaths(reg, sample, syst, true);
                std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(paths);
                if (!matrix) {
                    exit(EXIT_FAILURE);
                }

                const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
                const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
                if (nRecoBins != reg->fNumberUnfoldingRecoBins) {
                    WriteErrorStatus("TRExFit::ProcessUnfoldingSystematics", "Number of reco bins do not match the number of reco bins for the migration matrix in region: " + reg->fName);
                    exit(EXIT_FAILURE);
                }
                if (nTruthBins != fNumberUnfoldingTruthBins) {
                    WriteErrorStatus("TRExFit::ProcessUnfoldingSystematics", "Number of truth bins do not match the number of truth bins for the migration matrix: " + reg->fName);
                    exit(EXIT_FAILURE);
                }
                UnfoldingTools::NormalizeMatrix(matrix.get(), !horizontal);

                if (fPlotSystematicMigrations) {
                    PlotMigrationResponse(matrix.get(), true, reg->fName, syst->GetName());
                }

                manager->SetMigrationMatrix(matrix.get(), false);
                if (fPlotSystematicMigrations) {
                    file->cd();
                    matrix->Write((reg->fName + "_" + syst->GetName() + "_migration").c_str());
                }
            }

            // selectio eff now
            {
                const std::vector<std::string>& paths = FullSelectionEffPaths(reg, sample, syst, true);
                std::unique_ptr<TH1> eff = Common::CombineHistosFromFullPaths(paths);
                if (!eff) {
                    exit(EXIT_FAILURE);
                }

                const int nbins = eff->GetNbinsX();
                if (nbins != fNumberUnfoldingTruthBins) {
                    WriteErrorStatus("TRExFit::ProcessUnfoldingSystematics", "Number of efficiency selection bins doesnt match the number of truth bins");
                    exit(EXIT_FAILURE);
                }

                manager->SetSelectionEfficiency(eff.get());
            }

            if (fHasAcceptance || syst->GetHasAcceptance() || reg->fHasAcceptance) {
                const std::vector<std::string>& paths = FullAcceptancePaths(reg, sample, syst, true);
                std::unique_ptr<TH1> acc = Common::CombineHistosFromFullPaths(paths);
                if (!acc) {
                    exit(EXIT_FAILURE);
                }

                const int nbins = acc->GetNbinsX();
                if (nbins != reg->fNumberUnfoldingRecoBins) {
                    WriteErrorStatus("TRExFit::ProcessUnfoldingSystematics", "Number of acceptance bins doesnt match the number of reco bins in region " + reg->fName);
                    exit(EXIT_FAILURE);
                }

                manager->SetAcceptance(acc.get());
            }
            manager->CalculateResponseMatrix(true);
            if (fPlotSystematicMigrations) {
                PlotMigrationResponse(manager->GetResponseMatrix(), false, reg->fName, syst->GetName());
                file->cd();
                manager->GetResponseMatrix()->Write((reg->fName + "_" + syst->GetName() + "_response").c_str());
            }
        }
        manager->FoldTruth();

        // create the folder structure
        const std::string folderName = isUp ? syst->GetName() + "_Up" : syst->GetName() + "_Down";
        TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(folderName.c_str()));
        if (!dir) {
            file->cd();
            file->mkdir(folderName.c_str());
        }

        const std::string histoName = reg->fName + "_" + sample->GetName();
        manager->WriteFoldedToHisto(file, folderName, histoName);
    };

    if (syst->fHasUpVariation) {
        ProcessOneVariation(true);
    }
    if (syst->fHasDownVariation) {
        ProcessOneVariation(false);
    }
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullResponseMatrixPaths(const Region* reg,
                                                          const UnfoldingSample* smp,
                                                          const UnfoldingSystematic* syst,
                                                          const bool isUp) const {
    // protection against nullptr
    if(!reg) {
        WriteErrorStatus("TRExFit::FullResponseMatrixPaths","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        WriteErrorStatus("TRExFit::FullResponseMatrixPaths","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fResponseMatrixPathsUp.size()  >0) paths = syst->fResponseMatrixPathsUp;
            if(syst->fResponseMatrixFilesUp.size()  >0) files = syst->fResponseMatrixFilesUp;
            if(syst->fResponseMatrixNamesUp.size()  >0) names = syst->fResponseMatrixNamesUp;
        } else {
            if(syst->fResponseMatrixPathsDown.size()>0) paths = syst->fResponseMatrixPathsDown;
            if(syst->fResponseMatrixFilesDown.size()>0) files = syst->fResponseMatrixFilesDown;
            if(syst->fResponseMatrixNamesDown.size()>0) names = syst->fResponseMatrixNamesDown;
        }
    }
    if(paths.size()==0 && smp->fResponseMatrixPaths.size()>0) paths = smp->fResponseMatrixPaths;
    if(files.size()==0 && smp->fResponseMatrixFiles.size()>0) files = smp->fResponseMatrixFiles;
    if(names.size()==0 && smp->fResponseMatrixNames.size()>0) names = smp->fResponseMatrixNames;

    if(paths.size()==0 && reg->fResponseMatrixPaths.size()>0) paths = reg->fResponseMatrixPaths;
    if(files.size()==0 && reg->fResponseMatrixFiles.size()>0) files = reg->fResponseMatrixFiles;
    if(names.size()==0 && reg->fResponseMatrixNames.size()>0) names = reg->fResponseMatrixNames;

    if(paths.size()==0 && fResponseMatrixPaths.size()>0) paths = fResponseMatrixPaths;
    if(files.size()==0 && fResponseMatrixFiles.size()>0) files = fResponseMatrixFiles;
    if(names.size()==0 && fResponseMatrixNames.size()>0) names = fResponseMatrixNames;

    if(syst) {
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixPathSuffs,smp->fResponseMatrixPathSuffs), syst->fResponseMatrixPathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixPathSuffs,smp->fResponseMatrixPathSuffs), syst->fResponseMatrixPathSuffsDown);
    }
    else{
        pathSuffs = Common::CombinePathSufs(reg->fResponseMatrixPathSuffs, smp->fResponseMatrixPathSuffs);
    }

    if(syst) {
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixFileSuffs, smp->fResponseMatrixFileSuffs), syst->fResponseMatrixFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixFileSuffs, smp->fResponseMatrixFileSuffs), syst->fResponseMatrixFileSuffsDown);
    }
    else{
        fileSuffs = Common::CombinePathSufs(reg->fResponseMatrixFileSuffs, smp->fResponseMatrixFileSuffs);
    }

    if(syst) {
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixNameSuffs, smp->fResponseMatrixNameSuffs), syst->fResponseMatrixNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixNameSuffs, smp->fResponseMatrixNameSuffs), syst->fResponseMatrixNameSuffsDown);
    }
    else{
      nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixNameSuffs, smp->fResponseMatrixNameSuffs), fResponseMatrixNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullMigrationMatrixPaths(const Region* reg,
                                                           const UnfoldingSample* smp,
                                                           const UnfoldingSystematic* syst,
                                                           const bool isUp) const {
    // protection against nullptr
    if(!reg) {
        WriteErrorStatus("TRExFit::FullMigrationMatrixPaths","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        WriteErrorStatus("TRExFit::FullMigrationMatrixPaths","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fMigrationPathsUp.size()  >0) paths = syst->fMigrationPathsUp;
            if(syst->fMigrationFilesUp.size()  >0) files = syst->fMigrationFilesUp;
            if(syst->fMigrationNamesUp.size()  >0) names = syst->fMigrationNamesUp;
        } else {
            if(syst->fMigrationPathsDown.size()>0) paths = syst->fMigrationPathsDown;
            if(syst->fMigrationFilesDown.size()>0) files = syst->fMigrationFilesDown;
            if(syst->fMigrationNamesDown.size()>0) names = syst->fMigrationNamesDown;
        }
    }
    if(paths.size()==0 && smp->fMigrationPaths.size()>0) paths = smp->fMigrationPaths;
    if(files.size()==0 && smp->fMigrationFiles.size()>0) files = smp->fMigrationFiles;
    if(names.size()==0 && smp->fMigrationNames.size()>0) names = smp->fMigrationNames;

    if(paths.size()==0 && reg->fMigrationPaths.size()>0) paths = reg->fMigrationPaths;
    if(files.size()==0 && reg->fMigrationFiles.size()>0) files = reg->fMigrationFiles;
    if(names.size()==0 && reg->fMigrationNames.size()>0) names = reg->fMigrationNames;

    if(paths.size()==0 && fMigrationPaths.size()>0) paths = fMigrationPaths;
    if(files.size()==0 && fMigrationFiles.size()>0) files = fMigrationFiles;
    if(names.size()==0 && fMigrationNames.size()>0) names = fMigrationNames;

    if(syst) {
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationPathSuffs,smp->fMigrationPathSuffs), syst->fMigrationPathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationPathSuffs,smp->fMigrationPathSuffs), syst->fMigrationPathSuffsDown);
    }
    else{
        pathSuffs = Common::CombinePathSufs(reg->fMigrationPathSuffs, smp->fMigrationPathSuffs);
    }

    if(syst) {
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationFileSuffs, smp->fMigrationFileSuffs), syst->fMigrationFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationFileSuffs, smp->fMigrationFileSuffs), syst->fMigrationFileSuffsDown);
    }
    else{
        fileSuffs = Common::CombinePathSufs(reg->fMigrationFileSuffs, smp->fMigrationFileSuffs);
    }

    if(syst) {
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationNameSuffs, smp->fMigrationNameSuffs), syst->fMigrationNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationNameSuffs, smp->fMigrationNameSuffs), syst->fMigrationNameSuffsDown);
    }
    else{
      nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationNameSuffs, smp->fMigrationNameSuffs), fMigrationNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullAcceptancePaths(const Region* reg,
                                                      const UnfoldingSample* smp,
                                                      const UnfoldingSystematic* syst,
                                                      const bool isUp) const {
    // protection against nullptr
    if(!reg) {
        WriteErrorStatus("TRExFit::FullAcceptancePaths","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        WriteErrorStatus("TRExFit::FullAcceptancePaths","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fAcceptancePathsUp.size()  >0) paths = syst->fAcceptancePathsUp;
            if(syst->fAcceptanceFilesUp.size()  >0) files = syst->fAcceptanceFilesUp;
            if(syst->fAcceptanceNamesUp.size()  >0) names = syst->fAcceptanceNamesUp;
        } else {
            if(syst->fAcceptancePathsDown.size()>0) paths = syst->fAcceptancePathsDown;
            if(syst->fAcceptanceFilesDown.size()>0) files = syst->fAcceptanceFilesDown;
            if(syst->fAcceptanceNamesDown.size()>0) names = syst->fAcceptanceNamesDown;
        }
    }
    if(paths.size()==0 && smp->fAcceptancePaths.size()>0) paths = smp->fAcceptancePaths;
    if(files.size()==0 && smp->fAcceptanceFiles.size()>0) files = smp->fAcceptanceFiles;
    if(names.size()==0 && smp->fAcceptanceNames.size()>0) names = smp->fAcceptanceNames;

    if(paths.size()==0 && reg->fAcceptancePaths.size()>0) paths = reg->fAcceptancePaths;
    if(files.size()==0 && reg->fAcceptanceFiles.size()>0) files = reg->fAcceptanceFiles;
    if(names.size()==0 && reg->fAcceptanceNames.size()>0) names = reg->fAcceptanceNames;

    if(paths.size()==0 && fAcceptancePaths.size()>0) paths = fAcceptancePaths;
    if(files.size()==0 && fAcceptanceFiles.size()>0) files = fAcceptanceFiles;
    if(names.size()==0 && fAcceptanceNames.size()>0) names = fAcceptanceNames;

    if(syst) {
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptancePathSuffs,smp->fAcceptancePathSuffs), syst->fAcceptancePathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptancePathSuffs,smp->fAcceptancePathSuffs), syst->fAcceptancePathSuffsDown);
    }
    else{
        pathSuffs = Common::CombinePathSufs(reg->fAcceptancePathSuffs, smp->fAcceptancePathSuffs);
    }

    if(syst) {
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceFileSuffs, smp->fAcceptanceFileSuffs), syst->fAcceptanceFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceFileSuffs, smp->fAcceptanceFileSuffs), syst->fAcceptanceFileSuffsDown);
    }
    else{
        fileSuffs = Common::CombinePathSufs(reg->fAcceptanceFileSuffs, smp->fAcceptanceFileSuffs);
    }

    if(syst) {
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceNameSuffs, smp->fAcceptanceNameSuffs), syst->fAcceptanceNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceNameSuffs, smp->fAcceptanceNameSuffs), syst->fAcceptanceNameSuffsDown);
    }
    else{
      nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceNameSuffs, smp->fAcceptanceNameSuffs), fAcceptanceNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullSelectionEffPaths(const Region* reg,
                                                        const UnfoldingSample* smp,
                                                        const UnfoldingSystematic* syst,
                                                        const bool isUp) const {
    // protection against nullptr
    if(!reg) {
        WriteErrorStatus("TRExFit::FullSelectionEffPaths","Null pointer for Region.");
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        WriteErrorStatus("TRExFit::FullSelectionEffPaths","Null pointer for Sample.");
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fSelectionEffPathsUp.size()  >0) paths = syst->fSelectionEffPathsUp;
            if(syst->fSelectionEffFilesUp.size()  >0) files = syst->fSelectionEffFilesUp;
            if(syst->fSelectionEffNamesUp.size()  >0) names = syst->fSelectionEffNamesUp;
        } else {
            if(syst->fSelectionEffPathsDown.size()>0) paths = syst->fSelectionEffPathsDown;
            if(syst->fSelectionEffFilesDown.size()>0) files = syst->fSelectionEffFilesDown;
            if(syst->fSelectionEffNamesDown.size()>0) names = syst->fSelectionEffNamesDown;
        }
    }
    if(paths.size()==0 && smp->fSelectionEffPaths.size()>0) paths = smp->fSelectionEffPaths;
    if(files.size()==0 && smp->fSelectionEffFiles.size()>0) files = smp->fSelectionEffFiles;
    if(names.size()==0 && smp->fSelectionEffNames.size()>0) names = smp->fSelectionEffNames;

    if(paths.size()==0 && reg->fSelectionEffPaths.size()>0) paths = reg->fSelectionEffPaths;
    if(files.size()==0 && reg->fSelectionEffFiles.size()>0) files = reg->fSelectionEffFiles;
    if(names.size()==0 && reg->fSelectionEffNames.size()>0) names = reg->fSelectionEffNames;

    if(paths.size()==0 && fSelectionEffPaths.size()>0) paths = fSelectionEffPaths;
    if(files.size()==0 && fSelectionEffFiles.size()>0) files = fSelectionEffFiles;
    if(names.size()==0 && fSelectionEffNames.size()>0) names = fSelectionEffNames;

    if(syst) {
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffPathSuffs,smp->fSelectionEffPathSuffs), syst->fSelectionEffPathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffPathSuffs,smp->fSelectionEffPathSuffs), syst->fSelectionEffPathSuffsDown);
    }
    else{
        pathSuffs = Common::CombinePathSufs(reg->fSelectionEffPathSuffs, smp->fSelectionEffPathSuffs);
    }

    if(syst) {
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffFileSuffs, smp->fSelectionEffFileSuffs), syst->fSelectionEffFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffFileSuffs, smp->fSelectionEffFileSuffs), syst->fSelectionEffFileSuffsDown);
    }
    else{
        fileSuffs = Common::CombinePathSufs(reg->fSelectionEffFileSuffs, smp->fSelectionEffFileSuffs);
    }

    if(syst) {
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffNameSuffs, smp->fSelectionEffNameSuffs), syst->fSelectionEffNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffNameSuffs, smp->fSelectionEffNameSuffs), syst->fSelectionEffNameSuffsDown);
    }
    else{
      nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffNameSuffs, smp->fSelectionEffNameSuffs), fSelectionEffNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
void TRExFit::PlotUnfold(TH1D* data,
                         TGraphAsymmErrors* total) const {

    std::vector<std::unique_ptr<TH1> > mc;
    std::vector<std::string> legendNames;
    for (const auto& isample : fTruthSamples) {
        if (!isample->GetUseForPlotting()) continue;
        mc.emplace_back(std::move(isample->GetHisto(this)));
        mc.back()->SetLineColor(isample->GetLineColor());
        mc.back()->SetLineStyle(isample->GetLineStyle());
        legendNames.emplace_back(isample->GetTitle());
    }

    if (mc.empty()) {
        WriteWarningStatus("TRExFit::PlotUnfold", "No MC samples set for plotting. Will not create the final plots");
        return;
    }
    
    if(fUnfoldNormXSec){
        for (auto& itruth : mc) {
            itruth->Scale(1./itruth->Integral());
        }
    }

    if (fUnfoldingDivideByBinWidth) {
        Common::ScaleByBinWidth(data);
        Common::ScaleByBinWidth(total);
        for (auto& imc : mc) {
            Common::ScaleByBinWidth(imc.get());
        }
    }

    TCanvas c("","",600,600);
    TPad pad1("pad1","pad1",0.0, 0.3, 1.0, 1.00);
    TPad pad2("pad2","pad2", 0.0, 0.010, 1.0, 0.3);
    pad1.SetBottomMargin(0.001);
    pad1.SetBorderMode(0);
    pad2.SetBottomMargin(0.5);
    pad1.SetTicks(1,1);
    pad2.SetTicks(1,1);
    pad1.Draw();
    pad2.Draw();

    if (fUnfoldingLogX) {
        pad1.SetLogx();
        pad2.SetLogx();
    }
    if (fUnfoldingLogY) {
        pad1.SetLogy();
    }

    pad1.cd();
    total->SetMarkerStyle(20);
    total->SetMarkerSize(1.2);
    total->SetLineColor(kBlack);
    total->SetLineStyle(1);

    total->GetYaxis()->SetLabelSize(0.05);
    total->GetYaxis()->SetLabelFont(42);
    total->GetYaxis()->SetTitleFont(42);
    total->GetYaxis()->SetTitleSize(0.07);
    total->GetYaxis()->SetTitleOffset(1.1);

    std::unique_ptr<TH1> h_dummy(static_cast<TH1*>(mc[0]->Clone()));
    const double corr = fUnfoldingScaleRangeY > 0 ? fUnfoldingScaleRangeY : (fUnfoldingLogY ? 1e6 : 1.5);
    h_dummy->GetYaxis()->SetRangeUser(0.0001, corr*h_dummy->GetMaximum());
    h_dummy->GetYaxis()->SetTitle(fUnfoldingTitleY.c_str());
    h_dummy->GetYaxis()->SetTitleOffset(fUnfoldingTitleOffsetY*h_dummy->GetYaxis()->GetTitleOffset());
    h_dummy->SetLineWidth(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->DrawClone("HIST");

    std::unique_ptr<TGraphAsymmErrors> error(static_cast<TGraphAsymmErrors*>(total->Clone()));
    error->SetFillStyle(1001);
    error->SetFillColor(kGray);
    error->SetLineColor(kGray);
    error->Draw("E2 same");

    for (auto& itruth : mc) {
        itruth->SetLineWidth(2);
        itruth->Draw("HIST same");
    }

    total->Draw("pX SAME");

    float legH = (mc.size()+2)*0.07;
    TLegend leg(0.55, 0.9-legH, 0.9, 0.9);
    if (fAlternativeAsimovTruthSample == "") {
        leg.AddEntry(total, "Unfolded data", "p");
    } else {
        leg.AddEntry(total, "Unfolded pseudo-data", "p");
    }
    for (std::size_t imc = 0; imc < mc.size(); ++imc) {
        leg.AddEntry(mc.at(imc).get(), legendNames.at(imc).c_str(), "l");
    }
    leg.AddEntry(error.get(), "Total uncertainty","f");

    leg.SetFillColor(0);
    leg.SetLineColor(0);
    leg.SetBorderSize(0);
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetTextSize(gStyle->GetTextSize());
    leg.Draw("SAME");

    if (fAtlasLabel != "none") ATLASLabel(0.2,0.87,fAtlasLabel.c_str());
    myText(0.2,0.8,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    pad1.RedrawAxis();

    // Plot ratio
    pad2.cd();
    h_dummy->GetXaxis()->SetTitleOffset(fUnfoldingTitleOffsetX*h_dummy->GetXaxis()->GetTitleOffset());
    h_dummy->GetYaxis()->SetTitleOffset(fUnfoldingTitleOffsetY*0.55*h_dummy->GetYaxis()->GetTitleOffset());
    h_dummy->GetYaxis()->SetRangeUser(fUnfoldingRatioYmin,fUnfoldingRatioYmax);
    h_dummy->GetYaxis()->SetTitle("#frac{Prediction}{Data}");
    h_dummy->GetYaxis()->SetNdivisions(505);
    h_dummy->GetXaxis()->SetTitle(fUnfoldingTitleX.c_str());
    h_dummy->Draw("HIST");

    // Plot the error band
    std::unique_ptr<TGraphAsymmErrors> band = Common::GetRatioBand(total, data);
    band->SetFillStyle(1001);
    band->SetFillColor(kGray);
    band->SetLineColor(kGray);
    band->SetMarkerStyle(0);
    band->SetLineWidth(2);
    band->Draw("E2 SAME");

    // Plot the theory predictions
    std::vector<std::unique_ptr<TH1D> > ratios;
    for (const auto& itruth : mc) {
        ratios.emplace_back(static_cast<TH1D*>(itruth->Clone()));
        ratios.back()->Divide(data);
        ratios.back()->Draw("HIST same");
    }

    pad2.RedrawAxis();

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fName+"/UnfoldedData."+ format).c_str());
    }
}

//__________________________________________________________________________________
//
void TRExFit::PlotMigrationResponse(const TH2* matrix,
                                    const bool isMigration,
                                    const std::string& regionName,
                                    const std::string& systematicName) const {

    std::unique_ptr<TH2D> m(static_cast<TH2D*>(matrix->Clone()));

    gStyle->SetPalette(87);
    if (isMigration) gStyle->SetPaintTextFormat("1.2f");
    else gStyle->SetPaintTextFormat("2.1f");
    TCanvas c("","",600,600);
    c.cd();

    TPad pad("","",0,0.0,1,1);
    pad.SetRightMargin(0.15);
    pad.SetLeftMargin(0.15);
    pad.SetTopMargin(0.15);
    pad.SetBottomMargin(0.15);
    pad.Draw();
    pad.cd();

    if (fMigrationLogX) {
        pad.SetLogx();
    }
    if (fMigrationLogY) {
        pad.SetLogy();
    }

    m->SetMarkerSize(850);
    m->SetMarkerColor(kBlack);
    m->GetXaxis()->SetTitle(fMigrationTitleX.c_str());
    m->GetXaxis()->SetTitleOffset(fMigrationTitleOffsetX * m->GetXaxis()->GetTitleOffset());
    m->GetYaxis()->SetTitleOffset(fMigrationTitleOffsetY * m->GetYaxis()->GetTitleOffset());
    m->GetYaxis()->SetTitle(fMigrationTitleY.c_str());
    m->GetXaxis()->SetNdivisions(505);
    m->GetZaxis()->SetTitleOffset(1.3);
    if (isMigration) {
        m->GetZaxis()->SetTitle("Migration");
        m->GetZaxis()->SetRangeUser(fMigrationZmin,fMigrationZmax);
    } else {
        m->GetZaxis()->SetTitle("Response");
        m->GetZaxis()->SetRangeUser(fResponseZmin,fResponseZmax);
    }

    c.SetGrid();
    if(fMigrationText) m->Draw("COLZ TEXT");
    else m->Draw("COLZ");

    if (fAtlasLabel != "none") ATLASLabel(0.03,0.92,fAtlasLabel.c_str());
    myText(0.68,0.92,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    c.RedrawAxis("g");

    if (systematicName != "") {
        gSystem->mkdir("Systematics");
        gSystem->mkdir(("Systematics/"+systematicName).c_str());
    }
    const std::string tmp = isMigration ? "migration" : "response";
    const std::string name = systematicName == "" ? tmp + "_" + regionName : "Systematics/" + systematicName+"/"+tmp + "_" + regionName;

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fName+"/" + name+ "."+ format).c_str());
    }
}

//__________________________________________________________________________________
//
void TRExFit::RunForceShape() {
    for (const auto& ireg : fRegions) {
        for (const auto& ismp : fSamples) {
            if(Common::FindInStringVector(ismp->fRegions, ireg->fName) < 0) continue;
            std::shared_ptr<SampleHist> sh = ireg->GetSampleHist(ismp->fName);
            if(!sh) continue;
            for (const auto& isyst : ismp->fSystematics) {
                if (isyst->fForceShape == HistoTools::FORCESHAPETYPE::NOSHAPE) continue;
                if(isyst->fRegions.size()>0 && Common::FindInStringVector(isyst->fRegions,ireg->fName)<0  ) continue;
                if(isyst->fExclude.size()>0 && Common::FindInStringVector(isyst->fExclude,ireg->fName)>=0 ) continue;
                if(isyst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(isyst->fExcludeRegionSample, ireg->fName, ismp->fName)>=0 ) continue;
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic(isyst->fName);
                if(!syh) continue;

                HistoTools::ForceShape(syh->fHistUp.get(), sh->fHist.get(), isyst->fForceShape);
                // Symmetrise
                for (int ibin = 1; ibin <= sh->fHist->GetNbinsX(); ++ibin) {
                    const double diff = syh->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                    syh->fHistDown->SetBinContent(ibin, sh->fHist->GetBinContent(ibin) - diff);
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
bool TRExFit::DoingMixedFitting() const {
    if(fFitNPValues.size()>0) return false;
    if(fFitNPValuesFromFitResults!="") return false;
    int dataType = -1;
    // loop on regions
    for (auto reg : fRegions) {
        if (reg->fRegionType == Region::VALIDATION) continue;
        if (dataType >= 0 && reg->fRegionDataType != dataType) return true;
        dataType = reg->fRegionDataType;
    }
    return false;
}

//__________________________________________________________________________________
//
std::vector < std:: string > TRExFit::ListRegionsToFit(int dataType) const {
    std::vector < std:: string > list;
    for (auto reg : fRegions) {
        if (reg->fRegionType == Region::VALIDATION) continue;
        if (dataType >= 0 && reg->fRegionDataType != dataType) continue;
        if (fFitRegion == CRONLY && reg->fRegionType != Region::CONTROL) continue;
        list.emplace_back( reg->fName );
    }
    return list;
}

//__________________________________________________________________________________
//
std::map < std::string, int > TRExFit::MapRegionDataTypes(const std::vector<std::string>& regionList,bool isBlind) const {
    std::map < std::string, int > dataTypes;
    for (auto name : regionList) {
        if (isBlind) {
            dataTypes[name] = Region::ASIMOVDATA;
            continue;
        }
        Region* reg = GetRegion(name);
        if (!reg) {
            WriteErrorStatus("TRExFit::MapRegionDataTtpe", "Trying to access a non-exisiting region: " + name + "...");
            exit(EXIT_FAILURE);
        }
        dataTypes[name] = reg->fRegionDataType;
    }
    return dataTypes;
}

//__________________________________________________________________________________
//
std::map < std::string, double > TRExFit::NPValuesFromFitResults(const std::string& fitResultsFile) {
    std::map < std::string, double > npValues;
    ReadFitResults(fitResultsFile);
    for(const auto& inp : fFitResults->fNuisPar) {
        npValues[inp->fName] = inp->fFitValue;
    }
    return npValues;
}

//__________________________________________________________________________________
//
void TRExFit::DrawToyPullPlot(std::vector<TH1D>& hist, TFile* out) const {

    std::vector<double> mean;
    std::vector<double> sigma;
    std::vector<double> meanError;
    std::vector<double> sigmaError;

    for (auto& ihist : hist) {
        mean.emplace_back(ihist.GetMean());
        sigma.emplace_back(ihist.GetRMS());
        meanError.emplace_back(ihist.GetMeanError());
        sigmaError.emplace_back(ihist.GetRMSError());
        
        out->cd();
        ihist.Write();
    }

    TH1D h_mean("","",hist.size(), 0, hist.size());
    TH1D h_sigma("","",hist.size(), 0, hist.size());

    for (std::size_t i = 0; i < hist.size(); ++i) {
        h_mean.SetBinContent(i+1, mean.at(i));
        h_mean.SetBinError(i+1, meanError.at(i));
        h_sigma.SetBinContent(i+1, sigma.at(i));
        h_sigma.SetBinError(i+1, sigmaError.at(i));

        h_mean.GetXaxis()->SetBinLabel(i+1, ("bin " + std::to_string(i+1)).c_str());
    }

    h_mean.GetYaxis()->SetRangeUser(-1, 3);
    h_mean.GetYaxis()->SetTitle("#frac{#mu_{fitted}-1}{#sigma_{fitted}}");
    h_mean.LabelsOption("v");
    h_mean.SetLineColor(kBlue);
    h_mean.SetMarkerColor(kBlue);
    h_sigma.SetLineColor(kRed);
    h_sigma.SetMarkerColor(kRed);
    
    TCanvas c("c","c",800,600);
    h_mean.Draw("PE");
    h_sigma.Draw("PE SAME");

    TLine l_mean(0, 0, hist.size(), 0);
    TLine l_sigma(0, 1, hist.size(), 1);

    l_mean.SetLineColor(kBlue);
    l_mean.SetLineWidth(2);
    l_mean.SetLineStyle(2);
    l_sigma.SetLineColor(kRed);
    l_sigma.SetLineWidth(2);
    l_sigma.SetLineStyle(2);

    l_mean.Draw("same");
    l_sigma.Draw("same");

    TLegend leg(0.7, 0.8, 0.9, 0.9);
    leg.SetFillStyle(0.);
    leg.SetBorderSize(0.);
    leg.AddEntry(&h_mean, "mean", "p");
    leg.AddEntry(&h_sigma, "sigma", "p");
    leg.Draw("same");
        
    if (fAtlasLabel != "none") ATLASLabel(0.3,0.85,fAtlasLabel.c_str());

    gPad->RedrawAxis();

    out->cd();
    h_mean.Write("pull_mean");
    h_sigma.Write("pull_sigma");

    for(const auto& format : TRExFitter::IMAGEFORMAT) {
        c.SaveAs((fName+"/Toys/PullPlot."+format).c_str());
    }
}

//__________________________________________________________________________________
//
void TRExFit::FixUnfoldingExpressions() {
    if (!fUnfoldNormXSec) return;
    // Get truth histogram
    std::unique_ptr<TFile> input(TFile::Open((fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
    if (!input) {
        WriteErrorStatus("TRExFit::FixUnfoldingExpressions", "Cannot read file from " + fName + "/UnfoldingHistograms/FoldedHistograms.root");
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<TH1D> truth(dynamic_cast<TH1D*>(input->Get("truth_distribution")));
    if (!truth) {
        WriteErrorStatus("TRExFit::FixUnfoldingExpressions", "Cannot read the truth distribution");
        exit(EXIT_FAILURE);
    }
    truth->SetDirectory(nullptr);
    const double N = truth->Integral();
    // Get expression of last bin norm-factor
    const std::string NFname = "Bin_" + Common::IntToFixLenStr(fUnfoldNormXSecBinN) + "_mu";
    for (const auto& inorm : fNormFactors) {
        if (inorm->fName==NFname) {
            TString expression(inorm->fExpression.first);
            for (int i = 0; i < fNumberUnfoldingTruthBins; ++i) {
                // replace all "Nj/N" with truth->bincontent/truth->integral
                expression.ReplaceAll("N"+std::to_string(i+1)+"/N",std::to_string(truth->GetBinContent(i+1)/N));
            }
            inorm->fExpression.first = expression.Data();
            // title will contain the expression FIXME
            inorm->fTitle = inorm->fExpression.first;
            TRExFitter::SYSTMAP[inorm->fName] = inorm->fExpression.first;
        }
    }
    input->Close();
}
