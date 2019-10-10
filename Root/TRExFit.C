// Class include
#include "TRExFitter/TRExFit.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExPlot.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/PruningUtil.h"

// CommonStatTiils includes
#include "CommonStatTools/runSig.h"
#include "CommonStatTools/runAsymptoticsCLs.h"

//Roofit headers
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooMinimizer.h"

//HistFactory headers
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"

// ATLAS stuff
#include "AtlasUtils/AtlasStyle.h"
#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

// ROOT includes
#include "Math/DistFunc.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFormula.h"
#include "TGaxis.h"
#include "TGraph2D.h"
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
TRExFit::TRExFit(std::string name){
    fSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
    fDir = "";
    fName = name;
    fInputName = name;
    fLabel = "";
    fCmeLabel = "13 TeV";
    fLumiLabel = "3.2 fb^{-1}";

    fNRegions = 0;
    fNSamples = 0;
    fNSyst = 0;

    fPOI = "";
    fPOIunit = "";
    fUseStatErr = false;
    fStatErrThres = 0.05;
    fUseGammaPulls = false;

    fLumi = 1.;
    fLumiErr = 0.000001;
    fLumiScale = 1.;

    fThresholdSystPruning_Normalisation = -1;
    fThresholdSystPruning_Shape = -1;
    fThresholdSystLarge = - 1;

    fNtuplePaths.clear();
    fNtupleFiles.clear();
    fNtupleNames.clear();
    fMCweight = "1";
    fSelection = "1";

    fHistoPaths.clear();
    fHistoFiles.clear();
    fHistoNames.clear();
    fHistoNamesNominal.clear();

    fFitResults = 0;

    fWithPullTables = false;

    fRegions.clear();
    fSamples.clear();
    fSystematics.clear();

    fIntCode_overall = 4;
    fIntCode_shape = 0;

    fConfig = new ConfigParser();

    fInputType = HIST;

    fSuffix = "";
    fSaveSuffix = "";

    fUpdate = false;
    fKeepPruning = false;

    fBlindingThreshold = -1;
    fBlindingType = SOVERB;

    fRankingMaxNP = 10;
    fRankingOnly = "all";
    fRankingPlot = "Merge";
    fAtlasLabel = "Internal";

    fStatOnly = false;
    fGammasInStatOnly = false;
    fStatOnlyFit = false;
    fFixNPforStatOnlyFit = false;
    fSystDataPlot_upFrame = false;

    fSummaryPlotRegions.clear();
    fSummaryPlotLabels.clear();

    fSummaryPlotValidationRegions.clear();
    fSummaryPlotValidationLabels.clear();

    fYmin = 0;
    fYmax = 0;

    fFitResultsFile = "";

    fDoSummaryPlot = true;
    fDoMergedPlot = false;
    fDoTables = true;
    fDoSignalRegionsPlot = true;
    fDoPieChartPlot = true;

    fSummaryPrefix = "";

    fGroupedImpactCategory = "all";

    //
    // Fit caracteristics
    //
    fFitType = UNDEFINED;
    fFitRegion = CRSR;
    fFitRegionsToFit.clear();
    fFitNPValues.clear();
    fFitPOIAsimov = 0;
    fFitIsBlind = false;
    fUseRnd = false;
    fRndRange = 0.1;
    fRndSeed = -999;
    fVarNameLH.clear();
    fLHscanMin = 999999;
    fLHscanMax = -999999;
    fLHscanSteps = 30;
    fParal2D = false;
    fParal2Dstep = -1;
    fLHscanMinY = 999999;
    fLHscanMaxY = -999999;
    fLHscanStepsY = 30;
    fVarNameMinos.clear();
    fVarNameHide.clear();
    fWorkspaceFileName = "";
    fDoGroupedSystImpactTable = false;
    fSubCategoryImpactMap.clear();

    //
    // Limit type
    //
    fLimitType = ASYMPTOTIC;
    fLimitIsBlind = false;
    fLimitPOIAsimov = 0;
    fSignalInjection = false;
    fSignalInjectionValue = 0;
    fLimitParamName = "parameter";
    fLimitParamValue = 0;
    fLimitOutputPrefixName = "myLimit";
    fLimitsConfidence = 0.95;

    //
    // Significance parameters
    //
    fSignificanceIsBlind = false;
    fSignificancePOIAsimov = 0;
    fSignificanceParamName = "parameter";
    fSignificanceParamValue = 0;
    fSignificanceOutputPrefixName = "mySignificance";

    fImageFormat = "png";
    TRExFitter::IMAGEFORMAT.clear();
    TRExFitter::IMAGEFORMAT.push_back("png");

    //
    fSystematics.clear();
    fSystematicNames.clear();
    fNSyst = 0;
    //
    fNormFactors.clear();
    fNormFactorNames.clear();
    fNNorm = 0;
    //
    fShapeFactors.clear();
    fShapeFactorNames.clear();
    fNShape = 0;

    fCleanTables = false;
    fSystCategoryTables = false;

    fRegionGroups.clear();

    // Increase the limit for formula evaluations
    ROOT::v5::TFormula::SetMaxima(100000,1000,1000000);

    fKeepPrefitBlindedBins = false;
    fBlindedBins = nullptr;

    fRatioYmax = 1.5;
    fRatioYmin = 0.5;
    fRatioYmaxPostFit = 1.5;
    fRatioYminPostFit = 0.5;
    fRatioYtitle = "";
    fRatioType = "DATA/MC";

    fCustomAsimov = "";
    fTableOptions = "STANDALONE";

    fGetGoodnessOfFit = false;
    fGetChi2 = 0; // 0: no, 1: stat-only, 2: with syst

    fCustomFunctions.clear();
    fSuppressNegativeBinWarnings = false;

    fMorphParams.clear();
    fTemplateInterpolationOption = TRExFit::LINEAR;

    fBootstrap = "";
    fBootstrapSyst = "";
    fBootstrapIdx = -1;

    fDecorrSysts.clear();
    fDecorrSuff = "_decor";

    fDoNonProfileFit = false;
    fNonProfileFitSystThreshold = 0;
    fFitToys = 0;
    fToysHistoMin = 9999;
    fToysHistoMax = -9999;
    fToysHistoNbins = 50;
    fToysPseudodataNP = "";
    fToysPseudodataNPShift = 1.;
    fSmoothMorphingTemplates = "";

    fPOIPrecision = 2;

    fRankingPOIName = "#mu";

    fUseATLASRoundingTxt = false;
    fUseATLASRoundingTex = false;

    fuseGammasForCorr = false;
    fPropagateSystsForMorphing = false;
    fPruningType = SEPARATESAMPLE;

    fLabelX = -1;
    fLabelY = -1;
    fLegendX1 = -1;
    fLegendX2 = -1;
    fLegendY = -1;

    fLabelXSummary = -1;
    fLabelYSummary = -1;
    fLegendX1Summary = -1;
    fLegendX2Summary = -1;
    fLegendYSummary = -1;

    fLabelXMerge = -1;
    fLabelYMerge = -1;
    fLegendX1Merge = -1;
    fLegendX2Merge = -1;
    fLegendYMerge = -1;

    fLegendNColumns = 2;
    fLegendNColumnsSummary = 3;
    fLegendNColumnsMerge = 3;

    fShowRatioPad = true;

    fExcludeFromMorphing = "";

    fSaturatedModel = false;

    fDebugNev = -1;
}

//__________________________________________________________________________________
//
TRExFit::~TRExFit(){
    if(fFitResults) delete fFitResults;

    for(unsigned int i =0 ; i < fRegions.size(); ++i){
        if(fRegions[i]){
            delete fRegions[i];
        }
    }
    fRegions.clear();

    for(unsigned int i =0 ; i < fSamples.size(); ++i){
        if(fSamples[i]){
            delete fSamples[i];
        }
    }
    fSamples.clear();
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
Sample* TRExFit::NewSample(const std::string& name,int type){
    fSamples.push_back(new Sample(name,type));
    //
    fNSamples ++;
    return fSamples[fNSamples-1];
}

//__________________________________________________________________________________
//
Systematic* TRExFit::NewSystematic(const std::string& name){
    fSystematics.push_back(new Systematic(name));
    fNSyst ++;
    return fSystematics[fNSyst-1];
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
                fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(fSmoothOption, syst, false);
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
                        fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(fSmoothOption, fSystematics.at(i_syst)->fName, true);
                    } else {
                        // check if the syst has been smoothed already
                        if (std::find(usedSysts.begin(), usedSysts.end(), i_syst) != usedSysts.end()) continue;
                        // Need to apply special smoothing
                        // smooth the reference sample
                        SampleHist *sh = GetSampleHistFromName(fRegions[i_ch], fSystematics.at(i_syst)->fReferenceSmoothing);
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
                                sh->SmoothSyst(fSmoothOption, fSystematics.at(i_syst)->fName, true);

                                // save the smoothed histograms
                                nominal_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(fRegions[i_ch]->fSampleHists[i_sample]->fHist->Clone()));
                                systIndex = GetSystIndex(fRegions[i_ch]->fSampleHists[i_sample], fSystematics.at(i_syst)->fName);
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
                            delete fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistUp;
                            delete fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistDown;
                            fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistUp = CopySmoothedHisto(fRegions[i_ch]->fSampleHists[i_sample],nominal_cpy.get(),up_cpy.get(),down_cpy.get(),true);
                            fRegions[i_ch]->fSampleHists[i_sample]->fSyst[systIndex]->fHistDown = CopySmoothedHisto(fRegions[i_ch]->fSampleHists[i_sample],nominal_cpy.get(),up_cpy.get(),down_cpy.get(),false);
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
    std::string fileName = "";
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    //
    if(singleOutputFile){
        if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos" + fSaveSuffix + ".root";
        else                 fileName = fName + "/Histograms/" + fInputName + "_histos" + fSaveSuffix + ".root";
        // Bootstrap
        if(fBootstrap!="" && fBootstrapIdx>=0){
            fileName = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fBootstrapIdx));
        }
        WriteInfoStatus("TRExFit::CreateRootFiles","-------------------------------------------");
        WriteInfoStatus("TRExFit::CreateRootFiles","Creating/updating file " + fileName + " ...");
        if(recreate) fFiles.push_back(new TFile(fileName.c_str(),"RECREATE"));
        else         fFiles.push_back(new TFile(fileName.c_str(),"UPDATE"));
        TRExFitter::TFILEMAP.insert(std::make_pair(fileName,fFiles[fFiles.size()-1]));
    }
    else{
        for(int i_ch=0;i_ch<fNRegions;i_ch++){
            if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
            else                 fileName = fName + "/Histograms/" + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
            // Bootstrap
            if(fBootstrap!="" && fBootstrapIdx>=0){
                fileName = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fBootstrapIdx));
            }
            WriteInfoStatus("TRExFit::CreateRootFiles","-------------------------------------------");
            WriteInfoStatus("TRExFit::CreateRootFiles","Creating/updating file " + fileName + " ...");
            if(recreate) fFiles.push_back(new TFile(fileName.c_str(),"RECREATE"));
            else         fFiles.push_back(new TFile(fileName.c_str(),"UPDATE"));
            TRExFitter::TFILEMAP.insert(std::make_pair(fileName,fFiles[fFiles.size()-1]));
        }
    }
}

//__________________________________________________________________________________
// fill files with all the histograms
void TRExFit::WriteHistos(bool reWriteOrig) const{
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    SampleHist* sh;
    std::string fileName = "";
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        //
        if(singleOutputFile) fileName = fFiles[0]   ->GetName();
        else                 fileName = fFiles[i_ch]->GetName();
        WriteInfoStatus("TRExFit::WriteHistos","-------------------------------------------");
        WriteInfoStatus("TRExFit::WriteHistos","Writing histograms to file " + fileName + " ...");
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if(!sh){
                WriteDebugStatus("TRExFit::WriteHistos", "SampleHist[" + std::to_string(i_smp) + "] for sample " + fSamples[i_smp]->fName + " not there.");
                continue;
            }
            // set file and histo names for nominal
            sh->fHistoName = sh->fHist->GetName();
            sh->fFileName = fileName;
            // set file and histo names for systematics
            for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
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
        for(auto sh : reg->fSampleHists){
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
    SampleHist* sh;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
            sh = fRegions[i_ch]->fSampleHists[i_smp];
            sh->DrawSystPlot("all");
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots for sombined samples
void TRExFit::DrawSystPlotsSumSamples() const{
    WriteInfoStatus("TRExFit::DrawSystPlotsSumSamples", "-------------------------------------------");
    WriteInfoStatus("TRExFit::DrawSystPlotsSumSamples", "Drawing combined plots of syst effects on data...");
    TH1* h_dataCopy = nullptr;
    for(auto reg : fRegions){
        SampleHist* hist = new SampleHist();
        bool empty = true;
        std::set<std::string> systNames;
        for(int i_regSmp=0; i_regSmp<reg->fNSamples; i_regSmp++){
            for(int i_smSyst=0; i_smSyst<reg->fSampleHists[i_regSmp]->fNSyst; i_smSyst++){
                systNames.insert(reg->fSampleHists[i_regSmp]->fSyst[i_smSyst]->fName);
            }
        }
        for(int i_smp=0;i_smp<reg->fNSamples;i_smp++){
            if(reg->fSampleHists[i_smp]->fSample->fType==Sample::DATA) h_dataCopy=(TH1*)reg->fSampleHists[i_smp]->fHist->Clone();
            else if(reg->fSampleHists[i_smp]->fSample->fType==Sample::GHOST) continue;
            else {
                double scale = GetNominalMorphScale(reg->fSampleHists[i_smp]);
                if(empty){
                    hist->CloneSampleHist(reg->fSampleHists[i_smp],systNames, scale);
                    hist->fName = reg->fName + "_Combined";
                    empty=false;
                } else {
                    hist->SampleHistAdd(reg->fSampleHists[i_smp], scale);
                }
            }
        }
        hist->DrawSystPlot("all", h_dataCopy, true, fSystDataPlot_upFrame);
        delete hist;
    }
}

//__________________________________________________________________________________
// for each region, add a SampleHist for each Sample in the Fit, reading from ntuples
void TRExFit::ReadNtuples(){
    WriteInfoStatus("TRExFit::ReadNtuples", "-------------------------------------------");
    WriteInfoStatus("TRExFit::ReadNtuples", "Reading ntuples...");
    TH1D* h = nullptr;
    TH1D* hUp = nullptr;
    TH1D* hDown = nullptr;
    std::string variable;
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;
    std::vector<std::string> empty;
    SampleHist *sh;
    //
    // Import custom functions from .C files
    //
    for(auto file : fCustomFunctions){
        WriteInfoStatus("TRExFit::ReadNtuples", "  Loading function from " + file + " ...");
        gROOT->ProcessLineSync((".L "+file+"+").c_str());
    }
    //
    // Loop on regions
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        WriteInfoStatus("TRExFit::ReadNtuples", "  Region region " + fRegions[i_ch]->fName + " ...");
        //
        if(TRExFitter::SPLITHISTOFILES) fFiles[i_ch]->cd();
        //
        if(fRegions[i_ch]->fBinTransfo != "") ComputeBinning(i_ch);
        if(fRegions[i_ch]->fCorrVar1 != ""){
            if(fRegions[i_ch]->fCorrVar2 == ""){
                WriteWarningStatus("TRExFit::ReadNtuples", "Only first correlation variable defined, do not read region : " + fRegions[i_ch]->fName);
                continue;
            }
            WriteDebugStatus("TRExFit::ReadNtuples", "calling the function 'DefineVariable(i_ch)'");
            DefineVariable(i_ch);
        }
        else if(fRegions[i_ch]->fCorrVar2 != ""){
            WriteWarningStatus("TRExFit::ReadNtuples", "Only second correlation variable defined, do not read region : " + fRegions[i_ch]->fName);
            continue;
        }
        //
        // Loop on samples
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            WriteInfoStatus("TRExFit::ReadNtuples","    Reading sample " + fSamples[i_smp]->fName);
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            // set variables, selection, weight and paths
            variable      = Variable(       fRegions[i_ch],fSamples[i_smp]);
            fullSelection = FullSelection(  fRegions[i_ch],fSamples[i_smp]);
            fullMCweight  = FullWeight(     fRegions[i_ch],fSamples[i_smp]);
            fullPaths     = FullNtuplePaths(fRegions[i_ch],fSamples[i_smp]);
            //
            h = nullptr;
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                TH1D* htmp = nullptr;
                if(fRegions[i_ch]->fHistoBins){
                    htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                 variable, fRegions[i_ch]->fHistoNBinsRebin, fRegions[i_ch]->fHistoBins,
                                                 fullSelection, fullMCweight, fDebugNev);
                }
                else{
                    htmp = HistFromNtuple( fullPaths[i_path],
                                           variable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                           fullSelection, fullMCweight, fDebugNev);
                    //Pre-processing of histograms (rebinning, lumi scaling)
                    if(fRegions[i_ch]->fHistoNBinsRebin != -1){
                        htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin);
                    }
                }
                //
                if(fSamples[i_smp]->fNormalizedByTheory && fSamples[i_smp]->fType!=Sample::DATA) htmp -> Scale(fLumi);
                //
                if(fSamples[i_smp]->fLumiScales.size()>i_path)  htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                //
                if(i_path==0) h = (TH1D*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                delete htmp;
            }
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );
            //
            // Importing the histogram in TRExFitter
            sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
            sh->fHist_orig = h_orig;
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            //
            //  -----------------------------------
            //
            // read norm factors
            for(int i_norm=0;i_norm<fSamples[i_smp]->fNNorm;i_norm++){
                NormFactor *nf = fSamples[i_smp]->fNormFactors[i_norm];
                //
                // eventually skip norm factor / region combination
                if( nf->fRegions.size()>0 && FindInStringVector(nf->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( nf->fExclude.size()>0 && FindInStringVector(nf->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadNtuples", "Adding norm " + nf->fName);
                //
                sh->AddNormFactor( nf );
            }

            //
            //  -----------------------------------
            //
            // read shape factors
            for(int i_shape=0;i_shape<fSamples[i_smp]->fNShape;i_shape++){
                ShapeFactor *sf = fSamples[i_smp]->fShapeFactors[i_shape];
                //
                // eventually skip shape factor / region combination
                if( sf->fRegions.size()>0 && FindInStringVector(sf->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( sf->fExclude.size()>0 && FindInStringVector(sf->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadNtuples", "Adding shape " + sf->fName);
                //
                sh->AddShapeFactor( sf );
            }

            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                Systematic * syst = fSamples[i_smp]->fSystematics[i_syst];
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadNtuples", "Adding syst " + syst->fName);
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL){
                    SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddOverallSyst(syst->fName,syst->fOverallUp,syst->fOverallDown);
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
                    SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddStatSyst(syst->fName,syst->fBins[0]);
                    syh->fSystematic = syst;
                    continue;
                }
                // else ...
                if(FindInStringVector(syst->fDummyForSamples,smp->fName)>=0){
                    WriteInfoStatus("TRExFit::ReadNtuples", "Systematic " + syst->fName + " set as dummy for sample " + smp->fName + " (region " + reg->fName + ")");
                    hUp   = (TH1D*)sh->fHist->Clone(Form("h_%s_%s_%sUp",  reg->fName.c_str(),smp->fName.c_str(),syst->fStoredName.c_str()));
                    hDown = (TH1D*)sh->fHist->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),smp->fName.c_str(),syst->fStoredName.c_str()));
                    SystematicHist *syh = sh->AddHistoSyst(syst->fName,hUp,hDown);
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
                    if(GetSample(syst->fReferenceSample)!=nullptr) smp = GetSample(syst->fReferenceSample);
                    else{
//                         WriteWarningStatus("TRExFit::ReadNtuples", "Reference Sample "+syst->fReferenceSample+" for systematic "+syst->fName+", sample "+fSamples[i_smp]->fName+", region "+fRegions[i_ch]->fName+" not found. Ignoring.");
                        WriteErrorStatus("TRExFit::ReadNtuples", "Reference sample: " + syst->fReferenceSample + " does not exist for region: " + reg->fName + ". Please check this!");
                        WriteErrorStatus("TRExFit::ReadNtuples", "This probably means that you run over a specific sample, you need to run over the reference sample as well.");
                        WriteErrorStatus("TRExFit::ReadNtuples", "Ignoring SeparateSample setting.");
                    }
                }
                //
                // Up
                //
                hUp = nullptr;
                if(syst->fHasUpVariation){
                    // set variables, selection, weight and paths
                    fullMCweight  = FullWeight(     fRegions[i_ch],fSamples[i_smp],syst,true);
                    fullPaths     = FullNtuplePaths(fRegions[i_ch],fSamples[i_smp],syst,true);
                    WriteDebugStatus("TRExFit::ReadNtuples", "  Syst Up full weight: " + fullMCweight);
                    for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                        TH1D* htmp = nullptr;
                        if(reg->fHistoBins){
                            htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                        variable, reg->fHistoNBinsRebin, reg->fHistoBins,
                                                        fullSelection, fullMCweight, fDebugNev);
                        }
                        else{
                            htmp = HistFromNtuple( fullPaths[i_path],
                                                  variable, reg->fNbins, reg->fXmin, reg->fXmax,
                                                  fullSelection, fullMCweight, fDebugNev);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1){
                                htmp->Rebin(reg->fHistoNBinsRebin);
                            }
                        }
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        
                        //
                        // Importing histogram in TRExFitter
                        if(i_path==0){
                            hUp = (TH1D*)htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hUp->Add(htmp);
                        delete htmp;
                    } // end loop over files
                    
                    
                    // BW
                    // pulled this out of the file loop to apply it only to the fully constructed histogram insead of file by file
                        
                    // obtain relative variation and apply it to proper sample
                    // & try to keep also the same total relative variation
                    if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr){
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                        TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                           
                        // Protection added: fix empty bins before starting to divide and multiply 
                            
                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin< hUp->GetNbinsX()+2;i_bin++) if(hUp ->GetBinContent(i_bin)<=1e-6) hUp ->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) hUp ->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                        //
                        double relVar   = hUp->Integral(0,hUp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                            
                        // get copies with no error
                        auto hrefTmp = GetHistCopyNoError(href);
                        auto hnomTmp = GetHistCopyNoError(hnom);
                        hUp->Divide(   hrefTmp.get() );
                        hUp->Multiply( hnomTmp.get() );
                        double newVar   = hUp->Integral(0,hUp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                        if( syst->fKeepReferenceOverallVar && TMath::Abs(relVar-1) > 0.0001 && TMath::Abs(newVar) > 0.0001) hUp->Scale( relVar / newVar );
                    }
                    // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                    else if (syst->fReferenceSample!="" && syst->fReferenceSample!=fSamples[i_smp]->fName && syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr) {
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                        TH1* href_up = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistUp;
                        TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                        
                        // Protection added: fix empty bins before starting to divide and multiply
                        
                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin< hUp->GetNbinsX()+2;i_bin++) if( hUp->GetBinContent(i_bin)<=1e-6) hUp->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) hUp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6

                        // Formula: UpHisto = [1+(up-nom)/nom-(DataUp-Data)/Data]*nom = up+nom+DataUp/Data*nom
                        TH1* href_up_Tmp = (TH1*)href_up->Clone(Form("%s_Tmp", href_up->GetName()));
                        // get copies with no error
                        auto hrefTmp = GetHistCopyNoError(href);
                        auto hnomTmp = GetHistCopyNoError(hnom);
                        href_up_Tmp->Divide(hrefTmp.get());
                        href_up_Tmp->Multiply(hnomTmp.get());
                        hUp->Add(hnomTmp.get());
                        auto href_up_TmpNoError = GetHistCopyNoError(href_up_Tmp);
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
                    fullMCweight  = FullWeight(     fRegions[i_ch],fSamples[i_smp],syst,false);
                    fullPaths     = FullNtuplePaths(fRegions[i_ch],fSamples[i_smp],syst,false);
                    for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                        TH1D* htmp = nullptr;
                        if(reg->fHistoBins){
                            htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                        variable, reg->fHistoNBinsRebin, reg->fHistoBins,
                                                        fullSelection, fullMCweight, fDebugNev);
                        }
                        else{
                            htmp = HistFromNtuple( fullPaths[i_path],
                                                  variable, reg->fNbins, reg->fXmin, reg->fXmax,
                                                  fullSelection, fullMCweight, fDebugNev);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1){
                                htmp->Rebin(reg->fHistoNBinsRebin);
                            }
                        }
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);

                        //
                        // Importing histogram in TRExFitter
                        if(i_path==0){
                            hDown = (TH1D*)htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
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
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                        TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                        
                        // Protection added: fix empty bins before starting to divide and multiply
                        
                        for(int i_bin=0;i_bin< href->GetNbinsX()+2;i_bin++) if(href ->GetBinContent(i_bin)<=1e-6) href ->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin<hDown->GetNbinsX()+2;i_bin++) if(hDown->GetBinContent(i_bin)<=1e-6) hDown->SetBinContent(i_bin,1e-6);
                        for(int i_bin=0;i_bin< href->GetNbinsX()+2;i_bin++) if(href ->GetBinContent(i_bin)<=1e-6) hDown->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                        //
                        double relVar   = hDown->Integral(0,hDown->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                        hDown->Divide(   href );
                        hDown->Multiply( hnom );
                        double newVar   = hDown->Integral(0,hDown->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                        if( syst->fKeepReferenceOverallVar && TMath::Abs(relVar-1) > 0.0001 && TMath::Abs(newVar-1) > 0.0001) hDown->Scale( relVar / newVar );
                    }
                    // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                    else if (syst->fReferenceSample!="" && syst->fReferenceSample!=fSamples[i_smp]->fName && syst->fSubtractRefSampleVar && reg->GetSampleHist(syst->fReferenceSample)!=nullptr) {
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                        TH1* href_down = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistDown;
                        TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                        
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
                if(hUp==nullptr)   hUp   = (TH1D*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                if(hDown==nullptr) hDown = (TH1D*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                //
                SystematicHist *syh = sh->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                syh->fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                syh->fScaleUp = fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0)
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName]!=0)
                        syh->fScaleUp *= fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName];
                syh->fScaleDown = fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0)
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName]!=0)
                        syh->fScaleDown *= fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName];
            }
            // END SYST LOOP
            //
        }
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
            if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            SampleHist *sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            TH1* h_orig = (TH1*)sh->fHist_orig;
            TH1* h = nullptr;
            if(h_orig!=nullptr) h = (TH1*)h_orig->Clone(sh->fHist->GetName());
            sh->fHist = h;
            if(sh->fHist==nullptr) continue;
            sh->fHist->SetLineColor(linecolor);
            sh->fHist->SetFillColor(fillcolor);
            //
            // loop on systematics
            for(auto syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                // skip also separate gamma systs
                if(syst->fName.find("stat_")!=std::string::npos) continue;
                //
                // get the original syst histograms & reset the syst histograms
                SystematicHist *syh = sh->GetSystematic( syst->fName );
                //
                if(syh==nullptr) continue;
                TH1* hUp_orig   = syh->fHistUp_orig;
                TH1* hDown_orig = syh->fHistDown_orig;
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
                if(hUp_orig!=nullptr)   hUp   = (TH1*)hUp_orig  ->Clone( syh->fHistUp  ->GetName() );
                if(hDown_orig!=nullptr) hDown = (TH1*)hDown_orig->Clone( syh->fHistDown->GetName() );
                syh->fHistUp   = hUp;
                syh->fHistDown = hDown;
            }
        }
        //
        // 2. Rebin
        for(auto smp : fSamples){
            //
            // eventually skip sample / region combination
            if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            SampleHist *sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            //
            // Rebinning (FIXME: better to introduce a method Region::Rebin() ?)
            if(reg->fHistoNBinsRebinPost>0){
                WriteDebugStatus("TRExFit::CorrectHistograms", "rebinning " + smp->fName + " to " + std::to_string(reg->fHistoNBinsRebinPost) + " bins.");
                sh->fHist = sh->fHist->Rebin(reg->fHistoNBinsRebinPost,"",reg->fHistoBinsPost);
                for(auto syh : sh->fSyst){
                    WriteDebugStatus("TRExFit::CorrectHistograms", "  systematic " + syh->fName + " to " + std::to_string(reg->fHistoNBinsRebinPost) + " bins.");
                    if(syh==nullptr) continue;
                    if(syh->fSystematic->fSampleUp==""   && syh->fSystematic->fHasUpVariation   && syh->fHistUp!=nullptr)   syh->fHistUp   = syh->fHistUp  ->Rebin(reg->fHistoNBinsRebinPost,"",reg->fHistoBinsPost);
                    else                                                                                                    syh->fHistUp   = (TH1*)sh->fHist->Clone(syh->fHistUp->GetName());
                    if(syh->fSystematic->fSampleDown=="" && syh->fSystematic->fHasDownVariation && syh->fHistDown!=nullptr) syh->fHistDown = syh->fHistDown->Rebin(reg->fHistoNBinsRebinPost,"",reg->fHistoBinsPost);
                    else                                                                                                    syh->fHistDown = (TH1*)sh->fHist->Clone(syh->fHistDown->GetName());
                }
                //
                // rebin also separate-gamma hists!
                if(smp->fSeparateGammas){
                    SystematicHist *syh = sh->GetSystematic( "stat_"+smp->fName );
                    if(syh==nullptr) continue;
                    if(syh->fHistUp!=nullptr)   syh->fHistUp  ->Rebin(reg->fHistoNBinsRebinPost,"",reg->fHistoBinsPost);
                    if(syh->fHistDown!=nullptr) syh->fHistDown->Rebin(reg->fHistoNBinsRebinPost,"",reg->fHistoBinsPost);
                }
            }
        }

        // Randomize MC (before add/multiply/scale)
        if(TRExFitter::OPTION["RandomizeMC"]!=0){
            gRandom->SetSeed(TRExFitter::OPTION["RandomizeMC"]);
            for(auto smp : fSamples){
                //
                // eventually skip sample / region combination
                if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                //
                SampleHist *sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                if(sh->fHist==nullptr) continue;
                //
                if(smp->fUseMCStat){
                    TH1* hTmp = sh->fHist;
                    for(int i_bin=1;i_bin<=hTmp->GetNbinsX();i_bin++){
                        hTmp->SetBinContent(i_bin,gRandom->Poisson( hTmp->GetBinContent(i_bin) ));
                    }
                    for(auto syh : sh->fSyst){
                        for(int i_ud=0;i_ud<2;i_ud++){
                            if(i_ud==0) hTmp = syh->fHistUp;
                            else        hTmp = syh->fHistDown;
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
            if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            SampleHist *sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            //
            // Subtraction / Addition of sample
            for(auto sample : smp->fSubtractSamples){
                WriteDebugStatus("TRExFit::CorrectHistograms"," subtracting sample " + sample + " from sample " + smp->fName);
                SampleHist *smph0 = reg->GetSampleHist(sample);
                if(smph0!=nullptr) sh->Add(smph0,-1);
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+sample+" not found ...");
            }
            for(auto sample : smp->fAddSamples){
                WriteDebugStatus("TRExFit::CorrectHistograms", "adding sample " + sample + " to sample " + smp->fName);
                SampleHist *smph0 = reg->GetSampleHist(sample);
                if(smph0!=nullptr) sh->Add(smph0);
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+sample+" not found ...");
            }
            // Division & Multiplication by other samples
            if(smp->fMultiplyBy!=""){
                WriteDebugStatus("TRExFit::CorrectHistograms", "multiplying " + smp->fName  + " by sample " + smp->fMultiplyBy);
                SampleHist *smph0 = reg->GetSampleHist(smp->fMultiplyBy);
                if(smph0!=nullptr) sh->Multiply(smph0);
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+smp->fMultiplyBy+" not found ...");
            }
            if(smp->fDivideBy!=""){
                WriteDebugStatus("TRExFit::CorrectHistograms", "dividing " + smp->fName  + " by sample " + smp->fDivideBy + " from sample " + smp->fName);
                SampleHist *smph0 = reg->GetSampleHist(smp->fDivideBy);
                if(smph0!=nullptr) sh->Divide(smph0);
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+smp->fDivideBy+" not found ...");
            }
            // Norm to sample
            if(smp->fNormToSample!=""){
                WriteDebugStatus("TRExFit::CorrectHistograms", "normalizing " + smp->fName  + " to sample " + smp->fNormToSample);
                SampleHist *smph0 = reg->GetSampleHist(smp->fNormToSample);
                if(smph0!=nullptr) sh->Scale(smph0->fHist->Integral()/sh->fHist->Integral());
                else WriteWarningStatus("TRExFit::CorrectHistograms","Sample Hist of sample "+smp->fNormToSample+" not found ...");
            }

            //
            // For SampleUp / SampleDown
            for(auto syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                // get the original syst histograms & reset the syst histograms
                SystematicHist *syh = sh->GetSystematic( syst->fName );
                //
                // if syst defined with SampleUp / SampleDown
                if( syst->fSampleUp != "" || syst->fSampleDown != "" ){
                    bool isDummy = ( syst->fDummyForSamples.size()>0 && FindInStringVector(syst->fDummyForSamples,smp->fName)>=0 );
                    TH1 *h_up   = sh->fHist;
                    if(syst->fSampleUp   !="" && !isDummy){
                        if(reg->GetSampleHist(syst->fSampleUp  )){
                            h_up   = reg->GetSampleHist(syst->fSampleUp  )->fHist;
                        }
                    }
                    TH1 *h_down = sh->fHist;
                    if(syst->fSampleDown !="" && !isDummy){
                        if(reg->GetSampleHist(syst->fSampleDown)){
                            h_down = reg->GetSampleHist(syst->fSampleDown)->fHist;
                        }
                    }
                    syh = sh->AddHistoSyst(syst->fName,h_up,h_down);
                    syh->fSystematic = syst;
                }
            }

            //
            // Save to _preSmooth histograms (to be shown in syst plots) at this point
            sh->fHist_preSmooth = (TH1*)sh->fHist->Clone(Form("%s_preSmooth",sh->fHist->GetName()));
            for(auto syh : sh->fSyst){
                if(syh!=nullptr){
                    if(syh->fHistUp!=nullptr)   syh->fHistUp_preSmooth   = (TH1*)syh->fHistUp->Clone(  Form("%s_preSmooth",syh->fHistUp->GetName()  ));
                    else                        syh->fHistUp_preSmooth   = (TH1*)sh->fHist_preSmooth->Clone();
                    if(syh->fHistDown!=nullptr) syh->fHistDown_preSmooth = (TH1*)syh->fHistDown->Clone(Form("%s_preSmooth",syh->fHistDown->GetName()));
                    else                        syh->fHistDown_preSmooth = (TH1*)sh->fHist_preSmooth->Clone();
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
                    SmoothHistogramTtres( sh->fHist );
                } else {
                    isFlat = SmoothHistogram( sh->fHist );
                }
                h_correction->Divide( h0 );
            }

            //
            // Systematics
            for(auto syst : smp->fSystematics){
                if(syst==nullptr) continue;
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                SystematicHist *syh = sh->GetSystematic( syst->fName );
                if(syh==nullptr) continue;
                TH1* hUp   = syh->fHistUp;
                TH1* hDown = syh->fHistDown;
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
                    if(hUp!=nullptr  ) SmoothHistogram( hUp  , isFlat );
                    if(hDown!=nullptr) SmoothHistogram( hDown, isFlat );
                }

                //
                // Histogram smoothing, Symmetrisation, Massaging...
                if(!reg->fSkipSmoothing) syh -> fSmoothType = syst -> fSmoothType;
                else                                syh -> fSmoothType = 0;
                syh -> fSymmetrisationType = syst -> fSymmetrisationType;

            }  // end syst loop
            //
            // Histograms checking
            for(auto syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                if( sh->GetSystematic( syst->fName )==nullptr ) continue;
                //
                HistoTools::CheckHistograms( reg->GetSampleHist(smp->fName)->fHist /*nominal*/,
                                            sh->GetSystematic( syst->fName ) /*systematic*/,
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
        for(auto sh : reg->fSampleHists){
            if(sh->fHist==nullptr) continue;
            for(auto syst : fSystematics){
                if(  FindInStringVector(syst->fDropNormIn, reg->fName)>=0
                  || FindInStringVector(syst->fDropNormIn, sh->fSample->fName)>=0
                  || FindInStringVector(syst->fDropNormIn, "all")>=0
                  ){
                    SystematicHist* syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    if(sh->fHist->Integral()!=0){
                        WriteDebugStatus("TRExFit::CorrectHistograms", "  Normalising syst " + syst->fName + " for sample " + sh->fSample->fName);
                        if(syh->fHistUp  !=nullptr) syh->fHistUp  ->Scale(sh->fHist->Integral()/syh->fHistUp  ->Integral());
                        if(syh->fHistDown!=nullptr) syh->fHistDown->Scale(sh->fHist->Integral()/syh->fHistDown->Integral());
                    }
                }
            }
        }
    }

    // drop shape part of systematic according to fDropShapeIn
    for(auto reg : fRegions){
        for(auto sh : reg->fSampleHists){
            if(sh->fHist==nullptr) continue;
            for(auto syst : fSystematics){
                if(  FindInStringVector(syst->fDropShapeIn, reg->fName)>=0
                  || FindInStringVector(syst->fDropShapeIn, sh->fSample->fName)>=0
                  || FindInStringVector(syst->fDropShapeIn, "all")>=0
                  ){
                    SystematicHist* syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    WriteDebugStatus("TRExFit::CorrectHistograms", "  Removing shape component of syst " + syst->fName + " for sample " + sh->fSample->fName);
                    if(syh->fHistUp != nullptr) {
                        const double ratioUp = syh->fHistUp->Integral()/sh->fHist->Integral();
                        delete syh->fHistUp;
                        syh->fHistUp = static_cast<TH1*>(sh->fHist->Clone());
                        syh->fHistUp->Scale(ratioUp);
                    }
                    if(syh->fHistDown != nullptr) {
                        const double ratioDown = syh->fHistDown->Integral()/sh->fHist->Integral();
                        delete syh->fHistDown;
                        syh->fHistDown = static_cast<TH1*>(sh->fHist->Clone());
                        syh->fHistDown->Scale(ratioDown);
                    }
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
                    if(FindInStringVector(subSamples,smp->fName)<0) continue;
                    SampleHist *sh = reg->GetSampleHist(smp->fName);
                    if(sh==nullptr) continue;
                    SystematicHist *syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    yieldNominal += sh ->fHist    ->Integral();
                    yieldUp      += syh->fHistUp  ->Integral();
                    yieldDown    += syh->fHistDown->Integral();
                }
                // scale each syst variation
                for(auto smp : fSamples){
                    if(FindInStringVector(subSamples,smp->fName)<0) continue;
                    SampleHist *sh = reg->GetSampleHist(smp->fName);
                    if(sh==nullptr) continue;
                    SystematicHist *syh = sh->GetSystematic(syst->fName);
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
                SampleHist *shNominal = nullptr;
                for(auto sh : reg->fSampleHists){
                    if(!sh->fSample->fIsMorph[par]) continue;
                    if(sh->fSample->fMorphValue[par]==nominalValue){ // FIXME: eventually add something to flag a sample as nominal for morphing
                        shNominal = sh;
                        break;
                    }
                }
                // loop on all other samples
                for(auto sh : reg->fSampleHists){
                    if(!sh->fSample->fIsMorph[par]) continue;
                    if(sh!=shNominal){
                        for(auto syh : shNominal->fSyst){
                            Systematic * syst = syh->fSystematic;
                            if(syst->fIsNormOnly){
                                SystematicHist* syhNew = sh->AddOverallSyst(syst->fName,syst->fOverallUp,syst->fOverallDown);
                                syhNew->fSystematic = syst;
                            }
                            else{
                                TH1* hUpNew   = (TH1*)syh->fHistUp->Clone();
                                TH1* hDownNew = (TH1*)syh->fHistDown->Clone();
                                hUpNew->Divide(shNominal->fHist);
                                hUpNew->Multiply(sh->fHist);
                                hDownNew->Divide(shNominal->fHist);
                                hDownNew->Multiply(sh->fHist);
                                SystematicHist* syhNew = sh->AddHistoSyst(syst->fName,hUpNew,hDownNew);
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
                if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                SampleHist *sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                sh->fSample->fUseSystematics = true;
                //
                SampleHist *shReference = reg->GetSampleHist(smp->fSystFromSample);
                for(auto syh : shReference->fSyst){
                    Systematic * syst = syh->fSystematic;
                    if(syst->fIsNormOnly){
                        SystematicHist* syhNew = sh->AddOverallSyst(syst->fName,syst->fOverallUp,syst->fOverallDown);
                        syhNew->fSystematic = syst;
                    }
                    else{
                        TH1* hUpNew   = (TH1*)syh->fHistUp->Clone();
                        TH1* hDownNew = (TH1*)syh->fHistDown->Clone();
                        hUpNew->Divide(shReference->fHist);
                        hUpNew->Multiply(sh->fHist);
                        hDownNew->Divide(shReference->fHist);
                        hDownNew->Multiply(sh->fHist);
                        SystematicHist* syhNew = sh->AddHistoSyst(syst->fName,hUpNew,hDownNew);
                        syhNew->fSystematic = syst;
                    }
                }
            }
        }
    }

    //
    // set the hasData flag
    bool hasData = false;
    for(auto smp : fSamples){
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
                    TH1 *hdata = reg->fData->fHist;
                    if(TRExFitter::OPTION["PoissonizeData"]>0) gRandom->SetSeed(TRExFitter::OPTION["PoissonizeData"]);
                    for(int i_bin=1;i_bin<=hdata->GetNbinsX();i_bin++){
                        hdata->SetBinContent(i_bin,gRandom->Poisson( hdata->GetBinContent(i_bin) ));
                        hdata->SetBinError(i_bin,sqrt(hdata->GetBinContent(i_bin)));
                    }
                }
            }
        }
    }

    //
    // Drop bins
    for(auto reg : fRegions){
        if(reg->fDropBins.size()!=0){
            for(auto smp : fSamples){
                // eventually skip sample / region combination
                if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                SampleHist *sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                DropBins(sh->fHist,reg->fDropBins);
                for(auto syst : smp->fSystematics){
                    // eventually skip systematic / region combination
                    if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                    if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                    if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                    SystematicHist *syh = sh->GetSystematic( syst->fName );
                    if(syh==nullptr) continue;
                    DropBins(syh->fHistUp,  reg->fDropBins);
                    DropBins(syh->fHistDown,reg->fDropBins);
                    if(syh->fHistShapeUp)   DropBins(syh->fHistShapeUp,  reg->fDropBins);
                    if(syh->fHistShapeDown) DropBins(syh->fHistShapeDown,reg->fDropBins);
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::ReadHistograms(){
    TH1D* hUp = nullptr;
    TH1D* hDown = nullptr;

    //
    // Loop on regions and samples
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        WriteInfoStatus("TRExFit::ReadHistograms", "  Region " + fRegions[i_ch]->fName + " ...");
        //
        if(TRExFitter::SPLITHISTOFILES) fFiles[i_ch]->cd();
        //
        if(fRegions[i_ch]->fBinTransfo != "") ComputeBinning(i_ch);
        // first we must read the DATA samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType!=Sample::DATA) continue;
            WriteDebugStatus("TRExFit::ReadHistograms", "  Reading DATA sample " + fSamples[i_smp]->fName);
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            std::vector<std::string> fullPaths = FullHistogramPaths(fRegions[i_ch],fSamples[i_smp]);

            TH1D* h = ReadSingleHistogram(fullPaths, nullptr, i_ch, i_smp, true, false); // is nominal and not MC
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );
            //
            // Importing the histogram in TRExFitter
            SampleHist *sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
            sh->fHist_orig = h_orig;
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            // in fact DATA can be used for systs that have SubtractRefSampleVar: TRUE
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                Systematic *syst = fSamples[i_smp]->fSystematics[i_syst];
                // only relevant for systs that have this sample as reference
                if (!syst->fSubtractRefSampleVar || syst->fReferenceSample != fSamples[i_smp]->fName) continue;

                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadHistograms", "Adding syst " + syst->fName);
                //

                //
                // Up
                //
                hUp = nullptr;
                if(syst->fHasUpVariation){
                    fullPaths = FullHistogramPaths(fRegions[i_ch],fSamples[i_smp],syst,true);
                    hUp = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, true, false); // is up variation and not MC
                }
                //
                // Down
                //
                hDown = nullptr;
                if(syst->fHasDownVariation){
                    fullPaths = FullHistogramPaths(fRegions[i_ch],fSamples[i_smp],syst,false);
                    hDown = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, false, false); // is down variation and not MC
                }
                //
                if(hUp==nullptr){
                    hUp   = (TH1D*)fRegions[i_ch]->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                }
                if(hDown==nullptr){
                    hDown = (TH1D*)fRegions[i_ch]->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                }
                //
                SystematicHist *syh = sh->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                syh->fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                syh->fScaleUp = fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0){
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[fRegions[i_ch]->fName]!=0){
                        syh->fScaleUp *= fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[fRegions[i_ch]->fName];
                    }
                }
                syh->fScaleDown = fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0){
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[fRegions[i_ch]->fName]!=0){
                        syh->fScaleDown *= fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[fRegions[i_ch]->fName];
                    }
                }
            }
        }

        // then we can read the other samples
        std::set < std::string > files_names;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            WriteDebugStatus("TRExFit::ReadHistograms", "  Reading " + fSamples[i_smp]->fName);
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            std::vector<std::string>fullPaths = FullHistogramPaths(fRegions[i_ch],fSamples[i_smp]);
            for (const auto& ipath : fullPaths){
                files_names.insert(ipath);
            }
            TH1D* h = ReadSingleHistogram(fullPaths, nullptr, i_ch, i_smp, true, true); // is MC
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );
            //
            // Importing the histogram in TRExFitter
            SampleHist *sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
            sh->fHist_orig = h_orig;
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            // end here no systematics allowed (e.g. generally for GHOST samples)
            if (!fSamples[i_smp]->fUseSystematics) continue;

            //
            //  -----------------------------------
            //
            // read norm factors
            for(int i_norm=0;i_norm<fSamples[i_smp]->fNNorm;i_norm++){
                NormFactor *nf = fSamples[i_smp]->fNormFactors[i_norm];
                //
                // eventually skip systematic / region combination
                if( nf->fRegions.size()>0 && FindInStringVector(nf->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( nf->fExclude.size()>0 && FindInStringVector(nf->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadHistograms", "Adding norm " + nf->fName);
                //
                sh->AddNormFactor( nf );
            }

            //
            //  -----------------------------------
            //
            // read shape factors
            for(int i_shape=0;i_shape<fSamples[i_smp]->fNShape;i_shape++){
                ShapeFactor *sf = fSamples[i_smp]->fShapeFactors[i_shape];
                //
                // eventually skip systematic / region combination
                if( sf->fRegions.size()>0 && FindInStringVector(sf->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( sf->fExclude.size()>0 && FindInStringVector(sf->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadHistograms", "Adding shape " + sf->fName);
                //
                sh->AddShapeFactor( sf );
            }

            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                Systematic *syst = fSamples[i_smp]->fSystematics[i_syst];
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                WriteDebugStatus("TRExFit::ReadHistograms", "Adding syst " + syst->fName);
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL){
                    SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddOverallSyst(syst->fName,syst->fOverallUp,syst->fOverallDown);
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
                // else ...
                //
                if(syst->fReferenceSample!="") smp = GetSample(syst->fReferenceSample);
                //
                // Up
                //
                hUp = nullptr;
                if(syst->fHasUpVariation){
                    fullPaths     = FullHistogramPaths(fRegions[i_ch],fSamples[i_smp],syst,true);
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
                    fullPaths     = FullHistogramPaths(fRegions[i_ch],fSamples[i_smp],syst,false);
                    for (const auto& ipath : fullPaths){
                        files_names.insert(ipath);
                    }
                    hDown = ReadSingleHistogram(fullPaths, syst, i_ch, i_smp, false, true); // isUp and isMC
                }
                //
                if(hUp==nullptr)   hUp   = (TH1D*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                if(hDown==nullptr) hDown = (TH1D*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                //
                SystematicHist *syh = sh->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                syh->fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                syh->fScaleUp = fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0)
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName]!=0)
                        syh->fScaleUp *= fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[reg->fName];
                syh->fScaleDown = fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0)
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName]!=0)
                        syh->fScaleDown *= fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[reg->fName];
            }
            //closing the files for this sample
            CloseFiles( files_names );
            files_names.clear();
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::ReadHistos(/*string fileName*/){
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
        if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos.root";
        else                 fileName = fName + "/Histograms/" + fInputName + "_histos.root";
        // Bootstrap
        if(fBootstrap!="" && fBootstrapIdx>=0){
            fileName = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fBootstrapIdx));
        }
        WriteInfoStatus("TRExFit::ReadHistos", "-----------------------------");
        WriteInfoStatus("TRExFit::ReadHistos", "Reading histograms from file " + fileName + " ...");
    }
    //
    std::vector< TH2F* > histPrun;
    std::unique_ptr<TFile> filePrun = nullptr;
    if( fKeepPruning ){
        filePrun = std::unique_ptr<TFile>(new TFile( (fName+"/Pruning.root").c_str() ));
        if(!filePrun) fKeepPruning = false;
    }
    //
    // when we multply/divide by or subtract/add other samples, need to add systematics on the other samples
    Systematic * tmpsyst;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(!fSamples[i_smp]->fUseSystematics) continue;
        if(fSamples[i_smp]->fDivideBy!=""){
            Sample* smp = GetSample(fSamples[i_smp]->fDivideBy);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("TRExFit::ReadHistos", " The sample " + fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("TRExFit::ReadHistos", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst])[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fNSyst;j_syst++){
                       if(fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if( FindInStringVector(fSystematics[j_syst]->fSamples,fSamples[i_smp]->fName)<0 ) fSystematics[j_syst]->fSamples.push_back(fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
        if(fSamples[i_smp]->fMultiplyBy!=""){
            Sample* smp = GetSample(fSamples[i_smp]->fMultiplyBy);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("TRExFit::ReadHistos", " The sample " + fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("TRExFit::ReadHistos", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst])[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fNSyst;j_syst++){
                       if(fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if( FindInStringVector(fSystematics[j_syst]->fSamples,fSamples[i_smp]->fName)<0 ) fSystematics[j_syst]->fSamples.push_back(fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
        for(auto sample : fSamples[i_smp]->fSubtractSamples){
            Sample* smp = GetSample(sample);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("TRExFit::ReadHistos", " The sample " + fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("TRExFit::ReadHistos", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst])[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fNSyst;j_syst++){
                       if(fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if( FindInStringVector(fSystematics[j_syst]->fSamples,fSamples[i_smp]->fName)<0 ) fSystematics[j_syst]->fSamples.push_back(fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
        for(auto sample : fSamples[i_smp]->fAddSamples){
            Sample* smp = GetSample(sample);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systNPName = smp->fSystematics[i_syst]->fNuisanceParameter;
                if(!fSamples[i_smp]->HasNuisanceParameter(systNPName)){
                    WriteDebugStatus("TRExFit::ReadHistos", " The sample " + fSamples[i_smp]->fName + " doesn't have natively NP "+ systNPName);
                    WriteDebugStatus("TRExFit::ReadHistos", "                Inheriting it from "+smp->fName);
                    tmpsyst = new Systematic((smp->fSystematics[i_syst])[0]);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    fSamples[i_smp]->AddSystematic(tmpsyst);
                    for(int j_syst=0;j_syst<fNSyst;j_syst++){
                       if(fSystematics[j_syst]->fName==tmpsyst->fName) {
                          if( FindInStringVector(fSystematics[j_syst]->fSamples,fSamples[i_smp]->fName)<0 ) fSystematics[j_syst]->fSamples.push_back(fSamples[i_smp]->fName);
                       }
                    }
                }
            }
        }
    }
    //
    // Syst for morphing samples inherited from nominal sample
    if (fPropagateSystsForMorphing){
        for(auto par : fMorphParams){
            double nominalValue = 0.;
            for(auto norm : fNormFactors){
                if(norm->fName==par) nominalValue = norm->fNominal;
            }
            Sample *smpNominal = nullptr;
            for(auto smp : fSamples){
                if(!smp->fIsMorph[par]) continue;
                if(smp->fMorphValue[par]==nominalValue){ // FIXME: eventually add something to flag a sample as nominal for morphing
                    smpNominal = smp;
                    break;
                }
            }
            for(auto smp : fSamples){
                if(!smp->fIsMorph[par]) continue;
                if(smp==smpNominal) continue;
                for(auto syst : smpNominal->fSystematics){
                    smp->AddSystematic(syst);
                }
            }
        }
    }
    //
    // fSystFromSample
    for(auto smp : fSamples){
        if(smp->fSystFromSample!=""){
            Sample *smpReference = nullptr;
            for(auto smp2 : fSamples){
                if(smp2->fName==smp->fName) continue;
                if(smp2->fName==smp->fSystFromSample){
                    smpReference = smp2;
                    break;
                }
            }
            if(smpReference!=nullptr){
                for(auto syst : smpReference->fSystematics){
                    smp->AddSystematic(syst);
                }
            }
        }
    }
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        if( fKeepPruning ){
            histPrun.push_back( (TH2F*)filePrun->Get( Form("h_prun_%s_toSave", fRegions[i_ch]->fName.c_str()) ) );
        }
        regionName = fRegions[i_ch]->fName;
        WriteDebugStatus("TRExFit::ReadHistos","  Reading region " + regionName);
        //
        if(!singleOutputFile){
            if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + regionName + "_histos.root";
            else                 fileName = fName + "/Histograms/" + fInputName + "_" + regionName + "_histos.root";
            // Bootstrap
            if(fBootstrap!="" && fBootstrapIdx>=0){
                if( fBootstrapSyst == "" )
                        fileName = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fBootstrapIdx));
                else
                        fileNameBootstrap = ReplaceString(fileName,"_histos.root",Form("_histos__%d.root",fBootstrapIdx));
            }
            WriteInfoStatus("TRExFit::ReadHistos", "-----------------------------");
            WriteInfoStatus("TRExFit::ReadHistos", "Reading histograms from file " + fileName + " ...");
        }
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,regionName)<0 && fSamples[i_smp]->fName.find("customAsimov_")==std::string::npos ) continue;
            //
            sampleName = fSamples[i_smp]->fName;
            WriteDebugStatus("TRExFit::ReadHistos", "    Reading sample " + sampleName);
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp],regionName+"_"+sampleName,fileName);
            sh = fRegions[i_ch]->GetSampleHist(sampleName);
            if(sh==nullptr) continue;
            //
            // separate gammas -> Add systematic
            if(fSamples[i_smp]->fSeparateGammas){
                std::string systName = "stat_"+fSamples[i_smp]->fName;
                std::string systStoredName = systName;
                WriteDebugStatus("TRExFit::ReadHistos", "adding separate gammas as SHAPE systematic " + systName);
                SystematicHist *syh_tmp = sh->AddHistoSyst(systName,
                                                        Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                        Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                        0
                                                        );
                if(syh_tmp==nullptr){
                    WriteWarningStatus("TRExFit::ReadHistos", "No histogram found for separate gamma, but may be you will create it right now.");
                }
                else{
                    Systematic *gamma = nullptr;
                    if(FindInStringVector(fSystematicNames,systName)>=0) gamma = fSystematics[FindInStringVector(fSystematicNames,systName)];  //GetSystematic(systName);
                    if(gamma==nullptr) gamma = NewSystematic(systName);
                    WriteDebugStatus("TRExFit::ReadHistos", "adding separate gammas as SHAPE systematic " + systName);
                    gamma->fType = Systematic::SHAPE;
                    gamma->fRegions.clear();
                    gamma->fRegions.push_back(fRegions[i_ch]->fName);
                    syh_tmp->fSystematic = gamma;
                    gamma->fNuisanceParameter = gamma->fName;
                    TRExFitter::NPMAP[gamma->fName] = gamma->fNuisanceParameter;
                }
            }
            //
            // norm factors
            for(int i_norm=0;i_norm<fSamples[i_smp]->fNNorm;i_norm++){
                //
                // eventually skip norm factor / region combination
                if( fSamples[i_smp]->fNormFactors[i_norm]->fRegions.size()>0 && FindInStringVector(fSamples[i_smp]->fNormFactors[i_norm]->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( fSamples[i_smp]->fNormFactors[i_norm]->fExclude.size()>0 && FindInStringVector(fSamples[i_smp]->fNormFactors[i_norm]->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                //
                normName = fSamples[i_smp]->fNormFactors[i_norm]->fName;
                WriteDebugStatus("TRExFit::ReadHistos", "      Reading norm " + normName);
                // norm only
                sh->AddNormFactor(fSamples[i_smp]->fNormFactors[i_norm]);
            }
            //
            // shape factors
            for(int i_shape=0;i_shape<fSamples[i_smp]->fNShape;i_shape++){
                //
                // eventually skip shape factor / region combination
                if( fSamples[i_smp]->fShapeFactors[i_shape]->fRegions.size()>0 && FindInStringVector(fSamples[i_smp]->fShapeFactors[i_shape]->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( fSamples[i_smp]->fShapeFactors[i_shape]->fExclude.size()>0 && FindInStringVector(fSamples[i_smp]->fShapeFactors[i_shape]->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                //
                shapeName = fSamples[i_smp]->fShapeFactors[i_shape]->fName;
                WriteDebugStatus("TRExFit::ReadHistos", "      Reading shape " + shapeName);
                // shape only
                sh->AddShapeFactor(fSamples[i_smp]->fShapeFactors[i_shape]);
            }
            //
            // systematics
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                //
                // eventually skip systematic / region combination
                if( fSamples[i_smp]->fSystematics[i_syst]->fRegions.size()>0 && FindInStringVector(fSamples[i_smp]->fSystematics[i_syst]->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( fSamples[i_smp]->fSystematics[i_syst]->fExclude.size()>0 && FindInStringVector(fSamples[i_smp]->fSystematics[i_syst]->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                if( fSamples[i_smp]->fSystematics[i_syst]->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(fSamples[i_smp]->fSystematics[i_syst]->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                std::string systName       = fSamples[i_smp]->fSystematics[i_syst]->fName;
                std::string systStoredName = fSamples[i_smp]->fSystematics[i_syst]->fStoredName; // if no StoredName specified in the config, this should be == fName
                //
                // eventually skip systematics if pruned
                int xbin,ybin,bin;
                int binContent = 0;
                if( fKeepPruning && histPrun[i_ch]!=nullptr ){
                    xbin = histPrun[i_ch]->GetXaxis()->FindBin( sampleName.c_str() ); // sample
                    ybin = histPrun[i_ch]->GetYaxis()->FindBin( systName.c_str() ); // syst
                    bin = histPrun[i_ch]->GetBin(xbin,ybin);
                    binContent = histPrun[i_ch]->GetBinContent(bin);
                    if( binContent <= -4 || binContent == -1 || binContent >= 3 ){
                        WriteDebugStatus("TRExFit::ReadHistos", "SKIPPING systematic " + systName);
                        continue;
                    }
                    //{kBlack,6,kBlue, kGray, 8, kYellow, kOrange-3, kRed}
                }
                WriteDebugStatus("TRExFit::ReadHistos", "      Reading syst " + systName);
                // norm only
                if(fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    if( fKeepPruning ){
                        if( binContent == -2 || binContent == 2 ) continue;
                    }
                    syh = sh->AddOverallSyst(systName,
                         fSamples[i_smp]->fSystematics[i_syst]->fOverallUp,
                         fSamples[i_smp]->fSystematics[i_syst]->fOverallDown);
                }
                // histo syst
                else{
                    int pruned = 0;
                    if( fKeepPruning ){
                        if(binContent==1 || binContent==-2) pruned = 1;
                        if(binContent==2 || binContent==-3) pruned = 2;
                    }
                    if(fBootstrap!="" && fBootstrapIdx>=0 && fBootstrapSyst == systName ){
                        syh = sh->AddHistoSyst(systName,
                                          Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileNameBootstrap,
                                          Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileNameBootstrap,
                                          pruned
                                          );
                    }
                    else{
                        syh = sh->AddHistoSyst(systName,
                                          Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                          Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                          pruned
                                          );
                    }
                    if(syh==nullptr){
                        if (!pruned) WriteWarningStatus("TRExFit::ReadHistos", "No syst histo found for syst " + systName + ", sample " + sampleName + ", region " + regionName);
                        continue;
                    }
                }
                // for both
                syh->fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                syh->fHistoNameShapeUp   = Form("%s_%s_%s_Shape_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                syh->fHistoNameShapeDown = Form("%s_%s_%s_Shape_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                if(fBootstrap!="" && fBootstrapIdx>=0 && fBootstrapSyst == systName ){
                    syh->fFileNameShapeUp    = fileNameBootstrap;
                    syh->fFileNameShapeDown  = fileNameBootstrap;
                }
                else{
                    syh->fFileNameShapeUp    = fileName;
                    syh->fFileNameShapeDown  = fileName;
                }
                syh->fScaleUp = fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0){
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName]!=0){
                        syh->fScaleUp *= fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName];
                    }
                }
                syh->fScaleDown = fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0){
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName]!=0){
                        syh->fScaleDown *= fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName];
                    }
                }
                //
                if(fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
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

//__________________________________________________________________________________
//
void TRExFit::CloseInputFiles(){
    //
    // Close all input files
    for(auto it : TRExFitter::TFILEMAP){
        TDirectory *dir = gDirectory;
        TFile *f = it.second;
        if(f!=nullptr)
        dir->cd();
        f->Close();
        delete f;
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
        for(auto reg : fRegions){
            TH1* hTot = reg->GetTotHist(true);
            double totPred = hTot->Integral();
            double totData = reg->fData->fHist->Integral();
            double totToScale = 0;
            std::vector<SampleHist*> shToScale;
            for(auto sh : reg->fSampleHists){
                if(sh->fHist==nullptr) continue;
                if(sh->fSample->fType==Sample::GHOST){
                    WriteWarningStatus("TRExFit::CorrectHistograms","Requested to scale to data a GHOST sample, " + sh->fSample->fName + ". Skipping this sample.");
                    continue;
                }
                if(FindInStringVector(fScaleSamplesToData,sh->fSample->fName)>=0){
                    shToScale.emplace_back(sh);
                    double morph_scale = GetNominalMorphScale(sh);
                    totToScale += morph_scale*sh->fHist->Integral();
                }
            }
            if(totToScale<=0 || shToScale.size()==0) continue;
            double scale = (totData-(totPred-totToScale))/totToScale;
            for(auto sh : shToScale){
                WriteInfoStatus("TRExFit::CorrectHistograms","Scaling sample " + sh->fSample->fName + " by " + std::to_string(scale) + " in region " + reg->fName);
                sh->Scale(scale);
            }

        }
    }
    else if(fScaleSamplesToData.size()>0){
        for(auto reg : fRegions){
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
        TRExPlot *p = nullptr;
        fRegions[i_ch]->fUseStatErr = fUseStatErr;
        fRegions[i_ch]->fATLASlabel = fAtlasLabel;
        //
        if(fCustomAsimov!=""){
            std::string name = "customAsimov_"+fCustomAsimov;
            SampleHist* cash = fRegions[i_ch]->GetSampleHist(name);
            if(cash==nullptr){
                WriteWarningStatus("TRExFit::DrawAndSaveAll", "No Custom Asimov " + fCustomAsimov + " available. Taking regular Asimov.");
            }
            else{
                std::string s = cash->fHist->GetName();
                WriteDebugStatus("TRExFit::DrawAndSaveAll", "  Adding Custom-Asimov Data: " + s);
                fRegions[i_ch]->fData = cash;
            }
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
            for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
                p->SaveAs(     (fName+"/Plots/"+fRegions[i_ch]->fName+"_postFit"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format] ).c_str());

            if(fWithPullTables){
                pullTex << "\\hline\\hline\n" << std::endl;
                pullTex << "\\end{tabular}"  << std::endl;
                pullTex << "\\end{document}" << std::endl;
                pullTex.close();
            }
            delete p;
        }
        else{
            if(fRegions[i_ch]->fRegionDataType==Region::ASIMOVDATA) p = fRegions[i_ch]->DrawPreFit(fPrePostFitCanvasSize, opt+" blind");
            else                                                    p = fRegions[i_ch]->DrawPreFit(fPrePostFitCanvasSize, opt);
            // this line to fix the y-axis maximum getting doubled in some cases (FIXME)
            if((fRegions[i_ch]->fYmin==0) && (fRegions[i_ch]->fYmax==0) && (fRegions[i_ch]->fYmaxScale==0)){
                if(!fRegions[i_ch]->fLogScale) p->h_dummy->GetYaxis()->SetRangeUser(p->h_dummy->GetYaxis()->GetXmin(),p->h_dummy->GetMaximum());
                else                           p->h_dummy->GetYaxis()->SetRangeUser(1                                ,p->h_dummy->GetMaximum());
            }
            for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
                p->SaveAs(     (fName+"/Plots/"+fRegions[i_ch]->fName+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format] ).c_str());
            }
            if( !fKeepPrefitBlindedBins ) delete p;
        }
    }
}

//__________________________________________________________________________________
//
TRExPlot* TRExFit::DrawSummary(std::string opt, TRExPlot* prefit_plot) {
    WriteInfoStatus("TRExFit::DrawSummary", "-------------------------------------------");
    WriteInfoStatus("TRExFit::DrawSummary", "Building Summary Plot...");
    gSystem->mkdir(fName.c_str(),true);
    const bool isPostFit = opt.find("post")!=std::string::npos;
    const bool checkVR = opt.find("valid")!=std::string::npos;
    if(TRExFitter::POISSONIZE) opt += " poissonize";
    // build one bin per region
    TH1D* h_data = 0;
    TH1D* h_sig[MAXsamples];
    TH1D* h_bkg[MAXsamples];
    TH1D *h_tot;
    TGraphAsymmErrors *g_err;
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
    TH1* h; // to store varius histograms temporary
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
        SampleHist *sh = nullptr;
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
            h_sig[Nsig] = new TH1D(name.c_str(),title.c_str(), Nbin,0,Nbin);
            std::string temp_string = h_sig[Nsig]->GetTitle();
            WriteDebugStatus("TRExFit::DrawSummary", "Adding Signal: " + temp_string);
            h_sig[Nsig]->SetLineColor(lineColor);
            h_sig[Nsig]->SetFillColor(fillColor);
            h_sig[Nsig]->SetLineWidth(lineWidth);
            for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=nullptr){
                    if(isPostFit)  h = (TH1D*)sh->fHist_postFit->Clone(); // Michele
                    else           h = (TH1D*)sh->fHist->Clone(); // Michele
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
                h_sig[Nsig]->SetBinContent( i_bin,integral );
                h_sig[Nsig]->SetBinError( i_bin,intErr );
            }
            Nsig++;
        }
        else if(fSamples[i_smp]->fType==Sample::BACKGROUND){
            h_bkg[Nbkg] = new TH1D(name.c_str(),title.c_str(), Nbin,0,Nbin);
            std::string temp_string = h_bkg[Nbkg]->GetTitle();
            WriteDebugStatus("TRExFit::DrawSummary", "Adding Bkg:    " + temp_string);
            h_bkg[Nbkg]->SetLineColor(lineColor);
            h_bkg[Nbkg]->SetFillColor(fillColor);
            h_bkg[Nbkg]->SetLineWidth(lineWidth);
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=nullptr){
                    if(isPostFit)  h = (TH1D*)sh->fHist_postFit->Clone(); // Michele
                    else           h = (TH1D*)sh->fHist->Clone(); // Michele
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
                    //
                    // this becuase MC stat is taken into account by the gammas
                    if( (isPostFit && fUseGammaPulls) || !fUseStatErr || (!sh->fSample->fUseMCStat && !sh->fSample->fSeparateGammas))
                        intErr = 0.;
                }
                else{
                    integral = 0.;
                    intErr = 0.;
                }
                h_bkg[Nbkg]->SetBinContent( i_bin,integral );
                h_bkg[Nbkg]->SetBinError( i_bin,intErr );
            }
            Nbkg++;
        }
        else if(fSamples[i_smp]->fType==Sample::DATA){
            h_data = new TH1D(name.c_str(),title.c_str(), Nbin,0,Nbin);
            std::string temp_string = h_data->GetTitle();
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
    TRExPlot *p;
    //
    if (fSummaryCanvasSize.size() == 0){
        p = new TRExPlot(fInputName+"_summary",900,700,TRExFitter::NORATIO);
    } else {
        p = new TRExPlot(fInputName+"_summary",fSummaryCanvasSize.at(0),fSummaryCanvasSize.at(1),TRExFitter::NORATIO);
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
    if(fBlindingThreshold>=0){
        p->SetBinBlinding(true,fBlindingThreshold,fBlindingType);
        if(isPostFit && fKeepPrefitBlindedBins && fBlindedBins) p->SetBinBlinding(true,fBlindedBins,fBlindingType);
    }
    //
    if(h_data) p->SetData(h_data, h_data->GetTitle());
    for(int i=0;i<Nsig;i++){
        if(TRExFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    h_sig[i],h_sig[i]->GetTitle());
        if(TRExFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(h_sig[i],h_sig[i]->GetTitle());
        if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(h_sig[i],h_sig[i]->GetTitle());
    }
    for(int i=0;i<Nbkg;i++){
        p->AddBackground(h_bkg[i],h_bkg[i]->GetTitle());
    }

    if( TRExFitter::PREFITONPOSTFIT && isPostFit) {
      p->h_tot_bkg_prefit = (TH1*)prefit_plot->GetTotBkg()->Clone("h_tot_bkg_prefit");
    }

    //
    // Build tot
    //
    h_tot = new TH1D("h_Tot_summary","h_Tot_summary", Nbin,0,Nbin);

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
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    TH1* h_tmp_Up;
    TH1* h_tmp_Down;
    std::vector<std::string> systNames;
    std::vector<std::string> npNames;
    int i_np = -1;
    // actual systematics
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(isPostFit && fSystematics[i_syst]->fType == Systematic::SHAPE) continue;
        std::string systName = fSystematics[i_syst]->fName;
        std::string systNuisPar = systName;
        systNames.push_back( systName );
        if(fSystematics[i_syst]!=nullptr)
            systNuisPar = fSystematics[i_syst]->fNuisanceParameter;
        if(FindInStringVector(npNames,systNuisPar)<0){
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
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                }
                else{
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp_postFit[syst_idx];
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown_postFit[syst_idx];
                }
            }
            else{
                if(syst_idx<0){
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTot;
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTot;
                }
                else{
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp[syst_idx];
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown[syst_idx];
                }
            }
            if(i_bin==1){
                h_up.  push_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,systName.c_str()), Form("h_Tot_%s_Up_TMP",  systName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1D(Form("h_Tot_%s_Down_TMP",systName.c_str()), Form("h_Tot_%s_Down_TMP",systName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
            //
            // look for other syst with the same np
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(j_syst==syst_idx) continue;
                if(systNuisPar==TRExFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
                    TH1* h_tmp = nullptr;
                    if(isPostFit){
                        h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp_postFit[j_syst];
                        h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown_postFit[j_syst];
                        h_tmp      = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                    }
                    else{
                        h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp[j_syst];
                        h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown[j_syst];
                        h_tmp      = fRegions[regionVec[i_bin-1]]->fTot;
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
                    h_tmp_Up   = region->fTot_postFit;
                    h_tmp_Down = region->fTot_postFit;
                }
                else{
                    h_tmp_Up   = region->fTotUp_postFit[syst_idx];
                    h_tmp_Down = region->fTotDown_postFit[syst_idx];
                }
                h_up.  push_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,gammaName.c_str()), Form("h_Tot_%s_Up_TMP",  gammaName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1D(Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Nbin,0,Nbin) );
                h_up[i_np]  ->SetBinContent( i_ch,h_tmp_Up  ->Integral() );
                h_down[i_np]->SetBinContent( i_ch,h_tmp_Down->Integral() );
            }
        }
    }
    // add the sample-specific gammas (only if post-fit)
    if(isPostFit && fUseGammaPulls){
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
                        h_tmp_Up   = region->fTot_postFit;
                        h_tmp_Down = region->fTot_postFit;
                    }
                    else{
                        h_tmp_Up   = region->fTotUp_postFit[syst_idx];
                        h_tmp_Down = region->fTotDown_postFit[syst_idx];
                    }
                    h_up.  push_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,gammaName.c_str()), Form("h_Tot_%s_Up_TMP",  gammaName.c_str()), Nbin,0,Nbin) );
                    h_down.push_back( new TH1D(Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Form("h_Tot_%s_Down_TMP",gammaName.c_str()), Nbin,0,Nbin) );
                    h_up[i_np]  ->SetBinContent( i_ch,h_tmp_Up  ->Integral() );
                    h_down[i_np]->SetBinContent( i_ch,h_tmp_Down->Integral() );
                }
            }
        }
    }
    // add the norm factors
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        std::string normName = fNormFactors[i_norm]->fName;
        if(FindInStringVector(npNames,normName)<0){
            npNames.push_back(normName);
            i_np++;
        }
        else
            continue;
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
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                }
                else{
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp_postFit[syst_idx];
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown_postFit[syst_idx];
                }
            }
            else{
                if(syst_idx<0){
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTot;
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTot;
                }
                else{
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp[syst_idx];
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown[syst_idx];
                }
            }
            if(i_bin==1){
                h_up.  push_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,normName.c_str()), Form("h_Tot_%s_Up_TMP",  normName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1D(Form("h_Tot_%s_Down_TMP",normName.c_str()), Form("h_Tot_%s_Down_TMP",normName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
        }
    }
    //
    if(isPostFit)  g_err = BuildTotError( h_tot, h_up, h_down, npNames, fFitResults->fCorrMatrix );
    else           g_err = BuildTotError( h_tot, h_up, h_down, npNames );
    //
    p->SetTotBkg(h_tot);
    p->BlindData();
    p->SetTotBkgAsym(g_err);
    //
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        p->SetBinLabel(i_bin,fRegions[regionVec[i_bin-1]]->fShortLabel.c_str());
    }
    p->Draw(opt);
    if(!isPostFit) fBlindedBins = p->h_blinding;
    //
    if(divisionVec.size()>0){
        p->pad0->cd();
        TLine line(0,1,0,1);
        line.SetNDC(0);
        line.SetLineStyle(7);
        line.SetLineColor(kBlack);
        line.SetLineWidth(2);
        line.DrawLine(divisionVec[0],((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMinimum(),
                      divisionVec[0],pow(((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(),0.73) );
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
                if(divisionVec.size()<ii) break;
                double xmax = Nbin;
                double xmin = 0.;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                double xpos = xmin + 0.5*(xmax - xmin);
                double ypos = pow(((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(), 0.61 );
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
                if(divisionVec.size()<ii) break;
                double xmax = Nbin;
                double xmin = 0.;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                double xpos = xmin + 0.5*(xmax - xmin);
                double ypos = pow(((TH1D*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(), 0.61 );
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
    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
        if(fSummaryPrefix!=""){
            if(isPostFit)  p->SaveAs((fName+"/Plots/"+fSummaryPrefix+"_Summary_postFit"+(checkVR?"_VR":"")+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
            else           p->SaveAs((fName+"/Plots/"+fSummaryPrefix+"_Summary"        +(checkVR?"_VR":"")+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
        }
        else{
            if(isPostFit)  p->SaveAs((fName+"/Plots/Summary_postFit"+(checkVR?"_VR":"")+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
            else           p->SaveAs((fName+"/Plots/Summary"        +(checkVR?"_VR":"")+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
        }
    }
    //
    for(int i_syst=0;i_syst<(int)h_up.size();i_syst++){
        delete h_up[i_syst];
        delete h_down[i_syst];
    }
    h_up.clear();
    h_down.clear();
    //
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
    double ymaxTmp = -1.;
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
        ymaxTmp = 0;
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
            for(auto sampleHist : region->fSampleHists){
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
                            ScaleNominal(sampleHist, h_tmp);
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
            hDataVec[i_channel]->SetBinError(i_bin,sqrt(hDataVec[i_channel]->GetBinContent(i_bin)));
        }
        hDataVec[i_channel]->Scale( scale );
        for(auto hVec : hSignalVec)     hVec[i_channel]->Scale( scale );
        for(auto hVec : hBackgroundVec) hVec[i_channel]->Scale( scale );
    }
    //
    // merge and plot them
    TRExPlot *p;
    int cWidthMerge = 1200;
    int cHeightMerge = 600;
    if(fMergeCanvasSize.size()>1){
        cWidthMerge = fMergeCanvasSize.at(0);
        cHeightMerge = fMergeCanvasSize.at(1);
    }
    p = new TRExPlot(fInputName+"_merge",cWidthMerge,cHeightMerge,TRExFitter::NORATIO);
    //
    p->SetData(MergeHistograms(hDataVec),"");
    for(unsigned int i_sig=0;i_sig<hSignalVec.size();i_sig++){
        if(TRExFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    MergeHistograms(hSignalVec[i_sig]),"");
        if(TRExFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(MergeHistograms(hSignalVec[i_sig]),"");
        if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(MergeHistograms(hSignalVec[i_sig]),"");

    }
    for(unsigned int i_bkg=0;i_bkg<hBackgroundVec.size();i_bkg++) p->AddBackground(MergeHistograms(hBackgroundVec[i_bkg]),"");
    p->SetTotBkg(MergeHistograms(hTotVec));
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
    if(fBlindingThreshold>=0){
        p->SetBinBlinding(true,fBlindingThreshold,fBlindingType);
//         if(isPostFit && fKeepPrefitBlindedBins && fBlindedBins) p->SetBinBlinding(true,fBlindedBins,fBlindingType); // FIXME
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
    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
        p->SaveAs((saveName+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
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
    TH1D* h_smp[MAXsamples];
    TH1D *h_tot;
    TGraphAsymmErrors *g_err[MAXsamples];
    TGraphAsymmErrors *g_err_tot;
    //
    std::string name;
    std::string title;
    //
    double intErr; // to store the integral error
    TH1* h0; // to store varius histograms temporary
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
    SampleHist *sh = nullptr;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        name = fSamples[i_smp]->fName;
        title = fSamples[i_smp]->fTitle;
        //
        int idx = FindInStringVector(titleVec,title);
        if(idx>=0){
            idxVec.push_back(idx);
        }
        else{
            idxVec.push_back(i_smp);
            h_smp[idxVec[i_smp]] = new TH1D(("h_"+name).c_str(),title.c_str(), Nbin,0,Nbin);
        }
        for(unsigned int i_bin=1;i_bin<=regionVec.size();i_bin++){
            sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
            if(sh!=nullptr){
                if(isPostFit && fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::GHOST)
                    h0 = sh->fHist_postFit;
                else
                    h0 = sh->fHist;
                double tmpErr = h_smp[idxVec[i_smp]]->GetBinError(i_bin); // Michele -> get the error before adding content to bin, to avoid ROOT automatically increasing it!
                double scale = 1.;
                if (!isPostFit){
                    scale = GetNominalMorphScale(sh);
                }
                h_smp[idxVec[i_smp]]->AddBinContent( i_bin,scale*h0->IntegralAndError(1,h0->GetNbinsX(),intErr) );
                intErr*=scale;
                if( (isPostFit && fUseGammaPulls) || !fUseStatErr || (!sh->fSample->fUseMCStat && !sh->fSample->fSeparateGammas))
                    h_smp[idxVec[i_smp]]->SetBinError(i_bin,0.);
                else
                    h_smp[idxVec[i_smp]]->SetBinError(i_bin, sqrt( pow(tmpErr,2) + pow(intErr,2) ) );
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
        std::vector< TH1* > h_up;
        std::vector< TH1* > h_down;
        TH1* h_tmp_Up;
        TH1* h_tmp_Down;
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
                            h_tmp_Up   = sh->fHist_postFit;
                            h_tmp_Down = sh->fHist_postFit;
                        }
                        else{
                            h_tmp_Up   = sh->GetSystematic(systName)->fHistUp_postFit;
                            h_tmp_Down = sh->GetSystematic(systName)->fHistDown_postFit;
                        }
                    }
                    else {
                        if(syst_idx<0 || sh->GetSystematic(systName)==nullptr){
                            h_tmp_Up   = sh->fHist;
                            h_tmp_Down = sh->fHist;
                        }
                        else{
                            h_tmp_Up   = sh->GetSystematic(systName)->fHistUp;
                            h_tmp_Down = sh->GetSystematic(systName)->fHistDown;
                        }
                    }
                }
                else {
                    h_tmp_Up   = new TH1D(Form("h_DUMMY_%s_up_%i",  systName.c_str(),i_bin-1),"h_dummy",1,0,1);
                    h_tmp_Down = new TH1D(Form("h_DUMMY_%s_down_%i",systName.c_str(),i_bin-1),"h_dummy",1,0,1);
                }
                if(i_bin==1){
                    h_up.  push_back( new TH1D(Form("h_%s_%s_Up_TMP",  name.c_str(),systName.c_str()),Form("h_%s_%s_Up_TMP",  name.c_str(),systName.c_str()), Nbin,0,Nbin) );
                    h_down.push_back( new TH1D(Form("h_%s_%s_Down_TMP",name.c_str(),systName.c_str()),Form("h_%s_%s_Down_TMP",name.c_str(),systName.c_str()), Nbin,0,Nbin) );
                }
                double scale = 1.;
                if (!isPostFit){
                    scale = GetNominalMorphScale(sh);
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
                                h_tmp_Up   = sh->fHist_postFit;
                                h_tmp_Down = sh->fHist_postFit;
                            }
                            else{
                                h_tmp_Up   = sh->GetSystematic(systName)->fHistUp_postFit;
                                h_tmp_Down = sh->GetSystematic(systName)->fHistDown_postFit;
                            }
                        }
                        else{
                            if(syst_idx<0 || sh->GetSystematic(systName)==nullptr){
                                h_tmp_Up   = sh->fHist;
                                h_tmp_Down = sh->fHist;
                            }
                            else{
                                h_tmp_Up   = sh->GetSystematic(systName)->fHistUp;
                                h_tmp_Down = sh->GetSystematic(systName)->fHistDown;
                            }
                        }
                        double morph_scale = 1.;
                        if (!isPostFit) {
                            morph_scale = GetNominalMorphScale(sh);
                        }
                        h_up[i_np]  ->AddBinContent( i_bin,(h_tmp_Up  ->Integral(1,h_tmp_Up->GetNbinsX()))*morph_scale );
                        h_down[i_np]->AddBinContent( i_bin,(h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()))*morph_scale );
                    }
                }
            }
        }
        //
        if(isPostFit)  g_err[i_smp] = BuildTotError( h_smp[i_smp], h_up, h_down, npNames, fFitResults->fCorrMatrix );
        else           g_err[i_smp] = BuildTotError( h_smp[i_smp], h_up, h_down, npNames );
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
        if(fSamples[i_smp]->fType==Sample::DATA) texout << "\\hline " << std::endl;
        if(fSamples[i_smp]->fTexTitle!="") texout << "  " << fSamples[i_smp]->fTexTitle << "  ";
        else                               texout << "  " << fSamples[i_smp]->fTitle << "  ";
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
            double mean = h_smp[i_smp]->GetBinContent(i_bin);
            double uncertainty = ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
            double mean_rounded = mean;
            double uncertainty_rounded = uncertainty;
            int n = -1; // this will contain the number of decimal places
            if (fUseATLASRoundingTxt || fUseATLASRoundingTex){
                n = ApplyATLASrounding(mean_rounded, uncertainty_rounded);
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
    }
    //
    // Build tot
    //
    h_tot = new TH1D("h_Tot_","h_Tot", Nbin,0,Nbin);
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot_postFit->IntegralAndError(1,fRegions[regionVec[i_bin-1]]->fTot_postFit->GetNbinsX(),intErr) );
        else          h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot->IntegralAndError(        1,fRegions[regionVec[i_bin-1]]->fTot->GetNbinsX(),        intErr) );
        h_tot->SetBinError( i_bin, intErr );
    }
    //
    //   Build error band
    // build the vectors of variations
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    TH1* h_tmp_Up;
    TH1* h_tmp_Down;
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
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                }
                else{
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp_postFit[syst_idx];
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown_postFit[syst_idx];
                }
            }
            else{
                if(syst_idx<0){
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTot;
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTot;
                }
                else{
                    h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp[syst_idx];
                    h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown[syst_idx];
                }
            }
            if(i_bin==1){
                h_up.  push_back( new TH1D(Form("h_Tot_%s_Up_TMP"  ,systName.c_str()), Form("h_Tot_%s_Up_TMP",  systName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1D(Form("h_Tot_%s_Down_TMP",systName.c_str()), Form("h_Tot_%s_Down_TMP",systName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
            //
            // look for other syst with the same np
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(j_syst==syst_idx) continue;
                if(systNuisPar==TRExFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
                    TH1* h_tmp = nullptr;
                    if(isPostFit){
                        h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp_postFit[j_syst];
                        h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown_postFit[j_syst];
                        h_tmp      = fRegions[regionVec[i_bin-1]]->fTot_postFit;
                    }
                    else{
                        h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp[j_syst];
                        h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown[j_syst];
                        h_tmp      = fRegions[regionVec[i_bin-1]]->fTot;
                    }
                    h_up[i_np]  ->AddBinContent( i_bin,h_tmp_Up  ->Integral()-h_tmp->Integral() );
                    h_down[i_np]->AddBinContent( i_bin,h_tmp_Down->Integral()-h_tmp->Integral() );
                }
            }
        }
    }
    //
    if(isPostFit)  g_err_tot = BuildTotError( h_tot, h_up, h_down, npNames, fFitResults->fCorrMatrix );
    else           g_err_tot = BuildTotError( h_tot, h_up, h_down, npNames );
    //
    if(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) out << " | Total | ";
    else                                                    out << " | Tot.Bkg. | ";
    texout << "\\hline " << std::endl;
    if(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) texout << "  Total ";
    else                                                    texout << "  Total background ";
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        double mean = h_tot->GetBinContent(i_bin);
        double uncertainty = ( g_err_tot->GetErrorYhigh(i_bin-1) + g_err_tot->GetErrorYlow(i_bin-1) )/2.;
        double mean_rounded = mean;
        double uncertainty_rounded = uncertainty;
        int n = -1; // this will contain the number of decimal places
        if (fUseATLASRoundingTxt || fUseATLASRoundingTex){
            n = ApplyATLASrounding(mean_rounded, uncertainty_rounded);
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

    //
    // Print data
    if( !fFitIsBlind ){
        texout << "\\hline " << std::endl;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if( fSamples[i_smp]->fType!=Sample::DATA  ) continue;
            if(idxVec[i_smp]!=i_smp) continue;
            //
            // print values
            out << " | " << fSamples[i_smp]->fTitle << " | ";
            if(fSamples[i_smp]->fTexTitle!="") texout << "  " << fSamples[i_smp]->fTexTitle << "  ";
            else                               texout << "  " << fSamples[i_smp]->fTitle << "  ";
            for(int i_bin=1;i_bin<=Nbin;i_bin++){
                texout << " & ";
                out << h_smp[i_smp]->GetBinContent(i_bin);
                texout << Form("%.0f",h_smp[i_smp]->GetBinContent(i_bin));
                out << " | ";
            }
            out << std::endl;
            texout << " \\\\ ";
            texout << std::endl;
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
    //
    for(int i_syst=0;i_syst<(int)h_up.size();i_syst++){
        delete h_up[i_syst];
        delete h_down[i_syst];
    }
    h_up.clear();
    h_down.clear();
    //
    if(fCleanTables){
        std::string shellcommand = "cat "+fName+"/Tables/Yields"+suffix+".tex|sed -e \"s/\\#/ /g\" > "+fName+"/Tables/Yields";
        if(isPostFit) shellcommand += "_postFit";
        shellcommand += suffix+"_clean.tex";
        gSystem->Exec(shellcommand.c_str());
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
    if(TRExFitter::OPTION["FourTopStyle"]!=0) H0 = 75; // height of the top label pad
    if(TRExFitter::OPTION["FourTopStyle"]!=0) Hp = 200;
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
                const double scale = GetNominalMorphScale(regions[i]->fSig[i_sig]);
                S[i] += scale * regions[i]->fSig[i_sig]->fHist->Integral();
            }
        }
        for(int i_bkg=0;i_bkg<regions[i]->fNBkg;i_bkg++){
            if(regions[i]->fBkg[i_bkg]!=nullptr) {
                const double scale = GetNominalMorphScale(regions[i]->fBkg[i_bkg]);
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
        h.back().SetBinContent(2,S[i]/sqrt(B[i]));
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
        yMax = TMath::Max(yMax,h.back().GetMaximum());
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
    for(std::size_t i_format=0;i_format<TRExFitter::IMAGEFORMAT.size();i_format++) {
        c.SaveAs((fName+"/SignalRegions"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
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
    for(std::size_t i_format=0;i_format<TRExFitter::IMAGEFORMAT.size(); ++i_format){
        c.SaveAs((fName+"/PieChart" + fSuffix + ( isPostFit ? "_postFit" : "" ) + "."+TRExFitter::IMAGEFORMAT[i_format]).c_str());
    }
}

//__________________________________________________________________________________
// called before w in case of CustomAsimov
void TRExFit::CreateCustomAsimov() const{
    WriteDebugStatus("TRExFit::CreateCustomAsimov", "Running CreateCustomAsimov");
    // get a list of all CustomAsimov to create
    std::vector<std::string> customAsimovList;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fAsimovReplacementFor.first!="" && FindInStringVector(customAsimovList,fSamples[i_smp]->fAsimovReplacementFor.first)<0)
            customAsimovList.push_back(fSamples[i_smp]->fAsimovReplacementFor.first);
    }
    //
    // fill a different CustomAsimov data-set for each element in the list
    for(auto customAsimov : customAsimovList){
        WriteDebugStatus("TRExFit::CreateCustomAsimov", "CustomAsimov: " + customAsimov);
        Sample *ca = GetSample("customAsimov_"+customAsimov);
        // create a new data sample taking the nominal S and B
        for(int i_ch=0;i_ch<fNRegions;i_ch++){
            Region *reg = fRegions[i_ch];
            SampleHist *cash = reg->SetSampleHist(ca,(TH1*)reg->fData->fHist->Clone());
            cash->fHist_orig->SetName( Form("%s_orig",cash->fHist->GetName()) ); // fix the name
            cash->fHist->Scale(0.);
            //
            std::vector<std::string> smpToExclude;
            for(int i_smp=0;i_smp<fNSamples;i_smp++){
                SampleHist* h = reg->GetSampleHist(fSamples[i_smp]->fName);
                if( h==0 ) continue;
                if( h->fSample->fType==Sample::DATA ) continue;
                if( h->fSample->fType==Sample::GHOST ){
                    if( h->fSample->fAsimovReplacementFor.first!=customAsimov ) continue;
                    if( h->fSample->fAsimovReplacementFor.second!="" ) smpToExclude.push_back(h->fSample->fAsimovReplacementFor.second);
                }
                if( FindInStringVector( smpToExclude,fSamples[i_smp]->fName )>=0 ) continue;
                //
                // bug-fix: change normalisation factors to nominal value!
                double factor = 1.;
                for(auto norm : fSamples[i_smp]->fNormFactors){
                    WriteDebugStatus("TRExFit::CreateCustomAsimov", "setting norm factor to " + std::to_string(norm->fNominal));
                    factor *= norm->fNominal;
                }
                //
                cash->fHist->Add(h->fHist,factor);
            }
            cash->fHist->Sumw2(false);
        }
    }
}

//__________________________________________________________________________________
// turn to RooStat::HistFactory
void TRExFit::ToRooStat(bool makeWorkspace, bool exportOnly){

    WriteInfoStatus("TRExFit::ToRooStat", "-------------------------------------------");
    WriteInfoStatus("TRExFit::ToRooStat", "Exporting to RooStats...");

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);

    //Suffix used for the regular bin transformed histogram
    const std::string suffix_regularBinning = "_regBin";

    RooStats::HistFactory::Measurement meas((fInputName+fSuffix).c_str(), (fInputName+fSuffix).c_str());
    if(fBootstrap!="" && fBootstrapIdx>=0)
        meas.SetOutputFilePrefix((fName+"/RooStats/"+fBootstrapSyst+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName).c_str());
    else
        meas.SetOutputFilePrefix((fName+"/RooStats/"+fInputName).c_str());
    meas.SetExportOnly(exportOnly);
    meas.SetPOI(fPOI.c_str());
    meas.SetLumi(fLumiScale);
    if(fLumiErr==0){
        meas.AddConstantParam("Lumi");
        meas.SetLumiRelErr(0.1);
    } else {
        meas.SetLumiRelErr(fLumiErr);
    }

    for(int i_ch=0;i_ch<fNRegions;i_ch++){

        if(fRegions[i_ch]->fRegionType==Region::VALIDATION) continue;

        WriteDebugStatus("TRExFit::ToRooStat", "Adding Channel: " + fRegions[i_ch]->fName);
        RooStats::HistFactory::Channel chan(fRegions[i_ch]->fName.c_str());

        //Checks if a data sample exists
        bool hasData = false;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType==Sample::DATA){
                hasData = true;
                break;
            }
        }
        if(fCustomAsimov!=""){
            std::string name = "customAsimov_"+fCustomAsimov;
            SampleHist* cash = fRegions[i_ch]->GetSampleHist(name);
            if(cash==nullptr){
                if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
                WriteWarningStatus("TRExFit::ToRooStat", "No Custom Asimov " + fCustomAsimov + " available. Taking regular Asimov.");
                if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
            }
            else{
                std::string temp_string = cash->fHist->GetName();
                WriteDebugStatus("TRExFit::ToRooStat", "  Adding Custom-Asimov Data: " + temp_string);
                chan.SetData(cash->fHistoName+suffix_regularBinning, cash->fFileName);
            }
        }
        else if(hasData){
            std::string temp_string = fRegions[i_ch]->fData->fHist->GetName();
            WriteDebugStatus("TRExFit::ToRooStat", "  Adding Data: " + temp_string);
            chan.SetData(fRegions[i_ch]->fData->fHistoName+suffix_regularBinning, fRegions[i_ch]->fData->fFileName);
        } else {
            chan.SetData("", "");
        }

        // fStatErrCons is upper case after config reading if the MCstatThreshold option is used, otherwise it defaults to "Poisson"
        // HistFactory expects the constraint not in all uppercase, but in form "Poisson"/"Gaussian" instead
        if(fStatErrCons=="Poisson" || fStatErrCons=="POISSON") chan.SetStatErrorConfig(fStatErrThres, "Poisson");
        else if(fStatErrCons=="GAUSSIAN")                      chan.SetStatErrorConfig(fStatErrThres, "Gaussian");

        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            SampleHist* h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if( h != nullptr && h->fSample->fType!=Sample::DATA && h->fSample->fType!=Sample::GHOST ){
                WriteDebugStatus("TRExFit::ToRooStat", "  Adding Sample: " + fSamples[i_smp]->fName);
                RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());
                if(fUseStatErr && fSamples[i_smp]->fUseMCStat) sample.ActivateStatError();
                sample.SetHistoName(h->fHistoName+suffix_regularBinning);
                sample.SetInputFile(h->fFileName);
                sample.SetNormalizeByTheory(fSamples[i_smp]->fNormalizedByTheory);
                // norm factors
                for(int i_norm=0;i_norm<h->fNNorm;i_norm++){
                    WriteDebugStatus("TRExFit::ToRooStat", "    Adding NormFactor: " + h->fNormFactors[i_norm]->fName + ", " + std::to_string(h->fNormFactors[i_norm]->fNominal));
                    sample.AddNormFactor( h->fNormFactors[i_norm]->fName,
                                         h->fNormFactors[i_norm]->fNominal,
                                         h->fNormFactors[i_norm]->fMin,
                                         h->fNormFactors[i_norm]->fMax);
                    if (h->fNormFactors[i_norm]->fConst) meas.AddConstantParam( h->fNormFactors[i_norm]->fName );
                    if (fStatOnly && fFixNPforStatOnlyFit && h->fNormFactors[i_norm]->fName!=fPOI)
                        meas.AddConstantParam( h->fNormFactors[i_norm]->fName );
                }
                // shape factors
                for(int i_shape=0;i_shape<h->fNShape;i_shape++){
                    WriteDebugStatus("TRExFit::ToRooStat", "    Adding ShapeFactor: " + h->fShapeFactors[i_shape]->fName + ", " + std::to_string(h->fShapeFactors[i_shape]->fNominal));
                    sample.AddShapeFactor( h->fShapeFactors[i_shape]->fName );
                    if (h->fShapeFactors[i_shape]->fConst
                        || (fStatOnly && fFixNPforStatOnlyFit && h->fShapeFactors[i_shape]->fName!=fPOI)
                    ){
                        for(int i_bin=0;i_bin<h->fShapeFactors[i_shape]->fNbins;i_bin++){
                            meas.AddConstantParam( "gamma_" + h->fShapeFactors[i_shape]->fName + "_bin_" + std::to_string(i_bin) );
                        }
                    }
                }
                // systematics
                if(!fStatOnly){
                    for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
                        // add normalization part
                        WriteDebugStatus("TRExFit::ToRooStat", "    Adding Systematic: " + h->fSyst[i_syst]->fName);
                        if ( h->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
                            std::string npName = "shape_";
                            npName += h->fSyst[i_syst]->fSystematic->fNuisanceParameter+"_";
                            std::string regionName = fRegions[i_ch]->fName;
                            if(h->fSyst[i_syst]->fSystematic->fNuisanceParameter.find("stat_")!=std::string::npos){
                                // see if there are regions to correlate with others
                                for(auto set : h->fSample->fCorrelateGammasInRegions){
                                    for(unsigned int i_reg=0;i_reg<set.size();i_reg++){
                                        if(i_reg!=0 && regionName==set[i_reg]){
                                            regionName = set[0];
                                            break;
                                        }
                                    }
                                }
                                // eventually correlate MC stat with other samples
                                if(h->fSample->fCorrelateGammasWithSample!=""){
                                    npName = "shape_stat_"+h->fSample->fCorrelateGammasWithSample+"_";
                                }
                            }
                            npName += regionName;
                            sample.AddShapeSys( npName, RooStats::HistFactory::Constraint::Poisson,
                                                (h->fSyst[i_syst]->fHistoNameUp+"_Var"+suffix_regularBinning),
                                                h->fSyst[i_syst]->fFileNameUp, "");
                        }
                        else{
                            if ( !h->fSyst[i_syst]->fSystematic->fIsShapeOnly  &&
                                !h->fSyst[i_syst]->fNormPruned  &&
                                !h->fSyst[i_syst]->fBadNorm
                              ) {
                                if(h->fSyst[i_syst]->fSystematic->fNuisanceParameter=="ttXsec"){
                                    WriteDebugStatus("TRExFit::ToRooStat", "Syst norm up: " + std::to_string(h->fSyst[i_syst]->fNormUp));
                                }
                                sample.AddOverallSys( h->fSyst[i_syst]->fSystematic->fNuisanceParameter,
                                                      1+h->fSyst[i_syst]->fNormDown,
                                                      1+h->fSyst[i_syst]->fNormUp   );
                            }
                            // eventually add shape part
                            if ( h->fSyst[i_syst]->fIsShape  &&
                                !h->fSyst[i_syst]->fSystematic->fIsNormOnly  &&
                                !h->fSyst[i_syst]->fShapePruned  &&
                                !h->fSyst[i_syst]->fBadShape
                              ){
                                sample.AddHistoSys( h->fSyst[i_syst]->fSystematic->fNuisanceParameter,
                                                    h->fSyst[i_syst]->fHistoNameShapeDown+suffix_regularBinning, h->fSyst[i_syst]->fFileNameShapeDown, "",
                                                    h->fSyst[i_syst]->fHistoNameShapeUp+suffix_regularBinning,   h->fSyst[i_syst]->fFileNameShapeUp,   ""  );
                            }
                        }
                    }
                }
                else{
                    sample.AddOverallSys( "Dummy",1,1 );
                }
                chan.AddSample(sample);
            }
        }
        meas.AddChannel(chan);
    }
    // Experimental: turn off constraints for given systematics
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSystematics[i_syst]->fIsFreeParameter) meas.AddUniformSyst(fSystematics[i_syst]->fName.c_str());
    }
    // morphing
    for(const TRExFit::TemplateWeight& itemp : fTemplateWeightVec){
        std::string normName = "morph_"+itemp.name+"_"+ReplaceString(std::to_string(itemp.value),"-","m");
        WriteDebugStatus("TRExFit::ToRooStat", "Morphing: normName: " + normName);
        meas.AddPreprocessFunction(normName, itemp.function, itemp.range);
    }
    for(auto nf : fNormFactors){
        if(nf->fExpression.first!=""){
            meas.AddPreprocessFunction(nf->fName,nf->fExpression.first,nf->fExpression.second);
        }
    }
    //
    if(fBootstrap!="" && fBootstrapIdx>=0)
        meas.PrintXML((fName+"/RooStats/"+fBootstrapSyst+"_BSId"+Form("%d",fBootstrapIdx)+"/").c_str());
    else
        meas.PrintXML((fName+"/RooStats/").c_str());
    CloseInputFiles();
    meas.CollectHistograms();
    meas.PrintTree();

    if(makeWorkspace) RooStats::HistFactory::MakeModelAndMeasurementFast(meas);

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
}

//__________________________________________________________________________________
//
void TRExFit::SystPruning() const{
    WriteInfoStatus("TRExFit::SystPruning", "------------------------------------------------------");
    WriteInfoStatus("TRExFit::SystPruning", "Apply Systematics Pruning ...");
    if(fSystematics.size()==0 || fStatOnly){
        WriteInfoStatus("TRExFit::SystPruning", "No systematics => No Pruning applied.");
        return;
    }

    PruningUtil *pu = new PruningUtil();
    pu->SetStrategy((int)fPruningType);
    pu->SetThresholdNorm(fThresholdSystPruning_Normalisation);
    pu->SetThresholdShape(fThresholdSystPruning_Shape);
    pu->SetThresholdIsLarge(fThresholdSystLarge);

    for(auto reg : fRegions){
        // if want to skip validation regions from pruning, add a condition here
        reg->SystPruning(pu);
    }

    // it also writes the txt file actually
    DrawPruningPlot();
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
    std::vector< Sample* > samplesVec;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        samplesVec.push_back(fSamples[i_smp]);
        nSmp++;
    }
    // make a list of non-gamma systematics only
    std::vector<Systematic*> nonGammaSystematics;
    std::vector<std::string> uniqueSyst;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSystematics[i_syst]->fType!=Systematic::SHAPE){
            nonGammaSystematics.push_back(fSystematics[i_syst]);

            // fill names of unique systs
            if (std::find(uniqueSyst.begin(), uniqueSyst.end(), fSystematics[i_syst]->fName) == uniqueSyst.end()){
                uniqueSyst.emplace_back(fSystematics[i_syst]->fName);
            }
        }
    }
    const size_t NnonGammaSyst = nonGammaSystematics.size();
    if(NnonGammaSyst==0){
        WriteInfoStatus("TRExFit::DrawPruningPlot", "No non-gamma systematics found => No Pruning plot generated.");
        return;
    }
    //
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        if(fRegions[i_reg]->fRegionType==Region::VALIDATION) continue;

        out << "In Region : " << fRegions[i_reg]->fName << std::endl ;
        histPrun.emplace_back( std::move(std::unique_ptr<TH2F>(new TH2F (Form("h_prun_%s", fRegions[i_reg]->fName.c_str()  ),fRegions[i_reg]->fShortLabel.c_str(),nSmp,0,nSmp, uniqueSyst.size(),0,uniqueSyst.size()))));
        histPrun.back()->SetDirectory(0);

        for(int i_smp=0;i_smp<nSmp;i_smp++){
            out << " -> In Sample : " << samplesVec[i_smp]->fName << std::endl;

            for(std::size_t uniqueIndex = 0; uniqueIndex < uniqueSyst.size(); ++uniqueIndex){
               histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -1 );
            }

            SampleHist *sh = fRegions[i_reg]->GetSampleHist(samplesVec[i_smp]->fName);
            if (sh == nullptr) continue;

            for(size_t i_syst=0;i_syst<NnonGammaSyst;i_syst++){
                // find the corresponding index of unique syst
                auto it = std::find(uniqueSyst.begin(), uniqueSyst.end(), nonGammaSystematics.at(i_syst)->fName);
                const std::size_t uniqueIndex = std::distance(uniqueSyst.begin(), it);
                out << " --->>  " << nonGammaSystematics[i_syst]->fName << "     " ;
                if( (FindInStringVector(nonGammaSystematics[i_syst]->fSamples,samplesVec[i_smp]->fName)>=0 || nonGammaSystematics[i_syst]->fSamples[0] == "all")
                    && sh->HasSyst(nonGammaSystematics[i_syst]->fName)
                ){
                    SystematicHist *syh = sh->GetSystematic(nonGammaSystematics[i_syst]->fName);
                    histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 0 );
                    //
                    if(syh->fShapePruned && syh->fNormPruned) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 3 );
                    else if(syh->fShapePruned) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 1 );
                    else if(syh->fNormPruned) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 2 );
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
    int regionSize = 20*fNSamples;
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
    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++) {
        c.SaveAs( (fName+"/Pruning"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }

    //
    // Save prunign hist for future usage
    std::unique_ptr<TFile> filePrun = nullptr;
    // - checking if Pruning.root exists
    // if yes
    if(!gSystem->AccessPathName( (fName+"/Pruning.root").c_str() )){
        // ...
        filePrun = std::unique_ptr<TFile>( new TFile( (fName+"/Pruning.root").c_str() ));
    }
    else{
        filePrun = std::unique_ptr<TFile> (new TFile( (fName+"/Pruning.root").c_str(),"RECREATE" ));
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

    //
    //Checks if a data sample exists
    //
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    RooDataSet* data = nullptr;
    RooWorkspace* ws = nullptr;

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
        TFile *rootFile = new TFile(fWorkspaceFileName.c_str(),"read");
        ws = (RooWorkspace*) rootFile->Get("combined");
        if(!ws){
            WriteErrorStatus("TRExFit::Fit", "The workspace (\"combined\") cannot be found in file " + fWorkspaceFileName + ". Please check !");
            exit(EXIT_FAILURE);
        }
        if(!fFitIsBlind && hasData) data = (RooDataSet*)ws->data("obsData");
        else                        data = (RooDataSet*)ws->data("asimovData");
    }
    //
    // Otherwise go on with normal fit
    //
    else{
        if (!isLHscanOnly){
            WriteInfoStatus("TRExFit::Fit","");
            WriteInfoStatus("TRExFit::Fit","-------------------------------------------");
            WriteInfoStatus("TRExFit::Fit","Performing nominal fit...");
        }
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std:: string > regionsToFit;
        std::map < std::string, int > regionDataType;
        for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
            bool isToFit = false;

            if ( fFitRegion == CRONLY ) {
                if( fRegions[i_ch] -> fRegionType == Region::CONTROL ){
                    isToFit = true;
                }
            } else if ( fFitRegion == CRSR ){
                if( fRegions[i_ch] -> fRegionType == Region::CONTROL || fRegions[i_ch] -> fRegionType == Region::SIGNAL ){
                    isToFit = true;
                }
            }
            if ( ! isToFit ){
                for (unsigned int iReg = 0; iReg < fFitRegionsToFit.size(); ++iReg ){
                    if( fFitRegionsToFit[iReg] == fRegions[i_ch] -> fName ){
                        isToFit = true;
                        break;
                    }
                }
            }
            //
            if(isToFit){
                regionsToFit.push_back( fRegions[i_ch] -> fName );
                Region::DataType dataType;
                if(fFitIsBlind || !hasData){
                    dataType = Region::ASIMOVDATA;
                }
                else{
                    dataType = fRegions[i_ch] -> fRegionDataType;
                }
                regionDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , dataType) );
            }
        }
        //
        // Creating the combined model with the regions to fit only
        //
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        ws = PerformWorkspaceCombination( regionsToFit );
        if (!ws){
            WriteErrorStatus("TRExFit::Fit","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        //
        // If needed (only if needed), create a RooDataset object
        //
        data = DumpData( ws, regionDataType, fFitNPValues, fFitPOIAsimov );
        //
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    }

    // Calls the PerformFit() function to actually do the fit
    //
    if (!isLHscanOnly) PerformFit( ws, data, fFitType, true, TRExFitter::DEBUGLEVEL);

    //
    // Toys
    //
    if(fFitToys>0 && !isLHscanOnly){
        RunToys(ws);
    }

    //
    // Fit result on Asimov with shifted systematics
    //
    if(fDoNonProfileFit && !isLHscanOnly){
        //
        std::map<std::string,double> systGroups;
        std::vector<std::string> systGroupNames;
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
            out.open((fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts.txt").c_str());
            tex.open((fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts.tex").c_str());
            tex2.open((fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts_grouped.tex").c_str());
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
        std::map < std::string, double > npValues;
        // nominal fit on Asimov
        RooStats::ModelConfig *mc = (RooStats::ModelConfig*)ws -> obj("ModelConfig");
        ws->saveSnapshot("InitialStateModelGlob",   *mc->GetGlobalObservables());
        ws->saveSnapshot("InitialStateModelNuis",   *mc->GetNuisanceParameters());
        WriteInfoStatus("TRExFit::Fit","Fitting nominal Asimov...");
        data = DumpData( ws, regionDataType, fFitNPValues, fFitPOIAsimov );
        npValues = PerformFit( ws, data, fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
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
                for(auto sh : reg->fSampleHists){
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
                    if(hTot==nullptr && sh->fHist!=nullptr) hTot = (TH1*)sh->fHist->Clone("h_tot");
                    else if(            sh->fHist!=nullptr) hTot->Add(sh->fHist);
                }
                if(hTot==nullptr) continue;
                for(int i_bin=1;i_bin<=hTot->GetNbinsX();i_bin++){
                    double statErr = hTot->GetBinError(i_bin)/hTot->GetBinContent(i_bin);
                    std::string gammaName = "gamma_stat_"+reg->fName+"_bin_"+std::to_string(i_bin-1);
                    // up
                    npVal = fFitNPValues;
                    npVal[gammaName] = 1+statErr;
                    WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1+statErr));
                    data = DumpData( ws, regionDataType, npVal, fFitPOIAsimov );
                    npVal[gammaName] = 1;
                    ws->loadSnapshot("InitialStateModelGlob");
                    ws->loadSnapshot("InitialStateModelNuis");
                    npValues = PerformFit( ws, data, fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                    MCstatUp = sqrt(pow(MCstatUp,2)+pow(npValues[fPOI]-nominalPOIval,2));
                    // down
                    npVal = fFitNPValues;
                    npVal[gammaName] = 1-statErr;
                    WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1-statErr));
                    data = DumpData( ws, regionDataType, npVal, fFitPOIAsimov );
                    npVal[gammaName] = 1;
                    ws->loadSnapshot("InitialStateModelGlob");
                    ws->loadSnapshot("InitialStateModelNuis");
                    npValues = PerformFit( ws, data, fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                    MCstatDo = -sqrt(pow(MCstatDo,2)+pow(npValues[fPOI]-nominalPOIval,2));
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
                for(auto sh : reg->fSampleHists){
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
                        data = DumpData( ws, regionDataType, npVal, fFitPOIAsimov );
                        npVal[gammaName] = 1;
                        ws->loadSnapshot("InitialStateModelGlob");
                        ws->loadSnapshot("InitialStateModelNuis");
                        npValues = PerformFit( ws, data, fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                        MCstatUpSample[sh->fSample->fName] = sqrt(pow(MCstatUpSample[sh->fSample->fName],2) + pow(npValues[fPOI]-nominalPOIval,2));
                        npVal = fFitNPValues;
                        npVal[gammaName] = 1-statErr;
                        WriteDebugStatus("TRExFit::Fit","Setting "+gammaName+" to "+std::to_string(1-statErr));
                        data = DumpData( ws, regionDataType, npVal, fFitPOIAsimov );
                        npVal[gammaName] = 1;
                        ws->loadSnapshot("InitialStateModelGlob");
                        ws->loadSnapshot("InitialStateModelNuis");
                        npValues = PerformFit( ws, data, fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                        MCstatDoSample[sh->fSample->fName] = -sqrt(pow(MCstatDoSample[sh->fSample->fName],2) + pow(npValues[fPOI]-nominalPOIval,2));
                        smpTexTitle[sh->fSample->fName] = sh->fSample->fTexTitle;
                    }
                    MCstatUpSample[sh->fSample->fName]*=sh->fSample->fMCstatScale;
                    MCstatDoSample[sh->fSample->fName]*=sh->fSample->fMCstatScale;
                }
            }
        }
        //
        std::map < std::string, double > newPOIvalUp;
        std::map < std::string, double > newPOIvalDo;
        std::vector < std::string > npList;
        for(auto syst : fSystematics){
            if(FindInStringVector(npList,syst->fNuisanceParameter)<0){
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
                    data = DumpData( ws, regionDataType, npVal, fFitPOIAsimov );
                    // - again a stat-only fit to that Asimov
                    fFitFixedNPs[syst->fNuisanceParameter] = 0;
                    fFitFixedNPs["alpha_"+syst->fNuisanceParameter] = 0;
                    npValues = PerformFit( ws, data, fFitType, false, TRExFitter::DEBUGLEVEL<2 ? 0 : TRExFitter::DEBUGLEVEL);
                    double newPOIval = npValues[fPOI];
                    if(ud==0) newPOIvalUp[syst->fNuisanceParameter] = newPOIval;
                    if(ud==1) newPOIvalDo[syst->fNuisanceParameter] = newPOIval;
                }
                std::string category = syst->fCategory;
                if(syst->fSubCategory!="") category = syst->fSubCategory;
                if(FindInStringVector(systGroupNames,category)<0) systGroupNames.push_back(category);
                if(fabs(newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval)>fNonProfileFitSystThreshold
                || fabs(newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval)>fNonProfileFitSystThreshold)
                    systGroups[category] = sqrt(pow(systGroups[category],2)
                        +pow((fabs(newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval)+fabs(newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval))/2,2));
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
        totUp = sqrt(pow(totUp,2)+pow(MCstatUp,2));
        totDo = sqrt(pow(totDo,2)+pow(MCstatDo,2));
        // MC stat for separate gamma samples
        for(auto sepGammaPair : MCstatUpSample){
            std::string smpName = sepGammaPair.first;
            std::cout << "Stat." << smpName << "\t" << MCstatUpSample[smpName] << "\t" << MCstatDoSample[smpName] << std::endl;
            out       << "Stat." << smpName << "\t" << MCstatUpSample[smpName] << "\t" << MCstatDoSample[smpName] << std::endl;
            tex       << "  Stat (" << smpTexTitle[smpName] << ") & $+" << Form("%.2f",MCstatUpSample[smpName]) << "$ / $" << Form("%.2f",MCstatDoSample[smpName]) << "$ \\\\" << std::endl;
            tex       << "\\hline" << std::endl;
            totUp = sqrt(pow(totUp,2)+pow(MCstatUpSample[smpName],2));
            totDo = sqrt(pow(totDo,2)+pow(MCstatDoSample[smpName],2));
        }
        std::cout << "-----------------------------------" << std::endl;
        // systematics
        for(auto syst : fSystematics){
            if(syst->fName.find("stat_")!=std::string::npos && syst->fType==Systematic::SHAPE) continue;
            if(FindInStringVector(npList,syst->fNuisanceParameter)<0){
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
                        if(ud==0 && valUp>0) totUp = sqrt(pow(totUp,2)+pow(valUp,2));
                        if(ud==0 && valUp<0) totDo = sqrt(pow(totDo,2)+pow(valUp,2));
                        if(ud==1 && valDo>0) totUp = sqrt(pow(totUp,2)+pow(valDo,2));
                        if(ud==1 && valDo<0) totDo = sqrt(pow(totDo,2)+pow(valDo,2));
                    }
                }
                std::cout << std::endl;
                out       << std::endl;
                tex       << "\\\\" << std::endl;
            }
        }
        totUp = sqrt(pow(totUp,2)+pow(fNonProfileFitSystThreshold,2));
        totDo = sqrt(pow(totDo,2)+pow(fNonProfileFitSystThreshold,2));
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "TotalSystematic\t" << totUp << "\t-" << totDo << std::endl;
        out       << "TotalSystematic\t" << totUp << "\t-" << totDo << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "TotalStat+Syst\t" << sqrt(pow(totUp,2)+pow(statUp,2)) << "\t-" << sqrt(pow(totDo,2)+pow(statDo,2)) << std::endl;
        out       << "TotalStat+Syst\t" << sqrt(pow(totUp,2)+pow(statUp,2)) << "\t-" << sqrt(pow(totDo,2)+pow(statDo,2)) << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        out.close();
        tex << "\\hline" << std::endl;
        tex << "  Total systematics & $+" << Form("%.2f",totUp) << "$ / $-" << Form("%.2f",totDo) << "$ \\\\" << std::endl;
        tex << "\\hline" << std::endl;
        tex << "  Total stat+syst & $+"   << Form("%.2f",sqrt(pow(totUp,2)+pow(statUp,2)))     << "$ / $-" << Form("%.2f",sqrt(pow(totDo,2)+pow(statDo,2)))     << "$ \\\\" << std::endl;
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
        std::cout << "TotalSystematic\t" << (fabs(totUp)+fabs(totDo))/2. << std::endl;
        std::cout << "TotalStat+Syst\t" << (sqrt(pow(totUp,2)+pow(statUp,2))+sqrt(pow(totDo,2)+pow(statDo,2)))/2. << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        tex2 << "\\hline" << std::endl;
        tex2 << "Total systematic uncertainty & " << Form("$%.2f$",(fabs(totUp)+fabs(totDo))/2.) << " \\\\" << std::endl;
        tex2 << "\\hline" << std::endl;
        tex2 << "Total & "   << Form("$%.2f$",(sqrt(pow(totUp,2)+pow(statUp,2))+sqrt(pow(totDo,2)+pow(statDo,2)))/2.) << " \\\\" << std::endl;
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
                    GetLikelihoodScan( ws, it->first, data);
                }
            }
            else{
                for(unsigned int i=0; i<fVarNameLH.size(); ++i){
                    GetLikelihoodScan( ws, fVarNameLH[i], data);
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
                GetLikelihoodScan( ws, it->first, data);
            }
        } else {
            GetLikelihoodScan( ws, fVarNameLH[0], data);
        }
    }

    // run 2D likelihood scan
    if(fVarName2DLH.size()>0){
        for (const auto & ipair : fVarName2DLH) {
            Get2DLikelihoodScan( ws, ipair, data);
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
    for( const std::pair < std::string, int > dataType : regionDataType ){
        WriteDebugStatus("TRExFit::DumpData", "       - Region: " + dataType.first + "       DataType: " + std::to_string(dataType.second));
    }
    if(npValues.size()){
        WriteDebugStatus("TRExFit::DumpData", "    * Injected NP values ");
        for ( const std::pair < std::string, double > npValue : npValues ){
            WriteDebugStatus("TRExFit::DumpData", "       - NP: " + npValue.first + "       Value: " + std::to_string(npValue.second));
        }
    }
    else {
        WriteDebugStatus("TRExFit::DumpData", "    * No NP values injected ");
    }
    WriteDebugStatus("TRExFit::DumpData", "    * POI value: " + std::to_string(poiValue) );

    RooStats::ModelConfig *mc = (RooStats::ModelConfig*)ws -> obj("ModelConfig");

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

    //Setting binned likelihood option
    RooFIter rfiter = ws->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = rfiter.next())) {
        if (arg->IsA() == RooRealSumPdf::Class()) {
            arg->setAttribute("BinnedLikelihood");
            std::string temp_string = arg->GetName();
            WriteDebugStatus("TRExFit::DumpData", "Activating binned likelihood attribute for " + temp_string);
        }
    }

    //Creating a set
    const char* weightName="weightVar";
    RooArgSet obsAndWeight;
    obsAndWeight.add(*mc->GetObservables());

    RooRealVar* weightVar = NULL;
    if ( !(weightVar = ws->var(weightName)) ){
        ws->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
        weightVar = ws->var(weightName);
    }
    obsAndWeight.add(*ws->var(weightName));
    ws->defineSet("obsAndWeight",obsAndWeight);

    //
    // Getting observed data (in case some regions are unblinded)
    //
    RooDataSet* realData = (RooDataSet*)ws -> data("obsData");

    //
    // Set some parameters for the Asimov production
    //     |-> Values of NPs
    //     |-> Values of POI
    //

    //-- POI
    RooRealVar * poi = (RooRealVar*) mc->GetParametersOfInterest()->first();
    if (!poi){
        if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
        WriteErrorStatus("TRExFit::DumpData", "Cannot find POI in workspace, exiting...");
        exit(EXIT_FAILURE);
    }
    poi -> setVal(poiValue);

    //-- Nuisance parameters
    RooRealVar* var(nullptr);
    TIterator *npIterator = mc -> GetNuisanceParameters() -> createIterator();
    while( (var = (RooRealVar*) npIterator->Next()) ){
        std::map < std::string, double >::const_iterator it_npValue = npValues.find( var -> GetName() );
        if( it_npValue != npValues.end() ){
            var -> setVal(it_npValue -> second);
        }
    }

    //Looping over regions
    std::map<std::string, RooDataSet*> asimovDataMap;
    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
    RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;
    int iFrame = 0;
    int i = 0;
    while( (tt = (RooCatType*) iter -> Next()) ) {

        channelCat->setIndex(i);
        iFrame++;
        i++;

        //Check the type of data to store for this region !
        int dataType = Region::ASIMOVDATA;//default is AsimovData
        std::map < std::string, int >::const_iterator it_dataType = regionDataType.find( channelCat->getLabel() );
        if( it_dataType == regionDataType.end() ){
            std::string temp_string = channelCat->getLabel();
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
            RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
            double expectedEvents = pdftmp->expectedEvents(*obstmp);
            double thisNorm = 0;

            for(int jj=0; jj<thisObs->numBins(); ++jj){
                thisObs->setBin(jj);
                thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
                if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
            }
            obsDataUnbinned->Print();
            if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
                exit(1);
            }
            asimovDataMap[std::string(channelCat->getLabel())] = obsDataUnbinned;

        } else if(dataType==Region::REALDATA) {
            RooAbsData *datatmp = realData->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),tt->GetName()));
            asimovDataMap[std::string(channelCat->getLabel())] = (RooDataSet*)datatmp;
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
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        if (fNormFactors[i_norm]->fName == fPOI) continue;
        NPnames. emplace_back( fNormFactors[i_norm]->fName);
        NPvalues.emplace_back( fNormFactors[i_norm]->fNominal);
    }
    //
    // Fit configuration (SPLUSB or BONLY)
    //
    FittingTool *fitTool = new FittingTool();
    fitTool -> SetDebug(debugLevel);
    if(fitType==BONLY){
        fitTool -> ValPOI(0.);
        fitTool -> ConstPOI(true);
    } else if(fitType==SPLUSB){
        fitTool -> ValPOI(fFitPOIAsimov);
        fitTool -> ConstPOI(false);
    }
    fitTool -> SetNPs( NPnames,NPvalues );
    fitTool -> SetRandomNP(fRndRange, fUseRnd, fRndSeed);
    if(fStatOnly){
        if(!fGammasInStatOnly) fitTool -> NoGammas();
        fitTool -> NoSystematics();
    }

    //
    // Fit starting from custom point
    if(fFitResultsFile!=""){
        ReadFitResults(fFitResultsFile);
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(unsigned int i_np=0;i_np<fFitResults->fNuisPar.size();i_np++){
            npNames.push_back(  fFitResults->fNuisPar[i_np]->fName );
            npValues.push_back( fFitResults->fNuisPar[i_np]->fFitValue );
        }
        fitTool -> SetNPs( npNames,npValues );
    }

    //
    // Set Minos
    if(fVarNameMinos.size()>0){
        WriteDebugStatus("TRExFit::PerformFit", "Setting the variables to use MINOS with");
        fitTool -> UseMinos(fVarNameMinos);
    }

    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    if(inputData){
        data = inputData;
    } else {
        WriteWarningStatus("TRExFit::PerformFit", "You didn't provide inputData => will use the observed data !");
        data = (RooDataSet*)ws->data("obsData");
        if(data==nullptr){
            WriteWarningStatus("TRExFit::PerformFit", "No observedData found => will use the Asimov data !");
            data = (RooDataSet*)ws->data("asimovData");
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
        for(unsigned int i_np=0;i_np<fFitResults->fNuisPar.size();i_np++){
            if(!fFixNPforStatOnlyFit && FindInStringVector(fNormFactorNames,fFitResults->fNuisPar[i_np]->fName)>=0) continue;
            if(!fFixNPforStatOnlyFit && FindInStringVector(fShapeFactorNames,fFitResults->fNuisPar[i_np]->fName)>=0) continue;
            npNames.push_back(  fFitResults->fNuisPar[i_np]->fName );
            npValues.push_back( fFitResults->fNuisPar[i_np]->fFitValue );
        }
        fitTool -> FixNPs(npNames,npValues);
    }

    // FixNP
    if(fFitFixedNPs.size()>0){
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& nuisParToFix : fFitFixedNPs){
            npNames.push_back( nuisParToFix.first );
            npValues.push_back( nuisParToFix.second );
        }
        fitTool -> FixNPs(npNames,npValues);
    }

    // Tikhonov regularization (for unfolding)
    RooArgList l;
    std::vector<double> tauVec;
    for(auto nf : fNormFactors){
        if(nf->fTau!=0){
            l.add(*ws->var(nf->fName.c_str()));
            tauVec.push_back( nf->fTau );
        }
    }
    if(tauVec.size()>0){
        TMatrixDSym cov(tauVec.size());
        for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
            cov(i_tau,i_tau) = pow(1./tauVec[i_tau],2);
        }
        RooConstVar nominalValue("1","1",1);
        RooArgList nominal;
        for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++) nominal.add(nominalValue);
        RooMultiVarGaussian r("regularization","regularization",l,nominal,cov);
        ws->import(r);
        ws->defineSet("myConstraints","regularization");
        simPdf->setStringAttribute("externalConstraints","myConstraints");
        //
        const RooArgSet* externalConstraints = 0;
        if(simPdf->getStringAttribute("externalConstraints")){
            WriteInfoStatus("TRExFit::PerformFit",Form("Building NLL with external constraints %s",simPdf->getStringAttribute("externalConstraints")));
            externalConstraints = ws->set(simPdf->getStringAttribute("externalConstraints"));
            fitTool->SetExternalConstraints( externalConstraints );
        }
    }

    // save snapshot before fit
    ws->saveSnapshot("snapshot_BeforeFit_POI", *(mc->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_BeforeFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_BeforeFit_GO" , *(mc->GetGlobalObservables())    );

    //
    // Get initial ikelihood value from Asimov
    double nll0 = 0.;
    if (fBlindedParameters.size() > 0) std::cout.setstate(std::ios_base::failbit);
    if(fGetGoodnessOfFit) nll0 = fitTool -> FitPDF( mc, simPdf, (RooDataSet*)ws->data("asimovData"), false, true );

    // save snapshot before fit
    ws->saveSnapshot("snapshot_AfterFit_POI", *(mc->GetParametersOfInterest()) );
    ws->saveSnapshot("snapshot_AfterFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(mc->GetGlobalObservables())    );

    //
    // Get number of degrees of freedom
    // - number of bins
    int ndof = inputData->numEntries();
    // - minus number of free & non-constant parameters
    int nNF = 0;
    for(int i_nf=0;i_nf<fNNorm;i_nf++){
        if(fNormFactors[i_nf]->fConst) continue;
        if(fFitType==BONLY && fPOI==fNormFactors[i_nf]->fName) continue;
        // skip if it's a morphing parameter
        if(fNormFactors[i_nf]->fName.find("morph_")!=std::string::npos) continue;
        // skip if it has an "Expression"
        if(fNormFactors[i_nf]->fExpression.first!="") continue;
        // skip if not in the ws (e.g. because assigned to a sample or region not present in the fit)
        if(!ws->obj(fNormFactors[i_nf]->fName.c_str())) continue;
        nNF++;
    }
    ndof -= nNF;

    // Performs the fit
    fitTool -> MinimType("Minuit2");
    double nll = fitTool -> FitPDF( mc, simPdf, data );
    if (debugLevel < 1 && fBlindedParameters.size() == 0) std::cout.clear();
    if(save){
        if(fBootstrap!="" && fBootstrapIdx>=0){
            gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
            if(fStatOnlyFit) fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+"_statOnly.txt", fBlindedParameters);
            else             fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".txt", fBlindedParameters);
        }
        else{
            gSystem -> mkdir((fName+"/Fits/").c_str(),true);
            if(fStatOnlyFit) fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_statOnly.txt", fBlindedParameters);
            else             fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+".txt", fBlindedParameters);
        }
    }
    result = fitTool -> ExportFitResultInMap();
    if (fBlindedParameters.size() > 0) std::cout.clear();

    // If SaturatedModel used, superseed Asimov-based GOF
    if(fSaturatedModel && fGetGoodnessOfFit && !fDoGroupedSystImpactTable){
        ws->loadSnapshot("snapshot_BeforeFit_POI");
        ws->loadSnapshot("snapshot_BeforeFit_GO");
        ws->loadSnapshot("snapshot_BeforeFit_NP");
        //
        // perform fit to saturated model and store resulting nll as nll0
        nll0 = fitTool -> FitPDF( mc, simPdf, data, false, false, true );
        // could be removed, but kept for debugging FIXME
        if(save){
            if(fStatOnlyFit) fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_saturatedModel_statOnly.txt", fBlindedParameters);
            else             fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_saturatedModel.txt", fBlindedParameters);
        }
    }

    //
    // Goodness of fit
    if(fGetGoodnessOfFit && !fDoGroupedSystImpactTable){
        double deltaNLL = nll-nll0;
        double prob = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndof);
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
            gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
            outNameGroupedImpact = fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+"GroupedImpact"+fSuffix;
        }
        if(fGroupedImpactCategory!="all") outNameGroupedImpact += "_"+fGroupedImpactCategory;
        outNameGroupedImpact += ".txt";

        ProduceSystSubCategoryMap();                        // fill fSubCategoryImpactMap first
        fitTool -> SetSystMap( fSubCategoryImpactMap );     // hand over the map to the FittingTool
        fitTool -> GetGroupedImpact( mc, simPdf, data, ws, fGroupedImpactCategory, outNameGroupedImpact);
    }

    delete fitTool;
    delete fFitResults;
    return result;
}

//__________________________________________________________________________________
//
RooWorkspace* TRExFit::PerformWorkspaceCombination( std::vector < std::string > &regionsToFit ) const{

    //
    // Definition of the fit regions
    //
    std::vector < RooWorkspace* > vec_ws;
    std::vector < std::string > vec_chName;
    RooStats::HistFactory::Measurement *measurement = 0;
    //
    // Take the measurement from the combined workspace, to be sure to have all the systematics (even the ones which are not there in the first region)
    std::unique_ptr<TFile> rootFileCombined;
    if(fBootstrap!="" && fBootstrapIdx>=0)
        rootFileCombined = std::make_unique<TFile>( (fName+"/RooStats/"+fBootstrapSyst+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read");
    else
        rootFileCombined = std::make_unique<TFile>( (fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read");
    if(rootFileCombined!=nullptr) measurement = (RooStats::HistFactory::Measurement*) rootFileCombined -> Get( (fInputName+fSuffix).c_str());
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        bool isToFit = false;
        for(unsigned int iRegion = 0; iRegion < regionsToFit.size(); ++iRegion){
            if(fRegions[i_ch] -> fName == regionsToFit[iRegion]){
                isToFit = true;
                break;
            }
        }
        if(isToFit){
            std::string fileName = fName+"/RooStats/"+fInputName+"_"+fRegions[i_ch]->fName+"_"+fInputName+fSuffix+"_model.root";
            if(fBootstrap!="" && fBootstrapIdx>=0) fileName = fName+"/RooStats/"+fBootstrapSyst+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName+"_"+fRegions[i_ch]->fName+"_"+fInputName+fSuffix+"_model.root";
            TFile *rootFile = new TFile(fileName.c_str(),"read");
            RooWorkspace* m_ws = (RooWorkspace*) rootFile->Get((fRegions[i_ch]->fName).c_str());
            if(!m_ws){
                WriteErrorStatus("TRExFit::PerformWorkspaceCombination", "The workspace (\"" + fRegions[i_ch] -> fName + "\") cannot be found in file " + fileName + ". Please check !");
            }
            vec_ws.push_back(m_ws);
            vec_chName.push_back(fRegions[i_ch] -> fName);
            // if failed to get the measurement from the combined ws, take it from the first region
            if(!measurement){
                measurement = (RooStats::HistFactory::Measurement*) rootFile -> Get( (fInputName+fSuffix).c_str());
            }
        }
    }

    //
    // Create the HistoToWorkspaceFactoryFast object to perform safely the combination
    //
    if(!measurement){
        WriteErrorStatus("TRExFit::PerformWorkspaceCombination", "The measurement object has not been retrieved ! Please check.");
        return 0;
    }
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    RooStats::HistFactory::HistoToWorkspaceFactoryFast factory(*measurement);

    // Creating the combined model
    RooWorkspace* ws = factory.MakeCombinedModel( vec_chName, vec_ws );

    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws, *measurement );
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    rootFileCombined->Close();

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
    ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
    if(fFitResults){
        fFitResults->fNuisParToHide = fVarNameHide;
        std::set < std::string > npCategories;
        for(unsigned int i=0;i<fSystematics.size();i++){
            npCategories.insert(fSystematics[i]->fCategory);
        }
        for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
            if(!fStatOnly && !fStatOnlyFit){
                fFitResults->DrawNPPulls(fName+"/NuisPar"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format],"all",fNormFactors, fBlindedParameters);
                fFitResults->DrawGammaPulls(fName+"/Gammas"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format], fBlindedParameters);
            }
            fFitResults->DrawNormFactors(fName+"/NormFactors"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format],fNormFactors, fBlindedParameters);
        }
        if(npCategories.size()>1 && !fStatOnly && !fStatOnlyFit){
            for( const std::string cat : npCategories ){
                std::string cat_for_name = cat;
                std::replace( cat_for_name.begin(), cat_for_name.end(), ' ', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '#', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '{', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '}', '_');
                for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
                  fFitResults->DrawNPPulls(fName+"/NuisPar_"+cat_for_name+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format],cat,fNormFactors, fBlindedParameters);
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
        for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
            fFitResults->DrawCorrelationMatrix(fName+"/CorrMatrix"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format],
                                               fuseGammasForCorr, TRExFitter::CORRELATIONTHRESHOLD);
    }
}

//__________________________________________________________________________________
//
void TRExFit::GetLimit(){

    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    //
    // If a workspace file name is specified, do simple limit
    //
    int sigDebug = 3 - TRExFitter::DEBUGLEVEL;
    if (sigDebug < 0) sigDebug = 0;
    if(fWorkspaceFileName!=""){
        std::string dataName = "obsData";
        if(!hasData || fLimitIsBlind) dataName = "asimovData";
        runAsymptoticsCLs(fWorkspaceFileName.c_str(), "combined", "ModelConfig", dataName.c_str(), fLimitParamName.c_str(), fLimitParamValue, (fLimitOutputPrefixName+fSuffix).c_str(), (fName+"/Limits/").c_str(), fLimitIsBlind, fLimitsConfidence, "asimovData_0", fSignalInjection, fSignalInjectionValue, sigDebug);
    }
    else{
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std:: string > regionsForFit;
        std::vector < std::string > regionsForLimit;
        std::map < std::string, int > regionsForFitDataType;
        std::map < std::string, int > regionsForLimitDataType;
        bool onlyUseRealData = true;
        for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
            if( fRegions[i_ch] -> fRegionType == Region::VALIDATION ) continue;
            if( hasData && fRegions[i_ch] -> fRegionDataType == Region::REALDATA && !fLimitIsBlind ){
                Region::DataType dataType = fRegions[i_ch] -> fRegionDataType;
                regionsForFit.push_back( fRegions[i_ch] -> fName );
                regionsForFitDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , dataType) );
            }
            regionsForLimit.push_back(fRegions[i_ch] -> fName);
            regionsForLimitDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , (fLimitIsBlind || !hasData) ? Region::ASIMOVDATA : fRegions[i_ch] -> fRegionDataType) );
            if(fLimitIsBlind || !hasData || fRegions[i_ch] -> fRegionDataType == Region::ASIMOVDATA){
                onlyUseRealData = false;
            }
        }

        std::map < std::string, double > npValues;
        RooDataSet* data = 0;

        if(regionsForFit.size()>0 && !onlyUseRealData){
            //
            // Creates a combined workspace with the regions to be used *in the fit*
            //
            WriteInfoStatus("TRExFit::GetLimit","Creating ws for regions with real data only...");
            RooWorkspace* ws_forFit = PerformWorkspaceCombination( regionsForFit );
            if (!ws_forFit){
                WriteErrorStatus("TRExFit::GetLimit","Cannot retrieve the workspace, exiting!");
                exit(EXIT_FAILURE);
            }

            //
            // Calls the PerformFit() function to actually do the fit
            //
            WriteInfoStatus("TRExFit::GetLimit","Performing a fit in regions with real data only...");
            npValues = PerformFit( ws_forFit, data, FitType::BONLY, false, TRExFitter::DEBUGLEVEL);
            WriteInfoStatus("TRExFit::GetLimit","Now will use the fit results to create the Asimov in the regions without real data!");
        }

        //
        // Create the final asimov dataset for limit setting
        //
        RooWorkspace* ws_forLimit = PerformWorkspaceCombination( regionsForLimit );
        if (!ws_forLimit){
            WriteErrorStatus("TRExFit::GetLimit","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        data = DumpData( ws_forLimit, regionsForLimitDataType, npValues, npValues.find(fPOI)==npValues.end() ? fLimitPOIAsimov : npValues[fPOI] );

        //
        // Set all saturated model factors to constant
        RooRealVar* var = nullptr;
        RooArgSet vars = ws_forLimit->allVars();
        TIterator* it = vars.createIterator();
        while( (var = (RooRealVar*) it->Next()) ){
            std::string name = var->GetName();
            if(name.find("saturated_model_sf_")!=std::string::npos){
                WriteInfoStatus("TRExFit::GetLimit","Fixing parameter " + name );
                var->setConstant( 1 );
            }
        }

        //
        // Gets the measurement object in the original combined workspace (created with the "w" command)
        //
        const std::string originalCombinedFile = fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root";
        TFile *f_origin = new TFile(originalCombinedFile.c_str(), "read");
        RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get((fInputName+fSuffix).c_str());
        TString outputName = f_origin->GetName();
        f_origin -> Close();

        //
        // Creating the rootfile used as input for the limit setting :-)
        //
        outputName = outputName.ReplaceAll(".root","_forLimits.root");
        TFile *f_clone = new TFile( outputName, "recreate" );
        ws_forLimit -> import(*data,Rename("ttHFitterData"));
        originalMeasurement -> Write();
        ws_forLimit -> Write();
        f_clone -> Close();
        std::string outputName_s = static_cast<std::string> (outputName);
        runAsymptoticsCLs(outputName_s.c_str(), "combined", "ModelConfig", "ttHFitterData", fLimitParamName.c_str(), fLimitParamValue, (fLimitOutputPrefixName+fSuffix).c_str(), (fName+"/Limits/").c_str(), fLimitIsBlind, fLimitsConfidence, "asimovData_0", fSignalInjection, fSignalInjectionValue, sigDebug);
    }
}

//__________________________________________________________________________________
//
void TRExFit::GetSignificance(){

    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    //
    // If a workspace file name is specified, do simple significance
    //
    int sigDebug = 3 - TRExFitter::DEBUGLEVEL;
    if (sigDebug < 0) sigDebug = 0;
    if(fWorkspaceFileName!=""){
        std::string dataName = "obsData";
        if(!hasData || fSignificanceIsBlind) dataName = "asimovData";
        runSig(fWorkspaceFileName.c_str(), "combined", "ModelConfig", dataName.c_str(), fSignificanceParamName.c_str(), fSignificanceParamValue, fSignificanceOutputPrefixName.c_str(), (fName+"/Significance").c_str(), fSignificanceIsBlind, "asimovData_1", "conditionalGlobs_1", "nominalGlobs", false, fSignificancePOIAsimov, sigDebug);
    }
    else{
        //
        // Fills a vector of regions to consider for fit
        //
        std::vector < std:: string > regionsForFit;
        std::vector < std::string > regionsForSign;
        std::map < std::string, int > regionsForFitDataType;
        std::map < std::string, int > regionsForSignDataType;
        bool onlyUseRealData = true;
        for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
            if( fRegions[i_ch] -> fRegionType == Region::VALIDATION ) continue;
            if( hasData && fRegions[i_ch] -> fRegionDataType == Region::REALDATA && !fSignificanceIsBlind){
                Region::DataType dataType = fRegions[i_ch] -> fRegionDataType;
                regionsForFit.push_back( fRegions[i_ch] -> fName );
                regionsForFitDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , dataType) );
            }
            regionsForSign.push_back(fRegions[i_ch] -> fName);
            regionsForSignDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , (!hasData || fSignificanceIsBlind) ? Region::ASIMOVDATA : fRegions[i_ch] -> fRegionDataType) );
            if(fSignificanceIsBlind || !hasData || fRegions[i_ch] -> fRegionDataType == Region::ASIMOVDATA){
                onlyUseRealData = false;
            }
        }

        std::map < std::string, double > npValues;
        RooDataSet* data = 0;
        if(regionsForFit.size()>0 && !onlyUseRealData){
            //
            // Creates a combined workspace with the regions to be used *in the fit*
            //
            WriteInfoStatus("TRExFit::GetSignificance","Creating ws for regions with real data only...");
            RooWorkspace* ws_forFit = PerformWorkspaceCombination( regionsForFit );
            if (!ws_forFit){
                WriteErrorStatus("TRExFit::GetSignificance","Cannot retrieve the workspace, exiting!");
                exit(EXIT_FAILURE);
            }

            //
            // Calls the PerformFit() function to actually do the fit
            //
            WriteInfoStatus("TRExFit::GetSignificance","Performing a fit in regions with real data only...");
            npValues = PerformFit( ws_forFit, data, FitType::BONLY, false, TRExFitter::DEBUGLEVEL);
            WriteInfoStatus("TRExFit::GetSignificance","Now will use the fit results to create the Asimov in the regions without real data!");
        }

        //
        // Create the final asimov dataset for limit setting
        //
        RooWorkspace* ws_forSignificance = PerformWorkspaceCombination( regionsForSign );
        if (!ws_forSignificance){
            WriteErrorStatus("TRExFit::GetSignificance","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        data = DumpData( ws_forSignificance, regionsForSignDataType, npValues, npValues.find(fPOI)==npValues.end() ? fSignificancePOIAsimov : npValues[fPOI] );

        //
        // Set all saturated model factors to constant
        RooRealVar* var = nullptr;
        RooArgSet vars = ws_forSignificance->allVars();
        TIterator* it = vars.createIterator();
        while( (var = (RooRealVar*) it->Next()) ){
            std::string name = var->GetName();
            if(name.find("saturated_model_sf_")!=std::string::npos){
                WriteInfoStatus("TRExFit::GetSignificance","Fixing parameter " + name );
                var->setConstant( 1 );
            }
        }

        //
        // Gets the measurement object in the original combined workspace (created with the "w" command)
        //
        const std::string originalCombinedFile = fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root";
        TFile *f_origin = new TFile(originalCombinedFile.c_str(), "read");
        RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get((fInputName + fSuffix).c_str());
        TString outputName = f_origin->GetName();
        f_origin -> Close();

        //
        // Creating the rootfile used as input for the limit setting :-)
        //
        outputName = outputName.ReplaceAll(".root","_forSignificance.root");
        TFile *f_clone = new TFile( outputName, "recreate" );
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
    delete fFitResults;
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
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if (std::find(tmp.begin(), tmp.end(), fSystematics[i_syst]->fName) == tmp.end()){
            WriteInfoStatus("TRExFit::PrintConfigSummary"," "+fSystematics[i_syst]->fName);
            tmp.emplace_back(fSystematics[i_syst]->fName);
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
Sample* TRExFit::GetSample(const std::string& name) const{
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
            bkg->Add(fRegions[i_ch]->fBkg[i_bkg]->fHist); // add the rest
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

        if(fAtlasLabel!="none") ATLASLabelNew(0.20,0.84+0.04,(char*)(fAtlasLabel+"  Simulation").c_str(), kBlack, gStyle->GetTextSize());

        std::ostringstream SEP;
        SEP.precision(3);
        SEP << "Separation: " << GetSeparation(sig.get(),bkg.get())*100 << "%";
        myText(0.55,0.73,1,SEP.str().c_str());

        for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
            dummy3.SaveAs((fName+"/Plots/Separation/"+fRegions[i_ch]->fName+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format] ).c_str());

    }// regions

   return;
}

//____________________________________________________________________________________
//
void TRExFit::ProduceNPRanking( std::string NPnames/*="all"*/ ){

    if(fFitType==BONLY){
        WriteErrorStatus("TRExFit::ProduceNPRanking", "For ranking plots, the SPLUSB FitType is needed.");
        abort();
    }

    //
    // List of systematics to check
    //
    std::vector< std::string > nuisPars;
    std::vector< bool > isNF;
    std::vector<std::string> systNames_unique;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(NPnames=="all" || NPnames==fSystematics[i_syst]->fNuisanceParameter ||
            ( atoi(NPnames.c_str())==i_syst && (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0) )
            ){
            if(fSystematics[i_syst]->fType == Systematic::SHAPE) continue;
            if (std::find(systNames_unique.begin(), systNames_unique.end(), fSystematics[i_syst]->fNuisanceParameter) == systNames_unique.end())
                systNames_unique.push_back(fSystematics[i_syst]->fNuisanceParameter);
            else continue;
            nuisPars.push_back( fSystematics[i_syst]->fNuisanceParameter );
            isNF.push_back( false );
        }
    }
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        if(fPOI==fNormFactors[i_norm]->fName) continue;
        if(NPnames=="all" || NPnames==fNormFactors[i_norm]->fName ||
            ( atoi(NPnames.c_str())-fNSyst==i_norm && (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0) )
            ){
            nuisPars.push_back( fNormFactors[i_norm]->fName );
            isNF.push_back( true );
        }
    }

    //
    // Text files containing information necessary for drawing of ranking plot
    //
    std::string outName = fName+"/Fits/NPRanking"+fSuffix;
    if(fBootstrap!="" && fBootstrapIdx>=0){
        gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
        outName = fName+"/Fits/"+fBootstrapSyst+Form("_BSId%d/",fBootstrapIdx)+"NPRanking"+fSuffix;
    }
    if(NPnames!="all") outName += "_"+NPnames;
    outName += ".txt";
    std::ofstream outName_file(outName.c_str());
    //
    double central;
    double up;
    double down;
    double muhat;
    double dMuUp, dMuDown;
    std::map< std::string,double > muVarUp;
    std::map< std::string,double > muVarDown;
    std::map< std::string,double > muVarNomUp;
    std::map< std::string,double > muVarNomDown;

    //
    // Fills a vector of regions to consider for fit
    //
    std::vector < std:: string > regionsToFit;
    std::map < std::string, int > regionDataType;
    for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
        bool isToFit = false;

        if ( fFitRegion == CRONLY ) {
            if( fRegions[i_ch] -> fRegionType == Region::CONTROL ){
                isToFit = true;
            }
        } else if ( fFitRegion == CRSR ){
            if( fRegions[i_ch] -> fRegionType == Region::CONTROL || fRegions[i_ch] -> fRegionType == Region::SIGNAL ){
                isToFit = true;
            }
        }
        if ( ! isToFit ){
            for (unsigned int iReg = 0; iReg < fFitRegionsToFit.size(); ++iReg ){
                if( fFitRegionsToFit[iReg] == fRegions[i_ch] -> fName ){
                    isToFit = true;
                    break;
                }
            }
        }

        if(isToFit){
            regionsToFit.push_back( fRegions[i_ch] -> fName );
            Region::DataType dataType;
            if(fFitIsBlind){
                dataType = Region::ASIMOVDATA;
            } else {
                dataType = fRegions[i_ch] -> fRegionDataType;
            }
            regionDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , dataType) );
        }
    }

    //
    // Creating the combined model
    //
    std::unique_ptr<TFile> customWSfile = nullptr;
    RooWorkspace* ws = nullptr;
    if (fWorkspaceFileName!="") { // has custom worspace
        customWSfile = std::make_unique<TFile>(fWorkspaceFileName.c_str(),"read");
        ws = static_cast<RooWorkspace*>(customWSfile->Get("combined"));
    } else {
        ws = PerformWorkspaceCombination( regionsToFit );
    }
    if (!ws){
        WriteErrorStatus("TRExFit::ProduceNPRanking","Cannot retrieve the workspace, exiting!");
        exit(EXIT_FAILURE);
    }

    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig *mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());
    RooDataSet* data = nullptr;

    if (fWorkspaceFileName!=""){
        bool hasData = false;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType==Sample::DATA){
                hasData = true;
                break;
            }
        }

        if(!fFitIsBlind && hasData) data = static_cast<RooDataSet*>(ws->data("obsData"));
        else                        data = static_cast<RooDataSet*>(ws->data("asimovData"));
    } else {
        data = DumpData( ws, regionDataType, fFitNPValues, fFitPOIAsimov );
    }

    if (!mc || !simPdf || !data){
        WriteErrorStatus("TRExFit::ProduceNPRanking","At least one of the objects that is needed to run ranking is not present");
        exit(EXIT_FAILURE);
    }

    // Loop on NPs to find gammas and add to the list to be ranked
    if(NPnames=="all" || NPnames.find("gamma")!=std::string::npos || (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0)){
        RooRealVar* var = NULL;
        RooArgSet* nuis = (RooArgSet*) mc->GetNuisanceParameters();
        if(nuis){
            TIterator* it2 = nuis->createIterator();
            int i_gamma = 0;
            while( (var = (RooRealVar*) it2->Next()) ){
                std::string np = var->GetName();
                if(np.find("gamma")!=std::string::npos){
                    // add the nuisance parameter to the list nuisPars if it's there in the ws
                    // remove "gamma"...
                    if(np==NPnames || (atoi(NPnames.c_str())-fNSyst-fNNorm==i_gamma && (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0)) || NPnames=="all"){
                        nuisPars.push_back(ReplaceString(np,"gamma_",""));
                        isNF.push_back( true );
                        if(NPnames!="all") break;
                    }
                    i_gamma++;
                }
            }
        }
    }

    //
    // Create snapshot to keep inital values
    //
    ws -> saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));

    //
    // Initialize the FittingTool object
    //
    FittingTool *fitTool = new FittingTool();
    fitTool -> SetDebug(TRExFitter::DEBUGLEVEL);
    fitTool -> ValPOI(fFitPOIAsimov);
    fitTool -> ConstPOI(false);
    if(fStatOnly){
        fitTool -> NoGammas();
        fitTool -> NoSystematics();
    }

    // Set initial NP to random value if specified
    fitTool -> SetRandomNP(fRndRange, fUseRnd, fRndSeed);

    ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
    {
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(int i_norm=0;i_norm<fNNorm;i_norm++){
            if (fNormFactors[i_norm]->fName == fPOI) continue;
            npNames. emplace_back( fNormFactors[i_norm]->fName);
            npValues.emplace_back( fNormFactors[i_norm]->fNominal);
        }
        fitTool -> SetNPs( npNames,npValues );
    }

    muhat = fFitResults -> GetNuisParValue( fPOI );

    for(unsigned int i=0;i<nuisPars.size();i++){
        //
        // Getting the postfit values of the nuisance parameter
        central = fFitResults -> GetNuisParValue(   nuisPars[i] );
        up      = fFitResults -> GetNuisParErrUp(   nuisPars[i] );
        down    = fFitResults -> GetNuisParErrDown( nuisPars[i] );
        //// Thomas : We should be careful with changing naming convention compared to RooFit !!
        // TRExFitter store gammas names as stat_Reg_bin_i (i.e. remove the gamma_ at the beginning)
        // Now there is no real identifier in the NP name to state if it is a gamma or not and add back gamma_ except this _bin_
        if( (nuisPars[i].find("_bin_")!=std::string::npos) ){
            nuisPars[i] = "gamma_" + nuisPars[i];
        }
        outName_file <<  nuisPars[i] << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";
        //
        // Experimental: reduce the range of ranking
        if(TRExFitter::OPTION["ReduceRanking"]!=0){
            up   *= TRExFitter::OPTION["ReduceRanking"];
            down *= TRExFitter::OPTION["ReduceRanking"];
        }
        //
        // Set the NP to its post-fit *up* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool -> ResetFixedNP();
        // fix NPs that are fixed in the config
        // this has nothing to do with fixing NPs for the ranking
        // this is just needed to be compatible with the normal fit
        if(fFitFixedNPs.size()>0){
            for(const auto& nuisParToFix : fFitFixedNPs){
                fitTool -> FixNP(nuisParToFix.first,nuisParToFix.second);
            }
        }

        fitTool -> FixNP( nuisPars[i], central + TMath::Abs(up  ) );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
        //
        // Set the NP to its post-fit *down* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool -> ResetFixedNP();
        fitTool -> FixNP( nuisPars[i], central - TMath::Abs(down) );
        if(fFitFixedNPs.size()>0){
            for(const auto& nuisParToFix : fFitFixedNPs){
                fitTool -> FixNP(nuisParToFix.first,nuisParToFix.second);
            }
        }
        fitTool -> FitPDF( mc, simPdf, data );
        muVarDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        //
        dMuUp   = muVarUp[nuisPars[i]]-muhat;
        dMuDown = muVarDown[nuisPars[i]]-muhat;
        //
        // Experimental: reduce the range of ranking
        if(TRExFitter::OPTION["ReduceRanking"]!=0){
            dMuUp   /= TRExFitter::OPTION["ReduceRanking"];
            dMuDown /= TRExFitter::OPTION["ReduceRanking"];
        }
        //
        outName_file << dMuUp << "   " << dMuDown << "  ";

        if(isNF[i]){
            muVarNomUp[   nuisPars[i] ] = muhat;
            muVarNomDown[ nuisPars[i] ] = muhat;
        }
        else{
            up   = 1.;
            down = 1.;
            //
            // Experimental: reduce the range of ranking
            if(TRExFitter::OPTION["ReduceRanking"]!=0){
                up   *= TRExFitter::OPTION["ReduceRanking"];
                down *= TRExFitter::OPTION["ReduceRanking"];
            }
            //
            // Set the NP to its pre-fit *up* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool -> ResetFixedNP();
            fitTool -> FixNP( nuisPars[i], central + TMath::Abs(up  ) );
            fitTool -> FitPDF( mc, simPdf, data );
            if(fFitFixedNPs.size()>0){
                for(const auto& nuisParToFix : fFitFixedNPs){
                    fitTool -> FixNP(nuisParToFix.first,nuisParToFix.second);
                }
            }
            muVarNomUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
            //
            // Set the NP to its pre-fit *down* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool -> ResetFixedNP();
            fitTool -> FixNP( nuisPars[i], central - TMath::Abs(down) );
            if(fFitFixedNPs.size()>0){
                for(const auto& nuisParToFix : fFitFixedNPs){
                    fitTool -> FixNP(nuisParToFix.first,nuisParToFix.second);
                }
            }
            fitTool -> FitPDF( mc, simPdf, data );
            //
            muVarNomDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        }
        dMuUp   = muVarNomUp[nuisPars[i]]-muhat;
        dMuDown = muVarNomDown[nuisPars[i]]-muhat;
        //
        // Experimental: reduce the range of ranking
        if(TRExFitter::OPTION["ReduceRanking"]!=0){
            dMuUp   /= TRExFitter::OPTION["ReduceRanking"];
            dMuDown /= TRExFitter::OPTION["ReduceRanking"];
        }
        //
       outName_file << dMuUp << "   " << dMuDown << " " << std::endl;

    }
    outName_file.close();
    ws->loadSnapshot("tmp_snapshot");
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
void TRExFit::PlotNPRanking(bool flagSysts, bool flagGammas) const{
    //
    std::string fileToRead = fName+"/Fits/NPRanking"+fSuffix+".txt";
    //
    // trick to merge the ranking outputs produced in parallel:
    std::string cmd = " if [[ `ls "+fName+"/Fits/NPRanking"+fSuffix+"_*` != \"\" ]] ; then";
    cmd       += " if [[ ! -f "+fName+"/Fits/NPRanking"+fSuffix+".txt ]] ; then";
    cmd       += " cat "+fName+"/Fits/NPRanking"+fSuffix+"_* > "+fileToRead+" ; ";
    cmd       += " fi ;";
    cmd       += " fi ;";
    gSystem->Exec(cmd.c_str());
    //
    unsigned int maxNP = fRankingMaxNP;
    //
    std::string paramname;
    double nuiphat;
    double nuiperrhi;
    double nuiperrlo;
    double PoiUp;
    double PoiDown;
    double PoiNomUp;
    double PoiNomDown;
    std::vector<std::string> parname;
    std::vector<double> nuhat;
    std::vector<double> nuerrhi;
    std::vector<double> nuerrlo;
    std::vector<double> poiup;
    std::vector<double> poidown;
    std::vector<double> poinomup;
    std::vector<double> poinomdown;
    std::vector<double> number;

    std::ifstream fin( fileToRead.c_str() );
    fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    std::string temp_string = "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"";
    if (paramname=="Luminosity"){
        WriteErrorStatus("TRExFit::PlotNPRanking", temp_string);
        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    }
    while (!fin.eof()){
        if(paramname.find("gamma")!=std::string::npos && !flagGammas){
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            if (paramname=="Luminosity"){
                WriteErrorStatus("TRExFit::PlotNPRanking", temp_string);
                fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            }
            continue;
        }
        if(paramname.find("gamma")==std::string::npos && !flagSysts){
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            if (paramname=="Luminosity"){
                WriteErrorStatus("TRExFit::PlotNPRanking", temp_string);
                fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            }
            continue;
        }
        parname.push_back(paramname);
        nuhat.push_back(nuiphat);
        nuerrhi.push_back(nuiperrhi);
        nuerrlo.push_back(nuiperrlo);
        poiup.push_back(PoiUp);
        poidown.push_back(PoiDown);
        poinomup.push_back(PoiNomUp);
        poinomdown.push_back(PoiNomDown);
        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        if (paramname=="Luminosity"){
            WriteErrorStatus("TRExFit::PlotNPRanking", temp_string);
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        }
    }

    unsigned int SIZE = parname.size();
    WriteDebugStatus("TRExFit::PlotNPRanking", "NP ordering...");
    number.push_back(0.5);
    for (unsigned int i=1;i<SIZE;i++){
        number.push_back(i+0.5);
        double sumi = 0.0;
        int index=-1;
        sumi += TMath::Max( TMath::Abs(poiup[i]),TMath::Abs(poidown[i]) );
        for (unsigned int j=1;j<=i;j++){
            double sumii = 0.0;
            sumii += TMath::Max( TMath::Abs(poiup[i-j]),TMath::Abs(poidown[i-j]) );
            if (sumi<sumii){
                if (index==-1){
                    std::swap(poiup[i],poiup[i-j]);
                    std::swap(poidown[i],poidown[i-j]);
                    std::swap(poinomup[i],poinomup[i-j]);
                    std::swap(poinomdown[i],poinomdown[i-j]);
                    std::swap(nuhat[i],nuhat[i-j]);
                    std::swap(nuerrhi[i],nuerrhi[i-j]);
                    std::swap(nuerrlo[i],nuerrlo[i-j]);
                    std::swap(parname[i],parname[i-j]);
                    index=i-j;
                }
                else{
                    std::swap(poiup[index],poiup[i-j]);
                    std::swap(poidown[index],poidown[i-j]);
                    std::swap(poinomup[index],poinomup[i-j]);
                    std::swap(poinomdown[index],poinomdown[i-j]);
                    std::swap(nuhat[index],nuhat[i-j]);
                    std::swap(nuerrhi[index],nuerrhi[i-j]);
                    std::swap(nuerrlo[index],nuerrlo[i-j]);
                    std::swap(parname[index],parname[i-j]);
                    index=i-j;
                }
            }
            else{
                break;
            }
        }
    }
    number.push_back(parname.size()-0.5);

    double poimax = 0;
    for (unsigned int i=0;i<SIZE;i++) {
        poimax = TMath::Max(poimax,TMath::Max( TMath::Abs(poiup[i]),TMath::Abs(poidown[i]) ));
        poimax = TMath::Max(poimax,TMath::Max( TMath::Abs(poinomup[i]),TMath::Abs(poinomdown[i]) ));
        nuerrlo[i] = TMath::Abs(nuerrlo[i]);
    }
    poimax *= 1.2;

    for (unsigned int i=0;i<SIZE;i++) {
        poiup[i]     *= (2./poimax);
        poidown[i]   *= (2./poimax);
        poinomup[i]  *= (2./poimax);
        poinomdown[i]*= (2./poimax);
    }

    // Restrict to the first N
    if(SIZE>maxNP) SIZE = maxNP;

    // Graphical part - rewritten taking DrawPulls in TRExFitter
    double lineHeight  =  30.;
    double offsetUp    =  60.; // external
    double offsetDown  =  60.;
    double offsetUp1   = 100.; // internal
    double offsetDown1 =  15.;
    int offset = offsetUp + offsetDown + offsetUp1 + offsetDown1;
    int newHeight = offset + SIZE*lineHeight;

    double xmin = -2.;
    double xmax =  2.;
    double max  =  0.;

    TGraphAsymmErrors g{};
    TGraphAsymmErrors g1{};
    TGraphAsymmErrors g2{};
    TGraphAsymmErrors g1a{};
    TGraphAsymmErrors g2a{};

    int idx = 0;
    std::vector< std::string > Names;
    std::string parTitle;

    for(unsigned int i = parname.size()-SIZE; i<parname.size(); ++i){
        g.SetPoint(idx, nuhat[i],  idx+0.5);
        g.SetPointEXhigh(      idx, nuerrhi[i]);
        g.SetPointEXlow(       idx, nuerrlo[i]);

        g1.SetPoint(      idx, 0.,idx+0.5);
        g1.SetPointEXhigh(idx, poiup[i]);
        g1.SetPointEXlow( idx, 0.);
        g1.SetPointEYhigh(idx, 0.4);
        g1.SetPointEYlow( idx, 0.4);

        g2.SetPoint(      idx, 0.,idx+0.5);
        g2.SetPointEXhigh(idx, poidown[i]);
        g2.SetPointEXlow( idx, 0.);
        g2.SetPointEYhigh(idx, 0.4);
        g2.SetPointEYlow( idx, 0.4);

        g1a.SetPoint(      idx, 0.,idx+0.5);
        g1a.SetPointEXhigh(idx, poinomup[i]);
        g1a.SetPointEXlow( idx, 0.);
        g1a.SetPointEYhigh(idx, 0.4);
        g1a.SetPointEYlow( idx, 0.4);

        g2a.SetPoint(      idx, 0.,idx+0.5);
        g2a.SetPointEXhigh(idx, poinomdown[i]);
        g2a.SetPointEXlow( idx, 0.);
        g2a.SetPointEYhigh(idx, 0.4);
        g2a.SetPointEYlow( idx, 0.4);

        if(parname[i].find("gamma")!=std::string::npos){
            // get name of the region
            std::vector<std::string> tmpVec = Vectorize(parname[i],'_');
            int nWords = tmpVec.size();
            std::string regName = tmpVec[2];
            for(int i_word=3;i_word<nWords-2;i_word++){
                regName += tmpVec[i_word];
            }
            // find the short label of this region
            std::string regTitle = regName;
            for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
                if(fRegions[i_ch]->fName==regName){
                    regTitle = fRegions[i_ch]->fShortLabel;
                    break;
                }
            }
            // build the title of the nuis par
            parTitle = "#gamma (" + regTitle + " bin " + tmpVec[nWords-1] + ")";
        }
        else parTitle = TRExFitter::SYSTMAP[ parname[i] ];

        Names.push_back(parTitle);

        idx ++;
        if(idx > max)  max = idx;
    }
    int newWidth = 600;
    if (fNPRankingCanvasSize.size() != 0){
        newWidth = fNPRankingCanvasSize.at(0);
        newHeight = fNPRankingCanvasSize.at(1);
    }
    TCanvas c("c","c",newWidth,newHeight);
    c.SetTicks(0,0);
    gPad->SetLeftMargin(0.4);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum( SIZE + offsetUp1/lineHeight   );
    h_dummy.SetMinimum(      - offsetDown1/lineHeight );
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);
    for(int i_bin=0;i_bin<h_dummy.GetNbinsX()+1;i_bin++){
        h_dummy.SetBinContent(i_bin,-10);
    }

    g1.SetFillColor(kAzure-4);
    g2.SetFillColor(kCyan);
    g1.SetLineColor(g1.GetFillColor());
    g2.SetLineColor(g2.GetFillColor());

    g1a.SetFillColor(kWhite);
    g2a.SetFillColor(kWhite);
    g1a.SetLineColor(kAzure-4);
    g2a.SetLineColor(kCyan);
    g1a.SetFillStyle(0);
    g2a.SetFillStyle(0);
    g1a.SetLineWidth(1);
    g2a.SetLineWidth(1);

    g.SetLineWidth(2);

    g1a.Draw("5 same");
    g2a.Draw("5 same");
    g1.Draw("2 same");
    g2.Draw("2 same");
    g.Draw("p same");

    TLatex systs{};
    systs.SetTextAlign(32);
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs.DrawLatex(xmin-0.1,i+0.5,Names[i].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    h_dummy.GetXaxis()->CenterTitle();
    h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy.GetXaxis()->SetTitleOffset(1.2);

    TGaxis axis_up( -2, SIZE + (offsetUp1)/lineHeight, 2, SIZE + (offsetUp1)/lineHeight, -poimax,poimax, 510, "-" );
    axis_up.SetLabelOffset( 0.01 );
    axis_up.SetLabelSize(   h_dummy.GetXaxis()->GetLabelSize() );
    axis_up.SetLabelFont(   gStyle->GetTextFont() );
    axis_up.Draw();
    axis_up.CenterTitle();
    axis_up.SetTitle(("#Delta"+fRankingPOIName).c_str());
    if(SIZE==20) axis_up.SetTitleOffset(1.5);
    axis_up.SetTitleSize(   h_dummy.GetXaxis()->GetLabelSize() );
    axis_up.SetTitleFont(   gStyle->GetTextFont() );

    TPad pad1("p1","Pad High",0,(newHeight-offsetUp-offsetUp1)/newHeight,0.4,1);
    pad1.Draw();

    pad1.cd();
    TLegend leg1(0.02,0.7,1,1.0,("Pre-fit impact on "+fRankingPOIName+":").c_str());
    leg1.SetFillStyle(0);
    leg1.SetBorderSize(0);
    leg1.SetMargin(0.25);
    leg1.SetNColumns(2);
    leg1.SetTextFont(gStyle->GetTextFont());
    leg1.SetTextSize(gStyle->GetTextSize());
    leg1.AddEntry(&g1a,"#theta = #hat{#theta}+#Delta#theta","f");
    leg1.AddEntry(&g2a,"#theta = #hat{#theta}-#Delta#theta","f");
    leg1.Draw();

    TLegend leg2(0.02,0.32,1,0.62,("Post-fit impact on "+fRankingPOIName+":").c_str());
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetMargin(0.25);
    leg2.SetNColumns(2);
    leg2.SetTextFont(gStyle->GetTextFont());
    leg2.SetTextSize(gStyle->GetTextSize());
    leg2.AddEntry(&g1,"#theta = #hat{#theta}+#Delta#hat{#theta}","f");
    leg2.AddEntry(&g2,"#theta = #hat{#theta}-#Delta#hat{#theta}","f");
    leg2.Draw();

    TLegend leg0(0.02,0.1,1,0.25);
    leg0.SetFillStyle(0);
    leg0.SetBorderSize(0);
    leg0.SetMargin(0.2);
    leg0.SetTextFont(gStyle->GetTextFont());
    leg0.SetTextSize(gStyle->GetTextSize());
    leg0.AddEntry(&g,"Nuis. Param. Pull","lp");
    leg0.Draw();

    c.cd();

    TLine l0(0,- offsetDown1/lineHeight,0,SIZE+0.5);// + offsetUp1/lineHeight);
    l0.SetLineStyle(kDashed);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    TLine l1 (-1,- offsetDown1/lineHeight,-1,SIZE+0.5);// + offsetUp1/lineHeight);
    l1.SetLineStyle(kDashed);
    l1.SetLineColor(kBlack);
    l1.Draw("same");
    TLine l2(1,- offsetDown1/lineHeight,1,SIZE+0.5);// + offsetUp1/lineHeight);
    l2.SetLineStyle(kDashed);
    l2.SetLineColor(kBlack);
    l2.Draw("same");

    if (fAtlasLabel!= "none") ATLASLabelNew(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.6*offsetUp1)/newHeight), fAtlasLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(       0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.3*offsetUp1)/newHeight), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    gPad->RedrawAxis();

    if(flagGammas && flagSysts){
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c.SaveAs( (fName+"/Ranking"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else if(flagGammas){
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c.SaveAs( (fName+"/RankingGammas"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else if(flagSysts){
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c.SaveAs( (fName+"/RankingSysts"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else{
        WriteWarningStatus("TRExFit::PlotNPRanking", "Your ranking plot felt in unknown category :s");
      for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        c.SaveAs( (fName+"/RankingUnknown"+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format]).c_str() );
    }
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
    for(auto syst : fSystematics){
        if(syst->fName!=syst->fNuisanceParameter){
            // if so, loop on other systematics to find one with name = nuisance parameter
            for(auto syst1 : fSystematics){
                if(syst->fName==syst1->fName) continue;
                if(!(syst->fNuisanceParameter==syst1->fName && syst1->fNuisanceParameter==syst1->fName)) continue;
                // now merge all SystematicHist in all regions
                WriteDebugStatus("TRExFit::MergeSystematics", "Found NP(syst) " + syst->fNuisanceParameter + "(" + syst->fName + ") = to syst name " + syst1->fName );
                for(auto reg : fRegions){
                    WriteDebugStatus("TRExFit::MergeSystematics", "Region: " + reg->fName);
                    for(auto sh : reg->fSampleHists){
                        SystematicHist *syh  = sh->GetSystematic(syst ->fName);
                        SystematicHist *syh1 = sh->GetSystematic(syst1->fName);
                        if(syh!=nullptr && syh1!=nullptr){
                            // FIXME...
                            // the issue here is that to combine uncertainties one has to act differently depending on the fact that the different sources come from a multiplication/division or not...
                            syh1 ->Add(syh);
                            syh1 ->Add(sh->fHist,-1);
                            WriteDebugStatus("TRExFit::MergeSystematics", "Adding syst of " + syh->fName +  " to " + syh1->fName);
                            WriteDebugStatus("TRExFit::MergeSystematics", "Setting to 0 all Up/Down of " +  syh->fName);
                            //
                            // set to zero the other syst
                            syh->fHistUp   = (TH1*)sh->fHist->Clone(syh->fHistUp  ->GetName());
                            syh->fHistDown = (TH1*)sh->fHist->Clone(syh->fHistDown->GetName());
                            syh->fHistShapeUp   = (TH1*)sh->fHist->Clone(syh->fHistShapeUp  ->GetName());
                            syh->fHistShapeDown = (TH1*)sh->fHist->Clone(syh->fHistShapeDown->GetName());
                            syh->fNormUp   = 0.;
                            syh->fNormDown = 0.;
                            syh->fNormPruned  = true;
                            syh->fShapePruned = true;
                        }
                    }
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
    std::vector<std::string> empty;
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
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[regIter]->fName)<0 ) continue;
            //
            fullSelection = FullSelection(  fRegions[regIter],fSamples[i_smp]);
            fullMCweight  = FullWeight(     fRegions[regIter],fSamples[i_smp]);
            fullPaths     = FullNtuplePaths(fRegions[regIter],fSamples[i_smp]);
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                int tmp_debugLevel=TRExFitter::DEBUGLEVEL;
                TRExFitter::SetDebugLevel(0);
                TH1D* htmp = nullptr;
                htmp = HistFromNtuple( fullPaths[i_path],
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
            fullPaths     = FullHistogramPaths(fRegions[regIter],fSamples[i_smp]);
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                int tmp_debugLevel=TRExFitter::DEBUGLEVEL;
                TRExFitter::SetDebugLevel(0);
                std::unique_ptr<TH1> htmp = HistFromFile( fullPaths[i_path] );
                if (!htmp) {
                    WriteErrorStatus("TRExFit::ReadHistograms", "Histo pointer is empty cannot continue running the code");
                    exit(EXIT_FAILURE);
                }
                TRExFitter::SetDebugLevel(tmp_debugLevel);
                //
                // Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[regIter]->fHistoBins){
                    const char *hname = htmp->GetName();
                    std::unique_ptr<TH1> tmp_copy(static_cast<TH1*>(htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin, "tmp_copy", fRegions[regIter]->fHistoBins)));
                    htmp.reset(tmp_copy.release());
                    htmp->SetName(hname);
                    if(TRExFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(htmp.get());
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
                sumSigL += nSigBin * log(1 + nSigBin / nBkgBin);
            }
            err2Bkg += pow(hbkg -> GetBinError(iBin), 2);
            //
            double err2RelBkg = 1;
            if (sumBkg != 0) {
                err2RelBkg = err2Bkg / pow(sumBkg, 2);
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
                pass = sqrt(err2Rel) < 1;
                // distance
                dist = fabs(err2Rel - 1);
            }
            else if(fRegions[regIter]->fBinTransfo == "TransfoF"){
                // "trafo F" with 5% bkg stat unc
                if (sumBkg != 0 && sumSigL != 0)
                  err2Rel = 1 / (sqrt(sumBkg / (nBkg / fRegions[regIter]->fTransfoFzBkg)) + sqrt(sumSigL / (1 / fRegions[regIter]->fTransfoFzSig)));
                else if (sumBkg != 0)
                  err2Rel = sqrt((nBkg / fRegions[regIter]->fTransfoFzBkg) / sumBkg);
                else if (sumSigL != 0)
                  err2Rel = sqrt((1 / fRegions[regIter]->fTransfoFzSig) / sumSigL);
                pass = sqrt(err2Rel) < 1 && sqrt(err2RelBkg) < 0.10;
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
void TRExFit::GetLikelihoodScan( RooWorkspace *ws, std::string varName, RooDataSet* data) const{
    WriteInfoStatus("TRExFit::GetLikelihoodScan", "Running likelihood scan for the parameter = " + varName);

    // shut-up RooFit!
    if(TRExFitter::DEBUGLEVEL<=1){
        if(TRExFitter::DEBUGLEVEL<=0) gErrorIgnoreLevel = kError;
        else if(TRExFitter::DEBUGLEVEL<=1) gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooMsgService::instance().getStream(1).removeTopic(Generation);
        RooMsgService::instance().getStream(1).removeTopic(Plotting);
        RooMsgService::instance().getStream(1).removeTopic(LinkStateMgmt);
        RooMsgService::instance().getStream(1).removeTopic(Eval);
        RooMsgService::instance().getStream(1).removeTopic(Caching);
        RooMsgService::instance().getStream(1).removeTopic(Optimization);
        RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
        RooMsgService::instance().getStream(1).removeTopic(InputArguments);
        RooMsgService::instance().getStream(1).removeTopic(Tracing);
        RooMsgService::instance().getStream(1).removeTopic(Contents);
        RooMsgService::instance().getStream(1).removeTopic(DataHandling);
        RooMsgService::instance().setStreamStatus(1,false);
    }

    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());

    bool isPoI = false;
    RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
    TString firstPOIname = (TString)firstPOI->GetName();
    if (firstPOIname.Contains(varName.c_str())) isPoI = true;

    RooRealVar* var = NULL;
    TString vname = "";
    std::string vname_s = "";
    bool foundSyst = false;
    Double_t minVal = -3;
    Double_t maxVal =  3;
    for(auto nf : fNormFactors){
        if(nf->fName == varName){
            minVal = nf->fMin;
            maxVal = nf->fMax;
        }
    }

    if (fLHscanMin < 99999) { // is actually set
        minVal = fLHscanMin;
    }

    if (fLHscanMax > -99999) { // is actually set
        maxVal = fLHscanMax;
    }

    if (isPoI){
        TIterator* it = mc->GetParametersOfInterest()->createIterator();
        while( (var = (RooRealVar*) it->Next()) ){
            vname=var->GetName();
            vname_s=var->GetName();
            if (vname == varName || vname == "alpha_"+varName) {
                WriteInfoStatus("TRExFit::GetLikelihoodScan", "GetLikelihoodScan for POI = " + vname_s);
                foundSyst=true;
                break;
            }
        }
    }
    else {
        TIterator* it = mc->GetNuisanceParameters()->createIterator();
        while( (var = (RooRealVar*) it->Next()) ){
        vname=var->GetName();
            vname_s=var->GetName();
            if (vname == varName || vname == "alpha_"+varName) {
                WriteInfoStatus("TRExFit::GetLikelihoodScan", "GetLikelihoodScan for NP = " + vname_s);
                foundSyst=true;
                break;
            }
        }
    }

    if(!foundSyst){
        WriteWarningStatus("TRExFit::GetLikelihoodScan", "systematic " + varName + " not found (most probably due to Pruning), skip LHscan !");
        return;
    }
    WriteInfoStatus("TRExFit::GetLikelihoodScan", "GetLikelihoodScan for parameter = " + vname_s);

    TCanvas can("NLLscan");

    RooAbsReal* nll = simPdf->createNLL(*data,
                                        Constrain(*mc->GetNuisanceParameters()),
                                        Offset(1),
                                        NumCPU(TRExFitter::NCPU, RooFit::Hybrid),
                                        RooFit::Optimize(kTRUE));

    std::vector<double> x(fLHscanSteps);
    std::vector<double> y(fLHscanSteps);
    RooMinimizer m(*nll); // get MINUIT interface of fit
    m.setErrorLevel(-1);
    m.setPrintLevel(-1);
    m.setStrategy(2); // set precision to high
    var->setConstant(kTRUE); // make POI constant in the fit
    double min = 9999999;
    for (int ipoint = 0; ipoint < fLHscanSteps; ++ipoint) {
        WriteInfoStatus("TRExFit::GetLikelihoodScan","Running LHscan for point " + std::to_string(ipoint+1) + " out of " + std::to_string(fLHscanSteps) + " points");
        x[ipoint] = minVal+ipoint*(maxVal-minVal)/(fLHscanSteps - 1);
        *var = x[ipoint]; // set POI
        m.migrad(); // minimize again with new posSigXsecOverSM value
        RooFitResult* r = m.save(); // save fit result
        y[ipoint] = r->minNll();
        if (y[ipoint] < min) min = y[ipoint];
    }
    var->setConstant(kFALSE); // make POI not constant otherwise we run in errors in the 2D scan

    for (auto & iY : y) {
        iY = iY - min;
    }

    TGraph graph(fLHscanSteps, &x[0], &y[0]);
    graph.Draw("ALP");
    graph.GetXaxis()->SetRangeUser(minVal,maxVal);

    // y axis
    graph.GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());

    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(vname);

    can.SetTitle(cname);
    can.SetName(cname);
    can.cd();

    TLatex tex{};
    tex.SetTextColor(kGray+2);

    TLine l1s(minVal,0.5,maxVal,0.5);
    l1s.SetLineStyle(kDashed);
    l1s.SetLineColor(kGray);
    l1s.SetLineWidth(2);
    if(graph.GetMaximum()>2){
        l1s.Draw();
        tex.DrawLatex(maxVal,0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }

    if(isPoI){
        if(graph.GetMaximum()>2){
            TLine l2s(minVal,2,maxVal,2);
            l2s.SetLineStyle(kDashed);
            l2s.SetLineColor(kGray);
            l2s.SetLineWidth(2);
            l2s.Draw();
            tex.DrawLatex(maxVal,2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
        }
        //
        if(graph.GetMaximum()>4.5){
            TLine l3s(minVal,4.5,maxVal,4.5);
            l3s.SetLineStyle(kDashed);
            l3s.SetLineColor(kGray);
            l3s.SetLineWidth(2);
            l3s.Draw();
            tex.DrawLatex(maxVal,4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
        }
        //
        TLine lv0(0,graph.GetMinimum(),0,graph.GetMaximum());
        lv0.Draw();
        //
        TLine lh0(minVal,0,maxVal,0);
        lh0.Draw();
    }

    TString LHDir("LHoodPlots/");
    system(TString("mkdir -vp ")+fName+"/"+LHDir);

    can.RedrawAxis();

    for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++)
        can.SaveAs( fName+"/"+LHDir+"NLLscan_"+varName+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format] );

    // write it to a ROOT file as well
    TFile *f = new TFile(fName+"/"+LHDir+"NLLscan_"+varName+fSuffix+"_curve.root","UPDATE");
    f->cd();
    graph.Write("LHscan",TObject::kOverwrite);
    f->Close();
    delete f;
}

//____________________________________________________________________________________
//
void TRExFit::Get2DLikelihoodScan( RooWorkspace *ws, const std::vector<std::string>& varNames, RooDataSet* data) const{
    if (varNames.size() != 2){
        WriteErrorStatus("TRExFit::Get2DLikelihoodScan", "Wrong number of parameters provided for 2D likelihood scan, returning");
        return;
    }
    WriteInfoStatus("TRExFit::Get2DLikelihoodScan", "Running 2D likelihood scan for the parameters = " + varNames.at(0) + " and " + varNames.at(1));

    // shut-up RooFit!
    if(TRExFitter::DEBUGLEVEL<=1){
        if(TRExFitter::DEBUGLEVEL<=0) gErrorIgnoreLevel = kError;
        else if(TRExFitter::DEBUGLEVEL<=1) gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooMsgService::instance().getStream(1).removeTopic(Generation);
        RooMsgService::instance().getStream(1).removeTopic(Plotting);
        RooMsgService::instance().getStream(1).removeTopic(LinkStateMgmt);
        RooMsgService::instance().getStream(1).removeTopic(Eval);
        RooMsgService::instance().getStream(1).removeTopic(Caching);
        RooMsgService::instance().getStream(1).removeTopic(Optimization);
        RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
        RooMsgService::instance().getStream(1).removeTopic(InputArguments);
        RooMsgService::instance().getStream(1).removeTopic(Tracing);
        RooMsgService::instance().getStream(1).removeTopic(Contents);
        RooMsgService::instance().getStream(1).removeTopic(DataHandling);
    }

    RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());



    //Vector for the two parameters
    RooRealVar* varX = nullptr;
    RooRealVar* varY = nullptr;
    //Get the parameters from the model
    TIterator* it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar* var_tmp = nullptr;
    TString vname = "";
    int count = 0;

    // iterate over NPs
    while ( (var_tmp = static_cast<RooRealVar*>(it->Next())) ){
        vname=var_tmp->GetName();
        if (vname == varNames.at(0) || vname == "alpha_"+varNames.at(0)){
            varX = var_tmp;
            count++;
        }
        if (vname == varNames.at(1) || vname == "alpha_"+varNames.at(1)){
            varY = var_tmp;
            count++;
        }
        if (count == 2) break;
    }

    // iterate over POIs
    if (count < 2){
        TIterator* it_POI = mc->GetParametersOfInterest()->createIterator();
        while ( (var_tmp = static_cast<RooRealVar*>(it_POI->Next())) ){
            vname=var_tmp->GetName();
            if (vname == varNames.at(0) || vname == "alpha_"+varNames.at(0)){
                varX = var_tmp;
                count++;
            }
            if (vname == varNames.at(1) || vname == "alpha_"+varNames.at(1)){
                varY = var_tmp;
                count++;
            }
            if (count == 2) break;
        }
    }
    if (count != 2) {
        WriteErrorStatus("TRExFit::Get2DLikelihoodScan","Did not find the two parameters you want to use in the 2D likelihood scan");
        return;
    }
    WriteInfoStatus("TRExFit::Get2DLikelihoodScan", "Setting up the NLL");

    //To set the boundaries
    Double_t minValX = varX->getMin();
    Double_t maxValX = varX->getMax();
    Double_t minValY = varY->getMin();
    Double_t maxValY = varY->getMax();

    if (fLHscanMin < 99999) { // is actually set
        minValX = fLHscanMin;
    }
    if (fLHscanMinY < 99999) { // is actually set
        minValY = fLHscanMinY;
    }
    if (fLHscanMax > -99999) { // is actually set
        maxValX = fLHscanMax;
    }
    if (fLHscanMaxY > -99999) { // is actually set
        maxValY = fLHscanMaxY;
    }

    unsigned int offset = 1;
    if (fParal2D) {
        // When we run in parrellel we cant set offset to 1
        // this caused problems with the offset between the different sup processes
        offset = 0;
    }
    RooAbsReal* nll = simPdf->createNLL(*data,
                                        Constrain(*mc->GetNuisanceParameters()),
                                        Offset(offset),
                                        NumCPU(TRExFitter::NCPU, RooFit::Hybrid),
                                        RooFit::Optimize(kTRUE));

    RooMinimizer m(*nll); // get MINUIT interface of fit
    m.setErrorLevel(-1);
    m.setPrintLevel(-1);
    m.setStrategy(2); // set precision to high
    //Set both POIs to constant
    varX->setConstant(kTRUE); // make POI constant in the fit
    varY->setConstant(kTRUE); // make POI constant in the fit

    //values for parameter1, parameter2 and the NLL value
    std::vector<double> x(fLHscanSteps);
    std::vector<double> y(fLHscanStepsY);
    std::vector<std::vector<double>> z(fLHscanSteps, std::vector<double>(fLHscanStepsY));

    double zmin = 9999999;

    //Actual scan
    WriteInfoStatus("TRExFit::Get2DLikelihoodScan", "Start of the 2D scan");
    for (int ipoint = 0; ipoint < fLHscanSteps; ++ipoint) {
        if (fParal2D && ipoint!=fParal2Dstep) // if you are parallelizing, only run the point corresponding to the one passed from command line
            continue;
        WriteInfoStatus("TRExFit::Get2DLikelihoodScan","Running LHscan for point " + std::to_string(ipoint+1) + " out of " + std::to_string(fLHscanSteps) + " points");
        // x[ipoint] = minValX + ipoint * (maxValX - minValX) / (fLHscanSteps);
        // We could alternatively use the line below to inlcude the max value in the scan
        x[ipoint] = minValX + ipoint * (maxValX - minValX) / (fLHscanSteps - 1);
        *(varX) = x[ipoint]; // set POI
        for (int jpoint = 0; jpoint < fLHscanStepsY; ++jpoint) {
            WriteInfoStatus("TRExFit::Get2DLikelihoodScan","Running LHscan for subpoint " + std::to_string(jpoint+1) + " out of " + std::to_string(fLHscanStepsY) + " points");
            // y[jpoint] = minValY + jpoint * (maxValY - minValY) / (fLHscanStepsY);
            // We could alternatively use the line below to inlcude the max value in the scan
            y[jpoint] = minValY + jpoint * (maxValY - minValY) / (fLHscanStepsY - 1);
            *(varY) = y[jpoint]; // set POI
            m.migrad(); // minimize again with new posSigXsecOverSM value
            RooFitResult* r = m.save(); // save fit result
            const double z_tmp = r->minNll();
            z[ipoint][jpoint] = z_tmp;

            // save the best values
            if (z_tmp < zmin) {
                zmin = z_tmp;
            }
        }
    }
    //Set both POIs not constant
    varX->setConstant(kTRUE); // make POI not constant after the fit
    varY->setConstant(kTRUE); // make POI not constant after the fit

    // end of scaning, now fill some plots


    // this is needed for potential blinding
    TRandom3 rand{};
    rand.SetSeed(1234567);
    const double rndNumber = rand.Uniform(5);
    bool blindVarX = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.at(0)) != fBlindedParameters.end();
    bool blindVarY = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.at(1)) != fBlindedParameters.end();
    if (blindVarX){
        minValX += rndNumber;
        maxValX += rndNumber;
        for (auto & iX : x) {
            iX+= rndNumber;
        }
    }
    if (blindVarY){
        minValY += rndNumber;
        maxValY += rndNumber;
        for (auto & iY : y) {
            iY+= rndNumber;
        }
    }

    // make plots
    TCanvas can("NLLscan_2D_");
    can.cd();

    TGraph2D graph(fLHscanSteps * fLHscanStepsY);

    TH2D h_nll("NLL", "NLL", fLHscanSteps, minValX, maxValX, fLHscanStepsY, minValY, maxValY);
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
        graph.GetXaxis()->SetRangeUser(minValX,maxValX);
        graph.GetYaxis()->SetRangeUser(minValY,maxValY);

        // y axis
        graph.GetXaxis()->SetTitle(varNames.at(0).c_str());
        graph.GetYaxis()->SetTitle(varNames.at(1).c_str());

        // Print the canvas
        for(int i_format=0;i_format<(int)TRExFitter::IMAGEFORMAT.size();i_format++){
            can.SaveAs( fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"."+TRExFitter::IMAGEFORMAT[i_format] );
        }

        // write it to a ROOT file as well
        std::unique_ptr<TFile> f = std::make_unique<TFile>(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"_curve.root","UPDATE");
        f->cd();
        graph.Write(("LHscan_2D_"+varNames.at(0)+"_"+varNames.at(1)).c_str(),TObject::kOverwrite);
        f->Close();
    }

    // Write histogram to Root file as well
    if (fParal2D) {
        std::ostringstream step_os;
        step_os << fParal2Dstep;
        std::string paral2Dstep_str=step_os.str();
        std::unique_ptr<TFile> f2 = std::make_unique<TFile>(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"_step"+paral2Dstep_str+fSuffix+"_histo.root","UPDATE");
        h_nll.Write("NLL",TObject::kOverwrite);
        f2->Close();
    } else {
        std::unique_ptr<TFile> f2 = std::make_unique<TFile>(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"_histo.root","UPDATE");
        h_nll.Write("NLL",TObject::kOverwrite);
        f2->Close();
    }
}

//____________________________________________________________________________________
//
void TRExFit::DefineVariable(int regIter){
    TH1::StatOverflows(true);  //////  What is the defaut in root for this ???
    WriteDebugStatus("TRExFit::DefineVariable", "//////// --------");
    WriteDebugStatus("TRExFit::DefineVariable", "// DEBUG CORR VAR");
    TH1* h1 = new TH1D("h1","h1",1,-2000.,1000.);
    TH1* h2 = new TH1D("h2","h2",1,-2000.,1000.);
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;
    std::vector<std::string> empty;

    // copy of NtupleReading function.
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        WriteDebugStatus("TRExFit::DefineVariable", "Processing sample : " + fSamples[i_smp]->fName);
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[regIter]->fName)<0 ) continue;
        WriteDebugStatus("TRExFit::DefineVariable", " -> is used in the considered region");
        //
        // set selection, weight and paths (no variables)
        fullSelection = FullSelection(  fRegions[regIter],fSamples[i_smp]);
        fullMCweight  = FullWeight(     fRegions[regIter],fSamples[i_smp]);
        fullPaths     = FullNtuplePaths(fRegions[regIter],fSamples[i_smp]);
        //
        for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
            WriteDebugStatus("TRExFit::DefineVariable", " -> Retrieving : " + fRegions[regIter]->fCorrVar1 +
                                                        " w/ weight " + fullMCweight + "*" + fullSelection  +
                                                        " from " +  fullPaths[i_path]);
            TH1* htmp1 = new TH1D("htmp1","htmp1",1,-2000.,1000.);
            TH1* htmp2 = new TH1D("htmp2","htmp2",1,-2000.,1000.);
            TChain *t = new TChain();
            t->Add(fullPaths[i_path].c_str());
            t->Draw( Form("%s>>htmp1",fRegions[regIter]->fCorrVar1.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            t->Draw( Form("%s>>htmp2",fRegions[regIter]->fCorrVar2.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            delete t;
            //
            if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp1 -> Scale(fLumi);
            if(fSamples[i_smp]->fLumiScales.size()>i_path)  htmp1 -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
            else if(fSamples[i_smp]->fLumiScales.size()==1) htmp1 -> Scale(fSamples[i_smp]->fLumiScales[0]);
            //
            if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp2 -> Scale(fLumi);
            if(fSamples[i_smp]->fLumiScales.size()>i_path)  htmp2 -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
            else if(fSamples[i_smp]->fLumiScales.size()==1) htmp2 -> Scale(fSamples[i_smp]->fLumiScales[0]);
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

    WriteDebugStatus("TRExFit::DefineVariable", "the new variable : ( ( (" + fRegions[regIter]->fCorrVar1 + ") - " + std::to_string(mean1) + " )*( (" +fRegions[regIter]->fCorrVar2 + ")-" + std::to_string(mean2) + " ) )/( " + std::to_string(rms1) + " * " + std::to_string(rms2) + ")");

    fRegions[regIter]->fVariable = Form("( ( (%s)-%f )*( (%s)-%f ) )/( %f * %f )",fRegions[regIter]->fCorrVar1.c_str(),mean1,fRegions[regIter]->fCorrVar2.c_str(),mean2,rms1,rms2);

    TH1::StatOverflows(false);  //////  What is the defaut in root for this ???

    delete h1;
    delete h2;
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
    double deltaXp = -1.; // |x(i+1)-x(i)|
    double deltaXm = -1.; // |x(i-1)-x(i)|
    std::string name;
    if (itemp < templatePair.size()){
        x_i = templatePair.at(itemp).first;
        name = templatePair.at(itemp).second;
    }
    else return fun;
    //
    if (opt == TRExFit::LINEAR){
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
    fun = ReplaceString(fun,"--","+");
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
    double k_init(80.);
    double k_left(k_init);
    double k_right(k_init);

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
    k_left = k_init/width_left;
    k_right = k_init/width_right;

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
    NormFactor *nf = nullptr;
    for(auto norm : fNormFactors){
        if(norm->fName == name) nf = norm;
    }
    // get one histogram per bin (per region)
    for(auto reg : fRegions){
        std::map<double,TH1*> hMap; // map (paramater-value,histogram)
        TH1* h_tmp = nullptr;
        int nTemplates = 0;
        double min = -999.;
        double max = -999.;
        for(auto sh : reg->fSampleHists){
            Sample* smp = sh->fSample;
            // if the sample has morphing
            if(smp->fIsMorph[name]){
                hMap[smp->fMorphValue[name]] = sh->fHist;
                if(smp->fMorphValue[name]<min) min = smp->fMorphValue[name];
                if(smp->fMorphValue[name]>max) max = smp->fMorphValue[name];
                nTemplates++;
                h_tmp = sh->fHist;
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
            for(auto format : TRExFitter::IMAGEFORMAT) c.SaveAs((fName+"/Morphing/g_"+name+"_"+reg->fName+"_bin"+std::to_string(i_bin)+"."+format).c_str());
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
   for(int i_syst=0;i_syst<fNSyst;i_syst++){
       if(fSystematics[i_syst]->fSubCategory=="Gammas" || fSystematics[i_syst]->fSubCategory=="FullSyst" || fSystematics[i_syst]->fSubCategory=="combine")
            WriteWarningStatus("TRExFit::ProduceSystSubCategoryMap"," use of \"Gammas\", \"FullSyst\" or \"combine\" as SubCategory names is not supported, you will likely run into issues");
       fSubCategoryImpactMap.insert(std::make_pair(("alpha_" + fSystematics[i_syst]->fNuisanceParameter).c_str(), fSystematics[i_syst]->fSubCategory));
   }

   // also add norm factors, no "alpha_" needed
   for(int i_nf=0;i_nf<fNNorm;i_nf++){
       if(fNormFactors[i_nf]->fSubCategory=="Gammas" || fNormFactors[i_nf]->fSubCategory=="FullSyst" || fNormFactors[i_nf]->fSubCategory=="combine")
            WriteWarningStatus("TRExFit::ProduceSystSubCategoryMap"," use of \"Gammas\", \"FullSyst\" or \"combine\" as SubCategory names is not supported, you will likely run into issues");
       if (fNormFactors[i_nf]->fName != fPOI) {
           fSubCategoryImpactMap.insert(std::make_pair(fNormFactors[i_nf]->fNuisanceParameter, fNormFactors[i_nf]->fSubCategory));
       }
   }
}

//____________________________________________________________________________________
// combine individual results from grouped impact evaluation into one table
void TRExFit::BuildGroupedImpactTable() const{
    WriteInfoStatus("TRExFit::BuildGroupedImpactTable", "merging grouped impact evaluations");
    std::string targetName = fName+"/Fits/GroupedImpact"+fSuffix+".txt";

    if(std::ifstream(targetName).good()){
        WriteWarningStatus("TRExFit::BuildGroupedImpactTable","file " + targetName + " already exists, will not overwrite");
    }
    else{
        std::string cmd = " if [[ `ls "+fName+"/Fits/GroupedImpact"+fSuffix+"_*` != \"\" ]] ; then";
        cmd            += " cat "+fName+"/Fits/GroupedImpact"+fSuffix+"_* > "+targetName+" ; ";
        cmd            += " fi ;";
        gSystem->Exec(cmd.c_str());
    }
}

//____________________________________________________________________________________
//
void TRExFit::RunToys(RooWorkspace* ws){
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
        for(auto reg : fRegions) regionDataType[reg->fName] = Region::ASIMOVDATA;
        for( int i_ch = 0; i_ch < fNRegions; i_ch++ ){
            if ( fFitRegion == CRONLY && fRegions[i_ch] -> fRegionType == Region::CONTROL )
                regionsToFit.push_back( fRegions[i_ch] -> fName );
            else if ( fFitRegion == CRSR && (fRegions[i_ch] -> fRegionType == Region::CONTROL || fRegions[i_ch] -> fRegionType == Region::SIGNAL) )
                regionsToFit.push_back( fRegions[i_ch] -> fName );
        }
        ws = PerformWorkspaceCombination( regionsToFit );
        //Setting binned likelihood option
        RooFIter rfiter = ws->components().fwdIterator();
        RooAbsArg* arg;
        while ((arg = rfiter.next())) {
            if (arg->IsA() == RooRealSumPdf::Class()) {
                arg->setAttribute("BinnedLikelihood");
                std::string temp_string = arg->GetName();
                WriteDebugStatus("TRExFit::DumpData", "Activating binned likelihood attribute for " + temp_string);
            }
        }
        if (!ws){
            WriteErrorStatus("TRExFit::RunToys","Cannot retrieve the workspace, exiting!");
            exit(EXIT_FAILURE);
        }
        // create map to store fit results
        // create histogram to store fitted POI values
        NormFactor POInf = *(fNormFactors[FindInStringVector(fNormFactorNames,fPOI)]);
        double min = POInf.fMin;
        double max = POInf.fMax;
        if (fToysHistoMin < 9900 && fToysHistoMax > -9000){
            min = fToysHistoMin;
            max = fToysHistoMax;
        }
        TH1D h_toys ("h_toys","h_toys",fToysHistoNbins,min,max);
        // get RooStats stuff
        RooStats::ModelConfig mc = *((RooStats::ModelConfig*)ws -> obj("ModelConfig"));
        RooSimultaneous simPdf = *((RooSimultaneous*)(mc.GetPdf()));
        RooAbsPdf *pdf = mc.GetPdf();
        const RooArgSet obsSet = *(mc.GetObservables());
        RooRealVar* poiVar = (RooRealVar*) (& ws->allVars()[fPOI.c_str()]);

        //Get the desired NP
        RooRealVar* NPtoShift = nullptr;
        if (fToysPseudodataNP != "") {
            NPtoShift = (RooRealVar*) (& ws->allVars()[fToysPseudodataNP.c_str()]);
        }

        //For the loop over NPs
        std::string varname{};

        //Create NLL only once
        RooDataSet *dummy = pdf->generate( obsSet, RooFit::Extended() );
        RooAbsReal* nll = simPdf.createNLL(*dummy,
                                           Constrain(*mc.GetNuisanceParameters()),
                                           Offset(1),
                                           NumCPU(1, RooFit::Hybrid),
                                           RooFit::Optimize(kTRUE));


        for(int i_toy=0;i_toy<fFitToys;i_toy++){

            if (fToysPseudodataNP != "") {
                TIterator* it = mc.GetNuisanceParameters()->createIterator();
                RooRealVar* var = nullptr;
                while( (var = (RooRealVar*) it->Next()) ){
                    varname = var->GetName();
                    if (varname.find("alpha_")!=std::string::npos) {
                        var->setConstant(1);
                        var->setVal(0);
                    } else {
                        var->setConstant(1);
                        var->setVal(1);
                    }
                }
                delete it;
                delete var;
            }

            // setting POI to constant, not to allow it to fluctuate in toy creation
            poiVar->setConstant(1);
            poiVar->setVal(fFitPOIAsimov);
            if (fToysPseudodataNP != "") {
                NPtoShift->setConstant(1);
                NPtoShift->setVal(fToysPseudodataNPShift);
            }

            WriteInfoStatus("TRExFit::RunToys","Generating toy n. " + std::to_string(i_toy+1) + " out of " + std::to_string(fFitToys) + " toys");
            RooDataSet *toyData = pdf->generate( obsSet, RooFit::Extended() );
            // re-set POI to free-floating, and to nominal value
            if (fToysPseudodataNP != "") {
                RooArgSet* nuis = (RooArgSet*) mc.GetNuisanceParameters();
                if (nuis){
                    RooRealVar* vartmp = nullptr;
                    TIterator* it2 = nuis->createIterator();
                    while( (vartmp = (RooRealVar*) it2->Next()) ){
                        std::string np = vartmp->GetName();
                        if (np.find("alpha_")!=std::string::npos) {
                            vartmp->setConstant(0);
                            vartmp->setVal(0);
                        }
                        else if( np.find("gamma_")!=std::string::npos ){
                            vartmp->setVal(1);
                            vartmp->setConstant(0);
                        }
                        else {  // for norm factors
                            vartmp->setVal( 1 );
                        }
                    }
                    delete it2;
                    delete vartmp;
                }
            }
            poiVar->setConstant(0);
            poiVar->setVal(POInf.fNominal);

            // NP is fixed constant for each fit, and to nominal value
            if (fToysPseudodataNP != "") {
                NPtoShift->setConstant(1);
                NPtoShift->setVal(0);
            }
            // extract POI from fit result and fill histogram
            WriteInfoStatus("TRExFit::RunToys","Fitting toy n. " + std::to_string(i_toy+1));
            //Set new dataset for NLL
            nll->setData(*toyData);
            RooMinimizer m(*nll); // get MINUIT interface of fit
            m.setErrorLevel(-1);
            m.setPrintLevel(-1);
            //m.setStrategy(2); // set precision to high
            m.migrad();
            RooFitResult* r = m.save(); // save fit result

            h_toys.Fill(poiVar->getVal());
            WriteInfoStatus("TRExFit::RunToys","Toy n. " + std::to_string(i_toy+1) + ", fitted value: " + std::to_string(poiVar->getVal()));

            delete r;
            delete toyData;
        }
        // plot, fit and save toy histogram
        TCanvas c("c","c",600,600);
        h_toys.Draw("E");
        TF1 g("g","gaus",POInf.fMin,POInf.fMax);
        g.SetLineColor(kRed);
        h_toys.Fit("g","RQ");
        g.Draw("same");
        h_toys.GetXaxis()->SetTitle(TRExFitter::SYSTMAP[fPOI].c_str());
        h_toys.GetYaxis()->SetTitle("Pseudo-experiements");
        myText(0.60,0.90,1,Form("Mean  = %.2f #pm %.2f",g.GetParameter(1), g.GetParError(1)));
        myText(0.60,0.85,1,Form("Sigma = %.2f #pm %.2f",g.GetParameter(2),g.GetParError(2)));
        myText(0.60,0.80,1,Form("#chi^{2}/ndf = %.2f / %d",g.GetChisquare(),g.GetNDF()));
        for(auto format : TRExFitter::IMAGEFORMAT) c.SaveAs((fName+"/Toys/ToysPlot."+format).c_str());
        fVarNameMinos = varMinosTmp; // retore Minos settings

        // Also create a ROOT file
        TFile *out = new TFile ((fName+"/Toys/Toys"+fSuffix+".root").c_str(), "RECREATE");
        out->cd();
        h_toys.Write();
        out->Close();
        delete out;
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
    if(!CheckExpression(variable)){
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
        selection = ReplaceString(selection,smp->fIgnoreSelection,"1");
    }
    // from Sample
    if(smp->fSelection!="" && smp->fSelection!="1"){
        if(selection!="") selection += " && ";
        selection += "("+smp->fSelection+")";
    }
    // check the final expression
    if(!CheckExpression(selection)){
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
        weight = ReplaceString(weight,smp->fIgnoreWeight,"1");
    }
    // from Sample (nominal...
    std::string sampleWeight = "";
    if(smp->fMCweight!="" && smp->fMCweight!="1"){
        sampleWeight = smp->fMCweight;
    }
    // ... and systematics)
    if(syst!=nullptr){
        if(syst->fIgnoreWeight!=""){
            weight = ReplaceString(weight,syst->fIgnoreWeight,"1");
            sampleWeight = ReplaceString(sampleWeight,syst->fIgnoreWeight,"1");
        }
        if(isUp){
            if(syst->fWeightUp!=""){
                sampleWeight = syst->fWeightUp;
            }
            else if(syst->fWeightSufUp!=""){
                if(sampleWeight!="") sampleWeight += " * ";
                sampleWeight += syst->fWeightSufUp;
            }
        }
        else{
            if(syst->fWeightDown!=""){
                sampleWeight = syst->fWeightDown;
            }
            else if(syst->fWeightSufDown!=""){
                if(sampleWeight!="") sampleWeight += " * ";
                sampleWeight += syst->fWeightSufDown;
            }
        }
    }
    if(sampleWeight!=""){
        if(weight!="") weight += " * ";
        weight += "("+sampleWeight+")";
    }
    // add Bootstrap weights
    if(fBootstrap!="" && fBootstrapIdx>=0){
        if(weight!="") weight += " * ";
        weight += "("+ReplaceString(fBootstrap,"BootstrapIdx",Form("%d",fBootstrapIdx))+")";
        gRandom->SetSeed(fBootstrapIdx);
    }
    // check the final expression
    WriteDebugStatus("TRExFit::FullWeight","Full weight expression : "+weight);
    if(!CheckExpression(weight)){
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
            if(isUp) pathSuffs = CombinePathSufs( CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), ToVec(syst->fNtuplePathSufUpRefSample) );
            else     pathSuffs = CombinePathSufs( CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), ToVec(syst->fNtuplePathSufDownRefSample) );
        }
        else{
            if(isUp) pathSuffs = CombinePathSufs( CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), ToVec(syst->fNtuplePathSufUp) );
            else     pathSuffs = CombinePathSufs( CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), ToVec(syst->fNtuplePathSufDown) );
        }
    }
    else{
        pathSuffs = CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) fileSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), ToVec(syst->fNtupleFileSufUpRefSample) );
            else     fileSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), ToVec(syst->fNtupleFileSufDownRefSample) );
        }
        else{
            if(isUp) fileSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), ToVec(syst->fNtupleFileSufUp) );
            else     fileSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), ToVec(syst->fNtupleFileSufDown) );
        }
    }
    else{
        fileSuffs = CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) nameSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), ToVec(syst->fNtupleNameSufUpRefSample) );
            else     nameSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), ToVec(syst->fNtupleNameSufDownRefSample) );
        }
        else{
            if(isUp) nameSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), ToVec(syst->fNtupleNameSufUp) );
            else     nameSuffs = CombinePathSufs( CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), ToVec(syst->fNtupleNameSufDown) );
        }
    }
    else{
        nameSuffs = CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs );
    }
    //
    // And finally put everything together
    fullPaths = CreatePathsList( paths,pathSuffs, files,fileSuffs, names,nameSuffs );
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
            if(isUp) pathSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), ToVec(syst->fHistoPathSufUpRefSample) );
            else     pathSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), ToVec(syst->fHistoPathSufDownRefSample) );
        }
        else{
            if(isUp) pathSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), ToVec(syst->fHistoPathSufUp) );
            else     pathSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs ), ToVec(syst->fHistoPathSufDown) );
        }
    }
    else{
        pathSuffs = CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) fileSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), ToVec(syst->fHistoFileSufUpRefSample) );
            else     fileSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), ToVec(syst->fHistoFileSufDownRefSample) );
        }
        else{
            if(isUp) fileSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), ToVec(syst->fHistoFileSufUp) );
            else     fileSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs ), ToVec(syst->fHistoFileSufDown) );
        }
    }
    else{
        fileSuffs = CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs );
    }
    //
    if(syst!=nullptr){
        if(isData && syst->fSubtractRefSampleVar && syst->fReferenceSample==smp->fName){
            if(isUp) nameSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), ToVec(syst->fHistoNameSufUpRefSample) );
            else     nameSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), ToVec(syst->fHistoNameSufDownRefSample) );
        }
        else{
            if(isUp) nameSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), ToVec(syst->fHistoNameSufUp) );
            else     nameSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs ), ToVec(syst->fHistoNameSufDown) );
        }
    }
    else{
      nameSuffs = CombinePathSufs( CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs), fHistoNamesNominal );
    }
    //
    // And finally put everything together
    fullPaths = CreatePathsList( paths,pathSuffs, files,fileSuffs, names,nameSuffs );
    return fullPaths;
}

//__________________________________________________________________________________
//
TH1D* TRExFit::ReadSingleHistogram(const std::vector<std::string>& fullPaths, Systematic* syst,
 int i_ch, int i_smp, bool isUp, bool isMC){
    TH1D* h = nullptr;
    for(unsigned int i_path = 0; i_path < fullPaths.size(); ++i_path){
        std::unique_ptr<TH1> htmp = HistFromFile( fullPaths.at(i_path) );
        if (!htmp) {
            WriteErrorStatus("TRExFit::ReadSingleHistogram", "Histo pointer is nullptr, cannot continue running the code");
            exit(EXIT_FAILURE);
        }
        //Pre-processing of histograms (rebinning, lumi scaling)
        if(fRegions[i_ch]->fHistoBins){
            const char *hname = htmp->GetName();
            std::unique_ptr<TH1> tmp_copy(static_cast<TH1D*>(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin, "tmp_copy", fRegions[i_ch]->fHistoBins)));
            htmp.reset(tmp_copy.release());
            htmp->SetName(hname);
            if(TRExFitter::MERGEUNDEROVERFLOW) MergeUnderOverFlow(htmp.get());
        }
        else if(fRegions[i_ch]->fHistoNBinsRebin != -1) {
            htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin);
        }

        if (isMC){
            if(fSamples[i_smp]->fNormalizedByTheory){
                htmp -> Scale(fLumi);
            }
        }

        if(fSamples[i_smp]->fLumiScales.size()>i_path){
             htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
        }
        else if(fSamples[i_smp]->fLumiScales.size()==1){
            htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
        }

        if (isMC && syst != nullptr){
            // obtain relative variation and apply it to proper sample
            // & try to keep also the same total relative variation
            if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                // check if the reference sample exists
                if (fRegions[i_ch]->GetSampleHist(syst->fReferenceSample) == nullptr){
                    WriteErrorStatus("TRExFit::ReadSingleHistogram", "Reference sample: " + syst->fReferenceSample + " does not exist for region: " + fRegions[i_ch]->fName + ". Please check this!");
                    exit(EXIT_FAILURE);
                }
                TH1* href = fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->fHist;
                TH1* hnom = fRegions[i_ch]->GetSampleHist( fSamples[i_smp]->fName )->fHist;
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
                if( syst->fKeepReferenceOverallVar && (TMath::Abs(relVar-1) > 0.0001) &&
                    (TMath::Abs(newVar-1) > 0.0001)){
                    htmp->Scale( relVar / newVar );
                }
            }
            // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
            else if (syst->fReferenceSample!="" && syst->fSubtractRefSampleVar) {
                // check if the reference sample exists
                if (fRegions[i_ch]->GetSampleHist(syst->fReferenceSample) == nullptr){
                    WriteErrorStatus("TRExFit::ReadSingleHistogram", "Reference sample: " + syst->fReferenceSample + " does not exist for region: " + fRegions[i_ch]->fName + ". Please check this!");
                    exit(EXIT_FAILURE);
                }
                TH1* href = fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->fHist;
                TH1* href_upDown = nullptr;
                if (isUp){
                    href_upDown = fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->
                        GetSystematic(syst->fName)->fHistUp;
                } else {
                    href_upDown = fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->
                        GetSystematic(syst->fName)->fHistDown;
                }
                TH1* hnom = fRegions[i_ch]->GetSampleHist( fSamples[i_smp]->fName )->fHist;
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
                h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),
                    fSamples[i_smp]->fName.c_str())));
            } else { // is syst
                if (isUp){ // up variation
                    h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s_%sUp",fRegions[i_ch]->fName.c_str(),
                        fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str())));
                } else { // down variation
                    h = static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s_%sDown",fRegions[i_ch]->fName.c_str(),
                        fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str())));
                }
            }
        }
        else{
            h->Add(htmp.get());
        }
    }

    return h;
}

//__________________________________________________________________________________
//
SampleHist* TRExFit::GetSampleHistFromName(const Region* const reg, const std::string& name) const{
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
    TH1* currentNominal = sh->fHist;

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
    for (int i = 0; i < sh->fNSyst; ++i){
        if (sh->fSyst[i]->fName == name){
            return i;
        }
    }

    return -1;
}
