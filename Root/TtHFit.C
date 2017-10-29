#include <cctype>
#include <iomanip>

//TtHFitter headers
#include "TtHFitter/FittingTool.h"
#include "TtHFitter/HistoTools.h"

//Roofit headers
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooRealSumPdf.h"

//HistFactory headers
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"

//Corresponding header
#include "TtHFitter/TtHFit.h"

#include "TFormula.h"

using namespace RooFit;

// -------------------------------------------------------------------------------------------------
// class TtHFit

//__________________________________________________________________________________
//
TtHFit::TtHFit(string name){
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
    fUseStatErr = false;
    fStatErrThres = 0.05;

    fLumi = 1.;
    fLumiErr = 0.000001;
    fLumiScale = 1.;

    fThresholdSystPruning_Normalisation = -1;
    fThresholdSystPruning_Shape = -1;
    fThresholdSystLarge = - 1;

    fNtuplePaths.clear();
    fMCweight = "1";
    fSelection = "1";
    fNtupleName = "";

    fHistoPaths.clear();
    fHistoName = "";

    fFitResults = 0;

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

    fRankingMaxNP = 10;
    fReduceNPforRanking = 0.;
    fRankingOnly = "all";
    fRankingPlot = "Merge";
    fAtlasLabel = "Internal";

    fStatOnly = false;
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
    fVarNameMinos.clear();
    fVarNameHide.clear();
    fWorkspaceFileName = "";

    //
    // Limit type
    //
    fLimitType = ASYMPTOTIC;
    fLimitIsBlind = false;
    fLimitPOIAsimov = 0;
    fSignalInjection = false;

    fImageFormat = "png";
    TtHFitter::IMAGEFORMAT.clear();
    TtHFitter::IMAGEFORMAT.push_back("png");

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
    fBlindedBins = 0x0;
    
    fRatioYmax = 1.5;
    fRatioYmin = 0.5;
    fRatioYmaxPostFit = 1.5;
    fRatioYminPostFit = 0.5;
    
    fCustomAsimov = "";
    fRandomPOISeed = -1;
    fTableOptions = "STANDALONE";
    
    fGetGoodnessOfFit = false;
    fGetChi2 = 0; // 0: no, 1: stat-only, 2: with syst
}

//__________________________________________________________________________________
//
TtHFit::~TtHFit(){
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
void TtHFit::SetPOI(string name){
    fPOI = name;
}

//__________________________________________________________________________________
//
void TtHFit::SetStatErrorConfig(bool useIt, float thres, string cons){
    fUseStatErr = useIt;
    fStatErrThres = thres;
    fStatErrCons = cons;
}

//__________________________________________________________________________________
//
void TtHFit::SetLumiErr(float err){
    fLumiErr = err;
}

//__________________________________________________________________________________
//
void TtHFit::SetLumi(const float lumi){
    fLumi = lumi;
}

//__________________________________________________________________________________
//
void TtHFit::SetFitType(FitType type){
    fFitType = type;
}

//__________________________________________________________________________________
//
void TtHFit::SetLimitType(LimitType type){
    fLimitType = type;
}

//__________________________________________________________________________________
//
std::string TtHFit::CheckName( const std::string &name ){
    if( isdigit( name.at(0) ) ){
        std::cerr << "\033[1;31m<!> ERROR in browsing name: " << name << ". A number has been detected at the first position of the name." << std::endl;
        std::cerr << "           This can lead to unexpected behaviours in HistFactory. Please change the name. " << std::endl;
        std::cout << "           The code is about to crash. \033[0m" << std::endl;
        abort();
    } else {
        return name;
    }
}

//__________________________________________________________________________________
//
void TtHFit::SetFitRegion(FitRegion region){
    fFitRegion = region;
}

//__________________________________________________________________________________
//
Sample* TtHFit::NewSample(string name,int type){
    fSamples.push_back(new Sample(name,type));
    //
    fNSamples ++;
    return fSamples[fNSamples-1];
}

//__________________________________________________________________________________
//
Systematic* TtHFit::NewSystematic(string name){
    fSystematics.push_back(new Systematic(name));
    fNSyst ++;
    return fSystematics[fNSyst-1];
}

//__________________________________________________________________________________
//
Region* TtHFit::NewRegion(string name){
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
    fRegions[fNRegions]->fKeepPrefitBlindedBins = fKeepPrefitBlindedBins;
    fRegions[fNRegions]->fRatioYmax = fRatioYmax;
    fRegions[fNRegions]->fRatioYmin = fRatioYmin;
    fRegions[fNRegions]->fRatioYmaxPostFit = fRatioYmaxPostFit;
    fRegions[fNRegions]->fRatioYminPostFit = fRatioYminPostFit;
    //
    fNRegions ++;
    return fRegions[fNRegions-1];
}

//__________________________________________________________________________________
//
void TtHFit::AddNtuplePath(string path){
    fNtuplePaths.push_back(path);
}

//__________________________________________________________________________________
//
void TtHFit::SetMCweight(string weight){
    fMCweight = weight;
}

//__________________________________________________________________________________
//
void TtHFit::SetSelection(string selection){
    fSelection = selection;
}

//__________________________________________________________________________________
//
void TtHFit::SetNtupleName(string name){
    fNtupleName = name;
}

//__________________________________________________________________________________
//
void TtHFit::SetNtupleFile(string name){
    fNtupleFile = name;
}

//__________________________________________________________________________________
//
void TtHFit::AddHistoPath(string path){
    fHistoPaths.push_back(path);
}

// ...

//__________________________________________________________________________________
// apply smoothing to systematics
void TtHFit::SmoothSystematics(string syst){
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Smoothing and/or Symmetrising Systematic Variations ..." << std::endl;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
//             if(fRegions[i_ch]->fSampleHists[i_smp]==0x0) continue;
//             if(fRegions[i_ch]->fSampleHists[i_smp]->fSample!=0x0){
//                 if(fRegions[i_ch]->fSampleHists[i_smp]->fSample->fType==Sample::DATA) continue;
//             }
            fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(syst);
        }
    }
}

//
// Try to split root file creation and histogram wiriting
//__________________________________________________________________________________
// create new root file(s)
void TtHFit::CreateRootFiles(){
    bool recreate = !fUpdate;
    gSystem->mkdir( fName.c_str());
    gSystem->mkdir( (fName + "/Histograms/").c_str() );
    string fileName = "";
    bool singleOutputFile = !TtHFitter::SPLITHISTOFILES;
    //
    if(singleOutputFile){
        if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos" + fSaveSuffix + ".root";
        else                 fileName = fName + "/Histograms/" + fInputName + "_histos" + fSaveSuffix + ".root";
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "Creating/updating file " << fileName << " ..." << std::endl;
        if(recreate) fFiles.push_back(new TFile(fileName.c_str(),"RECREATE"));
        else         fFiles.push_back(new TFile(fileName.c_str(),"UPDATE"));
    }
    else{
        for(int i_ch=0;i_ch<fNRegions;i_ch++){
            if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
            else                 fileName = fName + "/Histograms/" + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
            std::cout << "-------------------------------------------" << std::endl;
            std::cout << "Creating/updating file " << fileName << " ..." << std::endl;
            if(recreate) fFiles.push_back(new TFile(fileName.c_str(),"RECREATE"));
            else         fFiles.push_back(new TFile(fileName.c_str(),"UPDATE"));
        }
    }
}

//__________________________________________________________________________________
// fill files with all the histograms
void TtHFit::WriteHistos(){
    bool singleOutputFile = !TtHFitter::SPLITHISTOFILES;
    SampleHist* sh;
    string fileName = "";
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        //
        if(singleOutputFile) fileName = fFiles[0]   ->GetName();
        else                 fileName = fFiles[i_ch]->GetName();
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "Writing histograms to file " << fileName << " ..." << std::endl;
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if(sh == 0x0){
              if(TtHFitter::DEBUGLEVEL>0) std::cout << "SampleHist[" << i_smp << "] for sample " << fSamples[i_smp]->fName << " not there." << std::endl;
                continue;
            }
            // set file and histo names for nominal
            sh->fHistoName = sh->fHist->GetName();
            sh->fFileName = fileName;
            // set file and histo names for systematics
            for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
                if(sh->fSyst[i_syst]->fHistUp  ==0x0) continue;
                if(sh->fSyst[i_syst]->fHistDown==0x0) continue;
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
            if(singleOutputFile) sh->WriteToFile(fFiles[0]);
            else                 sh->WriteToFile(fFiles[i_ch]);
        }
    }
    std::cout << "-------------------------------------------" << endl;
}

// //__________________________________________________________________________________
// // create new root file with all the histograms
// void TtHFit::WriteHistos(/*string fileName*/){
//     bool recreate = !fUpdate;
//     gSystem->mkdir( fName.c_str());
//     gSystem->mkdir( (fName + "/Histograms/").c_str() );
//     string fileName = "";
//     TDirectory *dir = gDirectory;
//     TFile *f;
//     bool singleOutputFile = !TtHFitter::SPLITHISTOFILES;
//     //
//     if(singleOutputFile){
//         if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos" + fSaveSuffix + ".root";
//         else                 fileName = fName + "/Histograms/" + fInputName + "_histos" + fSaveSuffix + ".root";
//         std::cout << "-------------------------------------------" << std::endl;
//         std::cout << "Writing histograms to file " << fileName << " ..." << std::endl;
//         if(recreate){
//             f = new TFile(fileName.c_str(),"RECREATE");
//             f->~TFile();
//             dir->cd();
//         }
//     }
//     //
//     SampleHist* sh;
//     for(int i_ch=0;i_ch<fNRegions;i_ch++){
//         //
//         if(!singleOutputFile){
//             if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
//             else                 fileName = fName + "/Histograms/" + fInputName + "_" + fRegions[i_ch]->fName + "_histos" + fSaveSuffix + ".root";
//             std::cout << "-------------------------------------------" << std::endl;
//             std::cout << "Writing histograms to file " << fileName << " ..." << std::endl;
//             if(recreate){
//                 f = new TFile(fileName.c_str(),"RECREATE");
//                 f->~TFile();
//                 dir->cd();
//             }
//         }
//         //
//         for(int i_smp=0;i_smp<fNSamples;i_smp++){
//             sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
//             if(sh == 0x0){
//               if(TtHFitter::DEBUGLEVEL>0) std::cout << "SampleHist[" << i_smp << "] for sample " << fSamples[i_smp]->fName << " not there." << std::endl;
//                 continue;
//             }
//             // set file and histo names for nominal
//             sh->fHistoName = sh->fHist->GetName();
//             sh->fFileName = fileName;
//             // set file and histo names for systematics
//             for(int i_syst=0;i_syst<sh->fNSyst;i_syst++){
//                 sh->fSyst[i_syst]->fFileNameUp = fileName;
//                 sh->fSyst[i_syst]->fHistoNameUp = sh->fSyst[i_syst]->fHistUp->GetName();
//                 sh->fSyst[i_syst]->fFileNameDown = fileName;
//                 sh->fSyst[i_syst]->fHistoNameDown = sh->fSyst[i_syst]->fHistDown->GetName();
//                 if(sh->fSyst[i_syst]->fIsShape){
//                     sh->fSyst[i_syst]->fFileNameShapeUp = fileName;
//                     sh->fSyst[i_syst]->fHistoNameShapeUp = sh->fSyst[i_syst]->fHistShapeUp->GetName();
//                     sh->fSyst[i_syst]->fFileNameShapeDown = fileName;
//                     sh->fSyst[i_syst]->fHistoNameShapeDown = sh->fSyst[i_syst]->fHistShapeDown->GetName();
//                 }
//             }
//             sh->WriteToFile();
//         }
//     }
//     std::cout << "-------------------------------------------" << endl;
// }

//__________________________________________________________________________________
// Draw syst plots
void TtHFit::DrawSystPlots(){
    SampleHist* sh;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
//             sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            sh = fRegions[i_ch]->fSampleHists[i_smp];
            sh->DrawSystPlot("all");
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots
void TtHFit::DrawSystPlotsSumSamples(){
  TH1* h_dataCopy;
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    SampleHist* hist = new SampleHist();
    bool empty=true;
    std::set<std::string> systNames;
    for(int i_regSmp=0; i_regSmp<fRegions[i_ch]->fNSamples; i_regSmp++){
      for(int i_smSyst=0; i_smSyst<fRegions[i_ch]->fSampleHists[i_regSmp]->fNSyst; i_smSyst++){
        systNames.insert(fRegions[i_ch]->fSampleHists[i_regSmp]->fSyst[i_smSyst]->fName);
      }
    }
    for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
      if(fSamples[i_smp]->fType==Sample::DATA) h_dataCopy=(TH1*)fRegions[i_ch]->fSampleHists[i_smp]->fHist->Clone();
      else if(fSamples[i_smp]->fType==Sample::GHOST) continue;
      else if(empty){
        hist->CloneSampleHist(fRegions[i_ch]->fSampleHists[i_smp],systNames);
        hist->fName = fRegions[i_ch]->fName + "_Combined";
        empty=false;
      }
      else hist->SampleHistAdd(fRegions[i_ch]->fSampleHists[i_smp]);
    }
    hist->DrawSystPlot("all", h_dataCopy, true, fSystDataPlot_upFrame);
    delete hist;
  }
}

//__________________________________________________________________________________
// Build fit from config file
void TtHFit::ReadConfigFile(string fileName,string options){
    fConfig->ReadFile(fileName);
    ConfigSet *cs; // to store stuff later
    string param;
    std::vector< string > vec;
    int type;
    //
    // Read options (to skip stuff, or include only some regions, samples, systs...)
    // Syntax: .. .. Regions=ge4jge2b:Exclude=singleTop,wjets
    std::map< string,string > optMap; optMap.clear();
    std::vector< string > optVec;
    std::vector< string > onlyRegions; onlyRegions.clear();
    std::vector< string > onlySamples; onlySamples.clear();
    std::vector< string > onlySystematics; onlySystematics.clear();
    std::vector< string > toExclude; toExclude.clear();
    string onlySignal; onlySignal = "";


    //##########################################################
    //
    // COMMAND LINE options
    //
    //##########################################################
    if(options!=""){
        optVec = Vectorize(options,':');
        for(unsigned int i_opt=0;i_opt<optVec.size();i_opt++){
            std::vector< string > optPair;
            optPair = Vectorize(optVec[i_opt],'=');
            optMap[optPair[0]] = optPair[1];
        }
        //
        if(optMap["Regions"]!="")
            onlyRegions = Vectorize(optMap["Regions"],',');
        if(optMap["Samples"]!="")
            onlySamples = Vectorize(optMap["Samples"],',');
        if(optMap["Systematics"]!="")
            onlySystematics = Vectorize(optMap["Systematics"],',');
        if(optMap["Exclude"]!="")
            toExclude = Vectorize(optMap["Exclude"],',');
        if(optMap["Suffix"]!="")
            fSuffix = optMap["Suffix"]; // used for input & output  plots, txt files & workspaces - NOT for histograms file
        if(optMap["SaveSuffix"]!="")
            fSaveSuffix = optMap["SaveSuffix"]; // ... and this one for histograms file
        if(optMap["Update"]!="" && optMap["Update"]!="FALSE")
            fUpdate = true;
        if(optMap["StatOnly"]!="" && optMap["StatOnly"]!="FALSE")
            fStatOnly = true;
        if(optMap["StatOnlyFit"]!="" && optMap["StatOnlyFit"]!="FALSE")
            fStatOnlyFit = true;
        if(optMap["Ranking"]!="")
            fRankingOnly = optMap["Ranking"];
        if(optMap["Signal"]!="")
            onlySignal = optMap["Signal"];
        if(optMap["FitResults"]!="")
            fFitResultsFile = optMap["FitResults"];
        if(optMap["FitType"]!=""){
            if(optMap["FitType"]=="SPLUSB") SetFitType(SPLUSB);
            if(optMap["FitType"]=="BONLY")  SetFitType(BONLY);
        }
        //
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "Running options: " << std::endl;
        if(onlyRegions.size()>0){
            std::cout << "  Only these Regions: " << std::endl;
            for(int i=0;i<onlyRegions.size();i++){
                std::cout << "    " << onlyRegions[i] << std::endl;
            }
        }
        if(onlySamples.size()>0){
            std::cout << "  Only these Samples: " << std::endl;
            for(int i=0;i<onlySamples.size();i++){
                std::cout << "    " << onlySamples[i] << std::endl;
            }
        }
        if(onlySystematics.size()>0){
            std::cout << "  Only these Systematics: " << std::endl;
            for(int i=0;i<onlySystematics.size();i++){
                std::cout << "    " << onlySystematics[i] << std::endl;
            }
        }
        if(toExclude.size()>0){
            std::cout << "  Exclude: " << std::endl;
            for(int i=0;i<toExclude.size();i++){
                std::cout << "    " << toExclude[i] << std::endl;
            }
        }
        if(onlySignal!=""){
            std::cout << "  Only Signal: " << std::endl;
            std::cout << "    " << onlySignal << std::endl;
        }
    }

    //##########################################################
    //
    // JOB options
    //
    //##########################################################
    cs = fConfig->GetConfigSet("Job");

    fName = CheckName(cs->GetValue());
    fInputName = fName;

    param = cs->Get("OutputDir");
    if(param != ""){
      fDir = param;
      if(fDir.back() != '/') fDir += '/';
      fName = fDir + fName;
      gSystem->mkdir(fName.c_str(), true);
    }
    param = cs->Get("Label");  if(param!="") fLabel = param;
                               else          fLabel = fName;
    SetPOI(CheckName(cs->Get("POI")));
    param = cs->Get("ReadFrom");
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(      param=="HIST" || param=="HISTOGRAMS")  fInputType = 0;
        else if( param=="NTUP" || param=="NTUPLES" )    fInputType = 1;
        else{
            std::cerr << "ERROR: Invalid \"ReadFrom\" argument. Options: \"HIST\", \"NTUP\"" << std::endl;
            return;
        }
    if(fInputType==0){
        AddHistoPath( cs->Get("HistoPath") );
    }
    if(fInputType==1){
        SetNtupleFile( cs->Get("NtupleFile") );
        if(cs->Get("NtuplePath")!="") { AddNtuplePath( cs->Get("NtuplePath") ); }
        param = cs->Get("NtuplePaths");
        if( param != "" ){
            std::vector<string> paths = Vectorize( param,',' );
            for(int i=0;i<(int)paths.size();i++){
                AddNtuplePath( paths[i] );
            }
        }
        param = cs->Get("MCweight");  if(param!="") SetMCweight(param);
        param = cs->Get("Selection"); if(param!="") SetSelection(param);
        SetNtupleName( cs->Get("NtupleName") );
    }
    param = cs->Get("Lumi");              if( param != "" ) SetLumi( atof(param.c_str()) );
    param = cs->Get("LumiScale");         if( param != "" ){
        std::cout << "\033[1;33m<!> WARNING: \"LumiScale\" is only done for quick tests since it is inefficient. To normalize all the samples to the luminosity, use \"Lumi\" instead.\033[0m" << std::endl;
        fLumiScale = atof(param.c_str());
    }
    param = cs->Get("SystPruningShape");  if( param != "")  fThresholdSystPruning_Shape         = atof(param.c_str());
    param = cs->Get("SystPruningNorm");   if( param != "")  fThresholdSystPruning_Normalisation = atof(param.c_str());
    param = cs->Get("SystLarge");         if( param != "")  fThresholdSystLarge = atof(param.c_str());
    param = cs->Get("IntCodeOverall");    if( param != "")  fIntCode_overall  = atoi(param.c_str());
    param = cs->Get("IntCodeShape");      if( param != "")  fIntCode_shape    = atoi(param.c_str());
    param = cs->Get("MCstatThreshold");   if( param != ""){
        if(param=="NONE" || param=="None" || param=="none")  SetStatErrorConfig( false, 0. );
        else                                                 SetStatErrorConfig( true,  atof(param.c_str()) );
    }
    else{
        SetStatErrorConfig( true, 0. );
    }
    param = cs->Get("MCstatConstraint");  if( param != "")  fStatErrCons = param;
    param = cs->Get("DebugLevel");        if( param != "")  TtHFitter::SetDebugLevel( atoi(param.c_str()) );
    param = cs->Get("PlotOptions");       if( param != ""){
        vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "YIELDS") !=vec.end() )  TtHFitter::SHOWYIELDS     = true;
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TtHFitter::SHOWSTACKSIG   = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TtHFitter::SHOWNORMSIG    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TtHFitter::SHOWOVERLAYSIG = true;
        if( std::find(vec.begin(), vec.end(), "LEFT")   !=vec.end() )  TtHFitter::LEGENDLEFT     = true;
        if( std::find(vec.begin(), vec.end(), "CHI2")   !=vec.end() )  TtHFitter::SHOWCHI2       = true;
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TtHFitter::PREFITONPOSTFIT= true;
        // ...
    }
    param = cs->Get("PlotOptionsSummary");       if( param != ""){
        vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TtHFitter::SHOWSTACKSIG_SUMMARY   = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TtHFitter::SHOWNORMSIG_SUMMARY    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TtHFitter::SHOWOVERLAYSIG_SUMMARY = true;
    }
    else{
        TtHFitter::SHOWSTACKSIG_SUMMARY   = TtHFitter::SHOWSTACKSIG    ;
        TtHFitter::SHOWNORMSIG_SUMMARY    = TtHFitter::SHOWNORMSIG     ;
        TtHFitter::SHOWOVERLAYSIG_SUMMARY = TtHFitter::SHOWOVERLAYSIG  ;
    }
    param = cs->Get("TableOptions");       if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fTableOptions = param;
    }
    //
    param = cs->Get("SystControlPlots");  if( param != ""){
        if( param == "true" || param == "True" ||  param == "TRUE" ){
            TtHFitter::SYSTCONTROLPLOTS = true;
        } else {
            TtHFitter::SYSTCONTROLPLOTS = false;
        }
    }
    param = cs->Get("SystDataPlots");  if( param != "" ){
        if( param == "true" || param == "True" ||  param == "TRUE" ){
            TtHFitter::SYSTDATAPLOT = true;
            fSystDataPlot_upFrame=false;
        } else if( param == "fillUpFrame" ){
            TtHFitter::SYSTDATAPLOT = true;
            fSystDataPlot_upFrame=true;
        } else {
            TtHFitter::SYSTDATAPLOT = false;
            fSystDataPlot_upFrame=false;
        }
    }
    param = cs->Get("SystErrorBars");  if( param != ""){
        if( param == "true" || param == "True" ||  param == "TRUE" ){
            TtHFitter::SYSTERRORBARS = true;
        }
    }
    param = cs->Get("CorrelationThreshold"); if( param != ""){
        TtHFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }
    param = cs->Get("SignalRegionsPlot");  if(param != ""){
        fRegionsToPlot = Vectorize(param,',');
    }
    param = cs->Get("SummaryPlotRegions");  if(param != ""){
        fSummaryPlotRegions = Vectorize(param,',');
    }
    param = cs->Get("SummaryPlotLabels");  if(param != ""){
        fSummaryPlotLabels = Vectorize(param,',');
    }
    param = cs->Get("SummaryPlotValidationRegions");  if(param != ""){
        fSummaryPlotValidationRegions = Vectorize(param,',');
    }
    param = cs->Get("SummaryPlotValidationLabels");  if(param != ""){
        fSummaryPlotValidationLabels = Vectorize(param,',');
    }
    param = cs->Get("SummaryPlotYmin");  if(param != "") fYmin = atof(param.c_str());
    param = cs->Get("SummaryPlotYmax");  if(param != "") fYmax = atof(param.c_str());
    param = cs->Get("RatioYmin");  if(param != "") { fRatioYmin = atof(param.c_str()); fRatioYminPostFit = fRatioYmin; }
    param = cs->Get("RatioYmax");  if(param != "") { fRatioYmax = atof(param.c_str()); fRatioYmaxPostFit = fRatioYmax; }  
    param = cs->Get("RatioYminPostFit");  if(param != "") fRatioYminPostFit = atof(param.c_str());
    param = cs->Get("RatioYmaxPostFit");  if(param != "") fRatioYmaxPostFit = atof(param.c_str());    
    param = cs->Get("HistoChecks");  if(param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "NOCRASH" ){
            TtHFitter::HISTOCHECKCRASH = false;
        }
    }
    param = cs->Get("LumiLabel"); if( param != "") fLumiLabel = param;
    param = cs->Get("CmeLabel"); if( param != "") fCmeLabel = param;
    param = cs->Get("SplitHistoFiles");  if( param != ""){
        if( param == "true" || param == "True" ||  param == "TRUE" ){
            TtHFitter::SPLITHISTOFILES = true;
        } else {
            TtHFitter::SPLITHISTOFILES = false;
        }
    }
    param = cs->Get("BlindingThreshold");  if( param != ""){
        fBlindingThreshold = atof(param.c_str());
    }
    param = cs->Get("RankingMaxNP");  if( param != ""){
        fRankingMaxNP = atoi(param.c_str());
    }
    param = cs->Get("RankingPlot");  if( param != ""){
        fRankingPlot = param;
    }
    param = cs->Get("ReduceNPforRanking");  if( param != ""){
        fReduceNPforRanking = atof(param.c_str());
    }
    param = cs->Get("ImageFormat");  if( param != ""){
        fImageFormat = Vectorize(param,',')[0];
        TtHFitter::IMAGEFORMAT = Vectorize(param,',');
    }
    param = cs->Get("StatOnly");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fStatOnly = true;
        }
    }
    param = cs->Get("FixNPforStatOnly");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFixNPforStatOnlyFit = true;
        }
    }
    param = cs->Get("InputFolder");    if( param != "" ){
        fInputFolder = param;
    }
    param = cs->Get("InputName");    if( param != "" ){
        fInputName = param;
    }
    param = cs->Get("WorkspaceFileName");    if( param != "" ){
        fWorkspaceFileName = param;
    }
    param = cs->Get("KeepPruning");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fKeepPruning = true;
    }
    param = cs->Get("AtlasLabel"); if( param != "" ){
        fAtlasLabel = param;
    }
    param = cs->Get("CleanTables");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fCleanTables = true;
    }
    param = cs->Get("SystCategoryTables");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fSystCategoryTables = true;
    }
    param = cs->Get("Suffix"); if( param != "" ){
        fSuffix = param;
    }
    param = cs->Get("SaveSuffix"); if( param != "" ){
        fSaveSuffix = param;
    }
    param = cs->Get("HideNP"); if( param != "" ){
        fVarNameHide = Vectorize(param,',');
    }
    param = cs->Get("RegionGroups"); if( param != "" ) {
        vector<string> groups =   Vectorize(param,',');
        for(unsigned int i_gr=0;i_gr<groups.size();i_gr++)
        fRegionGroups.push_back(groups[i_gr]);
    }
    param = cs->Get("KeepPrefitBlindedBins");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ) fKeepPrefitBlindedBins = true;
    }
    param = cs->Get("CustomAsimov");    if( param != "" ){
//         std::transform(param.begin(), param.end(), param.begin(), ::toupper);
//         if( param == "TRUE" ) fCustomAsimov = true;
        fCustomAsimov = param;
    }
    param = cs->Get("RandomPOISeed");   if( param != "" ){
        int seed = atoi(param.c_str());
        if(seed>=0) fRandomPOISeed = seed;
    }
    param = cs->Get("GetChi2");   if( param != "" ){ // can be TRUE, SYST+STAT, STAT-ONLY... (if it contains STAT and no SYST => stat-only, ptherwise stat+syst)
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fGetChi2 = 2;
        }
        if( param.find("SYST")!=std::string::npos ){
            fGetChi2 = 2;
        }
        else if( param.find("STAT")!=std::string::npos ){
            fGetChi2 = 1;
        }
    }
    
    //
    // General options
    //
    cs = fConfig->GetConfigSet("Options");
    if(cs!=0x0){
        for(int i=0;i<cs->GetN();i++){
            if(cs->GetConfigValue(i) != ""){
                TtHFitter::OPTION[cs->GetConfigName(i)] = atof(cs->GetConfigValue(i).c_str());
            }
        }
    }

    //##########################################################
    //
    // FIT options
    //
    //##########################################################
    cs = fConfig->GetConfigSet("Fit");
    param = cs->Get("FitType");    if( param != "" && fFitType == UNDEFINED ){
        if( param == "SPLUSB" )
            SetFitType(TtHFit::SPLUSB);
        else if( param == "BONLY" )
            SetFitType(TtHFit::BONLY);
        else{
            std::cerr << "Unknown FitType argument : " << cs->Get("FitType") << std::endl;
            return;
        }
    }
    else if( fFitType == UNDEFINED ){
        std::cout << "TtHFit::INFO : Setting default fit Type SPLUSB" << std::endl;
        SetFitType(TtHFit::SPLUSB);
    }
    param = cs->Get("FitRegion");    if( param != "" ){
        if( param == "CRONLY" )
            SetFitRegion(TtHFit::CRONLY);
        else if( param == "CRSR" )
            SetFitRegion(TtHFit::CRSR);
        else{
            SetFitRegion(TtHFit::USERSPECIFIC);
            fFitRegionsToFit = Vectorize(param,',');
            if(fFitRegionsToFit.size()==0){
                std::cerr << "Unknown FitRegion argument : " << cs->Get("FitRegion") << std::endl;
                return;
            }
        }
    }
    param = cs->Get("FitBlind");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitIsBlind = true;
        } else if ( param == "FALSE" ){
            fFitIsBlind = false;
        }
    }
    param = cs->Get("POIAsimov");   if( param != "" ){ fFitPOIAsimov = atof(param.c_str()); };
    param = cs->Get("NPValues");    if( param != "" ){
        std::vector < std::string > temp_vec = Vectorize(param,',');
        for(unsigned int iNP = 0; iNP < temp_vec.size(); ++iNP){
            std::vector < std::string > np_value = Vectorize(temp_vec[iNP],':');
            if(np_value.size()==2){
                fFitNPValues.insert( std::pair < std::string, double >( np_value[0], atof(np_value[1].c_str()) ) );
            }
        }
    }
    param = cs->Get("FixNPs");    if( param != "" ){
        std::vector < std::string > temp_fixedNPs = Vectorize(param,',');
        for(unsigned int iNP = 0; iNP < temp_fixedNPs.size(); ++iNP){
            std::vector < std::string > fixed_nps = Vectorize(temp_fixedNPs[iNP],':');
            if(fixed_nps.size()==2){
                fFitFixedNPs.insert( std::pair < std::string, double >( fixed_nps[0], atof(fixed_nps[1].c_str()) ) );
            }
        }
    }
    param = cs->Get("doLHscan"); if( param != "" ){ fVarNameLH = Vectorize(param,','); }
    param = cs->Get("UseMinos"); if( param != "" ){ fVarNameMinos = Vectorize(param,','); }
    param = cs->Get("SetRandomInitialNPval");  if( param != ""){
        fUseRnd = true;
        fRndRange = atof(param.c_str());
    }
    param = cs->Get("SetRandomInitialNPvalSeed"); if( param != ""){
      fRndSeed = atol(param.c_str());
    }
    param = cs->Get("NumCPU"); if( param != "" ){ TtHFitter::NCPU = atoi( param.c_str()); }
    param = cs->Get("StatOnlyFit");    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fStatOnlyFit = true;
        }
    }
    param = cs->Get("GetGoodnessOfFit");   if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fGetGoodnessOfFit = true;
        }
    }

    //##########################################################
    //
    // Reads LIMIT parameters
    //
    //##########################################################
    cs = fConfig->GetConfigSet("Limit");
    if (cs) {
        param = cs->Get("LimitType");    if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param == "ASYMPTOTIC" )
                SetLimitType(TtHFit::ASYMPTOTIC);
            else if( param == "TOYS" )
                SetLimitType(TtHFit::TOYS);
            else{
                std::cerr << "Unknown LimitType argument : " << cs->Get("LimitType") << std::endl;
                return;
            }
        }
        param = cs->Get("LimitBlind");    if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param == "TRUE" ){
                fLimitIsBlind = true;
            } else if ( param == "FALSE" ){
                fLimitIsBlind = false;
            }
        }
        param = cs->Get("POIAsimov");  if( param != "" ){ fLimitPOIAsimov = atof(param.c_str()); }
        param = cs->Get("SignalInjection");  if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param == "TRUE" ){
                fSignalInjection = true;
            } else if ( param == "FALSE" ){
                fSignalInjection = false;
            }
        }
    }

    //##########################################################
    //
    // REGIONS options
    //
    //##########################################################
    int nReg = 0;
    Region *reg;
    std::vector<std::string> regNames;
    while(true){
        cs = fConfig->GetConfigSet("Region",nReg);
        if(cs==0x0) break;
        nReg++;
        if(onlyRegions.size()>0 && FindInStringVector(onlyRegions,cs->GetValue())<0) continue;
        if(toExclude.size()>0 && FindInStringVector(toExclude,cs->GetValue())>=0) continue;
        regNames.push_back( CheckName(cs->GetValue()) ); //why the CheckName is needed?? A: cs->GetValue() might have leading/trailing spaces...
        reg = NewRegion(CheckName(cs->GetValue()));
        reg->fGetChi2 = fGetChi2;
        reg->SetVariableTitle(cs->Get("VariableTitle"));
        reg->SetLabel(cs->Get("Label"),cs->Get("ShortLabel"));
        param = cs->Get("YaxisTitle"); if( param != "") reg->fYTitle = param;
        param = cs->Get("YmaxScale");  if(param!="") reg->fYmaxScale = atof(param.c_str());
        param = cs->Get("Ymin");       if(param!="") reg->fYmin = atof(param.c_str());
        param = cs->Get("Ymax");       if(param!="") reg->fYmax = atof(param.c_str());
        param = cs->Get("RatioYmin");  if(param!="") { reg->fRatioYmin = atof(param.c_str()); reg->fRatioYminPostFit = reg->fRatioYmin; }
        param = cs->Get("RatioYmax");  if(param!="") { reg->fRatioYmax = atof(param.c_str()); reg->fRatioYmaxPostFit = reg->fRatioYmax; }
        param = cs->Get("RatioYminPostFit");  if(param!="") reg->fRatioYminPostFit = atof(param.c_str());
        param = cs->Get("RatioYmaxPostFit");  if(param!="") reg->fRatioYmaxPostFit = atof(param.c_str());
        param = cs->Get("TexLabel");   if( param != "") reg->fTexLabel = param;
        param = cs->Get("LumiLabel");  if( param != "") reg->fLumiLabel = param;
        else reg->fLumiLabel = fLumiLabel;
        param = cs->Get("CmeLabel"); if( param != "") reg->fCmeLabel = param;
        else reg->fCmeLabel = fCmeLabel;
        param = cs->Get("LogScale"); if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") reg->fLogScale = true;
            if(param=="FALSE") reg->fLogScale = false;
        }
        param = cs->Get("Group"); if( param != "") reg->fGroup = param;
        if(fInputType==0){
            param = cs->Get("HistoFile"); if(param!="") reg->fHistoFiles.push_back( param );
            param = cs->Get("HistoName"); if(param!="") reg->SetHistoName( param );
            if(cs->Get("HistoPathSuff")!="") { reg->fHistoPathSuffs.clear(); reg->fHistoPathSuffs.push_back( Fix(cs->Get("HistoPathSuff")) ); }
            param = cs->Get("HistoPathSuffs");
            if( param != "" ){
                reg->fHistoPathSuffs.clear();
                std::vector<string> paths = Vectorize( param,',' );
                for(int i=0;i<(int)paths.size();i++){
                    reg->fHistoPathSuffs.push_back( Fix(paths[i]) );
                }
            }
        }
        else if(fInputType==1){
            vector<string> variable = Vectorize(cs->Get("Variable"),',');
            // fix variable vector if special functions are used
            if(variable[0].find("Alt$")!=string::npos || variable[0].find("MaxIf$")!=string::npos ||variable[0].find("MinIf$")!=string::npos ){
                variable[0]+=","+variable[1];
                variable.erase(variable.begin()+1);
            }
            vector<string> corrVar  = Vectorize(variable[0],'|');
            if(corrVar.size()==2){
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Have a correlation variable in reg " << regNames.back() << " : "
                                                      << corrVar[0] << " and " << corrVar[1] << std::endl;
                reg->SetVariable(  "corr_"+corrVar[0]+"_"+corrVar[1], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()), corrVar[0].c_str(), corrVar[1].c_str() );
            }
            else{
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Have a usual variable in reg " << regNames.back() << " : "
                                                      << variable[0] << " and size of corrVar=" << corrVar.size() << std::endl;
                reg->SetVariable(  variable[0], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()) );
            }
            //
            param = cs->Get("VariableForSample");    if( param != "" ){
                std::vector < std::string > temp_samplesAndVars = Vectorize(param,',');
                for(unsigned int ivar = 0; ivar < temp_samplesAndVars.size(); ++ivar){
                  std::vector < std::string > vars = Vectorize(temp_samplesAndVars[ivar],':');
                    if(vars.size()==2){
                        reg->SetAlternativeVariable(vars[1], vars[0]);
                    }
                }
            }
            //
            if(cs->Get("Selection")!="") reg->AddSelection( cs->Get("Selection") );
            param = cs->Get("NtupleName"); if(param!="") { reg->fNtupleNames.clear(); reg->fNtupleNames.push_back(param); }
            if(cs->Get("NtupleNameSuff")!="") { reg->fNtupleNameSuffs.clear(); reg->fNtupleNameSuffs.push_back( cs->Get("NtupleNameSuff") ); }
            param = cs->Get("NtupleNameSuffs");
            if( param != "" ){
                reg->fNtupleNameSuffs.clear();
                std::vector<string> paths = Vectorize( param,',' );
                for(int i=0;i<(int)paths.size();i++){
                    reg->fNtupleNameSuffs.push_back( paths[i] );
                }
            }
            reg->fMCweight = cs->Get("MCweight"); // this will override the global MCweight, if any
            if(cs->Get("NtuplePathSuff")!="") { reg->fNtuplePathSuffs.clear(); reg->fNtuplePathSuffs.push_back( cs->Get("NtuplePathSuff") ); }
            param = cs->Get("NtuplePathSuffs");
            if( param != "" ){
                reg->fNtuplePathSuffs.clear();
                std::vector<string> paths = Vectorize( param,',' );
                for(int i=0;i<(int)paths.size();i++){
                    reg->fNtuplePathSuffs.push_back( Fix(paths[i]) );
                }
            }
        }
        //Potential rebinning
        if(cs->Get("Rebin")!="") reg -> Rebin(atoi(cs->Get("Rebin").c_str()));
        if(cs->Get("Binning")!="" && cs->Get("Binning")!="-"){
            std::vector < string > vec_bins = Vectorize(cs->Get("Binning"), ',');
            if(vec_bins[0]=="AutoBin"){
                reg -> fBinTransfo = vec_bins[1];
                if(vec_bins[1]=="TransfoD"){
                    reg -> fTransfoDzSig=convertStoD(vec_bins[2]);
                    reg -> fTransfoDzBkg=convertStoD(vec_bins[3]);
		    if(vec_bins.size()>4){
		      for(unsigned int i_Bkgs=3; i_Bkgs<vec_bins.size(); ++i_Bkgs){
			reg -> fAutoBinBkgsInSig.push_back(vec_bins[i_Bkgs]);
		      }
		    }
                }
                else if(vec_bins[1]=="TransfoF"){
                    reg -> fTransfoFzSig=convertStoD(vec_bins[2]);
                    reg -> fTransfoFzBkg=convertStoD(vec_bins[3]);
		    if(vec_bins.size()>4){
		      for(unsigned int i_Bkgs=3; i_Bkgs<vec_bins.size(); ++i_Bkgs){
			reg -> fAutoBinBkgsInSig.push_back(vec_bins[i_Bkgs]);
		      }
		    }
                }
                else if(vec_bins[1]=="TransfoJ"){
                    if(vec_bins.size() > 2) reg -> fTransfoJpar1=convertStoD(vec_bins[2]);
                    else reg -> fTransfoJpar1 = 5.;
                    if(vec_bins.size() > 3) reg -> fTransfoJpar2=convertStoD(vec_bins[3]);
                    else reg -> fTransfoJpar2 = 1.;
                    if(vec_bins.size() > 4) reg -> fTransfoJpar3=convertStoD(vec_bins[4]);
                    else reg -> fTransfoJpar3 = 5.;
		    if(vec_bins.size()>5){
		      for(unsigned int i_Bkgs=4; i_Bkgs<vec_bins.size(); ++i_Bkgs){
			reg -> fAutoBinBkgsInSig.push_back(vec_bins[i_Bkgs]);
		      }
		    }
                }
                else{
                    std::cout<<" ERROR: Unknown transformation: "<<vec_bins[1]<<", try again" << std::endl;
                    exit(1);
                }
            }
            else{
                const int nBounds = vec_bins.size();
                double bins[nBounds];
                for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
                    bins[iBound] = atof(vec_bins[iBound].c_str());
                }
                reg -> SetBinning(nBounds-1,bins);
            }
        }
        if(cs->Get("BinWidth")!="") reg->fBinWidth = atof(cs->Get("BinWidth").c_str());
        if(cs->Get("Type")!=""){
            param = cs->Get("Type");
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="CONTROL" )     reg -> SetRegionType(Region::CONTROL);
            if( param=="VALIDATION" )  reg -> SetRegionType(Region::VALIDATION);
            if( param=="SIGNAL" )      reg -> SetRegionType(Region::SIGNAL);
        }
        if(cs->Get("DataType")!=""){
            param = cs->Get("DataType");
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="DATA" )     reg -> SetRegionDataType(Region::REALDATA);
            else if( param=="ASIMOV" )  reg -> SetRegionDataType(Region::ASIMOVDATA);
            else{
                std::cout << "<!> DataType is not recognised: " << param << std::endl;
            }
        }
        param = cs->Get("SkipSmoothing"); if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE")  reg->fSkipSmoothing = true;
        }
        param = cs->Get("DropBins");
        if( param != "" ){
            reg->fDropBins.clear();
            std::vector<std::string> s = Vectorize( param,',' );
            for(int i=0;i<(int)s.size();i++){
                reg->fDropBins.push_back(atoi(s[i].c_str()));
            }
        }
        param = cs->Get("BinLabels");
        if( param != "" ){
            reg->fBinLabels.clear();
            std::vector<std::string> vec_string = Vectorize( param,',' );
            for(std::string s: vec_string) {
              reg->fBinLabels.push_back(s);
            }
        }
    }

    //##########################################################
    //
    // SAMPLES options
    //
    //##########################################################
    int nSmp = 0;
    Sample *smp;
    NormFactor *nf = 0x0;
    ShapeFactor *sf = 0x0;
    while(true){
        cs = fConfig->GetConfigSet("Sample",nSmp);
        if(cs==0x0) break;
        nSmp++;
        if(onlySamples.size()>0 && FindInStringVector(onlySamples,cs->GetValue())<0) continue;
        if(toExclude.size()>0 && FindInStringVector(toExclude,cs->GetValue())>=0) continue;
        type = Sample::BACKGROUND;
        if(cs->Get("Type")=="signal" || cs->Get("Type")=="SIGNAL") type = Sample::SIGNAL;
        if(cs->Get("Type")=="data"   || cs->Get("Type")=="DATA")   type = Sample::DATA;
        if(cs->Get("Type")=="ghost"  || cs->Get("Type")=="GHOST")  type = Sample::GHOST;
        if(onlySignal!="" && type==Sample::SIGNAL && cs->GetValue()!=onlySignal) continue;
        smp = NewSample(CheckName(cs->GetValue()),type);
        smp->SetTitle(cs->Get("Title"));
        param = cs->Get("TexTitle"); if(param!="") smp->fTexTitle = param;
        param = cs->Get("Group"); if(param!="") smp->fGroup = param;
        if(fInputType==0){
            param = cs->Get("HistoFile"); if(param!="") smp->AddHistoFile( param );
            param = cs->Get("HistoName"); if(param!="") smp->fHistoNames.push_back( param );
            param = cs->Get("HistoPath"); if(param!="") smp->AddHistoPath( param );
        }
        if(fInputType==1){
            // ntuple files
            param = cs->Get("NtupleFile");
            if(param!="") smp->AddNtupleFile( param );
            param = cs->Get("NtupleFiles");
            if(param!=""){
                smp->fNtupleFiles = Vectorize( param ,',' );
            }
            param = cs->Get("NtupleName");
            if(param!="") smp->AddNtupleName( param );
            param = cs->Get("NtupleNames");
            if(param!=""){
                smp->fNtupleNames = Vectorize( param ,',' );
            }
            // ntuple paths
            param = cs->Get("NtuplePath");
            if(param!="") smp->AddNtuplePath( param );
            param = cs->Get("NtuplePaths");
            if(param!=""){
                smp->fNtuplePaths = Vectorize( param ,',' );
            }
            if(cs->Get("NtupleNameSuff")!="") { smp->fNtupleNameSuffs.clear(); smp->fNtupleNameSuffs.push_back( cs->Get("NtupleNameSuff") ); }
            param = cs->Get("NtupleNameSuffs");
            if( param != "" ){
                smp->fNtupleNameSuffs.clear();
                std::vector<string> paths = Vectorize( param,',' );
                for(int i=0;i<(int)paths.size();i++){
                    smp->fNtupleNameSuffs.push_back( paths[i] );
                }
            }
        }
        if(cs->Get("FillColor")!="")
            smp->SetFillColor(atoi(cs->Get("FillColor").c_str()));
        if(cs->Get("LineColor")!="")
            smp->SetLineColor(atoi(cs->Get("LineColor").c_str()));
        if(cs->Get("NormFactor")!=""){
            // check if the normfactor is called just with the name or with full definition
            if( Vectorize(cs->Get("NormFactor"),',').size()>1 ){
                bool isConst = false;
                if( Vectorize(cs->Get("NormFactor"),',').size()>4 &&
                  (    Vectorize(cs->Get("NormFactor"),',')[4] == "TRUE"
                    || Vectorize(cs->Get("NormFactor"),',')[4] == "true"
                    || Vectorize(cs->Get("NormFactor"),',')[4] == "True" ) ) isConst = true;
                nf = smp->AddNormFactor(
                    Vectorize(cs->Get("NormFactor"),',')[0],
                    atof(Vectorize(cs->Get("NormFactor"),',')[1].c_str()),
                    atof(Vectorize(cs->Get("NormFactor"),',')[2].c_str()),
                    atof(Vectorize(cs->Get("NormFactor"),',')[3].c_str()),
                    isConst
                );
            }
            else{
                nf = smp->AddNormFactor( Vectorize(cs->Get("NormFactor"),',')[0] );
            }
            if( FindInStringVector(fNormFactorNames,nf->fName)<0 ){
                fNormFactors.push_back( nf );
                fNormFactorNames.push_back( nf->fName );
                fNNorm++;
            }
        }
        if(cs->Get("ShapeFactor")!=""){
            // check if the shapefactor is called just with the name or with full definition
            if( Vectorize(cs->Get("ShapeFactor"),',').size()>1 ){
                bool isConst = false;
                if( Vectorize(cs->Get("ShapeFactor"),',').size()>4 &&
                  (    Vectorize(cs->Get("ShapeFactor"),',')[4] == "TRUE"
                    || Vectorize(cs->Get("ShapeFactor"),',')[4] == "true"
                    || Vectorize(cs->Get("ShapeFactor"),',')[4] == "True" ) ) isConst = true;
                sf = smp->AddShapeFactor(
                    Vectorize(cs->Get("ShapeFactor"),',')[0],
                    atof(Vectorize(cs->Get("ShapeFactor"),',')[1].c_str()),
                    atof(Vectorize(cs->Get("ShapeFactor"),',')[2].c_str()),
                    atof(Vectorize(cs->Get("ShapeFactor"),',')[3].c_str()),
                    isConst
                );
            }
            else{
                sf = smp->AddShapeFactor( Vectorize(cs->Get("ShapeFactor"),',')[0] );
            }
            if( FindInStringVector(fShapeFactorNames,sf->fName)<0 ){
                fShapeFactors.push_back( sf );
                fShapeFactorNames.push_back( sf->fName );
                fNShape++;
            }
        }
        if(cs->Get("NormalizedByTheory")!=""){
            param = cs->Get("NormalizedByTheory");
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") smp->NormalizedByTheory(false);
            else if(param=="TRUE") smp->NormalizedByTheory(true);
            else std::cout << "<!> NormalizedByTheory flag not recognized ... *" << param << "*" << std::endl;
        }
        if(fInputType==1){
            param = cs->Get("MCweight");
            if(param!="")  smp->SetMCweight( param );
            param = cs->Get("Selection");
            if(param!="")  smp->SetSelection( param );
        }
        // to specify only certain regions
        string regions_str = cs->Get("Regions");
        string exclude_str = cs->Get("Exclude");
        vector<string> regions = Vectorize(regions_str,',');
        vector<string> exclude = Vectorize(exclude_str,',');
        smp->fRegions.clear();
        for(int i_reg=0;i_reg<fNRegions;i_reg++){
            string regName = fRegions[i_reg]->fName;
            if( (regions_str=="" || regions_str=="all" || FindInStringVector(regions,regName)>=0)
                && FindInStringVector(exclude,regName)<0 ){
                smp->fRegions.push_back( fRegions[i_reg]->fName );
            }
        }
        // to scale sample by a factor (can be different for each input ntuple / histogram)
        // NOTE: be careful when speficifying more than one file and more than one path at the same time!!
        param = cs->Get("LumiScale");  if(param!="") smp->fLumiScales.push_back( atof(param.c_str()) );
        param = cs->Get("LumiScales"); if(param!=""){
            vector<string> lumiScales_str = Vectorize( param ,',' );
            for(unsigned int i=0;i<lumiScales_str.size();i++)
                smp->fLumiScales.push_back( atof(lumiScales_str[i].c_str()) );
        }
        // to skip global & region selection for this sample
        param = cs->Get("IgnoreSelection");
        if(param!=""){
            std::string param1 = param;
            std::transform(param1.begin(), param1.end(), param1.begin(), ::toupper);
            if(param1 == "TRUE") smp->fIgnoreSelection = "TRUE";
            else if(param1 == "FALSE") smp->fIgnoreSelection = "FALSE";
            else smp->fIgnoreSelection = param;
        }
        // to skip MC stat uncertainty for this sample
        param = cs->Get("UseMCstat");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") smp->fUseMCStat = false;
        }
        // to skip MC systematics for this sample
        param = cs->Get("UseSystematics");
        // set it to false for ghost samples and data and true for other samples
        if(type == Sample::GHOST || type == Sample::DATA) smp->fUseSystematics = false;
        else                                              smp->fUseSystematics = true;
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") smp->fUseSystematics = false;
            if(param == "TRUE" ) smp->fUseSystematics = true;
        }
        // divide and multiply by another sample
        smp->fDivideBy   = cs->Get("DivideBy");
        smp->fMultiplyBy = cs->Get("MultiplyBy");
        // add/subtract
        param = cs->Get("SubtractSample");
        if(param!=""){
            smp->fSubtractSamples.push_back(param);
        }
        param = cs->Get("SubtractSamples");
        if(param!=""){
            smp->fSubtractSamples = Vectorize(param,',');
        }
        param = cs->Get("AddSample");
        if(param!=""){
            smp->fAddSamples.push_back(param);
        }
        param = cs->Get("AddSamples");
        if(param!=""){
            smp->fAddSamples = Vectorize(param,',');
        }
        // allow smoothing of nominal histogram?
        param = cs->Get("Smooth");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE") smp->fSmooth = true;
        }
        // AsimovReplacementFor
        param = cs->Get("AsimovReplacementFor");
        if(param!=""){
            if(Vectorize(param,',').size()>1){
                smp->fAsimovReplacementFor.first  = Vectorize(param,',')[0];
                smp->fAsimovReplacementFor.second = Vectorize(param,',')[1];
            }
        }
        // separate gammas
        param = cs->Get("SeparateGammas");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE"){
                smp->fSeparateGammas = true;
                if(cs->Get("UseMCstat")=="") smp->fUseMCStat = false; // remove the usual gammas for this sample (only if no UseMCstat is specified!!)
            }
        }
        // in the form    CorrelateGammasInRegions: SR1:SR2,CR1:CR2:CR3
        param = cs->Get("CorrelateGammasInRegions");
        if(param!=""){
            std::vector<std::string> sets = Vectorize(param,',');
            for(auto set : sets){
                std::vector<std::string> regions = Vectorize(set,':');
                std::cout << "Correlating gammas for this sample in regions " << set << std::endl;
                smp->fCorrelateGammasInRegions.push_back(regions);
            }
        }
        // ...
    }
    //
    // build new samples if AsimovReplacementFor are specified
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fAsimovReplacementFor.first!=""){
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "Creating sample " << fSamples[i_smp]->fAsimovReplacementFor.first << std::endl;
            Sample *ca = NewSample("customAsimov_"+fSamples[i_smp]->fAsimovReplacementFor.first,Sample::GHOST);
            ca->SetTitle("Pseudo-Data ("+fSamples[i_smp]->fAsimovReplacementFor.first+")");
            ca->fUseSystematics = false;
        }
    }

    //##########################################################
    //
    // NORMFACTOR options
    //
    //##########################################################
    int nNorm = 0;
    NormFactor *norm;
    Sample *sam;

    while(true){
        cs = fConfig->GetConfigSet("NormFactor",nNorm);
        if(cs==0x0) break;
        nNorm++;
        if(toExclude.size()>0 && FindInStringVector(toExclude,cs->GetValue())>=0) continue;
        string samples_str = cs->Get("Samples");
        string regions_str = cs->Get("Regions");
        string exclude_str = cs->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        vector<string> samples = Vectorize(samples_str,',');
        vector<string> regions = Vectorize(regions_str,',');
        vector<string> exclude = Vectorize(exclude_str,',');
        norm = new NormFactor(CheckName(cs->GetValue()));
        TtHFitter::SYSTMAP[norm->fName] = norm->fName;
        if( FindInStringVector(fNormFactorNames,norm->fName)<0 ){
            fNormFactors.push_back( norm );
            fNormFactorNames.push_back( norm->fName );
            fNNorm++;
        }
        else{
            norm = fNormFactors[ FindInStringVector(fNormFactorNames,norm->fName) ];
        }
        if(cs->Get("NuisanceParameter")!=""){
            norm->fNuisanceParameter = cs->Get("NuisanceParameter");
            TtHFitter::NPMAP[norm->fName] = norm->fNuisanceParameter;
        }
        else{
            norm->fNuisanceParameter = norm->fName;
            TtHFitter::NPMAP[norm->fName] = norm->fName;
        }
        param = cs->Get("Constant");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") norm->fConst = true;
        }
        param = cs->Get("Category");
        if(param!=""){
            norm->fCategory = param;
        }
        param = cs->Get("Title"); if(param!=""){
            norm->fTitle = param;
            TtHFitter::SYSTMAP[norm->fName] = norm->fTitle;
        }
        param = cs->Get("TexTitle"); if(param!=""){
            TtHFitter::SYSTTEX[norm->fName] = param;
        }
        param = cs->Get("Min"); if(param!=""){ norm->fMin = atof(param.c_str()); }
        param = cs->Get("Max"); if(param!=""){ norm->fMax = atof(param.c_str()); }
        param = cs->Get("Nominal"); if(param!=""){ norm->fNominal = atof(param.c_str()); }
        //
        // save list of
        if(regions[0]!="all") norm->fRegions = regions;
        if(exclude[0]!="")    norm->fExclude = exclude;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sam = fSamples[i_smp];
            if(sam->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
               && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
            ){
                sam->AddNormFactor(norm);
            }
        }
        // ...
    }

    //##########################################################
    //
    // ShapeFactor options
    //
    //##########################################################
    int nShape = 0;
    ShapeFactor *shape;
    //    Sample *sam;

    while(true){
        cs = fConfig->GetConfigSet("ShapeFactor",nShape);
        if(cs==0x0) break;
        nShape++;
        if(toExclude.size()>0 && FindInStringVector(toExclude,cs->GetValue())>=0) continue;
        string samples_str = cs->Get("Samples");
        string regions_str = cs->Get("Regions");
        string exclude_str = cs->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        vector<string> samples = Vectorize(samples_str,',');
        vector<string> regions = Vectorize(regions_str,',');
        vector<string> exclude = Vectorize(exclude_str,',');
        shape = new ShapeFactor(CheckName(cs->GetValue()));
        if( FindInStringVector(fShapeFactorNames,shape->fName)<0 ){
            fShapeFactors.push_back( shape );
            fShapeFactorNames.push_back( shape->fName );
            fNShape++;
        }
        else{
            shape = fShapeFactors[ FindInStringVector(fShapeFactorNames,shape->fName) ];
        }
        if(cs->Get("NuisanceParameter")!=""){
            shape->fNuisanceParameter = cs->Get("NuisanceParameter");
            TtHFitter::NPMAP[shape->fName] = shape->fNuisanceParameter;
        }
        else{
            shape->fNuisanceParameter = shape->fName;
            TtHFitter::NPMAP[shape->fName] = shape->fName;
        }
        param = cs->Get("Constant");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") shape->fConst = true;
        }
        param = cs->Get("Category");
        if(param!=""){
            shape->fCategory = param;
        }
        param = cs->Get("Title"); if(param!=""){
            shape->fTitle = param;
            TtHFitter::SYSTMAP[shape->fName] = shape->fTitle;
        }
        param = cs->Get("TexTitle"); if(param!=""){
            TtHFitter::SYSTTEX[shape->fName] = param;
        }
        param = cs->Get("Min"); if(param!=""){ shape->fMin = atof(param.c_str()); }
        param = cs->Get("Max"); if(param!=""){ shape->fMax = atof(param.c_str()); }
        param = cs->Get("Nominal"); if(param!=""){ shape->fNominal = atof(param.c_str()); }
        //
        // save list of
        if(regions[0]!="all") shape->fRegions = regions;
        if(exclude[0]!="")    shape->fExclude = exclude;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sam = fSamples[i_smp];
            if(sam->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
               && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
            ){
                sam->AddShapeFactor(shape);
            }
        }
        // ...
    }


    //##########################################################
    //
    // SYSTEMATICS options
    //
    //##########################################################
    int nSys = 0;
    Systematic *sys;
//     Sample *sam;

    //Addition for StatOnly fit: dummy systematic for the significance computation and limit setting
    int typed=0;
    Systematic *sysd;
    if (fStatOnly) {
        typed = Systematic::OVERALL;
        sysd = new Systematic("Dummy",typed);
        sysd->fOverallUp   = 0.;
        sysd->fOverallDown = -0.;
        sysd->fScaleUp   = 1.;
        sysd->fScaleDown   = 1.;
        fSystematics.push_back( sysd );
        TtHFitter::SYSTMAP[sysd->fName] = "Dummy";
        fNSyst++;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sam = fSamples[i_smp];
            if(sam->fType == Sample::SIGNAL ) {
                sam->AddSystematic(sysd);
            }
        }
    }

    while(true){
        cs = fConfig->GetConfigSet("Systematic",nSys);
        if(cs==0x0) break;
        nSys++;
        if(onlySystematics.size()>0 && FindInStringVector(onlySystematics,cs->GetValue())<0) continue;
        if(toExclude.size()>0 && FindInStringVector(toExclude,cs->GetValue())>=0) continue;
        string samples_str = cs->Get("Samples");
        string regions_str = cs->Get("Regions");
        string exclude_str = cs->Get("Exclude");
        string excludeRegionSample_str = cs->Get("ExcludeRegionSample");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        vector<string> samples = Vectorize(samples_str,',');
        vector<string> regions = Vectorize(regions_str,',');
        vector<string> exclude = Vectorize(exclude_str,',');
        vector<string> excludeRegionSample = Vectorize(excludeRegionSample_str,',');
        type = Systematic::HISTO;
        if(cs->Get("Type")=="overall" || cs->Get("Type")=="OVERALL")
            type = Systematic::OVERALL;
        if(cs->Get("Type")=="shape" || cs->Get("Type")=="SHAPE")
            type = Systematic::SHAPE;
        if(cs->Get("Type")=="stat" || cs->Get("Type")=="STAT") {
            type = Systematic::STAT;
        }
        string decorrelate = cs->Get("Decorrelate");
        sys = new Systematic(CheckName(cs->GetValue()),type);
        TtHFitter::SYSTMAP[sys->fName] = sys->fTitle;
        if(cs->Get("Type")=="overall" || cs->Get("Type")=="OVERALL")
            sys->fIsNormOnly=true;
        if(cs->Get("Category")!=""){
            sys->fCategory = cs->Get("Category");
        }
        // Experimental
        if(cs->Get("IsFreeParameter")!="" && cs->Get("IsFreeParameter")!="FALSE" && cs->Get("IsFreeParameter")!="False" && cs->Get("IsFreeParameter")!="false"){
            sys->fIsFreeParameter = true;
        }
        // New: name to use when writing / reading the Histograms file
        if(cs->Get("StoredName")!=""){
            sys->fStoredName = cs->Get("StoredName");
        }
        //
        bool hasUp   = false;
        bool hasDown = false;
        if(type==Systematic::HISTO || type==Systematic::SHAPE){
            if(fInputType==0){
                if(cs->Get("HistoPathUp")!="")      { sys->fHistoPathsUp  .push_back(cs->Get("HistoPathUp"));   hasUp   = true; }
                if(cs->Get("HistoPathDown")!="")    { sys->fHistoPathsDown.push_back(cs->Get("HistoPathDown")); hasDown = true; }
                if(cs->Get("HistoPathSufUp")!="")   { sys->fHistoPathSufUp   = cs->Get("HistoPathSufUp");       hasUp   = true; }
                if(cs->Get("HistoPathSufDown")!="") { sys->fHistoPathSufDown = cs->Get("HistoPathSufDown");     hasDown = true; }
                if(cs->Get("HistoFileUp")!="")      { sys->fHistoFilesUp  .push_back(cs->Get("HistoFileUp"));   hasUp   = true; }
                if(cs->Get("HistoFileDown")!="")    { sys->fHistoFilesDown.push_back(cs->Get("HistoFileDown")); hasDown = true; }
                if(cs->Get("HistoFileSufUp")!="")   { sys->fHistoFileSufUp   = cs->Get("HistoFileSufUp");       hasUp   = true; }
                if(cs->Get("HistoFileSufDown")!="") { sys->fHistoFileSufDown = cs->Get("HistoFileSufDown");     hasDown = true; }
                if(cs->Get("HistoNameUp")!="")      { sys->fHistoNamesUp  .push_back(cs->Get("HistoNameUp"));   hasUp   = true; }
                if(cs->Get("HistoNameDown")!="")    { sys->fHistoNamesDown.push_back(cs->Get("HistoNameDown")); hasDown = true; }
                if(cs->Get("HistoNameSufUp")!="")   { sys->fHistoNameSufUp   = cs->Get("HistoNameSufUp");       hasUp   = true; }
                if(cs->Get("HistoNameSufDown")!="") { sys->fHistoNameSufDown = cs->Get("HistoNameSufDown");     hasDown = true; }
                // ...
            }
            else if(fInputType==1){
                if(cs->Get("NtuplePathUp")!="")      { sys->fNtuplePathsUp  .push_back(cs->Get("NtuplePathUp"));    hasUp   = true; }
                if(cs->Get("NtuplePathDown")!="")    { sys->fNtuplePathsDown.push_back(cs->Get("NtuplePathDown"));  hasDown = true; }
                if(cs->Get("NtuplePathsUp")!="")     { sys->fNtuplePathsUp   = Vectorize(cs->Get("NtuplePathsUp")  ,',');  hasUp   = true; }
                if(cs->Get("NtuplePathsDown")!="")   { sys->fNtuplePathsDown = Vectorize(cs->Get("NtuplePathsDown"),',');  hasDown = true; }
                if(cs->Get("NtuplePathSufUp")!="")   { sys->fNtuplePathSufUp   = cs->Get("NtuplePathSufUp");         hasUp   = true; }
                if(cs->Get("NtuplePathSufDown")!="") { sys->fNtuplePathSufDown = cs->Get("NtuplePathSufDown");       hasDown = true; }
                if(cs->Get("NtupleFileUp")!="")      { sys->fNtupleFilesUp   .push_back(cs->Get("NtupleFileUp"));    hasUp   = true; }
                if(cs->Get("NtupleFileDown")!="")    { sys->fNtupleFilesDown .push_back(cs->Get("NtupleFileDown"));  hasDown = true; }
                if(cs->Get("NtupleFilesUp")!="")     { sys->fNtupleFilesUp     = Vectorize(cs->Get("NtupleFilesUp"),  ',');  hasUp   = true; }
                if(cs->Get("NtupleFilesDown")!="")   { sys->fNtupleFilesDown   = Vectorize(cs->Get("NtupleFilesDown"),',');  hasDown = true; }
                if(cs->Get("NtupleFileSufUp")!="")   { sys->fNtupleFileSufUp   = cs->Get("NtupleFileSufUp");         hasUp   = true; }
                if(cs->Get("NtupleFileSufDown")!="") { sys->fNtupleFileSufDown = cs->Get("NtupleFileSufDown");       hasDown = true; }
                if(cs->Get("NtupleNameUp")!="")      { sys->fNtupleNamesUp  .push_back(cs->Get("NtupleNameUp"));     hasUp   = true; }
                if(cs->Get("NtupleNameDown")!="")    { sys->fNtupleNamesDown.push_back( cs->Get("NtupleNameDown"));  hasDown = true; }
                if(cs->Get("NtupleNamesUp")!="")     { sys->fNtupleNamesUp     = Vectorize(cs->Get("NtupleNamesUp"),  ','); hasUp   = true; }
                if(cs->Get("NtupleNamesDown")!="")   { sys->fNtupleNamesDown   = Vectorize(cs->Get("NtupleNamesDown"),','); hasDown = true; }
                if(cs->Get("NtupleNameSufUp")!="")   { sys->fNtupleNameSufUp   = cs->Get("NtupleNameSufUp");         hasUp   = true; }
                if(cs->Get("NtupleNameSufDown")!="") { sys->fNtupleNameSufDown = cs->Get("NtupleNameSufDown");       hasDown = true; }
                if(cs->Get("WeightUp")!="")          { sys->fWeightUp      = cs->Get("WeightUp");                    hasUp   = true; }
                if(cs->Get("WeightDown")!="")        { sys->fWeightDown    = cs->Get("WeightDown");                  hasDown = true; }
                if(cs->Get("WeightSufUp")!="")       { sys->fWeightSufUp   = cs->Get("WeightSufUp");                 hasUp   = true; }
                if(cs->Get("WeightSufDown")!="")     { sys->fWeightSufDown = cs->Get("WeightSufDown");               hasDown = true; }
                if(cs->Get("IgnoreWeight")!="")        sys->fIgnoreWeight  = cs->Get("IgnoreWeight");
                // ...
            }
            sys->fHasUpVariation   = hasUp  ;
            sys->fHasDownVariation = hasDown;
            //
            if(cs->Get("Symmetrisation")!=""){
                if(cs->Get("Symmetrisation")=="OneSided" || cs->Get("Symmetrisation")=="ONESIDED")
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEONESIDED;
                else if(cs->Get("Symmetrisation")=="TwoSided" || cs->Get("Symmetrisation")=="TWOSIDED")
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZETWOSIDED;
                else
                    std::cout << "Symetrisation scheme is not recognized ... " << std::endl;
            }
            if(cs->Get("Smoothing")!=""){
                sys->fSmoothType = atoi(cs->Get("Smoothing").c_str());
            }
            param = cs->Get("PreSmoothing");
            if(param!=""){
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if(param=="TRUE") sys->fPreSmoothing = true;
            }
            // ...
        }
        else if(type==Systematic::OVERALL){
            sys->fOverallUp   = atof( cs->Get("OverallUp").c_str() );
            sys->fOverallDown = atof( cs->Get("OverallDown").c_str() );
        }
        param = cs->Get("ScaleUp"); if(param!="")   {
            std::vector < std::string > temp_vec = Vectorize(param,',');
            if(temp_vec.size()==1 && Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleUp = atof(param.c_str());
            }
            else{
                for(unsigned int i_reg = 0; i_reg < temp_vec.size(); i_reg++){
                    std::vector < std::string > reg_value = Vectorize(temp_vec[i_reg],':');
                    if(reg_value.size()==2){
                        sys->fScaleUpRegions.insert( std::pair < std::string, float >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }
        param = cs->Get("ScaleDown"); if(param!="")   {
            std::vector < std::string > temp_vec = Vectorize(param,',');
            if(temp_vec.size()==1 && Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleDown = atof(param.c_str());
            }
            else{
                for(unsigned int i_reg = 0; i_reg < temp_vec.size(); i_reg++){
                    std::vector < std::string > reg_value = Vectorize(temp_vec[i_reg],':');
                    if(reg_value.size()==2){
                        sys->fScaleDownRegions.insert( std::pair < std::string, float >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }
        // a new way of defining systematics, just comparing directly with GHOST samples
        // --> this can be used only if this systematic is applied to a single sample
        param = cs->Get("SampleUp"); if(param!="")   {
            sys->fSampleUp = param;
        }
        param = cs->Get("SampleDown"); if(param!="")   {
            sys->fSampleDown = param;
        }
//         if(cs->Get("ScaleUp")  !="") sys->fScaleUp     = atof( cs->Get("ScaleUp")  .c_str() );
//         if(cs->Get("ScaleDown")!="") sys->fScaleDown   = atof( cs->Get("ScaleDown").c_str() );
        // this to obtain syst variation relatively to given sample
        param = cs->Get("ReferenceSample"); if(param!="") sys->fReferenceSample = param;
        param = cs->Get("KeepReferenceOverallVar");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "FALSE") sys->fKeepReferenceOverallVar = false;
            if(param == "TRUE" ) sys->fKeepReferenceOverallVar = true;
        }
        param = cs->Get("DropShapeIn"); if(param!=""){
            sys->fDropShapeIn = Vectorize(param,',');
        }
        param = cs->Get("DropNorm"); if(param!=""){
            sys->fDropNormIn = Vectorize(param,',');
        }
        param = cs->Get("KeepNormForSamples"); if(param!=""){
            sys->fKeepNormForSamples = Vectorize(param,',');
        }
        //
        if(regions[0]!="all") sys->fRegions = regions;
        if(exclude[0]!="")    sys->fExclude = exclude;
        for (unsigned int iERS=0; iERS<excludeRegionSample.size(); iERS++){
          vector<string> pair_ERS = Vectorize(excludeRegionSample.at(iERS),':');
          if (pair_ERS.size()==2){
            sys->fExcludeRegionSample.push_back(pair_ERS);
          }
        }
        //
        // save list of...
        //
        // Default
        if ( decorrelate == "" && type != Systematic::STAT) {
            fSystematics.push_back( sys );
            fNSyst++;
            if(cs->Get("NuisanceParameter")!=""){
                sys->fNuisanceParameter = cs->Get("NuisanceParameter");
                TtHFitter::NPMAP[sys->fName] = sys->fNuisanceParameter;
            }
            else{
                sys->fNuisanceParameter = sys->fName;
                TtHFitter::NPMAP[sys->fName] = sys->fName;
            }
            if(cs->Get("Title")!=""){
                sys->fTitle = cs->Get("Title");
                TtHFitter::SYSTMAP[sys->fName] = sys->fTitle;
            }
            param = cs->Get("TexTitle"); if(param!=""){
                TtHFitter::SYSTTEX[sys->fName] = param;
            }
            //
            // attach the syst to the proper samples
            for(int i_smp=0;i_smp<fNSamples;i_smp++){
                sam = fSamples[i_smp];
//                 if(sam->fType == Sample::DATA) continue;
                // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
                if(sam->fType == Sample::DATA){
                  if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
                    sam->AddSystematic(sys);
                  }
                  else continue;
                }
                if(!sam->fUseSystematics) continue;
                if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
                      && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
                ){
                    sam->AddSystematic(sys);
                }
            }
        }
        //
        // Decorrelate by region
        else if (decorrelate == "REGION" || type == Systematic::STAT)  {
            //
            // looping over the regions
            for(unsigned int i_reg=0; i_reg<regNames.size(); ++i_reg) {
                bool keepReg=false;
                if ( regions[0]=="all" ) keepReg=true;
                else {
                    for ( int i_GoodReg=0; i_GoodReg < (int)regions.size(); ++i_GoodReg ) {
                        if ( regNames[i_reg]==regions[i_GoodReg] ) keepReg=true;
                    }
                }
                for ( int i_BadReg=0; i_BadReg < (int)exclude.size(); ++i_BadReg ) {
                    if ( exclude[i_BadReg]==regNames[i_reg] ) keepReg=false;
                }

                if (!keepReg) {
                    std::cout << " IGNORING REGION: " << regNames[i_reg] << std::endl;
                    continue;
                }
                std::cout << " --> KEEPING IT!!! " << regNames[i_reg] << std::endl;

                if (type == Systematic::STAT) {
                  Region* reg = GetRegion(regNames[i_reg]);
                  unsigned int nbins = reg->fHistoNBinsRebin>0 ? reg->fHistoNBinsRebin : reg->fNbins;
                    std::cout << regNames[i_reg] << " " << nbins << std::endl;
                    // decorrelate by bin
                    for (int i_bin = 0; i_bin < nbins; i_bin++) {
                        Systematic* mySys= new Systematic(*sys);
                        mySys->fName=(mySys->fName)+"_"+regNames[i_reg]+"_bin"+std::to_string(i_bin);
                        std::vector<string> tmpReg;
                        tmpReg.push_back( regNames[i_reg] );
                        mySys->fRegions = tmpReg;
                        std::vector<int> tmpBin;
                        tmpBin.push_back( i_bin );
                        mySys->fBins = tmpBin;
                        fSystematics.push_back( mySys );
                        if(cs->Get("NuisanceParameter")!=""){
                            mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+regNames[i_reg]+"_bin"+std::to_string(i_bin);
                            TtHFitter::NPMAP[mySys->fName] = sys->fNuisanceParameter;
                        } else {
                            mySys->fNuisanceParameter = mySys->fName;
                            TtHFitter::NPMAP[mySys->fName] = mySys->fName;
                        }
                        if(cs->Get("Title")!=""){
                            mySys->fTitle = (sys->fTitle)+"_"+regNames[i_reg]+"_bin"+std::to_string(i_bin);
                            TtHFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
                        }
                        fNSyst++;
                        for (int i_smp=0;i_smp<fNSamples;i_smp++){
                            sam = fSamples[i_smp];
                            if(sam->fType == Sample::DATA) continue;
                            if(!sam->fUseSystematics) continue;
                            if((samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
                                && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
                                ){
                                sam->AddSystematic(mySys);
                            }
                        }
                    }
                } else {
                    //
                    // cloning the sys for each region
                    Systematic* mySys= new Systematic(*sys);
                    mySys->fName=(mySys->fName)+"_"+regNames[i_reg];
                    std::vector<string> tmpReg;
                    tmpReg.push_back( regNames[i_reg] );
                    mySys->fRegions = tmpReg;
                    fSystematics.push_back( mySys );
                    if(cs->Get("NuisanceParameter")!=""){
                        mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+regNames[i_reg];
                        TtHFitter::NPMAP[mySys->fName] = sys->fNuisanceParameter;
                    }
                    else{
                        mySys->fNuisanceParameter = mySys->fName;
                        TtHFitter::NPMAP[mySys->fName] = mySys->fName;
                    }
                    if(cs->Get("Title")!=""){
                        mySys->fTitle = (sys->fTitle)+"_"+regNames[i_reg];
                        TtHFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
                    }
                    fNSyst++;
                    //
                    for(int i_smp=0;i_smp<fNSamples;i_smp++){
                        sam = fSamples[i_smp];
                        // if(sam->fType == Sample::DATA) continue;
                        // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
                        if(sam->fType == Sample::DATA){
                          if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
                            sam->AddSystematic(mySys);
                          }
                          else continue;
                        }
                        if(!sam->fUseSystematics) continue;
                        if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
                              && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
                              ){
                          sam->AddSystematic(mySys);
                        }
                    }
                }
            }
            delete sys;
        }
        //
        // Decorrelate by Sample
        else if (decorrelate == "SAMPLE")  {
            //
            // (this is really messy)
            for(int i_smp=0;i_smp<fNSamples;i_smp++){
                sam = fSamples[i_smp];
//                 if(sam->fType == Sample::DATA) continue;
                // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
                if(sam->fType == Sample::DATA){
                  if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
                    sam->AddSystematic(sys);
                  }
                  else continue;
                }
                if(sam->fType == Sample::GHOST) continue;
                bool keepSam=false;
                if ( samples[0]=="all" ) keepSam=true;
                else {
                    if ( find(samples.begin(), samples.end(), sam->fName)!=samples.end() ) keepSam=true;
                }
                if ( find(exclude.begin(), exclude.end(), sam->fName)!=exclude.end() ) keepSam=false;
                if (!keepSam) {
                    std::cout << " IGNORING SAMPLE: " << sam->fName << std::endl;
                    continue;
                }
                std::cout << " --> KEEPING SAMPLE: " << sam->fName << std::endl;
                //
                // cloning the sys for each region
                Systematic* mySys= new Systematic(*sys);
                mySys->fName=(mySys->fName)+"_"+sam->fName;
                fSystematics.push_back( mySys );
                if(cs->Get("NuisanceParameter")!=""){
                    mySys->fNuisanceParameter = (sys->fNuisanceParameter)+"_"+sam->fName;
                    TtHFitter::NPMAP[mySys->fName] = sys->fNuisanceParameter;
                }
                else{
                    mySys->fNuisanceParameter = mySys->fName;
                    TtHFitter::NPMAP[mySys->fName] = mySys->fName;
                }
                if(cs->Get("Title")!=""){
                    mySys->fTitle = (sys->fTitle)+"_"+sam->fName;
                    TtHFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
                }
                fNSyst++;
                sam->AddSystematic(mySys);
            }
            delete sys;
        }
        //
        // Decorrelate  shape and normalisation
        else if (decorrelate == "SHAPEACC")  {
            //
            // cloning the sys
            Systematic* mySys1= new Systematic(*sys);
            mySys1->fName=(mySys1->fName)+"_Acc";
            fSystematics.push_back( mySys1 );
            mySys1->fIsNormOnly=true;
            if(cs->Get("NuisanceParameter")!=""){
                mySys1->fNuisanceParameter = (sys->fNuisanceParameter)+"_Acc";
                TtHFitter::NPMAP[mySys1->fName] = sys->fNuisanceParameter;
            }
            else{
                mySys1->fNuisanceParameter = mySys1->fName;
                TtHFitter::NPMAP[mySys1->fName] = mySys1->fName;
            }
            if(cs->Get("Title")!=""){
                mySys1->fTitle = (sys->fTitle)+"_Acc";
                TtHFitter::SYSTMAP[mySys1->fName] = mySys1->fTitle;
            }
            fNSyst++;
            //
            for(int i_smp=0;i_smp<fNSamples;i_smp++){
                sam = fSamples[i_smp];
//                 if(sam->fType == Sample::DATA) continue;
                // in principle, no syst on DATA, except if this syst has SubtractRefSampleVar: TRUE and this data sample is the ReferenceSample of that syst
                if(sam->fType == Sample::DATA){
                  if (sys->fSubtractRefSampleVar && sys->fReferenceSample == sam->fName) {
                    sam->AddSystematic(mySys1);
                  }
                  else continue;
                }
                if(!sam->fUseSystematics) continue;
                if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
                      && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
                      ){
                    sam->AddSystematic(mySys1);
                }
            }
            //
            if ( sys->fType!=Systematic::OVERALL ) {
                // cloning the sys
                Systematic* mySys2= new Systematic(*sys);
                mySys2->fName=(mySys2->fName)+"_Shape";
                mySys2->fIsNormOnly=false;
                mySys2->fIsShapeOnly=true;
                fSystematics.push_back( mySys2 );
                //
                if(cs->Get("NuisanceParameter")!=""){
                    mySys2->fNuisanceParameter = (sys->fNuisanceParameter)+"_Shape";
                    TtHFitter::NPMAP[mySys2->fName] = sys->fNuisanceParameter;
                }
                else{
                    mySys2->fNuisanceParameter = mySys2->fName;
                    TtHFitter::NPMAP[mySys2->fName] = mySys2->fName;
                }
                if(cs->Get("Title")!=""){
                    mySys2->fTitle = (sys->fTitle)+"_Shape";
                    TtHFitter::SYSTMAP[mySys2->fName] = mySys2->fTitle;
                }
                fNSyst++;
                //
                for(int i_smp=0;i_smp<fNSamples;i_smp++){
                    sam = fSamples[i_smp];
                    if(sam->fType == Sample::DATA) continue;
                    if(!sam->fUseSystematics) continue;
                    if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
                          && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
                          ){
                        sam->AddSystematic(mySys2);
                    }
                }
            }
            delete sys;
        }
        //
        // If not....
        else {
            std::cout << "decorrelate option: " << decorrelate  << "  not supported ... PLEASE USE ONLY: REGION, SAMPLE, SHAPEACC" <<  std::endl;
            return;
        }
        // New: for systeamtics which also vary Data (e.g. JER with Full NPs)
        // This will subtract linearly the relative variation on Data from each relative variation on MC
        param = cs->Get("SubtractRefSampleVar");
        if(param!=""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE" ) sys->fSubtractRefSampleVar = true;// default is false
        }
    }
    //
    // -- Post config-reading actions ---
    //
    // if StatOnly, also sets to OFF the MC stat
    if(fStatOnly){
        std::cout << "TtHFit::INFO: StatOnly option is setting to OFF the MC-stat (gammas) as well.\nTo keep them on use the command line option 'Systematics=NONE' or comment out all Systematics in config file." << std::endl;
        SetStatErrorConfig( false, 0. );
    }
}


//__________________________________________________________________________________
// for each region, add a SampleHist for each Sample in the Fit, reading from ntuples
void TtHFit::ReadNtuples(){
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Reading ntuples..." << std::endl;
    TH1F* h = 0x0;
    TH1F* hUp = 0x0;
    TH1F* hDown = 0x0;
    TH1F* htmp = 0x0;
    string fullSelection;
    string fullMCweight;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();
    SampleHist *sh;
    //
    // Loop on regions and samples
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        std::cout << "  Region " << fRegions[i_ch]->fName << " ..." << std::endl;
        //
        if(TtHFitter::SPLITHISTOFILES) fFiles[i_ch]->cd();
        //
        if(fRegions[i_ch]->fBinTransfo != "")  ComputeBining(i_ch);
        if(fRegions[i_ch]->fCorrVar1 != ""){
            if(fRegions[i_ch]->fCorrVar2 == ""){
                std::cout << "TtHFitter::WARNING : Only first correlation variable defined, do not read region : " << fRegions[i_ch]->fName << std::endl;
                continue;
            }
            std::cout << "calling the function" << std::endl;
            defineVariable(i_ch);
        }
        else if(fRegions[i_ch]->fCorrVar2 != ""){
            std::cout << "TtHFitter::WARNING : Only second correlation variable defined, do not read region : " << fRegions[i_ch]->fName << std::endl;
            continue;
        }

        // first loop on Data samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType!=Sample::DATA) continue;
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "  Reading DATA sample " << fSamples[i_smp]->fName << std::endl;
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            // set selection and weight
            fullSelection = "1";
            //             fSelection + " && " + fRegions[i_ch]->fSelection;
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSelection!="" && fSelection!="1")
                fullSelection += " && ("+fSelection+")";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fRegions[i_ch]->fSelection!="" && fRegions[i_ch]->fSelection!="1")
                fullSelection += " && ("+fRegions[i_ch]->fSelection+")";
            if(fSamples[i_smp]->fSelection!="" && fSamples[i_smp]->fSelection!="1")
                fullSelection += " && ("+fSamples[i_smp]->fSelection+")";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSamples[i_smp]->fIgnoreSelection!="FALSE" && fSamples[i_smp]->fIgnoreSelection!="")
                fullSelection = ReplaceString(fullSelection,fSamples[i_smp]->fIgnoreSelection,"1");
            //
            fullMCweight = "1";
            //
            //Check whether to use an alternative variable (e.g. in TRF) or not
            std::string variable;
            if (fRegions[i_ch]->UseAlternativeVariable(fSamples[i_smp]->fName)){
              variable = fRegions[i_ch]->fAlternativeVariables[fSamples[i_smp]->fName];
            }
            else{
              variable = fRegions[i_ch]->fVariable;
            }
            //
            // build a list of ntuples to read
            fullPaths.clear();
            vector<string> NtupleNames;
            for(unsigned int ns_ch=0; ns_ch<fRegions[i_ch]->fNtupleNames.size(); ++ns_ch){
                NtupleNames.push_back(fRegions[i_ch]->fNtupleNames.at(ns_ch));
            }
            for(unsigned int ns_smp=0; ns_smp<fSamples[i_smp]->fNtupleNames.size(); ++ns_smp){
                NtupleNames.push_back(fSamples[i_smp]->fNtupleNames.at(ns_smp));
            }
            vector<string> NtupleNameSuffs = CombinePathSufs( fSamples[i_smp]->fNtupleNameSuffs,
                                                             fRegions[i_ch]->fNtupleNameSuffs );
            fullPaths = CreatePathsList( fSamples[i_smp]->fNtuplePaths.size()>0 ? fSamples[i_smp]->fNtuplePaths : fNtuplePaths,
                                         fRegions[i_ch]->fNtuplePathSuffs,
                                         fSamples[i_smp]->fNtupleFiles.size()>0 ? fSamples[i_smp]->fNtupleFiles : ToVec(fNtupleFile), empty, // no ntuple file suffs for nominal (syst only)
                                         NtupleNames.size()>0 ? NtupleNames : ToVec( fNtupleName ),
                                         NtupleNameSuffs.size()>0 ? NtupleNameSuffs : empty  // NEW
                                         );
            htmp = 0x0;
            h = 0x0;
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                if(fRegions[i_ch]->fHistoBins){
                    htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                 variable, fRegions[i_ch]->fHistoNBinsRebin, fRegions[i_ch]->fHistoBins,
                                                 fullSelection, fullMCweight);
                }
                else{
                    htmp = HistFromNtuple( fullPaths[i_path],
                                           variable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                           fullSelection, fullMCweight);
                    //Pre-processing of histograms (rebinning, lumi scaling)
                    if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                }
                //
                if(fSamples[i_smp]->fLumiScales.size()>i_path)  htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                // scale stat uncertainty if sample is data-driven!
                if(!fSamples[i_smp]->fNormalizedByTheory){
                    float scale = 1.;
                    if(fSamples[i_smp]->fLumiScales.size()>i_path)  scale = fSamples[i_smp]->fLumiScales[i_path];
                    else if(fSamples[i_smp]->fLumiScales.size()==1) scale = fSamples[i_smp]->fLumiScales[0];
                    for(int i_bin=1;i_bin<=htmp->GetNbinsX();i_bin++){
                        htmp->SetBinError( i_bin, htmp->GetBinError( i_bin ) / sqrt(scale) );
                    }
                }
                //
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );
            //
            // Importing the histogram in TtHFitter
            sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
            sh->fHist_orig = h_orig;
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name
            
            // in fact DATA can be used for systs that have SubtractRefSampleVar: TRUE
            // so we need to get its systematics first
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                Systematic * syst = fSamples[i_smp]->fSystematics[i_syst];
                // only relevant for systs that have this sample as reference
                if (!syst->fSubtractRefSampleVar || syst->fReferenceSample != fSamples[i_smp]->fName) continue;

                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,fRegions[i_ch]->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,fRegions[i_ch]->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding DATA syst " << syst->fName << std::endl;
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                //
                // set selection
                fullSelection = "1";
                if(smp->fIgnoreSelection!="TRUE" && fSelection!="" && fSelection!="1")
                    fullSelection += " && ("+fSelection+")";
                if(smp->fIgnoreSelection!="TRUE" && reg->fSelection!="" && reg->fSelection!="1")
                    fullSelection += " && ("+reg->fSelection+")";
                if(smp->fSelection!="" && smp->fSelection!="1")
                    fullSelection += " && ("+smp->fSelection+")";
                if(smp->fIgnoreSelection!="TRUE" && smp->fIgnoreSelection!="FALSE" && smp->fIgnoreSelection!="")
                    fullSelection = ReplaceString(fullSelection,smp->fIgnoreSelection,"1");

                //
                // Up
                //
                hUp = 0x0;
                if(syst->fHasUpVariation){
                    fullMCweight = "1.";
                    fullPaths.clear();
                    vector<string> NtupleNameSuffsUp = CombinePathSufs( ToVec( syst->fNtupleNameSufUp ), reg->fNtupleNameSuffs );
                    vector<string> NtuplePaths       = fNtuplePaths;
                      if(smp->fNtuplePaths.size()>0)    NtuplePaths = smp->fNtuplePaths;
                      if(syst->fNtuplePathsUp.size()>0) NtuplePaths = syst->fNtuplePathsUp;
                    vector<string> NtuplePathSuffs   = CombinePathSufs( reg->fNtuplePathSuffs, ToVec( syst->fNtuplePathSufUp ) );
                    //
                    fullPaths = CreatePathsList(
                                                // path
                                                NtuplePaths,
                                                // path suf
                                                NtuplePathSuffs,
                                                // file
                                                syst->fNtupleFilesUp.size()==0 ?
                                                ( smp->fNtupleFiles.size()>0 ? smp->fNtupleFiles : ToVec(fNtupleFile) ) :
                                                syst->fNtupleFilesUp ,
                                                // file suf
                                                syst->fNtupleFileSufUp=="" ?
                                                empty :
                                                ToVec( syst->fNtupleFileSufUp ),
                                                // name
                                                syst->fNtupleNamesUp.size()==0 ?
                                                ( smp->fNtupleNames.size()==0 ? ToVec( fNtupleName ) : smp->fNtupleNames ) :
                                                syst->fNtupleNamesUp,
                                                // name suf
                                                NtupleNameSuffsUp.size()>0 ? NtupleNameSuffsUp : empty
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        if(reg->fHistoBins){
                            htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                        variable, reg->fHistoNBinsRebin, reg->fHistoBins,
                                                        fullSelection, fullMCweight);
                        }
                        else{
                            htmp = HistFromNtuple( fullPaths[i_path],
                                                  variable, reg->fNbins, reg->fXmin, reg->fXmax,
                                                  fullSelection, fullMCweight);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        }
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hUp = (TH1F*)htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hUp->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                // Down
                //
                hDown = 0x0;
                if(syst->fHasDownVariation){
                    fullMCweight = "1.";
                    fullPaths.clear();
                    vector<string> NtupleNameSuffsDown  = CombinePathSufs( ToVec( syst->fNtupleNameSufDown ), reg->fNtupleNameSuffs );
                    vector<string> NtuplePaths          = fNtuplePaths;
                      if(smp->fNtuplePaths.size()>0)      NtuplePaths = smp->fNtuplePaths;
                      if(syst->fNtuplePathsDown.size()>0) NtuplePaths = syst->fNtuplePathsDown;
                    vector<string> NtuplePathSuffs      = CombinePathSufs( reg->fNtuplePathSuffs, ToVec( syst->fNtuplePathSufDown ) );
                    //
                    fullPaths = CreatePathsList(
                                                // path
                                                NtuplePaths,
                                                // path suf
                                                NtuplePathSuffs,
                                                // file
                                                syst->fNtupleFilesDown.size()==0 ?
                                                ( smp->fNtupleFiles.size()>0 ? smp->fNtupleFiles : ToVec(fNtupleFile) ) :
                                                syst->fNtupleFilesDown ,
                                                // file suf
                                                syst->fNtupleFileSufDown=="" ?
                                                empty :
                                                ToVec( syst->fNtupleFileSufDown ),
                                                // name
                                                syst->fNtupleNamesDown.size()==0 ?
                                                ( smp->fNtupleNames.size()==0 ? ToVec( fNtupleName ) : smp->fNtupleNames ) :
                                                syst->fNtupleNamesDown,
                                                // name suf
                                                NtupleNameSuffsDown.size()>0 ? NtupleNameSuffsDown : empty
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        if(reg->fHistoBins){
                            htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                        variable, reg->fHistoNBinsRebin, reg->fHistoBins,
                                                        fullSelection, fullMCweight);
                        }
                        else{
                            htmp = HistFromNtuple( fullPaths[i_path],
                                                  variable, reg->fNbins, reg->fXmin, reg->fXmax,
                                                  fullSelection, fullMCweight);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        }
                        //
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hDown = (TH1F*)htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hDown->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                if(hUp==0x0)   hUp   = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                if(hDown==0x0) hDown = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
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
        }
        // Then loop again on non-data samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "  Reading " << fSamples[i_smp]->fName << std::endl;
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            // set selection and weight
            fullSelection = "1";
            //             fSelection + " && " + fRegions[i_ch]->fSelection;
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSelection!="" && fSelection!="1")
                fullSelection += " && ("+fSelection+")";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fRegions[i_ch]->fSelection!="" && fRegions[i_ch]->fSelection!="1")
                fullSelection += " && ("+fRegions[i_ch]->fSelection+")";
            if(fSamples[i_smp]->fSelection!="" && fSamples[i_smp]->fSelection!="1")
                fullSelection += " && ("+fSamples[i_smp]->fSelection+")";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSamples[i_smp]->fIgnoreSelection!="FALSE" && fSamples[i_smp]->fIgnoreSelection!="")
                fullSelection = ReplaceString(fullSelection,fSamples[i_smp]->fIgnoreSelection,"1");
            //
            if(!fSamples[i_smp]->fNormalizedByTheory){ // for data-driven bkg, use just the sample weight (FIXME)
                fullMCweight = fSamples[i_smp]->fMCweight;
            }
            else{
                fullMCweight = fMCweight + " * " + fSamples[i_smp]->fMCweight;
                if(fRegions[i_ch]->fMCweight!="") fullMCweight += " * " + fRegions[i_ch]->fMCweight;
            }
            //
            //Check whether to use an alternative variable (e.g. in TRF) or not
            std::string variable;
            if (fRegions[i_ch]->UseAlternativeVariable(fSamples[i_smp]->fName)){
              variable = fRegions[i_ch]->fAlternativeVariables[fSamples[i_smp]->fName];
            }
            else{
              variable = fRegions[i_ch]->fVariable;
            }
            //
            // build a list of ntuples to read
            fullPaths.clear();
            vector<string> NtupleNames;
            for(unsigned int ns_ch=0; ns_ch<fRegions[i_ch]->fNtupleNames.size(); ++ns_ch){
                NtupleNames.push_back(fRegions[i_ch]->fNtupleNames.at(ns_ch));
            }
            for(unsigned int ns_smp=0; ns_smp<fSamples[i_smp]->fNtupleNames.size(); ++ns_smp){
                NtupleNames.push_back(fSamples[i_smp]->fNtupleNames.at(ns_smp));
            }
            vector<string> NtupleNameSuffs = CombinePathSufs( fSamples[i_smp]->fNtupleNameSuffs,
                                                             fRegions[i_ch]->fNtupleNameSuffs );
            fullPaths = CreatePathsList( fSamples[i_smp]->fNtuplePaths.size()>0 ? fSamples[i_smp]->fNtuplePaths : fNtuplePaths,
                                         fRegions[i_ch]->fNtuplePathSuffs,
                                         fSamples[i_smp]->fNtupleFiles.size()>0 ? fSamples[i_smp]->fNtupleFiles : ToVec(fNtupleFile), empty, // no ntuple file suffs for nominal (syst only)
                                         NtupleNames.size()>0 ? NtupleNames : ToVec( fNtupleName ),
                                         NtupleNameSuffs.size()>0 ? NtupleNameSuffs : empty  // NEW
                                         );
            htmp = 0x0;
            h = 0x0;
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                if(fRegions[i_ch]->fHistoBins){
                    htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                 variable, fRegions[i_ch]->fHistoNBinsRebin, fRegions[i_ch]->fHistoBins,
                                                 fullSelection, fullMCweight);
                }
                else{
                    htmp = HistFromNtuple( fullPaths[i_path],
                                           variable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                           fullSelection, fullMCweight);
                    //Pre-processing of histograms (rebinning, lumi scaling)
                    if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                }
                //
                if(fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                //
                if(fSamples[i_smp]->fLumiScales.size()>i_path)  htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                //
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );
            //
            // Importing the histogram in TtHFitter
            sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
            sh->fHist_orig = h_orig;
            sh->fHist_orig->SetName( Form("%s_orig",sh->fHist->GetName()) ); // fix the name

            // end here if no systematics allowed (e.g. generally for GHOST)
            if(!fSamples[i_smp]->fUseSystematics) continue;

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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding norm " << nf->fName << std::endl;
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding shape " << sf->fName << std::endl;
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding syst " << syst->fName << std::endl;
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL){
//                     SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddOverallSyst(syst->fName,syst->fOverallUp*syst->fScaleUp,syst->fOverallDown*syst->fScaleDown);
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
                //
                if(syst->fReferenceSample!="") smp = GetSample(syst->fReferenceSample);
                //
                // set selection
                fullSelection = "1";
                if(smp->fIgnoreSelection!="TRUE" && fSelection!="" && fSelection!="1")
                    fullSelection += " && ("+fSelection+")";
                if(smp->fIgnoreSelection!="TRUE" && reg->fSelection!="" && reg->fSelection!="1")
                    fullSelection += " && ("+reg->fSelection+")";
                if(smp->fSelection!="" && smp->fSelection!="1")
                    fullSelection += " && ("+smp->fSelection+")";
                if(smp->fIgnoreSelection!="TRUE" && smp->fIgnoreSelection!="FALSE" && smp->fIgnoreSelection!="")
                    fullSelection = ReplaceString(fullSelection,smp->fIgnoreSelection,"1");
                //
                // Up
                //
                hUp = 0x0;
                if(syst->fHasUpVariation){
                    // Note: no need to change selection for systematics. If needed, can be done via weight (partially...)
                    fullMCweight = "1.";
                    if(fSamples[i_smp]->fNormalizedByTheory)  fullMCweight = fMCweight;
                    if(syst->fWeightUp!="")
                        fullMCweight += " * "+syst->fWeightUp;
                    else
                        fullMCweight += " * "+smp->fMCweight;
                    if(syst->fWeightSufUp!="")
                        fullMCweight += " * "+syst->fWeightSufUp;
                    if(reg->fMCweight!="" && fSamples[i_smp]->fNormalizedByTheory)
                        fullMCweight += " * "+reg->fMCweight;
                    if(syst->fIgnoreWeight!=""){
                        ReplaceString(fullMCweight, syst->fIgnoreWeight,"");
                        ReplaceString(fullMCweight,"*  *","*");
                        ReplaceString(fullMCweight,"* *","*");
                        ReplaceString(fullMCweight,"**","*");
                    }
                    if(TtHFitter::DEBUGLEVEL>0) std::cout << "  Syst Up full weight: " << fullMCweight << std::endl;
                    //
                    fullPaths.clear();
                    vector<string> NtupleNameSuffsUp = CombinePathSufs( ToVec( syst->fNtupleNameSufUp ), reg->fNtupleNameSuffs );
                    vector<string> NtuplePaths       = fNtuplePaths;
                      if(smp->fNtuplePaths.size()>0)    NtuplePaths = smp->fNtuplePaths;
                      if(syst->fNtuplePathsUp.size()>0) NtuplePaths = syst->fNtuplePathsUp;
                    vector<string> NtuplePathSuffs   = CombinePathSufs( reg->fNtuplePathSuffs, ToVec( syst->fNtuplePathSufUp ) );
                    //
                    fullPaths = CreatePathsList(
                                                // path
                                                NtuplePaths,
                                                // path suf
                                                NtuplePathSuffs,
                                                // file
                                                syst->fNtupleFilesUp.size()==0 ?
                                                ( smp->fNtupleFiles.size()>0 ? smp->fNtupleFiles : ToVec(fNtupleFile) ) :
                                                syst->fNtupleFilesUp ,
                                                // file suf
                                                syst->fNtupleFileSufUp=="" ?
                                                empty :
                                                ToVec( syst->fNtupleFileSufUp ),
                                                // name
                                                syst->fNtupleNamesUp.size()==0 ?
                                                ( smp->fNtupleNames.size()==0 ? ToVec( fNtupleName ) : smp->fNtupleNames ) :
                                                syst->fNtupleNamesUp,
                                                // name suf
                                                NtupleNameSuffsUp.size()>0 ? NtupleNameSuffsUp : empty
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        if(reg->fHistoBins){
                            htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                        variable, reg->fHistoNBinsRebin, reg->fHistoBins,
                                                        fullSelection, fullMCweight);
                        }
                        else{
                            htmp = HistFromNtuple( fullPaths[i_path],
                                                  variable, reg->fNbins, reg->fXmin, reg->fXmax,
                                                  fullSelection, fullMCweight);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        }
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        //
                        // obtain relative variation and apply it to proper sample
                        // & try to keep also the same total relative variation
                        if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            //
                            float relVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                            htmp->Divide(   href );
                            htmp->Multiply( hnom );
                            float newVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                            if( syst->fKeepReferenceOverallVar && TMath::Abs(relVar-1) > 0.0001 && TMath::Abs(newVar) > 0.0001) htmp->Scale( relVar / newVar );
                        }
                        // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                        else if (syst->fReferenceSample!="" && syst->fSubtractRefSampleVar) {
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* href_up = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistUp;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            
                            // Formula: UpHisto = [1+(up-nom)/nom-(DataUp-Data)/Data]*nom = up+nom+DataUp/Data*nom
                            TH1* href_up_Tmp = (TH1*)href_up->Clone(Form("%s_Tmp", href_up->GetName()));
                            href_up_Tmp->Divide(href);
                            href_up_Tmp->Multiply(hnom);
                            htmp->Add(hnom);
                            htmp->Add(href_up_Tmp,-1);
                            
                            delete href_up_Tmp;// it's a clone, and it's the purpose of clones to die
                        }
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hUp = (TH1F*)htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hUp->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                // Down
                //
                hDown = 0x0;
                if(syst->fHasDownVariation){
                    //
                    // Note: no need to change selection for systematics. If needed, can be done via weight (partially...)
                    fullMCweight = "1.";
                    if(fSamples[i_smp]->fNormalizedByTheory)  fullMCweight = fMCweight;
                    if(syst->fWeightDown!="")
                        fullMCweight += " * "+syst->fWeightDown;
                    else
                        fullMCweight += " * " +smp->fMCweight;
                    if(syst->fWeightSufDown!="")
                        fullMCweight += " * "+syst->fWeightSufDown;
                    if(reg->fMCweight!="" && fSamples[i_smp]->fNormalizedByTheory)
                        fullMCweight += " * "+reg->fMCweight;
                    if(syst->fIgnoreWeight!=""){
                        ReplaceString(fullMCweight, syst->fIgnoreWeight,"");
                        ReplaceString(fullMCweight,"*  *","*");
                        ReplaceString(fullMCweight,"* *","*");
                        ReplaceString(fullMCweight,"**","*");
                    }
                    //
                    fullPaths.clear();
                    vector<string> NtupleNameSuffsDown  = CombinePathSufs( ToVec( syst->fNtupleNameSufDown ), reg->fNtupleNameSuffs );
                    vector<string> NtuplePaths          = fNtuplePaths;
                      if(smp->fNtuplePaths.size()>0)      NtuplePaths = smp->fNtuplePaths;
                      if(syst->fNtuplePathsDown.size()>0) NtuplePaths = syst->fNtuplePathsDown;
                    vector<string> NtuplePathSuffs      = CombinePathSufs( reg->fNtuplePathSuffs, ToVec( syst->fNtuplePathSufDown ) );
                    //
                    fullPaths = CreatePathsList(
                                                // path
                                                NtuplePaths,
                                                // path suf
                                                NtuplePathSuffs,
                                                // file
                                                syst->fNtupleFilesDown.size()==0 ?
                                                ( smp->fNtupleFiles.size()>0 ? smp->fNtupleFiles : ToVec(fNtupleFile) ) :
                                                syst->fNtupleFilesDown ,
                                                // file suf
                                                syst->fNtupleFileSufDown=="" ?
                                                empty :
                                                ToVec( syst->fNtupleFileSufDown ),
                                                // name
                                                syst->fNtupleNamesDown.size()==0 ?
                                                ( smp->fNtupleNames.size()==0 ? ToVec( fNtupleName ) : smp->fNtupleNames ) :
                                                syst->fNtupleNamesDown,
                                                // name suf
                                                NtupleNameSuffsDown.size()>0 ? NtupleNameSuffsDown : empty
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        if(reg->fHistoBins){
                            htmp = HistFromNtupleBinArr( fullPaths[i_path],
                                                        variable, reg->fHistoNBinsRebin, reg->fHistoBins,
                                                        fullSelection, fullMCweight);
                        }
                        else{
                            htmp = HistFromNtuple( fullPaths[i_path],
                                                  variable, reg->fNbins, reg->fXmin, reg->fXmax,
                                                  fullSelection, fullMCweight);
                            // Pre-processing of histograms (rebinning, lumi scaling)
                            if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        }
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        //
                        // obtain relative variation and apply it to proper sample
                        // & try to keep also the same total relative variation
                        if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            //
                            float relVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                            htmp->Divide(   href );
                            htmp->Multiply( hnom );
                            float newVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                            if( syst->fKeepReferenceOverallVar && TMath::Abs(relVar-1) > 0.0001 && TMath::Abs(newVar-1) > 0.0001) htmp->Scale( relVar / newVar );
                        }
                        // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                        else if (syst->fReferenceSample!="" && syst->fSubtractRefSampleVar) {
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* href_down = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistDown;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            
                            // Formula: UpHisto = [1+(down-nom)/nom-(DataDown-Data)/Data]*nom = down+nom+DataDown/Data*nom
                            TH1* href_down_Tmp = (TH1*) href_down->Clone(Form("%s_Tmp", href_down->GetName()));
                            href_down_Tmp->Divide(href);
                            href_down_Tmp->Multiply(hnom);
                            htmp->Add(hnom);
                            htmp->Add(href_down_Tmp,-1);
                            
                            delete href_down_Tmp;// it's a clone, and it's the purpose of clones to die
                        }
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hDown = (TH1F*)htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hDown->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                if(hUp==0x0)   hUp   = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                if(hDown==0x0) hDown = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
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
    delete htmp;
}

//__________________________________________________________________________________
// this method takes care of rebinning, smoothing, fixing
void TtHFit::CorrectHistograms(){
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        Region *reg = fRegions[i_ch];
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            Sample *smp = fSamples[i_smp];
            //
            // eventually skip sample / region combination
            if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            SampleHist *sh = reg->GetSampleHist(smp->fName);
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            TH1* h_orig = (TH1*)sh->fHist_orig;
            TH1* h      = (TH1*)h_orig->Clone(sh->fHist->GetName());
            sh->fHist = h;

            //
            // Systematics
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                Systematic * syst = smp->fSystematics[i_syst];
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                // get the original syst histograms & reset the syst histograms
                SystematicHist *syh = sh->GetSystematic( syst->fName );
                //
                if(syh==0x0) continue;
                TH1* hUp_orig   = syh->fHistUp_orig;
                TH1* hDown_orig = syh->fHistDown_orig;
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL){
                    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
                        hUp_orig  ->SetBinContent(i_bin,h_orig->GetBinContent(i_bin)*(1.+syst->fOverallUp));
                        hDown_orig->SetBinContent(i_bin,h_orig->GetBinContent(i_bin)*(1.+syst->fOverallDown));
                    }
                }
                //
                TH1* hUp   = (TH1*)hUp_orig  ->Clone( syh->fHistUp  ->GetName() );
                TH1* hDown = (TH1*)hDown_orig->Clone( syh->fHistDown->GetName() );
                syh->fHistUp   = hUp;
                syh->fHistDown = hDown;
            }

            // ---> NEED TO MOVE TO READNTUPLES? -- FIXME
            // Subtraction / Addition of sample
            for(auto sample : fSamples[i_smp]->fSubtractSamples){
                std::cout << "INFO: subtracting sample " << sample << std::endl;
                SampleHist *smph0 = fRegions[i_ch]->GetSampleHist(sample);
                sh->Add(smph0,-1);
            }
            for(auto sample : fSamples[i_smp]->fAddSamples){
                std::cout << "INFO: adding sample " << sample << std::endl;
                SampleHist *smph0 = fRegions[i_ch]->GetSampleHist(sample);
                sh->Add(smph0);
            }
            // Division & Multiplication by other samples
            if(fSamples[i_smp]->fMultiplyBy!=""){
                std::cout << "INFO: multiplying by sample " << fSamples[i_smp]->fMultiplyBy << std::endl;
                SampleHist *smph0 = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fMultiplyBy);
                sh->Multiply(smph0);
            }
            if(fSamples[i_smp]->fDivideBy!=""){
                std::cout << "INFO: dividing by sample " << fSamples[i_smp]->fDivideBy << std::endl;
                SampleHist *smph0 = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fDivideBy);
                sh->Divide(smph0);
            }

            //
            // For SampleUp / SampleDown
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                Systematic * syst = smp->fSystematics[i_syst];
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                // get the original syst histograms & reset the syst histograms
                SystematicHist *syh = sh->GetSystematic( syst->fName );
                //
                // if syst defined with SampleUp / SampleDown
                if( syst->fSampleUp != "" || syst->fSampleDown != "" ){
                    TH1 *h_up   = syst->fSampleUp  !="" ? reg->GetSampleHist(syst->fSampleUp)  ->fHist : sh->fHist;
                    TH1 *h_down = syst->fSampleDown!="" ? reg->GetSampleHist(syst->fSampleDown)->fHist : sh->fHist;
                    syh = sh->AddHistoSyst(syst->fName,h_up,h_down);
                    syh->fSystematic = syst;
                }
            }
            
            //
            // Fix empty bins
            if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::SIGNAL){
                sh->FixEmptyBins();
            }

            //
            // Eventually smooth nominal histogram  (use with caution...)
            TH1* h_correction = 0x0;
            bool isFlat = false;
            if(smp->fSmooth && !reg->fSkipSmoothing){
//                 h_correction = SmoothHistogram( h );
                h_correction = (TH1*)h->Clone( Form("%s_corr",h->GetName()) );
                TH1* h0 = (TH1*)h->Clone( Form("%s_orig0",h->GetName()) );
                isFlat = SmoothHistogram( h );
                h_correction->Divide( h0 );
            }

            //
            // Systematics
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                Systematic * syst = smp->fSystematics[i_syst];
                if(syst==0x0) continue;
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                //
                SystematicHist *syh = sh->GetSystematic( syst->fName );
                if(syh==0x0) continue;
                TH1* hUp   = syh->fHistUp;
                TH1* hDown = syh->fHistDown;
                //
                // if Overall only, re-create it if smoothing was applied
                if(syst->fType==Systematic::OVERALL){
                    if(h_correction!=0x0 && fSamples[i_smp]->fSmooth){
                        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
                            hUp  ->SetBinContent(i_bin,h->GetBinContent(i_bin)*(1.+syst->fOverallUp));
                            hDown->SetBinContent(i_bin,h->GetBinContent(i_bin)*(1.+syst->fOverallDown));
                        }
                    }
                    continue;
                }
                //
                // correct according to the sample nominal smoothing
                if(h_correction!=0x0 && fSamples[i_smp]->fSmooth){
//                     if(hUp!=0x0  ) hUp  ->Multiply(h_correction);
//                     if(hDown!=0x0) hDown->Multiply(h_correction);
                    if(hUp!=0x0  ) SmoothHistogram( hUp  , isFlat );
                    if(hDown!=0x0) SmoothHistogram( hDown, isFlat );
                }

                //
                // Histogram smoothing, Symmetrisation, Massaging...
                if(!fRegions[i_ch]->fSkipSmoothing) syh -> fSmoothType = syst -> fSmoothType;
                else                                syh -> fSmoothType = 0;
                syh -> fSymmetrisationType = syst -> fSymmetrisationType;

            }  // end syst loop
            //
            // Histograms checking
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                Systematic * syst = smp->fSystematics[i_syst];
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                if( sh->GetSystematic( smp->fSystematics[i_syst]->fName )==0x0 ) continue;
                //
                HistoTools::CheckHistograms( fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist /*nominal*/,
                                            sh->GetSystematic( smp->fSystematics[i_syst]->fName ) /*systematic*/,
                                            fSamples[i_smp]->fType!=Sample::SIGNAL/*check bins with content=0*/,
                                            TtHFitter::HISTOCHECKCRASH /*cause crash if problem*/);
            }

            // set the fill color
            sh->fHist->SetFillColor(fillcolor);
            sh->fHist->SetLineColor(linecolor);
        } // end sample loop
        //
    } // end region loop
    //
    // Smooth systematics
    SmoothSystematics("all");
    
//     // NEW: scale systematics according to ScaleUp and ScaleDown
//     for(int i_ch=0;i_ch<fNRegions;i_ch++){
//         Region *reg = fRegions[i_ch];
//         //
//         for(int i_smp=0;i_smp<fNSamples;i_smp++){
//             Sample *smp = fSamples[i_smp];
//             //
//             // eventually skip sample / region combination
//             if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
//             //
//             SampleHist *sh = reg->GetSampleHist(smp->fName);
//             if(sh==0x0) continue;
//             //
//             // Systematics
//             for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
//                 Systematic *syst = smp->fSystematics[i_syst];
//                 if(syst==0x0) continue;
//                 //
//                 // eventually skip systematic / region combination
//                 if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
//                 if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
//                 if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
//                 //
//                 SystematicHist *syh = sh->GetSystematic( syst->fName );
//                 if(syh==0x0) continue;
//                 //
//                 if(syh->fScaleUp!=1.){
//                     TH1* h_tmp = (TH1*)syh->fHistUp->Clone();
//                     h_tmp->Add(sh->fHist,-1);
//                     h_tmp->Scale(syh->fScaleUp);
//                     h_tmp->Add(sh->fHist);
//                     syh->fHistUp = h_tmp;
//                     cout << syh->fNormUp << " ";
//                     syh->fNormUp *= syh->fScaleUp;
//                     cout << syh->fNormUp << endl;
//                 }
//                 //
//                 if(syh->fScaleDown!=1.){
//                     TH1* h_tmp = (TH1*)syh->fHistDown->Clone();
//                     h_tmp->Add(sh->fHist,-1);
//                     h_tmp->Scale(syh->fScaleDown);
//                     h_tmp->Add(sh->fHist);
//                     syh->fHistDown = h_tmp;
//                     syh->fNormDown *= syh->fScaleDown;
//                 }
//             }            
//         }
//     }
    
    
    //
    // NEW: artifificially set all systematics not to affect overall normalisation for sample or set of samples
    // (the form should be KeepNormForSamples: ttlight+ttc+ttb,wjets
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        Region *reg = fRegions[i_ch];
        for(int i_sys=0;i_sys<fNSyst;i_sys++){
            if(fSystematics[i_sys]->fKeepNormForSamples.size()==0) continue;
            for(int ii=0;ii<fSystematics[i_sys]->fKeepNormForSamples.size();ii++){
                std::vector<std::string> subSamples = Vectorize(fSystematics[i_sys]->fKeepNormForSamples[ii],'+');
                // get nominal yield and syst yields for this sum of samples
                float yieldNominal = 0.;
                float yieldUp = 0.;
                float yieldDown = 0.;
                for(int i_smp=0;i_smp<fNSamples;i_smp++){
                    if(FindInStringVector(subSamples,fSamples[i_smp]->fName)<0) continue;
                    SampleHist *sh = reg->GetSampleHist(fSamples[i_smp]->fName);
                    if(sh==0x0) continue;
                    SystematicHist *syh = sh->GetSystematic(fSystematics[i_sys]->fName);
                    if(syh==0x0) continue;
                    yieldNominal += sh ->fHist    ->Integral();
                    yieldUp      += syh->fHistUp  ->Integral();
                    yieldDown    += syh->fHistDown->Integral();
                }
                // scale each syst variation
                for(int i_smp=0;i_smp<fNSamples;i_smp++){
                    if(FindInStringVector(subSamples,fSamples[i_smp]->fName)<0) continue;
                    SampleHist *sh = reg->GetSampleHist(fSamples[i_smp]->fName);
                    if(sh==0x0) continue;
                    SystematicHist *syh = sh->GetSystematic(fSystematics[i_sys]->fName);
                    if(syh==0x0) continue;
                    if(TtHFitter::DEBUGLEVEL>0){
                        std::cout << "  Normalising syst " << fSystematics[i_sys]->fName << " for sample " << fSamples[i_smp]->fName;
                        std::cout << "\t scaling by " << yieldNominal/yieldUp << " (up), " << yieldNominal/yieldDown << " (down)" << std::endl;
                    }
                    syh->fHistUp  ->Scale(yieldNominal/yieldUp);
                    syh->fHistDown->Scale(yieldNominal/yieldDown);
                }
            }
        }
    }
    
    //
    // For data: restore original data sample (needed if bins were dropped previously)
    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }
    if(hasData){
        for(int i_ch=0;i_ch<fNRegions;i_ch++){
            Region *reg = fRegions[i_ch];
            if(reg->fData!=0x0){
                if(reg->fData->fHist_orig!=0x0 && reg->fData->fHist!=0x0) reg->fData->fHist = (TH1F*)reg->fData->fHist_orig->Clone(reg->fData->fHist->GetName());
            }
        }
    }
    
    //
    // Drop bins
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        Region *reg = fRegions[i_ch];
        if(reg->fDropBins.size()!=0){
            for(int i_smp=0;i_smp<fNSamples;i_smp++){
                Sample *smp = fSamples[i_smp];
                // eventually skip sample / region combination
                if( FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                SampleHist *sh = reg->GetSampleHist(smp->fName);
                if(sh==0x0) continue;
                DropBins(sh->fHist,reg->fDropBins);
                for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                    Systematic *syst = smp->fSystematics[i_syst];
                    // eventually skip systematic / region combination
                    if( syst->fRegions.size()>0 && FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                    if( syst->fExclude.size()>0 && FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                    if( syst->fExcludeRegionSample.size()>0 && FindInStringVectorOfVectors(syst->fExcludeRegionSample,fRegions[i_ch]->fName, fSamples[i_smp]->fName)>=0 ) continue;
                    SystematicHist *syh = sh->GetSystematic( syst->fName );
                    if(syh==0x0) continue;
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
void TtHFit::ReadHistograms(){
    TH1F* h = 0x0;
    TH1F* hUp = 0x0;
    TH1F* hDown = 0x0;
    TH1F* htmp = 0x0;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();
    SampleHist *sh;
 
    //
    // Loop on regions and samples
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        std::cout << "  Region " << fRegions[i_ch]->fName << " ..." << std::endl;
        //
        if(TtHFitter::SPLITHISTOFILES) fFiles[i_ch]->cd();
        //
        if(fRegions[i_ch]->fBinTransfo != "") ComputeBining(i_ch);
        // first we must read the DATA samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType!=Sample::DATA) continue;
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "  Reading DATA sample " << fSamples[i_smp]->fName << std::endl;
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            // build a list of histograms to read
            fullPaths.clear();
            std::vector<string> histoFiles;
            std::vector<string> histoNames;
            if(fSamples[i_smp]->fHistoFiles.size()>0)     histoFiles = fSamples[i_smp]->fHistoFiles;
            else if(fRegions[i_ch]->fHistoFiles.size()>0) histoFiles = fRegions[i_ch]->fHistoFiles;
            else                                          histoFiles = ToVec( fHistoFile );
            if(fSamples[i_smp]->fHistoNames.size()>0)     histoNames = fSamples[i_smp]->fHistoNames;
            else if(fRegions[i_ch]->fHistoNames.size()>0) histoNames = fRegions[i_ch]->fHistoNames;
            else                                          histoNames = ToVec( fHistoName );

            fullPaths = CreatePathsList( fHistoPaths, CombinePathSufs(fRegions[i_ch]->fHistoPathSuffs, fSamples[i_smp]->fHistoPaths),
                                         histoFiles, empty, // no histo file suffs for nominal (syst only)
                                         histoNames, empty  // same for histo name
                                        );
            htmp = 0x0;
            h = 0x0;
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                //Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[i_ch]->fHistoBins){
                    TH1F* htmp2 = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,"htmp2",fRegions[i_ch]->fHistoBins));
                    const char *hname = htmp->GetName();
                    htmp->~TH1F();
                    htmp = (TH1F*)htmp2->Clone();
                    delete htmp2;
                    htmp->SetName(hname);
                }
                else if(fRegions[i_ch]->fHistoNBinsRebin != -1) {
                    htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                }
                if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);

                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );
            //
            // Importing the histogram in TtHFitter
            sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding syst " << syst->fName << std::endl;
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                //
                // Up
                //
                hUp = 0x0;
                if(syst->fHasUpVariation){
                    fullPaths.clear();
                    fullPaths = CreatePathsList(
                                                // path
                                                fHistoPaths,
                                                // path suf
                                                CombinePathSufs(reg->fHistoPathSuffs,syst->fHistoPathsUp ),
                                                // file
                                                syst->fHistoFilesUp.size()==0 ?
                                                histoFiles :
                                                syst->fHistoFilesUp ,
                                                // file suf
                                                syst->fHistoFileSufUp=="" ?
                                                empty :
                                                ToVec( syst->fHistoFileSufUp ),
                                                // name
                                                syst->fHistoNamesUp.size()==0 ?
                                                histoNames :
                                                syst->fHistoNamesUp,
                                                // name suf
                                                syst->fHistoNameSufUp=="" ?
                                                empty :
                                                ToVec( syst->fHistoNameSufUp )
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                        // Pre-processing of histograms (rebinning, lumi scaling)
                        if(reg->fHistoBins){
                            TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                            const char *hname = htmp->GetName();
                            htmp->~TH1F();
                            htmp = (TH1F*)htmp2->Clone();
                            delete htmp2;
                            htmp->SetName(hname);
                        }
                        else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        //
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hUp = (TH1F*)htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hUp->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                // Down
                //
                hDown = 0x0;
                if(syst->fHasDownVariation){
                    fullPaths.clear();
                    fullPaths = CreatePathsList(
                                                // path
                                                fHistoPaths,
                                                // path suf
                                                CombinePathSufs(reg->fHistoPathSuffs, syst->fHistoPathsDown ),
                                                // file
                                                syst->fHistoFilesDown.size()==0 ?
                                                histoFiles :
                                                syst->fHistoFilesDown ,
                                                // file suf
                                                syst->fHistoFileSufDown=="" ?
                                                empty :
                                                ToVec( syst->fHistoFileSufDown ),
                                                // name
                                                syst->fHistoNamesDown.size()==0 ?
                                                histoNames :
                                                syst->fHistoNamesDown,
                                                // name suf
                                                syst->fHistoNameSufDown=="" ?
                                                empty :
                                                ToVec( syst->fHistoNameSufDown )
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        htmp = (TH1F*)HistFromFile( fullPaths[i_path] ) ;
                        // Pre-processing of histograms (rebinning, lumi scaling)
                        if(reg->fHistoBins){
                            TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                            const char *hname = htmp->GetName();
                            htmp->~TH1F();
                            htmp = (TH1F*)htmp2->Clone();
                            delete htmp2;
                            htmp->SetName(hname);
                        }
                        else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        //
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hDown = (TH1F*)htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hDown->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                if(hUp==0x0)   hUp   = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                if(hDown==0x0) hDown = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
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
        }

        // then we can read the other samples
        std::set < std::string > files_names;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "  Reading " << fSamples[i_smp]->fName << std::endl;
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[i_ch]->fName)<0 ) continue;
            //
            // read nominal
            //
            // build a list of histograms to read
            fullPaths.clear();
            std::vector<string> histoFiles;
            std::vector<string> histoNames;
            if(fSamples[i_smp]->fHistoFiles.size()>0)     histoFiles = fSamples[i_smp]->fHistoFiles;
            else if(fRegions[i_ch]->fHistoFiles.size()>0) histoFiles = fRegions[i_ch]->fHistoFiles;
            else                                          histoFiles = ToVec( fHistoFile );
            if(fSamples[i_smp]->fHistoNames.size()>0)     histoNames = fSamples[i_smp]->fHistoNames;
            else if(fRegions[i_ch]->fHistoNames.size()>0) histoNames = fRegions[i_ch]->fHistoNames;
            else                                          histoNames = ToVec( fHistoName );

            fullPaths = CreatePathsList( fHistoPaths, CombinePathSufs(fRegions[i_ch]->fHistoPathSuffs, fSamples[i_smp]->fHistoPaths),
                                         histoFiles, empty, // no histo file suffs for nominal (syst only)
                                         histoNames, empty  // same for histo name
                                        );
            htmp = 0x0;
            h = 0x0;
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                files_names.insert(fullPaths[i_path]);
                htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                //Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[i_ch]->fHistoBins){
                    TH1F* htmp2 = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,"htmp2",fRegions[i_ch]->fHistoBins));
                    const char *hname = htmp->GetName();
                    htmp->~TH1F();
                    htmp = (TH1F*)htmp2->Clone();
                    delete htmp2;
                    htmp->SetName(hname);
                }
                else if(fRegions[i_ch]->fHistoNBinsRebin != -1) {
                    htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                }

                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);

                if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);

                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            //
            // Save the original histogram
            TH1* h_orig = (TH1*)h->Clone( Form("%s_orig",h->GetName()) );

            //
            // Importing the histogram in TtHFitter
            sh = fRegions[i_ch]->SetSampleHist( fSamples[i_smp], h );
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding norm " << nf->fName << std::endl;
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding shape " << sf->fName << std::endl;
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding syst " << syst->fName << std::endl;
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                //
                // if Overall only ...
                if(syst->fType==Systematic::OVERALL){
//                     SystematicHist *syh = reg->GetSampleHist(smp->fName)->AddOverallSyst(syst->fName,syst->fOverallUp*syst->fScaleUp,syst->fOverallDown*syst->fScaleDown);
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
                hUp = 0x0;
                if(syst->fHasUpVariation){
                    fullPaths.clear();
                    fullPaths = CreatePathsList(
                                                // path
                                                fHistoPaths,
                                                // path suf
                                                CombinePathSufs(reg->fHistoPathSuffs,syst->fHistoPathsUp ),
                                                // file
                                                syst->fHistoFilesUp.size()==0 ?
                                                histoFiles :
                                                syst->fHistoFilesUp ,
                                                // file suf
                                                syst->fHistoFileSufUp=="" ?
                                                empty :
                                                ToVec( syst->fHistoFileSufUp ),
                                                // name
                                                syst->fHistoNamesUp.size()==0 ?
                                                histoNames :
                                                syst->fHistoNamesUp,
                                                // name suf
                                                syst->fHistoNameSufUp=="" ?
                                                empty :
                                                ToVec( syst->fHistoNameSufUp )
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        files_names.insert(fullPaths[i_path]);
                        htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                        // Pre-processing of histograms (rebinning, lumi scaling)
                        if(reg->fHistoBins){
                            TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                            const char *hname = htmp->GetName();
                            htmp->~TH1F();
                            htmp = (TH1F*)htmp2->Clone();
                            delete htmp2;
                            htmp->SetName(hname);
                        }
                        else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        //
                        // obtain relative variation and apply it to proper sample
                        // & try to keep also the same total relative variation
                        if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            //
                            float relVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                            htmp->Divide(   href );
                            htmp->Multiply( hnom );
                            float newVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                            if( syst->fKeepReferenceOverallVar && TMath::Abs(relVar-1) > 0.0001 && TMath::Abs(newVar-1) > 0.0001) htmp->Scale( relVar / newVar );
                        }
                        // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                        else if (syst->fReferenceSample!="" && syst->fSubtractRefSampleVar) {
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* href_up = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistUp;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            
                            // Formula: UpHisto = [1+(up-nom)/nom-(DataUp-Data)/Data]*nom = up+nom+DataUp/Data*nom
                            TH1* href_up_Tmp = (TH1*) href_up->Clone(Form("%s_Tmp", href_up->GetName()));
                            href_up_Tmp->Divide(href);
                            href_up_Tmp->Multiply(hnom);
                            htmp->Add(hnom);
                            htmp->Add(href_up_Tmp,-1);
                            
                            delete href_up_Tmp;// it's a clone, and it's the purpose of clones to die
                        }
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hUp = (TH1F*)htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hUp->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                // Down
                //
                hDown = 0x0;
                if(syst->fHasDownVariation){
                    fullPaths.clear();
                    fullPaths = CreatePathsList(
                                                // path
                                                fHistoPaths,
                                                // path suf
                                                CombinePathSufs(reg->fHistoPathSuffs, syst->fHistoPathsDown ),
                                                // file
                                                syst->fHistoFilesDown.size()==0 ?
                                                histoFiles :
                                                syst->fHistoFilesDown ,
                                                // file suf
                                                syst->fHistoFileSufDown=="" ?
                                                empty :
                                                ToVec( syst->fHistoFileSufDown ),
                                                // name
                                                syst->fHistoNamesDown.size()==0 ?
                                                histoNames :
                                                syst->fHistoNamesDown,
                                                // name suf
                                                syst->fHistoNameSufDown=="" ?
                                                empty :
                                                ToVec( syst->fHistoNameSufDown )
                                                );
                    for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                        files_names.insert(fullPaths[i_path]);
                        htmp = (TH1F*)HistFromFile( fullPaths[i_path] ) ;
                        // Pre-processing of histograms (rebinning, lumi scaling)
                        if(reg->fHistoBins){
                            TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                            const char *hname = htmp->GetName();
                            htmp->~TH1F();
                            htmp = (TH1F*)htmp2->Clone();
                            delete htmp2;
                            htmp->SetName(hname);
                        }
                        else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                        //
                        if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                        if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                        else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                        //
                        // obtain relative variation and apply it to proper sample
                        // & try to keep also the same total relative variation
                        if(syst->fReferenceSample!="" && !syst->fSubtractRefSampleVar){
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            //
                            float relVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                            htmp->Divide(   href );
                            htmp->Multiply( hnom );
                            float newVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                            if( syst->fKeepReferenceOverallVar && TMath::Abs(relVar-1) > 0.0001 && TMath::Abs(newVar-1) > 0.0001) htmp->Scale( relVar / newVar );
                        }
                        // new special case: we subtract from the relative uncertainty the relative uncertainty of another (data) sample
                        else if (syst->fReferenceSample!="" && syst->fSubtractRefSampleVar) {
                            TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                            TH1* href_down = reg->GetSampleHist(syst->fReferenceSample)->GetSystematic(syst->fName)->fHistDown;
                            TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<htmp->GetNbinsX()+2;i_bin++) if(htmp->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++) if(href->GetBinContent(i_bin)<=1e-6) htmp->SetBinContent(i_bin,1e-6); // this to avoid multiplying bins by 1e6
                            
                            // Formula: UpHisto = [1+(down-nom)/nom-(DataDown-Data)/Data]*nom = down+nom+DataDown/Data*nom
                            TH1* href_down_Tmp = (TH1*) href_down->Clone(Form("%s_Tmp", href_down->GetName()));
                            href_down_Tmp->Divide(href);
                            href_down_Tmp->Multiply(hnom);
                            htmp->Add(hnom);
                            htmp->Add(href_down_Tmp,-1);
                            
                            delete href_down_Tmp;// it's a clone, and it's the purpose of clones to die
                        }
                        //
                        // Importing histogram in TtHFitter
                        if(i_path==0){
                            hDown = (TH1F*)htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()));
                        }
                        else hDown->Add(htmp);
                        htmp->~TH1F();
                    }
                }
                //
                if(hUp==0x0)   hUp   = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                if(hDown==0x0) hDown = (TH1F*)reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
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
    delete htmp;
}

//__________________________________________________________________________________
//
void TtHFit::ReadHistos(/*string fileName*/){
    string fileName = "";
    TH1F* h;
    SampleHist *sh;
    SystematicHist *syh;
    string regionName;
    string sampleName;
    string normName;
    string shapeName;
    string systName;
    //
    bool singleOutputFile = !TtHFitter::SPLITHISTOFILES;
    if(singleOutputFile){
        if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos.root";
        else                 fileName = fName + "/Histograms/" + fInputName + "_histos.root";
        std::cout << "-----------------------------" << std::endl;
        std::cout << "Reading histograms from file " << fileName << " ..." << std::endl;
    }
    //
    vector< TH2F* > histPrun;
    TFile *filePrun = 0x0;
    if( fKeepPruning ){
//         filePrun = new TFile( (fName+"/Pruning"+fSuffix+".root").c_str() );
        filePrun = new TFile( (fName+"/Pruning.root").c_str() );
        if(!filePrun) fKeepPruning = false;
    }
    //
    // when we multply/divide by or subtract/add other samples, need to add systematics on the other samples
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fDivideBy!=""){
            Sample* smp = GetSample(fSamples[i_smp]->fDivideBy);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systName = smp->fSystematics[i_syst]->fName;
                if(!fSamples[i_smp]->HasSystematic(systName)){
                    fSamples[i_smp]->AddSystematic(smp->fSystematics[i_syst]);
                }
            }
        }
        if(fSamples[i_smp]->fMultiplyBy!=""){
            Sample* smp = GetSample(fSamples[i_smp]->fMultiplyBy);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systName = smp->fSystematics[i_syst]->fName;
                if(!fSamples[i_smp]->HasSystematic(systName)){
                    fSamples[i_smp]->AddSystematic(smp->fSystematics[i_syst]);
                }
            }
        }
        for(auto sample : fSamples[i_smp]->fSubtractSamples){
            Sample* smp = GetSample(sample);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systName = smp->fSystematics[i_syst]->fName;
                if(!fSamples[i_smp]->HasSystematic(systName)){
                    fSamples[i_smp]->AddSystematic(smp->fSystematics[i_syst]);
                }
            }
        }
        for(auto sample : fSamples[i_smp]->fAddSamples){
            Sample* smp = GetSample(sample);
            for(int i_syst=0;i_syst<smp->fNSyst;i_syst++){
                std::string systName = smp->fSystematics[i_syst]->fName;
                if(!fSamples[i_smp]->HasSystematic(systName)){
                    fSamples[i_smp]->AddSystematic(smp->fSystematics[i_syst]);
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
        if(TtHFitter::DEBUGLEVEL>0) std::cout << "  Reading region " << regionName << std::endl;
        //
        if(!singleOutputFile){
            if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + regionName + "_histos.root";
            else                 fileName = fName + "/Histograms/" + fInputName + "_" + regionName + "_histos.root";
            std::cout << "-----------------------------" << std::endl;
            std::cout << "Reading histograms from file " << fileName << " ..." << std::endl;
        }
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            //
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,regionName)<0 && fSamples[i_smp]->fName.find("customAsimov_")==std::string::npos ) continue;
            //
            sampleName = fSamples[i_smp]->fName;
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "    Reading sample " << sampleName << std::endl;
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp],regionName+"_"+sampleName,fileName);
            sh = fRegions[i_ch]->GetSampleHist(sampleName);
            //
            // separate gammas -> Add systematic
            if(fSamples[i_smp]->fSeparateGammas){
                std::string systName = "stat_"+fSamples[i_smp]->fName;
                std::string systStoredName = systName;
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "adding separate gammas as SHAPE systematic " << systName << std::endl;
                SystematicHist *syh = sh->AddHistoSyst(systName,
                                                        Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                        Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                        0
                                                        );
                if(syh==0x0){
                    std::cout << "TtHFit::WARNING: no histogram found for separate gamma, but may be you will create it right now." << std::endl;
                }
                else{
                    Systematic *gamma = 0x0;
                    if(FindInStringVector(fSystematicNames,systName)>=0) gamma = fSystematics[FindInStringVector(fSystematicNames,systName)];  //GetSystematic(systName);
                    if(gamma==0x0) gamma = NewSystematic(systName);
                    if(TtHFitter::DEBUGLEVEL>0) std::cout << "adding separate gammas as SHAPE systematic " << systName << std::endl;
                    gamma->fType = Systematic::SHAPE;
                    gamma->fRegions.clear();
                    gamma->fRegions.push_back(fRegions[i_ch]->fName);
                    syh->fSystematic = gamma;
                    gamma->fNuisanceParameter = gamma->fName;
                    TtHFitter::NPMAP[gamma->fName] = gamma->fNuisanceParameter;
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "      Reading norm " << normName << std::endl;
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
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "      Reading shape " << shapeName << std::endl;
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
                systName              = fSamples[i_smp]->fSystematics[i_syst]->fName;
                string systStoredName = fSamples[i_smp]->fSystematics[i_syst]->fStoredName; // if no StoredName specified in the config, this should be == fName
                //
                // eventually skip systematics if pruned
                int xbin,ybin,bin;
                int binContent = 0;
                if( fKeepPruning && histPrun[i_ch]!=0x0 ){
                    xbin = histPrun[i_ch]->GetXaxis()->FindBin( sampleName.c_str() ); // sample
                    ybin = histPrun[i_ch]->GetYaxis()->FindBin( systName.c_str() ); // syst
                    bin = histPrun[i_ch]->GetBin(xbin,ybin);
                    binContent = histPrun[i_ch]->GetBinContent(bin);
                    if( binContent <= -4 || binContent == -1 || binContent >= 3 ){
                        std::cout << "SKIPPING systematic " << systName << std::endl;
                        continue;
                    }
                    //{kBlack,6,kBlue, kGray, 8, kYellow, kOrange-3, kRed}
                }
                if(TtHFitter::DEBUGLEVEL>0) std::cout << "      Reading syst " << systName << std::endl;
                // norm only
                if(fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    if( fKeepPruning ){
                        if( binContent == -2 || binContent == 2 ) continue;
                    }
                    syh = sh->AddOverallSyst(systName,
//                          fSamples[i_smp]->fSystematics[i_syst]->fOverallUp*fSamples[i_smp]  ->fSystematics[i_syst]->fScaleUp,
//                          fSamples[i_smp]->fSystematics[i_syst]->fOverallDown*fSamples[i_smp]->fSystematics[i_syst]->fScaleDown);
                         fSamples[i_smp]->fSystematics[i_syst]->fOverallUp,
                         fSamples[i_smp]->fSystematics[i_syst]->fOverallDown);
//                     syh->fSystematic = fSamples[i_smp]->fSystematics[i_syst];
//                     syh->fScaleUp = fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
//                     if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0)
//                         if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName]!=0){
//                             syh->fScaleUp *= fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName];
//                         }
//                     syh->fScaleDown = fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
//                     if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0)
//                         if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName]!=0)
//                             syh->fScaleDown *= fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName];
                }
                // histo syst
                else{
                    int pruned = 0;
                    if( fKeepPruning ){
                        if(binContent==1 || binContent==-2) pruned = 1;
                        if(binContent==2 || binContent==-3) pruned = 2;
                    }
                        syh = sh->AddHistoSyst(systName,
                                               Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                               Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                               pruned
                                              );
                    }
                // for both
                    syh->fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                syh->fHistoNameShapeUp   = Form("%s_%s_%s_Shape_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                syh->fHistoNameShapeDown = Form("%s_%s_%s_Shape_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                    syh->fFileNameShapeUp    = fileName;
                    syh->fFileNameShapeDown  = fileName;
                syh->fScaleUp = fSamples[i_smp]->fSystematics[i_syst]->fScaleUp;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions.size()!=0)
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName]!=0){
                        syh->fScaleUp *= fSamples[i_smp]->fSystematics[i_syst]->fScaleUpRegions[regionName];
                    }
                syh->fScaleDown = fSamples[i_smp]->fSystematics[i_syst]->fScaleDown;
                if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions.size()!=0)
                    if(fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName]!=0)
                        syh->fScaleDown *= fSamples[i_smp]->fSystematics[i_syst]->fScaleDownRegions[regionName];
                //
                if(fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    syh->fNormUp   *= syh->fScaleUp;
                    syh->fNormDown *= syh->fScaleDown;
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::CloseInputFiles(){
    //
    // Close all input files
    for(auto it : TtHFitter::TFILEMAP){
        TDirectory *dir = gDirectory;
        TFile *f = it.second;
        if(f!=0x0)
        dir->cd();
        f->Close();
        delete f;
    }
}

//__________________________________________________________________________________
//
void TtHFit::DrawAndSaveAll(string opt){
    TthPlot *p;
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    bool isPostFit = opt.find("post")!=string::npos;
    if(isPostFit){
        if(fFitResultsFile!="")
            ReadFitResults(fFitResultsFile);
        else
            ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
        // scale signal by random number used in the ws creation
        if(fRandomPOISeed>=0){
            gRandom->SetSeed(fRandomPOISeed); float rnd = 1./pow(gRandom->Uniform(0,2),2);
            if(fFitResults->fNuisParIsThere[fPOI]){
                if(TtHFitter::DEBUGLEVEL>0)  std::cout << "Scaling POI " << fPOI << " by the secret random number XXX..." << std::endl;
                NuisParameter* np = fFitResults->fNuisPar[fFitResults->fNuisParIdx[fPOI]];
                np->fFitValue*=rnd;
                np->fPostFitUp*=rnd;
                np->fPostFitDown*=rnd;
            }
        }
    }
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->fUseStatErr = fUseStatErr;
        fRegions[i_ch]->fATLASlabel = fAtlasLabel;
        //
        if(fCustomAsimov!=""){
            std::string name = "customAsimov_"+fCustomAsimov;
            SampleHist* cash = fRegions[i_ch]->GetSampleHist(name);
            if(cash==0x0){
                std::cerr << "TtHFit::DrawAndSaveAll::ERROR: no Custom Asimov " << fCustomAsimov << " available. Taking regular Asimov." << std::endl;
            }
            else{
                if(TtHFitter::DEBUGLEVEL>0){
                    std::cout << "  Adding Custom-Asimov Data: " << cash->fHist->GetName() << std::endl;
                }
                fRegions[i_ch]->fData = cash;
            }
        }
        //
        if(isPostFit){
            gSystem->mkdir( (fName + "/Histograms/").c_str() );
            if(fRegions[i_ch]->fRegionDataType==Region::ASIMOVDATA) p = fRegions[i_ch]->DrawPostFit(fFitResults,opt+" blind");
            else                                                    p = fRegions[i_ch]->DrawPostFit(fFitResults,opt);
            for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
                p->SaveAs(     (fName+"/Plots/"+fRegions[i_ch]->fName+"_postFit"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format] ).c_str());
        }
        else{
            if(fRegions[i_ch]->fRegionDataType==Region::ASIMOVDATA) p = fRegions[i_ch]->DrawPreFit(opt+" blind");
            else                                                    p = fRegions[i_ch]->DrawPreFit(opt);
            for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
                p->SaveAs(     (fName+"/Plots/"+fRegions[i_ch]->fName+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format] ).c_str());
        }
    }
}

//__________________________________________________________________________________
//
TthPlot* TtHFit::DrawSummary(string opt, TthPlot* prefit_plot){
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Building Summary Plot..." << std::endl;
    gSystem->mkdir(fName.c_str(),true);
    const bool isPostFit = opt.find("post")!=string::npos;
    const bool checkVR = opt.find("valid")!=string::npos;
    // build one bin per region
    TH1F* h_data = 0;
    TH1F* h_sig[MAXsamples];
    TH1F* h_bkg[MAXsamples];
    TH1F *h_tot;
    TGraphAsymmErrors *g_err;
    int Nsig = 0;
    int Nbkg = 0;
    //
    string name;
    string title;
    int lineColor;
    int fillColor;
    int lineWidth;
    float integral;
    double intErr; // to store the integral error
    TH1* h; // to store varius histograms temporary
    //
    // Building region - bin correspondence
    //
    std::vector<int> regionVec;   regionVec.clear();
    std::vector<int> divisionVec; divisionVec.clear();
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
    if(Nbin<=0) return 0x0;
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        SampleHist *sh = 0x0;
        name = (fSamples[i_smp]->fName).c_str();
        title = fSamples[i_smp]->fTitle.c_str();
        if(fSamples[i_smp]->fGroup != "") title = fSamples[i_smp]->fGroup.c_str();
        // look for the first SampleHist defined for this sample
        for(int i_ch=0;i_ch<(int)regionVec.size();i_ch++){
            sh = fRegions[regionVec[i_ch]]->GetSampleHist( name );
            if(sh!=0x0) break;
        }
        // skip sample if no SampleHist found
        if(sh==0x0) continue;
        if(sh->fHist==0x0) continue;
        //
        lineColor = sh->fHist->GetLineColor();
        fillColor = sh->fHist->GetFillColor();
        lineWidth = sh->fHist->GetLineWidth();
        //
        if(fSamples[i_smp]->fType==Sample::SIGNAL){
            h_sig[Nsig] = new TH1F(name.c_str(),title.c_str(), Nbin,0,Nbin);
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding Signal: " << h_sig[Nsig]->GetTitle() << std::endl;
            h_sig[Nsig]->SetLineColor(lineColor);
            h_sig[Nsig]->SetFillColor(fillColor);
            h_sig[Nsig]->SetLineWidth(lineWidth);
            for(int i_bin=1;(int)i_bin<=regionVec.size();i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=0x0){
                    if(isPostFit)  h = (TH1F*)sh->fHist_postFit->Clone(); // Michele
                    else           h = (TH1F*)sh->fHist->Clone(); // Michele
                    //
                    if(!isPostFit){
                        // FIXME SF 
                        // scale it according to NormFactors
                        for(unsigned int i_nf=0;i_nf<sh->fSample->fNormFactors.size();i_nf++){
                            h->Scale(sh->fSample->fNormFactors[i_nf]->fNominal);
                            if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit::INFO: Scaling " << sh->fSample->fName << " by " << sh->fSample->fNormFactors[i_nf]->fNominal << std::endl;
                        }
                    }
                    //
                    integral = h->IntegralAndError(1,h->GetNbinsX(),intErr);
                }
                else{
                    integral = 0.;
                    intErr = 0.;
                }
                h_sig[Nsig]->SetBinContent( i_bin,integral );
                h_sig[Nsig]->SetBinError( i_bin,intErr );
            }
            Nsig++;
        }
        else if(fSamples[i_smp]->fType==Sample::BACKGROUND){
            h_bkg[Nbkg] = new TH1F(name.c_str(),title.c_str(), Nbin,0,Nbin);
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding Bkg:    " << h_bkg[Nbkg]->GetTitle() << std::endl;
            h_bkg[Nbkg]->SetLineColor(lineColor);
            h_bkg[Nbkg]->SetFillColor(fillColor);
            h_bkg[Nbkg]->SetLineWidth(lineWidth);
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=0x0){
                    if(isPostFit)  h = (TH1F*)sh->fHist_postFit->Clone(); // Michele
                    else           h = (TH1F*)sh->fHist->Clone(); // Michele
                    //
                    if(!isPostFit){
                        // FIXME SF 
                        // scale it according to NormFactors
                        for(unsigned int i_nf=0;i_nf<sh->fSample->fNormFactors.size();i_nf++){
                            h->Scale(sh->fSample->fNormFactors[i_nf]->fNominal);
                            if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit::INFO: Scaling " << sh->fSample->fName << " by " << sh->fSample->fNormFactors[i_nf]->fNominal << std::endl;
                        }
                    }
                    //
                    integral = h->IntegralAndError(1,h->GetNbinsX(),intErr);
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
            h_data = new TH1F(name.c_str(),title.c_str(), Nbin,0,Nbin);
            if(TtHFitter::DEBUGLEVEL>0) std::cout << "Adding Data:   " << h_data->GetTitle() << std::endl;
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                if(fRegions[regionVec[i_bin-1]]->fRegionDataType==Region::ASIMOVDATA)
                    h_data->SetBinContent( i_bin,0 );
                else
                    h_data->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fData->fHist->Integral() );
            }
        }
    }
    //
    TthPlot *p;
    //
    // For 4-top-style plots
    if(TtHFitter::OPTION["FourTopStyle"]>0){
//         p = new TthPlot(fName+"_summary",900.*TtHFitter::OPTION["CanvasWidth"]/600.,TtHFitter::OPTION["CanvasHeight"]);
        p = new TthPlot(fInputName+"_summary",TtHFitter::OPTION["CanvasWidthSummary"],TtHFitter::OPTION["CanvasHeight"]);
        if(fYmin!=0) p->fYmin = fYmin;
        else         p->fYmin = 1;
        if(fYmax!=0) p->fYmax = fYmax;
        else         p->SetYmaxScale(3);
        p->SetXaxis("",false);
        p->fLegendNColumns = TtHFitter::OPTION["LegendNColumnsSummary"];
        if(!checkVR && fSummaryPlotLabels.size()>0){
//             if(isPostFit) p->AddLabel("#font[52]{B-only fit}");
            if(isPostFit) p->AddLabel("#font[52]{Post-fit}");
            else          p->AddLabel("#font[52]{Pre-fit}");
        }
        else if(checkVR && fSummaryPlotValidationLabels.size()>0){
//                 if(isPostFit) p->AddLabel("#font[52]{B-only fit}");
                if(isPostFit) p->AddLabel("#font[52]{Post-fit}");
                else          p->AddLabel("#font[52]{Pre-fit}");
        }
        else{
            p->AddLabel(fLabel);
//             if(isPostFit) p->AddLabel("#font[52]{B-only fit}");
            if(isPostFit) p->AddLabel("#font[52]{Post-fit}");
            else          p->AddLabel("#font[52]{Pre-fit}");
        }
    }
    // TThbb style
    else if(TtHFitter::OPTION["TtHbbStyle"]>0){
        p = new TthPlot(fInputName+"_summary",TtHFitter::OPTION["CanvasWidthSummary"],TtHFitter::OPTION["CanvasHeight"]);
        if(fYmin!=0) p->fYmin = fYmin;
        else         p->fYmin = 1;
        if(fYmax!=0) p->fYmax = fYmax;
        else         p->SetYmaxScale(3);
        p->SetXaxis("",false);
        p->SetYaxis("Events / bin");
        p->fLegendNColumns = TtHFitter::OPTION["LegendNColumnsSummary"];
        p->AddLabel(fLabel);
        if(isPostFit) p->AddLabel("Post-Fit");
        else          p->AddLabel("Pre-Fit");
    }
    //
    // normal-/old-style plots
    else{
        p = new TthPlot(fInputName+"_summary",900,700);
        if(fYmin!=0) p->fYmin = fYmin;
        else         p->fYmin = 1;
        if(fYmax!=0) p->fYmax = fYmax;
        else         p->SetYmaxScale(2);
        p->SetXaxis("",false);
        p->AddLabel(fLabel);
        if(TtHFitter::OPTION["NoPrePostFit"]==0){
            if(isPostFit) p->AddLabel("Post-Fit");
            else          p->AddLabel("Pre-Fit");
        }
    }
    //
    if(isPostFit) p->fRatioYmax = fRatioYmaxPostFit;
    else          p->fRatioYmax = fRatioYmax;
    if(isPostFit) p->fRatioYmin = fRatioYminPostFit;
    else          p->fRatioYmin = fRatioYmin;
    //
    p->fATLASlabel = fAtlasLabel;
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    if(fBlindingThreshold>=0){
        p->SetBinBlinding(true,fBlindingThreshold);
        if(isPostFit && fKeepPrefitBlindedBins && fBlindedBins) p->SetBinBlinding(true,fBlindedBins);
    }
    //
    if(h_data) p->SetData(h_data, h_data->GetTitle());
    for(int i=0;i<Nsig;i++){
        if(TtHFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    h_sig[i],h_sig[i]->GetTitle());
//         if(TtHFitter::SHOWNORMSIG)    p->AddNormSignal(h_sig[i],((string)h_sig[i]->GetTitle())+"*");
        if(TtHFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(h_sig[i],((string)h_sig[i]->GetTitle()));
        if(TtHFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(h_sig[i],(h_sig[i]->GetTitle()));
    }
    for(int i=0;i<Nbkg;i++){
        p->AddBackground(h_bkg[i],h_bkg[i]->GetTitle());
    }


    if( TtHFitter::PREFITONPOSTFIT and isPostFit) {
      p->h_tot_bkg_prefit = (TH1*)prefit_plot->GetTotBkg()->Clone("h_tot_bkg_prefit");
    }
    
    //
    // Build tot
    //
    h_tot = new TH1F("h_Tot_summary","h_Tot_summary", Nbin,0,Nbin);

    for(int i_bin=1;i_bin<=Nbin;i_bin++){
//         if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot_postFit->Integral() );
//         else          h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot->Integral() );
//         h_tot->SetBinError( i_bin,0 );
        double mc_stat_err;
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot_postFit->IntegralAndError(1,fRegions[regionVec[i_bin-1]]->fTot_postFit->GetNbinsX(),mc_stat_err) );
        else          h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot->IntegralAndError(1,fRegions[regionVec[i_bin-1]]->fTot->GetNbinsX(),mc_stat_err) );
        h_tot->SetBinError( i_bin,mc_stat_err );
    }
    //
    //   Build error band
    // build the vectors of variations
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    TH1* h_tmp_Up;
    TH1* h_tmp_Down;
    std::vector<string> systNames;
    std::vector<string> npNames;
    systNames.clear();
    npNames.clear();
    int i_np = -1;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        string systName = fSystematics[i_syst]->fName;
        string systNuisPar = systName;
        systNames.push_back( systName );
        if(fSystematics[i_syst]!=0x0)
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
                if(systNuisPar==TtHFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
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
                h_up.  push_back( new TH1F(Form("h_Tot_%s_Up_TMP"  ,systName.c_str()), Form("h_Tot_%s_Up_TMP",  systName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1F(Form("h_Tot_%s_Down_TMP",systName.c_str()), Form("h_Tot_%s_Down_TMP",systName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_np]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np]->SetBinContent( i_bin,h_tmp_Down->Integral() );
            //
            // look for other syst with the same np
            for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                if(j_syst==syst_idx) continue;
                if(systNuisPar==TtHFitter::NPMAP[ fRegions[regionVec[i_bin-1]]->fSystNames[j_syst] ]){
                    TH1* h_tmp = 0x0;
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
    // FIXME SF 
    // add the norm factors
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        string normName = fNormFactors[i_norm]->fName;
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
                h_up.  push_back( new TH1F(Form("h_Tot_%s_Up_TMP"  ,normName.c_str()), Form("h_Tot_%s_Up_TMP",  normName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1F(Form("h_Tot_%s_Down_TMP",normName.c_str()), Form("h_Tot_%s_Down_TMP",normName.c_str()), Nbin,0,Nbin) );
            }
//             h_up[fNSyst+i_norm]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
//             h_down[fNSyst+i_norm]->SetBinContent( i_bin,h_tmp_Down->Integral() );
            h_up[i_np+1+i_norm]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_np+1+i_norm]->SetBinContent( i_bin,h_tmp_Down->Integral() );
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
        line.DrawLine(divisionVec[0],((TH1F*)p->pad0->GetPrimitive("h_dummy"))->GetMinimum(),
                      divisionVec[0],pow(((TH1F*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(),0.73) );
        p->pad1->cd();
        line.DrawLine(divisionVec[0],((TH1F*)p->pad1->GetPrimitive("h_dummy2"))->GetMinimum(),
                      divisionVec[0],((TH1F*)p->pad1->GetPrimitive("h_dummy2"))->GetMaximum() );
    }
    //
    p->pad0->cd();
    if(!checkVR){
        if(fSummaryPlotLabels.size()>0){
            TLatex tex;
            tex.SetNDC(0);
            tex.SetTextAlign(20);
            //
            for(int ii=0;ii<=(int)divisionVec.size();ii++){
                if(fSummaryPlotLabels.size()<ii+1) break;
                if(divisionVec.size()<ii) break;
                float xmax = Nbin;
                float xmin = 0;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                float xpos = xmin + 0.5*(xmax - xmin);
                float ypos = pow(((TH1F*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(), 0.61 );
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
            for(int ii=0;ii<=(int)divisionVec.size();ii++){
                if(fSummaryPlotValidationLabels.size()<ii+1) break;
                if(divisionVec.size()<ii) break;
                float xmax = Nbin;
                float xmin = 0;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                float xpos = xmin + 0.5*(xmax - xmin);
                float ypos = pow(((TH1F*)p->pad0->GetPrimitive("h_dummy"))->GetMaximum(), 0.61 );
                tex.DrawLatex(xpos,ypos,fSummaryPlotValidationLabels[ii].c_str());
            }
        }
    }
    //
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        if(TtHFitter::DEBUGLEVEL>0) std::cout << i_bin << ":\t" << h_tot->GetBinContent(i_bin) << "\t+" << g_err->GetErrorYhigh(i_bin-1) << "\t-" << g_err->GetErrorYlow(i_bin-1) << std::endl;
    }
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++){
        if(isPostFit)  p->SaveAs((fName+"/Plots/Summary_postFit"+(checkVR?"_VR":"")+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
        else           p->SaveAs((fName+"/Plots/Summary"        +(checkVR?"_VR":"")+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
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
void TtHFit::BuildYieldTable(string opt,string group){
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Building Yields Table..." << std::endl;
    bool isPostFit = opt.find("post")!=string::npos;
    ofstream out;
    ofstream texout;
    gSystem->mkdir(fName.c_str(),true);
    gSystem->mkdir((fName+"/Tables").c_str());
    string suffix = "";
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
    TH1F* h_smp[MAXsamples];
    TH1F *h_tot;
    TGraphAsymmErrors *g_err[MAXsamples];
    TGraphAsymmErrors *g_err_tot;
    int Nbkg = 0;
    //
    string name;
    string title;
    float err;
    //
    double intErr; // to store the integral error
    TH1* h0; // to store varius histograms temporary
    //
    // Building region - bin correspondence
    //
    std::vector<int> regionVec; regionVec.clear();
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
//         if(!checkVR && fRegions[i_ch]->fRegionType!=Region::VALIDATION){
            if(group!="" && fRegions[i_ch]->fGroup!=group) continue;
            regionVec.push_back(i_ch);
//         } else if(checkVR && fRegions[i_ch]->fRegionType==Region::VALIDATION){
//             regionVec.push_back(i_ch);
//         }
    }
    if(regionVec.size()==0) return;
    int Nbin = regionVec.size();
    //
    out << " |       | ";
    for(int i_bin=1;i_bin<=regionVec.size();i_bin++){
        out << fRegions[regionVec[i_bin-1]]->fLabel << " | ";
    }
    out << endl;
    if(fTableOptions.find("STANDALONE")!=string::npos){
        texout << "\\documentclass[10pt]{article}" << endl;
        texout << "\\usepackage{siunitx}" << endl;
        texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
        texout << "\\begin{document}" << endl;
    }
    if(fTableOptions.find("LANDSCAPE")!=string::npos){
        texout << "\\begin{landscape}" << endl;        
    }
    texout << "\\begin{table}[htbp]" << endl;
    texout << "\\begin{center}" << endl;
    if(fTableOptions.find("FOOTNOTESIZE")!=string::npos){
        texout << "\\footnotesize" << endl;        
    }
    texout << "\\begin{tabular}{|c" ;
    for(int i_bin=1;i_bin<=regionVec.size();i_bin++){
        texout << "|c";
    }
    texout << "|}" << endl;
    texout << "\\hline " << endl;
    for(int i_bin=1;i_bin<=regionVec.size();i_bin++){
//         texout << " & " << fRegions[regionVec[i_bin-1]]->fLabel ;
        if(fRegions[regionVec[i_bin-1]]->fTexLabel!="") texout << " & " << fRegions[regionVec[i_bin-1]]->fTexLabel ;
        else                                            texout << " & " << fRegions[regionVec[i_bin-1]]->fLabel ;
    }
    texout << "\\\\" << endl;
    texout << "\\hline " << endl;
    //
    std::vector< string > titleVec;
    std::vector< int > idxVec;
    SampleHist *sh = 0x0;
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
            h_smp[idxVec[i_smp]] = new TH1F(("h_"+name).c_str(),title.c_str(), Nbin,0,Nbin);
        }
        for(int i_bin=1;i_bin<=regionVec.size();i_bin++){
            sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
            if(sh!=0x0){
                if(isPostFit && fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::GHOST)
//                 if(isPostFit && fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fUseSystematics)
                    h0 = sh->fHist_postFit;
                else
                    h0 = sh->fHist;
                float tmpErr = h_smp[idxVec[i_smp]]->GetBinError(i_bin); // Michele -> get the error before adding content to bin, to avoid ROOT automatically increasing it!
                h_smp[idxVec[i_smp]]->AddBinContent( i_bin,h0->IntegralAndError(1,h0->GetNbinsX(),intErr) );
                if(!fUseStatErr || !sh->fSample->fUseMCStat) h_smp[idxVec[i_smp]]->SetBinError(i_bin,0.);
                else                                         h_smp[idxVec[i_smp]]->SetBinError(i_bin, sqrt( pow(tmpErr,2) + pow(intErr,2) ) );
            }
        }
        titleVec.push_back(title);
    }
    //
    // add tot uncertainty on each sample
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
//         if(!fSamples[i_smp]->fUseSystematics) continue;
        name = fSamples[i_smp]->fName;
        if(idxVec[i_smp]!=i_smp) continue;
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        // build the vectors of variations
        std::vector< TH1* > h_up;   h_up.clear();
        std::vector< TH1* > h_down; h_down.clear();
        TH1* h_tmp_Up;
        TH1* h_tmp_Down;
        std::vector<string> systNames;
        systNames.clear();
        for(int i_syst=0;i_syst<fNSyst;i_syst++){
            string systName = fSystematics[i_syst]->fName;
            systNames.push_back( systName );
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
                if(sh!=0x0){
                    if(isPostFit){
                        if(syst_idx<0){
                            h_tmp_Up   = sh->fHist_postFit;
                            h_tmp_Down = sh->fHist_postFit;
                        }
                        else{
                            h_tmp_Up   = sh->GetSystematic(systName)->fHistUp_postFit;
                            h_tmp_Down = sh->GetSystematic(systName)->fHistDown_postFit;
                        }
                    }
                    else {
                        if(syst_idx<0){
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
                    h_tmp_Up   = new TH1F(Form("h_DUMMY_%s_up_%i",  systName.c_str(),i_bin-1),"h_dummy",1,0,1);
                    h_tmp_Down = new TH1F(Form("h_DUMMY_%s_down_%i",systName.c_str(),i_bin-1),"h_dummy",1,0,1);
                }
                if(i_bin==1){
                    h_up.  push_back( new TH1F(Form("h_%s_%s_Up_TMP",  name.c_str(),systName.c_str()),Form("h_%s_%s_Up_TMP",  name.c_str(),systName.c_str()), Nbin,0,Nbin) );
                    h_down.push_back( new TH1F(Form("h_%s_%s_Down_TMP",name.c_str(),systName.c_str()),Form("h_%s_%s_Down_TMP",name.c_str(),systName.c_str()), Nbin,0,Nbin) );
                }
                h_up[i_syst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral(1,h_tmp_Up  ->GetNbinsX()) );
                h_down[i_syst]->SetBinContent( i_bin,h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()) );
                //
                // eventually add any other samples with the same title
                for(int j_smp=0;j_smp<fNSamples;j_smp++){
                    sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( fSamples[j_smp]->fName );
                    if(sh==0x0) continue;
                    if(idxVec[j_smp]==i_smp && i_smp!=j_smp){
                        if(isPostFit){
                            if(syst_idx<0){
                                h_tmp_Up   = sh->fHist_postFit;
                                h_tmp_Down = sh->fHist_postFit;
                            }
                            else{
                                h_tmp_Up   = sh->GetSystematic(systName)->fHistUp_postFit;
                                h_tmp_Down = sh->GetSystematic(systName)->fHistDown_postFit;
                            }
                        }
                        else{
                            if(syst_idx<0){
                                h_tmp_Up   = sh->fHist;
                                h_tmp_Down = sh->fHist;
                            }
                            else{
                                h_tmp_Up   = sh->GetSystematic(systName)->fHistUp;
                                h_tmp_Down = sh->GetSystematic(systName)->fHistDown;
                            }
                        }
                        h_up[i_syst]  ->AddBinContent( i_bin,h_tmp_Up  ->Integral(1,h_tmp_Up->GetNbinsX()) );
                        h_down[i_syst]->AddBinContent( i_bin,h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()) );
                    }
                }
            }
        }
	// FIXME SF 
        //
        // Only for post-fit, loop on norm factors as well
        if(isPostFit){
            for(int i_norm=0;i_norm<fNNorm;i_norm++){
                string normName = fNormFactors[i_norm]->fName;
                systNames.push_back( normName );
                for(int i_bin=1;i_bin<=Nbin;i_bin++){
                    sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                    //
                    // find the normfactor in the region
                    int syst_idx = -1;
                    for(int j_syst=0;j_syst<(int)fRegions[regionVec[i_bin-1]]->fSystNames.size();j_syst++){
                        if(normName==fRegions[regionVec[i_bin-1]]->fSystNames[j_syst]){
                            syst_idx = j_syst;
                        }
                    }
                    //
                    if(sh!=0x0){
                        if(isPostFit){
                            if(syst_idx<0){
                                h_tmp_Up   = sh->fHist_postFit;
                                h_tmp_Down = sh->fHist_postFit;
                            }
                            else{
                                h_tmp_Up   = sh->GetSystematic(normName)->fHistUp_postFit;
                                h_tmp_Down = sh->GetSystematic(normName)->fHistDown_postFit;
                            }
                        }
                        else {
                            if(syst_idx<0){
                                h_tmp_Up   = sh->fHist;
                                h_tmp_Down = sh->fHist;
                            }
                            else{
                                h_tmp_Up   = sh->GetSystematic(normName)->fHistUp;
                                h_tmp_Down = sh->GetSystematic(normName)->fHistDown;
                            }
                        }
                    }
                    else {
                        h_tmp_Up   = new TH1F(Form("h_DUMMY_%s_up_%i",  normName.c_str(),i_bin-1),"h_dummy",1,0,1);
                        h_tmp_Down = new TH1F(Form("h_DUMMY_%s_down_%i",normName.c_str(),i_bin-1),"h_dummy",1,0,1);
                    }
                    if(i_bin==1){
                        h_up.  push_back( new TH1F(Form("h_%s_%s_Up_TMP",  name.c_str(),normName.c_str()),Form("h_%s_%s_Up_TMP",  name.c_str(),normName.c_str()), Nbin,0,Nbin) );
                        h_down.push_back( new TH1F(Form("h_%s_%s_Down_TMP",name.c_str(),normName.c_str()),Form("h_%s_%s_Down_TMP",name.c_str(),normName.c_str()), Nbin,0,Nbin) );
                    }
                    h_up[i_norm+fNSyst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral(1,h_tmp_Up  ->GetNbinsX()) );
                    h_down[i_norm+fNSyst]->SetBinContent( i_bin,h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()) );
                    //
                    // eventually add any other samples with the same title
                    for(int j_smp=0;j_smp<fNSamples;j_smp++){
                        sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( fSamples[j_smp]->fName );
                        if(sh!=0){
                            if(idxVec[j_smp]==i_smp && i_smp!=j_smp){
                                if(isPostFit){
                                    if(syst_idx<0){
                                        h_tmp_Up   = sh->fHist_postFit;
                                        h_tmp_Down = sh->fHist_postFit;
                                    }
                                    else{
                                        h_tmp_Up   = sh->GetSystematic(normName)->fHistUp_postFit;
                                        h_tmp_Down = sh->GetSystematic(normName)->fHistDown_postFit;
                                    }
                                }
                                else{
                                    if(syst_idx<0){
                                        h_tmp_Up   = sh->fHist;
                                        h_tmp_Down = sh->fHist;
                                    }
                                    else{
                                        h_tmp_Up   = sh->GetSystematic(normName)->fHistUp;
                                        h_tmp_Down = sh->GetSystematic(normName)->fHistDown;
                                    }
                                }
                                h_up[i_norm+fNSyst]  ->AddBinContent( i_bin,h_tmp_Up  ->Integral(1,h_tmp_Up->GetNbinsX()) );
                                h_down[i_norm+fNSyst]->AddBinContent( i_bin,h_tmp_Down->Integral(1,h_tmp_Down->GetNbinsX()) );
                            }
                        }
                    }
                }
            }
        }
        //
        if(isPostFit)  g_err[i_smp] = BuildTotError( h_smp[i_smp], h_up, h_down, systNames, fFitResults->fCorrMatrix );
        else           g_err[i_smp] = BuildTotError( h_smp[i_smp], h_up, h_down, systNames );
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
//         texout << "  " << fSamples[i_smp]->fTitle << "  ";
        if(fSamples[i_smp]->fType==Sample::DATA) texout << "\\hline " << endl;
        if(fSamples[i_smp]->fTexTitle!="") texout << "  " << fSamples[i_smp]->fTexTitle << "  ";
        else                               texout << "  " << fSamples[i_smp]->fTitle << "  ";
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
            texout << " & ";
            out << h_smp[i_smp]->GetBinContent(i_bin);
//             texout << setprecision(3) << h_smp[i_smp]->GetBinContent(i_bin);
            texout << "\\num[round-mode=figures,round-precision=3]{";
            texout << h_smp[i_smp]->GetBinContent(i_bin);
            texout << "}";
            out << " pm ";
            texout << " $\\pm$ ";
            out << ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
//             texout << setprecision(3) << ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
            texout << "\\num[round-mode=figures,round-precision=3]{";
            texout << ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
            texout << "}";
            out << " | ";
        }
        out << endl;
        texout << " \\\\ ";
        texout << endl;
    }
    //
    // Build tot
    //
    h_tot = new TH1F("h_Tot_","h_Tot", Nbin,0,Nbin);
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
    std::vector<string> systNames;
    systNames.clear();
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        string systName = fSystematics[i_syst]->fName;
        systNames.push_back( systName );
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
                h_up.  push_back( new TH1F(Form("h_Tot_%s_Up_TMP"  ,systName.c_str()), Form("h_Tot_%s_Up_TMP",  systName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1F(Form("h_Tot_%s_Down_TMP",systName.c_str()), Form("h_Tot_%s_Down_TMP",systName.c_str()), Nbin,0,Nbin) );
            }
            h_up[i_syst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_syst]->SetBinContent( i_bin,h_tmp_Down->Integral() );
        }
    }
    // FIXME SF 
    // add the norm factors
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        string normName = fNormFactors[i_norm]->fName;
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
                h_up.  push_back( new TH1F(Form("h_Tot_%s_Up_TMP"  ,normName.c_str()), Form("h_Tot_%s_Up_TMP",  normName.c_str()), Nbin,0,Nbin) );
                h_down.push_back( new TH1F(Form("h_Tot_%s_Down_TMP",normName.c_str()), Form("h_Tot_%s_Down_TMP",normName.c_str()), Nbin,0,Nbin) );
            }
            h_up[fNSyst+i_norm]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[fNSyst+i_norm]->SetBinContent( i_bin,h_tmp_Down->Integral() );
        }
    }
    //
    if(isPostFit)  g_err_tot = BuildTotError( h_tot, h_up, h_down, systNames, fFitResults->fCorrMatrix );
    else           g_err_tot = BuildTotError( h_tot, h_up, h_down, systNames );
    //
    out << " | Total | ";
    texout << "\\hline " << endl;
    texout << "  Total ";
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        texout << " & ";
        out << h_tot->GetBinContent(i_bin);
//         texout << setprecision(3) << h_tot->GetBinContent(i_bin);
        texout << "\\num[round-mode=figures,round-precision=3]{";
        texout << h_tot->GetBinContent(i_bin);
        texout << "}";
        out << " pm ";
        texout << " $\\pm$ ";
        out << ( g_err_tot->GetErrorYhigh(i_bin-1) + g_err_tot->GetErrorYlow(i_bin-1) )/2.;
//         texout << setprecision(3) << g_err_tot->GetErrorYhigh(i_bin-1);
        texout << "\\num[round-mode=figures,round-precision=3]{";
        texout << ( g_err_tot->GetErrorYhigh(i_bin-1) + g_err_tot->GetErrorYlow(i_bin-1) )/2.;
        texout << "}";
        out << " | ";
    }
    out << endl;
    texout << " \\\\ ";
    texout << endl;

    //
    // Print data
    if( !fFitIsBlind ){
        texout << "\\hline " << endl;
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
                texout << h_smp[i_smp]->GetBinContent(i_bin);
                out << " | ";
            }
            out << endl;
            texout << " \\\\ ";
            texout << endl;
        }
    }
    //
    texout << "\\hline " << endl;
    texout << "\\end{tabular} " << endl;
    texout << "\\caption{Yields of the analysis} " << endl;
    texout << "\\end{center} " << endl;
    texout << "\\end{table} " << endl;
    if(fTableOptions.find("LANDSCAPE")!=string::npos){
        texout << "\\end{landscape}" << endl;        
    }
    if(fTableOptions.find("STANDALONE")!=string::npos){
        texout << "\\end{document}" << endl;
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
void TtHFit::DrawSignalRegionsPlot(int nCols,int nRows){
    std::vector< Region* > vRegions;
    vRegions.clear();
    if(fRegionsToPlot.size()>0){
        nCols = 1;
        nRows = 1;
        // first loop
        int nRegInRow = 0;
        for(unsigned int i=0;i<fRegionsToPlot.size();i++){
            if(TtHFitter::DEBUGLEVEL>0) std::cout << fRegionsToPlot[i] << std::endl;
            if(fRegionsToPlot[i].find("ENDL")!=string::npos){
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
void TtHFit::DrawSignalRegionsPlot(int nCols,int nRows, std::vector < Region* > &regions){
    gSystem->mkdir(fName.c_str(), true);
    float Hp = 250; // height of one mini-plot, in pixels
    float Wp = 200; // width of one mini-plot, in pixels
    float H0 = 100; // height of the top label pad
    if(TtHFitter::OPTION["FourTopStyle"]!=0) H0 = 75; // height of the top label pad
    if(TtHFitter::OPTION["FourTopStyle"]!=0) Hp = 200;
    if(TtHFitter::OPTION["SignalRegionSize"]!=0){
        Hp = TtHFitter::OPTION["SignalRegionSize"];
        Wp = (200./250.)*TtHFitter::OPTION["SignalRegionSize"];
    }
    float H = H0 + nRows*Hp; // tot height of the canvas
    float W = nCols*Wp; // tot width of the canvas
    if(TtHFitter::OPTION["FourTopStyle"]!=0) W += 50.; // FIXME
    else W += 0.1; // to fix eps format (why is this needed?)

    TCanvas *c = new TCanvas("c","c",W,H);
    TPad *pTop = new TPad("c0","c0",0,1-H0/H,1,1);
    pTop->Draw();
    pTop->cd();
    ATLASLabel(0.1/(W/200.),1.-0.3*(100./H0),fAtlasLabel.c_str());
    myText(    0.1/(W/200.),1.-0.6*(100./H0),1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(fLabel!="-") myText(    0.1/(W/200.),1.-0.9*(100./H0),1,Form("%s",fLabel.c_str()));

    TLegend *leg = 0x0;
    if(TtHFitter::OPTION["FourTopStyle"]!=0){
        leg = new TLegend(0.35,1.-0.6*(100./H0),1,1.-0.3*(100./H0));
        leg->SetNColumns(3);
        leg->SetFillStyle(0.);
        leg->SetBorderSize(0.);
        leg->SetTextSize( gStyle->GetTextSize() );
        leg->SetTextFont( gStyle->GetTextFont() );
        leg->Draw();
    }

    c->cd();

    TPad *pLeft = new TPad("c1","c1",0,0,0+(W-nCols*Wp)/W,1-H0/H);
    pLeft->Draw();
    pLeft->cd();
    TLatex *tex0 = new TLatex();
    tex0->SetNDC();
    tex0->SetTextAngle(90);
    tex0->SetTextAlign(23);
    tex0->DrawLatex(0.4,0.5,"S / #sqrt{ B }");

    c->cd();

//     TPad *pBottom = new TPad("c1","c1",0,0,1,1-H0/H);
    TPad *pBottom = new TPad("c1","c1",0+(W-nCols*Wp)/W,0,1,1-H0/H);
    pBottom->Draw();
    pBottom->cd();

    pBottom->Divide(nCols,nRows);
    int Nreg = nRows*nCols;
    if(Nreg>(int)regions.size()) Nreg = regions.size();
    TH1F* h[Nreg];
    float S[Nreg];
    float B[Nreg];
    double xbins[] = {0,0.1,0.9,1.0};
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(gStyle->GetTextSize());
    pBottom->cd(1);

    //
    // Get the values
    //
    for(int i=0;i<Nreg;i++){
        S[i] = 0.;
        B[i] = 0.;
        if(regions[i]==0x0) continue;
        for(int i_sig=0;i_sig<regions[i]->fNSig;i_sig++)
            if(regions[i]->fSig[i_sig]!=0x0)
                S[i] += regions[i]->fSig[i_sig]->fHist->Integral();
        for(int i_bkg=0;i_bkg<regions[i]->fNBkg;i_bkg++){
            if(regions[i]->fBkg[i_bkg]!=0x0)
                B[i] += regions[i]->fBkg[i_bkg]->fHist->Integral();
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
    for(int i=0;i<Nreg;i++){
        if(regions[i]==0x0) continue;
        pBottom->cd(i+1);
        if(TtHFitter::OPTION["LogSignalRegionPlot"]) gPad->SetLogy();
        if(TtHFitter::OPTION["FourTopStyle"]!=0){
            gPad->SetLeftMargin(0.1);
            gPad->SetRightMargin(0.);
        }
        string label = regions[i]->fShortLabel;
        h[i] = new TH1F(Form("h[%d]",i),label.c_str(),3,xbins);
        h[i]->SetBinContent(2,S[i]/sqrt(B[i]));
        if(TtHFitter::OPTION["FourTopStyle"]==0) h[i]->GetYaxis()->SetTitle("S / #sqrt{B}");
        h[i]->GetYaxis()->CenterTitle();
        //     h[i]->GetYaxis()->SetTitleSize(0.14);
//         h[i]->GetYaxis()->SetLabelOffset(1.5*h[i]->GetYaxis()->GetLabelOffset());
        h[i]->GetYaxis()->SetLabelOffset(1.5*h[i]->GetYaxis()->GetLabelOffset() / (Wp/200.));
        h[i]->GetYaxis()->SetTitleOffset(9*nRows/4. );
        if(Wp<200) h[i]->GetYaxis()->SetTitleOffset( h[i]->GetYaxis()->GetTitleOffset()*0.90 );
        h[i]->GetYaxis()->SetLabelSize( h[i]->GetYaxis()->GetLabelSize() * (Wp/200.) );
        if(TtHFitter::OPTION["FourTopStyle"]!=0) h[i]->GetYaxis()->SetLabelSize( h[i]->GetYaxis()->GetLabelSize() * 1.1 );
        h[i]->GetXaxis()->SetTickLength(0);
        if(TtHFitter::OPTION["LogSignalRegionPlot"]==0) h[i]->GetYaxis()->SetNdivisions(3);
//         if(TtHFitter::OPTION["LogSignalRegionPlot"]==0) //h[i]->GetYaxis()->SetMaxDigits(3);
        else //h[i]->GetYaxis()->SetNdivisions(2);
          TGaxis::SetMaxDigits(5);
        yMax = TMath::Max(yMax,h[i]->GetMaximum());
        h[i]->GetXaxis()->SetLabelSize(0);
        h[i]->SetLineWidth(1);
        h[i]->SetLineColor(kBlack);
        if(regions[i]->fRegionType==Region::SIGNAL)          h[i]->SetFillColor(kRed+1);
        else if(regions[i]->fRegionType==Region::VALIDATION) h[i]->SetFillColor(kGray);
        else                                                 h[i]->SetFillColor(kAzure-4);
        if(leg!=0x0){
            if(regions[i]->fRegionType==Region::CONTROL && !hasCR)    { leg->AddEntry(h[i],"Control Regions","f");    hasCR = true; }
            if(regions[i]->fRegionType==Region::VALIDATION && !hasVR) { leg->AddEntry(h[i],"Validation Regions","f"); hasVR = true; }
            if(regions[i]->fRegionType==Region::SIGNAL && !hasSR)     { leg->AddEntry(h[i],"Signal Regions","f");     hasSR = true; }
        }
        h[i]->Draw();
        gPad->SetLeftMargin( gPad->GetLeftMargin()*2.4 );
        gPad->SetRightMargin(gPad->GetRightMargin()*0.1);
        gPad->SetTicky(0);
        gPad->RedrawAxis();
        if(TtHFitter::OPTION["FourTopStyle"]==0) tex->DrawLatex(0.42,0.85,label.c_str());
        else                                     tex->DrawLatex(0.27,0.85,label.c_str());
        float SoB = S[i]/B[i];
        string SB = Form("%.1f%%",(100.*SoB));
//         if(Wp<200) SB = Form("#scale[0.75]{S/B=}%.1f%%",(100.*SoB));
        if(TtHFitter::OPTION["FourTopStyle"]!=0){
            if( (100.*SoB)<0.1 ){
//                 SB = Form("S/B = %.0e%%",(100.*SoB));
                SB = Form("%.0e%%",SoB);
                if(SB.find("0")!=string::npos)SB.replace(SB.find("0"), 1, "");
                if(SB.find("e")!=string::npos) SB.replace(SB.find("e"), 1, "#scale[0.75]{#times}10^{");
                if(SB.find("%")!=string::npos)SB.replace(SB.find("%"), 1, "}");
//                 SB = Form("S/B = %.0e%%",(100.*SoB));
            }
        }
        SB = "#scale[0.75]{S/B} = "+SB;
        if(TtHFitter::OPTION["FourTopStyle"]==0) tex->DrawLatex(0.42,0.72,SB.c_str());
        else                                     tex->DrawLatex(0.27,0.72,SB.c_str());
    }
    //
    for(int i=0;i<Nreg;i++){
        if(regions[i]==0x0) continue;
        if(TtHFitter::OPTION["LogSignalRegionPlot"]!=0){
            h[i]->SetMaximum(yMax*200);
            h[i]->SetMinimum(2e-4);
        }
        else{
            h[i]->SetMaximum(yMax*1.5);
            h[i]->SetMinimum(0.);
        }
    }
    //
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs((fName+"/SignalRegions"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str());

    //
    delete c;
}

//__________________________________________________________________________________
//
void TtHFit::DrawPieChartPlot(const std::string &opt, int nCols,int nRows){

    std::vector< Region* > vRegions;
    vRegions.clear();
    if(fRegionsToPlot.size()>0){
        nCols = 1;
        nRows = 1;
        // first loop
        int nRegInRow = 0;
        for(unsigned int i=0;i<fRegionsToPlot.size();i++){
            if(TtHFitter::DEBUGLEVEL>0) std::cout << fRegionsToPlot[i] << std::endl;
            if(fRegionsToPlot[i].find("ENDL")!=string::npos){
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
void TtHFit::DrawPieChartPlot(const std::string &opt, int nCols,int nRows, std::vector < Region* > &regions ){

    float Hp = 250; // height of one mini-plot, in pixels
    float Wp = 250; // width of one mini-plot, in pixels
    float H0 = 100; // height of the top label pad
    if(TtHFitter::OPTION["FourTopStyle"]>0) H0 = 75; // height of the top label pad

    if(TtHFitter::OPTION["PieChartSize"]!=0){
        Hp = TtHFitter::OPTION["PieChartSize"];
        Wp = TtHFitter::OPTION["PieChartSize"];
    }

    float H = H0 + nRows*Hp; // tot height of the canvas
    float W = nCols*Wp; // tot width of the canvas

    bool isPostFit = opt.find("post")!=string::npos;

    //
    // Create the canvas
    //
    TCanvas *c = new TCanvas("c","c",W,H);
    TPad *pTop = new TPad("c0","c0",0,1-H0/H,1,1);
    pTop->Draw();
    pTop->cd();

    if(TtHFitter::OPTION["FourTopStyle"]>0){
        ATLASLabel(0.1/(W/200.),1.-0.3*(100./H0),fAtlasLabel.c_str());
        myText(    0.1/(W/200.),1.-0.6*(100./H0),1,Form("#sqrt{s} = %s",fCmeLabel.c_str()));
        if(fLabel!="-") myText(    0.1/(W/200.),1.-0.9*(100./H0),1,Form("%s",fLabel.c_str()));
    }
    else{
        ATLASLabel(0.05 / (W/200),0.7,fAtlasLabel.c_str());
//         myText(    0.05 / (W/200),0.4,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
        myText(    0.05 / (W/200),0.4,1,Form("#sqrt{s} = %s",fCmeLabel.c_str()));
        if(fLabel!="-") myText(    0.05 / (W/200),0.1,1,Form("%s",fLabel.c_str()));
    }

    c->cd();
    TPad *pBottom = new TPad("c1","c1",0,0,1,1-H0/H);
    pBottom->Draw();
    pBottom->cd();
    pBottom->Divide(nCols,nRows);
    int Nreg = nRows*nCols;
    if(Nreg>(int)regions.size()) Nreg = regions.size();
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(gStyle->GetTextSize());
    pBottom->cd(1);

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

        if(regions[i]!=0x0){
            for(int i_bkg=regions[i]->fNBkg-1;i_bkg>=0;i_bkg--){
                if(regions[i]->fBkg[i_bkg]!=0x0){
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
    for(int i=0;i<Nreg;i++){
        if(regions[i]==0x0) continue;
        pBottom->cd(i+1);
        string label = regions[i]->fShortLabel;

        const int back_n = results[i].size();
        float values[back_n];
        int colors[back_n];
        for( unsigned int iTemp = 0; iTemp < back_n; ++iTemp ){
            values[iTemp] = 0.;
            colors[iTemp] = 0;
        }

        int count = 0;
        for ( std::pair < string, double > temp_pair : results[i] ){
            values[count] = temp_pair.second;
            colors[count] = results_color[i][temp_pair.first];
            count++;
        }

        TPie *pie = new TPie(("pie_"+label).c_str()," ",back_n, values, colors);
        pie -> SetRadius( pie -> GetRadius() * 0.8 );
        for( unsigned int iEntry = 0; iEntry < pie->GetEntries(); ++iEntry) pie -> SetEntryLabel(iEntry,"");
        pie -> Draw();
        tex->DrawLatex(0.1,0.85,label.c_str());
    }

    c -> cd();

    //
    // Adding the legend in the top panel
    //
    pTop->cd();
    TLegend *leg;
    if(TtHFitter::OPTION["FourTopStyle"]>0 || TtHFitter::OPTION["TtHbbStyle"]>0){
        leg = new TLegend(0.5,0.1,0.95,0.90);
        leg -> SetNColumns(3);
    }
    else{
        leg = new TLegend(0.7,0.1,0.95,0.90);
        if(map_for_legend.size()>4){
            leg -> SetNColumns(2);
        }
    }

    leg -> SetLineStyle(0);
    leg -> SetFillStyle(0);
    leg -> SetLineColor(0);
    leg -> SetBorderSize(0);
    leg -> SetTextFont( gStyle->GetTextFont() );
    leg -> SetTextSize( gStyle->GetTextSize() );

    std::vector<std::string> legVec; legVec.clear();
    for ( const std::pair < std::string, int > legend_entry : map_for_legend ) {
//         TH1F *dummy = new TH1F( ("legend_entry_" + legend_entry.first).c_str(), "",1,0,1);
//         dummy -> SetFillColor(legend_entry.second);
//         dummy -> SetLineColor(kBlack);
//         dummy -> SetLineWidth(1);
//         leg -> AddEntry(dummy,legend_entry.first.c_str(),"f");
        legVec.push_back(legend_entry.first);
    }
    for(int i_leg=legVec.size()-1;i_leg>=0;i_leg--){
        TH1F *dummy = new TH1F( ("legend_entry_" + legVec[i_leg]).c_str(), "",1,0,1);
        dummy -> SetFillColor(map_for_legend[legVec[i_leg]]);
        dummy -> SetLineColor(kBlack);
        dummy -> SetLineWidth(1);
        leg -> AddEntry(dummy,legVec[i_leg].c_str(),"f");
    }
    leg -> Draw();

    //
    // Stores the pie chart in the desired format
    //
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++){
        c->SaveAs((fName+"/PieChart" + fSuffix + ( isPostFit ? "_postFit" : "" ) + "."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
    }

    //
    delete c;
}

//__________________________________________________________________________________
// called before w in case of CustomAsimov
void TtHFit::CreateCustomAsimov(){
    cout << "Running CreateCustomAsimov" << endl;
    // get a list of all CustomAsimov to create
    std::vector<std::string> customAsimovList;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fAsimovReplacementFor.first!="" && FindInStringVector(customAsimovList,fSamples[i_smp]->fAsimovReplacementFor.first)<0)
            customAsimovList.push_back(fSamples[i_smp]->fAsimovReplacementFor.first);
    }
    //
    // fill a different CustomAsimov data-set for each element in the list
    for(auto customAsimov : customAsimovList){
        cout << customAsimov << endl;
        Sample *ca = GetSample("customAsimov_"+customAsimov);
        // create a new data sample taking the nominal S and B
        for(int i_ch=0;i_ch<fNRegions;i_ch++){
            Region *reg = fRegions[i_ch];
            SampleHist *cash = reg->SetSampleHist(ca,(TH1*)reg->fData->fHist->Clone());
            cash->fHist->Scale(0.);
            //
            std::vector<std::string> smpToExclude; smpToExclude.clear();
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
                float factor = 1.;
                for(auto norm : fSamples[i_smp]->fNormFactors){
                    cout << "setting norm factor to " << norm->fNominal << endl;
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
void TtHFit::ToRooStat(bool makeWorkspace, bool exportOnly){

    //Suffix used for the regular bin transformed histogram
    const std::string suffix_regularBinning = "_regBin";

    if(TtHFitter::DEBUGLEVEL>0){
        std::cout << "--------------------------------" << std::endl;
        std::cout << "|      Export to RooStat       |" << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }
    else{
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "Exporting to RooStats..." << std::endl;
    }

    RooStats::HistFactory::Measurement meas((fInputName+fSuffix).c_str(), (fInputName+fSuffix).c_str());
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

        if(TtHFitter::DEBUGLEVEL>0){
            std::cout << "Adding Channel: " << fRegions[i_ch]->fName << std::endl;
        }
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
            if(cash==0x0){
                std::cerr << "TtHFit::ToRooStat::ERROR: no Custom Asimov " << fCustomAsimov << " available. Taking regular Asimov." << std::endl;
            }
            else{
                if(TtHFitter::DEBUGLEVEL>0){
                    std::cout << "  Adding Custom-Asimov Data: " << cash->fHist->GetName() << std::endl;
                }
                chan.SetData(cash->fHistoName+suffix_regularBinning, cash->fFileName);
            }
        }
        else if(hasData){
            if(TtHFitter::DEBUGLEVEL>0){
                std::cout << "  Adding Data: " << fRegions[i_ch]->fData->fHist->GetName() << std::endl;
            }
            chan.SetData(fRegions[i_ch]->fData->fHistoName+suffix_regularBinning, fRegions[i_ch]->fData->fFileName);
        } else {
            chan.SetData("", "");
        }

        chan.SetStatErrorConfig(fStatErrThres,fStatErrCons.c_str()); // "Gaussian"
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            SampleHist* h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if( h != 0x0 && h->fSample->fType!=Sample::DATA && h->fSample->fType!=Sample::GHOST ){
                if(TtHFitter::DEBUGLEVEL>0){
                    std::cout << "  Adding Sample: " << fSamples[i_smp]->fName << std::endl;
                }
                RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());
                if(fUseStatErr && fSamples[i_smp]->fUseMCStat) sample.ActivateStatError();
                sample.SetHistoName(h->fHistoName+suffix_regularBinning);
                sample.SetInputFile(h->fFileName);
                sample.SetNormalizeByTheory(fSamples[i_smp]->fNormalizedByTheory);
                // norm factors
                for(int i_norm=0;i_norm<h->fNNorm;i_norm++){
                    if(TtHFitter::DEBUGLEVEL>0){
                        std::cout << "    Adding NormFactor: " << h->fNormFactors[i_norm]->fName << ", " << h->fNormFactors[i_norm]->fNominal << std::endl;
                    }
                    sample.AddNormFactor( h->fNormFactors[i_norm]->fName,
                                         h->fNormFactors[i_norm]->fNominal,
                                         h->fNormFactors[i_norm]->fMin,
                                         h->fNormFactors[i_norm]->fMax);
                    if (h->fNormFactors[i_norm]->fConst) meas.AddConstantParam( h->fNormFactors[i_norm]->fName );
                    if (fStatOnly && fFixNPforStatOnlyFit && h->fNormFactors[i_norm]->fName!=fPOI)
                        meas.AddConstantParam( h->fNormFactors[i_norm]->fName );
                    // radom POI tricK: if this norm factor is te POI, add another (fixed) norm factor equal to a random number generated from the given seed
                    if (fRandomPOISeed>=0){
                        if ( h->fNormFactors[i_norm]->fName == fPOI ){
                            if(TtHFitter::DEBUGLEVEL>0){
                                std::cout << "    Scaling POI by random secret number XXX ;)" << std::endl;
                            }
                            gRandom->SetSeed(fRandomPOISeed); float rnd = 1./pow(gRandom->Uniform(0,2),2);
                            sample.AddNormFactor( "POIScale",rnd,rnd-1,rnd+1);
                            meas.AddConstantParam( "POIScale" );
                        }
                    }
                }
                // shape factors
                for(int i_shape=0;i_shape<h->fNShape;i_shape++){
                    if(TtHFitter::DEBUGLEVEL>0){
                        std::cout << "    Adding ShapeFactor: " << h->fShapeFactors[i_shape]->fName << ", " << h->fShapeFactors[i_shape]->fNominal << std::endl;
                    }
                    sample.AddShapeFactor( h->fShapeFactors[i_shape]->fName );
                    if (h->fShapeFactors[i_shape]->fConst) meas.AddConstantParam( h->fShapeFactors[i_shape]->fName );
                    if (fStatOnly && fFixNPforStatOnlyFit && h->fShapeFactors[i_shape]->fName!=fPOI)
                        meas.AddConstantParam( h->fShapeFactors[i_shape]->fName );
                }
                // systematics
                if(!fStatOnly){
                    for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
                        // add normalization part
                        if(TtHFitter::DEBUGLEVEL>0){
                            std::cout << "    Adding Systematic: " << h->fSyst[i_syst]->fName << std::endl;
                        }
                        if ( h->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
                            std::string npName = "shape_";
                            npName += h->fSyst[i_syst]->fSystematic->fNuisanceParameter+"_";
                            std::string regionName = fRegions[i_ch]->fName;
                            if(h->fSyst[i_syst]->fSystematic->fNuisanceParameter.find("stat_")!=string::npos){
                                // see if there are regions to correlate with others
                                for(auto set : h->fSample->fCorrelateGammasInRegions){
                                    for(int i_reg=0;i_reg<set.size();i_reg++){
                                        if(i_reg!=0 && regionName==set[i_reg]){
                                            regionName = set[0];
                                            break;
                                        }
                                    }
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
                                if(h->fSyst[i_syst]->fSystematic->fNuisanceParameter=="ttXsec") cout << h->fSyst[i_syst]->fNormUp << endl;
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
//         if(fSystematics[i_syst]->fIsFreeParameter) meas.AddNoSyst(fSystematics[i_syst]->fName.c_str());
        if(fSystematics[i_syst]->fIsFreeParameter) meas.AddUniformSyst(fSystematics[i_syst]->fName.c_str());
    }
    //
    meas.PrintXML((fName+"/RooStats/").c_str());
    meas.CollectHistograms();
    meas.PrintTree();
    if(makeWorkspace) RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
}

//__________________________________________________________________________________
//
void TtHFit::DrawPruningPlot(){
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Drawing Pruning Plot ..." << std::endl;
    if(fSystematics.size()==0 || fStatOnly){
        std::cout << "TtHFit::INFO: Stat only fit => No Pruning plot generated." << std::endl;
        return;
    }
    //
    ofstream out;
    out.open((fName+"/PruningText.txt").c_str());
    out << "-------///////                 ///////-------" << endl ;
    out << "-------/////// IN PRUNING PLOT ///////-------" << endl ;
    out << "-------///////                 ///////-------" << endl ;
    //
    vector< TH2F* > histPrun;
    vector< TH2F* > histPrun_toSave;
    int iReg = 0;
    int nSmp = 0;
    vector< Sample* > samplesVec;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        samplesVec.push_back(fSamples[i_smp]);
        nSmp++;
    }
    //
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        if(fRegions[i_reg]->fRegionType!=Region::VALIDATION){
            out << "In Region : " << fRegions[i_reg]->fName << endl ;
            histPrun.push_back( new TH2F(Form("h_prun_%s", fRegions[i_reg]->fName.c_str()  ),fRegions[i_reg]->fShortLabel.c_str(),nSmp,0,nSmp, fNSyst,0,fNSyst) );
            histPrun[histPrun.size()-1]->SetDirectory(0);
            for(int i_smp=0;i_smp<nSmp;i_smp++){
                out << " -> In Sample : " << samplesVec[i_smp]->fName << endl ;
                for(int i_syst=0;i_syst<fNSyst;i_syst++){
                   histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), -1 );
                }
                SampleHist *sh = fRegions[i_reg]->GetSampleHist(samplesVec[i_smp]->fName);
                if(sh!=0x0){
                    for(int i_syst=0;i_syst<fNSyst;i_syst++){
                        out << " --->>  " << fSystematics[i_syst]->fName << "     " ;
                        if( sh->HasSyst(fSystematics[i_syst]->fName) ) {
                            SystematicHist *syh = sh->GetSystematic(fSystematics[i_syst]->fName);
                            histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 0 );
                            //
                            // set to 1 if shape pruned away
                            if( !sh->GetSystematic(fSystematics[i_syst]->fName)->fIsShape ||
                                (fThresholdSystPruning_Shape>-1 &&
                              ( !HistoTools::HasShape(sh->fHist, sh->GetSystematic(fSystematics[i_syst]->fName),fThresholdSystPruning_Shape)
                                || FindInStringVector( fSystematics[i_syst]->fDropShapeIn, fRegions[i_reg]->fName )>=0
                                || FindInStringVector( fSystematics[i_syst]->fDropShapeIn, samplesVec[i_smp]->fName )>=0
                                ) )
                              ) {
                                syh->fShapePruned = true;
                                histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 1 );
                            }
                            float normUp=TMath::Abs(sh->GetSystematic(fSystematics[i_syst]->fName)->fNormUp);
                            float normDo=TMath::Abs(sh->GetSystematic(fSystematics[i_syst]->fName)->fNormDown);
                            // set to 2 is normalization pruned away
                            if(    FindInStringVector( fSystematics[i_syst]->fDropNormIn, fRegions[i_reg]->fName )>=0
                                || FindInStringVector( fSystematics[i_syst]->fDropNormIn, samplesVec[i_smp]->fName )>=0
                                || FindInStringVector( fSystematics[i_syst]->fDropNormIn, "all" )>=0
                                || ( fThresholdSystPruning_Normalisation>-1 && normUp<fThresholdSystPruning_Normalisation && normDo<fThresholdSystPruning_Normalisation )
                                || (normUp!=normUp || normDo!=normDo)
                                ) {
                                syh->fNormPruned = true;
                                if(syh->fShapePruned || fSystematics[i_syst]->fIsNormOnly) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 3 );
                                else                  histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 2 );
                            }
                            //
                            // now check for crazy sys ....
                            if ( fThresholdSystLarge > -1 ) {
                                // first norm:
                                if ( normUp>fThresholdSystLarge || normDo>fThresholdSystLarge || normUp!=normUp || normDo!=normDo ) {
                                    syh->fBadNorm = true;
                                    histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst),-2);
                                }
                                //
                                // then shape
                                if ( sh->GetSystematic(fSystematics[i_syst]->fName)->fIsShape && ( HistoTools::HasShape(sh->fHist, sh->GetSystematic(fSystematics[i_syst]->fName),fThresholdSystLarge) ) ) {
                                    syh->fBadShape = true;
                                    if ( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )==-2 ) {
                                        histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst),-4);
                                    }
                                    else {
                                        histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst),-3);
                                    }
                                }
                            }
                        }
                        if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== -1 ) out << " is not present" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== 0 ) out << " is kept" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== 1 ) out << " is norm only" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== 2 ) out << " is shape only" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== 3 ) out << " is dropped" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== -2 ) out << " has bad norm" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== -3 ) out << " has bad shape" << endl;
                        else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst) )== -4 ) out << " is bad" << endl;
                    }
                }
            }
            //
            histPrun_toSave.push_back( (TH2F*)histPrun[iReg]->Clone(Form("%s_toSave",histPrun[iReg]->GetName())) );
            histPrun_toSave[iReg]->SetDirectory(0);
            //
            iReg++;
        }
    }
    //
    // draw the histograms
    int upSize = 50;
    int loSize = 150;
    int mainHeight = fNSyst*20;
    int leftSize = 250;
    int regionSize = 20*fNSamples;
    int separation = 10;
    int mainWidth = iReg*(regionSize+separation);
    //
    TCanvas *c = new TCanvas("c_pruning","Canvas - Pruning",leftSize+mainWidth,upSize+mainHeight+loSize);
    Int_t colors[] = {kBlack,6,kBlue, kGray, 8, kYellow, kOrange-3, kRed}; // #colors >= #levels - 1
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    TPad *pUp = new TPad("pUp","Pad High",0,(1.*loSize+mainHeight)/(upSize+mainHeight+loSize),1,1);
    pUp->Draw();
    c->cd();
//     TPad *pLo = new TPad("pLo","Pad Low",0,0,1,(1.*loSize)/(upSize+mainHeight+loSize));
    TPad *pReg[100];
    for(int i_reg=0;i_reg<(int)histPrun.size();i_reg++){
        c->cd();
//         pReg[i_reg] = new TPad(Form("pReg[%d]",i_reg),"Pad Region",
//                               (leftSize+1.*i_reg*(regionSize+separation))           /(leftSize+mainWidth),   (1.*loSize)           /(upSize+mainHeight+loSize),
//                               (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize) );
        if(i_reg==0){
            pReg[i_reg] = new TPad(Form("pReg[%d]",i_reg),"Pad Region",
                                  0,   0,
                                  (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize) );
            pReg[i_reg]->SetLeftMargin( (1.*leftSize) / (1.*leftSize+regionSize) );
        }
        else{
            pReg[i_reg] = new TPad(Form("pReg[%d]",i_reg),"Pad Region",
                                  (leftSize+1.*i_reg*(regionSize+separation))           /(leftSize+mainWidth),   0,
                                  (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize) );
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
                histPrun[i_reg]       ->GetYaxis()->SetBinLabel(i_bin,TtHFitter::SYSTMAP[fSystematics[i_bin-1]->fName].c_str());
            }
            else {
                histPrun[i_reg]->GetYaxis()->SetBinLabel(i_bin,"");
            }
            histPrun_toSave[i_reg]->GetYaxis()->SetBinLabel(i_bin,fSystematics[i_bin-1]->fName.c_str());
        }
        histPrun[i_reg]->Draw("COL");
        histPrun[i_reg]->GetYaxis()->SetLabelOffset(0.03);
//         histPrun[i_reg]->GetXaxis()->SetLabelOffset(0.05);
        gPad->SetTopMargin(0);
//         gPad->SetBottomMargin(0);
//         gPad->SetLeftMargin(0);
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
//         myText(    0.1,1.+(20./mainHeight) ,1,histPrun[i_reg]->GetTitle());
        //
        pUp->cd();
        myText((leftSize+1.*i_reg*(regionSize+separation))/(leftSize+mainWidth),0.1 ,1,histPrun[i_reg]->GetTitle());
    }
    c->cd();
    TPad *pLo = new TPad("pLo","Pad Low",0,0,(1.*leftSize)/(leftSize+mainWidth),(1.*loSize)/(upSize+mainHeight+loSize));
    pLo->Draw();
    //
    c->cd();
    pUp->cd();
    myText(0.01,0.5,1,fLabel.c_str());
    //
    pLo->cd();
//     TLegend *leg = new TLegend(0.005,0,(1.*leftSize)/(leftSize+mainWidth),0.95);
    TLegend *leg = new TLegend(0.005,0,0.95,0.95);
    TH1F* hGray   = new TH1F("hGray"  ,"hGray"  ,1,0,1);    hGray->SetFillColor(kGray);         hGray->SetLineWidth(0);
    TH1F* hYellow = new TH1F("hYellow","hYellow",1,0,1);    hYellow->SetFillColor(kYellow);     hYellow->SetLineWidth(0);
    TH1F* hOrange = new TH1F("hOrange","hOrange",1,0,1);    hOrange->SetFillColor(kOrange-3);   hOrange->SetLineWidth(0);
    TH1F* hRed    = new TH1F("hRed"   ,"hRed"   ,1,0,1);    hRed->SetFillColor(kRed);           hRed->SetLineWidth(0);
    TH1F* hGreen  = new TH1F("hGreen" ,"hGree"  ,1,0,1);    hGreen->SetFillColor(8);            hGreen->SetLineWidth(0);
    TH1F* hBlue   = new TH1F("hBlue"  ,"hBlue"  ,1,0,1);    hBlue->SetFillColor(kBlue);         hBlue->SetLineWidth(0);
    TH1F* hPurple = new TH1F("hPurple","hPurple",1,0,1);    hPurple->SetFillColor(6);           hPurple->SetLineWidth(0);
    TH1F* hBlack  = new TH1F("hBlack" ,"hBlack" ,1,0,1);    hBlack->SetFillColor(kBlack);       hBlack->SetLineWidth(0);
    //char sysLarg [10];
    //sprintf(sysLarg,"#geq %3.f %%",fThresholdSystLarge*100);
    string sysLarg="Dropped as >"+std::to_string((int)(fThresholdSystLarge*100))+"%";
    leg->SetBorderSize(0);
    leg->SetMargin(0.1);
    leg->SetFillStyle(0);
    leg->AddEntry(hGray,"Not present","f");
    leg->AddEntry(hGreen,"Kept","f");
    leg->AddEntry(hYellow, "Shape dropped","f");
    leg->AddEntry(hOrange, "Norm. dropped","f");
    leg->AddEntry(hRed, "Dropped","f");
    if (fThresholdSystLarge > -1) {
        leg->AddEntry(hBlue  , sysLarg.c_str() ,"f");
        leg->AddEntry(hPurple, "Bad shape" ,"f");
        leg->AddEntry(hBlack , "Bad shape & norm." ,"f");
    }
    leg->SetTextSize(0.85*gStyle->GetTextSize());
    leg->Draw();
    //
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/Pruning"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );

    //
    // Save prunign hist for future usage
    TFile *filePrun;
    // - checking if Pruning.root exists
    // if yes
    if(!gSystem->AccessPathName( (fName+"/Pruning.root").c_str() )){
        // ...
        filePrun = new TFile( (fName+"/Pruning.root").c_str() );
    }
    else{
//     if(!filePrun->IsOpen()){
        filePrun = new TFile( (fName+"/Pruning.root").c_str(),"RECREATE" );
//         filePrun = new TFile( (fName+"/Pruning"+fSuffix+".root").c_str(),"RECREATE" );
        for(int i_reg=0;i_reg<(int)histPrun.size();i_reg++){
            histPrun_toSave[i_reg]->Write("",TObject::kOverwrite);
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::Fit(){

    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    RooDataSet* data = 0x0;
    RooWorkspace* ws = 0x0;

    //
    // If there's a workspace specified, go on with simple fit, without looking for separate workspaces per region
    //
    if(fWorkspaceFileName!=""){
        TFile *rootFile = new TFile(fWorkspaceFileName.c_str(),"read");
        ws = (RooWorkspace*) rootFile->Get("combined");
        if(!ws){
            std::cout << "<!> Error in TtHFit::Fit: The workspace (\"combined\") cannot be found in file " << fWorkspaceFileName << ". Please check !" << std::endl;
        }
        if(!fFitIsBlind && hasData) data = (RooDataSet*)ws->data("obsData");
        else                        data = (RooDataSet*)ws->data("asimovData");
        PerformFit( ws, data, fFitType, true );
    }
    else{
        //
        // Fills a vector of regions to consider for fit
        //
//         bool simpleData = true; // this becomes false if some regions have data and others don't
        int previousDataType = -1;
        std::vector < std:: string > regionsToFit;
        std::map < std::string, int > regionDataType;
        for( unsigned int i_ch = 0; i_ch < fNRegions; i_ch++ ){
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
                } else {
                    dataType = fRegions[i_ch] -> fRegionDataType;
                }
                regionDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , dataType) );
                //
//                 if(previousDataType>=0 && previousDataType!=(int)dataType) simpleData = false;
                previousDataType = (int)dataType;
            }
        }
        //
        // Creating the combined model with the regions to fit only
        //
        ws = PerformWorkspaceCombination( regionsToFit );
        //
        // If needed (only if needed), create a RooDataset object
        //
//         if(simpleData){
//             if(!fFitIsBlind && hasData) data = (RooDataSet*)ws->data("obsData");
//             else                        data = (RooDataSet*)ws->data("asimovData");
//         }
//         else{
            data = DumpData( ws, regionDataType, fFitNPValues, fFitPOIAsimov );
//         }
        //
        // Calls the PerformFit() function to actually do the fit
        //
        PerformFit( ws, data, fFitType, true);
    }

    //
    // Calls the  function to create LH scan with respect to a parameter
    //
    if(fVarNameLH.size()>0){
        if (fVarNameLH[0]=="all"){
            for(map<string,string>::iterator it=TtHFitter::SYSTMAP.begin(); it!=TtHFitter::SYSTMAP.end(); ++it){
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

//__________________________________________________________________________________
//
RooDataSet* TtHFit::DumpData( RooWorkspace *ws,  std::map < std::string, int > &regionDataType, std::map < std::string, double > &npValues, const double poiValue  ){

    //
    // This function dumps a RooDataSet object using the input informations provided by the user
    //    |-> Used when testing Fit response (inject one NP in data and check fit result)
    //    |-> Used when using fit results in some regions to generate Asimov data in blinded regions
    //
    if(TtHFitter::DEBUGLEVEL>0){
        std::cout << "=> In TtHFit::DumpData(): Dumping data with the following parameters" << std::endl;
        std::cout << "    * Regions data type " << std::endl;
        for( const std::pair < std::string, int > dataType : regionDataType ){
            std::cout << "       - Region: " << dataType.first << "       DataType: " << dataType.second << std::endl;
        }
        if(npValues.size()){
            std::cout << "    * Injected NP values " << std::endl;
            for ( const std::pair < std::string, double > npValue : npValues ){
                std::cout << "       - NP: " << npValue.first << "       Value: " << npValue.second << std::endl;
            }
        } else {
            std::cout << "    * No NP values injected " << std::endl;
        }
        std::cout << "    * POI value: " << poiValue << std::endl;
    }

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
            cout << "TtHFit::INFO: Activating binned likelihood attribute for " << arg->GetName() << endl;
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
    poi -> setVal(poiValue);

    //-- Nuisance parameters
    if (!fStatOnly){
        RooRealVar* var(nullptr);
        TIterator *npIterator = mc -> GetNuisanceParameters() -> createIterator();
        while( (var = (RooRealVar*) npIterator->Next()) ){
          std::map < std::string, double >::const_iterator it_npValue = npValues.find( var -> GetName() );
          if( it_npValue != npValues.end() ){
            var -> setVal(it_npValue -> second);
          }
        }
    }

    //Looping over regions
    map<string, RooDataSet*> asimovDataMap;
    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
    RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;
    int nrIndices = 0;
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
            std::cout << "=> In TtHFit::DumpData(): the following region is not specified in the inputs to the function (" << channelCat->getLabel() << "): use Asimov" << std::endl;
            std::cout << "   This SHOULD NOT HAPPEN ! Please check if everying is fine !" << std::endl;
        } else {
            dataType = regionDataType[channelCat->getLabel()];
        }

        //A protection: if there is no real observed data, use only ASIMOV (but print a warning)
        if(dataType==Region::REALDATA && !realData){
            std::cout << "=> In TtHFit::DumpData(): you want real data for channel " << channelCat->getLabel() << " but none is available in the workspace. Using Asimov instead." << std::endl;
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
            asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;

        } else if(dataType==Region::REALDATA) {
            RooAbsData *datatmp = realData->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),tt->GetName()));
            asimovDataMap[string(channelCat->getLabel())] = (RooDataSet*)datatmp;
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

    return asimovData;
}

//__________________________________________________________________________________
//
std::map < std::string, double > TtHFit::PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType, bool save ){

    std::map < std::string, double > result;

    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////

    //
    // Fit configuration (SPLUSB or BONLY)
    //
    FittingTool *fitTool = new FittingTool();
    fitTool -> SetDebug(TtHFitter::DEBUGLEVEL);
    if(fitType==BONLY){
        fitTool -> ValPOI(0.);
        fitTool -> ConstPOI(true);
    } else if(fitType==SPLUSB){
        fitTool -> ValPOI(fFitPOIAsimov);
        fitTool -> ConstPOI(false);
    }
    fitTool -> SetRandomNP(fRndRange, fUseRnd, fRndSeed);
    if(fStatOnly){
      fitTool -> NoGammas();
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
        std::cout << "Setting the variables to use MINOS with" << std::endl;
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
    RooDataSet* data = 0x0;
    if(inputData){
        data = inputData;
    } else {
        std::cout << "In TtHFit::PerformFit() function: you didn't provide inputData => will use the observed data !" << std::endl;
        data = (RooDataSet*)ws->data("obsData");
        if(data==0x0){
            std::cout << "In TtHFit::PerformFit() function: no observedData found => will use the Asimov data !" << std::endl;
            data = (RooDataSet*)ws->data("asimovData");
        }
    }

    //
    // For stat-only fit on data:
    // - read fit resutls
    // - fix all NP to fitted ones before fitting
    if(fStatOnlyFit){
        std::cout << "Fitting stat-only: reading fit results from full fit from file: " << std::endl;
        std::cout << "  " << (fName+"/Fits/"+fInputName+fSuffix+".txt") << std::endl;
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

    // FixNP -- Roger
    if(fFitFixedNPs.size()>0){
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(auto nuisParToFix : fFitFixedNPs){
            npNames.push_back( nuisParToFix.first );
            npValues.push_back( nuisParToFix.second );
        }
        fitTool -> FixNPs(npNames,npValues);
    }

    //
    // Get initial ikelihood value from Asimov
    float nll0 = 0.;
    if(fGetGoodnessOfFit) nll0 = fitTool -> FitPDF( mc, simPdf, (RooDataSet*)ws->data("asimovData"), false, true );
    
    //
    // Get number of degrees of freedom
    // - number of bins
    int ndof = inputData->numEntries();
    // - minus number of free & non-constant parameters
    int nNF = 0;
    for(int i_nf=0;i_nf<fNNorm;i_nf++){
        if(fNormFactors[i_nf]->fConst) continue;
        if(fFitType==BONLY && fPOI==fNormFactors[i_nf]->fName) continue;
        nNF++;
    }
    ndof -= nNF;
    
//     RooArgSet* nuis = (RooArgSet*) mc->GetPdf()->getParameters(*data);
// //     RooArgSet* syst = (RooArgSet*) mc->GetPdf()->getAllConstraints(mc->GetObservables(),*nuis);
// //     RooArgSet* syst = (RooArgSet*) mc->GetConstraintParameters();
//     RooStats::RemoveConstantParameters(nuis);
// //     if(syst) RooStats::RemoveConstantParameters(syst);
//     cerr << ndof << endl;
//     cerr << nuis->getSize() << endl;
// //     if(syst) cerr << syst->getSize() << endl;
//     cerr << syst << endl;
//     ndof -= (nuis->getSize()-syst->getSize());
//     cerr << ndof << endl;
    
    // Performs the fit
    fitTool -> MinimType("Minuit2");
    float nll = fitTool -> FitPDF( mc, simPdf, data );
//     fitTool -> FitPDF( mc, simPdf, data, true );   // for fast fit
    if(save){
        gSystem -> mkdir((fName+"/Fits/").c_str(),true);
        if(fStatOnlyFit) fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_statOnly.txt");
        else             fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+".txt");
    }
    result = fitTool -> ExportFitResultInMap();

    //
    // Goodness of fit
    if(fGetGoodnessOfFit){
        cout << "----------------------- -------------------------- -----------------------" << endl;
        cout << "----------------------- GOODNESS OF FIT EVALUATION -----------------------" << endl;
        cout << "  NLL0        = " << nll0 << endl;
        cout << "  NLL         = " << nll  << endl;
        cout << "  ndof        = " << ndof << endl;
        float deltaNLL = nll-nll0;
        cout << "  dNLL        = " << deltaNLL << endl;
        cout << "  2dNLL/nof   = " << 2.*deltaNLL/ndof << endl;
        double prob = ROOT::Math::chisquared_cdf_c( 2* deltaNLL, ndof);
        cout << "  probability = " << prob << endl;
        cout << "----------------------- -------------------------- -----------------------" << endl;
        cout << "----------------------- -------------------------- -----------------------" << endl;
    }
    
    return result;
}

//__________________________________________________________________________________
//
RooWorkspace* TtHFit::PerformWorkspaceCombination( std::vector < std::string > &regionsToFit ){

    //
    // Definition of the fit regions
    //
    std::vector < RooWorkspace* > vec_ws;
    std::vector < std::string > vec_chName;
    RooStats::HistFactory::Measurement *measurement = 0;
    //
    // Take the measurement from the combined workspace, to be sure to have all the systematics (even the ones which are not there in the first region)
    TFile *rootFileCombined = new TFile( (fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read");
    if(rootFileCombined!=0x0) measurement = (RooStats::HistFactory::Measurement*) rootFileCombined -> Get( (fInputName+fSuffix).c_str());
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
            TFile *rootFile = new TFile(fileName.c_str(),"read");
            RooWorkspace* m_ws = (RooWorkspace*) rootFile->Get((fRegions[i_ch]->fName).c_str());
            if(!m_ws){
                std::cout << "<!> Error in TtHFit::PerformWorkspaceCombination: The workspace (\"" << fRegions[i_ch] -> fName << "\") cannot be found in file " << fileName << ". Please check !" << std::endl;
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
        std::cout << "<!> Error in TtHFit::PerformWorkspaceCombination() : The measurement object has not been retrieved ! Please check." << std::endl;
        return 0;
    }
    RooStats::HistFactory::HistoToWorkspaceFactoryFast factory(*measurement);

    // Creating the combined model
    RooWorkspace* ws = factory.MakeCombinedModel( vec_chName, vec_ws );

    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws, *measurement );

    return ws;
}

//__________________________________________________________________________________
//
void TtHFit::PlotFittedNP(){
    if(fStatOnly || fStatOnlyFit){
        std::cout << "TtHFit::INFO: Stat only fit => No NP Pull plots generated." << std::endl;
        return;
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
        for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++){
          fFitResults->DrawNPPulls(fName+"/NuisPar"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format],"all",fNormFactors);
          fFitResults->DrawGammaPulls(fName+"/Gammas"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]);
          fFitResults->DrawNormFactors(fName+"/NormFactors"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format],fNormFactors);
        }
        if(npCategories.size()>1){
            for( const std::string cat : npCategories ){
                std::string cat_for_name = cat;
                std::replace( cat_for_name.begin(), cat_for_name.end(), ' ', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '#', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '{', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '}', '_');
                for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++){
                  fFitResults->DrawNPPulls(fName+"/NuisPar_"+cat_for_name+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format],cat,fNormFactors);
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::PlotCorrelationMatrix(){
    if(fStatOnly || fStatOnlyFit){
        std::cout << "TtHFit::INFO: Stat only fit => No Correlation Matrix generated." << std::endl;
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TtHFitter::CORRELATIONTHRESHOLD)
    ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
    if(fFitResults){
        fFitResults->fNuisParToHide = fVarNameHide;
        for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
            fFitResults->DrawCorrelationMatrix(fName+"/CorrMatrix"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format],TtHFitter::CORRELATIONTHRESHOLD);
    }
}

//__________________________________________________________________________________
//
void TtHFit::GetLimit(){

    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    string cmd;

    //
    // If a workspace file name is specified, do simple limit
    //
    if(fWorkspaceFileName!=""){
        string dataName = "obsData";
        if(!hasData || fLimitIsBlind) dataName = "asimovData";
        if(fSignalInjection)
            cmd = "root -l -b -q 'runAsymptoticsCLs_inject.C+(\""+fWorkspaceFileName+"\",\"combined\",\"ModelConfig\",\""+dataName+"\",\"asimovData_0\",\""+fName+"/Limits/\",\""+fInputName+fSuffix+"\",0.95)'";
        else
            cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+fWorkspaceFileName+"\",\"combined\",\"ModelConfig\",\""+dataName+"\",\"asimovData_0\",\""+fName+"/Limits/\",\""+fInputName+fSuffix+"\",0.95)'";

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
        for( unsigned int i_ch = 0; i_ch < fNRegions; i_ch++ ){
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
            RooWorkspace* ws_forFit = PerformWorkspaceCombination( regionsForFit );

            //
            // Calls the PerformFit() function to actually do the fit
            //
            npValues = PerformFit( ws_forFit, data, FitType::BONLY);
        }

        //
        // Create the final asimov dataset for limit setting
        //
        RooWorkspace* ws_forLimit = PerformWorkspaceCombination( regionsForLimit );
        data = DumpData( ws_forLimit, regionsForLimitDataType, npValues, npValues.find(fPOI)==npValues.end() ? fLimitPOIAsimov : npValues[fPOI] );

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
        if(fSignalInjection)
            cmd = "root -l -b -q 'runAsymptoticsCLs_inject.C+(\""+(string)outputName+"\",\"combined\",\"ModelConfig\",\"ttHFitterData\",\"asimovData_0\",\""+fName+"/Limits/\",\""+fInputName+fSuffix+"\",0.95)'";
        else
            cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+(string)outputName+"\",\"combined\",\"ModelConfig\",\"ttHFitterData\",\"asimovData_0\",\""+fName+"/Limits/\",\""+fInputName+fSuffix+"\",0.95)'";
    }

    //
    // Finally computing the limit
    //
    gSystem->Exec(cmd.c_str());
}

//__________________________________________________________________________________
//
void TtHFit::GetSignificance(){

    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }

    string cmd;

    //
    // If a workspace file name is specified, do simple significance
    //
    if(fWorkspaceFileName!=""){
        string dataName = "obsData";
        if(!hasData || fFitIsBlind) dataName = "asimovData";
        cmd = "root -l -b -q 'runSig.C(\""+fWorkspaceFileName+"\",\"combined\",\"ModelConfig\",\""+dataName+"\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+fSuffix+"\",\""+fName+"/Significance\")'";
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
        for( unsigned int i_ch = 0; i_ch < fNRegions; i_ch++ ){
            if( fRegions[i_ch] -> fRegionType == Region::VALIDATION ) continue;
            if( hasData && fRegions[i_ch] -> fRegionDataType == Region::REALDATA && !fLimitIsBlind){
                Region::DataType dataType = fRegions[i_ch] -> fRegionDataType;
                regionsForFit.push_back( fRegions[i_ch] -> fName );
                regionsForFitDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , dataType) );
            }
            regionsForSign.push_back(fRegions[i_ch] -> fName);
            regionsForSignDataType.insert( std::pair < std::string, int >(fRegions[i_ch] -> fName , (!hasData || fLimitIsBlind) ? Region::ASIMOVDATA : fRegions[i_ch] -> fRegionDataType) );
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
            RooWorkspace* ws_forFit = PerformWorkspaceCombination( regionsForFit );

            //
            // Calls the PerformFit() function to actually do the fit
            //
            npValues = PerformFit( ws_forFit, data, FitType::BONLY);
        }

        //
        // Create the final asimov dataset for limit setting
        //
        RooWorkspace* ws_forSignificance = PerformWorkspaceCombination( regionsForSign );
        data = DumpData( ws_forSignificance, regionsForSignDataType, npValues, npValues.find(fPOI)==npValues.end() ? fLimitPOIAsimov : npValues[fPOI] );

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
        cmd = "root -l -b -q 'runSig.C(\""+(string)outputName+"\",\"combined\",\"ModelConfig\",\"ttHFitterData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fInputName+fSuffix+"\",\""+fName+"/Significance\")'";
    }

    gSystem->Exec(cmd.c_str());
}

//__________________________________________________________________________________
//
void TtHFit::ReadFitResults(string fileName){
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Reading fit results from file " << fileName << std::endl;
    fFitResults = new FitResults();
    if(fileName.find(".txt")!=string::npos){
        fFitResults->ReadFromTXT(fileName);
    }
    // make a list of systematics from all samples...
    // ...
    // assign to each NP in the FitResults a title, and a category according to the syst in the fitter
    for(unsigned int i=0;i<fFitResults->fNuisPar.size();i++){
        for(unsigned int j=0;j<fSystematics.size();j++){
            if(fSystematics[j]->fName == fFitResults->fNuisPar[i]->fName){
                fFitResults->fNuisPar[i]->fTitle = fSystematics[j]->fTitle;
                fFitResults->fNuisPar[i]->fCategory = fSystematics[j]->fCategory;
            }
        }
        for(unsigned int j=0;j<fNormFactors.size();j++){
            if(fNormFactors[j]->fName == fFitResults->fNuisPar[i]->fName){
                fFitResults->fNuisPar[i]->fTitle = fNormFactors[j]->fTitle;
                fFitResults->fNuisPar[i]->fCategory = fNormFactors[j]->fCategory;
            }
        }
	// FIXME SF probably there are several NPs associated to it 
        for(unsigned int j=0;j<fShapeFactors.size();j++){
            if(fShapeFactors[j]->fName == fFitResults->fNuisPar[i]->fName){
                fFitResults->fNuisPar[i]->fTitle = fShapeFactors[j]->fTitle;
                fFitResults->fNuisPar[i]->fCategory = fShapeFactors[j]->fCategory;
            }
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::Print(){
    std::cout << endl;
    std::cout << "-------------------------------------------" << endl;
    std::cout << "  TtHFit: " << fInputName << endl;
    std::cout << "      NtuplePaths ="; for(int i=0;i<(int)fNtuplePaths.size();i++) std::cout << " " << fNtuplePaths[i] << std::endl;
    std::cout << "      NtupleName  =";   std::cout << " " << fNtupleName << std::endl;
    std::cout << "      MCweight    =";   std::cout << " " << fMCweight   << std::endl;
    std::cout << "      Selection   =";   std::cout << " " << fSelection  << std::endl;
    std::cout << "      HistoPaths  ="; for(int i=0;i<(int)fHistoPaths.size();i++) std::cout << " " << fHistoPaths[i] << std::endl;
    std::cout << "      HistoName   =";   std::cout << " " << fHistoName << endl;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->Print();
    }
    std::cout << endl;
    std::cout << "-------------------------------------------" << std::endl;
}

//__________________________________________________________________________________
//
Region* TtHFit::GetRegion(string name){
    for(unsigned int i=0;i<fRegions.size();i++){
        if(fRegions[i]->fName == name) return fRegions[i];
    }
    return 0x0;
}

//__________________________________________________________________________________
//
Sample* TtHFit::GetSample(string name){
    for(unsigned int i=0;i<fSamples.size();i++){
        if(fSamples[i]->fName == name) return fSamples[i];
    }
    return 0x0;
}

//__________________________________________________________________________________
//
void TtHFit::DrawAndSaveSeparationPlots(){
//     std::cout << "In DrawAndSaveSeparationPlots:" << endl;

    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    gSystem->mkdir((fName+"/Plots/Separation").c_str());


    // loop over regions
    for(unsigned int i_ch=0; i_ch < fRegions.size(); i_ch++){
        // begin plotting
        TCanvas* dummy3 = new TCanvas("dummy3", "dummy3", 600,600);
        dummy3->cd();

        if(fRegions[i_ch]->fNSig==0){
            std::cout << "ERROR::TtHFit::DrawAndSaveSeparationPlots: No Signal Found" << std::endl;
            continue;
        }

        TH1F* sig = (TH1F*)fRegions[i_ch]->fSig[0]->fHist->Clone();

        TH1F* bkg = (TH1F*)fRegions[i_ch]->fBkg[0]->fHist->Clone(); // clone the first bkg
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

//         TLegend *legend3=new TLegend(0.48,0.72,0.94,0.87);
        TLegend *legend3=new TLegend(0.55,0.77,0.94,0.87);
//         legend3->SetTextFont(42);
//         legend3->SetTextSize(0.043);
        legend3->SetTextFont(gStyle->GetTextFont());
        legend3->SetTextSize(gStyle->GetTextSize());
        legend3->AddEntry(bkg, "Total background" , "l");
//         legend3->AddEntry(sig, "t#bar{t}H (m_{H} = 125 GeV)" , "l");
        legend3->AddEntry(sig, fRegions[i_ch]->fSig[0]->fSample->fTitle.c_str() , "l");
//         legend3->SetFillColor(0) ;
//         legend3->SetLineColor(0) ;
        legend3->SetFillStyle(0) ;
        legend3->SetBorderSize(0);

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

        legend3->Draw("same");

//         std::string identS = fRegions[i_ch]->fLabel;
//         TLatex ls;
//         ls.SetNDC();
//         ls.SetTextSize(0.03);
//         ls.SetTextColor(kBlack);
//         ls.SetTextFont(42);
//         ls.DrawLatex(0.20,0.73,identS.c_str());

//         TLatex ls2;
//         ls2.SetNDC();
//         ls2.SetTextSize(0.03);
//         ls2.SetTextColor(kBlack);
//         ls2.SetTextFont(62);
// //         ls2.DrawLatex(0.20,0.78,"Single lepton");
//         ls2.DrawLatex(0.20,0.78,fLabel.c_str());
        myText(0.20,0.78,1,fLabel.c_str());
        myText(0.20,0.73,1,fRegions[i_ch]->fLabel.c_str());

        std::string cme = fRegions[i_ch]->fCmeLabel;
        std::string lumi = fRegions[i_ch]->fLumiLabel;

//         TLatex ls3;
//         ls3.SetNDC();
//         ls3.SetTextSize(0.03);
//         ls3.SetTextColor(kBlack);
//         ls3.SetTextFont(42);
//         ls3.DrawLatex(0.20,0.83, Form("#sqrt{s} = %s, %s", cme.c_str(), lumi.c_str()));
        myText(0.20,0.83,1,Form("#sqrt{s} = %s, %s", cme.c_str(), lumi.c_str()));

//         ATLASLabelNew(0.20, 0.90,(char*)" Internal Simulation",kBlack, 0.03);
        if(fAtlasLabel!="none") ATLASLabelNew(0.20,0.84+0.04,(char*)(fAtlasLabel+"  Simulation").c_str(), kBlack, gStyle->GetTextSize());
        
//         TLatex ls4;
//         ls4.SetNDC();
//         ls4.SetTextSize(0.03);
//         ls4.SetTextColor(kBlack);
//         ls4.SetTextFont(42);
        ostringstream SEP;
        SEP.precision(3);
        SEP << "Separation: " << GetSeparation(sig,bkg)*100 << "%";
//         ls4.DrawLatex(0.20, 0.69, SEP.str().c_str());
        myText(0.55,0.73,1,SEP.str().c_str());

        for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
            dummy3->SaveAs((fName+"/Plots/Separation/"+fRegions[i_ch]->fName+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format] ).c_str());

    }// regions

   return;
}

//____________________________________________________________________________________
//
void TtHFit::ProduceNPRanking( string NPnames/*="all"*/ ){

    if(fFitType==BONLY){
        std::cerr << "\033[1;31m<!> ERROR in TtHFit::ProduceNPRanking(): For ranking plots, the SPLUSB FitType is needed.  \033[0m"<<std::endl;
        return;
    }

    //
    // List of systematics to check
    //
    std::vector< string > nuisPars;
    std::vector< bool > isNF;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(NPnames=="all" || NPnames==fSystematics[i_syst]->fName ||
            ( atoi(NPnames.c_str())==i_syst && (atoi(NPnames.c_str())>0 || strcmp(NPnames.c_str(),"0")==0 ) )
            ){
	    if(fSystematics[i_syst]->fType == Systematic::SHAPE) continue;
            nuisPars.push_back( fSystematics[i_syst]->fName );
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
    //Text files containing information necessary for drawing of ranking plot
    //
    string outName = fName+"/Fits/NPRanking"+fSuffix;
    if(NPnames!="all") outName += "_"+NPnames;
    outName += ".txt";
    ofstream outName_file(outName.c_str());
    //
    float central;
    float up;
    float down;
    float muhat;
    float dMuUp, dMuDown;
    std::map< string,float > muVarUp;
    std::map< string,float > muVarDown;
    std::map< string,float > muVarNomUp;
    std::map< string,float > muVarNomDown;

    //
    // Fills a vector of regions to consider for fit
    //
    std::vector < std:: string > regionsToFit;
    std::map < std::string, int > regionDataType;
    for( unsigned int i_ch = 0; i_ch < fNRegions; i_ch++ ){
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
    RooWorkspace* ws = PerformWorkspaceCombination( regionsToFit );

    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());
    RooDataSet* data = DumpData( ws, regionDataType, fFitNPValues, fFitPOIAsimov );

    // Loop on NPs to find gammas and add to the list to be ranked
//     if(NPnames=="all" || NPnames.find("gamma_stat")!=string::npos){
    if(NPnames=="all" || NPnames.find("gamma")!=string::npos){
        RooRealVar* var = NULL;
        RooArgSet* nuis = (RooArgSet*) mc->GetNuisanceParameters();
        if(nuis){
            TIterator* it2 = nuis->createIterator();
            while( (var = (RooRealVar*) it2->Next()) ){
                string np = var->GetName();
                if(np.find("gamma")!=string::npos){
            // 	    if(np==NPnames || NPnames=="all") nuisPars.push_back(np);
                    if(np==NPnames || NPnames=="all") nuisPars.push_back(ReplaceString(np,"gamma_",""));
                    isNF.push_back( true );
                    if(NPnames!="all") break;
                }
// 	  if(NPnames!="all"){
// 	    if(np==NPnames){
// 	      nuisPars.push_back(ReplaceString(np,"gamma_",""));
// 	      isNF.push_back( true );
// 	      break;
// 	    }
// 	  }
// 	  else if(np.find("gamma_stat")!=string::npos){
// 	    nuisPars.push_back(np);
// 	    isNF.push_back( true );
// 	  }
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
    fitTool -> SetDebug(TtHFitter::DEBUGLEVEL);
    fitTool -> ValPOI(1.);
    fitTool -> ConstPOI(false);
    if(fStatOnly){
      fitTool -> NoGammas();
      fitTool -> NoSystematics();
    }
    //
    // no point is setting Minos for ranking...
//     if(fVarNameMinos.size()>0){
//         std::cout << "Setting the variables to use MINOS with" << std::endl;
//         fitTool -> UseMinos(fVarNameMinos);
//     }

    ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".txt");
    muhat = fFitResults -> GetNuisParValue( fPOI );
    //if(!hasData) muhat = 1.;  // FIXME -> Loic: Do we actually need that ?

    for(unsigned int i=0;i<nuisPars.size();i++){
        //
        // Getting the postfit values of the nuisance parameter
        central = fFitResults -> GetNuisParValue(   nuisPars[i] );
        up      = fFitResults -> GetNuisParErrUp(   nuisPars[i] );
        down    = fFitResults -> GetNuisParErrDown( nuisPars[i] );
	//// Thomas : We should be careful with changing naming convention compared to RooFit !!
	// TtHFitter store gammas names as stat_Reg_bin_i (i.e. remove the gamma_ at the beginning)
	// Now there is no real identifier in the NP name to state if it is a gamma or not and add back gamma_ except this _bin_
	if( (NPnames=="all" && nuisPars[i].find("_bin_")!=string::npos)
	    || NPnames.find("gamma_stat")!=string::npos){
	  string tmpNuispar=nuisPars[i];
	  nuisPars[i]="gamma_"+tmpNuispar;
	}
        outName_file <<  nuisPars[i] << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";
        //
        // Experimental: reduce the range of ranking
        if(TtHFitter::OPTION["ReduceRanking"]!=0){
            up   *= TtHFitter::OPTION["ReduceRanking"];
            down *= TtHFitter::OPTION["ReduceRanking"];
        }
        //
        // Set the NP to its post-fit *up* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool -> ResetFixedNP();
        fitTool -> FixNP( nuisPars[i], central + TMath::Abs(up  ) );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
        //
        // Set the NP to its post-fit *down* variation and refit to get the fitted POI
        ws->loadSnapshot("tmp_snapshot");
        fitTool -> ResetFixedNP();
        fitTool -> FixNP( nuisPars[i], central - TMath::Abs(down) );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        //
        dMuUp   = muVarUp[nuisPars[i]]-muhat;
        dMuDown = muVarDown[nuisPars[i]]-muhat;
        //
        // Experimental: reduce the range of ranking
        if(TtHFitter::OPTION["ReduceRanking"]!=0){
            dMuUp   /= TtHFitter::OPTION["ReduceRanking"];
            dMuDown /= TtHFitter::OPTION["ReduceRanking"];
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
            if(TtHFitter::OPTION["ReduceRanking"]!=0){
                up   *= TtHFitter::OPTION["ReduceRanking"];
                down *= TtHFitter::OPTION["ReduceRanking"];
            }
            //
            // Set the NP to its pre-fit *up* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool -> ResetFixedNP();
            fitTool -> FixNP( nuisPars[i], central + TMath::Abs(up  ) );
            fitTool -> FitPDF( mc, simPdf, data );
            muVarNomUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
            //
            // Set the NP to its pre-fit *down* variation and refit to get the fitted POI (pre-fit impact on POI)
            ws->loadSnapshot("tmp_snapshot");
            fitTool -> ResetFixedNP();
            fitTool -> FixNP( nuisPars[i], central - TMath::Abs(down) );
            fitTool -> FitPDF( mc, simPdf, data );
            //
            muVarNomDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
            if(fReduceNPforRanking!=0){
                muVarNomUp[ nuisPars[i] ]   = muhat + (muVarNomUp[ nuisPars[i] ]   - muhat)/fReduceNPforRanking;
                muVarNomDown[ nuisPars[i] ] = muhat + (muVarNomDown[ nuisPars[i] ] - muhat)/fReduceNPforRanking;
            }
        }
        dMuUp   = muVarNomUp[nuisPars[i]]-muhat;
        dMuDown = muVarNomDown[nuisPars[i]]-muhat;
        //
        // Experimental: reduce the range of ranking
        if(TtHFitter::OPTION["ReduceRanking"]!=0){
            dMuUp   /= TtHFitter::OPTION["ReduceRanking"];
            dMuDown /= TtHFitter::OPTION["ReduceRanking"];
        }
        //
       outName_file << dMuUp << "   " << dMuDown << " "<<endl;

    }
    outName_file.close();
    ws->loadSnapshot("tmp_snapshot");

//     //
//     // Creating the rootfile
//     //
//     TFile *f_clone = new TFile( "ws_forRanking.root", "recreate" );
//     ws -> import(*data,Rename("ttHFitterData"));
//     ws -> Write();
//     f_clone -> Close();

}

//____________________________________________________________________________________
//
void TtHFit::PlotNPRankingManager(){
  if(fRankingPlot=="Merge"  || fRankingPlot=="all") PlotNPRanking(true,true);
  if(fRankingPlot=="Systs"  || fRankingPlot=="all") PlotNPRanking(true,false);
  if(fRankingPlot=="Gammas" || fRankingPlot=="all") PlotNPRanking(false,true);

}

//____________________________________________________________________________________
//
void TtHFit::PlotNPRanking(bool flagSysts, bool flagGammas){
    //
    string fileToRead = fName+"/Fits/NPRanking"+fSuffix+".txt";
    //
    // trick to merge the ranking outputs produced in parallel:
    string cmd = " if [[ `ls "+fName+"/Fits/NPRanking"+fSuffix+"_*` != \"\" ]] ; ";
    cmd       += " then rm "+fileToRead+" ; ";
    cmd       += " cat "+fName+"/Fits/NPRanking"+fSuffix+"_* > "+fileToRead+" ; ";
    cmd       += " fi ;";
    gSystem->Exec(cmd.c_str());
    //
    int maxNP = fRankingMaxNP;
    //
    string paramname;
    double nuiphat;
    double nuiperrhi;
    double nuiperrlo;
    double PoiUp;
    double PoiDown;
    double PoiNomUp;
    double PoiNomDown;
    std::vector<string> parname;
    std::vector<double> nuhat;
    std::vector<double> nuerrhi;
    std::vector<double> nuerrlo;
    std::vector<double> poiup;
    std::vector<double> poidown;
    std::vector<double> poinomup;
    std::vector<double> poinomdown;
    std::vector<double> number;

    ifstream fin( fileToRead.c_str() );
    fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    if (paramname=="Luminosity"){
        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    }
    while (!fin.eof()){
        if(paramname.find("gamma")!=string::npos && !flagGammas){
	  fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
	  if (paramname=="Luminosity"){
	    fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
	  }
	  continue;
	}
	if(paramname.find("gamma")==string::npos && !flagSysts){
	  fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
	  if (paramname=="Luminosity"){
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
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        }
    }

    unsigned int SIZE = parname.size();
    if(TtHFitter::DEBUGLEVEL>0) std::cout << "NP ordering..." << std::endl;
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
                    swap(poiup[i],poiup[i-j]);
                    swap(poidown[i],poidown[i-j]);
                    swap(poinomup[i],poinomup[i-j]);
                    swap(poinomdown[i],poinomdown[i-j]);
                    swap(nuhat[i],nuhat[i-j]);
                    swap(nuerrhi[i],nuerrhi[i-j]);
                    swap(nuerrlo[i],nuerrlo[i-j]);
                    swap(parname[i],parname[i-j]);
                    index=i-j;
                }
                else{
                    swap(poiup[index],poiup[i-j]);
                    swap(poidown[index],poidown[i-j]);
                    swap(poinomup[index],poinomup[i-j]);
                    swap(poinomdown[index],poinomdown[i-j]);
                    swap(nuhat[index],nuhat[i-j]);
                    swap(nuerrhi[index],nuerrhi[i-j]);
                    swap(nuerrlo[index],nuerrlo[i-j]);
                    swap(parname[index],parname[i-j]);
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
    for (int i=0;i<SIZE;i++) {
        poimax = TMath::Max(poimax,TMath::Max( TMath::Abs(poiup[i]),TMath::Abs(poidown[i]) ));
        poimax = TMath::Max(poimax,TMath::Max( TMath::Abs(poinomup[i]),TMath::Abs(poinomdown[i]) ));
        nuerrlo[i] = TMath::Abs(nuerrlo[i]);
    }
    poimax *= 1.2;

    for (int i=0;i<SIZE;i++) {
        poiup[i]     *= (2./poimax);
        poidown[i]   *= (2./poimax);
        poinomup[i]  *= (2./poimax);
        poinomdown[i]*= (2./poimax);
    }

    // Resttrict to the first N
    if(SIZE>maxNP) SIZE = maxNP;

    // Graphical part - rewritten taking DrawPulls in TtHFitter
    float lineHeight  =  30;
    float offsetUp    =  60; // external
    float offsetDown  =  60;
    float offsetUp1   = 100; // internal
    float offsetDown1 =  15;
    int offset = offsetUp + offsetDown + offsetUp1 + offsetDown1;
    int newHeight = offset + SIZE*lineHeight;

    float xmin = -2;
    float xmax =  2;
    float max  =  0;

    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    TGraphAsymmErrors *g1 = new TGraphAsymmErrors();
    TGraphAsymmErrors *g2 = new TGraphAsymmErrors();
    TGraphAsymmErrors *g1a = new TGraphAsymmErrors();
    TGraphAsymmErrors *g2a = new TGraphAsymmErrors();

    int idx = 0;
    std::vector< string > Names;
    Names.clear();
    string parTitle;

    for(unsigned int i = parname.size()-SIZE; i<parname.size(); ++i){
        g->SetPoint(      idx, nuhat[i],idx+0.5);
        g->SetPointEXhigh(idx, nuerrhi[i]);
        g->SetPointEXlow( idx, nuerrlo[i]);

        g1->SetPoint(      idx, 0.,idx+0.5);
        g1->SetPointEXhigh(idx, poiup[i]);
        g1->SetPointEXlow( idx, 0.);
        g1->SetPointEYhigh(idx, 0.4);
        g1->SetPointEYlow( idx, 0.4);

        g2->SetPoint(      idx, 0.,idx+0.5);
        g2->SetPointEXhigh(idx, poidown[i]);
        g2->SetPointEXlow( idx, 0.);
        g2->SetPointEYhigh(idx, 0.4);
        g2->SetPointEYlow( idx, 0.4);

        g1a->SetPoint(      idx, 0.,idx+0.5);
        g1a->SetPointEXhigh(idx, poinomup[i]);
        g1a->SetPointEXlow( idx, 0.);
        g1a->SetPointEYhigh(idx, 0.4);
        g1a->SetPointEYlow( idx, 0.4);

        g2a->SetPoint(      idx, 0.,idx+0.5);
        g2a->SetPointEXhigh(idx, poinomdown[i]);
        g2a->SetPointEXlow( idx, 0.);
        g2a->SetPointEYhigh(idx, 0.4);
        g2a->SetPointEYlow( idx, 0.4);

	if(parname[i].find("gamma")!=string::npos){
	  string tmpTitle=parname[i];
	  tmpTitle=ReplaceString(tmpTitle,"gamma_stat_","");
	  tmpTitle=ReplaceString(tmpTitle,"_"," ");
	  parTitle="#gamma ("+tmpTitle+")";
	}
        else parTitle = TtHFitter::SYSTMAP[ parname[i] ];

        Names.push_back(parTitle);

        idx ++;
        if(idx > max)  max = idx;
    }

    TCanvas *c = new TCanvas("c","c",600,newHeight);
    c->SetTicks(0,0);
    gPad->SetLeftMargin(0.4);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1F *h_dummy = new TH1F("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy->SetMaximum( SIZE + offsetUp1/lineHeight   );
    h_dummy->SetMinimum(      - offsetDown1/lineHeight );
    h_dummy->SetLineWidth(0);
    h_dummy->SetFillStyle(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->SetFillColor(kWhite);
    h_dummy->GetYaxis()->SetLabelSize(0);
    h_dummy->Draw();
    h_dummy->GetYaxis()->SetNdivisions(0);

    g1->SetFillColor(kAzure-4);
    g2->SetFillColor(kCyan);
    g1->SetLineColor(g1->GetFillColor());
    g2->SetLineColor(g2->GetFillColor());

    g1a->SetFillColor(kWhite);
    g2a->SetFillColor(kWhite);
    g1a->SetLineColor(kAzure-4);
    g2a->SetLineColor(kCyan);
    g1a->SetFillStyle(0);
    g2a->SetFillStyle(0);
    g1a->SetLineWidth(1);
    g2a->SetLineWidth(1);

    g->SetLineWidth(2);

    g1a->Draw("5 same");
    g2a->Draw("5 same");
    g1->Draw("2 same");
    g2->Draw("2 same");
    g->Draw("p same");

    TLatex *systs = new TLatex();
    systs->SetTextAlign(32);
    systs->SetTextSize( systs->GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs->DrawLatex(xmin-0.1,i+0.5,Names[i].c_str());
    }
    h_dummy->GetXaxis()->SetLabelSize( h_dummy->GetXaxis()->GetLabelSize()*0.9 );
    h_dummy->GetXaxis()->CenterTitle();
    h_dummy->GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy->GetXaxis()->SetTitleOffset(1.2);

    TGaxis *axis_up = new TGaxis( -2, SIZE + (offsetUp1)/lineHeight, 2, SIZE + (offsetUp1)/lineHeight, -poimax,poimax, 510, "-" );
    axis_up->SetLabelOffset( 0.01 );
    axis_up->SetLabelSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetLabelFont(   gStyle->GetTextFont() );
    axis_up->Draw();
    axis_up->CenterTitle();
    axis_up->SetTitle("#Delta#mu");
    if(SIZE==20) axis_up->SetTitleOffset(1.5);
    axis_up->SetTitleSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetTitleFont(   gStyle->GetTextFont() );

    TPad *pad1 = new TPad("p1","Pad High",0,(newHeight-offsetUp-offsetUp1)/newHeight,0.4,1);
    pad1->Draw();

    pad1->cd();
    TLegend *leg1 = new TLegend(0.02,0.7,1,1.0,"Pre-fit impact on #mu:");
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetMargin(0.33);
    leg1->SetNColumns(2);
    leg1->SetTextFont(gStyle->GetTextFont());
    leg1->SetTextSize(gStyle->GetTextSize());
    leg1->AddEntry(g1a,"#theta_{0}=+#Delta#theta","f");
    leg1->AddEntry(g2a,"#theta_{0}=-#Delta#theta","f");
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.02,0.32,1,0.62,"Post-fit impact on #mu:");
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.33);
    leg2->SetNColumns(2);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetTextSize(gStyle->GetTextSize());
    leg2->AddEntry(g1,"#theta_{0}=+#Delta#hat{#theta}","f");
    leg2->AddEntry(g2,"#theta_{0}=-#Delta#hat{#theta}","f");
    leg2->Draw();

    TLegend *leg0 = new TLegend(0.02,0.1,1,0.25);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetMargin(0.2);
    leg0->SetTextFont(gStyle->GetTextFont());
    leg0->SetTextSize(gStyle->GetTextSize());
    leg0->AddEntry(g,"Nuis. Param. Pull","lp");
    leg0->Draw();

    c->cd();

    TLine l0;
    TLine l1;
    TLine l2;
    l0 = TLine(0,- offsetDown1/lineHeight,0,SIZE+0.5);// + offsetUp1/lineHeight);
    l0.SetLineStyle(kDashed);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    l1 = TLine(-1,- offsetDown1/lineHeight,-1,SIZE+0.5);// + offsetUp1/lineHeight);
    l1.SetLineStyle(kDashed);
    l1.SetLineColor(kBlack);
    l1.Draw("same");
    l2 = TLine(1,- offsetDown1/lineHeight,1,SIZE+0.5);// + offsetUp1/lineHeight);
    l2.SetLineStyle(kDashed);
    l2.SetLineColor(kBlack);
    l2.Draw("same");

    ATLASLabelNew(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.6*offsetUp1)/newHeight), fAtlasLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(       0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.3*offsetUp1)/newHeight), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    gPad->RedrawAxis();

    if(flagGammas && flagSysts){
      for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/Ranking"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else if(flagGammas){
      for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/RankingGammas"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else if(flagSysts){
      for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/RankingSysts"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    else{
      std::cout << "WARNING::TtHFitter:Your ranking plot felt in unknown category :s" << endl;
      for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/RankingUnknown"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );
    }
    //
    delete c;
}

//____________________________________________________________________________________
//
void TtHFit::PrintSystTables(string opt){
    std::cout << "TtHFit::INFO: Printing syt tables" << std::endl;
    if(fCleanTables) opt += "clean";
    if(fSystCategoryTables) opt += "category";
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        fRegions[i_reg]->PrintSystTable(fFitResults,opt);
    }
}


//____________________________________________________________________________________
// this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter
void TtHFit::MergeSystematics(){
    // loop on systematics, see if any of them has name != nuisance
    for(auto syst : fSystematics){
        if(syst->fName!=syst->fNuisanceParameter){
            // if so, loop on other systematics to find one with name = nuisance parameter
            for(auto syst1 : fSystematics){
                if(syst->fName==syst1->fName) continue;
                if(!(syst->fNuisanceParameter==syst1->fName && syst1->fNuisanceParameter==syst1->fName)) continue;
                // now merge all SystematicHist in all regions
                for(auto reg : fRegions){
                    for(auto sh : reg->fSampleHists){
                        SystematicHist *syh  = sh->GetSystematic(syst ->fName);
                        SystematicHist *syh1 = sh->GetSystematic(syst1->fName);
                        if(syh!=0x0 && syh1!=0x0){
                            // FIXME...
                            // the issue here is that to combine uncertainties one has to act differently depending on the fact that the different sources come from a multiplication/division or not...
//                             if(sh->fSample->fMultiplyBy!="" && sh->fSample->fDivideBy!=""){
//                                 syh1 ->Multiply(syh);
//                                 syh1 ->Divide(sh->fHist);
//                             }
//                             else{
                                syh1 ->Add(syh);
                                syh1 ->Add(sh->fHist,-1);
//                             }
//                             //
//                             // up variation
//                             TH1* htmpUp   = (TH1*)syh->fHistUp->Clone();
//                             htmpUp->Add(sh->fHist,-1);
//                             htmpUp->Divide(sh->fHist);
//                             TH1* htmpUp1   = (TH1*)syh1->fHistUp->Clone();
//                             htmpUp1->Add(sh->fHist,-1);
//                             htmpUp1->Divide(sh->fHist);
//                             htmpUp->Add(htmpUp1);
//                             htmpUp->Multiply(sh->fHist);
//                             htmpUp->Add(sh->fHist);
//                             syh1->fHistUp = (TH1*)htmpUp->Clone();
//                             delete htmpUp;
//                             // down variation
//                             TH1* htmpDown   = (TH1*)syh->fHistDown->Clone();
//                             htmpDown->Add(sh->fHist,-1);
//                             htmpDown->Divide(sh->fHist);
//                             TH1* htmpDown1   = (TH1*)syh1->fHistDown->Clone();
//                             htmpDown1->Add(sh->fHist,-1);
//                             htmpDown1->Divide(sh->fHist);
//                             htmpDown->Add(htmpDown1);
//                             htmpDown->Multiply(sh->fHist);
//                             htmpDown->Add(sh->fHist);
//                             syh1->fHistDown = (TH1*)htmpDown->Clone();
//                             delete htmpDown;
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
void TtHFit::ComputeBining(int regIter){
    //
    //Creating histograms to rebin
    TH1F* hsig = 0x0;
    TH1F* hbkg = 0x0;
    TH1F* htmp = 0x0;
    bool nDefSig=true;
    bool nDefBkg=true;
    string fullSelection;
    string fullMCweight;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();
    bool bkgReg=false;
    bool flatBkg=false;
    if(fRegions[regIter]->fRegionType==Region::CONTROL) bkgReg=true;
    if(bkgReg && fRegions[regIter]->fTransfoDzSig<1e-3) flatBkg=true;
    //
    if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit DEBUG::ComputeBinning: Will compute binning with the following options: ";

    if(TtHFitter::DEBUGLEVEL>0 && (fRegions[regIter]->fTransfoDzSig>1e-3 || fRegions[regIter]->fTransfoDzBkg>1e-3) )
      std::cout << " TransfoD - zSig=" << fRegions[regIter]->fTransfoDzSig << " - zBkg=" << fRegions[regIter]->fTransfoDzBkg;
    if(TtHFitter::DEBUGLEVEL>0 && (fRegions[regIter]->fTransfoFzSig>1e-3 || fRegions[regIter]->fTransfoFzBkg>1e-3) )
      std::cout << " TransfoF - zSig=" << fRegions[regIter]->fTransfoFzSig << " - zBkg=" << fRegions[regIter]->fTransfoFzBkg;
    if(TtHFitter::DEBUGLEVEL>0 && (fRegions[regIter]->fTransfoJpar1>1e-3 || fRegions[regIter]->fTransfoJpar2>1e-3 || fRegions[regIter]->fTransfoJpar3>1e-3) )
      std::cout << " TransfoJ - z1=" << fRegions[regIter]->fTransfoJpar1 << " - z2=" << fRegions[regIter]->fTransfoJpar2 << " - z3=" << fRegions[regIter]->fTransfoJpar3;

    if(TtHFitter::DEBUGLEVEL>0 && bkgReg) std::cout << " - bkg reg";
    else if(TtHFitter::DEBUGLEVEL>0) std::cout << " - sig reg";

    if(TtHFitter::DEBUGLEVEL>0) std::cout << std::endl;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        //
        // using NTuples
        if(fInputType==1){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            if(fSamples[i_smp]->fType==Sample::GHOST) continue;
            if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[regIter]->fName)<0 ) continue;
            fullSelection = "1";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSelection!="" && fSelection!="1")
                fullSelection += " && ("+fSelection+")";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fRegions[regIter]->fSelection!="" && fRegions[regIter]->fSelection!="1")
                fullSelection += " && ("+fRegions[regIter]->fSelection+")";
            if(fSamples[i_smp]->fSelection!="" && fSamples[i_smp]->fSelection!="1")
                fullSelection += " && ("+fSamples[i_smp]->fSelection+")";
            if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSamples[i_smp]->fIgnoreSelection!="FALSE" && fSamples[i_smp]->fIgnoreSelection!="")
                fullSelection = ReplaceString(fullSelection,fSamples[i_smp]->fIgnoreSelection,"1");
            if(!fSamples[i_smp]->fNormalizedByTheory){
                fullMCweight = fSamples[i_smp]->fMCweight;
            }
            else{
                fullMCweight = fMCweight + " * " + fSamples[i_smp]->fMCweight;
                if(fRegions[regIter]->fMCweight!="") fullMCweight += " * " + fRegions[regIter]->fMCweight;
            }
            //
            // build a list of ntuples to read
            fullPaths.clear();
            vector<string> NtupleNames;
            for(unsigned int ns_ch=0; ns_ch<fRegions[regIter]->fNtupleNames.size(); ++ns_ch){
                NtupleNames.push_back(fRegions[regIter]->fNtupleNames.at(ns_ch));
            }
            for(unsigned int ns_smp=0; ns_smp<fSamples[i_smp]->fNtupleNames.size(); ++ns_smp){
                NtupleNames.push_back(fSamples[i_smp]->fNtupleNames.at(ns_smp));
            }
            vector<string> NtupleNameSuffs = CombinePathSufs( fSamples[i_smp]->fNtupleNameSuffs,
                                                              fRegions[regIter]->fNtupleNameSuffs );
            fullPaths = CreatePathsList( fSamples[i_smp]->fNtuplePaths.size()>0 ? fSamples[i_smp]->fNtuplePaths : fNtuplePaths,
                                        fRegions[regIter]->fNtuplePathSuffs,
                                        fSamples[i_smp]->fNtupleFiles.size()>0 ? fSamples[i_smp]->fNtupleFiles : ToVec(fNtupleFile), empty,
                                        NtupleNames.size()>0 ? NtupleNames : ToVec( fNtupleName ),
                                        NtupleNameSuffs.size()>0 ? NtupleNameSuffs : empty  // NEW
                                        );
              for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                  int tmp_debugLevel=TtHFitter::DEBUGLEVEL;
                  TtHFitter::SetDebugLevel(0);
                  htmp = HistFromNtuple( fullPaths[i_path],
                                        fRegions[regIter]->fVariable, 10000, fRegions[regIter]->fXmin, fRegions[regIter]->fXmax,
                                        fullSelection, fullMCweight);
                  TtHFitter::SetDebugLevel(tmp_debugLevel);
                  //
                  // Pre-processing of histograms (rebinning, lumi scaling)
                  if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                  //
                  if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                  else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                  //
                  // Importing the histogram in TtHFitter
                  if(fSamples[i_smp]->fType==Sample::SIGNAL){
                      if(nDefSig){
                          hsig = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
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
                            if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit DEBUG::ComputeBinning: Using "<<fSamples[i_smp]->fName<<" as signal" << std::endl;
                            if(nDefSig){
                                hsig = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                                nDefSig=false;
                            }
                            else hsig->Add(htmp);
                        }
                        else{
                            if(nDefBkg){
                              hbkg = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                              nDefBkg=false;
                            }
                            else hbkg->Add(htmp);
                        }
                    }
                    else{
                        if(nDefBkg){
                            hbkg = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp);
                    }
                }
                htmp->~TH1F();
            }
        }
        //
        // Input with hists
        else if(fInputType == 0){
            if(fSamples[i_smp]->fType==Sample::DATA) continue;
            if(fSamples[i_smp]->fType==Sample::GHOST) continue;
            //
            // build a list of histograms to read
            fullPaths.clear();
            std::vector<string> histoFiles;
            std::vector<string> histoNames;
            if(fSamples[i_smp]->fHistoFiles.size()>0)     histoFiles = fSamples[i_smp]->fHistoFiles;
            else if(fRegions[regIter]->fHistoFiles.size()>0) histoFiles = fRegions[regIter]->fHistoFiles;
            else                                          histoFiles = ToVec( fHistoFile );
            if(fSamples[i_smp]->fHistoNames.size()>0)     histoNames = fSamples[i_smp]->fHistoNames;
            else if(fRegions[regIter]->fHistoNames.size()>0) histoNames = fRegions[regIter]->fHistoNames;
            else                                          histoNames = ToVec( fHistoName );
            fullPaths = CreatePathsList( fHistoPaths, CombinePathSufs(fRegions[regIter]->fHistoPathSuffs, fSamples[i_smp]->fHistoPaths),
                        histoFiles, empty, // no histo file suffs for nominal (syst only)
                        histoNames, empty  // same for histo name
                        );

            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                int tmp_debugLevel=TtHFitter::DEBUGLEVEL;
                TtHFitter::SetDebugLevel(0);
                htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                TtHFitter::SetDebugLevel(tmp_debugLevel);
                //
                // Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[regIter]->fHistoBins){
                    TH1F* htmp2 = (TH1F*)(htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin,"htmp2",fRegions[regIter]->fHistoBins));
                    const char *hname = htmp->GetName();
                    htmp->~TH1F();
                    htmp = (TH1F*)htmp2->Clone();
                    delete htmp2;
                    htmp->SetName(hname);
                }
                else if(fRegions[regIter]->fHistoNBinsRebin != -1) {
                    htmp = (TH1F*)(htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin));
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
                    hsig = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
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
                    if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit DEBUG::ComputeBinning: Using "<<fSamples[i_smp]->fName<<" as signal" << std::endl;
                    if(nDefSig){
                        hsig = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                        nDefSig=false;
                    }
                    else hsig->Add(htmp);
                }
                else{
                    if(nDefBkg){
                        hbkg = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                        nDefBkg=false;
                    }
                    else hbkg->Add(htmp);
                }
            }
            else{
                if(nDefBkg){
                    hbkg = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                    nDefBkg=false;
                }
                else hbkg->Add(htmp);
            }
          }
          //
          htmp->~TH1F();
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
    vector<int> bins_vec;
    bins_vec.clear();
    //
    if (!hbkg || !hsig) {
        std::cout << "ERROR: please provide signal and background histograms!" << std::endl;
        gSystem -> Exit(1);
    }
    int iBin2 = 0;
    int nBins_Vec = hbkg -> GetNbinsX();
    int iBin = nBins_Vec; // skip overflow bin
    bins_vec.push_back(nBins_Vec + 1);
    float nBkg = hbkg -> Integral(1, nBins_Vec );
    float nSig = hsig -> Integral(1, nBins_Vec );
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
            float err2Rel = 1;
            if(fRegions[regIter]->fBinTransfo == "TransfoD"){
                // "trafo D"
                if (sumBkg != 0 && sumSig != 0)
                  err2Rel = 1 / (sumBkg / (nBkg / fRegions[regIter]->fTransfoDzBkg) + sumSig / (nSig / fRegions[regIter]->fTransfoDzSig));
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
                std::cout << "ERROR: transformation method '" << fRegions[regIter]->fBinTransfo << "' unknown, try again!" << std::endl;
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
    std::cout << "TtHFit::INFO: Your final binning from automatic binning function is:"  << std::endl;
    std::cout << "TtHFit::INFO: -->> - ";
    for(unsigned int i_bins=0; i_bins<bins_vec.size(); ++i_bins){
      std::cout << bins[i_bins] << " - ";
    }
    std::cout << std::endl;
    //
    delete hsig;
    delete hbkg;
    delete htmp;
    fRegions[regIter]->SetBinning(nBins-1, bins);
    delete[] bins;
}

//__________________________________________________________________________________
//
void TtHFit::GetLikelihoodScan( RooWorkspace *ws, string varName, RooDataSet* data){
  std::cout << "TtHFit::INFO: Running likelihood scan for the parameter = " << varName << std::endl;

  RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());

  bool isPoI = false;
  RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  TString firstPOIname = (TString)firstPOI->GetName();
  if (firstPOIname.Contains(varName.c_str())) isPoI = true;

  RooRealVar* var = NULL;
  TString vname = "";
  bool foundSyst = false;
  Double_t minVal = -3;
  Double_t maxVal =  3;

  if (isPoI){
    TIterator* it = mc->GetParametersOfInterest()->createIterator();
    while( (var = (RooRealVar*) it->Next()) ){
      vname=var->GetName();
      if (vname.Contains(varName.c_str())) { std::cout << "TtHFit::INFO: GetLikelihoodScan for POI = " << vname << std::endl; foundSyst=true; break; }
    }
  }
  else {
    TIterator* it = mc->GetNuisanceParameters()->createIterator();
    while( (var = (RooRealVar*) it->Next()) ){
      vname=var->GetName();
      if (vname.Contains(varName.c_str())) { std::cout << "TtHFit::INFO: GetLikelihoodScan for NP = " << vname << std::endl; foundSyst=true; break; }
    }
  }

  if(!foundSyst){
    std::cout << "TtHFit::WARNING: GetLikelihoodScan : systematic "<< varName <<" not found (most probably due to Pruning), skip LHscan !"<<std::endl;
    return;
  }
  std::cout << "TtHFit::INFO: GetLikelihoodScan for parameter = " << vname << std::endl;

  TF1* poly = new TF1("poly2","[0]+[1]*x+[2]*x*x",0,10);
  TCanvas* can = new TCanvas("NLLscan");

  RooAbsReal* nll = simPdf->createNLL(*data, Constrain(*mc->GetNuisanceParameters()), Offset(1), NumCPU(TtHFitter::NCPU, RooFit::Hybrid));

  TString tag("");
  RooAbsReal* pll = nll->createProfile(*var);
//   RooPlot* frameLH = var->frame(Title("-log(L) vs "+vname),Bins(30),Range(-1.5,3.5));
  RooPlot* frameLH = var->frame(Title("-log(L) vs "+vname),Bins(30),Range(minVal, maxVal));
  pll->plotOn(frameLH,RooFit::Precision(-1),LineColor(kRed), NumCPU(TtHFitter::NCPU));
  RooCurve* curve = frameLH->getCurve();
  curve->Draw();

  float val = var->getVal();
//   Double_t minVal = 0;
//   Double_t maxVal = 2;
//   if(val>1/5) {
//     minVal = val - 2; if(minVal<0) { minVal = 0; }
//     maxVal = val + 2;
//   }
  //Double_t minVal = -3.;
  //Double_t maxVal =  3.;
  frameLH->GetXaxis()->SetRangeUser(minVal,maxVal);

  // fit function
  for( int par(0); par<3; par++) { poly->SetParameter(par,0); }
  curve->Fit(poly,"RQN"); // R=range, Q=quiet, N=do not draw
  TString fitStr = Form("%5.2f + %5.2fx + %5.2fx^{2}", poly->GetParameter(0), poly->GetParameter(1), poly->GetParameter(2));
  TLatex* latex = new TLatex();
  latex->SetNDC(); // latex->SetTextSize(0.055);
//   latex->SetTextSize(gStyle->GetTextSize());
  latex->SetTextAlign(32);
  latex->SetText(0.925,0.925, fitStr);

  // y axis
  //frameLH->updateYAxis(minVal,maxVal,"");
//   frameLH->GetYaxis()->SetRangeUser(minVal,5.0);
  //frameLH->GetYaxis()->SetRangeUser(0,5.0);
  frameLH->GetYaxis()->SetTitle("#Delta [-Log(L)]");
  if(TtHFitter::NPMAP[varName]!="") frameLH->GetXaxis()->SetTitle(TtHFitter::NPMAP[varName].c_str());

  TString cname="";
  cname.Append("NLLscan_");
  cname.Append(vname);

  can->SetTitle(cname);
  can->SetName(cname);
  can->cd();
  frameLH->Draw();
//   latex->Draw("same");

  TLine *l1s = new TLine(minVal,0.5,maxVal,0.5);
  l1s->SetLineStyle(kDotted);
  l1s->Draw();

  TString LHDir("LHoodPlots/");
  system(TString("mkdir -vp ")+fName+"/"+LHDir);

  for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
      can->SaveAs( fName+"/"+LHDir+"NLLscan_"+varName+"."+TtHFitter::IMAGEFORMAT[i_format] );
}

//____________________________________________________________________________________
//
void TtHFit::defineVariable(int regIter){
    TH1::StatOverflows(true);  //////  What is the defaut in root for this ???
    if(TtHFitter::DEBUGLEVEL>0) std::cout << "//////// --------" << std::endl
                                          << "// DEBUG CORR VAR" << std::endl;
    TH1* h1 = new TH1D("h1","h1",1,-2000.,1000.);
    TH1* h2 = new TH1D("h2","h2",1,-2000.,1000.);
    string fullSelection;
    string fullMCweight;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();

    // copy of NtupleReading function.
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(TtHFitter::DEBUGLEVEL>0) std::cout<<"Processing sample : "<<fSamples[i_smp]->fName<<std::endl;
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        if( FindInStringVector(fSamples[i_smp]->fRegions,fRegions[regIter]->fName)<0 ) continue;
        if(TtHFitter::DEBUGLEVEL>0) std::cout<<" -> is used in the considered region"<<std::endl;
        fullSelection = "1";
        if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSelection!="" && fSelection!="1")
            fullSelection += " && ("+fSelection+")";
        if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fRegions[regIter]->fSelection!="" && fRegions[regIter]->fSelection!="1")
            fullSelection += " && ("+fRegions[regIter]->fSelection+")";
        if(fSamples[i_smp]->fSelection!="" && fSamples[i_smp]->fSelection!="1")
            fullSelection += " && ("+fSamples[i_smp]->fSelection+")";
        if(fSamples[i_smp]->fIgnoreSelection!="TRUE" && fSamples[i_smp]->fIgnoreSelection!="FALSE" && fSamples[i_smp]->fIgnoreSelection!="")
            fullSelection = ReplaceString(fullSelection,fSamples[i_smp]->fIgnoreSelection,"1");
        //
        if(!fSamples[i_smp]->fNormalizedByTheory){ // for data-driven bkg, use just the sample weight (FIXME)
            fullMCweight = fSamples[i_smp]->fMCweight;
        }
        else{
            fullMCweight = fMCweight + " * " + fSamples[i_smp]->fMCweight;
            if(fRegions[regIter]->fMCweight!="") fullMCweight += " * " + fRegions[regIter]->fMCweight;
        }
        fullPaths.clear();
        vector<string> NtupleNames;
        for(unsigned int ns_ch=0; ns_ch<fRegions[regIter]->fNtupleNames.size(); ++ns_ch){
            NtupleNames.push_back(fRegions[regIter]->fNtupleNames.at(ns_ch));
        }
        for(unsigned int ns_smp=0; ns_smp<fSamples[i_smp]->fNtupleNames.size(); ++ns_smp){
            NtupleNames.push_back(fSamples[i_smp]->fNtupleNames.at(ns_smp));
        }
        vector<string> NtupleNameSuffs = CombinePathSufs( fSamples[i_smp]->fNtupleNameSuffs,
                                                          fRegions[regIter]->fNtupleNameSuffs );
        fullPaths = CreatePathsList( fSamples[i_smp]->fNtuplePaths.size()>0 ? fSamples[i_smp]->fNtuplePaths : fNtuplePaths,
                                    fRegions[regIter]->fNtuplePathSuffs,
                                    fSamples[i_smp]->fNtupleFiles.size()>0 ? fSamples[i_smp]->fNtupleFiles : ToVec(fNtupleFile), empty, // no ntuple file suffs for nominal (syst only)
                                    NtupleNames.size()>0 ? NtupleNames : ToVec( fNtupleName ),
                                    NtupleNameSuffs.size()>0 ? NtupleNameSuffs : empty  // NEW
                                    );
        //
        for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
            if(TtHFitter::DEBUGLEVEL>0) std::cout <<" -> Retriving : "<<Form("%s",fRegions[regIter]->fCorrVar1.c_str())
                                                  <<" w/ weight "<<Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str())
                                                  <<" from "<< fullPaths[i_path].c_str() << std::endl;
            TH1* htmp1 = new TH1D("htmp1","htmp1",1,-2000.,1000.);
            TH1* htmp2 = new TH1D("htmp2","htmp2",1,-2000.,1000.);
            TChain *t = new TChain();
            t->Add(fullPaths[i_path].c_str());
            t->Draw( Form("%s>>htmp1",fRegions[regIter]->fCorrVar1.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            t->Draw( Form("%s>>htmp2",fRegions[regIter]->fCorrVar2.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
//             t->~TChain();
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
//             htmp1->~TH1();
//             htmp2->~TH1();
            delete htmp1;
            delete htmp2;
        }
    }

    double mean1 = h1->GetMean();
    double rms1 = h1->GetRMS();
    double mean2 = h2->GetMean();
    double rms2 = h2->GetRMS();

    if(TtHFitter::DEBUGLEVEL>0) std::cout<<"the new variable : "<<Form("( ( (%s)-%f )*( (%s)-%f ) )/( %f * %f )",fRegions[regIter]->fCorrVar1.c_str(),mean1,fRegions[regIter]->fCorrVar2.c_str(),mean2,rms1,rms2)<<std::endl;

    fRegions[regIter]->fVariable = Form("( ( (%s)-%f )*( (%s)-%f ) )/( %f * %f )",fRegions[regIter]->fCorrVar1.c_str(),mean1,fRegions[regIter]->fCorrVar2.c_str(),mean2,rms1,rms2);

    TH1::StatOverflows(false);  //////  What is the defaut in root for this ???

    delete h1;
    delete h2;
}
