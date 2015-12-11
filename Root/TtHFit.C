#include <cctype>

//TtHFitter headers
#include "TtHFitter/FittingTool.h"
#include "TtHFitter/HistoTools.h"

//Roofit headers
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"

//HistFactory headers
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"

//Corresponding header
#include "TtHFitter/TtHFit.h"

using namespace RooFit;

// -------------------------------------------------------------------------------------------------
// class TtHFit

//__________________________________________________________________________________
//
TtHFit::TtHFit(string name){
    fName = name;
    fLabel = "";
    fCmeLabel = "8 TeV";
    fLumiLabel = "20.3 fb^{-1}";
    
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
    
    fUpdate = false;
    
    fBlindingThreshold = -1;
    
    fRankingMaxNP = 10;
    fRankingOnly = "all";
    
    fStatOnly = false;
    
    //
    // Fit caracteristics
    //
    fFitType = SPLUSB;
    fFitRegion = CRSR;
    fFitRegionsToFit.clear();
    fFitNPValues.clear();
    fFitPOIAsimov = 0;
    fFitIsBlind = false;
    
    //
    // Limit type
    //
    fLimitType = ASYMPTOTIC;
    fLimitIsBlind = false;
    fLimitPOIAsimov = 0;
    
    fImageFormat = "png";
    TtHFitter::IMAGEFORMAT.clear();
    TtHFitter::IMAGEFORMAT.push_back("png");
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
    fRegions[fNRegions]->fFitLabel = fLabel;
    fRegions[fNRegions]->fFitType = fFitType;
    fRegions[fNRegions]->fPOI = fPOI;
    fRegions[fNRegions]->fIntCode_overall = fIntCode_overall;
    fRegions[fNRegions]->fIntCode_shape   = fIntCode_shape;
    fRegions[fNRegions]->fLumiScale = fLumiScale;
    fRegions[fNRegions]->fBlindingThreshold = fBlindingThreshold;
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
    cout << "-------------------------------------------" << endl;
    cout << "Smoothing and/or Symmetrising Systematic Variations ..." << endl;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
            fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(syst);
        }
    }
}

//__________________________________________________________________________________
// create new root file with all the histograms
void TtHFit::WriteHistos(/*string fileName*/){
    bool recreate = !fUpdate;
    gSystem->mkdir( fName.c_str() );
    gSystem->mkdir( (fName + "/Histograms/").c_str() );
    string fileName = "";
    TDirectory *dir = gDirectory;
    TFile *f;
    bool singleOutputFile = !TtHFitter::SPLITHISTOFILES;
    //
    if(singleOutputFile){
//         fileName = fName + "/Histograms/" + fName + "_histos"+fSaveSuf+".root";
        fileName = fName + "/Histograms/" + fName + "_histos"+fSuffix+".root";
        cout << "-------------------------------------------" << endl;
        cout << "Writing histograms to file " << fileName << " ..." << endl;
        if(recreate){
            f = new TFile(fileName.c_str(),"RECREATE");
            f->~TFile();
            dir->cd();
        }
    }
    //
    SampleHist* h;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        //
        if(!singleOutputFile){
//             fileName = fName + "/Histograms/" + fName + "_" + fRegions[i_ch]->fName + "_histos"+fSaveSuf+".root";
            fileName = fName + "/Histograms/" + fName + "_" + fRegions[i_ch]->fName + "_histos"+fSuffix+".root";
            cout << "-------------------------------------------" << endl;
            cout << "Writing histograms to file " << fileName << " ..." << endl;
            if(recreate){
                f = new TFile(fileName.c_str(),"RECREATE");
                f->~TFile();
                dir->cd();
            }
        }
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if(h == 0x0){
                if(TtHFitter::DEBUGLEVEL>0) cout << "SampleHist[" << i_smp << "] not there." << endl;
                continue;
            }
            // set file and histo names for nominal
            h->fHistoName = h->fHist->GetName();
            h->fFileName = fileName;
            // set file and histo names for systematics
            for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
                h->fSyst[i_syst]->fFileNameUp = fileName;
                h->fSyst[i_syst]->fHistoNameUp = h->fSyst[i_syst]->fHistUp->GetName();
                h->fSyst[i_syst]->fFileNameDown = fileName;
                h->fSyst[i_syst]->fHistoNameDown = h->fSyst[i_syst]->fHistDown->GetName();
                if(h->fSyst[i_syst]->fIsShape){
                    h->fSyst[i_syst]->fFileNameShapeUp = fileName;
                    h->fSyst[i_syst]->fHistoNameShapeUp = h->fSyst[i_syst]->fHistShapeUp->GetName();
                    h->fSyst[i_syst]->fFileNameShapeDown = fileName;
                    h->fSyst[i_syst]->fHistoNameShapeDown = h->fSyst[i_syst]->fHistShapeDown->GetName();
                }
            }
            h->WriteToFile();
        }
    }
    cout << "-------------------------------------------" << endl;
}

//__________________________________________________________________________________
// Draw syst plots
void TtHFit::DrawSystPlots(){
    SampleHist* h;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
//             h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            h = fRegions[i_ch]->fSampleHists[i_smp];
            h->DrawSystPlot("all");
        }
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
//         optVec = Vectorize(options,';');
//         optVec = Vectorize(options,'-');
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
        if(optMap["Update"]!="" && optMap["Update"]!="FALSE")
            fUpdate = true;
        if(optMap["Ranking"]!="")
            fRankingOnly = optMap["Ranking"];
        if(optMap["Signal"]!="")
            onlySignal = optMap["Signal"];
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
    param = cs->Get("LumiScale");         if( param != "" ) fLumiScale = atof(param.c_str());
    param = cs->Get("SystPruningShape");  if( param != "")  fThresholdSystPruning_Shape         = atof(param.c_str());
    param = cs->Get("SystPruningNorm");   if( param != "")  fThresholdSystPruning_Normalisation = atof(param.c_str());
    param = cs->Get("IntCodeOverall");    if( param != "")  fIntCode_overall  = atoi(param.c_str());
    param = cs->Get("IntCodeShape");      if( param != "")  fIntCode_shape    = atoi(param.c_str());
    param = cs->Get("MCstatThreshold");   if( param != "")  SetStatErrorConfig( true,  atof(param.c_str()) );
                                          else              SetStatErrorConfig( false, 0.0 );
    param = cs->Get("DebugLevel");        if( param != "")  TtHFitter::SetDebugLevel( atoi(param.c_str()) );
    param = cs->Get("PlotOptions");       if( param != ""){
        vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "YIELDS")!=vec.end() )   TtHFitter::SHOWYIELDS = true;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TtHFitter::SHOWNORMSIG = true;
        // ...
    }
    param = cs->Get("SystControlPlots");  if( param != ""){
        if( param == "true" || param == "True" ||  param == "TRUE" ){
            TtHFitter::SYSTCONTROLPLOTS = true;
        } else {
            TtHFitter::SYSTCONTROLPLOTS = false;
        }
    }
    param = cs->Get("CorrelationThreshold"); if( param != ""){
        TtHFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }
    param = cs->Get("SignalRegionsPlot");  if(param != ""){
        fRegionsToPlot = Vectorize(param,',');
    }
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
    param = cs->Get("ImageFormat");  if( param != ""){
        fImageFormat = Vectorize(param,',')[0];
        TtHFitter::IMAGEFORMAT = Vectorize(param,',');
    }
    
    //##########################################################
    //
    // FIT options
    //
    //##########################################################
    cs = fConfig->GetConfigSet("Fit");
    param = cs->Get("FitType");    if( param != "" ){
        if( param == "SPLUSB" )
            SetFitType(TtHFit::SPLUSB);
        else if( param == "BONLY" )
            SetFitType(TtHFit::BONLY);
        else{
            std::cerr << "Unknown FitType argument : " << cs->Get("FitType") << std::endl;
            return;
        }
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
    param = cs->Get("StatOnly");    if( param != "" ){
         std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fStatOnly = true;
        } else if ( param == "FALSE" ){
            fStatOnly = false;
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
        param = cs->Get("POIAsimov");  if( param != "" ){ fLimitPOIAsimov = atof(param.c_str()); };
    }

    //##########################################################
    //
    // REGIONS options
    //
    //##########################################################
    int nReg = 0;
    Region *reg;
    while(true){
        cs = fConfig->GetConfigSet("Region",nReg);
        if(cs==0x0) break;
        nReg++;
        if(onlyRegions.size()>0 && FindInStringVector(onlyRegions,cs->GetValue())<0) continue;
        if(toExclude.size()>0 && FindInStringVector(toExclude,cs->GetValue())>=0) continue;
        reg = NewRegion(CheckName(cs->GetValue()));
        reg->SetVariableTitle(cs->Get("VariableTitle"));
        reg->SetLabel(cs->Get("Label"),cs->Get("ShortLabel"));
        param = cs->Get("LumiLabel"); if( param != "") reg->fLumiLabel = param; 
        else reg->fLumiLabel = fLumiLabel;
        param = cs->Get("CmeLabel"); if( param != "") reg->fCmeLabel = param;
        else reg->fCmeLabel = fCmeLabel;
        param = cs->Get("LogScale"); if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="TRUE") reg->fLogScale = true;
            if(param=="FALSE") reg->fLogScale = false;
        }
        if(fInputType==0){
            param = cs->Get("HistoFile"); if(param!="") reg->fHistoFiles.push_back( param );
            param = cs->Get("HistoName"); if(param!="") reg->SetHistoName( param );
        }
        else if(fInputType==1){
            vector<string> variable = Vectorize(cs->Get("Variable"),',');
            reg->SetVariable(  variable[0], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()) );
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
//             reg->AddMCweight(  cs->Get("MCweight") );
            reg->fMCweight = cs->Get("MCweight"); // this will override the global MCweight, if any
            if(cs->Get("NtuplePathSuff")!="") { reg->fNtuplePathSuffs.clear(); reg->fNtuplePathSuffs.push_back( cs->Get("NtuplePathSuff") ); }
            param = cs->Get("NtuplePathSuffs");
            if( param != "" ){
                reg->fNtuplePathSuffs.clear();
                std::vector<string> paths = Vectorize( param,',' );
                for(int i=0;i<(int)paths.size();i++){
//                     reg->fNtuplePathSuffs.push_back( paths[i] );
                    reg->fNtuplePathSuffs.push_back( Fix(paths[i]) );
                }
            }
        }
        //Potential rebinning
        if(cs->Get("Rebin")!="") reg -> Rebin(atoi(cs->Get("Rebin").c_str()));
        if(cs->Get("Binning")!=""){
            std::vector < string > vec_bins = Vectorize(cs->Get("Binning"), ',');
            const int nBounds = vec_bins.size();
            double bins[nBounds];
            for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
                bins[iBound] = atof(vec_bins[iBound].c_str());
            }
            reg -> SetBinning(nBounds-1,bins);
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
    }
    
    //##########################################################
    //
    // SAMPLES options
    //
    //##########################################################
    int nSmp = 0;
    Sample *smp;
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
        param = cs->Get("Group"); if(param!="") smp->fGroup = param;
        if(fInputType==0){
            param = cs->Get("HistoFile"); if(param!="") smp->AddHistoFile( param );
            param = cs->Get("HistoName"); if(param!="") smp->fHistoNames.push_back( param );
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
        if(cs->Get("NormFactor")!="")
            smp->AddNormFactor(
                Vectorize(cs->Get("NormFactor"),',')[0],
                atof(Vectorize(cs->Get("NormFactor"),',')[1].c_str()),
                atof(Vectorize(cs->Get("NormFactor"),',')[2].c_str()),
                atof(Vectorize(cs->Get("NormFactor"),',')[3].c_str())
            );
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
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "TRUE") smp->fIgnoreSelection = true;
        }
//         param = cs->Get("IsGhost");
//         if(param!=""){
//             std::transform(param.begin(), param.end(), param.begin(), ::toupper);
//             if(param == "TRUE") smp->fIsGhost = true;
//         }
        // ...
    }
    
    //##########################################################
    //
    // SYSTEMATICS options
    //
    //##########################################################
    int nSys = 0;
    Systematic *sys;
    Sample *sam;


    //Addition for StatOnly fit: dummy systematic for the significance computation and limit setting
    int typed=0;
    Systematic *sysd; 
    if (fStatOnly) {
      typed = Systematic::OVERALL;
      sysd = new Systematic("Dummy",typed);
      sysd->fOverallUp   = 0.;
      sysd->fOverallDown = -0.;
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
        string exclude_str = cs->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        vector<string> samples = Vectorize(samples_str,',');
        vector<string> exclude = Vectorize(exclude_str,',');
        type = Systematic::HISTO;
        if(cs->Get("Type")=="overall" || cs->Get("Type")=="OVERALL")
            type = Systematic::OVERALL;
        sys = new Systematic(CheckName(cs->GetValue()),type);
        fSystematics.push_back( sys );
        fNSyst++;
        if(cs->Get("Title")!=""){
            sys->fTitle = cs->Get("Title");
            TtHFitter::SYSTMAP[sys->fName] = sys->fTitle;
        }
        if(cs->Get("Category")!=""){
            sys->fCategory = cs->Get("Category");
        }
        
        if(type==Systematic::HISTO){
            if(fInputType==0){
                if(cs->Get("HistoPathUp")!="")      sys->fHistoPathsUp  .push_back(cs->Get("HistoPathUp"));
                if(cs->Get("HistoPathDown")!="")    sys->fHistoPathsDown.push_back(cs->Get("HistoPathDown"));
                if(cs->Get("HistoPathSufUp")!="")   sys->fHistoPathSufUp   = cs->Get("HistoPathSufUp");
                if(cs->Get("HistoPathSufDown")!="") sys->fHistoPathSufDown = cs->Get("HistoPathSufDown");
                if(cs->Get("HistoFileUp")!="")      sys->fHistoFilesUp  .push_back(cs->Get("HistoFileUp"));
                if(cs->Get("HistoFileDown")!="")    sys->fHistoFilesDown.push_back(cs->Get("HistoFileDown"));
                if(cs->Get("HistoFileSufUp")!="")   sys->fHistoFileSufUp   = cs->Get("HistoFileSufUp");
                if(cs->Get("HistoFileSufDown")!="") sys->fHistoFileSufDown = cs->Get("HistoFileSufDown");
                if(cs->Get("HistoNameUp")!="")      sys->fHistoNamesUp  .push_back(cs->Get("HistoNameUp"));
                if(cs->Get("HistoNameDown")!="")    sys->fHistoNamesDown.push_back(cs->Get("HistoNameDown"));
                if(cs->Get("HistoNameSufUp")!="")   sys->fHistoNameSufUp   = cs->Get("HistoNameSufUp");
                if(cs->Get("HistoNameSufDown")!="") sys->fHistoNameSufDown = cs->Get("HistoNameSufDown");
                // ...
            }
            else if(fInputType==1){
                if(cs->Get("NtuplePathUp")!="")      sys->fNtuplePathsUp  .push_back(cs->Get("NtuplePathsUp"));
                if(cs->Get("NtuplePathDown")!="")    sys->fNtuplePathsDown.push_back( cs->Get("NtuplePathsDown"));
                if(cs->Get("NtuplePathSufUp")!="")   sys->fNtuplePathSufUp   = cs->Get("NtuplePathSufUp");
                if(cs->Get("NtuplePathSufDown")!="") sys->fNtuplePathSufDown = cs->Get("NtuplePathSufDown");
                if(cs->Get("NtupleFileUp")!="")      sys->fNtupleFilesUp  .push_back(cs->Get("NtupleFileUp"));
                if(cs->Get("NtupleFileDown")!="")    sys->fNtupleFilesDown.push_back( cs->Get("NtupleFileDown"));
                if(cs->Get("NtupleFileSufUp")!="")   sys->fNtupleFileSufUp   = cs->Get("NtupleFileSufUp");
                if(cs->Get("NtupleFileSufDown")!="") sys->fNtupleFileSufDown = cs->Get("NtupleFileSufDown");
                if(cs->Get("NtupleNameUp")!="")      sys->fNtupleNamesUp  .push_back(cs->Get("NtupleNameUp"));
                if(cs->Get("NtupleNameDown")!="")    sys->fNtupleNamesDown.push_back( cs->Get("NtupleNameDown"));
                if(cs->Get("NtupleNameSufUp")!="")   sys->fNtupleNameSufUp   = cs->Get("NtupleNameSufUp");
                if(cs->Get("NtupleNameSufDown")!="") sys->fNtupleNameSufDown = cs->Get("NtupleNameSufDown");
                if(cs->Get("WeightUp")!="")          sys->fWeightUp      = cs->Get("WeightUp");
                if(cs->Get("WeightDown")!="")        sys->fWeightDown    = cs->Get("WeightDown");
                if(cs->Get("WeightSufUp")!="")       sys->fWeightSufUp   = cs->Get("WeightSufUp");
                if(cs->Get("WeightSufDown")!="")     sys->fWeightSufDown = cs->Get("WeightSufDown");
                if(cs->Get("IgnoreWeight")!="")      sys->fIgnoreWeight  = cs->Get("IgnoreWeight");
                // ...
            }
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
            // ...
        }
        else if(type==Systematic::OVERALL){
            sys->fOverallUp   = atof( cs->Get("OverallUp").c_str() );
            sys->fOverallDown = atof( cs->Get("OverallDown").c_str() );
        }
        // this to obtain syst variation relatively to given sample
        param = cs->Get("ReferenceSample"); if(param!="") sys->fReferenceSample = param;
        // attach the syst to the proper samples
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sam = fSamples[i_smp];
            if(sam->fType == Sample::DATA) continue;
            if(sam->fType == Sample::GHOST) continue;
            if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
               && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
            ){
                sam->AddSystematic(sys);
            }
        }
        // ...
    }
}


//__________________________________________________________________________________
// for each region, add a SampleHist for each Sample in the Fit, reading from ntuples
void TtHFit::ReadNtuples(){
    cout << "-------------------------------------------" << endl;
    cout << "Reading ntuples..." << endl;
    TH1F* h = 0x0;
    TH1F* hUp = 0x0;
    TH1F* hDown = 0x0;
    TH1F* htmp = 0x0;
    //   string ntupleFullPath;
    string fullSelection;
    string fullMCweight;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();
    //
    // loop on regions and samples
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        cout << "  Region " << fRegions[i_ch]->fName << " ..." << endl;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
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
            if(!fSamples[i_smp]->fIgnoreSelection && fSelection!="" && fSelection!="1")
                fullSelection += " && "+fSelection;
            if(!fSamples[i_smp]->fIgnoreSelection && fRegions[i_ch]->fSelection!="" && fRegions[i_ch]->fSelection!="1")
                fullSelection += " && "+fRegions[i_ch]->fSelection;
            if(fSamples[i_smp]->fSelection!="" && fSamples[i_smp]->fSelection!="1")
                fullSelection += " && "+fSamples[i_smp]->fSelection;
            //
            if(fSamples[i_smp]->fType==Sample::DATA) fullMCweight = "1";
            else if(!fSamples[i_smp]->fNormalizedByTheory){ // for data-driven bkg, use just the sample weight (FIXME)
                fullMCweight = fSamples[i_smp]->fMCweight;
            }
            else{
                fullMCweight = fMCweight + " * " + fSamples[i_smp]->fMCweight;
                if(fRegions[i_ch]->fMCweight!="") fullMCweight += " * " + fRegions[i_ch]->fMCweight;
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
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                htmp = HistFromNtuple( fullPaths[i_path],
                                      fRegions[i_ch]->fVariable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                      fullSelection, fullMCweight);
                
                //Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[i_ch]->fHistoBins){
                    TH1F* htmp2 = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,"htmp2",fRegions[i_ch]->fHistoBins));
                    const char *hname = htmp->GetName();
                    htmp->~TH1F();
                    htmp = htmp2;
                    htmp->SetName(hname);
                }
                else if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                
                if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                
                //Importing the histogram in TtHFitter
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp], h );
            std::map < int, bool > applyCorrection; applyCorrection.clear();
            
            if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::SIGNAL){
                for(unsigned int iBin = 1; iBin <= fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetNbinsX(); ++iBin ){
                    double content = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetBinContent(iBin);
                    if( content<=0 ){
                        std::cout << "WARNING: Checking your nominal histogram for sample " << fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fName << ": negative/null bin " << iBin << " content ! Trying to fix it." << std::endl;
                        fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->SetBinContent(iBin,1e-06);
                        applyCorrection.insert( std::pair < int, bool > (iBin, true) );
                    } else {
                        applyCorrection.insert( std::pair < int, bool > (iBin, false) );
                    }
                }
            }
            
            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                // if not Overall only...
                if(fSamples[i_smp]->fSystematics[i_syst]->fType==Systematic::OVERALL)
                    continue;
                if(TtHFitter::DEBUGLEVEL>0) cout << "Adding syst " << fSamples[i_smp]->fSystematics[i_syst]->fName << endl;
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                Systematic *syst = smp->fSystematics[i_syst];
                if(syst->fReferenceSample!="") smp = GetSample(syst->fReferenceSample);
                //
                // set selection
                fullSelection = "1";
                if(!smp->fIgnoreSelection && fSelection!="" && fSelection!="1")
                    fullSelection += " && "+fSelection;
                if(!smp->fIgnoreSelection && reg->fSelection!="" && reg->fSelection!="1")
                    fullSelection += " && "+reg->fSelection;
                if(smp->fSelection!="" && smp->fSelection!="1")
                    fullSelection += " && "+smp->fSelection;
                
                //
                // Up
                //
                // Note: no need to change selection for systematics. If needed, can be done via weight (partially...)
                fullMCweight = fMCweight;
                if(syst->fWeightUp!="")
                    fullMCweight += " * "+syst->fWeightUp;
                else
                    fullMCweight += " * "+smp->fMCweight;
                if(syst->fWeightSufUp!="")
                    fullMCweight += " * "+syst->fWeightSufUp;
                if(reg->fMCweight!="")
                    fullMCweight += " * "+reg->fMCweight;
                if(syst->fIgnoreWeight!=""){
                    ReplaceString(fullMCweight, syst->fIgnoreWeight,"");
                    ReplaceString(fullMCweight,"*  *","*");
                    ReplaceString(fullMCweight,"* *","*");
                    ReplaceString(fullMCweight,"**","*");
                }
                vector<string> s = CombinePathSufs(
                                                   reg->fNtuplePathSuffs,
                                                   syst->fNtuplePathsUp );
                //
                fullPaths.clear();
                vector<string> NtupleNameSuffsUp = CombinePathSufs( ToVec( syst->fNtupleNameSufUp ),
                                                                   reg->fNtupleNameSuffs );
                fullPaths = CreatePathsList(
                                            // path
                                            smp->fNtuplePaths.size()>0 ? smp->fNtuplePaths : fNtuplePaths,
                                            // path suf
                                            CombinePathSufs(reg->fNtuplePathSuffs,syst->fNtuplePathsUp ),
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
                    htmp = HistFromNtuple( fullPaths[i_path],
                                          reg->fVariable, reg->fNbins, reg->fXmin, reg->fXmax,
                                          fullSelection, fullMCweight);
                    //Pre-processing of histograms (rebinning, lumi scaling)
                    if(reg->fHistoBins){
                        //                         htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,htmp->GetName(),reg->fHistoBins));
                        TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                        const char *hname = htmp->GetName();
                        htmp->~TH1F();
                        htmp = htmp2;
                        htmp->SetName(hname);
                    }
                    else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                    
                    if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                    
                    if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                    else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                    
//                     // obtain relative variation and apply it to proper sample
//                     if(syst->fReferenceSample!=""){
//                         TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
//                         htmp->Divide(href);
//                         htmp->Multiply( reg->GetSampleHist( fSamples[i_smp]->fName )->fHist );
//                     }
                    // obtain relative variation and apply it to proper sample
                    // & try to keep also the same total relative variation
                    if(syst->fReferenceSample!=""){
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                        TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                        float relVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                        htmp->Divide(   href );
                        htmp->Multiply( hnom );
                        float newVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                        if(relVar > 0.0001 && newVar > 0.0001) htmp->Scale( relVar / newVar );
                    }
                    
                    //Importing histogram in TtHFitter
                    if(i_path==0){
                        hUp = (TH1F*)htmp->Clone(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fName.c_str()));
                    }
                    else hUp->Add(htmp);
                    htmp->~TH1F();
                }
                //
                // Down
                //
                // Note: no need to change selection for systematics. If needed, can be done via weight (partially...)
                fullMCweight = fMCweight;
                if(syst->fWeightDown!="")
                    fullMCweight += " * "+syst->fWeightDown;
                else
                    fullMCweight += " * " +smp->fMCweight;
                if(syst->fWeightSufDown!="")
                    fullMCweight += " * "+syst->fWeightSufDown;
                if(reg->fMCweight!="")
                    fullMCweight += " * "+reg->fMCweight;
                if(syst->fIgnoreWeight!=""){
                    ReplaceString(fullMCweight, syst->fIgnoreWeight,"");
                    ReplaceString(fullMCweight,"*  *","*");
                    ReplaceString(fullMCweight,"* *","*");
                    ReplaceString(fullMCweight,"**","*");
                }
                
                //
                fullPaths.clear();
                vector<string> NtupleNameSuffsDown = CombinePathSufs( ToVec( syst->fNtupleNameSufDown ),
                                                                     reg->fNtupleNameSuffs );
                fullPaths = CreatePathsList(
                                            // path
                                            smp->fNtuplePaths.size()>0 ? smp->fNtuplePaths : fNtuplePaths,
                                            // path suf
                                            CombinePathSufs(
                                                            reg->fNtuplePathSuffs,
                                                            syst->fNtuplePathsDown ),
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
                    htmp = HistFromNtuple( fullPaths[i_path],
                                          reg->fVariable, reg->fNbins, reg->fXmin, reg->fXmax,
                                          fullSelection, fullMCweight);
                    //Pre-processing of histograms (rebinning, lumi scaling)
                    if(reg->fHistoBins){
                        //                         htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,htmp->GetName(),reg->fHistoBins));
                        TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                        const char *hname = htmp->GetName();
                        htmp->~TH1F();
                        htmp = htmp2;
                        htmp->SetName(hname);
                    }
                    else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                    
                    if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                    
                    if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                    else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                    
//                     // obtain relative variation and apply it to proper sample
//                     if(syst->fReferenceSample!=""){
//                         TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
//                         htmp->Divide(href);
//                         htmp->Multiply( reg->GetSampleHist( fSamples[i_smp]->fName )->fHist );
//                     }
                    // obtain relative variation and apply it to proper sample
                    // & try to keep also the same total relative variation
                    if(syst->fReferenceSample!=""){
                        TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                        TH1* hnom = reg->GetSampleHist( fSamples[i_smp]->fName )->fHist;
                        float relVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / href->Integral(0,href->GetNbinsX()+1);
                        htmp->Divide(   href );
                        htmp->Multiply( hnom );
                        float newVar   = htmp->Integral(0,htmp->GetNbinsX()+1) / hnom->Integral(0,hnom->GetNbinsX()+1);
                        if(relVar > 0.0001 && newVar > 0.0001) htmp->Scale( relVar / newVar );
                    }
                    
                    //Importing histogram in TtHFitter
                    if(i_path==0){
                        hDown = (TH1F*)htmp->Clone(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fName.c_str()));
                        //                         hDown->SetName(Form("h_%s_%s_%sDown",reg->fName.c_str(),smp->fName.c_str(),syst->fName.c_str()));
                    }
                    else hDown->Add(htmp);
                    htmp->~TH1F();
                }
                //
                // Histogram smoothing, Symmetrisation, Massaging...
                SystematicHist *sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                sh -> fSmoothType = fSamples[i_smp]->fSystematics[i_syst] -> fSmoothType;
                sh -> fSymmetrisationType = fSamples[i_smp]->fSystematics[i_syst] -> fSymmetrisationType;
                sh -> fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                
                //
                // Histograms checking
                //
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::SIGNAL){
                    for(unsigned int iBin = 1; iBin <= fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetNbinsX(); ++iBin ){
                        if( applyCorrection[iBin]){
                            sh -> fHistUp   -> SetBinContent(iBin,1e-06);
                            sh -> fHistDown -> SetBinContent(iBin,1e-06);
                        }
                    }
                }
                HistoTools::CheckHistograms( fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist /*nominal*/,
                                            sh /*systematic*/,
                                            fSamples[i_smp]->fType!=Sample::SIGNAL/*check bins with content=0*/,
                                            TtHFitter::HISTOCHECKCRASH /*cause crash if problem*/);
            }
        }
    }
    delete htmp;
}

//__________________________________________________________________________________
//
void TtHFit::ReadHistograms(){
    TH1F* h;
    TH1F* hUp;
    TH1F* hDown;
    TH1F* htmp;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();
    //
    // loop on regions and samples
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(TtHFitter::DEBUGLEVEL>0) cout << "Reading " << fSamples[i_smp]->fName << endl;
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

            fullPaths = CreatePathsList( fHistoPaths, fRegions[i_ch]->fHistoPathSuffs,
                                        histoFiles, empty, // no histo file suffs for nominal (syst only)
                                        histoNames, empty  // same for histo name
                                        );
            
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                
                htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                
                //Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[i_ch]->fHistoBins){
                    TH1F* htmp2 = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,"htmp2",fRegions[i_ch]->fHistoBins));
                    const char *hname = htmp->GetName();
                    htmp->~TH1F();
                    htmp = htmp2;
                    htmp->SetName(hname);
                } else if(fRegions[i_ch]->fHistoNBinsRebin != -1) {
                    htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                }
                
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                
                if(fSamples[i_smp]->fLumiScales.size()>i_path) htmp -> Scale(fSamples[i_smp]->fLumiScales[i_path]);
                else if(fSamples[i_smp]->fLumiScales.size()==1) htmp -> Scale(fSamples[i_smp]->fLumiScales[0]);
                
                //Importing the histogram in TtHFitter
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
                
            }
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp], h );
            
            std::map < int, bool > applyCorrection;
            if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::SIGNAL){
                for(unsigned int iBin = 1; iBin <= fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetNbinsX(); ++iBin ){
                    double content = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetBinContent(iBin);
                    if( content<=0 ){
                        std::cout << "WARNING: Checking your nominal histogram for sample " << fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fName << ": negative/null content ! Trying to fix it." << std::endl;
                        fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->SetBinContent(iBin,1e-06);
                        applyCorrection.insert( std::pair < int, bool > (iBin, true) );
                    } else {
                        applyCorrection.insert( std::pair < int, bool > (iBin, false) );
                    }
                }
            }
            
            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                
                if(TtHFitter::DEBUGLEVEL>0) cout << "Adding syst " << fSamples[i_smp]->fSystematics[i_syst]->fName << endl;
                //
                Region *reg = fRegions[i_ch];
                Sample *smp = fSamples[i_smp];
                Systematic *syst = smp->fSystematics[i_syst];
                if(syst->fReferenceSample!="") smp = GetSample(syst->fReferenceSample);

                //
                // Up
                //
                // For histo syst:
                if(syst->fType==Systematic::HISTO){
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
                      //Pre-processing of histograms (rebinning, lumi scaling)
                      if(reg->fHistoBins){
                          TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                          const char *hname = htmp->GetName();
                          htmp->~TH1F();
                          htmp = htmp2;
                          htmp->SetName(hname);
                      }
                      else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                      
                      if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                      
                      if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                      else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                      
                      // obtain relative variation and apply it to proper sample
                      if(syst->fReferenceSample!=""){
                          TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                          htmp->Divide(href);
                          htmp->Multiply( reg->GetSampleHist( fSamples[i_smp]->fName )->fHist );
                      }
                      
                      //Importing histogram in TtHFitter
                      if(i_path==0) hUp = (TH1F*)htmp->Clone();
                      else hUp->Add(htmp);
                      htmp->~TH1F();
                  }
                }
                // For Overall syst
                else{
                  hUp = (TH1F*)h->Clone();
                  hUp->Scale(1+syst->fOverallUp);
                }
                hUp->SetName(Form("h_%s_%s_%sUp",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fName.c_str()));
                //
                // Down
                //
                // For histo syst:
                if(syst->fType==Systematic::HISTO){
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
                      //Pre-processing of histograms (rebinning, lumi scaling)
                      if(reg->fHistoBins){
                          TH1F* htmp2 = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin,"htmp2",reg->fHistoBins));
                          const char *hname = htmp->GetName();
                          htmp->~TH1F();
                          htmp = htmp2;
                          htmp->SetName(hname);
                      }
                      else if(reg->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(reg->fHistoNBinsRebin));
                      
                      if(smp->fType!=Sample::DATA && smp->fNormalizedByTheory) htmp -> Scale(fLumi);
                      
                      if(smp->fLumiScales.size()>i_path) htmp -> Scale(smp->fLumiScales[i_path]);
                      else if(smp->fLumiScales.size()==1) htmp -> Scale(smp->fLumiScales[0]);
                      
                      // obtain relative variation and apply it to proper sample
                      if(syst->fReferenceSample!=""){
                          TH1* href = reg->GetSampleHist(syst->fReferenceSample)->fHist;
                          htmp->Divide(href);
                          htmp->Multiply( reg->GetSampleHist( fSamples[i_smp]->fName )->fHist );
                      }
                      
                      //Importing histogram in TtHFitter
                      if(i_path==0) hDown = (TH1F*)htmp->Clone();
                      else hDown->Add(htmp);
                      htmp->~TH1F();
                  }
                }
                // For Overall syst
                else{
                  hDown = (TH1F*)h->Clone();
                  hDown->Scale(1+syst->fOverallDown);
                }
                hDown->SetName(Form("h_%s_%s_%sDown",reg->fName.c_str(),fSamples[i_smp]->fName.c_str(),syst->fName.c_str()));
                
                // 
                // Histogram smoothing, Symmetrisation, Massaging...
                //
                SystematicHist *sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                sh -> fSmoothType = fSamples[i_smp]->fSystematics[i_syst] -> fSmoothType;
                sh -> fSymmetrisationType = fSamples[i_smp]->fSystematics[i_syst] -> fSymmetrisationType;
                sh -> fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                
                //
                // Histograms checking
                //
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::SIGNAL){
                    for(unsigned int iBin = 1; iBin <= fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetNbinsX(); ++iBin ){
                        if( applyCorrection[iBin]){
                            sh -> fHistUp   -> SetBinContent(iBin,1e-06);
                            sh -> fHistDown -> SetBinContent(iBin,1e-06);
                        }
                    }
                }
                HistoTools::CheckHistograms( fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist /*nominal*/,
                                            sh /*systematic*/,
                                            fSamples[i_smp]->fType!=Sample::SIGNAL/*check bins with content=0*/,
                                            TtHFitter::HISTOCHECKCRASH /*cause crash if problem*/);
            }
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
    string systName;
    //
    bool singleOutputFile = !TtHFitter::SPLITHISTOFILES;
    if(singleOutputFile){
        fileName = fName + "/Histograms/" + fName + "_histos.root";
        cout << "-----------------------------" << endl;
        cout << "Reading histograms from file " << fileName << " ..." << endl;
    }
    //
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        regionName = fRegions[i_ch]->fName;
        if(TtHFitter::DEBUGLEVEL>0) cout << "  Reading region " << regionName << endl;
        //
        if(!singleOutputFile){
            fileName = fName + "/Histograms/" + fName + "_" + regionName + "_histos.root";
            cout << "-----------------------------" << endl;
            cout << "Reading histograms from file " << fileName << " ..." << endl;
        }
        //
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            // 
            // eventually skip sample / region combination
            //
            if( FindInStringVector(fSamples[i_smp]->fRegions,regionName)<0 ) continue;
            //
            sampleName = fSamples[i_smp]->fName;
            if(TtHFitter::DEBUGLEVEL>0) cout << "    Reading sample " << sampleName << endl;
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp],regionName+"_"+sampleName,fileName);
            sh = fRegions[i_ch]->GetSampleHist(sampleName);
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                systName = fSamples[i_smp]->fSystematics[i_syst]->fName;
                if(TtHFitter::DEBUGLEVEL>0) cout << "      Reading syst " << systName << endl;
                // norm only
                if(fSamples[i_smp]->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    sh->AddOverallSyst(systName,fSamples[i_smp]->fSystematics[i_syst]->fOverallUp,fSamples[i_smp]->fSystematics[i_syst]->fOverallDown);
                }
                // histo syst
                else{
                    syh = sh->AddHistoSyst(systName,
                                           Form("%s_%s_%s_Up",regionName.c_str(),sampleName.c_str(),systName.c_str()),   fileName,
                                           Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systName.c_str()), fileName);
                    syh->fHistoNameShapeUp   = Form("%s_%s_%s_Shape_Up",regionName.c_str(),sampleName.c_str(),systName.c_str());
                    syh->fHistoNameShapeDown = Form("%s_%s_%s_Shape_Down",regionName.c_str(),sampleName.c_str(),systName.c_str());
                    syh->fFileNameShapeUp    = fileName;
                    syh->fFileNameShapeDown  = fileName;
                }
            }
        }
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
        ReadFitResults(fName+"/Fits/"+fName+fSuffix+".txt");
    }
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->fUseStatErr = fUseStatErr;
        if(isPostFit){
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
TthPlot* TtHFit::DrawSummary(string opt){
    cout << "-------------------------------------------" << endl;
    cout << "Building Summary Plot..." << endl;
    gSystem->mkdir(fName.c_str());
    bool isPostFit = opt.find("post")!=string::npos;
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
    std::vector<int> regionVec; regionVec.clear();
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        if(fRegions[i_ch]->fRegionType!=Region::VALIDATION){
            regionVec.push_back(i_ch);
        }
    }
    int Nbin = (int)regionVec.size();
    if(Nbin<=0) return 0x0;
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        SampleHist *sh = 0x0;
        name = fSamples[i_smp]->fName.c_str();
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
            if(TtHFitter::DEBUGLEVEL>0) cout << "Adding Signal: " << h_sig[Nsig]->GetTitle() << endl;
            h_sig[Nsig]->SetLineColor(lineColor);
            h_sig[Nsig]->SetFillColor(fillColor);
            h_sig[Nsig]->SetLineWidth(lineWidth);
            for(int i_bin=1;(int)i_bin<=regionVec.size();i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=0x0){
                    if(isPostFit)  h = sh->fHist_postFit;
                    else           h = sh->fHist;
                    //
                    if(!isPostFit){
                        // scale it according to NormFactors
                        for(unsigned int i_nf=0;i_nf<sh->fSample->fNormFactors.size();i_nf++){
                            h->Scale(sh->fSample->fNormFactors[i_nf]->fNominal);
                            if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit::INFO: Scaling " << sh->fSample->fName << " by " << sh->fSample->fNormFactors[i_nf]->fNominal << std::endl;
                        }
                    }
                    //
                    integral = h->IntegralAndError(0,h->GetNbinsX()+1,intErr);
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
            if(TtHFitter::DEBUGLEVEL>0) cout << "Adding Bkg:    " << h_bkg[Nbkg]->GetTitle() << endl;
            h_bkg[Nbkg]->SetLineColor(lineColor);
            h_bkg[Nbkg]->SetFillColor(fillColor);
            h_bkg[Nbkg]->SetLineWidth(lineWidth);
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                sh = fRegions[regionVec[i_bin-1]]->GetSampleHist( name );
                if(sh!=0x0){
                    if(isPostFit)  h = sh->fHist_postFit;
                    else           h = sh->fHist;
                    //
                    if(!isPostFit){
                        // scale it according to NormFactors
                        for(unsigned int i_nf=0;i_nf<sh->fSample->fNormFactors.size();i_nf++){
                            h->Scale(sh->fSample->fNormFactors[i_nf]->fNominal);
                            if(TtHFitter::DEBUGLEVEL>0) std::cout << "TtHFit::INFO: Scaling " << sh->fSample->fName << " by " << sh->fSample->fNormFactors[i_nf]->fNominal << std::endl;
                        }
                    }
                    //
                    integral = h->IntegralAndError(0,h->GetNbinsX()+1,intErr);
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
            if(TtHFitter::DEBUGLEVEL>0) cout << "Adding Data:   " << h_data->GetTitle() << endl;
            for(int i_bin=1;i_bin<=(int)regionVec.size();i_bin++){
                if(fRegions[regionVec[i_bin-1]]->fRegionDataType==Region::ASIMOVDATA)
                    h_data->SetBinContent( i_bin,0 );
                else
                    h_data->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fData->fHist->Integral() );
            }
        }
    } 
    //
    TthPlot *p = new TthPlot(fName+"_summary",900,700);
    p->fYmin = 1;
    p->SetXaxis("",false);
    p->AddLabel(fLabel);
    if(isPostFit) p->AddLabel("Post-Fit");
    else          p->AddLabel("Pre-Fit");
    p->fATLASlabel = "Internal";
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    if(fBlindingThreshold>=0) p->SetBinBlinding(true,fBlindingThreshold);
    //
    if(h_data) p->SetData(h_data, h_data->GetTitle());
//     if(h_sig) p->AddSignal(h_sig,h_sig->GetTitle());
//     if(TtHFitter::SHOWNORMSIG) p->AddNormSignal(h_sig,((string)h_sig->GetTitle())+"(norm)");
    for(int i=0;i<Nsig;i++){
        p->AddSignal(h_sig[i],h_sig[i]->GetTitle());
        if(TtHFitter::SHOWNORMSIG) p->AddNormSignal(h_sig[i],((string)h_sig[i]->GetTitle())+"(norm)");
    }
    for(int i=0;i<Nbkg;i++)
        p->AddBackground(h_bkg[i],h_bkg[i]->GetTitle());
    //
    // Build tot
    h_tot = new TH1F("h_Tot_summary","h_Tot_summary", Nbin,0,Nbin);
    
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot_postFit->Integral() );
        else          h_tot->SetBinContent( i_bin,fRegions[regionVec[i_bin-1]]->fTot->Integral() );
        h_tot->SetBinError( i_bin,0 );
    }

    //
    //   Build error band
    // build the vectors of variations
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    TH1* h_tmp_Up;
    TH1* h_tmp_Down;
    for(int i_syst=0;i_syst<(int)fRegions[0]->fSystNames.size();i_syst++){
        for(int i_bin=1;i_bin<=Nbin;i_bin++){
	  if(isPostFit){
                h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp_postFit[i_syst];
                h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown_postFit[i_syst];
            }
            else{
                h_tmp_Up   = fRegions[regionVec[i_bin-1]]->fTotUp[i_syst];
                h_tmp_Down = fRegions[regionVec[i_bin-1]]->fTotDown[i_syst];
            }
            if(i_bin==1){
                h_up.  push_back( new TH1F(Form("%s_TMP",h_tmp_Up->GetName()),  h_tmp_Up->GetTitle(),   Nbin,0,Nbin) );
                h_down.push_back( new TH1F(Form("%s_TMP",h_tmp_Down->GetName()),h_tmp_Down->GetTitle(), Nbin,0,Nbin) );
            }
            h_up[i_syst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral() );
            h_down[i_syst]->SetBinContent( i_bin,h_tmp_Down->Integral() );
        }
    }
    if(isPostFit)  g_err = BuildTotError( h_tot, h_up, h_down, fRegions[0]->fSystNames, fFitResults->fCorrMatrix );
    else           g_err = BuildTotError( h_tot, h_up, h_down, fRegions[0]->fSystNames );
    //
    p->SetTotBkg(h_tot);
    p->SetTotBkgAsym(g_err);
    
    //
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        p->SetBinLabel(i_bin,fRegions[regionVec[i_bin-1]]->fShortLabel.c_str());
    }
    p->Draw(opt);
    //
    for(int i_bin=1;i_bin<=Nbin;i_bin++){
        if(TtHFitter::DEBUGLEVEL>0) cout << i_bin << ":\t" << h_tot->GetBinContent(i_bin) << "\t+" << g_err->GetErrorYhigh(i_bin-1) << "\t-" << g_err->GetErrorYlow(i_bin-1) << endl;
    }
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++){
	if(isPostFit)  p->SaveAs((fName+"/Plots/Summary_postFit"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
	else           p->SaveAs((fName+"/Plots/Summary"        +fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
    }
    //
    return p;
}

//__________________________________________________________________________________
//
void TtHFit::BuildYieldTable(string opt){
    cout << "-------------------------------------------" << endl;
    cout << "Building Yields Table..." << endl;
    bool isPostFit = opt.find("post")!=string::npos;
    ofstream out;
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Tables").c_str());
    if(!isPostFit)  out.open((fName+"/Tables/Yields.txt").c_str());
    else            out.open((fName+"/Tables/Yields_postFit.txt").c_str());
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
    out << " |       | ";
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        out << fRegions[i_bin-1]->fLabel << " | ";
    }
    out << endl;
    //
    std::vector< string > titleVec;
    std::vector< int > idxVec;
    SampleHist *sh = 0x0;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
//         if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        name = fSamples[i_smp]->fName;
        title = fSamples[i_smp]->fTitle;
        //
        int idx = FindInStringVector(titleVec,title);
        if(idx>=0){
            idxVec.push_back(idx);
        }
        else{
            idxVec.push_back(i_smp);
            h_smp[idxVec[i_smp]] = new TH1F(("h_"+name).c_str(),title.c_str(), fNRegions,0,fNRegions);
        }
        for(int i_bin=1;i_bin<=fNRegions;i_bin++){
            sh = fRegions[i_bin-1]->GetSampleHist( name );
            if(sh!=0x0){
                if(isPostFit && fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fType!=Sample::GHOST)
                    h0 = sh->fHist_postFit;
                else
                    h0 = sh->fHist;
                float tmpErr = h_smp[idxVec[i_smp]]->GetBinError(i_bin); // Michele -> get the error before adding content to bin, to avoid ROOT automatically increasing it!
                h_smp[idxVec[i_smp]]->AddBinContent( i_bin,h0->IntegralAndError(0,h0->GetNbinsX()+1,intErr) );
                h_smp[idxVec[i_smp]]->SetBinError(   i_bin, sqrt( pow(tmpErr,2) + pow(intErr,2) ) );
            }
        }
        titleVec.push_back(title);
    }
    //
    // add tot uncertainty on each sample
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        name = fSamples[i_smp]->fName;
        if(idxVec[i_smp]!=i_smp) continue;
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        // build the vectors of variations
        std::vector< TH1* > h_up;   h_up.clear();
        std::vector< TH1* > h_down; h_down.clear();
        TH1* h_tmp_Up;
        TH1* h_tmp_Down;
        for(int i_syst=0;i_syst<(int)fRegions[0]->fSystNames.size();i_syst++){
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                sh = fRegions[i_bin-1]->GetSampleHist( name );
                if(sh!=0x0){
                    if(isPostFit){
                        h_tmp_Up   = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp_postFit;
                        h_tmp_Down = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown_postFit;
                    } else {
                        h_tmp_Up   = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp;
                        h_tmp_Down = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown;
                    }
                } else {
                    h_tmp_Up   = new TH1F(Form("h_DUMMY_%s_up_%i",  fRegions[0]->fSystNames[i_syst].c_str(),i_bin-1),"h_dummy",1,0,1);
                    h_tmp_Down = new TH1F(Form("h_DUMMY_%s_down_%i",fRegions[0]->fSystNames[i_syst].c_str(),i_bin-1),"h_dummy",1,0,1);
                }
                if(i_bin==1){
                    h_up.  push_back( new TH1F(Form("%s_TMP",h_tmp_Up->GetName()),  h_tmp_Up->GetTitle(),   fNRegions,0,fNRegions) );
                    h_down.push_back( new TH1F(Form("%s_TMP",h_tmp_Down->GetName()),h_tmp_Down->GetTitle(), fNRegions,0,fNRegions) );
                }
                h_up[i_syst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral(0,h_tmp_Up->GetNbinsX()+1) );
                h_down[i_syst]->SetBinContent( i_bin,h_tmp_Down->Integral(0,h_tmp_Down->GetNbinsX()+1) );
                // eventually add any other samples with the same title
                for(int j_smp=0;j_smp<fNSamples;j_smp++){
                    sh = fRegions[i_bin-1]->GetSampleHist( fSamples[j_smp]->fName );
                    if(idxVec[j_smp]==i_smp && i_smp!=j_smp){
                        if(isPostFit){
                            h_tmp_Up   = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp_postFit;
                            h_tmp_Down = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown_postFit;
                        }
                        else{
                            h_tmp_Up   = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp;
                            h_tmp_Down = sh->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown;
                        }
                        h_up[i_syst]  ->AddBinContent( i_bin,h_tmp_Up  ->Integral(0,h_tmp_Up->GetNbinsX()+1) );
                        h_down[i_syst]->AddBinContent( i_bin,h_tmp_Down->Integral(0,h_tmp_Down->GetNbinsX()+1) );
                    }
                }
            }
        }
        //
        //
        if(isPostFit)  g_err[i_smp] = BuildTotError( h_smp[i_smp], h_up, h_down, fRegions[0]->fSystNames, fFitResults->fCorrMatrix );
        else           g_err[i_smp] = BuildTotError( h_smp[i_smp], h_up, h_down, fRegions[0]->fSystNames );
    }
    
    //
    // Print samples except data
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if( fSamples[i_smp]->fType==Sample::DATA  ) continue;
        if( fSamples[i_smp]->fType==Sample::GHOST ) continue;
        if( fSamples[i_smp]->fType==Sample::SIGNAL && (fFitType==FitType::BONLY && isPostFit) ) continue;
        if(idxVec[i_smp]!=i_smp) continue;
        //
        // print values
        out << " | " << fSamples[i_smp]->fTitle << " | ";
        for(int i_bin=1;i_bin<=fNRegions;i_bin++){
            out << h_smp[i_smp]->GetBinContent(i_bin);
            out << " pm ";
            out << ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
            out << " | ";
        }
        out << endl;
    }
    
    //
    // Build tot
    //
    h_tot = new TH1F("h_Tot_","h_Tot", fNRegions,0,fNRegions);
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
//         if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot_postFit->Integral(0,fRegions[i_bin-1]->fTot_postFit->GetNbinsX()+1) );
//         else          h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot->Integral(        0,fRegions[i_bin-1]->fTot->GetNbinsX()+1) );
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot_postFit->IntegralAndError(0,fRegions[i_bin-1]->fTot_postFit->GetNbinsX()+1,intErr) );
        else          h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot->IntegralAndError(        0,fRegions[i_bin-1]->fTot->GetNbinsX()+1,        intErr) );
        h_tot->SetBinError( i_bin, intErr );
    }
    
    
    //
    //   Build error band
    // build the vectors of variations
    std::vector< TH1* > h_up;
    std::vector< TH1* > h_down;
    TH1* h_tmp_Up;
    TH1* h_tmp_Down;
    for(int i_syst=0;i_syst<(int)fRegions[0]->fSystNames.size();i_syst++){
        for(int i_bin=1;i_bin<=fNRegions;i_bin++){
            if(isPostFit){
                h_tmp_Up   = fRegions[i_bin-1]->fTotUp_postFit[i_syst];
                h_tmp_Down = fRegions[i_bin-1]->fTotDown_postFit[i_syst];
            }
            else{
                h_tmp_Up   = fRegions[i_bin-1]->fTotUp[i_syst];
                h_tmp_Down = fRegions[i_bin-1]->fTotDown[i_syst];
            }
            if(i_bin==1){
                h_up.  push_back( new TH1F(Form("h_%s_TMP",h_tmp_Up->GetName()),  h_tmp_Up->GetTitle(),   fNRegions,0,fNRegions) );
                h_down.push_back( new TH1F(Form("h_%s_TMP",h_tmp_Down->GetName()),h_tmp_Down->GetTitle(), fNRegions,0,fNRegions) );
            }
            h_up[i_syst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral(0,h_tmp_Up->GetNbinsX()+1) );
            h_down[i_syst]->SetBinContent( i_bin,h_tmp_Down->Integral(0,h_tmp_Down->GetNbinsX()+1) );
        }
    }
    if(isPostFit)  g_err_tot = BuildTotError( h_tot, h_up, h_down, fRegions[0]->fSystNames, fFitResults->fCorrMatrix );
    else           g_err_tot = BuildTotError( h_tot, h_up, h_down, fRegions[0]->fSystNames );
    //
    out << " | Total | ";
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        out << h_tot->GetBinContent(i_bin);
        out << " pm ";
        out << g_err_tot->GetErrorYhigh(i_bin-1);
        out << " | ";
    }
    out << endl;
    
    //
    // Print the data at last
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType!=Sample::DATA) continue;
        if(idxVec[i_smp]!=i_smp) continue;
        out << " | " << fSamples[i_smp]->fTitle << " | ";
        for(int i_bin=1;i_bin<=fNRegions;i_bin++){
            out << h_smp[i_smp]->GetBinContent(i_bin);
            out << " | ";
        }
    }
    out << endl;
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
            if(TtHFitter::DEBUGLEVEL>0) cout << fRegionsToPlot[i] << endl;
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
    gSystem->mkdir(fName.c_str());
    float Hp = 250; // height of one mini-plot, in pixels
    float Wp = 200; // width of one mini-plot, in pixels
    float H0 = 100; // height of the top label pad
    float H = H0 + nRows*Hp; // tot height of the canvas
    float W = nCols*Wp; // tot width of the canvas
    
    TCanvas *c = new TCanvas("c","c",W,H);
    TPad *pTop = new TPad("c0","c0",0,1-H0/H,1,1);
    pTop->Draw();
    pTop->cd();
    ATLASLabel(0.1,0.7,(char*)"Internal");
    myText(    0.1,0.4,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    myText(    0.1,0.1,1,Form("%s",fLabel.c_str()));
    
    c->cd();
    TPad *pBottom = new TPad("c1","c1",0,0,1,1-H0/H);
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
        if(regions[i]->fNSig > 0){
            if(regions[i]->fSig[0]!=0x0 && regions[i]->fSig[0]->fHist!=0x0)
                S[i] = regions[i]->fSig[0]->fHist->Integral();
        }
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
    for(int i=0;i<Nreg;i++){
        if(regions[i]==0x0) continue;
        pBottom->cd(i+1);
        string label = regions[i]->fShortLabel;
        h[i] = new TH1F(Form("h[%d]",i),label.c_str(),3,xbins);
        h[i]->SetBinContent(2,S[i]/sqrt(B[i]));
        h[i]->GetYaxis()->SetTitle("S / #sqrt{B}");
        h[i]->GetYaxis()->CenterTitle();
        //     h[i]->GetYaxis()->SetTitleSize(0.14);
        h[i]->GetYaxis()->SetLabelOffset(1.5*h[i]->GetYaxis()->GetLabelOffset());
        h[i]->GetYaxis()->SetTitleOffset(9*nRows/4.);
        //     h[i]->GetYaxis()->SetLabelSize(0.12);
        h[i]->GetXaxis()->SetTickLength(0);
        h[i]->GetYaxis()->SetNdivisions(3);
        yMax = TMath::Max(yMax,h[i]->GetMaximum());
        h[i]->GetXaxis()->SetLabelSize(0);
        h[i]->SetLineWidth(1);
        h[i]->SetLineColor(kBlack);
        if(regions[i]->fRegionType==Region::SIGNAL)          h[i]->SetFillColor(kRed+1);
        else if(regions[i]->fRegionType==Region::VALIDATION) h[i]->SetFillColor(kGray);
        else                                                 h[i]->SetFillColor(kAzure-4);
        h[i]->Draw();
        gPad->SetLeftMargin( gPad->GetLeftMargin()*2.4 );
        gPad->SetRightMargin(gPad->GetRightMargin()*0.1);
        gPad->SetTicky(0);
        gPad->RedrawAxis();
        tex->DrawLatex(0.4,0.85,label.c_str());
        float SoB = S[i]/B[i];
        string SB = Form("S/B = %.1f%%",(100.*SoB));
        tex->DrawLatex(0.4,0.72,SB.c_str());
    }
    //
    for(int i=0;i<Nreg;i++){
        if(regions[i]==0x0) continue;
        h[i]->SetMaximum(yMax*1.5);
    }
    //
//     c->SaveAs((fName+"/SignalRegions"+fSaveSuf+"."+fImageFormat).c_str());
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs((fName+"/SignalRegions"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
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
            if(TtHFitter::DEBUGLEVEL>0) cout << fRegionsToPlot[i] << endl;
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
    
    gSystem->mkdir((fName+"/PieChart").c_str());
    
    float Hp = 250; // height of one mini-plot, in pixels
    float Wp = 250; // width of one mini-plot, in pixels
    float H0 = 100; // height of the top label pad
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
    ATLASLabel(0.05,0.7,(char*)"Internal");
    myText(    0.05,0.4,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    myText(    0.05,0.1,1,Form("%s",fLabel.c_str()));
    
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
        if(regions[i]==0x0) continue;
        std::map < std::string, double > temp_map_for_region;
        std::map < std::string, int > temp_map_for_region_color;
        
//         for(int i_bkg=0;i_bkg<regions[i]->fNBkg;i_bkg++){
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
    TLegend *leg = new TLegend(0.7,0.1,0.95,0.90);
    if(map_for_legend.size()>4){
        leg -> SetNColumns(2);
    }
    leg -> SetLineStyle(0);
    leg -> SetFillStyle(0);
    leg -> SetLineColor(0);
    leg -> SetBorderSize(0);
    
    for ( const std::pair < std::string, int > legend_entry : map_for_legend ) {
        TH1F *dummy = new TH1F( ("legend_entry_" + legend_entry.first).c_str(), "",1,0,1);
        dummy -> SetFillColor(legend_entry.second);
        dummy -> SetLineColor(kBlack);
        dummy -> SetLineWidth(1);
        leg -> AddEntry(dummy,legend_entry.first.c_str(),"f");
    }
    leg -> Draw();
//     c->SaveAs((fName+"/Plots/PieChart" + fSuffix + ( isPostFit ? "_postFit" : "" ) + "."+fImageFormat).c_str());
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs((fName+"/PieChart" + fSuffix + ( isPostFit ? "_postFit" : "" ) + "."+TtHFitter::IMAGEFORMAT[i_format]).c_str());
}

//__________________________________________________________________________________
// turn to RooStat::HistFactory
void TtHFit::ToRooStat(bool makeWorkspace, bool exportOnly){
    
    //Suffix used for the regular bin transformed histogram
    const std::string suffix_regularBinning = "_regBin";
    
    if(TtHFitter::DEBUGLEVEL>0){
        cout << "--------------------------------" << endl;
        cout << "|      Export to RooStat       |" << endl;
        cout << "--------------------------------" << endl;
    }
    else{
        cout << "-------------------------------------------" << endl;
        cout << "Exporting to RooStats..." << endl;
    }
    
    RooStats::HistFactory::Measurement meas((fName+fSuffix).c_str(), (fName+fSuffix).c_str());
    meas.SetOutputFilePrefix((fName+"/RooStats/"+fName).c_str());
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
            cout << "Adding Channel: " << fRegions[i_ch]->fName << endl;
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
        if(hasData){
            if(TtHFitter::DEBUGLEVEL>0){
                cout << "  Adding Data: " << fRegions[i_ch]->fData->fHist->GetName() << endl;
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
                    cout << "  Adding Sample: " << fSamples[i_smp]->fName << endl;
                }
                RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());
                if(fUseStatErr) sample.ActivateStatError();
                sample.SetHistoName(h->fHistoName+suffix_regularBinning);
                sample.SetInputFile(h->fFileName);
                sample.SetNormalizeByTheory(fSamples[i_smp]->fNormalizedByTheory);
                // norm factors
                for(int i_norm=0;i_norm<h->fNNorm;i_norm++){
                    if(TtHFitter::DEBUGLEVEL>0){
                        cout << "    Adding NormFactor: " << h->fNormFactors[i_norm]->fName << endl;
                    }
                    sample.AddNormFactor( h->fNormFactors[i_norm]->fName,
                                         h->fNormFactors[i_norm]->fNominal,
                                         h->fNormFactors[i_norm]->fMin,
                                         h->fNormFactors[i_norm]->fMax  );
                }
                // systematics
                if(!fStatOnly){
                    for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
//                         if(fStatOnly)
//                             if(h->fSyst[i_syst]->fSystematic)
//                                 if(h->fSyst[i_syst]->fSystematic->fName!="Dummy") continue;
//                         std::cout << "OU?" << std::endl;
                        // add normalization part
                        if(TtHFitter::DEBUGLEVEL>0){
                            cout << "    Adding Systematic: " << h->fSyst[i_syst]->fName << endl;
                        }
                        
                        if(
                          (fThresholdSystPruning_Normalisation>-1 && (TMath::Abs(h->fSyst[i_syst]->fNormUp)>fThresholdSystPruning_Normalisation || TMath::Abs(h->fSyst[i_syst]->fNormDown)>fThresholdSystPruning_Normalisation)) ||
                            (fThresholdSystPruning_Normalisation==-1)
                          ){
                            sample.AddOverallSys( h->fSyst[i_syst]->fName,
                                                1+h->fSyst[i_syst]->fNormDown,
                                                1+h->fSyst[i_syst]->fNormUp   );
                        }
                        // eventually add shape part
                        if( h->fSyst[i_syst]->fIsShape && (fThresholdSystPruning_Shape==-1 || HistoTools::HasShape(h->fHist, h->fSyst[i_syst],fThresholdSystPruning_Shape) ) ){
                            sample.AddHistoSys( h->fSyst[i_syst]->fName,
                                              h->fSyst[i_syst]->fHistoNameShapeDown+suffix_regularBinning, h->fSyst[i_syst]->fFileNameShapeDown, "",
                                              h->fSyst[i_syst]->fHistoNameShapeUp+suffix_regularBinning,   h->fSyst[i_syst]->fFileNameShapeUp,   ""  );
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
    meas.PrintXML((fName+"/RooStats/").c_str());
    meas.CollectHistograms();
    meas.PrintTree();
    if(makeWorkspace) RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
}

//__________________________________________________________________________________
//
void TtHFit::DrawPruningPlot(){
    cout << "------------------------------------------------------" << endl;
    cout << "Drawing Pruning Plot ..." << endl;
    if(fSystematics.size()==0 || fStatOnly){
        std::cout << "TtHFit::INFO: Stat only fit => No Pruning plot generated." << std::endl;
//         cout << "... No systematics. Skipping." << endl;
        return;
    }
    vector< TH2F* > histPrun;
    int iReg = 0;
    int nSmp = 0;
    vector< Sample* > samplesVec;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        if(fSamples[i_smp]->fType==Sample::GHOST) continue;
        samplesVec.push_back(fSamples[i_smp]);
        nSmp++;
    }
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        if(fRegions[i_reg]->fRegionType!=Region::VALIDATION){
            histPrun.push_back( 
                new TH2F(Form("h_prun_%s", fRegions[i_reg]->fName.c_str()  ),
                         fRegions[i_reg]->fShortLabel.c_str(),
                         nSmp,0,nSmp, fNSyst,0,fNSyst
                        )
            );
            for(int i_smp=0;i_smp<nSmp;i_smp++){
                for(int i_syst=0;i_syst<fNSyst;i_syst++){
                    histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), -1 );
                }
                SampleHist *sh = fRegions[i_reg]->GetSampleHist(samplesVec[i_smp]->fName);
                if(sh!=0x0){
                    for(int i_syst=0;i_syst<fNSyst;i_syst++){
                        if(sh->HasSyst(fSystematics[i_syst]->fName)){
                            histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 0 );
                            // set to 1 if shape pruned away
                            if(sh->GetSystematic(fSystematics[i_syst]->fName)->fIsShape && 
                               fThresholdSystPruning_Shape>-1 && 
                               !HistoTools::HasShape(sh->fHist, sh->GetSystematic(fSystematics[i_syst]->fName),fThresholdSystPruning_Shape)
                              ){
                                histPrun[iReg]->AddBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 1 );
                            }
                            // set to 2 is normalization pruned away
                            if(fThresholdSystPruning_Normalisation>-1 && 
                               TMath::Abs(sh->GetSystematic(fSystematics[i_syst]->fName)->fNormUp)<fThresholdSystPruning_Normalisation &&
                               TMath::Abs(sh->GetSystematic(fSystematics[i_syst]->fName)->fNormDown)<fThresholdSystPruning_Normalisation
                              ){
                                histPrun[iReg]->AddBinContent( histPrun[iReg]->FindBin(i_smp,i_syst), 2 );
                            }
                        }
                    }
                }
            }
            iReg++;
        }
    }
    //
    // draw the histograms
    TCanvas *c = new TCanvas("c_pruning","Canvas - Pruning",200*(1+iReg),20*(fNSyst)+100+50);
    Int_t colors[] = {kGray, kGreen, kYellow, kOrange-3, kRed}; // #colors >= #levels - 1
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    c->Divide(1+iReg);
    for(int i_reg=0;i_reg<(int)histPrun.size();i_reg++){
        c->cd(i_reg+2);
        gPad->SetGridy();
        histPrun[i_reg]->Draw("COL");
        for(int i_bin=1;i_bin<=histPrun[i_reg]->GetNbinsX();i_bin++){
            histPrun[i_reg]->GetXaxis()->SetBinLabel(i_bin,samplesVec[i_bin-1]->fTitle.c_str());
        }
        for(int i_bin=1;i_bin<=histPrun[i_reg]->GetNbinsY();i_bin++){
            if(i_reg==0)
                histPrun[i_reg]->GetYaxis()->SetBinLabel(i_bin,TtHFitter::SYSTMAP[fSystematics[i_bin-1]->fName].c_str());
            else
                histPrun[i_reg]->GetYaxis()->SetBinLabel(i_bin,"");
        }
        histPrun[i_reg]->GetYaxis()->SetLabelOffset(0.04);
        gPad->SetBottomMargin(100./(20*(fNSyst)+100+50));
        gPad->SetTopMargin(50./(20*(fNSyst)+100+50));
        gPad->SetLeftMargin(0);
        gPad->SetRightMargin(0);
        histPrun[i_reg]->GetXaxis()->LabelsOption("v");
        histPrun[i_reg]->GetXaxis()->SetLabelSize( histPrun[i_reg]->GetXaxis()->GetLabelSize()*0.75 );
        histPrun[i_reg]->GetYaxis()->SetLabelSize( histPrun[i_reg]->GetYaxis()->GetLabelSize()*0.75 );
        gPad->SetTickx(0);
        gPad->SetTicky(0);
        histPrun[i_reg]->SetMinimum(-1);
        histPrun[i_reg]->SetMaximum(3);
        histPrun[i_reg]->GetYaxis()->SetTickLength(0);
        histPrun[i_reg]->GetXaxis()->SetTickLength(0);  
        gPad->SetGrid();
        myText(    0.1,1.-40./(20.*(fNSyst)+100.+50.),1,histPrun[i_reg]->GetTitle());
    }
    c->cd(1);
    myText(0.1,1.-10./(20.*(fNSyst)+100.+50.),1,fLabel.c_str());
    TLegend *leg = new TLegend(0.05,90./(20*(fNSyst)+100+50),0.95,0);
    TH1F* hGray = new TH1F("hGray","hGray",1,0,1);       hGray->SetFillColor(kGray);       hGray->SetLineWidth(0);
    TH1F* hYellow = new TH1F("hYellow","hYellow",1,0,1); hYellow->SetFillColor(kYellow);   hYellow->SetLineWidth(0);
    TH1F* hOrange = new TH1F("hOrange","hOrange",1,0,1); hOrange->SetFillColor(kOrange-3); hOrange->SetLineWidth(0);
    TH1F* hRed = new TH1F("hRed","hRed",1,0,1);          hRed->SetFillColor(kRed);         hRed->SetLineWidth(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hGray,"Not present","f");
//     leg->AddEntry(hGreen,"Kept","f");
    leg->AddEntry(hYellow, "Shape dropped","f");
    leg->AddEntry(hOrange, "Norm. dropped","f");
    leg->AddEntry(hRed, "Dropped","f");
    leg->SetTextSize(0.85*gStyle->GetTextSize());
    leg->Draw();
    //
//     c->SaveAs( (fName+"/Pruning"+fSaveSuf+"."+fImageFormat).c_str() );
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/Pruning"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );
}

//__________________________________________________________________________________
//
void TtHFit::Fit(){
    
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
    // Creating the combined model with the regions to fit only
    //
    RooWorkspace* ws = PerformWorkspaceCombination( regionsToFit );
    
    //
    // If needed, create a RooDataset object
    //
    RooDataSet* data = DumpData( ws, regionDataType, fFitNPValues, fFitPOIAsimov );
    
    //
    // Calls the PerformFit() function to actually do the fit
    //
    PerformFit( ws, regionsToFit, data);

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
    if (!fStatOnly){
      ws->saveSnapshot("InitialStateModelNuis",   *mc->GetNuisanceParameters());
    }

    //Be sure to take the initial values of the NP
    ws->loadSnapshot("InitialStateModelGlob");
    if (!fStatOnly){
      ws->loadSnapshot("InitialStateModelNuis");
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
std::map < std::string, double > TtHFit::PerformFit( RooWorkspace *ws, std::vector < std::string > &regionsToFit, RooDataSet* inputData, FitType fitType ){
    
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
    if(fFitType==BONLY || fitType==BONLY){
        fitTool -> ValPOI(0.);
        fitTool -> ConstPOI(true);
    } else if(fFitType==SPLUSB){
        fitTool -> ValPOI(1.);
        fitTool -> ConstPOI(false);
    }
    
    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());
    
    //
    // Creates the data object
    //
    RooDataSet* data = 0;
    if(inputData){
        data = inputData;
    } else {
        std::cout << "In TtHFit::PerformFit() function: you didn't provide inputData => will use the observed data !" << std::endl;
        data = (RooDataSet*)ws->data("obsData");
    }
    
    // Performs the fit
    gSystem -> mkdir((fName+"/Fits/").c_str(),true);
    fitTool -> MinimType("Minuit2");
    fitTool -> FitPDF( mc, simPdf, data );
    if(fitType==FitType::SPLUSB)fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fName+fSuffix+".txt");
    result = fitTool -> ExportFitResultInMap();
    
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
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        bool isToFit = false;
        for(unsigned int iRegion = 0; iRegion < regionsToFit.size(); ++iRegion){
            if(fRegions[i_ch] -> fName == regionsToFit[iRegion]){
                isToFit = true;
                break;
            }
        }
        if(isToFit){
            std::string fileName = fName+"/RooStats/"+fName+"_"+fRegions[i_ch]->fName+"_"+fName+fSuffix+"_model.root";
            TFile *rootFile = new TFile(fileName.c_str(),"read");
            RooWorkspace* m_ws = (RooWorkspace*) rootFile->Get((fRegions[i_ch]->fName).c_str());
            if(!m_ws){
                std::cout << "<!> Error in TtHFit::PerformWorkspaceCombination: The workspace (\"" << fRegions[i_ch] -> fName << "\") cannot be found in file " << fileName << ". Please check !" << std::endl;
            }
            vec_ws.push_back(m_ws);
            vec_chName.push_back(fRegions[i_ch] -> fName);
            if(!measurement){
                measurement = (RooStats::HistFactory::Measurement*) rootFile -> Get( (fName+fSuffix).c_str());
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
    if(fStatOnly){
        std::cout << "TtHFit::INFO: Stat only fit => No NP Pull plots generated." << std::endl;
        return;
    }
    //
    // plot the NP fit pull plot
    //
    ReadFitResults(fName+"/Fits/"+fName+fSuffix+".txt");
    if(fFitResults){
        std::set < std::string > npCategories;
        for(unsigned int i=0;i<fSystematics.size();i++){
            npCategories.insert(fSystematics[i]->fCategory);
        }
        for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
	    fFitResults->DrawPulls(fName+"/NuisPar"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format],"all");
        if(npCategories.size()>1){
            for( const std::string cat : npCategories ){
                std::string cat_for_name = cat;
                std::replace( cat_for_name.begin(), cat_for_name.end(), ' ', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '#', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '{', '_');
                std::replace( cat_for_name.begin(), cat_for_name.end(), '}', '_');
                for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
		    fFitResults->DrawPulls(fName+"/NuisPar_"+cat_for_name+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format],cat);
            }
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::PlotCorrelationMatrix(){
    if(fStatOnly){
        std::cout << "TtHFit::INFO: Stat only fit => No Correlation Matrix generated." << std::endl;
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TtHFitter::CORRELATIONTHRESHOLD)
    ReadFitResults(fName+"/Fits/"+fName+fSuffix+".txt");
    if(fFitResults){
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
        npValues = PerformFit( ws_forFit, regionsForFit, data, FitType::BONLY);
    }
    
    //
    // Create the final asimov dataset for limit setting
    //
    RooWorkspace* ws_forLimit = PerformWorkspaceCombination( regionsForLimit );
    data = DumpData( ws_forLimit, regionsForLimitDataType, npValues, npValues.find(fPOI)==npValues.end() ? fLimitPOIAsimov : npValues[fPOI] );
    
    //
    // Gets the measurement object in the original combined workspace (created with the "w" command)
    //
    const std::string originalCombinedFile = fName+"/RooStats/"+fName+"_combined_"+fName+fSuffix+"_model.root";
    TFile *f_origin = new TFile(originalCombinedFile.c_str(), "read");
//     RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get(fName.c_str());
    RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get((fName+fSuffix).c_str());
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
    string cmd;
    cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+(string)outputName+"\",\"combined\",\"ModelConfig\",\"ttHFitterData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+fSuffix+"\",0.95)'";
    
    //
    // Finally computing the limit
    //
    gSystem->Exec(cmd.c_str());
}


////__________________________________________________________________________________
////
//void TtHFit::GetLimit(){
//
//    //Checks if a data sample exists
//    bool hasData = false;
//    for(int i_smp=0;i_smp<fNSamples;i_smp++){
//        if(fSamples[i_smp]->fType==Sample::DATA){
//            hasData = true;
//            break;
//        }
//    }
////     string workspace = fName+"/RooStats/"+fName+"_combined_"+fName+"_model.root";
//    string workspace = fName+"/RooStats/"+fName+fSuffix+"_combined_"+fName+"_model.root";
//    if(hasData){
////         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\")'";
////         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_0\",\"./limits/\",\""+fName+"\",0.95)'";
////         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+workspace+"\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+fSaveSuf+"\",0.95)'";
//        string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+workspace+"\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+fSuffix+"\",0.95)'";
//        gSystem->Exec(cmd.c_str());
//    } else {
////         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_0\",\"./limits/\",\""+fName+"_blind\",0.95)'";
////         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+workspace+"\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+fSaveSuf+"\",0.95)'";
//        string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+workspace+"\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+fSuffix+"\",0.95)'";
//        gSystem->Exec(cmd.c_str());
//    }
//}

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
    
    //
    // Fills a vector of regions to consider for fit
    //
    std::vector < std:: string > regionsForFit;
    std::vector < std::string > regionsForSign;
    std::map < std::string, int > regionsForFitDataType;
    std::map < std::string, int > regionsForSignDataType;
    bool onlyUseRealData = true;
    for( unsigned int i_ch = 0; i_ch < fNRegions; i_ch++ ){
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
        npValues = PerformFit( ws_forFit, regionsForFit, data, FitType::BONLY);
    }
    
    //
    // Create the final asimov dataset for limit setting
    //
    RooWorkspace* ws_forSignificance = PerformWorkspaceCombination( regionsForSign );
    data = DumpData( ws_forSignificance, regionsForSignDataType, npValues, npValues.find(fPOI)==npValues.end() ? fLimitPOIAsimov : npValues[fPOI] );
    
    //
    // Gets the measurement object in the original combined workspace (created with the "w" command)
    //
    const std::string originalCombinedFile = fName+"/RooStats/"+fName+"_combined_"+fName+fSuffix+"_model.root";
    TFile *f_origin = new TFile(originalCombinedFile.c_str(), "read");
    RooStats::HistFactory::Measurement *originalMeasurement = (RooStats::HistFactory::Measurement*)f_origin -> Get(fName.c_str());
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
    string cmd = "root -l -b -q 'runSig.C(\""+(string)outputName+"\",\"combined\",\"ModelConfig\",\"ttHFitterData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\""+fName+"/Significance\")'";
    gSystem->Exec(cmd.c_str());
}

////__________________________________________________________________________________
////
//void TtHFit::GetSignificance(){
//    
//    //Checks if a data sample exists
//    bool hasData = false;
//    for(int i_smp=0;i_smp<fNSamples;i_smp++){
//        if(fSamples[i_smp]->fType==Sample::DATA){
//            hasData = true;
//            break;
//        }
//    }
////     string workspace = fName+"/RooStats/"+fName+"_combined_"+fName+"_model.root";
//    string workspace = fName+"/RooStats/"+fName+fSuffix+"_combined_"+fName+"_model.root";
//    if(hasData){
////         string cmd = "root -l -b -q 'runSig.C(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\"significance\")'";
//        string cmd = "root -l -b -q 'runSig.C(\""+workspace+"\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\""+fName+"/Significance\")'";
//        gSystem->Exec(cmd.c_str());
//        
//    } else {
////         string cmd = "root -l -b -q 'runSig.C(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\"significance\")'";
//        string cmd = "root -l -b -q 'runSig.C(\""+workspace+"\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\""+fName+"/Significance\")'";
//        gSystem->Exec(cmd.c_str());
//    }
//}

//__________________________________________________________________________________
//
void TtHFit::ReadFitResults(string fileName){
    cout << "------------------------------------------------------" << endl;
    cout << "Reading fit results from file " << fileName << endl;
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
    }
}

//__________________________________________________________________________________
//
void TtHFit::Print(){
    cout << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  TtHFit: " << fName << endl;
    cout << "      NtuplePaths ="; for(int i=0;i<(int)fNtuplePaths.size();i++) cout << " " << fNtuplePaths[i] << endl;
    cout << "      NtupleName  =";   cout << " " << fNtupleName << endl;
    cout << "      MCweight    =";   cout << " " << fMCweight << endl;
    cout << "      Selection   =";   cout << " " << fSelection << endl;
    cout << "      HistoPaths  ="; for(int i=0;i<(int)fHistoPaths.size();i++) cout << " " << fHistoPaths[i] << endl;
    cout << "      HistoName   =";   cout << " " << fHistoName << endl;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->Print();
    }
    cout << endl;
    cout << "-------------------------------------------" << endl;
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
    cout << "In DrawAndSaveSeparationPlots:" << endl; 
 
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    gSystem->mkdir((fName+"/Plots/Separation").c_str());

 
    // loop over regions
    for(unsigned int i_ch=0; i_ch < fRegions.size(); i_ch++)      {
      //cout << fRegions[i_ch]->fSig->fSample->fName << endl;
      //for(int i_bkg=0; i_bkg< fRegions[i_ch] -> fNBkg; i_bkg++){
      //      cout << fRegions[i_ch]->fBkg[i_bkg]->fSample->fName << endl; // will be used to discriminate other bkgs
      //  }

      // begin plotting
      TCanvas* dummy3 = new TCanvas("dummy3", "dummy3", 600,600);
      dummy3->cd();

      if(fRegions[i_ch]->fNSig==0){
          std::cout << "ERROR::TtHFit::DrawAndSaveSeparationPlots: No Signal Found" << std::endl;
          return;
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

      TLegend *legend3=new TLegend(0.48,0.72,0.94,0.87);
      legend3->SetTextFont(42);
      legend3->SetTextSize(0.043);
      legend3->AddEntry(bkg, "Total background" , "l");
      legend3->AddEntry(sig, "t#bar{t}H (m_{H} = 125 GeV)" , "l");
      legend3->SetFillColor(0) ;
      legend3->SetLineColor(0) ;
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

      std::string identS = fRegions[i_ch]->fLabel;      
      TLatex ls;
      ls.SetNDC();
      ls.SetTextSize(0.03);
      ls.SetTextColor(kBlack);
      ls.SetTextFont(42);
      ls.DrawLatex(0.20,0.73,identS.c_str());

      TLatex ls2;
      ls2.SetNDC();
      ls2.SetTextSize(0.03);
      ls2.SetTextColor(kBlack);
      ls2.SetTextFont(62);
      ls2.DrawLatex(0.20,0.78,"Single lepton");

      std::string cme = fRegions[i_ch]->fCmeLabel;
      std::string lumi = fRegions[i_ch]->fLumiLabel;

      TLatex ls3;
      ls3.SetNDC();
      ls3.SetTextSize(0.03);
      ls3.SetTextColor(kBlack);
      ls3.SetTextFont(42);
      ls3.DrawLatex(0.20,0.83, Form("#sqrt{s} = %s, %s", cme.c_str(), lumi.c_str()));

      ATLASLabelNew(0.20, 0.90,(char*)" Internal Simulation",kBlack, 0.03);

      TLatex ls4;
      ls4.SetNDC();
      ls4.SetTextSize(0.03);
      ls4.SetTextColor(kBlack);
      ls4.SetTextFont(42);
      ostringstream SEP;
      SEP.precision(3);
      SEP << "Separation: " << GetSeparation(sig,bkg)*100 << "%";
      ls4.DrawLatex(0.20, 0.69, SEP.str().c_str());
  
//       dummy3->SaveAs((fName+"/Plots/Separation/"+fRegions[i_ch]->fName+fSaveSuf+"_sep."+fImageFormat ).c_str());
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
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(NPnames=="all" || NPnames==fSystematics[i_syst]->fName || atoi(NPnames.c_str())==i_syst )
            nuisPars.push_back( fSystematics[i_syst]->fName );
    }
    
    //
    //Text files containing information necessary for drawing of ranking plot
    //     string outName = fName+"/Fits/NPRanking"+fSaveSuf;
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
    
    //
    // Initialize the FittingTool object
    //
    FittingTool *fitTool = new FittingTool();
    fitTool -> ValPOI(1.);
    fitTool -> ConstPOI(false);
    ReadFitResults(fName+"/Fits/"+fName+fSuffix+".txt");
    muhat = fFitResults -> GetNuisParValue( fPOI );
    //if(!hasData) muhat = 1.;  // FIXME -> Loic: Do we actually need that ?
    
    for(unsigned int i=0;i<nuisPars.size();i++){
        
        //Getting the postfit values of the nuisance parameter
        central = fFitResults -> GetNuisParValue(   nuisPars[i] );
        up      = fFitResults -> GetNuisParErrUp(   nuisPars[i] );
        down    = fFitResults -> GetNuisParErrDown( nuisPars[i] );
        outName_file <<  nuisPars[i] << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";
        
        //Set the NP to its post-fit *up* variation and refit to get the fitted POI
        fitTool -> FixNP( nuisPars[i], central + TMath::Abs(up  ) );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
        
        //Set the NP to its post-fit *down* variation and refit to get the fitted POI
        fitTool -> FixNP( nuisPars[i], central - TMath::Abs(down) );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        
        outName_file << muVarUp[nuisPars[i]]-muhat << "   " <<  muVarDown[nuisPars[i]]-muhat<< "  ";
        
        
        //Set the NP to its pre-fit *up* variation and refit to get the fitted POI (pre-fit impact on POI)
        fitTool -> FixNP( nuisPars[i], central + 1. );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarNomUp[ nuisPars[i] ]   = (fitTool -> ExportFitResultInMap())[ fPOI ];
        
        //Set the NP to its pre-fit *down* variation and refit to get the fitted POI (pre-fit impact on POI)
        fitTool -> FixNP( nuisPars[i], central - 1. );
        fitTool -> FitPDF( mc, simPdf, data );
        muVarNomDown[ nuisPars[i] ] = (fitTool -> ExportFitResultInMap())[ fPOI ];
        
        outName_file << muVarNomUp[nuisPars[i]]-muhat << "   " <<  muVarNomDown[nuisPars[i]]-muhat<< " "<<endl;
        
    }
    outName_file.close();

}

//____________________________________________________________________________________
//
void TtHFit::PlotNPRanking(){
    //
//     string fileToRead = fName+"/Fits/NPRanking"+fLoadSuf+".txt";
    string fileToRead = fName+"/Fits/NPRanking"+fSuffix+".txt";
    //
    // trick to merge the ranking outputs produced in parallel:
//     string cmd = " if [[ `ls "+fName+"/Fits/NPRanking"+fLoadSuf+"_*` != \"\" ]] ; ";
    string cmd = " if [[ `ls "+fName+"/Fits/NPRanking"+fSuffix+"_*` != \"\" ]] ; ";
    cmd       += " then rm "+fileToRead+" ; ";
//     cmd       += " cat "+fName+"/Fits/NPRanking"+fLoadSuf+"_* > "+fileToRead+" ; ";
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
    if(TtHFitter::DEBUGLEVEL>0) cout << "NP ordering..." << endl;
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
//     string npToExclude[] = {"SigXsecOverSM","gamma_","stat_"};
//     bool brazilian = true;
//     bool grayLines = false;
    
    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    TGraphAsymmErrors *g1 = new TGraphAsymmErrors();
    TGraphAsymmErrors *g2 = new TGraphAsymmErrors();
    TGraphAsymmErrors *g1a = new TGraphAsymmErrors();
    TGraphAsymmErrors *g2a = new TGraphAsymmErrors();
    
//     NuisParameter *par;
    int idx = 0;
    std::vector< string > Names;
    Names.clear();
    string parTitle;
    
//     for(unsigned int i = 0; i<parname.size(); ++i){
    for(unsigned int i = parname.size()-SIZE; i<parname.size(); ++i){
//         par = fNuisPar[i];
//         bool skip = false;
//         for(int ii=0; ii<sizeof(npToExclude)/sizeof(string); ii++){
//             if(par->fName.find(npToExclude[ii])!=string::npos){
//                 skip = true;
//                 continue;
//             }
//         }
//         if(skip) continue;
        
//         g->SetPoint(idx,par->fFitValue,idx+0.5);
//         g->SetPointEXhigh(idx, par->fPostFitUp);
//         g->SetPointEXlow( idx,-par->fPostFitDown);

//         idx++;
        
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
        
//         Poidown, Poiup
//         Poinomdown, Poinomup
        
//         parTitle = systMap[parname[i]];
           parTitle = TtHFitter::SYSTMAP[ parname[i] ];
//         h2->GetYaxis()->SetBinLabel(idx+1,parTitle.c_str());
        
//         Names.push_back(par->fName);
        Names.push_back(parTitle);
        
        idx ++;
        if(idx > max)  max = idx;      
    }

    TCanvas *c = new TCanvas("c","c",600,newHeight);
    c->SetTicks(0,0);
//     gPad->SetLeftMargin(0.33);
    gPad->SetLeftMargin(0.4);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);
    
    TH1F *h_dummy = new TH1F("h_dummy","h_dummy",10,xmin,xmax);
//     h_dummy->SetMaximum( SIZE );
//     h_dummy->SetMinimum( 0 );
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
    
//     TGaxis *axis_up = new TGaxis(-poimax,g->GetYaxis()->GetXmin()+0.05,poimax,SIZE+0.05,g->GetYaxis()->GetXmin(),SIZE+0.05,10);
//                Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax,
//                Double_t wmin, Double_t wmax, Int_t ndiv, Option_t *chopt,
//                Double_t gridlength
    TGaxis *axis_up = new TGaxis( -2, SIZE + (offsetUp1)/lineHeight, 2, SIZE + (offsetUp1)/lineHeight, -poimax,poimax, 510, "-" );
    axis_up->SetLabelOffset( 0.01 );
    axis_up->SetLabelSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetLabelFont(   gStyle->GetTextFont() );
    axis_up->Draw();
    axis_up->CenterTitle();
    axis_up->SetTitle("#Delta#mu");
    axis_up->SetTitleSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetTitleFont(   gStyle->GetTextFont() );
    
//     TLegend *leg = new TLegend(0.02,0.7,0.3,0.95);
//     leg->SetFillStyle(0);
//     leg->SetBorderSize(0);
//     leg->SetMargin(0.25);
//     leg->SetTextFont(gStyle->GetTextFont());
//     leg->SetTextSize(gStyle->GetTextSize());
//     leg->AddEntry(g,"(#hat{#theta}-#theta_{0})/#Delta#theta","lp");
//     leg->AddEntry(g1a,"#Delta#mu for #theta_{0}=+#Delta#theta","f");
//     leg->AddEntry(g2a,"#Delta#mu for #theta_{0}=-#Delta#theta","f");
//     leg->AddEntry(g1,"#Delta#mu for #theta_{0}=+#Delta#hat{#theta}","f");
//     leg->AddEntry(g2,"#Delta#mu for #theta_{0}=-#Delta#hat{#theta}","f");
//     leg->Draw();

//     TPad (const char *name, const char *title, Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Color_t color=-1, Short_t bordersize=-1, Short_t bordermode=-2)
//     TPad *pad0 = gPad;
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
    
//     ATLASLabelNew(0.42,0.80*(1.*newHeight/(offset + 6.*lineHeight)), (char*)"Internal", kBlack, gStyle->GetTextSize());
//     myText(       0.42,0.72*(1.*newHeight/(offset + 6.*lineHeight)), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    ATLASLabelNew(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.6*offsetUp1)/newHeight), (char*)"Internal", kBlack, gStyle->GetTextSize());
    myText(       0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.3*offsetUp1)/newHeight), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    
    gPad->RedrawAxis();
    
    for(int i_format=0;i_format<(int)TtHFitter::IMAGEFORMAT.size();i_format++)
        c->SaveAs( (fName+"/Ranking"+fSuffix+"."+TtHFitter::IMAGEFORMAT[i_format]).c_str() );
}

//____________________________________________________________________________________
//
void TtHFit::PrintSystTables(){
    std::cout << "TtHFit::INFO: Printing syt tables" << std::endl;
    for(int i_reg=0;i_reg<fNRegions;i_reg++){
        fRegions[i_reg]->PrintSystTable();
    }
}
