#include "TtHFitter/TtHFit.h"
#include "TtHFitter/HistoTools.h"

// -------------------------------------------------------------------------------------------------
// class TtHFit

//__________________________________________________________________________________
//
TtHFit::TtHFit(string name){
    fName = name;
    fLabel = "";
//     fResultsFolder = "results/";
//     fResultsFolder = fName+"/RooStats/";
    fFitType = CONTROLSIGNAL;

//     gSystem->mkdir(fName.c_str());
    
    fNRegions = 0;
    fNSamples = 0;
    fNSyst = 0;
    
    fPOI = "";
    fUseStatErr = false;
    fStatErrThres = 0.05;
    
    fLumi = 1.;
    fLumiErr = 0.000001;
    
    fThresholdSystPruning_Normalisation = -1;
    fThresholdSystPruning_Shape = -1;
    
    fNtuplePaths.clear();
    fMCweight = "";
    fSelection = "";
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
    
    fHistoCheckCrash = true;
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
void TtHFit::AddHistoPath(string path){
    fHistoPaths.push_back(path);
}

// ...

//__________________________________________________________________________________
// apply smoothing to systematics
void TtHFit::SmoothSystematics(string syst){
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
            fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(syst);
        }
    }
}

//__________________________________________________________________________________
// create new root file with all the histograms
void TtHFit::WriteHistos(string fileName,bool recreate){
  gSystem->mkdir( fName.c_str() );
//     if(fileName=="") fileName = fName + "_histos.root";
    if(fileName==""){
        fileName = fName + "/Histograms/" + fName + "_histos.root";
        gSystem->mkdir( (fName + "/Histograms/").c_str() );
    }
    cout << "-----------------------------" << endl;
    cout << "Writing histograms to file " << fileName << " ..." << endl;
    TDirectory *dir = gDirectory;
    TFile *f;
    if(recreate){
        f = new TFile(fileName.c_str(),"RECREATE");
        f->~TFile();
        dir->cd();
    }
    //
    SampleHist* h;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
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
//             h->DrawSystPlot();
            h->WriteToFile();
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots
void TtHFit::DrawSystPlots(){
    SampleHist* h;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            h->DrawSystPlot("all",TtHFitter::SYSTCONTROLPLOTS);
        }
    }
}

//__________________________________________________________________________________
// Build fit from config file
void TtHFit::ReadConfigFile(string fileName){
    fConfig->ReadFile(fileName);
    ConfigSet *cs; // to store stuff later
    string param;
    std::vector< string > vec;
    int type;
    //
    // set fit
    cs = fConfig->GetConfigSet("Fit");
    fName = cs->GetValue();
    param = cs->Get("Label");  if(param!="") fLabel = param;
                               else          fLabel = fName;
    SetPOI(cs->Get("POI"));
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
        if(cs->Get("NtuplePath")!="") { AddNtuplePath( cs->Get("NtuplePath") ); }
        param = cs->Get("NtuplePaths");
        if( param != "" ){
            std::vector<string> paths = Vectorize( param,',' );
            for(int i=0;i<(int)paths.size();i++){
                AddNtuplePath( paths[i] );
            }
        }
        SetMCweight(   cs->Get("MCweight")   );
        SetSelection(  cs->Get("Selection")  );
        SetNtupleName( cs->Get("NtupleName") );
    }
    param = cs->Get("LumiScale");  if( param != "" ) SetLumi( atof(param.c_str()) );
    param = cs->Get("FitType");    if( param != "" ){
        if(     param == "ControlSignalRegion" || param == "CONTROLSIGNAL")
            SetFitType(TtHFit::CONTROLSIGNAL);
        else if(param == "ControlRegion"       || param == "CONTROL")
            SetFitType(TtHFit::CONTROL);
        else{
            std::cerr << "Unknown FitType argument : " << cs->Get("FitType") << std::endl;
            return;
        }
    }
    param = cs->Get("SystPruningShape");  if( param != "")  fThresholdSystPruning_Shape         = atof(param.c_str());
    param = cs->Get("SystPruningNorm");   if( param != "")  fThresholdSystPruning_Normalisation = atof(param.c_str());
    param = cs->Get("IntCodeOverall");    if( param != "")  fIntCode_overall  = atoi(param.c_str());
    param = cs->Get("IntCodeShape");      if( param != "")  fIntCode_shape    = atoi(param.c_str());
    param = cs->Get("MCstatThreshold");   if( param != "")  SetStatErrorConfig( true,  atof(param.c_str()) );
                                          else              SetStatErrorConfig( false, 0.0 );
    param = cs->Get("DebugLevel");        if( param != "")  TtHFitter::SetDebugLevel( atoi(param.c_str()) );
    param = cs->Get("PlotOptions");       if( param != ""){
        vec = Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "YIELDS")!=vec.end() )  TtHFitter::SHOWYIELDS = true;
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
//         fRegionsToPlot.clear();
//         vec = Vectorize(param,',');
//         for(unsigned int i=0;i<vec.size();i++){
//         }
        fRegionsToPlot = Vectorize(param,',');
    }
    param = cs->Get("HistoChecks");  if(param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "NOCRASH" ){
            fHistoCheckCrash = false;
        }
    }
    
    //
    // set regions
    int nReg = 0;
    Region *reg;
    while(true){
        cs = fConfig->GetConfigSet("Region",nReg);
        if(cs==0x0) break;
        reg = NewRegion(cs->GetValue());
        reg->SetVariableTitle(cs->Get("VariableTitle"));
        reg->SetLabel(cs->Get("Label"),cs->Get("ShortLabel"));
        if(fInputType==0){
            param = cs->Get("HistoFile"); if(param!="") reg->fHistoFiles.push_back( param );
            param = cs->Get("HistoName"); if(param!="") reg->SetHistoName( param );
        }
        else if(fInputType==1){
          vector<string> variable = Vectorize(cs->Get("Variable"),',');
          reg->SetVariable(  variable[0], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()) );
          reg->AddSelection( cs->Get("Selection") );
//           reg->AddMCweight(  cs->Get("MCweight") );
          reg->fMCweight = cs->Get("MCweight"); // this will override the global MCweight, if any
          if(cs->Get("NtuplePathSuff")!="") { reg->fNtuplePathSuffs.clear(); reg->fNtuplePathSuffs.push_back( cs->Get("NtuplePathSuff") ); }
          param = cs->Get("NtuplePathSuffs");
          if( param != "" ){
              reg->fNtuplePathSuffs.clear();
              std::vector<string> paths = Vectorize( param,',' );
              for(int i=0;i<(int)paths.size();i++){
                  reg->fNtuplePathSuffs.push_back( paths[i] );
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
      if(cs->Get("Type")!=""){
          param = cs->Get("Type");
          std::transform(param.begin(), param.end(), param.begin(), ::toupper);
          if( param=="CONTROL" )     reg -> SetRegionType(Region::CONTROL);
          if( param=="VALIDATION" )  reg -> SetRegionType(Region::VALIDATION);
          if( param=="SIGNAL" )      reg -> SetRegionType(Region::SIGNAL);
      }
      nReg++;
    }
    //
    // set samples 
    int nSmp = 0;
    Sample *smp;
    while(true){
        cs = fConfig->GetConfigSet("Sample",nSmp);
        if(cs==0x0) break;
        type = Sample::BACKGROUND;
        if(cs->Get("Type")=="signal" || cs->Get("Type")=="SIGNAL") type = Sample::SIGNAL;
        if(cs->Get("Type")=="data"   || cs->Get("Type")=="DATA")   type = Sample::DATA;
        smp = NewSample(cs->GetValue(),type);
        smp->SetTitle(cs->Get("Title"));
        if(fInputType==0){
            param = cs->Get("HistoFile"); if(param!="") smp->AddHistoFile( param );
            param = cs->Get("HistoName"); if(param!="") smp->fHistoNames.push_back( param );
        }
        if(fInputType==1){
            param = cs->Get("NtupleFile");
            if(param!="") smp->AddNtupleFile( param );
            param = cs->Get("NtupleFiles");
            if(param!=""){
                smp->fNtupleFiles = Vectorize( param ,',' );
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
        }
        // ...
        nSmp++;
    }
    //
    // set systs
    int nSys = 0;
    Systematic *sys;
    Sample *sam;
    while(true){
        cs = fConfig->GetConfigSet("Systematic",nSys);
        if(cs==0x0) break;
        string samples_str = cs->Get("Samples");
        string exclude_str = cs->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        vector<string> samples = Vectorize(samples_str,',');
        vector<string> exclude = Vectorize(exclude_str,',');
        type = Systematic::HISTO;
        if(cs->Get("Type")=="overall" || cs->Get("Type")=="OVERALL")
            type = Systematic::OVERALL;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sam = fSamples[i_smp];
            if(sam->fType == Sample::DATA) continue;
            if(   (samples[0]=="all" || find(samples.begin(), samples.end(), sam->fName)!=samples.end() )
               && (exclude[0]==""    || find(exclude.begin(), exclude.end(), sam->fName)==exclude.end() )
            ){
                sys = sam->AddSystematic(cs->GetValue(),type);
                if(cs->Get("Title")!="") sys->fTitle = cs->Get("Title");
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
//                         if(cs->Get("NtupleFilesUp")!="")   sys->fNtupleFilesUp   = Vectorize( cs->Get("NtupleFilesUp"),  ',' );
//                         if(cs->Get("NtupleFilesDown")!="") sys->fNtupleFilesDown = Vectorize( cs->Get("NtupleFilesDown"),',' );
                        if(cs->Get("NtuplePathUp")!="")      sys->fNtuplePathsUp  .push_back(cs->Get("NtuplePathsUp"));
                        if(cs->Get("NtuplePathDown")!="")    sys->fNtuplePathsDown.push_back( cs->Get("NtuplePathsDown"));
                        if(cs->Get("NtuplePathSufUp")!="")   sys->fNtuplePathSufUp   = cs->Get("NtuplePathSufUp");
                        if(cs->Get("NtuplePathSufDown")!="") sys->fNtuplePathSufDown = cs->Get("NtuplePathSufDown");
                        if(cs->Get("NtupleFileUp")!="")      sys->fNtupleFilesUp  .push_back(cs->Get("NtupleFilesUp"));
                        if(cs->Get("NtupleFileDown")!="")    sys->fNtupleFilesDown.push_back( cs->Get("NtupleFilesDown"));
                        if(cs->Get("NtupleFileSufUp")!="")   sys->fNtupleFileSufUp   = cs->Get("NtupleFileSufUp");
                        if(cs->Get("NtupleFileSufDown")!="") sys->fNtupleFileSufDown = cs->Get("NtupleFileSufDown");
                        if(cs->Get("NtupleNameUp")!="")      sys->fNtupleNamesUp  .push_back(cs->Get("NtupleNamesUp"));
                        if(cs->Get("NtupleNameDown")!="")    sys->fNtupleNamesDown.push_back( cs->Get("NtupleNamesDown"));
                        if(cs->Get("NtupleNameSufUp")!="")   sys->fNtupleNameSufUp   = cs->Get("NtupleNameSufUp");
                        if(cs->Get("NtupleNameSufDown")!="") sys->fNtupleNameSufDown = cs->Get("NtupleNameSufDown");
                        if(cs->Get("WeightUp")!="")          sys->fWeightUp      = cs->Get("WeightUp");
                        if(cs->Get("WeightDown")!="")        sys->fWeightDown    = cs->Get("WeightDown");
                        if(cs->Get("WeightSufUp")!="")       sys->fWeightSufUp   = cs->Get("WeightSufUp");
                        if(cs->Get("WeightSufDown")!="")     sys->fWeightSufDown = cs->Get("WeightSufDown");
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
            }
        }
        // ...
        nSys++;
    }
}


//__________________________________________________________________________________
// for each region, add a SampleHist for each Sample in the Fit, reading from ntuples
void TtHFit::ReadNtuples(){
    TH1F* h;
    TH1F* hUp;
    TH1F* hDown;
    TH1F* htmp;
    //   string ntupleFullPath;
    string fullSelection;
    string fullMCweight;
    vector<string> fullPaths;
    vector<string> empty; empty.clear();
    //
    // loop on regions and samples
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            //
            // read nominal
            //
            // set selection and weight
            fullSelection = fSelection + " && " + fRegions[i_ch]->fSelection;
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
            // Common::CreatePathsList( vector<string> paths, vector<string> pathSufs,
            //                          vector<string> files, vector<string> filesSuf,
            //                          vector<string> names, vector<string> namesSuf);
            fullPaths = CreatePathsList( fNtuplePaths, fRegions[i_ch]->fNtuplePathSuffs,
                                        fSamples[i_smp]->fNtupleFiles, empty, // no ntuple file suffs for nominal (syst only)
                                        ToVec( fNtupleName ), empty  // same for ntuple name
                                        );
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                htmp = HistFromNtuple( fullPaths[i_path],
                                      fRegions[i_ch]->fVariable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                      fullSelection, fullMCweight);
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp], h );
            
            std::map < int, bool > applyCorrection;
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
                // Up
                //
                // Note: no need to change selection for systematics. If needed, can be done via weight (partially...)
                fullMCweight = fMCweight;
                if(fSamples[i_smp]->fSystematics[i_syst]->fWeightUp!="")
                    fullMCweight += " * "+fSamples[i_smp]->fSystematics[i_syst]->fWeightUp;
                else
                    fullMCweight += " * " + fSamples[i_smp]->fMCweight;
                if(fSamples[i_smp]->fSystematics[i_syst]->fWeightSufUp!="")
                    fullMCweight += " * "+fSamples[i_smp]->fSystematics[i_syst]->fWeightSufUp;
                if(fRegions[i_ch]->fMCweight!="")
                    fullMCweight += " * " + fRegions[i_ch]->fMCweight;
                vector<string> s = CombinePathSufs(
                                                   fRegions[i_ch]->fNtuplePathSuffs,
                                                   fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsUp );
                //
                fullPaths.clear();
                fullPaths = CreatePathsList(
                                            // path
                                            fNtuplePaths,
                                            // path suf
                                            CombinePathSufs(
                                                            fRegions[i_ch]->fNtuplePathSuffs,
                                                            fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsUp ),
                                            // file
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesUp.size()==0 ?
                                            fSamples[i_smp]->fNtupleFiles :
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesUp ,
                                            // file suf
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleFileSufUp=="" ?
                                            empty :
                                            ToVec( fSamples[i_smp]->fSystematics[i_syst]->fNtupleFileSufUp ),
                                            // name
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleNamesUp.size()==0 ?
                                            ToVec( fNtupleName ) :
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleNamesUp,
                                            // name suf
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameSufUp=="" ?
                                            empty :
                                            ToVec( fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameSufUp )
                                            );
                for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                    htmp = HistFromNtuple( fullPaths[i_path],
                                          fRegions[i_ch]->fVariable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                          fullSelection, fullMCweight);
                    if(i_path==0) hUp = (TH1F*)htmp->Clone();
                    else hUp->Add(htmp);
                    htmp->~TH1F();
                }
                hUp->SetName(Form("h_%s_%s_%sUp",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
                //
                // Down
                //
                // Note: no need to change selection for systematics. If needed, can be done via weight (partially...)
                fullMCweight = fMCweight;
                if(fSamples[i_smp]->fSystematics[i_syst]->fWeightDown!="")
                    fullMCweight += " * "+fSamples[i_smp]->fSystematics[i_syst]->fWeightDown;
                else
                    fullMCweight += " * " + fSamples[i_smp]->fMCweight;
                if(fSamples[i_smp]->fSystematics[i_syst]->fWeightSufDown!="")
                    fullMCweight += " * "+fSamples[i_smp]->fSystematics[i_syst]->fWeightSufDown;
                if(fRegions[i_ch]->fMCweight!="")
                    fullMCweight += " * " + fRegions[i_ch]->fMCweight;
                //
                fullPaths.clear();
                fullPaths = CreatePathsList(
                                            // path
                                            fNtuplePaths,
                                            // path suf
                                            CombinePathSufs(
                                                            fRegions[i_ch]->fNtuplePathSuffs,
                                                            fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsDown ),
                                            // file
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesDown.size()==0 ?
                                            fSamples[i_smp]->fNtupleFiles :
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleFilesDown ,
                                            // file suf
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleFileSufDown=="" ?
                                            empty :
                                            ToVec( fSamples[i_smp]->fSystematics[i_syst]->fNtupleFileSufDown ),
                                            // name
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleNamesDown.size()==0 ?
                                            ToVec( fNtupleName ) :
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleNamesDown,
                                            // name suf
                                            fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameSufDown=="" ?
                                            empty :
                                            ToVec( fSamples[i_smp]->fSystematics[i_syst]->fNtupleNameSufDown )
                                            );
                for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                    htmp = HistFromNtuple( fullPaths[i_path],
                                          fRegions[i_ch]->fVariable, fRegions[i_ch]->fNbins, fRegions[i_ch]->fXmin, fRegions[i_ch]->fXmax,
                                          fullSelection, fullMCweight);
                    if(i_path==0) hDown = (TH1F*)htmp->Clone();
                    else hDown->Add(htmp);
                    htmp->~TH1F();
                }
                hDown->SetName(Form("h_%s_%s_%sDown",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
                //
                // Histogram smoothing, Symmetrisation, Massaging...
                SystematicHist *sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                sh -> fSmoothType = fSamples[i_smp]->fSystematics[i_syst] -> fSmoothType;
                sh -> fSymmetrisationType = fSamples[i_smp]->fSystematics[i_syst] -> fSymmetrisationType;
                sh -> fSystematic = fSamples[i_smp]->fSystematics[i_syst];
                
                //
                // Histograms checking
                //
                for(unsigned int iBin = 1; iBin <= fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetNbinsX(); ++iBin ){
                    if( applyCorrection[iBin]){
                        sh -> fHistUp   -> SetBinContent(iBin,1e-06);
                        sh -> fHistDown -> SetBinContent(iBin,1e-06);
                    }
                }
                HistoTools::CheckHistograms( fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist /*nominal*/,
                                            sh /*systematic*/, fHistoCheckCrash /*cause crash if problem*/);
                
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
                if(fRegions[i_ch]->fHistoBins) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,htmp->GetName(),fRegions[i_ch]->fHistoBins));
                else if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                
                if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                
                //Importing the histogram in TtHFitter
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp], h );
            
            std::map < int, bool > applyCorrection;
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
            
            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                
                if(TtHFitter::DEBUGLEVEL>0) cout << "Adding syst " << fSamples[i_smp]->fSystematics[i_syst]->fName << endl;
                //
                // Up
                //
                // For histo syst:
                if(fSamples[i_smp]->fSystematics[i_syst]->fType==Systematic::HISTO){
                  fullPaths.clear();
                  fullPaths = CreatePathsList(
                                              // path
                                              fHistoPaths,
                                              // path suf
                                              CombinePathSufs(
                                                              fRegions[i_ch]->fHistoPathSuffs,
                                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoPathsUp ),
                                              // file
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoFilesUp.size()==0 ?
//                                               fSamples[i_smp]->fHistoFiles :
                                              histoFiles :
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoFilesUp ,
                                              // file suf
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufUp=="" ?
                                              empty :
                                              ToVec( fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufUp ),
                                              // name
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoNamesUp.size()==0 ?
//                                               ToVec( fRegions[i_ch]->fHistoName ) :
                                              histoNames :
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoNamesUp,
                                              // name suf
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoNameSufUp=="" ?
                                              empty :
                                              ToVec( fSamples[i_smp]->fSystematics[i_syst]->fHistoNameSufUp )
                                              );
                  for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                      htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                      //Pre-processing of histograms (rebinning, lumi scaling)
                      if(fRegions[i_ch]->fHistoBins) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,htmp->GetName(),fRegions[i_ch]->fHistoBins));
                      else if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                      
                      if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                      
                      //Importing histogram in TtHFitter
                      if(i_path==0) hUp = (TH1F*)htmp->Clone();
                      else hUp->Add(htmp);
                      htmp->~TH1F();
                  }
                }
                // For Overall syst
                else{
                  hUp = (TH1F*)h->Clone();
                  hUp->Scale(1+fSamples[i_smp]->fSystematics[i_syst]->fOverallUp);
                }
                hUp->SetName(Form("h_%s_%s_%sUp",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
                //
                // Down
                //
                // For histo syst:
                if(fSamples[i_smp]->fSystematics[i_syst]->fType==Systematic::HISTO){
                  fullPaths.clear();
                  fullPaths = CreatePathsList(
                                              // path
                                              fHistoPaths,
                                              // path suf
                                              CombinePathSufs(
                                                              fRegions[i_ch]->fHistoPathSuffs,
                                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoPathsDown ),
                                              // file
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoFilesDown.size()==0 ?
//                                               fSamples[i_smp]->fHistoFiles :
                                              histoFiles :
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoFilesDown ,
                                              // file suf
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufDown=="" ?
                                              empty :
                                              ToVec( fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufDown ),
                                              // name
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoNamesDown.size()==0 ?
//                                               ToVec( fRegions[i_ch]->fHistoName ) :
                                              histoNames :
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoNamesDown,
                                              // name suf
                                              fSamples[i_smp]->fSystematics[i_syst]->fHistoNameSufDown=="" ?
                                              empty :
                                              ToVec( fSamples[i_smp]->fSystematics[i_syst]->fHistoNameSufDown )
                                              );
                  for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                      htmp = (TH1F*)HistFromFile( fullPaths[i_path] ) ;
                      //Pre-processing of histograms (rebinning, lumi scaling)
                      if(fRegions[i_ch]->fHistoBins) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,htmp->GetName(),fRegions[i_ch]->fHistoBins));
                      else if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                      
                      if(fSamples[i_smp]->fType!=Sample::DATA && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                      
                      //Importing histogram in TtHFitter
                      if(i_path==0) hDown = (TH1F*)htmp->Clone();
                      else hDown->Add(htmp);
                      htmp->~TH1F();
                  }
                }
                // For Overall syst
                else{
                  hDown = (TH1F*)h->Clone();
                  hDown->Scale(1+fSamples[i_smp]->fSystematics[i_syst]->fOverallDown);
                }
                hDown->SetName(Form("h_%s_%s_%sDown",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
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
                for(unsigned int iBin = 1; iBin <= fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist->GetNbinsX(); ++iBin ){
                    if( applyCorrection[iBin]){
                        sh -> fHistUp   -> SetBinContent(iBin,1e-06);
                        sh -> fHistDown -> SetBinContent(iBin,1e-06);
                    }
                }
                HistoTools::CheckHistograms( fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->fHist /*nominal*/,
                                                sh /*systematic*/, fHistoCheckCrash /*cause crash if problem*/);
            }
        }
    }
    delete htmp;
}

//__________________________________________________________________________________
//
void TtHFit::ReadHistos(string fileName){
//     if(fileName=="") fileName = fName + "_histos.root";
    if(fileName=="") fileName = fName + "/Histograms/" + fName + "_histos.root";
    cout << "-----------------------------" << endl;
    cout << "Reading histograms from file " << fileName << " ..." << endl;
    TH1F* h;
    SampleHist *sh;
    SystematicHist *syh;
    string regionName;
    string sampleName;
    string systName;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        regionName = fRegions[i_ch]->fName;
        if(TtHFitter::DEBUGLEVEL>0) cout << "  Reading region " << regionName << endl;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sampleName = fSamples[i_smp]->fName;
            if(TtHFitter::DEBUGLEVEL>0) cout << "    Reading sample " << sampleName << endl;
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp],regionName+"_"+sampleName,fileName);
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                systName = fSamples[i_smp]->fSystematics[i_syst]->fName;
                if(TtHFitter::DEBUGLEVEL>0) cout << "      Reading syst " << systName << endl;
                sh = fRegions[i_ch]->GetSampleHist(sampleName);
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
        if(fFitType==CONTROL)            ReadFitResults(fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_conditionnal_mu0.txt");
        else if(fFitType==CONTROLSIGNAL) ReadFitResults(fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt");
    }
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->fUseStatErr = fUseStatErr;
        if(isPostFit){
            p = fRegions[i_ch]->DrawPostFit(fFitResults,opt);
            p->SaveAs(     (fName+"/Plots/"+fRegions[i_ch]->fName+"_postFit.png" ).c_str());
//             p->WriteToFile((fName+"/"+fRegions[i_ch]->fName+"_postFit.root").c_str());
        }
        else{
            p = fRegions[i_ch]->DrawPreFit(opt);
            p->SaveAs(     (fName+"/Plots/"+fRegions[i_ch]->fName+".png" ).c_str());            
//             p->WriteToFile((fName+"/"+fRegions[i_ch]->fName+".root").c_str());
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
    TH1F* h_sig = 0;
    TH1F* h_data = 0;
    TH1F* h_bkg[MAXsamples];
    TH1F *h_tot;
    TGraphAsymmErrors *g_err;
    int Nbkg = 0;
    //
    string name;
    string title;
    int lineColor;
    int fillColor;
    int lineWidth;
    double intErr; // to store the integral error
    TH1* h; // to store varius histograms temporary
    //
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        name = fSamples[i_smp]->fName.c_str();
        title = fSamples[i_smp]->fTitle.c_str();
        lineColor = fRegions[0]->fSampleHists[i_smp]->fHist->GetLineColor();
        fillColor = fRegions[0]->fSampleHists[i_smp]->fHist->GetFillColor();
        lineWidth = fRegions[0]->fSampleHists[i_smp]->fHist->GetLineWidth();
        //
        if(fSamples[i_smp]->fType==Sample::SIGNAL){
            h_sig = new TH1F(name.c_str(),title.c_str(), fNRegions,0,fNRegions);
            if(TtHFitter::DEBUGLEVEL>0) cout << "Adding Signal: " << h_sig->GetTitle() << endl;
            h_sig->SetLineColor(lineColor);
            h_sig->SetFillColor(fillColor);
            h_sig->SetLineWidth(lineWidth);
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                if(isPostFit)  h = fRegions[i_bin-1]->fSampleHists[i_smp]->fHist_postFit;
                else           h = fRegions[i_bin-1]->fSampleHists[i_smp]->fHist;
                h_sig->SetBinContent( i_bin,h->IntegralAndError(0,h->GetNbinsX()+1,intErr) );
                h_sig->SetBinError( i_bin,intErr );
            }
        }
        else if(fSamples[i_smp]->fType==Sample::BACKGROUND){
            h_bkg[Nbkg] = new TH1F(name.c_str(),title.c_str(), fNRegions,0,fNRegions);
            if(TtHFitter::DEBUGLEVEL>0) cout << "Adding Bkg:    " << h_bkg[Nbkg]->GetTitle() << endl;
            h_bkg[Nbkg]->SetLineColor(lineColor);
            h_bkg[Nbkg]->SetFillColor(fillColor);
            h_bkg[Nbkg]->SetLineWidth(lineWidth);
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                if(isPostFit)  h = fRegions[i_bin-1]->fSampleHists[i_smp]->fHist_postFit;
                else           h = fRegions[i_bin-1]->fSampleHists[i_smp]->fHist;
                h_bkg[Nbkg]->SetBinContent( i_bin,h->IntegralAndError(0,h->GetNbinsX()+1,intErr) );
                h_bkg[Nbkg]->SetBinError( i_bin,intErr );
            }
            Nbkg++;
        }
        else if(fSamples[i_smp]->fType==Sample::DATA){
            h_data = new TH1F(name.c_str(),title.c_str(), fNRegions,0,fNRegions);
            if(TtHFitter::DEBUGLEVEL>0) cout << "Adding Data:   " << h_data->GetTitle() << endl;
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                h_data->SetBinContent( i_bin,fRegions[i_bin-1]->fData->fHist->Integral() );
            }
        }
    }
    //
    TthPlot *p = new TthPlot(fName+"_summary",900,700);
    p->fShowYields = TtHFitter::SHOWYIELDS;
    p->fYmin = 10;
    p->SetXaxis("",false);
    p->AddLabel(fLabel);
    if(isPostFit) p->AddLabel("Post-Fit");
    else          p->AddLabel("Pre-Fit");
    p->fATLASlabel = "Internal";
    //
    if(h_data) p->SetData(h_data, h_data->GetTitle());
    if(h_sig) p->AddSignal(h_sig,h_sig->GetTitle());
    //   p->AddNormSignal(h_sig,((string)h_sig->GetTitle())+"(norm)");
    for(int i=0;i<Nbkg;i++)
        p->AddBackground(h_bkg[i],h_bkg[i]->GetTitle());
    //
    // Build tot
    h_tot = new TH1F("h_Tot_summary","h_Tot_summary", fNRegions,0,fNRegions);
    
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot_postFit->Integral() );
        else          h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot->Integral() );
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
                h_up.  push_back( new TH1F(Form("%s_TMP",h_tmp_Up->GetName()),  h_tmp_Up->GetTitle(),   fNRegions,0,fNRegions) );
                h_down.push_back( new TH1F(Form("%s_TMP",h_tmp_Down->GetName()),h_tmp_Down->GetTitle(), fNRegions,0,fNRegions) );
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
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        p->SetBinLabel(i_bin,fRegions[i_bin-1]->fShortLabel.c_str());
    }
    p->Draw(opt);
    //
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        if(TtHFitter::DEBUGLEVEL>0) cout << i_bin << ":\t" << h_tot->GetBinContent(i_bin) << "\t+" << g_err->GetErrorYhigh(i_bin-1) << "\t-" << g_err->GetErrorYlow(i_bin-1) << endl;
    }
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    if(isPostFit)  p->SaveAs((fName+"/Plots/Summary_postFit.png").c_str());
    else           p->SaveAs((fName+"/Plots/Summary.png").c_str());
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
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        name = fSamples[i_smp]->fName;//.c_str();
        title = fSamples[i_smp]->fTitle;//.c_str();
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
            if(isPostFit && fSamples[i_smp]->fType!=Sample::DATA)
                h0 = fRegions[i_bin-1]->fSampleHists[i_smp]->fHist_postFit;
            else
                h0 = fRegions[i_bin-1]->fSampleHists[i_smp]->fHist;
            h_smp[idxVec[i_smp]]->AddBinContent( i_bin,h0->IntegralAndError(0,h0->GetNbinsX()+1,intErr) );
//             h_smp[idxVec[i_smp]]->SetBinError( i_bin, sqrt( pow(h_smp[idxVec[i_smp]]->GetBinError(i_bin),2) + pow(intErr,2) ) );
        }
        titleVec.push_back(title);
    }
    //
    // add tot uncertainty on each sample
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(idxVec[i_smp]!=i_smp) continue;
        if(fSamples[i_smp]->fType==Sample::DATA) continue;
        // build the vectors of variations
        std::vector< TH1* > h_up;   h_up.clear();
        std::vector< TH1* > h_down; h_down.clear();
        TH1* h_tmp_Up;
        TH1* h_tmp_Down;
        for(int i_syst=0;i_syst<(int)fRegions[0]->fSystNames.size();i_syst++){
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                if(isPostFit){
//                     h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[i_smp]->fSyst[i_syst]->fHistUp_postFit;
//                     h_tmp_Down = fRegions[i_bin-1]->fSampleHists[i_smp]->fSyst[i_syst]->fHistDown_postFit;
                    h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[i_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp_postFit;
                    h_tmp_Down = fRegions[i_bin-1]->fSampleHists[i_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown_postFit;
                }
                else{
//                     h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[i_smp]->fSyst[i_syst]->fHistUp;
//                     h_tmp_Down = fRegions[i_bin-1]->fSampleHists[i_smp]->fSyst[i_syst]->fHistDown;
                    h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[i_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp;
                    h_tmp_Down = fRegions[i_bin-1]->fSampleHists[i_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown;
                }
                if(i_bin==1){
                    h_up.  push_back( new TH1F(Form("%s_TMP",h_tmp_Up->GetName()),  h_tmp_Up->GetTitle(),   fNRegions,0,fNRegions) );
                    h_down.push_back( new TH1F(Form("%s_TMP",h_tmp_Down->GetName()),h_tmp_Down->GetTitle(), fNRegions,0,fNRegions) );
                }
                h_up[i_syst]  ->SetBinContent( i_bin,h_tmp_Up  ->Integral(0,h_tmp_Up->GetNbinsX()+1) );
                h_down[i_syst]->SetBinContent( i_bin,h_tmp_Down->Integral(0,h_tmp_Down->GetNbinsX()+1) );
                // eventually add any other samples with the same title
                for(int j_smp=0;j_smp<fNSamples;j_smp++){
                    if(idxVec[j_smp]==i_smp && i_smp!=j_smp){
                        if(isPostFit){
//                             h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[j_smp]->fSyst[i_syst]->fHistUp_postFit;
//                             h_tmp_Down = fRegions[i_bin-1]->fSampleHists[j_smp]->fSyst[i_syst]->fHistDown_postFit;
                            h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[j_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp_postFit;
                            h_tmp_Down = fRegions[i_bin-1]->fSampleHists[j_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown_postFit;
                        }
                        else{
//                             h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[j_smp]->fSyst[i_syst]->fHistUp;
//                             h_tmp_Down = fRegions[i_bin-1]->fSampleHists[j_smp]->fSyst[i_syst]->fHistDown;
                            h_tmp_Up   = fRegions[i_bin-1]->fSampleHists[j_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistUp;
                            h_tmp_Down = fRegions[i_bin-1]->fSampleHists[j_smp]->GetSystematic(fRegions[0]->fSystNames[i_syst])->fHistDown;
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
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(idxVec[i_smp]!=i_smp) continue;
        //
        // print values
        out << " | " << fSamples[i_smp]->fTitle << " | ";
        for(int i_bin=1;i_bin<=fNRegions;i_bin++){
            out << h_smp[i_smp]->GetBinContent(i_bin);
            out << " +/- ";
            if(fSamples[i_smp]->fType!=Sample::DATA){
                out << ( g_err[i_smp]->GetErrorYhigh(i_bin-1) + g_err[i_smp]->GetErrorYlow(i_bin-1) )/2.;
            }
            else{
//                 out << h_smp[i_smp]->GetBinError(i_bin);
                out << sqrt(h_smp[i_smp]->GetBinContent(i_bin) );
            }
            out << " | ";
        }
        out << endl;
    }
    //
    // Build tot
    h_tot = new TH1F("h_Tot_","h_Tot", fNRegions,0,fNRegions);
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        if(isPostFit) h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot_postFit->Integral(0,fRegions[i_bin-1]->fTot_postFit->GetNbinsX()+1) );
        else          h_tot->SetBinContent( i_bin,fRegions[i_bin-1]->fTot->Integral(        0,fRegions[i_bin-1]->fTot->GetNbinsX()+1) );
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
        out << " +/- ";
        out << g_err_tot->GetErrorYhigh(i_bin-1);
        out << " | ";
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
//     TCanvas *c = new TCanvas("c","c",200*nCols,100+250*nRows);
//     TPad *pTop = new TPad("c0","c0",0,1-100./(100.+150*nCols),1,1);
    TCanvas *c = new TCanvas("c","c",W,H);
    TPad *pTop = new TPad("c0","c0",0,1-H0/H,1,1);
    pTop->Draw();
    pTop->cd();
    ATLASLabel(0.1,0.7,(char*)"Internal");
    myText(    0.1,0.4,1,Form("#sqrt{s} = 8 TeV, 20.3 fb^{-1}"));
    myText(    0.1,0.1,1,Form("%s",fLabel.c_str()));
    //
    c->cd();
    TPad *pBottom = new TPad("c1","c1",0,0,1,1-H0/H);
    pBottom->Draw();
    pBottom->cd();
    pBottom->Divide(nCols,nRows);
    int Nreg = nRows*nCols;
    if(Nreg>fNRegions) Nreg = fNRegions;
    TH1F* h[Nreg];
    float S[Nreg];
    float B[Nreg];
    double xbins[] = {0,0.1,0.9,1.0};
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(gStyle->GetTextSize());
    pBottom->cd(1);
    
    //
    // get the values
    for(int i=0;i<Nreg;i++){
        S[i] = regions[i]->fSig->fHist->Integral();
        B[i] = 0.;
        for(int i_bkg=0;i_bkg<regions[i]->fNBkg;i_bkg++){
            B[i] += regions[i]->fBkg[i_bkg]->fHist->Integral();
        }
    }
    //
    double yMax = 0;
    //
    for(int i=0;i<Nreg;i++){
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
        h[i]->SetMaximum(yMax*1.5);
    }
    //
    c->SaveAs((fName+"/SignalRegions.png").c_str());
}

//__________________________________________________________________________________
//
void TtHFit::DrawPieChartPlot(){
    // still to implement...
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
    else
        cout << "Exporting to RooStats..." << endl;
    
    RooStats::HistFactory::Measurement meas(fName.c_str(), fName.c_str());
    //   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    meas.SetOutputFilePrefix((fName+"/RooStats/"+fName).c_str());
    meas.SetExportOnly(exportOnly);
    meas.SetPOI(fPOI.c_str());
    meas.SetLumi(1.0);
    if(fLumiErr==0){
        meas.AddConstantParam("Lumi");
        meas.SetLumiRelErr(0.1);
    } else {
        meas.SetLumiRelErr(fLumiErr);
    }
    
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        
        if(fFitType==CONTROL){
            if(fRegions[i_ch]->fRegionType==Region::SIGNAL || fRegions[i_ch]->fRegionType==Region::VALIDATION) continue;
        } else if(fFitType==CONTROLSIGNAL){
            if(fRegions[i_ch]->fRegionType==Region::VALIDATION) continue;
        }
        
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
            if( h != 0x0 && h->fSample->fType!=Sample::DATA){
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
                for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
                    // add normalization part
                    if(TtHFitter::DEBUGLEVEL>0){
                        cout << "    Adding Systematic: " << h->fSyst[i_syst]->fName << endl;
                    }
                    
                    if(
                       (fThresholdSystPruning_Normalisation>-1 && (TMath::Abs(h->fSyst[i_syst]->fNormDown)>fThresholdSystPruning_Normalisation || TMath::Abs(h->fSyst[i_syst]->fNormDown)>fThresholdSystPruning_Normalisation)) ||
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
void TtHFit::Fit(){
    // PlotHistosBeforeFit=0,
    // PlotMorphingControlPlots=1, 
    // PlotHistosAfterFitEachSubChannel=2, 
    // PlotHistosAfterFitGlobal=3, 
    // PlotsNuisanceParametersVSmu=4, 
    // PlotsStatisticalTest=5
    int algo = 3;
//       int algo = 0;
//     string workspace = "results/"+fName+"_combined_"+fName+"_model.root";
    string workspace = fName+"/RooStats/"+fName+"_combined_"+fName+"_model.root";
    string outputDir = fName+"/FitResults/";
    
    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==Sample::DATA){
            hasData = true;
            break;
        }
    }
    if(hasData){
        if(fFitType==CONTROL){
            string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 1,\"%s\",\"./%s/\",\"combined\",\"ModelConfig\",\"obsData\")'",
                      algo,workspace.c_str(),outputDir.c_str());
            gSystem->Exec(cmd.c_str());
        } else if(fFitType==CONTROLSIGNAL){
            string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 0,\"%s\",\"./%s/\",\"combined\",\"ModelConfig\",\"obsData\")'",
                              algo,workspace.c_str(),outputDir.c_str());
            gSystem->Exec(cmd.c_str());
        }
    } else {
        if(fFitType==CONTROL){
            string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 1,\"%s\",\"./%s/\",\"combined\",\"ModelConfig\",\"asimovData\")'",
                              algo,workspace.c_str(),outputDir.c_str());
            gSystem->Exec(cmd.c_str());
        } else if(fFitType==CONTROLSIGNAL){
            string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 0,\"%s\",\"./%s/\",\"combined\",\"ModelConfig\",\"asimovData\")'",
                              algo,workspace.c_str(),outputDir.c_str());
            gSystem->Exec(cmd.c_str());
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::PlotFittedNP(){    
//     // plot the NP fit plot
//     string cmd = "python plotNP.py";
// //     cmd += " --outFile "+fName+"/NuisPar.png";
// //     if(fFitType==CONTROL) cmd += " xcheckResults/"+fName+"/TextFileFitResult/GlobalFit_fitres_conditionnal_mu0.txt";
// //     else if(fFitType==CONTROLSIGNAL) cmd += " xcheckResults/"+fName+"/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt";
//     if(fFitType==CONTROL)            cmd += " --outFile "+fName+"/NuisPar_Bonly.png";
//     else if(fFitType==CONTROLSIGNAL) cmd += " --outFile "+fName+"/NuisPar_SplusB.png";
//     if(fFitType==CONTROL)            cmd += " "+fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_conditionnal_mu0.txt";
//     else if(fFitType==CONTROLSIGNAL) cmd += " "+fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt";
//     gSystem->Exec(cmd.c_str());
    //
    if(fFitType==CONTROL)            ReadFitResults(fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_conditionnal_mu0.txt");
    else if(fFitType==CONTROLSIGNAL) ReadFitResults(fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt");
    if(fFitResults){
        if(fFitType==CONTROL)            fFitResults->DrawPulls(fName+"/NuisPar_Bonly.png");
        else if(fFitType==CONTROLSIGNAL) fFitResults->DrawPulls(fName+"/NuisPar_SplusB.png");
    }
}

//__________________________________________________________________________________
//
void TtHFit::PlotCorrelationMatrix(){
    if(fFitType==CONTROL)            ReadFitResults(fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_conditionnal_mu0.txt");
    else if(fFitType==CONTROLSIGNAL) ReadFitResults(fName+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt");
    if(fFitResults){
//         fFitResults->DrawCorrelationMatrix(fName+"/Plots/",TtHFitter::CORRELATIONTHRESHOLD);
        fFitResults->DrawCorrelationMatrix(fName+"/",TtHFitter::CORRELATIONTHRESHOLD);
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
    string workspace = fName+"/RooStats/"+fName+"_combined_"+fName+"_model.root";
    if(hasData){
//         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\")'";
//         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_0\",\"./limits/\",\""+fName+"\",0.95)'";
        string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+workspace+"\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+"\",0.95)'";
        gSystem->Exec(cmd.c_str());
    } else {
//         string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_0\",\"./limits/\",\""+fName+"_blind\",0.95)'";
        string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+workspace+"\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+"\",0.95)'";
        gSystem->Exec(cmd.c_str());
    }
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
    string workspace = fName+"/RooStats/"+fName+"_combined_"+fName+"_model.root";
    if(hasData){
//         string cmd = "root -l -b -q 'runSig.C(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\"significance\")'";
        string cmd = "root -l -b -q 'runSig.C(\""+workspace+"\",\"combined\",\"ModelConfig\",\"obsData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\""+fName+"/Significance\")'";
        gSystem->Exec(cmd.c_str());
        
    } else {
//         string cmd = "root -l -b -q 'runSig.C(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\"significance\")'";
        string cmd = "root -l -b -q 'runSig.C(\""+workspace+"\",\"combined\",\"ModelConfig\",\"asimovData\",\"asimovData_1\",\"conditionalGlobs_1\",\"nominalGlobs\",\""+fName+"\",\""+fName+"/Significance\")'";
        gSystem->Exec(cmd.c_str());
    }
}

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
    // assign to each NP in the FitResults a title, according to the syst in the fitter
    for(unsigned int i=0;i<fFitResults->fNuisPar.size();i++){
        for(unsigned int j=0;j<fSystematics.size();j++){
            if(fSystematics[j]->fName == fFitResults->fNuisPar[i]->fName){
                fFitResults->fNuisPar[i]->fTitle = fSystematics[j]->fTitle;
            }
            cout << endl;
        }
    }
}

//__________________________________________________________________________________
//
void TtHFit::Print(){
    cout << endl;
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
}

//__________________________________________________________________________________
//
Region* TtHFit::GetRegion(string name){
    for(unsigned int i=0;i<fRegions.size();i++){
        if(fRegions[i]->fName == name) return fRegions[i];
    }
    return 0x0;
}
