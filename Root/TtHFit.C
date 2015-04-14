#include "TtHFitter/TtHFit.h"
#include "TtHFitter/HistoTools.h"

// -------------------------------------------------------------------------------------------------
// class TtHFit

TtHFit::TtHFit(string name){
    fName = name;
    fResultsFolder = "results/";
    fFitType = ControlSignalRegion;

    gSystem->mkdir(name.c_str());
    
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
}

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

void TtHFit::SetPOI(string name){
    fPOI = name;
}

void TtHFit::SetStatErrorConfig(bool useIt, float thres, string cons){
    fUseStatErr = useIt;
    fStatErrThres = thres;
    fStatErrCons = cons;
}

void TtHFit::SetLumiErr(float err){
    fLumiErr = err;
}

void TtHFit::SetLumi(const float lumi){
    fLumi = lumi;
}

void TtHFit::SetFitType(FitType type){
    fFitType = type;
}

Sample* TtHFit::NewSample(string name,int type){
    //fSamples[fNSamples] = new Sample(name,type);
    fSamples.push_back(new Sample(name,type));
    // propagate stuff
    //   for(int i_path=0;i_path<(int)fNtuplePaths.size();i_path++){
    //     fSamples[fNSamples]->AddNtuplePath(fNtuplePaths[i_path]);
    //   }
    //   if(fMCweight!="") fSamples[fNSamples]->SetMCweight(fMCweight);
    //   if(fNtupleName!="") fSamples[fNSamples]->fNtupleNames.push_back(fNtupleName);
    //
    fNSamples ++;
    return fSamples[fNSamples-1];
}

Systematic* TtHFit::NewSystematic(string name){
    //fSystematics[fNSyst] = new Systematic(name);
    fSystematics.push_back(new Systematic(name));
    fNSyst ++;
    return fSystematics[fNSyst-1];
}

Region* TtHFit::NewRegion(string name){
    //fRegions[fNRegions] = new Region(name);
    fRegions.push_back(new Region(name));
    
    //   // propagate stuff
    //   for(int i_smp=0;i_smp<fNSamples;i_smp++){
    //     fRegions[fNRegions]->AddSample(fSamples[i_smp]);
    //   }
    //   for(int i_syst=0;i_syst<fNSyst;i_syst++){
    //     fRegions[fNRegions]->AddSystematic(fSystematics[i_syst]);
    //   }
    //   if(fSelection!="") fRegions[fNRegions]->AddSelection(fSelection);
    //
    fRegions[fNRegions]->fFitName = fName;
    fNRegions ++;
    return fRegions[fNRegions-1];
}

// ntuple stuff
void TtHFit::AddNtuplePath(string path){
    fNtuplePaths.push_back(path);
}

void TtHFit::SetMCweight(string weight){
    fMCweight = weight;
}

void TtHFit::SetSelection(string selection){
    fSelection = selection;
}

void TtHFit::SetNtupleName(string name){
    fNtupleName = name;
}

// histogram stuff
void TtHFit::AddHistoPath(string path){
    fHistoPaths.push_back(path);
}

// ...

// apply smoothing to systematics
void TtHFit::SmoothSystematics(string syst){
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        for(int i_smp=0;i_smp<fRegions[i_ch]->fNSamples;i_smp++){
            fRegions[i_ch]->fSampleHists[i_smp]->SmoothSyst(syst);
        }
    }
}

// create new root file with all the histograms
void TtHFit::WriteHistos(string fileName,bool recreate){
    if(fileName=="") fileName = fName + "_histos.root";
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
                cout << "SampleHist[" << i_smp << "] not there." << endl;
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
            h->DrawSystPlot();
            h->WriteToFile();
        }
    }
}

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
            if(fSamples[i_smp]->fType==SampleType::Data) fullMCweight = "1";
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
            //
            // fix the bin contents FIXME (to avoid fit problems)
            if(fSamples[i_smp]->fType!=SampleType::Data) fRegions[i_ch]->fSampleHists[i_smp]->FixEmptyBins();
            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                // if not Overall only...
                if(fSamples[i_smp]->fSystematics[i_syst]->fType==SystType::Overall)
                    continue;
                cout << "Adding syst " << fSamples[i_smp]->fSystematics[i_syst]->fName << endl;
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
                //         cout << fullSelection << endl;
                //         cout << fullMCweight << endl;
                vector<string> s = CombinePathSufs(
                                                   fRegions[i_ch]->fNtuplePathSuffs,
                                                   fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsUp );
                //         cout << s[0] << endl;
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
                //
                fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,
                                                                                    hUp,hDown);
            }
        }
    }
    delete htmp;
}

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
            cout << "Reading " << fSamples[i_smp]->fName << endl;
            //
            // read nominal
            //
            // build a list of histograms to read
            fullPaths.clear();
            // Common::CreatePathsList( vector<string> paths, vector<string> pathSufs,
            //                          vector<string> files, vector<string> filesSuf,
            //                          vector<string> names, vector<string> namesSuf);
            fullPaths = CreatePathsList( fHistoPaths, fRegions[i_ch]->fHistoPathSuffs,
                                        fSamples[i_smp]->fHistoFiles, empty, // no histo file suffs for nominal (syst only)
                                        ToVec( fRegions[i_ch]->fHistoName ), empty  // same for histo name
                                        );
            for(int i_path=0;i_path<(int)fullPaths.size();i_path++){
                //         cout << " " << fullPaths[i_path] << endl;
                //         cout << "  " << fRegions[i_ch]->fHistoName << endl;
                //TH1* HistFromFile(string fileName,string histoName)
                //         htmp = (TH1F*)HistFromFile( fullPaths[i_path],fRegions[i_ch]->fHistoName);
                htmp = (TH1F*)HistFromFile( fullPaths[i_path] );
                
                //Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[i_ch]->fHistoBins) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin,htmp->GetName(),fRegions[i_ch]->fHistoBins));
                else if(fRegions[i_ch]->fHistoNBinsRebin != -1) htmp = (TH1F*)(htmp->Rebin(fRegions[i_ch]->fHistoNBinsRebin));
                
                if(fSamples[i_smp]->fType!=SampleType::Data && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                
                //Importing the histogram in TtHFitter
                if(i_path==0) h = (TH1F*)htmp->Clone(Form("h_%s_%s",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str()));
                else h->Add(htmp);
                htmp->~TH1F();
            }
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp], h );
            
            //
            // fix the bin contents FIXME (to avoid fit problems)
            if(fSamples[i_smp]->fType!=SampleType::Data) fRegions[i_ch]->fSampleHists[i_smp]->FixEmptyBins();
            //
            //  -----------------------------------
            //
            // read systematics (Shape and Histo)
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                // if not Overall only...
                if(fSamples[i_smp]->fSystematics[i_syst]->fType==SystType::Overall)
                    continue;
                cout << "Adding syst " << fSamples[i_smp]->fSystematics[i_syst]->fName << endl;
                //
                // Up
                //
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
                                            fSamples[i_smp]->fHistoFiles :
                                            fSamples[i_smp]->fSystematics[i_syst]->fHistoFilesUp ,
                                            // file suf
                                            fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufUp=="" ?
                                            empty :
                                            ToVec( fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufUp ),
                                            // name
                                            fSamples[i_smp]->fSystematics[i_syst]->fHistoNamesUp.size()==0 ?
                                            ToVec( fRegions[i_ch]->fHistoName ) :
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
                    
                    if(fSamples[i_smp]->fType!=SampleType::Data && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                    
                    //Importing histogram in TtHFitter
                    if(i_path==0) hUp = (TH1F*)htmp->Clone();
                    else hUp->Add(htmp);
                    htmp->~TH1F();
                }
                hUp->SetName(Form("h_%s_%s_%sUp",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
                //
                // Down
                //
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
                                            fSamples[i_smp]->fHistoFiles :
                                            fSamples[i_smp]->fSystematics[i_syst]->fHistoFilesDown ,
                                            // file suf
                                            fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufDown=="" ?
                                            empty :
                                            ToVec( fSamples[i_smp]->fSystematics[i_syst]->fHistoFileSufDown ),
                                            // name
                                            fSamples[i_smp]->fSystematics[i_syst]->fHistoNamesDown.size()==0 ?
                                            ToVec( fRegions[i_ch]->fHistoName ) :
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
                    
                    if(fSamples[i_smp]->fType!=SampleType::Data && fSamples[i_smp]->fNormalizedByTheory) htmp -> Scale(fLumi);
                    
                    //Importing histogram in TtHFitter
                    if(i_path==0) hDown = (TH1F*)htmp->Clone();
                    else hDown->Add(htmp);
                    htmp->~TH1F();
                }
                hDown->SetName(Form("h_%s_%s_%sDown",fRegions[i_ch]->fName.c_str(),fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fSystematics[i_syst]->fName.c_str()));
                
                SystematicHist *sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName)->AddHistoSyst(fSamples[i_smp]->fSystematics[i_syst]->fName,hUp,hDown);
                sh -> fSmoothType = fSamples[i_smp]->fSystematics[i_syst] -> fSmoothType;
                sh -> fSymmetrisationType = fSamples[i_smp]->fSystematics[i_syst] -> fSymmetrisationType;
                
            }
        }
    }
    delete htmp;
}

void TtHFit::ReadHistos(string fileName){
    if(fileName=="") fileName = fName + "_histos.root";
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
        cout << "  Reading region " << regionName << endl;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            sampleName = fSamples[i_smp]->fName;
            cout << "    Reading sample " << sampleName << endl;
            fRegions[i_ch]->SetSampleHist(fSamples[i_smp],regionName+"_"+sampleName,fileName);
            for(int i_syst=0;i_syst<fSamples[i_smp]->fNSyst;i_syst++){
                systName = fSamples[i_smp]->fSystematics[i_syst]->fName;
                cout << "      Reading syst " << systName << endl;
                sh = fRegions[i_ch]->GetSampleHist(sampleName);
                // norm only
                if(fSamples[i_smp]->fSystematics[i_syst]->fType == SystType::Overall){
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
    cout << "-----------------------------" << endl;
}

// old
void TtHFit::ReadAll(bool readNtuples,string fileName){
    //   for(int i_ch=0;i_ch<fNRegions;i_ch++){
    //     fRegions[i_ch]->SetAllSamples(readNtuples,fileName);
    //   }
    //   if(!readNtuples){
    //     SampleHist* h;
    //     SystematicHist* sh;
    //     for(int i_ch=0;i_ch<fNRegions;i_ch++){
    //       for(int i_smp=0;i_smp<fNSamples;i_smp++){
    //         h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
    //         h->fHistoName = h->fHist->GetName();
    //         h->fFileName = fileName;
    //         if( h != 0x0 && !h->fIsData ){
    //           for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
    //             sh = h->fSyst[i_syst];
    //             sh->fNormUp = sh->fHistUp->Integral() / h->fHist->Integral();
    //             sh->fNormDown = sh->fHistDown->Integral() / h->fHist->Integral();
    //             if(h->fSyst[i_syst]->fIsShape){
    // //                 h->fSystUp[i_syst]->Scale( 1./h->fSystNormUp[i_syst] );
    // //                 h->fSystDown[i_syst]->Scale( 1./h->fSystNormDown[i_syst] );
    // //                 h->fSystUp[i_syst]->Write("",TObject::kOverwrite);
    // //                 h->fSystDown[i_syst]->Write("",TObject::kOverwrite);
    //               h->fSyst[i_syst]->fHistoNameShapeUp = h->fSyst[i_syst]->fHistShapeUp->GetName();
    //               h->fSyst[i_syst]->fFileNameShapeUp = fileName;
    //               h->fSyst[i_syst]->fHistoNameShapeDown = h->fSyst[i_syst]->fHistShapeDown->GetName();
    //               h->fSyst[i_syst]->fFileNameShapeDown = fileName;
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
}

void TtHFit::DrawAndSaveAll(string opt){
    gSystem->mkdir(fName.c_str());
    bool isPostFit = opt.find("post")!=string::npos;
    if(isPostFit){
        ReadFitResults("xcheckResults/"+fName+"/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt");
    }
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->fUseStatErr = fUseStatErr;
        if(isPostFit){
            fRegions[i_ch]->DrawPostFit(fFitResults,opt) -> SaveAs((fName+"/"+fRegions[i_ch]->fName+"_postFit.png").c_str());
            
        }
        else{
            fRegions[i_ch]->DrawPreFit(opt)              -> SaveAs((fName+"/"+fRegions[i_ch]->fName+".png").c_str());            
        }
    }
    //   DrawSummary(opt+" log")->SaveAs("Summary.png");
}

TthPlot* TtHFit::DrawSummary(string opt){
    bool isPostFit = opt.find("post")!=string::npos;
    // build one bin per region
    TH1F* h_sig = 0;
    TH1F* h_data = 0;
    TH1F* h_bkg[MAXsamples];
    int Nbkg = 0;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        
        if(fSamples[i_smp]->fType==SampleType::Signal){
            h_sig = new TH1F(fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fTitle.c_str(), fNRegions,0,fNRegions);
            cout << "Adding Signal: " << h_sig->GetTitle() << endl;
            h_sig->SetLineColor(fRegions[0]->fSampleHists[i_smp]->fHist->GetLineColor());
            h_sig->SetFillColor(fRegions[0]->fSampleHists[i_smp]->fHist->GetFillColor());
            h_sig->SetLineWidth(fRegions[0]->fSampleHists[i_smp]->fHist->GetLineWidth());
            
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                h_sig->SetBinContent( i_bin,fRegions[i_bin-1]->fSampleHists[i_smp]->fHist->Integral() );
            }
        }
        else if(fSamples[i_smp]->fType==SampleType::Background){
            h_bkg[Nbkg] = new TH1F(fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fTitle.c_str(), fNRegions,0,fNRegions);
            cout << "Adding Bkg:    " << h_bkg[Nbkg]->GetTitle() << endl;
            h_bkg[Nbkg]->SetLineColor(fRegions[0]->fSampleHists[i_smp]->fHist->GetLineColor());
            h_bkg[Nbkg]->SetFillColor(fRegions[0]->fSampleHists[i_smp]->fHist->GetFillColor());
            h_bkg[Nbkg]->SetLineWidth(fRegions[0]->fSampleHists[i_smp]->fHist->GetLineWidth());
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                h_bkg[Nbkg]->SetBinContent( i_bin,fRegions[i_bin-1]->fSampleHists[i_smp]->fHist->Integral() );
            }
            Nbkg++;
        }
        else if(fSamples[i_smp]->fType==SampleType::Data){
            h_data = new TH1F(fSamples[i_smp]->fName.c_str(),fSamples[i_smp]->fTitle.c_str(), fNRegions,0,fNRegions);
            cout << "Adding Data:   " << h_data->GetTitle() << endl;
            for(int i_bin=1;i_bin<=fNRegions;i_bin++){
                h_data->SetBinContent( i_bin,fRegions[i_bin-1]->fData->fHist->Integral() );
            }
        }
    }
    
    //
    TthPlot *p = new TthPlot(fName+"_summary",900,700);
    p->SetXaxis("",false);
    p->SetChannel("Single Lepton");
    p->fATLASlabel = "Internal";
    
    //
    if(h_data){
        p->SetData(h_data, h_data->GetTitle());
    }
    
    p->AddSignal(h_sig,h_sig->GetTitle());
    
    //   p->AddNormSignal(h_sig,((string)h_sig->GetTitle())+"(norm)");
    for(int i=0;i<Nbkg;i++)
        p->AddBackground(h_bkg[i],h_bkg[i]->GetTitle());
    //
    //   BuildPreFitErrorHist();
    //
    //   p->SetTotBkg((TH1*)fTot);
    //   p->SetTotBkgAsym(fErr);
    for(int i_bin=1;i_bin<=fNRegions;i_bin++){
        p->SetBinLabel(i_bin,fRegions[i_bin-1]->fShortLabel.c_str());
    }
    p->Draw(opt);
    gSystem->mkdir(fName.c_str());
    p->SaveAs((fName+"/Summary.png").c_str());
    return p;
}

//void TtHFit::DrawSignalRegionsPlot(int nCols,int nRows,Region *regions[MAXregions]){
void TtHFit::DrawSignalRegionsPlot(int nCols,int nRows){
    DrawSignalRegionsPlot(nCols,nRows,fRegions);
}

void TtHFit::DrawSignalRegionsPlot(int nCols,int nRows, std::vector < Region* > &regions){
    
    gSystem->mkdir(fName.c_str());
    TCanvas *c = new TCanvas("c","c",200*nCols,100+250*nRows);
    //   c->SetTopMargin(100/(100+150*nCols));
    TPad *pTop = new TPad("c0","c0",0,1-100./(100.+150*nCols),1,1);
    pTop->Draw();
    ATLASLabel(0.1,0.90,(char*)"Internal");
    myText(    0.1,0.85,1,Form("#sqrt{s} = 8 TeV, 20.3 fb^{-1}"));
    myText(    0.1,0.80,1,Form("Single Lepton"));
    //
    TPad *pBottom = new TPad("c1","c1",0,0,1,1-100./(100.+150*nCols));
    pBottom->Draw();
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
    for(int i=0;i<Nreg;i++){
        pBottom->cd(i+1);
        string label = regions[i]->fShortLabel;
        gPad->SetLeftMargin( gPad->GetLeftMargin()*1.5 );
        h[i] = new TH1F(Form("h[%d]",i),label.c_str(),3,xbins);
        h[i]->SetBinContent(2,S[i]/sqrt(B[i]));
        h[i]->GetYaxis()->SetTitle("S / #sqrt{B}");
        h[i]->GetYaxis()->CenterTitle();
        //     h[i]->GetYaxis()->SetTitleSize(0.14);
        h[i]->GetYaxis()->SetLabelOffset(1.5*h[i]->GetYaxis()->GetLabelOffset());
        h[i]->GetYaxis()->SetTitleOffset(3.5);
        //     h[i]->GetYaxis()->SetLabelSize(0.12);
        h[i]->GetXaxis()->SetTickLength(0);
        h[i]->GetYaxis()->SetNdivisions(3);
        h[i]->SetMaximum(0.2);
        h[i]->GetXaxis()->SetLabelSize(0);
        h[i]->SetLineWidth(1);
        h[i]->SetLineColor(kBlack);
        if(i==Nreg-1) h[i]->SetFillColor(kRed+1);
        else          h[i]->SetFillColor(kAzure-4);
        h[i]->Draw();
        gPad->SetLeftMargin(gPad->GetLeftMargin()*1.25);
        gPad->SetTicky(0);
        gPad->RedrawAxis();
        tex->DrawLatex(0.35,0.85,label.c_str());
        float SoB = S[i]/B[i];
        string SB = Form("S/B = %.1f%%",(100.*SoB));
        tex->DrawLatex(0.35,0.72,SB.c_str());
    }
    //
    c->SaveAs((fName+"/SignalRegions.png").c_str());
}

void TtHFit::DrawPieChartPlot(){
    // still to implement...
}


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
    meas.SetOutputFilePrefix((fResultsFolder+"/"+fName).c_str());//"results/myMeasurement");
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
        
        if(fFitType==ControlRegion){
            if(fRegions[i_ch]->fRegionType==Region::SIGNAL || fRegions[i_ch]->fRegionType==Region::VALIDATION) continue;
        } else if(fFitType==ControlSignalRegion){
            if(fRegions[i_ch]->fRegionType==Region::VALIDATION) continue;
        }
        
        if(TtHFitter::DEBUGLEVEL>0){
            cout << "Adding Channel: " << fRegions[i_ch]->fName << endl;
        }
        RooStats::HistFactory::Channel chan(fRegions[i_ch]->fName.c_str());

        //Checks if a data sample exists
        bool hasData = false;
        for(int i_smp=0;i_smp<fNSamples;i_smp++){
            if(fSamples[i_smp]->fType==SampleType::Data){
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
            if( h != 0x0 && h->fSample->fType!=SampleType::Data){
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
    meas.PrintXML();
    meas.CollectHistograms();
    meas.PrintTree();
    if(makeWorkspace) RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
}


void TtHFit::Fit(){
    // PlotHistosBeforeFit=0,
    // PlotMorphingControlPlots=1, 
    // PlotHistosAfterFitEachSubChannel=2, 
    // PlotHistosAfterFitGlobal=3, 
    // PlotsNuisanceParametersVSmu=4, 
    // PlotsStatisticalTest=5
    int algo = 3;
    //     int algo = 0;
    string workspace = "results/"+fName+"_combined_"+fName+"_model.root";
    
    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==SampleType::Data){
            hasData = true;
            break;
        }
    }
    if(hasData){
        string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 0,\"%s\",\"./xcheckResults/%s/\",\"combined\",\"ModelConfig\",\"obsData\")'",
                      algo,workspace.c_str(),fName.c_str());
        gSystem->Exec(cmd.c_str());
    } else {
        string cmd = Form("root -l -b -q 'FitCrossCheckForLimits.C+(%d, 0, 1, 0,\"%s\",\"./xcheckResults/%s/\",\"combined\",\"ModelConfig\",\"asimovData\")'",
                          algo,workspace.c_str(),fName.c_str());
        gSystem->Exec(cmd.c_str());
    }
}

void TtHFit::PlotFittedNP(){
    // plot the NP fit plot
    string cmd = "python plotNP.py";
    cmd += " --outFile "+fName+"/NuisPar.png";
    cmd += " xcheckResults/"+fName+"/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt";
    gSystem->Exec(cmd.c_str());
}

void TtHFit::GetLimit(){
    
    //Checks if a data sample exists
    bool hasData = false;
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
        if(fSamples[i_smp]->fType==SampleType::Data){
            hasData = true;
            break;
        }
    }
    if(hasData){
        string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"obsData\")'";
        gSystem->Exec(cmd.c_str());
    } else {
        string cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\"results/"+fName+"_combined_"+fName+"_model.root\",\"combined\",\"ModelConfig\",\"asimovData\")'";
        gSystem->Exec(cmd.c_str());
    }
}

void TtHFit::ReadFitResults(string fileName){
    fFitResults = new FitResults();
    if(fileName.find(".txt")!=string::npos)
        fFitResults->ReadFromTXT(fileName);
}


void TtHFit::Print(){
    cout << endl;
    cout << "  TtHFit: " << fName << endl;
    cout << "      ntuplePaths ="; for(int i=0;i<(int)fNtuplePaths.size();i++) cout << " " << fNtuplePaths[i] << endl;
    cout << "      ntupleName  =";   cout << " " << fNtupleName << endl;
    cout << "      MCweight    =";   cout << " " << fMCweight << endl;
    cout << "      selection   =";   cout << " " << fSelection << endl;
    cout << "      histoPaths  ="; for(int i=0;i<(int)fHistoPaths.size();i++) cout << " " << fHistoPaths[i] << endl;
    cout << "      histoName   =";   cout << " " << fHistoName << endl;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
        fRegions[i_ch]->Print();
    }
    cout << endl;
}
