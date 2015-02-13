#include "TtHFitter/TtHFit.h"
#include "TSystem.h"

// -------------------------------------------------------------------------------------------------
// class TtHFit

TtHFit::TtHFit(string name){
  fNRegions = 0;
  fNSamples = 0;
  fNSyst = 0;
  fLumiErr = 0.000001;
  fName = name;
  fResultsFolder = "results/";
  fSelection = "";
  fNtuplePaths.clear();
  fPOI = "";
  fUseStatErr = false;
  fStatErrThres = 0.05;
}
TtHFit::~TtHFit(){}

void TtHFit::SetPOI(string name){
  fPOI = name;
}

void TtHFit::SetStatErrorConfig(bool useIt, float thres){
  fUseStatErr = useIt;
  fStatErrThres = thres;
}

void TtHFit::SetLumiErr(float err){
  fLumiErr = err;
}
  
Sample* TtHFit::NewSample(string name,int type){
  fSamples[fNSamples] = new Sample(name,type);
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
  fSystematics[fNSyst] = new Systematic(name);
  fNSyst ++;
  return fSystematics[fNSyst-1];
}

Region* TtHFit::NewRegion(string name){
  fRegions[fNRegions] = new Region(name);
//   // propagate stuff
//   for(int i_smp=0;i_smp<fNSamples;i_smp++){
//     fRegions[fNRegions]->AddSample(fSamples[i_smp]);
//   }
//   for(int i_syst=0;i_syst<fNSyst;i_syst++){
//     fRegions[fNRegions]->AddSystematic(fSystematics[i_syst]);
//   }
//   if(fSelection!="") fRegions[fNRegions]->AddSelection(fSelection);
  //
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
        cout << fullSelection << endl;
        cout << fullMCweight << endl;
        vector<string> s = CombinePathSufs(
                                        fRegions[i_ch]->fNtuplePathSuffs,
                                        fSamples[i_smp]->fSystematics[i_syst]->fNtuplePathsUp );
                                        cout << "ok" << endl;
        cout << s[0] << endl;
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
  bool isPostFit = opt.find("post")!=string::npos;
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    fRegions[i_ch]->fUseStatErr = fUseStatErr;
    if(isPostFit) fRegions[i_ch]->DrawPostFit(fFitResults,opt) -> SaveAs((fRegions[i_ch]->fName+"_postFit.png").c_str());
    else          fRegions[i_ch]->DrawPreFit(opt)              -> SaveAs((fRegions[i_ch]->fName+".png").c_str());
  }
}

// turn to RooStat::HistFactory
void TtHFit::ToRooStat(bool makeWorkspace, bool exportOnly){
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
  }
  else{
    meas.SetLumiRelErr(fLumiErr);
  }
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    if(TtHFitter::DEBUGLEVEL>0){
      cout << "Adding Channel: " << fRegions[i_ch]->fName << endl;
    }
    RooStats::HistFactory::Channel chan(fRegions[i_ch]->fName.c_str());
    if(TtHFitter::DEBUGLEVEL>0){
      cout << "  Adding Data: " << fRegions[i_ch]->fData->fHist->fName << endl;
    }
    chan.SetData(fRegions[i_ch]->fData->fHistoName, fRegions[i_ch]->fData->fFileName);
    chan.SetStatErrorConfig(fStatErrThres,"Gaussian");
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
      SampleHist* h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
      if( h != 0x0 && h->fSample->fType!=SampleType::Data){
        if(TtHFitter::DEBUGLEVEL>0){
          cout << "  Adding Sample: " << fSamples[i_smp]->fName << endl;
        }
        RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());
        if(fUseStatErr) sample.ActivateStatError();
        sample.SetHistoName(h->fHistoName);
        sample.SetInputFile(h->fFileName);
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
          sample.AddOverallSys( h->fSyst[i_syst]->fName,
                                1+h->fSyst[i_syst]->fNormDown,
                                1+h->fSyst[i_syst]->fNormUp   );
          // eventually add shape part
          if(h->fSyst[i_syst]->fIsShape){
            sample.AddHistoSys( h->fSyst[i_syst]->fName,
                                h->fSyst[i_syst]->fHistoNameShapeDown, h->fSyst[i_syst]->fFileNameShapeDown, "",
                                h->fSyst[i_syst]->fHistoNameShapeUp,   h->fSyst[i_syst]->fFileNameShapeUp,   ""  );
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
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    fRegions[i_ch]->Print();
  }
  cout << endl;
}
