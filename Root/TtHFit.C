#include "TtHFitter/TtHFit.h"

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
}
TtHFit::~TtHFit(){}

void TtHFit::SetPOI(string name){
  fPOI = name;
}

void TtHFit::SetLumiErr(float err){
  fLumiErr = err;
}
  
Sample* TtHFit::NewSample(string name,int type){
  fSamples[fNSamples] = new Sample(name,type);
  // propagate stuff
  for(int i_path=0;i_path<(int)fNtuplePaths.size();i_path++){
    fSamples[fNSamples]->AddNtuplePath(fNtuplePaths[i_path]);
  }
  if(fMCweight!="") fSamples[fNSamples]->SetMCweight(fMCweight);
  if(fNtupleName!="") fSamples[fNSamples]->fNtupleNames.push_back(fNtupleName);
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

void TtHFit::ReadAll(bool readNtuples,string fileName){
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    fRegions[i_ch]->SetAllSamples(readNtuples,fileName);
  }
  if(!readNtuples){
    SampleHist* h;
    SystematicHist* sh;
    for(int i_ch=0;i_ch<fNRegions;i_ch++){
      for(int i_smp=0;i_smp<fNSamples;i_smp++){
        h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
        h->fHistoName = h->fHist->GetName();
        h->fFileName = fileName;
        if( h != 0x0 && !h->fIsData ){
          for(int i_syst=0;i_syst<h->fNSyst;i_syst++){
            sh = h->fSyst[i_syst];
            sh->fNormUp = sh->fHistUp->Integral() / h->fHist->Integral();
            sh->fNormDown = sh->fHistDown->Integral() / h->fHist->Integral();
            if(h->fSyst[i_syst]->fIsShape){
//                 h->fSystUp[i_syst]->Scale( 1./h->fSystNormUp[i_syst] );
//                 h->fSystDown[i_syst]->Scale( 1./h->fSystNormDown[i_syst] );
//                 h->fSystUp[i_syst]->Write("",TObject::kOverwrite);
//                 h->fSystDown[i_syst]->Write("",TObject::kOverwrite);
              h->fSyst[i_syst]->fHistoNameShapeUp = h->fSyst[i_syst]->fHistShapeUp->GetName();
              h->fSyst[i_syst]->fFileNameShapeUp = fileName;
              h->fSyst[i_syst]->fHistoNameShapeDown = h->fSyst[i_syst]->fHistShapeDown->GetName();
              h->fSyst[i_syst]->fFileNameShapeDown = fileName;
            }
          }
        }
      }
    }
  }
}

void TtHFit::DrawAndSaveAll(){
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    fRegions[i_ch]->DrawPreFit() -> SaveAs((fRegions[i_ch]->fName+".png").c_str());
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
  meas.SetLumiRelErr(fLumiErr);
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    if(TtHFitter::DEBUGLEVEL>0){
      cout << "Adding Channel: " << fRegions[i_ch]->fName << endl;
    }
    RooStats::HistFactory::Channel chan(fRegions[i_ch]->fName.c_str());
    if(TtHFitter::DEBUGLEVEL>0){
      cout << "  Adding Data: " << fRegions[i_ch]->fData->fHist->fName << endl;
    }
    chan.SetData(fRegions[i_ch]->fData->fHistoName, fRegions[i_ch]->fData->fFileName);
    for(int i_smp=0;i_smp<fNSamples;i_smp++){
      SampleHist* h = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
      if( h != 0x0 && h->fSample->fType!=SampleType::Data){
        if(TtHFitter::DEBUGLEVEL>0){
          cout << "  Adding Sample: " << fSamples[i_smp]->fName << endl;
        }
        RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());
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


void TtHFit::Print(){
  cout << endl;
  cout << "  TtHFit: " << fName << endl;
  for(int i_ch=0;i_ch<fNRegions;i_ch++){
    fRegions[i_ch]->Print();
  }
  cout << endl;
}
