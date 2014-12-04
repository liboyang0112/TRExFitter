#include "TtHFitter/Common.h"

// -------------------------------------------------------------------------------------------------
// FUNCTIONS

TH1F* HistFromNtuple(string ntuple, string variable, int nbin, float xmin, float xmax, string selection, string weight){
  TH1F* h = new TH1F("h","h",nbin,xmin,xmax);
  cout << "  Extracting histogram from  " << ntuple << "  ..." << endl;
  TChain *t = new TChain();
  t->Add(ntuple.c_str());
  h->Sumw2();
  t->Draw( Form("%s>>h",variable.c_str()), Form("(%s)*(%s)",weight.c_str(),selection.c_str()), "goff");
  t->~TChain();
  return h;
}

TH1* HistFromFile(string fileName,string histoName){
  if(fileName=="") return 0x0;
  if(histoName=="") return 0x0;
  TH1* h;
  TDirectory *dir = gDirectory;
  TFile *f = new TFile(fileName.c_str());
  h = (TH1*)f->Get(fileName.c_str())->Clone();
  dir->cd();
}

void WriteHistToFile(TH1* h,string fileName,string option){
  TDirectory *dir = gDirectory;
  TFile *f = new TFile(fileName.c_str(),option.c_str());
  h->Write("",TObject::kOverwrite);
  f->~TFile();
  dir->cd();
}

void TtHFitter::SetDebugLevel(int level){
  DEBUGLEVEL = level;
}
