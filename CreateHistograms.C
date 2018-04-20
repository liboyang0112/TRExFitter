#include "TRandom.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSystem.h"

void CreateHistograms(){

  gSystem->mkdir("ExampleInputs");
  
  //
  // create signal histos
  //
  TFile *f_sig = new TFile("ExampleInputs/sig.root","RECREATE");
  TH1F *h_sig = new TH1F("h_sig","Signal",4,0,400);
  float yield_sig[] = {0.1,0.5,1.5,3.1};
  for(int i_bin=1;i_bin<=sizeof(yield_sig)/sizeof(float);i_bin++){
    h_sig->SetBinContent(i_bin,yield_sig[i_bin-1]);
  }
  h_sig->Write("HTj",TObject::kOverwrite);
  //in the VR now
  TH1F *h_sig_VR = (TH1F*) h_sig -> Clone("HTj_VR");
  h_sig_VR -> SetDirectory(0);
  h_sig_VR -> Scale(0.1);
  h_sig_VR -> Write("HTj_VR",TObject::kOverwrite);

  //
  // create sig systs
  //
  TH1F *h_sig_jesUp = new TH1F("h_sig_jesUp","Signal jesUp",4,0,400);
  float yield_sig_jesUp[] = {0.11,0.51,1.52,3.13};
  for(int i_bin=1;i_bin<=sizeof(yield_sig_jesUp)/sizeof(float);i_bin++){
    h_sig_jesUp->SetBinContent(i_bin,yield_sig_jesUp[i_bin-1]);
  }
  h_sig_jesUp->Write("HTj_jesUp",TObject::kOverwrite);
  //in the VR now
  TH1F *h_sig_jesUp_VR = (TH1F*) h_sig_jesUp -> Clone("h_sig_jesUp_VR");
  h_sig_jesUp_VR -> SetDirectory(0);
  h_sig_jesUp_VR -> Scale(0.1);
  h_sig_jesUp -> Write("HTj_VR_jesUp",TObject::kOverwrite);

  TH1F *h_sig_jesDown = new TH1F("h_sig_jesDown","Signal jesDown",4,0,400);
  float yield_sig_jesDown[] = {0.09,0.49,1.48,3.06};
  for(int i_bin=1;i_bin<=sizeof(yield_sig_jesDown)/sizeof(float);i_bin++){
    h_sig_jesDown->SetBinContent(i_bin,yield_sig_jesDown[i_bin-1]);
  }
  h_sig_jesDown->Write("HTj_jesDown",TObject::kOverwrite);
  //in the VR now
  TH1F *h_sig_jesDown_VR = (TH1F*) h_sig_jesDown -> Clone("h_sig_jesDown_VR");
  h_sig_jesDown_VR -> SetDirectory(0);
  h_sig_jesDown_VR -> Scale(0.1);
  h_sig_jesDown_VR -> Write("HTj_VR_jesDown",TObject::kOverwrite);

  //
  // create bkg1 histos
  //
  TFile *f_bkg1 = new TFile("ExampleInputs/bkg1.root","RECREATE");
  TH1F *h_bkg1 = new TH1F("h_bkg1","Background1",4,0,400);
  float yield_bkg1[] = {100.,75.,20.,10.};
  float stat_bkg1[] = {1.,1.,2.,0.5};
  for(int i_bin=1;i_bin<=sizeof(yield_bkg1)/sizeof(float);i_bin++){
    h_bkg1->SetBinContent(i_bin,yield_bkg1[i_bin-1]);
    h_bkg1->SetBinError(i_bin,stat_bkg1[i_bin-1]);
  }
  h_bkg1->Write("HTj",TObject::kOverwrite);
  //in the VR now
  TH1F *h_bkg1_VR = (TH1F*) h_bkg1 -> Clone("HTj_VR");
  h_bkg1_VR -> SetDirectory(0);
  h_bkg1_VR -> Scale(3.);
  h_bkg1_VR -> Write("HTj_VR",TObject::kOverwrite);

  //
  // create bkg1 systs
  //
  TH1F *h_bkg1_jesUp = new TH1F("h_bkg1_jesUp","Background1 jesUp",4,0,400);
  float yield_bkg1_jesUp[] = {90.,70.,30.,30.};
  for(int i_bin=1;i_bin<=sizeof(yield_bkg1_jesUp)/sizeof(float);i_bin++){
    h_bkg1_jesUp->SetBinContent(i_bin,yield_bkg1_jesUp[i_bin-1]);
    h_bkg1_jesUp->SetBinError(i_bin,stat_bkg1[i_bin-1]);
  }
  h_bkg1_jesUp->Write("HTj_jesUp",TObject::kOverwrite);
  //in the VR now
  TH1F *h_bkg1_jesUp_VR = (TH1F*) h_bkg1_jesUp -> Clone("h_bkg1_VR_jesUp");
  h_bkg1_jesUp_VR -> SetDirectory(0);
  h_bkg1_jesUp_VR -> Scale(3.);
  h_bkg1_jesUp_VR -> Write("HTj_VR_jesUp",TObject::kOverwrite);


  TH1F *h_bkg1_jesDown = new TH1F("h_bkg1_jesDown","Background1 jesDown",4,0,400);
  float yield_bkg1_jesDown[] = {110.,80.,15.,5.};
  for(int i_bin=1;i_bin<=sizeof(yield_bkg1_jesDown)/sizeof(float);i_bin++){
    h_bkg1_jesDown->SetBinContent(i_bin,yield_bkg1_jesDown[i_bin-1]);
    h_bkg1_jesDown->SetBinError(i_bin,stat_bkg1[i_bin-1]);
  }
  h_bkg1_jesDown->Write("HTj_jesDown",TObject::kOverwrite);
  //in the VR now
  TH1F *h_bkg1_jesDown_VR = (TH1F*) h_bkg1_jesDown -> Clone("h_bkg1_VR_jesDown");
  h_bkg1_jesDown_VR -> SetDirectory(0);
  h_bkg1_jesDown_VR -> Scale(3.);
  h_bkg1_jesDown_VR -> Write("HTj_VR_jesDown",TObject::kOverwrite);
  
  //
  // create bkg2 histos
  //
  TFile *f_bkg2 = new TFile("ExampleInputs/bkg2.root","RECREATE");
  TH1F *h_bkg2 = new TH1F("h_bkg2","Background2",4,0,400);
  float yield_bkg2[] = {10.,30.,50.,15.};
  for(int i_bin=1;i_bin<=sizeof(yield_bkg2)/sizeof(float);i_bin++){
    h_bkg2->SetBinContent(i_bin,yield_bkg2[i_bin-1]);
  }
  h_bkg2->Write("HTj",TObject::kOverwrite);
  //in the VR now
  TH1F *h_bkg2_VR = (TH1F*) h_bkg2 -> Clone("HTj_VR");
  h_bkg2_VR -> SetDirectory(0);
  h_bkg2_VR -> Scale(2.);
  h_bkg2_VR -> Write("HTj_VR",TObject::kOverwrite);
  
  //
  // create data histos
  //
  TFile *f_data = new TFile("ExampleInputs/data.root","RECREATE");
  TH1F *h_data = new TH1F("h_data","Data",4,0,400);
  for(int i_bin=1;i_bin<=sizeof(yield_sig)/sizeof(float);i_bin++){
    h_data->SetBinContent(i_bin,gRandom->Poisson( yield_sig[i_bin-1]+yield_bkg1[i_bin-1]+yield_bkg2[i_bin-1] ));
  }
  h_data->Write("HTj",TObject::kOverwrite);

  TH1F *h_data_VR = new TH1F("h_data_VR","Data",4,0,400);
  for(int i_bin=1;i_bin<=sizeof(yield_sig)/sizeof(float);i_bin++){
    const double background_prediction =  yield_sig[i_bin-1]*0.1 + yield_bkg1[i_bin-1]*3.0 + yield_bkg2[i_bin-1]*2.0;
    const int yield_data = gRandom->Poisson( background_prediction );
    h_data_VR->SetBinContent(i_bin, yield_data);
  }
  h_data_VR->Write("HTj_VR",TObject::kOverwrite);
}
