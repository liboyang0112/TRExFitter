#include "TtHFitter/TthPlot.h"

#include "AtlasLabels.C"

using namespace std;

TthPlot::TthPlot(){
  fName = "c";
  Init();
}
TthPlot::TthPlot(string name){
  fName = name;
  Init();
}
void TthPlot::Init(){
  c = new TCanvas(fName.c_str(),fName.c_str(),500,600);
  //
  pad0 = new TPad("pad0","pad0",0,0.28,1,1,0,0,0);
  pad0->SetTickx(false);
  pad0->SetTicky(false);
  pad0->SetTopMargin(0.05);
  pad0->SetBottomMargin(0);
  pad0->SetLeftMargin(0.14);
  pad0->SetRightMargin(0.05);
  pad0->SetFrameBorderMode(0);
  //
  pad1 = new TPad("pad1","pad1",0,0,1,0.28,0,0,0);
  pad1->SetTickx(false);
  pad1->SetTicky(false);
  pad1->SetTopMargin(0.0);
  pad1->SetBottomMargin(0.37);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.05);
  pad1->SetFrameBorderMode(0);
  //
  pad1->Draw();
  pad0->Draw();
  pad0->cd();
  h_stack = new THStack("h_stack","h_stack");
  h_tot = 0;
  xtitle = "Variable [GeV]";
  ytitle = "Events";
//   leg  = new TLegend(0.54,0.93-2*(0.05),0.80,0.93);
//   leg1 = new TLegend(0.80,0.93-2*(0.05),0.94,0.93);
  leg_title = "e/#mu+jets";
  lumi = "4.7 fb^{-1}";
  yMaxScale = 3.;
  Chi2prob = -1;
  KSprob = -1;
  //
  h_signal = 0x0;
  h_normsig = 0x0;
  //
  fIsNjet = false;
}

void TthPlot::SetChannel(string name){
  leg_title = name;
}

void TthPlot::SetLumi(string name){
  lumi = name;
}

void TthPlot::SetXaxis(string name,bool isNjet){
  xtitle = name;
  fIsNjet = isNjet;
}

void TthPlot::SetYaxis(string name){
  ytitle = name;
}

void TthPlot::SetYmaxScale(float scale){
  yMaxScale = scale;
}

void TthPlot::SetData(TH1* h,string name){
  h_data = (TH1*)h;
  data_name = name;
}

void TthPlot::AddSignal(TH1* h,string name){
  h_signal = h;
  sample_name.push_back(name);
}

void TthPlot::AddNormSignal(TH1* h,string name){
  h_normsig = h;
  sample_name.push_back(name);
}

void TthPlot::AddBackground(TH1* h,string name){
  h_stack->Add(h);
  if(h_tot==0x0) h_tot = (TH1*)h->Clone();
  else h_tot->Add(h);
  h_bkg[sample_name.size()] = h;
  sample_name.push_back(name);
}

void TthPlot::SetTotBkg(TH1* h){
  h_tot = h;
  g_tot = new TGraphAsymmErrors(h);
}

void TthPlot::SetTotBkgAsym(TGraphAsymmErrors* g){
  g_tot = g;
  for(int i=1;i<h_tot->GetNbinsX()+1;i++){
    h_tot->SetBinContent(i,g_tot->GetY()[i-1]);
  }
}

void TthPlot::SetChi2KS(float chi2,float ks){
  Chi2prob = chi2;
  KSprob = ks;
}

void TthPlot::Draw(string options){
  pad0->cd();
  if(options.find("log")!=string::npos) pad0->SetLogy();
  //
  h_data->Draw("E");
  if(h_signal!=0x0) h_stack->Add(h_signal);
  h_stack->Draw("HISTsame");
  //
  g_tot->SetFillStyle(3354);
  g_tot->SetFillColor(kBlue-7);
  g_tot->SetLineColor(kWhite);
  g_tot->SetLineWidth(0);
  g_tot->SetMarkerSize(0);
  g_tot->Draw("sameE2");
  //
  if(h_normsig!=0x0){
    h_signal = (TH1*)h_normsig->Clone();
    h_normsig->Scale(h_tot->Integral()/h_normsig->Integral());
    h_normsig->SetLineColor(h_normsig->GetFillColor());
    h_normsig->SetFillColor(0);
    h_normsig->SetFillStyle(0);
    h_normsig->SetLineStyle(2);
    h_normsig->Draw("HISTsame");
  }
  //
  gStyle->SetEndErrorSize(4.);
  h_data->Draw("sameE1");

  h_data->GetXaxis()->SetTitle(xtitle.c_str());
  h_data->GetYaxis()->SetTitle(ytitle.c_str());
  h_data->GetYaxis()->SetTitleSize(0.05);
  if(fIsNjet){
    for(int i_bin=1;i_bin<h_data->GetNbinsX()+1;i_bin++){
      int nj = (int)h_data->GetXaxis()->GetBinCenter(i_bin);
      if(i_bin<h_data->GetNbinsX()) h_data->GetXaxis()->SetBinLabel( i_bin,Form("%d",nj) );
      else                          h_data->GetXaxis()->SetBinLabel( i_bin,Form("#geq%d",nj) );
    }
//     h_data->GetXaxis()->SetLabelSize(0.05);
  }
  h_data->GetYaxis()->SetTitleOffset(1.5);
  h_data->SetMarkerSize(1.2);
  h_data->SetMinimum(0);
//   h_data->SetMaximum(yMaxScale*TMath::Max(h_data->GetMaximum(),h_tot->GetMaximum()));
  h_data->SetMaximum(yMaxScale*h_tot->GetMaximum());
  pad0->RedrawAxis();
  
  ATLASLabel(0.2,0.89,"Internal");
  myText(0.2,0.8,1,"#sqrt{s} = 8 TeV, 20.3 fb^{-1}",0.045);
  myText(0.2,0.7,1,Form("%s",leg_title.c_str()),0.045);
  
  leg  = new TLegend(0.54,0.93-(sample_name.size()+2)*0.06,0.80,0.93);
  leg1 = new TLegend(0.80,leg->GetY1(),0.94,leg->GetY2());
  
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextAlign(32);
  leg->SetTextFont(42);
  leg->SetTextSize(0.042);
  leg->SetMargin(0.22);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextAlign(32);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.042);
  leg1->SetMargin(0.);
  leg->AddEntry(h_data,data_name.c_str(),"lep");
  leg1->AddEntry((TObject*)0,Form("%.1f",h_data->Integral()),"");
  leg->AddEntry(h_signal,sample_name[0].c_str(),"f");
  leg1->AddEntry((TObject*)0,Form("%.1f",h_signal->Integral()),"");
  for(int i_smp=sample_name.size()-1;i_smp>0;i_smp--){
    leg->AddEntry(h_bkg[i_smp], sample_name[i_smp].c_str(),"f");
    leg1->AddEntry((TObject*)0,Form("%.1f",h_bkg[i_smp]->Integral()),"");
  }
  leg->AddEntry((TObject*)0,"Total","");
  leg1->AddEntry((TObject*)0,Form("%.1f",h_tot->Integral()),"");
// //   leg1->AddEntry((TObject*)0,Form("%.1f",g_tot->Integral()),"");
// //   leg->AddEntry(h_tot,"Uncertainty","f");
  leg->AddEntry(g_tot,"Uncertainty","f");
  leg1->AddEntry((TObject*)0," ","");
  leg->Draw();
  leg1->Draw();

  
  // subpad
  TLine *hline;
  pad1->cd();
  pad1->GetFrame()->SetY1(2);
  TH1 *h_ratio = (TH1*)h_data->Clone("h_ratio");
  TH1 *h_ratio2 = (TH1*)h_tot->Clone("h_ratio2");
  TGraphAsymmErrors *g_ratio2 = (TGraphAsymmErrors*)g_tot->Clone("g_ratio2");
  TH1 *h_tot_nosyst = (TH1*)h_tot->Clone("h_tot_nosyst");
  for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){
    h_tot_nosyst->SetBinError(i_bin,0);
  }
  h_ratio->SetTitle("Data/MC");
  h_ratio->GetYaxis()->SetTitle("Data / MC");
  h_ratio->GetYaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetTitleOffset(0.65);
  h_ratio->GetYaxis()->SetLabelSize(0.12); // 0.04
  h_ratio->Divide(h_tot_nosyst);
  for(int i_bin=1;i_bin<h_tot_nosyst->GetNbinsX()+1;i_bin++){
    g_ratio2->SetPoint(i_bin-1,g_ratio2->GetX()[i_bin-1],g_ratio2->GetY()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
    g_ratio2->SetPointEXlow(i_bin-1,g_ratio2->GetEXlow()[i_bin-1]);
    g_ratio2->SetPointEXhigh(i_bin-1,g_ratio2->GetEXhigh()[i_bin-1]);
    g_ratio2->SetPointEYlow(i_bin-1,g_ratio2->GetEYlow()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
    g_ratio2->SetPointEYhigh(i_bin-1,g_ratio2->GetEYhigh()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
  }
  hline = new TLine(h_ratio->GetXaxis()->GetXmin(),1,h_ratio->GetXaxis()->GetXmax(),1);
  hline->SetLineColor(kRed);
  hline->SetLineWidth(2);
  hline->SetLineStyle(2);
  h_ratio->SetMarkerStyle(24);
  h_ratio->SetMarkerSize(0.8);
  gStyle->SetEndErrorSize(4.);
  h_ratio->GetYaxis()->CenterTitle();
  h_ratio->GetYaxis()->SetNdivisions(504,false);
  h_ratio->Draw("0E1");
  hline->Draw();
  h_ratio->SetMinimum(0.);
  h_ratio->SetMaximum(2.);
  if(options.find("prefit")!=string::npos){
    h_ratio->SetMinimum(0.00);
    h_ratio->SetMaximum(2.00);
  }
  else{
    h_ratio->SetMinimum(0.50);
    h_ratio->SetMaximum(1.50);
  }
  h_ratio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
  h_ratio->GetXaxis()->SetTitleSize(0.14);
  h_ratio->GetXaxis()->SetTitleOffset(1.);
  h_data->GetXaxis()->SetTitle("");
  h_ratio->GetXaxis()->SetLabelSize(0.14);
  if(fIsNjet) h_ratio->GetXaxis()->SetLabelSize(0.2);
  h_data->GetXaxis()->SetLabelSize(0);
  gPad->RedrawAxis();
  // to hide the upper limit
  TLine line(0.01,1,0.1,1);
  line.SetLineColor(kWhite);
  line.SetLineWidth(20);
  line.DrawLineNDC(0.07,1,0.135,1);
  g_ratio2->Draw("sameE2");
  //

  // by Valerio
  //// putting the stuff per Rick's requests.
  for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){
    if (i_bin==0 || i_bin>h_tot_nosyst->GetNbinsX()) continue;
    float val=h_ratio->GetBinContent(i_bin);    
    int isUp=0; //1==up, 0==nothing, -1==down
    if (options.find("prefit")!=string::npos) {
      if ( val>=2 ) isUp=1;
    } else {
      if ( val<0.5 ) isUp=-1;
      else if (val>1.5 ) isUp=1;
      if (val==0) isUp=0;
    }
    if (isUp!=0) {
      TArrow *arrow;
      if (isUp==1) arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),1.41, h_ratio->GetXaxis()->GetBinCenter(i_bin),1.45,0.030,"|>");
      else         arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),0.54, h_ratio->GetXaxis()->GetBinCenter(i_bin),0.51,0.030,"|>");
      arrow->SetFillColor(2);
      arrow->SetFillStyle(1001);
      arrow->SetLineColor(2);
      arrow->SetAngle(40);
      arrow->Draw();

      TLine *fix=0;
      TLine *fixUp=0;
      TLine *fixDo=0;
      cout << " bin: " << i_bin << "  with val: " << val << " and lower: " << val-h_ratio->GetBinError(i_bin) <<  endl;
      cout << " bin: " << i_bin << "  with val: " << val << " and upper: " << val+h_ratio->GetBinError(i_bin) <<  endl;
      if (val-h_ratio->GetBinError(i_bin)<1.5) {
        float x = h_ratio->GetXaxis()->GetBinCenter(i_bin);
        float y0 = TMath::Max(val-h_ratio->GetBinError(i_bin),0.50);
        float y1 = TMath::Min(val+h_ratio->GetBinError(i_bin),1.50);
        fix = new TLine( x, y0, x, y1 );
        fixUp = new TLine( x-0.07, y1, x+0.07, y1 );
        fixDo = new TLine( x-0.07, y0, x+0.07, y0 );
      }
      if (fix!=0) {
        fix->SetLineColor(1);
        fix->SetLineWidth(1);
        fix->Draw("SAME");
        fixUp->SetLineColor(1);
        fixUp->SetLineWidth(1);
        fixUp->Draw("SAME");
        fixDo->SetLineColor(1);
        fixDo->SetLineWidth(1);
        fixDo->Draw("SAME");
      }
    }
  }
  // ---

  pad1->cd();
  TLatex *KSlab = new TLatex();
  KSlab->SetNDC(1);
  KSlab->SetTextFont(42);
  KSlab->SetTextSize(0.1);
//  if(KSprob!=0)
//    KSlab->DrawLatex(0.6,0.9,Form("KS = %.2f",KSprob));
//  if(KStest == "chi2")
//    KSlab->DrawLatex(0.6,0.9,Form("#chi^{2} = %.2f",Chi2));
  if(KSprob >= 0 && Chi2prob >= 0)
    KSlab->DrawLatex(0.15,0.9,Form("#chi^{2} prob = %.2f,   KS prob = %.2f",Chi2prob,KSprob));
  //
  pad0->cd();
}

void TthPlot::SaveAs(string name){
  c->SaveAs(name.c_str());
}

void TthPlot::WriteToFile(string name){
  TFile *f = new TFile(name.c_str(),"RECREATE");
  f->cd();
  h_data->Write("h_Data",TObject::kOverwrite);
  h_tot->Write("h_TotBkg",TObject::kOverwrite);
  for(int i_smp=sample_name.size()-1;i_smp>0;i_smp--){
    h_bkg[i_smp]->Write("",TObject::kOverwrite);
  }
  h_signal->Write("h_ttH",TObject::kOverwrite);
  f->Close();
  f->~TFile();
  delete f;
}

TCanvas* TthPlot::GetCanvas(){
  return c;
}
