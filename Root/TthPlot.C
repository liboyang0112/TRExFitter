#include "TtHFitter/TthPlot.h"

#include "AtlasLabels.C"
#include "AtlasUtils.C"

using namespace std;

// TthPlot::TthPlot(){
//   fName = "c";
//   Init();
// }
TthPlot::TthPlot(string name,int canvasWidth,int canvasHeight){
  fName = name;
//   Init();
// }
// void TthPlot::Init(int canvasWidth,int canvasHeight){
//   c = new TCanvas(fName.c_str(),fName.c_str(),500,600);
  c = new TCanvas(fName.c_str(),fName.c_str(),canvasWidth,canvasHeight);
  //
  pad0 = new TPad("pad0","pad0",0,0.20,1,1,0,0,0);
  pad0->SetTickx(false);
  pad0->SetTicky(false);
  pad0->SetTopMargin(0.05);
  pad0->SetBottomMargin(0.1);
  pad0->SetLeftMargin(0.14);
  pad0->SetRightMargin(0.05);
  pad0->SetFrameBorderMode(0);
  pad0->SetFillStyle(0);
  //
  pad1 = new TPad("pad1","pad1",0,0,1,0.28,0,0,0);
  pad1->SetTickx(false);
  pad1->SetTicky(false);
  pad1->SetTopMargin(0.0);
  pad1->SetBottomMargin(0.37);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.05);
  pad1->SetFrameBorderMode(0);
  pad1->SetFillStyle(0);
  //
  pad1->Draw();
  pad0->Draw();
  pad0->cd();
  h_stack = new THStack("h_stack","h_stack");
  h_tot = 0x0;
  g_tot = 0x0;
  xtitle = "Variable [GeV]";
  ytitle = "Events";
//   leg  = new TLegend(0.54,0.93-2*(0.05),0.80,0.93);
//   leg1 = new TLegend(0.80,0.93-2*(0.05),0.94,0.93);
  leg_title = "e/#mu+jets";
  fLumi = "20.3 fb^{-1}";
  fCME = "8 TeV";
  fATLASlabel = "none";
  yMaxScale = 2.;
  Chi2prob = -1;
  KSprob = -1;
  //
  h_data = 0x0;
  g_data = 0x0;
  h_signal = 0x0;
  h_normsig = 0x0;
  //
  fIsNjet = false;
  //
  for(int i_bin=0;i_bin<MAXbins;i_bin++)
    fBinLabel[i_bin] = "";
}

void TthPlot::SetChannel(string name){
  leg_title = name;
}

void TthPlot::SetLumi(string name){
  fLumi = name;
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

void TthPlot::SetBinLabel(int bin,string name){
  fBinLabel[bin] = name;
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
  h_normsig = (TH1*)h->Clone(name.c_str());
//   sample_name.push_back(name);
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
  gStyle->SetEndErrorSize(4.);
  if(g_tot==0x0) g_tot = new TGraphAsymmErrors(h_tot);
  pad0->cd();
  TH1* h_dummy = (TH1*)h_tot->Clone("h_dummy");
  h_dummy->Scale(0);
  h_dummy->Draw("HIST");
  //
  if(options.find("log")!=string::npos) pad0->SetLogy();
  //
  bool hasData = true;
  if(h_data){
//     h_data->Draw("E");
//     h_data->SetMarkerSize(1.2);
    h_data->SetMarkerSize(1.4);
    h_data->SetLineWidth(2);
    // build asym data
    g_data = poissonize(h_data);
    g_data->SetMarkerSize(h_data->GetMarkerSize());
    g_data->SetMarkerColor(h_data->GetMarkerColor());
    g_data->SetMarkerStyle(h_data->GetMarkerStyle());
    g_data->SetLineWidth(h_data->GetLineWidth());
  }
  else{
    hasData = false;
    h_data = (TH1F*)h_tot->Clone("dummyData");
    h_data->SetTitle("Asimov Data");
//     h_data->SetLineWidth(0);
//     h_data->SetMarkerSize(0);
//     h_data->SetLineColor(kWhite);
//     h_data->SetFillStyle(0);
//     h_data->Draw("HIST same");
  }
  if(h_signal!=0x0) h_stack->Add(h_signal);
  h_stack->Draw("HIST same");
  //
  g_tot->SetFillStyle(3354);
  g_tot->SetFillColor(kBlue-7);
  g_tot->SetLineColor(kWhite);
  g_tot->SetLineWidth(0);
  g_tot->SetMarkerSize(0);
  g_tot->Draw("sameE2");
  //
  if(h_normsig!=0x0){
//     h_signal = (TH1*)h_normsig->Clone();
    h_normsig->Scale(h_tot->Integral()/h_normsig->Integral());
    h_normsig->SetLineColor(h_normsig->GetFillColor());
    h_normsig->SetFillColor(0);
    h_normsig->SetFillStyle(0);
    h_normsig->SetLineStyle(2);
    h_normsig->Draw("HISTsame");
  }
  //
//   if(hasData) h_data->Draw("sameE1");
  if(hasData) g_data->Draw("Ep1 same");

  h_dummy->GetXaxis()->SetTitle(xtitle.c_str());
  h_dummy->GetYaxis()->SetTitle(ytitle.c_str());
  if(fIsNjet){
    for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
      int nj = (int)h_dummy->GetXaxis()->GetBinCenter(i_bin);
      if(i_bin<h_dummy->GetNbinsX()) h_dummy->GetXaxis()->SetBinLabel( i_bin,Form("%d",nj) );
      else                           h_dummy->GetXaxis()->SetBinLabel( i_bin,Form("#geq%d",nj) );
    }
  }
  else{
    for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
      if(fBinLabel[i_bin]!="") h_dummy->GetXaxis()->SetBinLabel( i_bin, fBinLabel[i_bin].c_str());
    }
  }
  h_dummy->GetYaxis()->SetTitleOffset(2);
  if(options.find("log")==string::npos){
    h_dummy->SetMinimum(0);
    if(hasData) h_dummy->SetMaximum(yMaxScale*TMath::Max(h_tot->GetMaximum(),h_data->GetMaximum()+GC_up(h_data->GetMaximum())));
    else        h_dummy->SetMaximum(yMaxScale*h_tot->GetMaximum());
  }
  else{
    h_dummy->SetMaximum(h_tot->GetMaximum()*pow(10,yMaxScale));
    h_dummy->SetMinimum(1.);
  }
  pad0->RedrawAxis();
  
//   if(fATLASlabel!="none") ATLASLabel(0.18,0.85,(char*)fATLASlabel.c_str());
//   myText(0.18,0.8,1,Form("#sqrt{s} = %s, %s",fCME.c_str(),fLumi.c_str()));//,0.045);
//   myText(0.18,0.75,1,Form("%s",leg_title.c_str()));//,0.045);
  if(fATLASlabel!="none") ATLASLabel(0.18,0.85+0.04,(char*)fATLASlabel.c_str());
  myText(0.18,0.8+0.04,1,Form("#sqrt{s} = %s, %s",fCME.c_str(),fLumi.c_str()));//,0.045);
  myText(0.18,0.75+0.04,1,Form("%s",leg_title.c_str()));//,0.045);
  
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
//   bool drawRatio = true;
  TLine *hline;
  pad1->cd();
  pad1->GetFrame()->SetY1(2);
  //
  TH1* h_dummy2 = (TH1*)h_tot->Clone("h_dummy2");
  h_dummy2->Scale(0);
  h_dummy2->Draw("HIST");
  //
  TH1 *h_ratio;
  TGraphAsymmErrors *g_ratio;
//   if(h_data)
  h_ratio = (TH1*)h_data->Clone("h_ratio");
  //
  h_dummy2->GetYaxis()->SetTitleOffset(1.4*h_dummy->GetTitleOffset());
//   else {
//     h_ratio = (TH1*)h_tot->Clone("h_ratio");
//     h_ratio->Scale(0.01);
//   }
  TH1 *h_ratio2 = (TH1*)h_tot->Clone("h_ratio2");
  TGraphAsymmErrors *g_ratio2 = (TGraphAsymmErrors*)g_tot->Clone("g_ratio2");
  TH1 *h_tot_nosyst = (TH1*)h_tot->Clone("h_tot_nosyst");
  for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){
    h_tot_nosyst->SetBinError(i_bin,0);
  }
  h_dummy2->SetTitle("Data/MC");
  h_dummy2->GetYaxis()->CenterTitle();
  h_dummy2->GetYaxis()->SetTitle("Data / Pred.");
  h_dummy2->GetYaxis()->SetLabelSize(0.8*h_ratio->GetYaxis()->GetLabelSize());
  h_dummy2->GetYaxis()->SetLabelOffset(0.02);
  h_dummy2->GetYaxis()->SetNdivisions(504,false);
//   h_dummy2->GetYaxis()->SetNdivisions(303,true);
  gStyle->SetEndErrorSize(4.);
  h_ratio->Divide(h_tot_nosyst);
    h_ratio->SetMarkerStyle(24);
    h_ratio->SetMarkerSize(1.4);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineWidth(2);
  g_ratio = histToGraph(h_ratio);
  for(int i_bin=1;i_bin<=h_ratio->GetNbinsX();i_bin++){
    g_ratio->SetPointEYhigh( i_bin-1,g_data->GetErrorYhigh(i_bin-1)/h_tot->GetBinContent(i_bin) );
    g_ratio->SetPointEYlow(  i_bin-1,g_data->GetErrorYlow(i_bin-1) /h_tot->GetBinContent(i_bin) );
  }
  for(int i_bin=1;i_bin<h_tot_nosyst->GetNbinsX()+1;i_bin++){
    g_ratio2->SetPoint(i_bin-1,g_ratio2->GetX()[i_bin-1],g_ratio2->GetY()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
    g_ratio2->SetPointEXlow(i_bin-1,g_ratio2->GetEXlow()[i_bin-1]);
    g_ratio2->SetPointEXhigh(i_bin-1,g_ratio2->GetEXhigh()[i_bin-1]);
    g_ratio2->SetPointEYlow(i_bin-1,g_ratio2->GetEYlow()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
    g_ratio2->SetPointEYhigh(i_bin-1,g_ratio2->GetEYhigh()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
  }
  hline = new TLine(h_dummy2->GetXaxis()->GetXmin(),1,h_dummy2->GetXaxis()->GetXmax(),1);
  hline->SetLineColor(kRed);
  hline->SetLineWidth(2);
  hline->SetLineStyle(2);
//   h_ratio->Draw("0E1 same");
  g_ratio->Draw("Ep1 same");
  hline->Draw();
  h_dummy2->SetMinimum(0.);
  h_dummy2->SetMaximum(2.);
  if(options.find("prefit")!=string::npos){
    h_dummy2->SetMinimum(0.00);
    h_dummy2->SetMaximum(2.00);
  }
  else{
    h_dummy2->SetMinimum(0.50);
    h_dummy2->SetMaximum(1.50);
//     h_ratio->SetMinimum(0.501);
//     h_ratio->SetMaximum(1.499);
  }
  h_dummy2->GetXaxis()->SetTitle(h_dummy->GetXaxis()->GetTitle());
//   h_ratio->GetXaxis()->SetTitleSize(0.14);
  h_dummy2->GetXaxis()->SetTitleOffset(4.);
  h_dummy->GetXaxis()->SetTitle("");
//   h_ratio->GetXaxis()->SetLabelSize(0.14);
//   if(fIsNjet) h_ratio->GetXaxis()->SetLabelSize(0.2);
  h_dummy->GetXaxis()->SetLabelSize(0);
  for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
//     if(((string)h_dummy->GetXaxis()->GetBinLabel(i_bin))!="") h_ratio->GetXaxis()->SetBinLabel( i_bin, h_dummy->GetXaxis()->GetBinLabel(i_bin));
    if(((string)h_dummy->GetXaxis()->GetBinLabel(i_bin))!="") h_dummy2->GetXaxis()->SetBinLabel( i_bin, h_dummy->GetXaxis()->GetBinLabel(i_bin));
  }
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
//       if (isUp==1) arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),1.41, h_ratio->GetXaxis()->GetBinCenter(i_bin),1.45,0.030,"|>");
//       else         arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),0.54, h_ratio->GetXaxis()->GetBinCenter(i_bin),0.51,0.030,"|>");
//       arrow->SetFillColor(2);
//       arrow->SetFillStyle(1001);
//       arrow->SetLineColor(2);
      if (isUp==1) arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),1.45, h_ratio->GetXaxis()->GetBinCenter(i_bin),1.5,0.030,"|>");
      else         arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),0.55, h_ratio->GetXaxis()->GetBinCenter(i_bin),0.5,0.030,"|>");
      arrow->SetFillColor(10);
      arrow->SetFillStyle(1001);
      arrow->SetLineColor(kBlue-7);
      arrow->SetLineWidth(2);
      arrow->SetAngle(40);
//       arrow->Draw();
      //
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
        fix->SetLineColor(h_ratio->GetLineColor());
        fix->SetLineWidth(h_ratio->GetLineWidth());
        fix->Draw("SAME");
        fixUp->SetLineColor(1);
        fixUp->SetLineWidth(1);
        fixUp->Draw("SAME");
        fixDo->SetLineColor(1);
        fixDo->SetLineWidth(1);
        fixDo->Draw("SAME");
      }
      arrow->Draw();
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

// function to get asymmetric error bars for hists (Used in WZ observation)
double GC_up(double data) {
  if (data == 0 ) return 0;
  return 0.5*TMath::ChisquareQuantile(1.-0.1586555,2.*(data+1))-data;
}
double GC_down(double data) {
  if (data == 0 ) return 0;
  return data-0.5*TMath::ChisquareQuantile(0.1586555,2.*data);
}
TGraphAsymmErrors* poissonize(TH1 *h) {
  vector<int> points_to_remove;
  TGraphAsymmErrors* gr= new TGraphAsymmErrors(h);
  for (UInt_t i=0; i< (UInt_t)gr->GetN(); i++) {
    double content = (gr->GetY())[i];
//     gr->SetPointError(i,0.5*h->GetBinWidth(i),0.5*h->GetBinWidth(i),GC_down(content),GC_up(content));
    gr->SetPointError(i,0.499*h->GetBinWidth(i),0.5*h->GetBinWidth(i),GC_down(content),GC_up(content));
    if(content==0){
      gr->RemovePoint(i);
      i--;
    }
  }
  return gr;
}

TGraphAsymmErrors* histToGraph(TH1* h){
  TGraphAsymmErrors* gr= new TGraphAsymmErrors(h);
  for (UInt_t i=0; i< (UInt_t)gr->GetN(); i++) {
    gr->SetPointError(i,0.499*h->GetBinWidth(i),0.5*h->GetBinWidth(i),0,0);
  }
  gr->SetMarkerStyle(h->GetMarkerStyle());
  gr->SetMarkerSize(h->GetMarkerSize());
  gr->SetMarkerColor(h->GetMarkerColor());
  gr->SetLineWidth(h->GetLineWidth());
  gr->SetLineColor(h->GetLineColor());
  gr->SetLineStyle(h->GetLineStyle());
  return gr;
}
