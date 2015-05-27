#include "AtlasLabels.C"
#include "AtlasUtils.C"

void PlotLimits(){
  
  float xmax = 120;
  string process = "t#bar{t}t#bar{t} (SM)";
  
//   // Fit titles
//   vector<string> titles;
//   titles.push_back("Extended l+jets");
//   titles.push_back("CMS-like l+jets");
//   titles.push_back("Base l+jets");
//   // Fit names (the same as the .root file names!!)
//   vector<string> names;
//   names.push_back("New4t");
//   names.push_back("CMS4t");
//   names.push_back("IFAE4t");

  // ---
  
  // Fit titles
  vector<string> titles;
  titles.push_back("20 GeV & |#eta|<4.5");
  titles.push_back("|#eta|<4.5 jets");
  titles.push_back("20 GeV jets");
  titles.push_back("Default jets");
  // Fit names (the same as the .root file names!!)
  vector<string> names;
  names.push_back("Test_AllJets");
  names.push_back("Test_ForwJets");
  names.push_back("Test_LowPtJets");
  names.push_back("Test_Default");
//   xmax = 110;
  
  // ---
  
//   process = "t#bar{t}H(b#bar{b})";
//   xmax = 10;
//   // Fit titles
//   vector<string> titles;
//   titles.push_back("Default");
//   titles.push_back("20 GeV jets");
//   // Fit names (the same as the .root file names!!)
//   vector<string> names;
//   names.push_back("ttH_Default");
//   names.push_back("ttH_LowPtJets");

  // ---
  
  bool showObs = false;
  
  int N = names.size();
  
  float ymin = -0.5;
  float ymax = N-0.5;
//   float xmax = 120;
  
  TCanvas *c = new TCanvas("c","c",600,400);
  
  TGraphErrors *g_obs = new TGraphErrors(N);
  TGraphErrors *g_exp = new TGraphErrors(N);
  TGraphAsymmErrors *g_1s = new TGraphAsymmErrors(N);
  TGraphAsymmErrors *g_2s = new TGraphAsymmErrors(N);
  
  TH1F* h_dummy = new TH1F("h_dummy","h_dummy",1,0,xmax);
  h_dummy->Draw();
  h_dummy->SetMinimum(ymin);
  h_dummy->SetMaximum(ymax);
  h_dummy->SetLineColor(kWhite);
  
  int Ndiv = N+1;
  
  h_dummy->GetYaxis()->Set(N,ymin,ymax);
  h_dummy->GetYaxis()->SetNdivisions(Ndiv);
  
  TFile *f;
  TH1* h;
  
  // get values
  for(int i=0;i<N;i++){
    h_dummy->GetYaxis()->SetBinLabel(i+1,titles[i].c_str());
    
//     f = new TFile(Form("root-files/limits/%s.root",names[i].c_str()) );
    f = new TFile(Form("%s/Limits/%s.root",names[i].c_str(),names[i].c_str()) );
    h = (TH1*)f->Get("limit");
    
    g_obs->SetPoint(i,h->GetBinContent(1),i);
    g_exp->SetPoint(i,h->GetBinContent(2),i);
    g_1s->SetPoint(i,h->GetBinContent(2),i);
    g_2s->SetPoint(i,h->GetBinContent(2),i);
    g_obs->SetPointError(i,0,0.5);
    g_exp->SetPointError(i,0,0.5);
    g_1s->SetPointError(i,h->GetBinContent(2)-h->GetBinContent(5),h->GetBinContent(4)-h->GetBinContent(2),0.5,0.5);
    g_2s->SetPointError(i,h->GetBinContent(2)-h->GetBinContent(6),h->GetBinContent(3)-h->GetBinContent(2),0.5,0.5);
  }
    
  g_obs->SetLineWidth(2);
  g_exp->SetLineWidth(2);
  g_exp->SetLineStyle(2);
  g_1s->SetFillColor(kGreen);
  g_1s->SetLineWidth(2);
  g_1s->SetLineStyle(2);
  g_2s->SetFillColor(kYellow);
  g_2s->SetLineColor(kYellow);
  g_2s->SetLineWidth(0);
  
  g_2s->SetMarkerSize(0);
  g_1s->SetMarkerSize(0);
  g_exp->SetMarkerSize(0);
  g_obs->SetMarkerSize(0);
  
  g_2s->Draw("E2 same");
  g_1s->Draw("E2 same");
  g_exp->Draw("E same");
  if(showObs) g_obs->Draw("E same");

  TLine *l_SM = new TLine(1,-0.5,1,N-0.5);
  l_SM->SetLineWidth(2);
  l_SM->SetLineColor(kGray);
  l_SM->Draw("same");
  
  c->RedrawAxis();

  gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
  gPad->SetBottomMargin( 1.15*gPad->GetBottomMargin() );
  gPad->SetTopMargin( 1.8*gPad->GetTopMargin() );
  h_dummy->GetXaxis()->SetTitle("95% CL limit on #sigma/#sigma_{SM}");

  ATLASLabel(0.2,0.93,"    Internal",kBlack);
  myText(0.50,0.93,kBlack,process.c_str());
  myText(0.65,0.93,kBlack,"#sqrt{s} = 8 TeV, 20.3 fb^{-1}");
  
  TLegend *leg;
  if(showObs) leg = new TLegend(0.60,0.2,0.95,0.40);
  else        leg = new TLegend(0.60,0.2,0.95,0.35);
  leg->SetTextSize(gStyle->GetTextSize());
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g_1s,"Expected #pm 1 #sigma","lf");
  leg->AddEntry(g_2s,"Expected #pm 2 #sigma","f");
  if(showObs) leg->AddEntry(g_obs,"Observed","l");
  leg->Draw();
  
  myText(0.75,0.4,kBlack,"Stat. only");
  
  c->SaveAs("Limits.png");
}
