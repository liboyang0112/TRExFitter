// C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <map>

// Root
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"

void PlotRanking(){
  
    map<string,string> systMap;
    systMap[""] = "";
    
    systMap["BTag_B_NP1"]        = "b-tag 1";
    systMap["BTag_B_NP2"]        = "b-tag 2";
    systMap["BTag_B_NP3"]        = "b-tag 3";
    systMap["JER"]               = "JER";
    systMap["JES_Scenario1_NP1"] = "JES 1";
    systMap["JES_Scenario1_NP2"] = "JES 2";
    systMap["JES_Scenario1_NP3"] = "JES 3";
    systMap["stXsec"]            = "single top x-sec";
    systMap["ttXsec"]            = "t#bar{t} x-sec";
    systMap["tt_Rad_HiLow"]      = "t#bar{t} Radiation";
    systMap["tt_Shower"]         = "t#bar{t} Parton Shower";
    
    ifstream fin("NPR_mu.txt");
    string paramname;
    double nuiphat;
    double nuiperrhi;
    double nuiperrlo;
    double PoiUp;
    double PoiDown;
    double PoiNomUp;
    double PoiNomDown;
    vector<string> parname;
    vector<double> nuhat;
    vector<double> nuerrhi;
    vector<double> nuerrlo;
    vector<double> poiup;
    vector<double> poidown;
    vector<double> poinomup;
    vector<double> poinomdown;
    vector<double> number;

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
    cout<<"Beginning of NP ordering "<<endl;
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
    cout<<"End of ordering"<<endl;
    number.push_back(parname.size()-0.5);

    double poimax = 0;
    for (int i=0;i<SIZE;i++) {
        poimax=TMath::Max(poimax,TMath::Max( TMath::Abs(poiup[i]),TMath::Abs(poidown[i]) ));
        cout << poiup[i] << " " << poimax << endl;
        nuerrlo[i] = TMath::Abs(nuerrlo[i]);
    }
    poimax *= 1.2;

    for (int i=0;i<SIZE;i++) {
        poiup[i]     *= (2./poimax);
        poidown[i]   *= (2./poimax);
        poinomup[i]  *= (2./poimax);
        poinomdown[i]*= (2./poimax);
    }
  

  
    // Graphical part - rewritten taking DrawPulls in TtHFitter
    
    float lineHeight = 20;
    float offsetUp = 60; // external
    float offsetDown = 60;
    float offsetUp1 = 100; // internal
    float offsetDown1 = 10;
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
    
//     NuisParameter *par;
    int idx = 0;
    std::vector< string > Names;
    Names.clear();
    string parTitle;
    
    for(unsigned int i = 0; i<parname.size(); ++i){
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

        idx = i;
        
        g->SetPoint(      idx, nuhat[idx],idx+0.5);
        g->SetPointEXhigh(idx, nuerrhi[idx]);
        g->SetPointEXlow( idx, nuerrlo[idx]);
        
        g1->SetPoint(      idx, 0.,idx+0.5);
        g1->SetPointEXhigh(idx, poiup[idx]);
        g1->SetPointEXlow( idx, 0.);
        g1->SetPointEYhigh(idx, 0.4);
        g1->SetPointEYlow( idx, 0.4);
        
        g2->SetPoint(      idx, 0.,idx+0.5);
        g2->SetPointEXhigh(idx, poidown[idx]);
        g2->SetPointEXlow( idx, 0.);
        g2->SetPointEYhigh(idx, 0.4);
        g2->SetPointEYlow( idx, 0.4);
        
//         Poidown, Poiup
//         Poinomdown, Poinomup
        
        parTitle = systMap[parname[i]];
//         h2->GetYaxis()->SetBinLabel(idx+1,parTitle.c_str());
        
//         Names.push_back(par->fName);
        Names.push_back(parTitle);
    
        idx ++;
        if(idx > max)  max = idx;      
    }

    TCanvas *c = new TCanvas("c","c",600,newHeight);
    c->SetTicks(0,0);
    gPad->SetLeftMargin(0.33);
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

//     TLine l0;
//     TBox b1, b2;
//     if(brazilian){
//         l0 = TLine(0,0,0,max);
//         l0.SetLineStyle(7);
//         l0.SetLineColor(kBlack);
//         b1 = TBox(-1,0,1,max);
//         b2 = TBox(-2,0,2,max);
//         b1.SetFillColor(kGreen);
//         b2.SetFillColor(kYellow);
//         b2.Draw("same");
//         b1.Draw("same");
//         l0.Draw("same");
//     }
    
//     g1->SetFillColor(kAzure-5);
//     g2->SetFillColor(kAzure+6);
    g1->SetFillColor(kAzure-4);
    g2->SetFillColor(kCyan);
    g1->SetLineColor(g1->GetFillColor());
    g2->SetLineColor(g2->GetFillColor());
    
    g->SetLineWidth(2);
    
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
    axis_up->SetLabelOffset( 0. );
    axis_up->SetLabelSize(   h_dummy->GetXaxis()->GetLabelSize() );
    axis_up->SetLabelFont(   gStyle->GetTextFont() );
    axis_up->Draw();
    axis_up->CenterTitle();
    axis_up->SetTitle("#Delta#mu");
    
    TLegend *leg = new TLegend(0.02,0.7,0.3,0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.2);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->AddEntry(g,"(#hat{#theta}-#theta_{0})/#Delta#theta","lep");
    leg->AddEntry(g1,"#Delta#mu for #theta_{0}=+#Delta#hat{#theta}","f");
    leg->AddEntry(g2,"#Delta#mu for #theta_{0}=-#Delta#hat{#theta}","f");
    leg->Draw();
    
    gPad->RedrawAxis();
    
//     c->SaveAs((path+"/NuisPar.png").c_str());
//     c->SaveAs(path.c_str());
//     delete c;
    
    
    // END
  
  
  
  
//     // Start graphical part
//     TCanvas* c1 = new TCanvas("c1","c1",1024,1448);
//     TPad *pad1 = new TPad("pad1", "pad1", 0  , 0  , 1.0 , 1.0  , 0);
//     pad1->SetLeftMargin(0.40);
//     pad1->SetRightMargin(0.05);
//     pad1->SetBottomMargin(0.09);
//     pad1->SetTopMargin(0.10);
//     pad1->Draw();
//     pad1->cd();
//     pad1->SetTicks(0,1);
//     float markerSize = 1.;
//     
// //     TGraphAsymmErrors* gr2 = new TGraphAsymmErrors(parname.size(), Nuhat, Number, One, One, NULL, NULL);
// //     gr2->SetLineColor(kRed);
// //     gr2->SetMarkerColor(kBlack);
// //     gr2->SetMarkerStyle(20);
// //     gr2->SetLineStyle(1);
// //     gr2->SetLineWidth(2);
// //     gr2->SetMarkerSize(markerSize);
// //     double xlimup=TMath::Max(gr2->GetXaxis()->GetXmax()+0.1,borderHi/poimax+0.1);
// //     double xlimdown=TMath::Min(gr2->GetXaxis()->GetXmin()-0.1,borderLo/poimax-0.1);
// //     gr2->GetXaxis()->SetLimits(xlimdown,xlimup);
// //     gr2->GetHistogram()->SetMaximum(SIZE+0.15);
// //     gr2->GetHistogram()->SetMinimum(0.);
// //     gr2->GetYaxis()->SetBinLabel(1,"");
// //     gr2->GetYaxis()->SetTickLength(0);
//  
//     TGraphAsymmErrors* gr = new TGraphAsymmErrors(parname.size(), Nuhat, Number, Nuerrlo, Nuerrhi, NULL, NULL);
//     gr->SetLineColor(kBlack);
//     gr->SetMarkerColor(kBlack);
//     gr->SetMarkerStyle(20);
//     gr->SetLineStyle(1);
//     gr->SetLineWidth(2);
//     gr->SetMarkerSize(markerSize);
//     gr->GetXaxis()->SetLimits(xlimdown,xlimup);
//     gr->GetHistogram()->SetMaximum(SIZE+0.15);
//     gr->GetHistogram()->SetMinimum(0.);
//     gr->GetYaxis()->SetTickLength(0);
//  
//     TGraphAsymmErrors* gr3 = new TGraphAsymmErrors(parname.size(), Zero, Number, Poidown, Poiup, Height, Height);
//     gr3->SetFillColor(kBlue);
//     gr3->SetFillStyle(3354);
//     gr3->GetXaxis()->SetLimits(xlimdown,xlimup);
//     gr3->GetHistogram()->SetMaximum(SIZE+0.15);
//     gr3->GetHistogram()->SetMinimum(0.);
//     gr3->GetYaxis()->SetTickLength(0);
//     gr3->SetMarkerSize(0);
// 
//     TGraphAsymmErrors* gr4 = new TGraphAsymmErrors(parname.size(), Zero, Number, Poinomdown, Poinomup, Height, Height);
//     gr4->SetFillColor(kYellow-7);
//     gr4->SetMarkerSize(0);
//     gr4->GetXaxis()->SetLimits(xlimdown,xlimup);
//     gr4->GetHistogram()->SetMaximum(SIZE+0.15);
//     gr4->GetHistogram()->SetMinimum(0.);
//     gr4->GetYaxis()->SetTickLength(0.);
//     gr4->GetYaxis()->SetBinLabel(1,"");
//  
//     TLine* l0=new TLine(0,0.,0,parname.size()+0.1);
//     l0->SetLineStyle(2);
//     l0->SetLineColor(40);
//     TLine* l1=new TLine(1,0.,1,parname.size()+0.1);
//     l1->SetLineStyle(2);
//     l1->SetLineColor(40);
//     TLine* l2=new TLine(-1,0.,-1,parname.size()+0.1);
//     l2->SetLineStyle(2);
//     l2->SetLineColor(40);
// 
//     gr4->Draw("P2A");
//     gr3->Draw("P2same");
//     gr2->Draw("P");
//     gr->Draw("p");
//     l0->Draw();
//     l1->Draw();
//     l2->Draw();
// 
//     TH2F *h2 = new TH2F("h2", "", 1,xlimdown, xlimup, parname.size(),0+0.05,parname.size()+0.05);
//     string parTitle;
//     for (unsigned int i=0;i<parname.size();i++){
// //       h2->GetYaxis()->SetBinLabel(i+1,parname[i].c_str());
//       parTitle = systMap[parname[i]];
//       h2->GetYaxis()->SetBinLabel(i+1,parTitle.c_str());
//     }
//     h2->GetYaxis()->SetTickLength(0);
//     axis2->ImportAxisAttributes(h2->GetYaxis());
//     axis2->SetTickSize(0.);
//     axis2->Draw();
//     TGaxis *axis=new TGaxis(xlimdown,SIZE+0.15,xlimup,SIZE+0.15,xlimdown*poimax,xlimup*poimax,510,"-");
//     axis->SetLabelColor(kBlack);
//     axis->Draw();
// 
//     TLegend* leg = new TLegend(0.65,0.1,0.95,0.29);
//     leg->AddEntry(gr, "Pull","lp");
//     leg->AddEntry(gr2, "1 standard deviation","l");
//     leg->AddEntry(gr4, "Prefit Impact on #hat#mu", "f");
//     leg->AddEntry(gr3, "Postfit Impact on #hat#mu", "f");    
//     leg->Draw();

//     MainDirRank->cd();
//     c1->Write();
    //c1->Print(dirName+"/NPRanking.eps");
    //c1->Print(dirName+"/NPRanking.png");
//     c->Print("NPRanking.eps");
    c->Print("NPRanking.png");
//     c1->Close();
//     gROOT->cd();
//     cout<<"             "<<endl;
//     cout<<"             "<<endl;
//     swatch.Print();
}
