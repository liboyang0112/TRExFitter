// Class include
#include "TRExFitter/TRExPlot.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/TRExFit.h"

// ATLAS stuff
#include "AtlasUtils/AtlasStyle.h"
#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

// ROOT includes
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

// c++ includes
#include <algorithm>
#include <iostream>

using namespace std;

//_____________________________________________________________________________
//
TRExPlot::TRExPlot(std::string name,int canvasWidth,int canvasHeight,bool hideRatioPad) :
    fName(name),
    h_data(nullptr),
    g_data(nullptr),
    h_stack(new THStack("h_stack","h_stack")),
    h_tot(nullptr),
    g_tot(nullptr),
    h_tot_bkg_prefit(nullptr),
    h_dummy(nullptr),
    c(new TCanvas(fName.c_str(),fName.c_str(),canvasWidth,canvasHeight)),
    leg(nullptr),
    leg1(nullptr),
    pad0(nullptr),
    pad1(nullptr),
    xtitle("Variable [GeV]"),
    ytitle("Events"),
    fDataName("Data"),
    fLumi("XXX fb^{-1}"),
    fCME("13 TeV"),
    fATLASlabel("none"),
    yMaxScale(2.),
    NDF(-1),
    Chi2val(-1),
    Chi2prob(-1),
    KSprob(-1),
    fYmax(0),
    fYmin(0),
    fRatioYmax(2.),
    fRatioYmin(0.),
    fBinWidth(-1),
    fIsNjet(false),
    fShowYields(false),
    fLumiScale(1.),
    fLegendNColumns(2),
    fRatioYtitle(""),
    fRatioType("DATA/MC"),
    fLabelX(-1),
    fLabelY(-1),
    fLegendX1(-1),
    fLegendX2(-1),
    fLegendY(-1),
    fBlinding(nullptr) {

    if(hideRatioPad){
        pad0 = new TPad("pad0","pad0",0,0,1,1,0,0,0);
    }
    else{
        pad0 = new TPad("pad0","pad0",0,0.20,1,1,0,0,0);
    }
    pad0->SetTicks(1,1);
    pad0->SetTopMargin(0.05*(700./canvasHeight));
    if(hideRatioPad){
        pad0->SetBottomMargin(0.14*(600./canvasHeight));
    }
    else{
        pad0->SetBottomMargin(0.1);
    }
    pad0->SetLeftMargin(0.14*(600./canvasWidth));
    pad0->SetRightMargin(0.05*(600./canvasWidth));
    pad0->SetFrameBorderMode(0);
    pad0->SetFillStyle(0);
    //
    if(hideRatioPad){
        pad1 = nullptr;
    }
    else{
        pad1 = new TPad("pad1","pad1",0,0,1,0.28,0,0,0);
        pad1->SetTicks(1,1);
        pad1->SetTopMargin(0.0);
        pad1->SetBottomMargin(0.37*(700./canvasHeight));
        pad1->SetLeftMargin(0.14*(600./canvasWidth));
        pad1->SetRightMargin(0.05*(600./canvasWidth));
        pad1->SetFrameBorderMode(0);
        pad1->SetFillStyle(0);
    }
    if(pad1!=nullptr) pad1->Draw();
    pad0->Draw();
    pad0->cd();
}

//_____________________________________________________________________________
//
TRExPlot::~TRExPlot(){
    delete h_stack;
    delete h_tot_bkg_prefit;
    delete h_dummy;
    delete leg;
    delete leg1;
    delete pad0;
    delete pad1;
    delete c;
}

//_____________________________________________________________________________
//
void TRExPlot::SetChannel(const std::string& name){
    fLabels.clear();
    fLabels.push_back(name);
}

//_____________________________________________________________________________
//
void TRExPlot::AddLabel(const std::string& name){
    fLabels.push_back(name);
}

//_____________________________________________________________________________
//
void TRExPlot::SetLumi(const std::string& name){
    fLumi = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetLumiScale(double scale){
    fLumiScale = scale;
}

//_____________________________________________________________________________
//
void TRExPlot::SetCME(const std::string& name){
    fCME = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetXaxis(const std::string& name,bool isNjet){
    xtitle = name;
    fIsNjet = isNjet;
}

//_____________________________________________________________________________
//
void TRExPlot::SetYaxis(const std::string& name){
    ytitle = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetYmaxScale(double scale){
    yMaxScale = scale;
}

//_____________________________________________________________________________
//
void TRExPlot::ResizeBinLabel(const int n) {
    fBinLabel.resize(n);
}

//_____________________________________________________________________________
//
void TRExPlot::SetBinLabel(int bin, const std::string& name){
    fBinLabel[bin] = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetBinWidth(double width){
    fBinWidth = width;
}

//_____________________________________________________________________________
//
void TRExPlot::SetData(TH1* h,std::string name){
    h_data.reset(static_cast<TH1*>(h->Clone()));
    h_data->SetDirectory(nullptr);
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    fDataName = name;
}

//_____________________________________________________________________________
//
void TRExPlot::AddSignal(TH1* h,std::string name){
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    unsigned int idx = std::find(fSigNames.begin(),fSigNames.end(),name) - fSigNames.begin();
    if(idx<fSigNames.size()){
        h_signal[idx]->Add(h,fLumiScale);
    }
    else{
        h_signal.push_back((TH1*)h->Clone());
        h_signal[fSigNames.size()]->Scale(fLumiScale);
        fSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::AddNormSignal(TH1* h,std::string name){
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    unsigned int idx = std::find(fNormSigNames.begin(),fNormSigNames.end(),name) - fNormSigNames.begin();
    if(idx<fNormSigNames.size()){
        h_normsig[idx]->Add(h,fLumiScale);
    }
    else{
        h_normsig.push_back((TH1*)h->Clone());
        h_normsig[fNormSigNames.size()]->Scale(fLumiScale);
        fNormSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::AddOverSignal(TH1* h,std::string name){
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    unsigned int idx = std::find(fOverSigNames.begin(),fOverSigNames.end(),name) - fOverSigNames.begin();
    if(idx<fOverSigNames.size()){
        h_oversig[idx]->Add(h,fLumiScale);
    }
    else{
        h_oversig.push_back((TH1*)h->Clone());
        h_oversig[fOverSigNames.size()]->Scale(fLumiScale);
        fOverSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::AddBackground(TH1* h,std::string name){
    if(h_tot==nullptr) {
        h_tot.reset(static_cast<TH1*>(h->Clone()));
        h_tot->SetDirectory(nullptr);
    }
    else h_tot->Add(h);
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    //
    unsigned int idx = std::find(fBkgNames.begin(),fBkgNames.end(),name) - fBkgNames.begin();
    if(idx<fBkgNames.size()){
        h_bkg[idx]->Add(h,fLumiScale);
    }
    else{
        h_bkg.push_back((TH1*)h->Clone());
        h_bkg[fBkgNames.size()]->Scale(fLumiScale);
        fBkgNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SetTotBkg(TH1* h){
    h_tot.reset(static_cast<TH1*>(h->Clone()));
    h_tot->Scale(fLumiScale);
    h_tot->SetDirectory(nullptr);
    g_tot = std::make_unique<TGraphAsymmErrors>(h);
    for(int i=0;i<g_tot->GetN();i++){
        g_tot->GetY()[i]      *= fLumiScale;
        g_tot->GetEYlow()[i]  *= fLumiScale;
        g_tot->GetEYhigh()[i] *= fLumiScale;
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SetTotBkgAsym(TGraphAsymmErrors* g){
    g_tot.reset(static_cast<TGraphAsymmErrors*>(g->Clone()));
    for(int i=0;i<g_tot->GetN();i++){
        g_tot->GetY()[i] *= fLumiScale;
        g_tot->GetEYlow()[i]  *= fLumiScale;
        g_tot->GetEYhigh()[i] *= fLumiScale;
    }
    for(int i=1;i<h_tot->GetNbinsX()+1;i++){
        h_tot->SetBinContent(i,g_tot->GetY()[i-1]);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SetChi2KS(double chi2prob,double ksprob,double chi2val,int ndf){
    Chi2prob = chi2prob;
    KSprob = ksprob;
    Chi2val = chi2val;
    NDF = ndf;
}

//_____________________________________________________________________________
//
void TRExPlot::BlindData(){
    //
    // Eventually blind bins
    
    if(h_data && ((fSigNames.size() > 0) || (fOverSigNames.size() > 0)) && h_tot) {
        if(fBlinding) {
            Common::BlindDataHisto(h_data.get(),fBlinding.get());
        } else{
            fBlinding = Common::BlindDataHisto(h_data.get(), fBlindedBins);
            // if more than one signal:
            if(fSigNames.size() > 1) {
                for(std::size_t i_sig = 1; i_sig < fSigNames.size(); ++i_sig) {
                    auto tmp = Common::BlindDataHisto(h_data.get(), fBlindedBins);
                    fBlinding->Add(tmp.get());
                    fBlinding->Scale(2.);
                }
            }
            // now blind if oversig is used
            if(fOverSigNames.size() > 0) {
                for(std::size_t i_sig = 0; i_sig < fOverSigNames.size(); ++i_sig) {
                    auto tmp = Common::BlindDataHisto(h_data.get(), fBlindedBins);
                    fBlinding->Add(tmp.get());
                    fBlinding->Scale(2.);
                }
            }
        }
    } else {
        WriteWarningStatus("TRExPlot::BlindData", "Either h_data, h_signal, h_tot not defined.");
        WriteWarningStatus("TRExPlot::BlindData", "Blinding not possible. Skipped.");
    }
}

//_____________________________________________________________________________
//
TH1* TRExPlot::GetTotBkg() const{
    TH1* h = (TH1*)h_tot->Clone("h_tot_bkg");
    for (unsigned int i=0; i<fSigNames.size(); ++i) {
      h->Add( h_signal[i], -1);
    }
    h->SetDirectory(nullptr);
    return h;
}

//_____________________________________________________________________________
//
void TRExPlot::Draw(std::string options){

    /////////////////////////
    //
    // Main function of the class
    // ==========================
    //   It takes the data, background, signal to perform the full comparison (stack, ratio plot, ...)
    //
    /////////////////////////

    //
    // Draws an empty histogram to reserve the upper pad and set style
    //
    gStyle->SetEndErrorSize(0);
    pad0->cd();
    h_dummy = (TH1*)h_tot->Clone("h_dummy");
    h_dummy->SetDirectory(nullptr);
    h_dummy->Scale(0);
    if(pad0->GetWw() > pad0->GetWh()){
        h_dummy->GetYaxis()->SetTickLength(0.01);
        h_dummy->GetXaxis()->SetTickLength(0.02);
    }
    if (fXaxisRange.size() > 1){
        h_dummy->GetXaxis()->SetRangeUser(fXaxisRange.at(0), fXaxisRange.at(1));
    }
    h_dummy->Draw("HIST");
    if(options.find("log")!=std::string::npos) pad0->SetLogy();

    if(g_tot==nullptr) g_tot = std::make_unique<TGraphAsymmErrors>(h_tot.get());

    //
    // Determines if the data is real (and computes the poisson uncertainty) or not
    //
    bool hasData = true;
    if(h_data){
        h_data->SetMarkerSize(1.4);
        h_data->SetLineWidth(2);
        // build asym data
        if(options.find("poiss")!=std::string::npos) g_data = poissonize(h_data.get());
        else                                         g_data = histToGraph(h_data.get());
    }
    else{
        hasData = false;
        h_data.reset(static_cast<TH1D*>(h_tot->Clone("dummyData")));//tajes data = total
        h_data->SetTitle("Asimov Data");
        h_data->SetDirectory(nullptr);
        g_data = std::make_unique<TGraphAsymmErrors>(h_data.get());
    }

    //
    // Add Bkg's to the stack
    //
    for(int i_smp=fBkgNames.size()-1;i_smp>=0;i_smp--){
        h_bkg[i_smp]->SetLineWidth(1);
        h_stack->Add(h_bkg[i_smp]);
    }

    //
    // Eventually add Signal(s)
    //
    for(int i_smp=fSigNames.size()-1;i_smp>=0;i_smp--){
        h_signal[i_smp]->SetLineWidth(1);
        h_stack->Add(h_signal[i_smp]);
    }

    //
    // Draw
    //
    h_stack->Draw("HIST same");

    if( TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit ) {
      h_tot_bkg_prefit->SetFillColor(0);
      h_tot_bkg_prefit->SetLineStyle(kDashed);
      h_tot_bkg_prefit->SetLineColor(kBlue);
      h_tot_bkg_prefit->SetLineWidth(2);
      h_tot_bkg_prefit->Draw("HIST same");
    }

    //
    // Total error bands style setting
    //
    g_tot->SetFillStyle(3354);
    if(TRExFitter::OPTION["SystFillStyle"]!=0) g_tot->SetFillStyle(TRExFitter::OPTION["SystFillStyle"]);
    g_tot->SetFillColor(kBlue-7);
    g_tot->SetLineColor(kWhite);
    g_tot->SetLineWidth(0);
    g_tot->SetMarkerSize(0);
    g_tot->Draw("sameE2");

    //
    // Draw a normalized signal distribution
    //
    std::vector<double> signalScale;
    signalScale.resize(fNormSigNames.size());
    for(int i_smp=fNormSigNames.size()-1;i_smp>=0;i_smp--){
        if (std::fabs(h_normsig[i_smp]->Integral()) < 1e-10) {
            // division by zero
            WriteWarningStatus("TRExPlot::Draw", " --- Signal " + fNormSigNames[i_smp] + " has integral equal to zero - cannot scale, returning scale = 1");
            signalScale[i_smp] = 1;
        } else {
            signalScale[i_smp] = h_tot->Integral()/h_normsig[i_smp]->Integral();
            WriteInfoStatus("TRExPlot::Draw", "--- Signal " + fNormSigNames[i_smp] + " scaled by " + std::to_string(signalScale[i_smp]));
        }
        h_normsig[i_smp]->Scale(signalScale[i_smp]);
        h_normsig[i_smp]->SetLineColor(h_normsig[i_smp]->GetFillColor());
        h_normsig[i_smp]->SetFillColor(0);
        h_normsig[i_smp]->SetFillStyle(0);
        h_normsig[i_smp]->SetLineStyle(2);
        h_normsig[i_smp]->SetLineWidth(2);
        h_normsig[i_smp]->Draw("HISTsame");
    }

    //
    // Draw a overlayed signal distribution
    //
    //if(h_oversig!=nullptr){
        for(int i_smp=fOverSigNames.size()-1;i_smp>=0;i_smp--){
            h_oversig[i_smp]->SetLineColor(h_oversig[i_smp]->GetFillColor());
            h_oversig[i_smp]->SetFillColor(0);
            h_oversig[i_smp]->SetFillStyle(0);
            h_oversig[i_smp]->SetLineStyle(2);
            h_oversig[i_smp]->SetLineWidth(2);
            h_oversig[i_smp]->Draw("HISTsame");
        }
    //}

    //
    // Draw data (if it is real data of course)
    //
    if(hasData) g_data->Draw("Ep1 same");

    //
    // Draw blinding markers
    //
    TH1D* h_blind = nullptr;
    if(fBlinding) {
        h_blind = static_cast<TH1D*>(fBlinding->Clone("h_blind"));
        h_blind->SetLineWidth(0);
        h_blind->SetLineColor(kGray);
        h_blind->SetFillColor(kGray);
        h_blind->SetFillStyle(3345);
        h_blind->Draw("same HIST");
    }

    //
    // Axes labelling and style
    //
    h_dummy->GetXaxis()->SetTitle(xtitle.c_str());
    h_dummy->GetYaxis()->SetTitle(ytitle.c_str());
    if(fIsNjet){
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            int nj = (int)h_dummy->GetXaxis()->GetBinCenter(i_bin);
            if(i_bin<h_dummy->GetNbinsX()) h_dummy->GetXaxis()->SetBinLabel( i_bin,Form("%d",nj) );
            else                           h_dummy->GetXaxis()->SetBinLabel( i_bin,Form("#geq%d",nj) );
        }
    }
    else if ((int)fBinLabel.size() > h_dummy->GetNbinsX()) {
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            if(fBinLabel[i_bin]!="") h_dummy->GetXaxis()->SetBinLabel( i_bin, fBinLabel[i_bin].c_str());
        }
    }
    if(fBinLabel.size() > 1 && fBinLabel[1]!="") h_dummy->GetXaxis()->LabelsOption("d");
    double offset = 2.*pad0->GetWh()/pad0->GetWw();
    h_dummy->GetYaxis()->SetTitleOffset( offset );
    h_dummy->GetXaxis()->SetTitleOffset( 2 );

    //
    // Fix / redraw axis
    //
    pad0->RedrawAxis();

    double textHeight = 0.05*(672./pad0->GetWh());
    if(pad1==nullptr) textHeight *= 0.8;

    //
    // ATLAS labels
    //
    double labelX = 0.18*(600./pad0->GetWw());
    if(fLabelX>=0){
        labelX = fLabelX;
    }
    // was 0.84-textHeight+0.04
    double labelY = 1-0.08*(700./c->GetWh());
    if(fLabelY>=0){
        labelY = fLabelY;
    }
    // scale it down to give space to text (in this way one can set the same value for LabelY and LegendY and have them vertically alligned)
    labelY -= textHeight - 0.015;

    if(fATLASlabel!="none") ATLASLabel(labelX,labelY,(char*)fATLASlabel.c_str());
    myText(labelX,labelY-textHeight,1,Form("#sqrt{s} = %s, %s",fCME.c_str(),fLumi.c_str()));//,0.045);
    for(unsigned int i_lab=0;i_lab<fLabels.size();i_lab++){
        myText(labelX,labelY-(i_lab+2)*textHeight,1,Form("%s",fLabels[i_lab].c_str()));//,0.045);
    }

    double legX1 = 1-fLegendNColumns*0.2*(600./pad0->GetWw())-0.1*(600./pad0->GetWw());
    if(fLegendX1>=0){
        legX1 = fLegendX1;
    }
    double legX2 = 1-0.1*(600./pad0->GetWw());
    if(fLegendX2>=0){
        legX2 = fLegendX2;
    }
    double legXmid = legX1+0.5*(legX2-legX1);

    double legY = 1-0.08*(700./c->GetWh());
    if(fLegendY>=0){
        legY = fLegendY;
    }

    if(fShowYields){
        legXmid = legX1+0.6*(legX2-legX1);
        leg  = new TLegend(legX1,legY-(fBkgNames.size()+fSigNames.size()+2+(hasData))*textHeight, legXmid,legY);
        leg1 = new TLegend(legXmid,leg->GetY1(), legX2,leg->GetY2());
        //
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        if(!TRExFitter::LEGENDLEFT) leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetMargin(0.22);
        leg1->SetFillStyle(0);
        leg1->SetBorderSize(0);
        leg1->SetTextAlign(32);
        leg1->SetTextFont(gStyle->GetTextFont());
        leg1->SetTextSize(gStyle->GetTextSize());
        leg1->SetMargin(0.);

        // Add data in the legend if real data are here
        if(hasData){
            leg->AddEntry(h_data.get(),fDataName.c_str(),"lep");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_data->Integral()),"");
        }

        // Signal and background legends
        for(unsigned int i_smp=0;i_smp<fSigNames.size();i_smp++){
            leg->AddEntry(h_signal[i_smp], fSigNames[i_smp].c_str(),"f");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_signal[i_smp]->Integral()),"");
        }
        for(unsigned int i_smp=0;i_smp<fBkgNames.size();i_smp++){
            leg->AddEntry(h_bkg[i_smp], fBkgNames[i_smp].c_str(),"f");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_bkg[i_smp]->Integral()),"");
        }
        leg->AddEntry((TObject*)0,"Total","");
        leg1->AddEntry((TObject*)0,Form("%.1f",h_tot->Integral()),"");
        if(TRExFitter::OPTION["TRExbbStyle"]!=0) leg->AddEntry(g_tot.get(),"Total unc.","f");
        else leg->AddEntry(g_tot.get(),"Uncertainty","f");
        leg1->AddEntry((TObject*)0," ","");

        if(TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit) {
            leg->AddEntry(h_tot_bkg_prefit,"Pre-Fit Bkgd.","l");
            leg1->AddEntry((TObject*)0," ","");
        }

        leg->Draw();
        leg1->Draw();
    }
    else{
        int Nrows = fBkgNames.size()+fSigNames.size()+fNormSigNames.size()+fOverSigNames.size();
        if(hasData) Nrows ++;
        Nrows ++; // for "Uncertainty"
        double legHeight = ((Nrows+fLegendNColumns-1)/fLegendNColumns)*textHeight;
        leg  = new TLegend(legX1,legY-legHeight, legX2,legY);
        leg->SetNColumns(fLegendNColumns);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        if(TRExFitter::LEGENDRIGHT) leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetMargin(fLegendNColumns*((0.85*pad0->GetWh())/(1.*pad0->GetWw()))*textHeight/(legX2-legX1));
        //
        // Draws data in the legend only is real data
        if(hasData){
            if(TRExFitter::REMOVEXERRORS) leg->AddEntry(h_data.get(),fDataName.c_str(),"ep");
            else                          leg->AddEntry(h_data.get(),fDataName.c_str(),"lep");
        }
        //
        // Signal and background legend
        for(unsigned int i_smp=0;i_smp<fSigNames.size();i_smp++)     leg->AddEntry(h_signal[i_smp], fSigNames[i_smp].c_str(),"f");
        if(TRExFitter::OPTION["TRExbbStyle"]==0){
            for(unsigned int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" *").c_str(),"l");
            for(unsigned int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }
        for(unsigned int i_smp=0;i_smp<fBkgNames.size();i_smp++)     leg->AddEntry(h_bkg[i_smp], fBkgNames[i_smp].c_str(),"f");
        if(TRExFitter::OPTION["TRExbbStyle"]!=0) leg->AddEntry(g_tot.get(),"Total unc.","f");
        else leg->AddEntry(g_tot.get(),"Uncertainty","f");
        if(TRExFitter::OPTION["TRExbbStyle"]!=0){
            for(unsigned int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" (norm)").c_str(),"l");
            for(unsigned int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }
        //
        if(TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit) leg->AddEntry(h_tot_bkg_prefit,"Pre-Fit Bkgd.","l");
        //
        leg->Draw();
        //
        if(TRExFitter::OPTION["TRExbbStyle"]==0 && fNormSigNames.size()>0){
            myText(legX1,0.96,  1,"*: normalised to total Bkg.");
        }
    }

    //
    // Ratio pad: drawing dummy histogram
    //
    if(pad1!=nullptr){
        pad1->cd();
        pad1->GetFrame()->SetY1(2);
        TH1* h_dummy2 = (TH1*)h_tot->Clone("h_dummy2");
        h_dummy2->Scale(0);
        if(pad0->GetWw() > pad0->GetWh()) h_dummy2->GetYaxis()->SetTickLength(0.01);
        h_dummy2->Draw("HIST");
        h_dummy2->GetYaxis()->SetTitleOffset(1.*h_dummy->GetYaxis()->GetTitleOffset());
        if (fXaxisRange.size() > 1){
            h_dummy2->GetXaxis()->SetRangeUser(fXaxisRange.at(0),fXaxisRange.at(1));
        }

        //
        // Initialising the ratios
        //    h_ratio: is the real Data/MC ratio
        //    h_ratio2: is a MC/MC ratio to plot the uncertainty band
        //
        if(fRatioType=="S/B" || fRatioType=="S/SQRT(B)" || fRatioType=="S/SQRT(S+B)"){ // ...
            if(fSigNames.size()>0)          h_ratio.reset(static_cast<TH1*>(h_signal[0] ->Clone("h_ratio")));
            else if(fNormSigNames.size()>0){
                h_ratio.reset(static_cast<TH1*>(h_normsig[0]->Clone("h_ratio")));
                h_ratio->Scale(1./signalScale[0]);
            }
            else if(fOverSigNames.size()>0) h_ratio.reset(static_cast<TH1*>(h_oversig[0]->Clone("h_ratio")));
            else{
                h_ratio.reset(static_cast<TH1*>(h_tot->Clone("h_ratio")));
                h_ratio->Scale(0);
            }
        }
        else{
            h_ratio.reset(static_cast<TH1*>(h_data->Clone("h_ratio")));
        }

        // in case of S/B,.. and several signal samples, build other ratios
        std::vector<TH1*> h_addRatioVec;
        if(fRatioType=="S/B" || fRatioType=="S/SQRT(B)" || fRatioType=="S/SQRT(S+B)"){ // ...
            if(fSigNames.size()>1){
                for(unsigned int i_sig=1;i_sig<fSigNames.size();i_sig++){
                    h_addRatioVec.push_back((TH1*)h_signal[i_sig] ->Clone(Form("h_ratio_%d",i_sig)));
                }
            }
            else if(fNormSigNames.size()>1){
                for(unsigned int i_sig=1;i_sig<fNormSigNames.size();i_sig++){
                    h_addRatioVec.push_back((TH1*)h_normsig[i_sig]->Clone(Form("h_ratio_%d",i_sig)));
                    h_addRatioVec[i_sig-1]->Scale(1./signalScale[i_sig]);
                }
            }
            else if(fOverSigNames.size()>1){
                for(unsigned int i_sig=1;i_sig<fOverSigNames.size();i_sig++){
                    h_addRatioVec.push_back((TH1*)h_oversig[i_sig]->Clone(Form("h_ratio_%d",i_sig)));
                }
            }
        }

        h_tot_nosyst.reset(static_cast<TH1*>(h_tot->Clone("h_tot_nosyst")));
        for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){
            h_tot_nosyst->SetBinError(i_bin,0);
        }
        TGraphAsymmErrors *g_ratio2 = (TGraphAsymmErrors*)g_tot->Clone("g_ratio2");

        //
        // Plots style
        //
        std::string ratioTitle = "";
        if(fRatioYtitle == ""){
            if(fRatioType == "DATA/MC")     ratioTitle = "Data / Pred.";
            if(fRatioType == "DATA/BKG")    ratioTitle = "Data / Bkg.";
            if(fRatioType == "S/B")         ratioTitle = "S / B";
            if(fRatioType == "S/SQRT(B)")   ratioTitle = "S / #sqrt{B}";
            if(fRatioType == "S/SQRT(S+B)") ratioTitle = "S / #sqrt{S+B}";
        }
        else{
            ratioTitle = fRatioYtitle;
        }
        if(ratioTitle == "-") ratioTitle = "";
        h_dummy2->GetYaxis()->SetTitle(ratioTitle.c_str());
        h_dummy2->GetYaxis()->SetNdivisions(504,false);
        gStyle->SetEndErrorSize(0);

        //
        // Compute Data/MC ratio
        //
        h_ratio->Divide(h_tot_nosyst.get());
        for(auto h_tmp : h_addRatioVec) h_tmp->Divide(h_tot_nosyst.get());
        if(TRExFitter::OPRATIO) h_ratio->SetMarkerStyle(24);
        else                    h_ratio->SetMarkerStyle(h_data->GetMarkerStyle());
        h_ratio->SetMarkerSize(1.4);
        h_ratio->SetMarkerColor(kBlack);
        h_ratio->SetLineWidth(2);
        g_ratio = histToGraph(h_ratio.get());
        for(int i_bin=1;i_bin<=h_ratio->GetNbinsX();i_bin++){
            //For the ratio plot, the error is just to illustrate the "poisson uncertainty on the data"
            if(TRExFitter::REMOVEXERRORS){
                g_ratio->SetPointEXhigh( i_bin-1, 0. );
                g_ratio->SetPointEXlow(  i_bin-1, 0. );
            }
            g_ratio->SetPointEYhigh( i_bin-1,g_data->GetErrorYhigh(i_bin-1)/h_tot->GetBinContent(i_bin) );
            g_ratio->SetPointEYlow(  i_bin-1,g_data->GetErrorYlow(i_bin-1) /h_tot->GetBinContent(i_bin) );
        }

        //
        // Compute the MC/MC ratio (for uncertainty band in the bottom pad)
        //
        for(int i_bin=1;i_bin<h_tot_nosyst->GetNbinsX()+1;i_bin++){
            g_ratio2->SetPoint(i_bin-1,g_ratio2->GetX()[i_bin-1],g_ratio2->GetY()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
            g_ratio2->SetPointEXlow(i_bin-1,g_ratio2->GetEXlow()[i_bin-1]);
            g_ratio2->SetPointEXhigh(i_bin-1,g_ratio2->GetEXhigh()[i_bin-1]);
            if(h_tot_nosyst->GetBinContent(i_bin)>1e-4){
                g_ratio2->SetPointEYlow(i_bin-1,g_ratio2->GetEYlow()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
                g_ratio2->SetPointEYhigh(i_bin-1,g_ratio2->GetEYhigh()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
            }
            else{
                g_ratio2->SetPointEYlow(i_bin-1,0.);
                g_ratio2->SetPointEYhigh(i_bin-1,0.);
            }
        }

        //
        // Now draws everything
        //
        hline = std::make_unique<TLine>(h_dummy2->GetXaxis()->GetXmin(),1,h_dummy2->GetXaxis()->GetXmax(),1);
        hline->SetLineColor(kBlack);
        hline->SetLineWidth(2);
        hline->SetLineStyle(2);
        hline->Draw();
        //
        h_dummy2->SetMinimum(fRatioYmin);
        h_dummy2->SetMaximum(fRatioYmax);
        //
        h_dummy2->GetXaxis()->SetTitle(h_dummy->GetXaxis()->GetTitle());
        h_dummy2->GetXaxis()->SetTitleOffset(5);
        //
        h_dummy->GetXaxis()->SetTitle("");
        h_dummy->GetXaxis()->SetLabelSize(0);

        if(fRatioType=="DATA/MC" || fRatioType=="DATA/BKG"){
            g_ratio2->Draw("sameE2");
        }

        bool customLabels = false;
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            if(((std::string)h_dummy->GetXaxis()->GetBinLabel(i_bin))!=""){
                h_dummy2->GetXaxis()->SetBinLabel( i_bin, h_dummy->GetXaxis()->GetBinLabel(i_bin));
                customLabels = true;
            }
        }

        if(fRatioType=="S/B" || fRatioType=="S/SQRT(B)" || fRatioType=="S/SQRT(S+B)"){ // ...
            h_ratio->SetFillStyle(0);
            h_ratio->SetLineColor(h_ratio->GetFillColor());
            h_ratio->Draw("HIST same");
            for(auto h_tmp : h_addRatioVec){
                h_tmp->SetFillStyle(0);
                h_tmp->SetLineColor(h_tmp->GetFillColor());
                h_tmp->Draw("HIST same");
            }
        }
        else if(hasData){
            g_ratio->Draw("pe0");
        }

        //
        // Mark blinded bins in ratio pad as  well
        //
        if(h_blind!=nullptr){
            TH1D* h_blindratio = (TH1D*)h_blind->Clone("h_blindratio");
            h_blindratio->Scale(2.);
            h_blindratio->Draw("HIST same");
        }

        if(fBinLabel.size() > 1 && fBinLabel[1]!="") h_dummy2->GetXaxis()->LabelsOption("d");
        h_dummy2->GetXaxis()->SetLabelOffset( h_dummy2->GetXaxis()->GetLabelOffset()+0.02 );
        if(customLabels && h_dummy->GetNbinsX()>10) h_dummy2->GetXaxis()->SetLabelSize(0.66*h_dummy2->GetXaxis()->GetLabelSize() );
        if(customLabels) h_dummy2->GetXaxis()->SetLabelOffset( h_dummy2->GetXaxis()->GetLabelOffset()+0.02 );
        gPad->RedrawAxis();

        h_dummy2->GetYaxis()->ChangeLabel(-1,-1,-1,-1,-1,-1," ");

        //
        // Add arrows when the ratio is beyond the limits of the ratio plot
        //
        if(fRatioType=="DATA/MC" || fRatioType=="DATA/BKG"){
            for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){

                if (i_bin==0 || i_bin>h_tot_nosyst->GetNbinsX()) continue; //skip under/overflow bins

                double val = h_ratio->GetBinContent(i_bin);

                double maxRange = h_dummy2->GetMaximum();
                double minRange = h_dummy2->GetMinimum();

                int isUp=0; //1==up, 0==nothing, -1==down
                if ( val<minRange ) isUp=-1;
                else if (val>maxRange ) isUp=1;
                if (val==0) isUp=0;

                if (isUp!=0) {
                    TArrow *arrow;
                    if (isUp==1) arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmax-0.05*(fRatioYmax-fRatioYmin), h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmax,0.030/(pad0->GetWw()/596.),"|>");
                    else         arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmin+0.05*(fRatioYmax-fRatioYmin), h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmin,0.030/(pad0->GetWw()/596.),"|>");
                    arrow->SetFillColor(10);
                    arrow->SetFillStyle(1001);
                    arrow->SetLineColor(kBlue-7);
                    arrow->SetLineWidth(2);
                    arrow->SetAngle(40);
                    arrow->Draw();
                }
            }
        }
        // ---

        pad1->cd();
        KSlab = std::make_unique<TLatex>();
        KSlab->SetNDC(1);
        KSlab->SetTextFont(42);
        KSlab->SetTextSize(0.1);
        std::string kslab = "";
        if(Chi2val >= 0)  kslab += Form("   #chi^{2}/ndf = %.1f",Chi2val);
        if(NDF >= 0)      kslab += Form(" / %d",NDF);
        if(Chi2prob >= 0) kslab += Form("  #chi^{2}prob = %.2f",Chi2prob);
        if(KSprob >= 0)   kslab += Form("  KS prob = %.2f",KSprob);
        KSlab->DrawLatex(0.15,0.9,kslab.c_str());
        //
        pad0->cd();
    }

    //
    // Set bin width and eventually divide larger bins by this bin width
    if(fBinWidth>0){
        for(unsigned int i_smp=0;i_smp<fSigNames.size();i_smp++)      SetHistBinWidth(h_signal[i_smp], fBinWidth);
        for(unsigned int i_smp=0;i_smp<fNormSigNames.size();i_smp++)  SetHistBinWidth(h_normsig[i_smp],fBinWidth);
        for(unsigned int i_smp=0;i_smp<fOverSigNames.size();i_smp++)  SetHistBinWidth(h_oversig[i_smp],fBinWidth);
        for(unsigned int i_smp=0;i_smp<fBkgNames.size();i_smp++)      SetHistBinWidth(h_bkg[i_smp],    fBinWidth);
        //
        if(h_tot) SetHistBinWidth(h_tot.get(),fBinWidth);
        if(g_tot) SetGraphBinWidth(g_tot.get(),fBinWidth);
        if(h_data) SetHistBinWidth(h_data.get(),fBinWidth);
        if(g_data) SetGraphBinWidth(g_data.get(),fBinWidth);
        // try to guess y axis label...
        if(ytitle=="Events"){
            if(xtitle.find("GeV")!=std::string::npos){
                if((int)fBinWidth==fBinWidth) ytitle = Form("Events / %.0f GeV",fBinWidth);
                else if((int)(fBinWidth*10)==(fBinWidth*10)) ytitle = Form("Events / %.1f GeV",fBinWidth);
                else if((int)(fBinWidth*100)==(fBinWidth*100)) ytitle = Form("Events / %.2f GeV",fBinWidth);
                else if((int)(fBinWidth*1000)==(fBinWidth*1000)) ytitle = Form("Events / %.3f GeV",fBinWidth);
                // ...
            }
            else{
                if((int)fBinWidth==fBinWidth) ytitle = Form("Events / %.0f",fBinWidth);
                else if((int)(fBinWidth*10)==(fBinWidth*10)) ytitle = Form("Events / %.1f",fBinWidth);
                else if((int)(fBinWidth*100)==(fBinWidth*100)) ytitle = Form("Events / %.2f",fBinWidth);
                else if((int)(fBinWidth*1000)==(fBinWidth*1000)) ytitle = Form("Events / %.3f",fBinWidth);
                // ...
            }
            h_dummy->GetYaxis()->SetTitle(ytitle.c_str());
        }
    }

    // turn off x-error bars
    if(TRExFitter::REMOVEXERRORS){
        for (int i=0; i < g_data->GetN(); i++) {
            g_data->SetPointEXlow(i,0);
            g_data->SetPointEXhigh(i,0);
        }
    }

    // Fix y max
    //
    double yMax = 0.;
    // take into account also total prediction uncertainty
    for(int i_bin=1;i_bin<h_tot->GetNbinsX()+1;i_bin++){
        double y = h_tot->GetBinContent(i_bin);
        if(y>yMax) yMax = y;
        if(hasData && h_data!=nullptr && g_data!=nullptr){
            if(h_data->Integral()>0 && h_data->GetBinContent(i_bin)>0 && g_data->GetY()[i_bin-1]>0 && g_data->GetEYhigh()[i_bin-1]>0){
                y = h_data->GetBinContent(i_bin)+g_data->GetEYhigh()[i_bin-1];
                if(y>yMax) yMax = y;
            }
        }
    }
    //
    if(options.find("log")==std::string::npos){
        if(fYmax!=0) h_dummy->SetMaximum(fYmax);
        else         h_dummy->SetMaximum(yMaxScale*yMax);
        if(fYmin>0)  h_dummy->SetMinimum(fYmin);
        else         h_dummy->SetMinimum(0.);
    }
    else{
        if(fYmax!=0) h_dummy->SetMaximum(fYmax);
        else         h_dummy->SetMaximum(yMax*pow(10,yMaxScale));
        if(fYmin>0)  h_dummy->SetMinimum(fYmin);
        else         h_dummy->SetMinimum(1.);
    }

    if(h_blind!=nullptr){
        h_blind->Scale(h_dummy->GetMaximum());
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SaveAs(const std::string& name) const{
    c->SaveAs(name.c_str());
}

//_____________________________________________________________________________
//
void TRExPlot::WriteToFile(const std::string& name) const{
    TDirectory *here = gDirectory;
    std::unique_ptr<TFile> f (TFile::Open(name.c_str(),"RECREATE"));
    if (!f) {
        WriteWarningStatus("TRExPlot::WriteToFile", "Cannot open file: " + name + ". Not writing into a file!");
        return;
    }
    f->cd();
    if(h_data) h_data->Write(Form("h_%s",fDataName.c_str()),TObject::kOverwrite);
    h_tot->Write("h_totErr",TObject::kOverwrite);
    if(g_tot) g_tot->Write("g_totErr",TObject::kOverwrite);
    for(int i_smp=fBkgNames.size()-1;i_smp>=0;i_smp--){
        h_bkg[i_smp]->Write(Form("h_%s",fBkgNames[i_smp].c_str()),TObject::kOverwrite);
    }
    for(int i_smp=fSigNames.size()-1;i_smp>=0;i_smp--){
        h_signal[i_smp]->Write(Form("h_%s",fSigNames[i_smp].c_str()),TObject::kOverwrite);
        if(h_normsig[i_smp]) h_normsig[i_smp]->Write(Form("h_%s_norm",fSigNames[i_smp].c_str()),TObject::kOverwrite);
    }
    here->cd();
    f->Close();
}

//_____________________________________________________________________________
//
void TRExPlot::SetBinBlinding(const std::vector<int>& bins){
    fBlindedBins = bins;
}

//_____________________________________________________________________________
//
const TH1D* TRExPlot::GetBlindingHisto() const {
    return fBlinding.get();
}

//_____________________________________________________________________________
//
const std::vector<int>& TRExPlot::GetBlindedBins() const {
    return fBlindedBins;
}

//_____________________________________________________________________________
//
void TRExPlot::SetBlindingHisto(const TH1* h) {
    fBlinding.reset(static_cast<TH1D*>(h->Clone()));
    if (fBlinding) fBlinding->SetDirectory(nullptr);
}

//_____________________________________________________________________________
// function to get asymmetric error bars for hists (Used in WZ observation)
double GC_up(double data) {
    if (data == 0 ) return 0;
    return 0.5*TMath::ChisquareQuantile(1.-0.1586555,2.*(data+1))-data;
}

//_____________________________________________________________________________
//
double GC_down(double data) {
    if (data == 0 ) return 0;
    return data-0.5*TMath::ChisquareQuantile(0.1586555,2.*data);
}

//_____________________________________________________________________________
//
std::unique_ptr<TGraphAsymmErrors> poissonize(const TH1 *h) {
    std::unique_ptr<TGraphAsymmErrors> gr = std::make_unique<TGraphAsymmErrors>(h);
    for (int i = 0; i < gr->GetN(); i++) {
        double content = gr->GetErrorYhigh(i) * gr->GetErrorYhigh(i); // this to fix the case of the merged plots, where histograms (even data) are scaled; so the actual content is the square of the stat. error (right?)
        gr->SetPointError(i,0.499*h->GetBinWidth(i+1),0.5*h->GetBinWidth(i+1),GC_down(content),GC_up(content));
        if(h->GetBinContent(i+1)==0){
            gr->SetPoint(i,gr->GetX()[i],-1);
            gr->SetPointError(i,0,0,0,0);
        }
    }
    gr->SetMarkerSize(h->GetMarkerSize());
    gr->SetMarkerColor(h->GetMarkerColor());
    gr->SetMarkerStyle(h->GetMarkerStyle());
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    gr->SetLineStyle(h->GetLineStyle());
    return gr;
}

//_____________________________________________________________________________
//
std::unique_ptr<TGraphAsymmErrors> histToGraph(const TH1* h){
    std::unique_ptr<TGraphAsymmErrors> gr = std::make_unique<TGraphAsymmErrors>(h);
    for (int i = 0; i < gr->GetN(); i++) {
        gr->SetPointEXlow(i,0.499*h->GetBinWidth(i+1));
        gr->SetPointEXhigh(i,0.5*h->GetBinWidth(i+1));
        if(h->GetBinContent(i+1)==0){
            gr->SetPoint(i,gr->GetX()[i],-1);
            gr->SetPointError(i,0,0,0,0);
        }
    }
    gr->SetMarkerStyle(h->GetMarkerStyle());
    gr->SetMarkerSize(h->GetMarkerSize());
    gr->SetMarkerColor(h->GetMarkerColor());
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    gr->SetLineStyle(h->GetLineStyle());
    return gr;
}

//_____________________________________________________________________________
//
void SetHistBinWidth(TH1* h,double width){
    static const double epsilon = 0.00000001;
    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
        if(std::abs(h->GetBinWidth(i_bin)-width)>epsilon){
            h->SetBinContent(i_bin,h->GetBinContent(i_bin)*width/h->GetBinWidth(i_bin));
            h->SetBinError(  i_bin,h->GetBinError(i_bin)  *width/h->GetBinWidth(i_bin));
        }
    }
}

//_____________________________________________________________________________
//
void SetGraphBinWidth(TGraphAsymmErrors* g,double width){
    static const double epsilon = 0.00000001;
    for(int i_bin=0;i_bin<g->GetN();i_bin++){
        const double w = g->GetErrorXhigh(i_bin)+g->GetErrorXlow(i_bin);
        if(std::abs(w-width)>epsilon){
            g->SetPoint(      i_bin,g->GetX()[i_bin], g->GetY()[i_bin]*width/w);
            g->SetPointEYhigh(i_bin,g->GetErrorYhigh(i_bin)*width/w);
            g->SetPointEYlow( i_bin,g->GetErrorYlow(i_bin) *width/w);
        }
    }
}
