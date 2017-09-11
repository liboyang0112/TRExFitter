#include "TtHFitter/TthPlot.h"

using namespace std;

//_____________________________________________________________________________
//
TthPlot::TthPlot(string name,int canvasWidth,int canvasHeight){
    fName = name;
    c = new TCanvas(fName.c_str(),fName.c_str(),canvasWidth,canvasHeight);
    //
    pad0 = new TPad("pad0","pad0",0,0.20,1,1,0,0,0);
//     pad0->SetTickx(false);
//     pad0->SetTicky(false);
    pad0->SetTicks(1,1);
    pad0->SetTopMargin(0.05);
    pad0->SetBottomMargin(0.1);
    pad0->SetLeftMargin(0.14);
    pad0->SetRightMargin(0.05);
    pad0->SetFrameBorderMode(0);
    pad0->SetFillStyle(0);
    //
    pad1 = new TPad("pad1","pad1",0,0,1,0.28,0,0,0);
//     pad1->SetTickx(false);
//     pad1->SetTicky(false);
    pad1->SetTicks(1,1);
    pad1->SetTopMargin(0.0);
    pad1->SetBottomMargin(0.37);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->SetFrameBorderMode(0);
    pad1->SetFillStyle(0);
    //
    if(canvasWidth>canvasHeight){ // FIXME
        pad0->SetLeftMargin(0.10);
        pad1->SetLeftMargin(0.10);
    }
    //
    pad1->Draw();
    pad0->Draw();
    pad0->cd();
    h_stack = new THStack("h_stack","h_stack");
    h_tot = 0x0;
    g_tot = 0x0;
    xtitle = "Variable [GeV]";
    ytitle = "Events";
    fLabels.clear();
    fLumi = "20.3 fb^{-1}";
    fCME = "8 TeV";
    fATLASlabel = "none";
    yMaxScale = 2.;
    Chi2val = -1;
    NDF = -1;
    Chi2prob = -1;
    KSprob = -1;
    //
    h_data = 0x0;
    g_data = 0x0;
    for(int i_smp=0;i_smp<MAXSAMPLES;i_smp++){
        h_bkg[i_smp]     = 0x0;
        h_signal[i_smp]  = 0x0;
        h_normsig[i_smp] = 0x0;
        h_oversig[i_smp] = 0x0;
    }

    //
    fIsNjet = false;
    fShowYields = false;
    //
    for(int i_bin=0;i_bin<MAXbins;i_bin++)
        fBinLabel[i_bin] = "";
    //
    fSigNames.clear();
    fNormSigNames.clear();
    fOverSigNames.clear();
    fBkgNames.clear();

    fBinWidth = -1;
    fLumiScale = 1.;
    fBlindingThreshold = -1; // if <0, no blinding
    fLegendNColumns = 0;
    
    fYmin = 0;
    fYmax = 0;
    fRatioYmin = 0.;
    fRatioYmax = 2.;
    
    h_blinding = 0x0;
    
    h_tot_bkg_prefit = 0x0;
}

//_____________________________________________________________________________
//
void TthPlot::SetChannel(string name){
    fLabels.clear();
    fLabels.push_back(name);
}

//_____________________________________________________________________________
//
void TthPlot::AddLabel(string name){
    fLabels.push_back(name);
}

//_____________________________________________________________________________
//
void TthPlot::SetLumi(string name){
    fLumi = name;
}

//_____________________________________________________________________________
//
void TthPlot::SetLumiScale(float scale){
    fLumiScale = scale;
}

//_____________________________________________________________________________
//
void TthPlot::SetCME(string name){
    fCME = name;
}

//_____________________________________________________________________________
//
void TthPlot::SetXaxis(string name,bool isNjet){
    xtitle = name;
    fIsNjet = isNjet;
}

//_____________________________________________________________________________
//
void TthPlot::SetYaxis(string name){
    ytitle = name;
}

//_____________________________________________________________________________
//
void TthPlot::SetYmaxScale(float scale){
    yMaxScale = scale;
}

//_____________________________________________________________________________
//
void TthPlot::SetBinLabel(int bin,string name){
    fBinLabel[bin] = name;
}

//_____________________________________________________________________________
//
void TthPlot::SetBinWidth(float width){
    fBinWidth = width;
}

//_____________________________________________________________________________
//
void TthPlot::SetData(TH1* h,string name){
    h_data = (TH1*)h->Clone();
    data_name = name;
}

//_____________________________________________________________________________
//
void TthPlot::AddSignal(TH1* h,string name){
    // if already there...
//     if(std::find(fSigNames.begin(),fSigNames.end(),name)!=fSigNames.end()){
//         h_signal[fSigNames.size()-1]->Add(h,fLumiScale);
//     }
    int idx = std::find(fSigNames.begin(),fSigNames.end(),name) - fSigNames.begin();
    if(idx<fSigNames.size()){
        h_signal[idx]->Add(h,fLumiScale);
    }
    else{
        h_signal[fSigNames.size()] = (TH1*)h->Clone();
        h_signal[fSigNames.size()]->Scale(fLumiScale);
        fSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TthPlot::AddNormSignal(TH1* h,string name){
    // if already there...
//     if(std::find(fNormSigNames.begin(),fNormSigNames.end(),name)!=fNormSigNames.end()){
//         h_normsig[fNormSigNames.size()-1]->Add(h,fLumiScale);
//     }
    int idx = std::find(fNormSigNames.begin(),fNormSigNames.end(),name) - fNormSigNames.begin();
    if(idx<fNormSigNames.size()){
        h_normsig[idx]->Add(h,fLumiScale);
    }
    else{
        h_normsig[fNormSigNames.size()] = (TH1*)h->Clone();
        h_normsig[fNormSigNames.size()]->Scale(fLumiScale);
        fNormSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TthPlot::AddOverSignal(TH1* h,string name){
    // if already there...
//     if(std::find(fOverSigNames.begin(),fOverSigNames.end(),name)!=fOverSigNames.end()){
//         h_oversig[fOverSigNames.size()-1]->Add(h,fLumiScale);
//     }
    int idx = std::find(fOverSigNames.begin(),fOverSigNames.end(),name) - fOverSigNames.begin();
    if(idx<fOverSigNames.size()){
        h_oversig[idx]->Add(h,fLumiScale);
    }
    else{
        h_oversig[fOverSigNames.size()] = (TH1*)h->Clone();
        h_oversig[fOverSigNames.size()]->Scale(fLumiScale);
        fOverSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TthPlot::AddBackground(TH1* h,string name){
    if(h_tot==0x0) h_tot = (TH1*)h->Clone();
    else h_tot->Add(h);
    // if already there...
//     if(std::find(fBkgNames.begin(),fBkgNames.end(),name)!=fBkgNames.end()){
//         h_bkg[fBkgNames.size()-1]->Add(h,fLumiScale);
//     }
//     int idx = -1;
//     for(int i=0;i<(int)fBkgNames.size();i++){
//         if(name==fBkgNames.at(ih)){
//             idx = i;
//             break;
//         }
//     }
//     if(idx>=0){
//         h_bkg[idx]->Add(h,fLumiScale);
//     }
    //
    int idx = std::find(fBkgNames.begin(),fBkgNames.end(),name) - fBkgNames.begin();
    if(idx<fBkgNames.size()){
        h_bkg[idx]->Add(h,fLumiScale);
    }
    else{
        h_bkg[fBkgNames.size()] = (TH1*)h->Clone();
        h_bkg[fBkgNames.size()]->Scale(fLumiScale);
        fBkgNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TthPlot::SetTotBkg(TH1* h){
    h_tot = (TH1*)h->Clone();
    h_tot->Scale(fLumiScale);
    g_tot = new TGraphAsymmErrors(h);
    for(int i=0;i<g_tot->GetN();i++){
        g_tot->GetY()[i]      *= fLumiScale;
//         g_tot->GetEXlow()[i]  *= fLumiScale;
//         g_tot->GetEXhigh()[i] *= fLumiScale;
        g_tot->GetEYlow()[i]  *= fLumiScale;
        g_tot->GetEYhigh()[i] *= fLumiScale;
    }
}

//_____________________________________________________________________________
//
void TthPlot::SetTotBkgAsym(TGraphAsymmErrors* g){
    g_tot = (TGraphAsymmErrors*)g->Clone();
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
void TthPlot::SetChi2KS(float chi2prob,float ksprob,float chi2val,int ndf){
    Chi2prob = chi2prob;
    KSprob = ksprob;
    Chi2val = chi2val;
    NDF = ndf;
}

//_____________________________________________________________________________
//
void TthPlot::BlindData(){
    //
    // Eventually blind bins
    //
    if(fBlindingThreshold>=0){
        if(h_data!=0x0 && fSigNames.size()>0 && h_tot!=0x0){
            if(h_blinding!=0x0){
                BlindDataHisto( h_data,h_blinding );
            }
            else{
                h_blinding = BlindDataHisto( h_data, h_tot, h_signal[0], fBlindingThreshold );
                // if more than one signal:
                if(fSigNames.size()>1){
                    for(unsigned int i_sig=1;i_sig<fSigNames.size();i_sig++){
                        h_blinding->Add( BlindDataHisto( h_data, h_tot, h_signal[i_sig], fBlindingThreshold ) );
                        h_blinding->Scale(2.);
                    }
                }
            }
        }
        else{
            std::cout << "TthPlot::WARNING: Either h_data (" << h_data << "), h_signal (" << h_signal << ") or h_tot (" << h_tot << ") not defined.";
            std::cout << " Blidning not possible. Skipped." << std::endl;
        }
    }
}

//_____________________________________________________________________________
//
TH1* TthPlot::GetTotBkg() const{
  TH1* h = (TH1*)h_tot->Clone("h_tot_bkg");
  for (int i=0; i<fSigNames.size(); ++i) {
    h->Add( h_signal[i], -1);
  }
  return h;
}

//_____________________________________________________________________________
//
void TthPlot::Draw(string options){

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
    gStyle->SetEndErrorSize(4.);
    pad0->cd();
    TH1* h_dummy = (TH1*)h_tot->Clone("h_dummy");
    h_dummy->Scale(0);
    if(pad0->GetWw() > pad0->GetWh()){
        h_dummy->GetYaxis()->SetTickLength(0.01);
        h_dummy->GetXaxis()->SetTickLength(0.02);
    }
    h_dummy->Draw("HIST");
    if(options.find("log")!=string::npos) pad0->SetLogy();

    if(g_tot==0x0) g_tot = new TGraphAsymmErrors(h_tot);

    //
    // Determines if the data is real (and computes the poisson uncertainty) or not
    //
    bool hasData = true;
    if(h_data){
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
        h_data = (TH1F*)h_tot->Clone("dummyData");//tajes data = total
        h_data->SetTitle("Asimov Data");
        g_data = new TGraphAsymmErrors(h_data);
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

    if( TtHFitter::PREFITONPOSTFIT ) {
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
    g_tot->SetFillColor(kBlue-7);
    g_tot->SetLineColor(kWhite);
    g_tot->SetLineWidth(0);
    g_tot->SetMarkerSize(0);
    g_tot->Draw("sameE2");

    //
    // Draw a normalized signal distribution
    //
    double signalScale = 1.;
    //if(h_normsig!=0x0){
        for(int i_smp=fNormSigNames.size()-1;i_smp>=0;i_smp--){
            signalScale = h_tot->Integral()/h_normsig[i_smp]->Integral();
            std::cout << "--- Signal " << fNormSigNames[i_smp] << " scaled by " <<signalScale << std::endl;
            h_normsig[i_smp]->Scale(signalScale);
            h_normsig[i_smp]->SetLineColor(h_normsig[i_smp]->GetFillColor());
            h_normsig[i_smp]->SetFillColor(0);
            h_normsig[i_smp]->SetFillStyle(0);
            h_normsig[i_smp]->SetLineStyle(2);
            h_normsig[i_smp]->SetLineWidth(2);
            h_normsig[i_smp]->Draw("HISTsame");
        }
    //}

    //
    // Draw a overlayed signal distribution
    //
    //if(h_oversig!=0x0){
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
    TH1F* h_blind = 0x0;
    if(h_blinding!=0x0){
        h_blind = (TH1F*)h_blinding->Clone("h_blind");
        h_blind->SetLineWidth(0);
        h_blind->SetLineColor(kGray);
        h_blind->SetFillColor(kGray);
        h_blind->SetFillStyle(3345);
//         h_blinding->SetFillColorAlpha(kWhite,0.75);
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
    else{
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            if(fBinLabel[i_bin]!="") h_dummy->GetXaxis()->SetBinLabel( i_bin, fBinLabel[i_bin].c_str());
        }
    }
    if(fBinLabel[1]!="") h_dummy->GetXaxis()->LabelsOption("d");
//     float offset = 2.3*(pad0->GetWh()/672.);
    float offset = 2.4*(pad0->GetWh()/672.);
    if(pad0->GetWw() > pad0->GetWh()) offset *= 0.8*596./pad0->GetWw();
    h_dummy->GetYaxis()->SetTitleOffset( offset );

    //
    // Fix / redraw axis
    //
    pad0->RedrawAxis();

    float textHeight = 0.05*(672./pad0->GetWh());
//     if(TtHFitter::OPTION["TtHbbStyle"]>0) textHeight *= 0.85;

    //
    // ATLAS labels
    //
    float labelX = 0.18;

    if(pad0->GetWw() > pad0->GetWh()) labelX = 0.12;

    if(fATLASlabel!="none") ATLASLabel(labelX,0.84+0.04,(char*)fATLASlabel.c_str());
    myText(labelX,0.84-textHeight+0.04,1,Form("#sqrt{s} = %s, %s",fCME.c_str(),fLumi.c_str()));//,0.045);
    for(unsigned int i_lab=0;i_lab<fLabels.size();i_lab++){
        myText(labelX,0.84-textHeight+0.04-(i_lab+1)*textHeight,1,Form("%s",fLabels[i_lab].c_str()));//,0.045);
    }

    float legX1 = 1-0.41*(596./pad0->GetWw())-0.08;
    if(TtHFitter::OPTION["FourTopStyle"]!=0 || TtHFitter::OPTION["TtHbbStyle"]!=0){
        legX1 = 1-0.5*(596./pad0->GetWw())-0.08;
    }
    float legX2 = 0.94;
    float legXmid = legX1+0.5*(legX2-legX1);

    if(fShowYields){
        legXmid = legX1+0.6*(legX2-legX1);
        leg  = new TLegend(legX1,0.93-(fBkgNames.size()+fSigNames.size()+2)*0.05, legXmid,0.93);
        leg1 = new TLegend(legXmid,leg->GetY1(), legX2,leg->GetY2());
        //
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize()*0.9);
        leg->SetMargin(0.22);
        leg1->SetFillStyle(0);
        leg1->SetBorderSize(0);
        leg1->SetTextAlign(32);
        leg1->SetTextFont(gStyle->GetTextFont());
        leg1->SetTextSize(gStyle->GetTextSize()*0.9);
        leg1->SetMargin(0.);

        if(hasData){//only add data in the legend if real data are here
            leg->AddEntry(h_data,data_name.c_str(),"lep");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_data->Integral()),"");
        }

        //Signal and background legends
        for(int i_smp=0;i_smp<fSigNames.size();i_smp++){
            leg->AddEntry(h_signal[i_smp], fSigNames[i_smp].c_str(),"f");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_signal[i_smp]->Integral()),"");
        }
        for(int i_smp=0;i_smp<fBkgNames.size();i_smp++){
            leg->AddEntry(h_bkg[i_smp], fBkgNames[i_smp].c_str(),"f");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_bkg[i_smp]->Integral()),"");
        }
        leg->AddEntry((TObject*)0,"Total","");
        leg1->AddEntry((TObject*)0,Form("%.1f",h_tot->Integral()),"");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0) leg->AddEntry(g_tot,"Total unc.","f");
        else leg->AddEntry(g_tot,"Uncertainty","f");
        leg1->AddEntry((TObject*)0," ","");

        if(TtHFitter::PREFITONPOSTFIT) {
          leg->AddEntry(h_tot_bkg_prefit,"Pre-Fit Bkgd.","l");
          leg1->AddEntry((TObject*)0," ","");
        }
        
        leg->Draw();
        leg1->Draw();
    }
    else if(fLegendNColumns==1){   //TtHFitter::OPTION["LegendNColumns"]==1){
        int Nrows = fBkgNames.size()+fSigNames.size()+fNormSigNames.size()+fOverSigNames.size();
        if(hasData) Nrows ++;
        Nrows ++; // for "Uncertainty"
        if(TtHFitter::OPTION["TtHbbStyle"]>0)
            leg  = new TLegend(legXmid+0.1*(legX2-legXmid),0.92-Nrows*textHeight*0.8, legX2,0.92);
        else
            leg  = new TLegend(legXmid+0.1*(legX2-legXmid),0.92-Nrows*textHeight, legX2,0.92);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);        if(TtHFitter::OPTION["TtHbbStyle"]>0)

        if(!TtHFitter::LEGENDLEFT) leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
          leg->SetTextSize(gStyle->GetTextSize());
//         leg->SetMargin(0.18);
        leg->SetMargin(0.22);

        //Draws data in the legend only is real data
        if(hasData)leg->AddEntry(h_data,data_name.c_str(),"lep");

        //Signal and background legend
        for(int i_smp=0;i_smp<fSigNames.size();i_smp++)     leg->AddEntry(h_signal[i_smp], fSigNames[i_smp].c_str(),"f");
//         for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], fNormSigNames[i_smp].c_str(),"f");
        if(TtHFitter::OPTION["TtHbbStyle"]==0){
            for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" *").c_str(),"l");
            for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }
        for(int i_smp=0;i_smp<fBkgNames.size();i_smp++)     leg->AddEntry(h_bkg[i_smp], fBkgNames[i_smp].c_str(),"f");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0) leg->AddEntry(g_tot,"Total unc.","f");
        else leg->AddEntry(g_tot,"Uncertainty","f");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0){
            for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" (norm.)").c_str(),"l");
            for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }
        
        if(TtHFitter::PREFITONPOSTFIT) leg->AddEntry(h_tot_bkg_prefit,"Pre-Fit Bkgd.","l");
        
        leg->Draw();
    }
    else if(fLegendNColumns==3){ //TtHFitter::OPTION["LegendNColumns"]==3){
        int Nrows = fBkgNames.size()+fSigNames.size()+fNormSigNames.size()+fOverSigNames.size();
        if(hasData) Nrows ++;
        Nrows ++; // for "Uncertainty"
        if(TtHFitter::OPTION["TtHbbStyle"]>0)
            leg  = new TLegend(0.4,0.92-((Nrows+2)/3)*textHeight, legX2,0.92);
        else
            leg  = new TLegend(0.4,0.92-((Nrows+2)/3)*textHeight, legX2,0.92);
        leg->SetNColumns(3);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        if(!TtHFitter::LEGENDLEFT) leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
          leg->SetTextSize(gStyle->GetTextSize());
        leg->SetMargin(0.20);

        //Draws data in the legend only is real data
        if(hasData)leg->AddEntry(h_data,data_name.c_str(),"lep");

        //Signal and background legend
        for(int i_smp=0;i_smp<fSigNames.size();i_smp++)     leg->AddEntry(h_signal[i_smp], fSigNames[i_smp].c_str(),"f");
        if(TtHFitter::OPTION["TtHbbStyle"]==0){
            for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" *").c_str(),"l");
            for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }
        for(int i_smp=0;i_smp<fBkgNames.size();i_smp++)     leg->AddEntry(h_bkg[i_smp], fBkgNames[i_smp].c_str(),"f");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0) leg->AddEntry(g_tot,"Total unc.","f");
        else leg->AddEntry(g_tot,"Uncertainty","f");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0){
            for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" (norm)").c_str(),"l");
            for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }

        if(TtHFitter::PREFITONPOSTFIT) leg->AddEntry(h_tot_bkg_prefit,"Pre-Fit Bkgd.","l");
        
        leg->Draw();
    }
    else{
        int Nrows = fBkgNames.size()+fSigNames.size()+fNormSigNames.size()+fOverSigNames.size();
        if(hasData) Nrows ++;
        Nrows ++; // for "Uncertainty"
        if(TtHFitter::OPTION["TtHbbStyle"]>0) legX1 = 0.43; // FIXME
        if(TtHFitter::OPTION["TtHbbStyle"]>0)
            leg  = new TLegend(legX1,0.8-((Nrows+1)/2)*0.05, legX2,0.8);
        else
            leg  = new TLegend(legX1,0.93-((Nrows+1)/2)*0.05, legX2,0.93);
        leg->SetNColumns(2);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(gStyle->GetTextFont());
        if(c->GetWw() > c->GetWh()) leg->SetTextSize(gStyle->GetTextSize());
        else                        leg->SetTextSize(gStyle->GetTextSize()*0.9);
        leg->SetMargin(0.22);

        //Draws data in the legend only is real data
        if(hasData)leg->AddEntry(h_data,data_name.c_str(),"lep");

        //Signal and background legend
        for(int i_smp=0;i_smp<fSigNames.size();i_smp++)     leg->AddEntry(h_signal[i_smp], fSigNames[i_smp].c_str(),"f");
        if(TtHFitter::OPTION["TtHbbStyle"]==0){
            for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" *").c_str(),"l");
            for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }
        for(int i_smp=0;i_smp<fBkgNames.size();i_smp++)     leg->AddEntry(h_bkg[i_smp], fBkgNames[i_smp].c_str(),"f");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0) leg->AddEntry(g_tot,"Total unc.","f");
        else leg->AddEntry(g_tot,"Uncertainty","f");
        if(TtHFitter::OPTION["TtHbbStyle"]!=0){
            for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp], (fNormSigNames[i_smp]+" (norm)").c_str(),"l");
            for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp], fOverSigNames[i_smp].c_str(),"l");
        }

        if(TtHFitter::PREFITONPOSTFIT) leg->AddEntry(h_tot_bkg_prefit,"Pre-Fit Bkgd.","l");
        
        leg->Draw();

        if(TtHFitter::OPTION["TtHbbStyle"]==0 && fNormSigNames.size()>0)
            myText(legX1,0.93-((Nrows+1)/2)*0.05 - 0.05,  1,"*: normalised to total Bkg.");
    }

    //
    // Ratio pad: drawing dummy histogram
    //
    pad1->cd();
    pad1->GetFrame()->SetY1(2);
    TH1* h_dummy2 = (TH1*)h_tot->Clone("h_dummy2");
    h_dummy2->Scale(0);
    if(pad0->GetWw() > pad0->GetWh()) h_dummy2->GetYaxis()->SetTickLength(0.01);
    h_dummy2->Draw("HIST");
    h_dummy2->GetYaxis()->SetTitleOffset(1.*h_dummy->GetYaxis()->GetTitleOffset());

    //
    // Initialising the ratios
    //    h_ratio: is the real Data/MC ratio
    //    h_ratio2: is a MC/MC ratio to plot the uncertainty band
    //
    TH1* h_ratio = 0x0;
    if(TtHFitter::OPTION["SoverBinRatio"]){
        if(fSigNames.size()>0) h_ratio = (TH1*)h_signal[0]->Clone("h_ratio");
        else if(fNormSigNames.size()>0) h_ratio = (TH1*)h_normsig[0]->Clone("h_ratio");
        else if(fOverSigNames.size()>0) h_ratio = (TH1*)h_oversig[0]->Clone("h_ratio");
        else                   h_ratio = (TH1*)h_tot      ->Clone("h_ratio");
    }
    else                                   h_ratio = (TH1*)h_data->Clone("h_ratio");

    TH1 *h_ratio2 = (TH1*)h_tot->Clone("h_ratio2");
    TH1 *h_tot_nosyst = (TH1*)h_tot->Clone("h_tot_nosyst");
    for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){
        h_tot_nosyst->SetBinError(i_bin,0);
    }
    TGraphAsymmErrors *g_ratio2 = (TGraphAsymmErrors*)g_tot->Clone("g_ratio2");

    //
    // Plots style
    //
    h_dummy2->SetTitle("Data/MC");
//     h_dummy2->GetYaxis()->CenterTitle();
    if(TtHFitter::OPTION["SoverBinRatio"]) h_dummy2->GetYaxis()->SetTitle("S / B");
    else                                   h_dummy2->GetYaxis()->SetTitle("Data / Pred. ");
    h_dummy2->GetYaxis()->SetLabelSize(0.8*h_ratio->GetYaxis()->GetLabelSize());
    if(pad0->GetWw() > pad0->GetWh()) h_dummy2->GetYaxis()->SetLabelOffset(0.01);
    else                              h_dummy2->GetYaxis()->SetLabelOffset(0.02);
    h_dummy2->GetYaxis()->SetNdivisions(504,false);
    gStyle->SetEndErrorSize(4.);
//     pad1 -> SetTicky();

    //
    // Compute Data/MC ratio
    //
    h_ratio->Divide(h_tot_nosyst);
    h_ratio->SetMarkerStyle(24);
    h_ratio->SetMarkerSize(1.4);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineWidth(2);
    TGraphAsymmErrors *g_ratio = histToGraph(h_ratio);
    for(int i_bin=1;i_bin<=h_ratio->GetNbinsX();i_bin++){
        //For the ratio plot, the error is just to illustrate the "poisson uncertainty on the data"
        if(h_data->GetBinContent(i_bin)<1 || h_ratio->GetBinContent(i_bin)<0.001){
            g_ratio->SetPointEYhigh( i_bin-1,0 );
            g_ratio->SetPointEYlow(  i_bin-1,0 );
            h_ratio->SetBinError( i_bin,0 );
        }
        else{
            g_ratio->SetPointEYhigh( i_bin-1,g_data->GetErrorYhigh(i_bin-1)/h_tot->GetBinContent(i_bin) );
            g_ratio->SetPointEYlow(  i_bin-1,g_data->GetErrorYlow(i_bin-1) /h_tot->GetBinContent(i_bin) );
        }
    }

    //
    // Compute the MC/MC ratio (for uncertainty band in the bottom pad)
    //
    for(int i_bin=1;i_bin<h_tot_nosyst->GetNbinsX()+1;i_bin++){
        g_ratio2->SetPoint(i_bin-1,g_ratio2->GetX()[i_bin-1],g_ratio2->GetY()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
        g_ratio2->SetPointEXlow(i_bin-1,g_ratio2->GetEXlow()[i_bin-1]);
        g_ratio2->SetPointEXhigh(i_bin-1,g_ratio2->GetEXhigh()[i_bin-1]);
        g_ratio2->SetPointEYlow(i_bin-1,g_ratio2->GetEYlow()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
        g_ratio2->SetPointEYhigh(i_bin-1,g_ratio2->GetEYhigh()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
    }

    //
    // Now draws everything
    //
    TLine *hline = new TLine(h_dummy2->GetXaxis()->GetXmin(),1,h_dummy2->GetXaxis()->GetXmax(),1);
//     hline->SetLineColor(kRed);
    hline->SetLineColor(kBlack);
    hline->SetLineWidth(2);
    hline->SetLineStyle(2);
    if(TtHFitter::OPTION["SoverBinRatio"]){
        h_ratio->SetFillStyle(0);
        h_ratio->SetLineColor(h_ratio->GetLineColor());
        h_ratio->Draw("HIST same");
    }
    else if(hasData){
        h_ratio->Draw("E0 same");
    }
    hline->Draw();
    //
    h_dummy2->SetMinimum(fRatioYmin);
    h_dummy2->SetMaximum(fRatioYmax);
    //
    h_dummy2->GetXaxis()->SetTitle(h_dummy->GetXaxis()->GetTitle());
//     h_dummy2->GetXaxis()->SetTitleOffset(5.);
//     h_dummy2->GetXaxis()->SetTitleOffset(5.*(pad0->GetWw()/596.));
    // FIXME
    h_dummy2->GetXaxis()->SetLabelSize( 0.9*h_dummy2->GetXaxis()->GetLabelSize() );
    h_dummy2->GetXaxis()->SetTitleOffset(5.05*(pad0->GetWw()/596.));
//     h_dummy2->GetXaxis()->SetTitleSize( h_dummy2->GetXaxis()->GetTitleSize() );
    //
    h_dummy->GetXaxis()->SetTitle("");
    h_dummy->GetXaxis()->SetLabelSize(0);

    g_ratio2->Draw("sameE2");

    bool customLabels = false;
    for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
        if(((string)h_dummy->GetXaxis()->GetBinLabel(i_bin))!=""){
            h_dummy2->GetXaxis()->SetBinLabel( i_bin, h_dummy->GetXaxis()->GetBinLabel(i_bin));
            customLabels = true;
        }
    }

    //
    // Mark blinded bins in ratio pad as  well
    //
    if(h_blind!=0x0){
        TH1F* h_blindratio = (TH1F*)h_blind->Clone("h_blindratio");
        h_blindratio->Scale(2.);
        h_blindratio->Draw("HIST same");
    }

    if(fBinLabel[1]!="") h_dummy2->GetXaxis()->LabelsOption("d");
    h_dummy2->GetXaxis()->SetLabelOffset( h_dummy2->GetXaxis()->GetLabelOffset()+0.02 );
    if(customLabels && h_dummy->GetNbinsX()>10) h_dummy2->GetXaxis()->SetLabelSize(0.66*h_dummy2->GetXaxis()->GetLabelSize() );
    if(customLabels) h_dummy2->GetXaxis()->SetLabelOffset( h_dummy2->GetXaxis()->GetLabelOffset()+0.02 );
    gPad->RedrawAxis();

    // to hide the upper limit (label) of the ratio plot
    TLine line(0.01,1,0.1,1);
    line.SetLineColor(kWhite);
    line.SetLineWidth(25);
    if(pad0->GetWw() > pad0->GetWh()) line.DrawLineNDC(0.05,1,0.090,1);
    else                              line.DrawLineNDC(0.07,1,0.135,1);

    //
    // Add arrows when the ratio is beyond the limits of the ratio plot
    //
    for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){

        if (i_bin==0 || i_bin>h_tot_nosyst->GetNbinsX()) continue; //skip under/overflow bins

        float val=h_ratio->GetBinContent(i_bin);

        double maxRange = h_dummy2->GetMaximum();
        double minRange = h_dummy2->GetMinimum();

        int isUp=0; //1==up, 0==nothing, -1==down
        if ( val<minRange ) isUp=-1;
        else if (val>maxRange ) isUp=1;
        if (val==0) isUp=0;

        if (isUp!=0) {
            TArrow *arrow;
//             if (isUp==1) arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),1.45, h_ratio->GetXaxis()->GetBinCenter(i_bin),1.5,0.030,"|>");
//             else         arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),0.55, h_ratio->GetXaxis()->GetBinCenter(i_bin),0.5,0.030,"|>");
            if (isUp==1) arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmax-0.05*(fRatioYmax-fRatioYmin), h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmax,0.030,"|>");
            else         arrow = new TArrow(h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmin+0.05*(fRatioYmax-fRatioYmin), h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmin,0.030,"|>");
            arrow->SetFillColor(10);
            arrow->SetFillStyle(1001);
            arrow->SetLineColor(kBlue-7);
            arrow->SetLineWidth(2);
            arrow->SetAngle(40);
            //       arrow->Draw();
            //
            // Note: the following is not needed if using option "0" when Drawing...
            //       TLine *fix=0;
            //       TLine *fixUp=0;
            //       TLine *fixDo=0;
            //       cout << " bin: " << i_bin << "  with val: " << val << " and lower: " << val-h_ratio->GetBinError(i_bin) <<  endl;
            //       cout << " bin: " << i_bin << "  with val: " << val << " and upper: " << val+h_ratio->GetBinError(i_bin) <<  endl;
            //       if (val-h_ratio->GetBinError(i_bin)<1.5) {
            //         float x = h_ratio->GetXaxis()->GetBinCenter(i_bin);
            //         float w = h_ratio->GetXaxis()->GetBinWidth(i_bin);
            //         float y0 = TMath::Max(val-h_ratio->GetBinError(i_bin),0.50);
            //         float y1 = TMath::Min(val+h_ratio->GetBinError(i_bin),1.50);
            //         fix = new TLine( x, y0, x, y1 );
            //         fixUp = new TLine( x-w*0.07, y1, x+w*0.07, y1 );
            //         fixDo = new TLine( x-w*0.07, y0, x+w*0.07, y0 );
            //       }
            //       if (fix!=0) {
            //         fix->SetLineColor(h_ratio->GetLineColor());
            //         fixUp->SetLineColor(h_ratio->GetLineColor());
            //         fixDo->SetLineColor(h_ratio->GetLineColor());
            //         fix->SetLineWidth(h_ratio->GetLineWidth());
            //         fixUp->SetLineWidth(h_ratio->GetLineWidth());
            //         fixDo->SetLineWidth(h_ratio->GetLineWidth());
            //         fix->Draw("SAME");
            //         fixUp->SetLineColor(1);
            //         fixUp->SetLineWidth(1);
            //         fixUp->Draw("SAME");
            //         fixDo->SetLineColor(1);
            //         fixDo->SetLineWidth(1);
            //         fixDo->Draw("SAME");
            //       }
            arrow->Draw();
        }
    }
    // ---

    pad1->cd();
    TLatex *KSlab = new TLatex();
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

    //
    // Set bin width and eventually divide larger bins by this bin width
    if(fBinWidth>0){
        for(int i_smp=0;i_smp<fSigNames.size();i_smp++)      SetHistBinWidth(h_signal[i_smp], fBinWidth);
        for(int i_smp=0;i_smp<fNormSigNames.size();i_smp++)  SetHistBinWidth(h_normsig[i_smp],fBinWidth);
        for(int i_smp=0;i_smp<fOverSigNames.size();i_smp++)  SetHistBinWidth(h_oversig[i_smp],fBinWidth);
        for(int i_smp=0;i_smp<fBkgNames.size();i_smp++)      SetHistBinWidth(h_bkg[i_smp],    fBinWidth);
        //
        if(h_tot) SetHistBinWidth(h_tot,fBinWidth);
        if(g_tot) SetGraphBinWidth(g_tot,fBinWidth);
        if(h_data) SetHistBinWidth(h_data,fBinWidth);
        if(g_data) SetGraphBinWidth(g_data,fBinWidth);
        // try to guess y axis label...
        if(ytitle=="Events"){
            if(xtitle.find("GeV")!=string::npos){
                if((int)fBinWidth==fBinWidth) ytitle = Form("Events / %.0f GeV",fBinWidth);
                else if((int)(fBinWidth*10)==(fBinWidth*10)) ytitle = Form("Events / %.1f GeV",fBinWidth);
                else if((int)(fBinWidth*100)==(fBinWidth*100)) ytitle = Form("Events / %.2f GeV",fBinWidth);
                // ...
            }
            else{
                ytitle = Form("Events / %.2f",fBinWidth);
            }
            h_dummy->GetYaxis()->SetTitle(ytitle.c_str());
        }
    }

    // Fix y max
    //
    float yMax = 0.;
    float y;
    // take into account also total prediction uncertainty
    for(int i_bin=1;i_bin<h_tot->GetNbinsX()+1;i_bin++){
//         y = h_tot->GetBinContent(i_bin)+g_tot->GetEYhigh()[i_bin-1]; // creating problems
        y = h_tot->GetBinContent(i_bin);
        if(y>yMax) yMax = y;
        if(hasData && h_data!=0x0 && g_data!=0x0){
            if(h_data->Integral()>0 && h_data->GetBinContent(i_bin)>0 && g_data->GetY()[i_bin-1]>0 && g_data->GetEYhigh()[i_bin-1]>0){
                y = h_data->GetBinContent(i_bin)+g_data->GetEYhigh()[i_bin-1];
                if(y>yMax) yMax = y;
            }
        }
    }
    //
    if(options.find("log")==string::npos){
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
    
    if(h_blind!=0x0){
        h_blind->Scale(h_dummy->GetMaximum());
    }

    //
    // eventually make y-axis labels smaller...
    if(pad0->GetWw()<596. && h_dummy->GetMaximum()>10000){
        h_dummy->GetYaxis()->SetLabelSize( h_dummy->GetYaxis()->GetLabelSize()*0.75 );
    }
    else if(pad0->GetWw()<596. && h_dummy->GetMaximum()>1000){
        h_dummy->GetYaxis()->SetLabelSize( h_dummy->GetYaxis()->GetLabelSize()*0.9 );
    }
    
//     if(TtHFitter::OPTION["TtHbbStyle"]==0 && fNormSigNames.size()>0)
//         myText(0.4,0.96,  1,"#scale[0.75]{*: signal normalised to total background}");
}

//_____________________________________________________________________________
//
void TthPlot::SaveAs(string name){
    c->SaveAs(name.c_str());
}

//_____________________________________________________________________________
//
void TthPlot::WriteToFile(string name){
    TDirectory *here = gDirectory;
    TFile *f = new TFile(name.c_str(),"RECREATE");
    f->cd();
    if(h_data) h_data->Write(Form("h_%s",data_name.c_str()),TObject::kOverwrite);
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
//     f->~TFile();
    delete f;
}

//_____________________________________________________________________________
//
TCanvas* TthPlot::GetCanvas(){
    return c;
}

//_____________________________________________________________________________
//
void TthPlot::SetBinBlinding(bool on,float threshold){
    fBlindingThreshold = threshold;
    if(!on) fBlindingThreshold = -1;
    std::cout << "TthPlot::INFO: Setting blinding threshold = " << fBlindingThreshold << std::endl;
}

//_____________________________________________________________________________
//
void TthPlot::SetBinBlinding(bool on,TH1F* h_blind){
    h_blinding = h_blind;
    if(!on) fBlindingThreshold = -1;
    std::cout << "TthPlot::INFO: Setting blinding bins:";
    for(int i_bin=1;i_bin<h_blinding->GetNbinsX()+1;i_bin++){
        std::cout << " " << h_blinding->GetBinContent(i_bin);
    }
    std::cout << std::endl;
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
TGraphAsymmErrors* poissonize(TH1 *h) {
    vector<int> points_to_remove;
    TGraphAsymmErrors* gr= new TGraphAsymmErrors(h);
    int hBinCounter = 1;
    for (UInt_t i=0; i< (UInt_t)gr->GetN(); i++) {
        double content = (gr->GetY())[i];
        gr->SetPointError(i,0.499*h->GetBinWidth(hBinCounter),0.5*h->GetBinWidth(hBinCounter),GC_down(content),GC_up(content));
        //     if(content==0){
        if(content<0.1){ // FIXME
            gr->RemovePoint(i);
            i--;
        }
        hBinCounter++;
    }
    return gr;
}

//_____________________________________________________________________________
//
TGraphAsymmErrors* histToGraph(TH1* h){
    TGraphAsymmErrors* gr= new TGraphAsymmErrors(h);
    for (UInt_t i=0; i< (UInt_t)gr->GetN(); i++) {
        gr->SetPointError(i,0.499*h->GetBinWidth(i+1),0.5*h->GetBinWidth(i+1),0,0);
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
void SetHistBinWidth(TH1* h,float width){
    float epsilon = 0.00000001;
    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
        if(TMath::Abs(h->GetBinWidth(i_bin)-width)>epsilon){
            h->SetBinContent(i_bin,h->GetBinContent(i_bin)*width/h->GetBinWidth(i_bin));
            h->SetBinError(  i_bin,h->GetBinError(i_bin)  *width/h->GetBinWidth(i_bin));
        }
    }
}

//_____________________________________________________________________________
//
void SetGraphBinWidth(TGraphAsymmErrors* g,float width){
    float epsilon = 0.00000001;
    float w;
    for(int i_bin=0;i_bin<g->GetN();i_bin++){
        w = g->GetErrorXhigh(i_bin)+g->GetErrorXlow(i_bin);
        if(TMath::Abs(w-width)>epsilon){
            g->SetPoint(      i_bin,g->GetX()[i_bin], g->GetY()[i_bin]*width/w);
            g->SetPointEYhigh(i_bin,g->GetErrorYhigh(i_bin)*width/w);
            g->SetPointEYlow( i_bin,g->GetErrorYlow(i_bin) *width/w);
        }
    }
}

