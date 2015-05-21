#include "TtHFitter/HistoTools.h"
#include "TtHFitter/SampleHist.h"

// -------------------------------------------------------------------------------------------------
// SampleHist

//_____________________________________________________________________________
//
SampleHist::SampleHist(Sample *sample,TH1 *hist){
    fSample = sample;
    fName = fSample->fName;
    //
    //Nominal histogram configuration
    fHist = (TH1*)hist->Clone(Form("h_%s",fName.c_str()));
    fHist->SetFillColor(fSample->fFillColor);
    fHist->SetLineColor(fSample->fLineColor);
    //
    fHist_postFit = 0x0;
    //
    fHistoName = "";
    fFileName = "";
    fFitName = "";
    fNSyst = 0;
    fNNorm = 0;
    fRegionName = "Region";
    fRegionLabel = "Region";
    fVariableTitle = "Variable";
    fSystSmoothed = false;
    // add overall systematics and normFactors from sample
    for(int i_syst=0;i_syst<sample->fNSyst;i_syst++){
        if(sample->fSystematics[i_syst]->fType == Systematic::OVERALL)
            AddOverallSyst(sample->fSystematics[i_syst]->fName,sample->fSystematics[i_syst]->fOverallUp,sample->fSystematics[i_syst]->fOverallDown);
    }
    for(int i_norm=0;i_norm<sample->fNNorm;i_norm++){
        AddNormFactor(sample->fNormFactors[i_norm]);
    }
}

//_____________________________________________________________________________
//
SampleHist::SampleHist(Sample *sample, string histoName, string fileName){
    fSample = sample;
    fHist = HistFromFile(fileName,histoName);
    fHist->SetFillColor(fSample->fFillColor);
    fHist->SetLineColor(fSample->fLineColor);
    fHist->SetLineWidth(1);
    fName = fSample->fName;
    fHistoName = histoName;
    fFileName = fileName;
    fHist_postFit = 0x0;
    fNSyst = 0;
    fNNorm = 0;
    fRegionName = "Region";
    fVariableTitle = "Variable";
    fSystSmoothed = false;
    // add overall systematics and normFactors from sample
    for(int i_syst=0;i_syst<sample->fNSyst;i_syst++){
        if(sample->fSystematics[i_syst]->fType == Systematic::OVERALL)
            AddOverallSyst(sample->fSystematics[i_syst]->fName,sample->fSystematics[i_syst]->fOverallUp,sample->fSystematics[i_syst]->fOverallDown);
    }
    for(int i_norm=0;i_norm<sample->fNNorm;i_norm++){
        AddNormFactor(sample->fNormFactors[i_norm]);
    }
}

//_____________________________________________________________________________
//
SampleHist::~SampleHist(){
    delete fHist;
    delete fHist_postFit;
    for(unsigned int i = 0; i<fSyst.size(); ++i){
        if(fSyst[i]){
            delete fSyst[i];
        }
    }
    fSyst.clear();
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddOverallSyst(string name,float up,float down){
    SystematicHist *sh;
    // try if it's already there...
    sh = GetSystematic(name);
    // ... and if not create a new one
    if(sh==0x0){
      sh = new SystematicHist(name);
      fSyst.push_back(sh);
      fNSyst ++;
    }
    //
    sh->fHistUp   = (TH1*)fHist->Clone(Form("%s_%s_Up",fHist->GetName(),name.c_str()));
    sh->fHistDown = (TH1*)fHist->Clone(Form("%s_%s_Down",fHist->GetName(),name.c_str()));
    sh->fHistUp->Scale(1.+up);
    sh->fHistDown->Scale(1.+down);
    sh->fIsOverall = true;
    sh->fIsShape   = false;
    sh->fNormUp   = up;
    sh->fNormDown = down;
    return sh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddHistoSyst(string name,TH1* h_up,TH1* h_down){
    SystematicHist *sh;
    // try if it's already there...
    sh = GetSystematic(name);
    // ... and if not create a new one
    if(sh==0x0){
      sh = new SystematicHist(name);
      fSyst.push_back(sh);
      fNSyst ++;
    }
    //
    sh->fHistUp   = (TH1*)h_up->Clone(Form("%s_%s_Up",fHist->GetName(),name.c_str()));
    sh->fHistDown = (TH1*)h_down->Clone(Form("%s_%s_Down",fHist->GetName(),name.c_str()));
    sh->fHistShapeUp = (TH1*)h_up->Clone(Form("%s_%s_Shape_Up",fHist->GetName(),name.c_str()));
    sh->fHistShapeDown = (TH1*)h_down->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),name.c_str()));
    sh->fHistShapeUp->Scale(fHist->Integral() / sh->fHistShapeUp->Integral());
    sh->fHistShapeDown->Scale(fHist->Integral() / sh->fHistShapeDown->Integral());
    sh->fIsOverall = true;
    sh->fIsShape   = true;
    sh->fNormUp   = ( sh->fHistUp->Integral() - fHist->Integral() ) / fHist->Integral();
    sh->fNormDown = ( sh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
    if(sh->fNormUp == 0 && sh->fNormDown == 0) sh->fIsOverall = false;
    return sh;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::AddHistoSyst(string name,string histoName_up, string fileName_up,string histoName_down, string fileName_down){
    SystematicHist *sh;
    // try if it's already there...
    sh = GetSystematic(name);
    // ... and if not create a new one
    if(sh==0x0){
      sh = new SystematicHist(name);
      fSyst.push_back(sh);
      fNSyst ++;
    }
    //
    sh->fFileNameUp = fileName_up;
    sh->fFileNameDown = fileName_down;
    sh->fHistoNameUp = histoName_up;
    sh->fHistoNameDown = histoName_down;
    sh->fHistUp = HistFromFile(sh->fFileNameUp,sh->fHistoNameUp);
    sh->fHistDown = HistFromFile(sh->fFileNameDown,sh->fHistoNameDown);
    sh->fHistShapeUp = (TH1*)sh->fHistUp->Clone(Form("%s_%s_Shape_Up",fHist->GetName(),name.c_str()));
    sh->fHistShapeDown = (TH1*)sh->fHistDown->Clone(Form("%s_%s_Shape_Down",fHist->GetName(),name.c_str()));
    sh->fHistShapeUp->Scale(fHist->Integral() / sh->fHistShapeUp->Integral());
    sh->fHistShapeDown->Scale(fHist->Integral() / sh->fHistShapeDown->Integral());
    sh->fIsOverall = true;
    sh->fIsShape   = true;
    sh->fNormUp   = ( sh->fHistUp->Integral() - fHist->Integral() ) / fHist->Integral();
    sh->fNormDown = ( sh->fHistDown->Integral() - fHist->Integral() ) / fHist->Integral();
    if(sh->fNormUp == 0 && sh->fNormDown == 0) sh->fIsOverall = false;
    return sh;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::AddNormFactor(NormFactor *normFactor){
    NormFactor *norm = GetNormFactor(normFactor->fName);
    if(norm==0x0){
        fNormFactors.push_back(normFactor);
        fNNorm ++;
    }
    else{
        norm = normFactor;
    }
    return norm;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::AddNormFactor(string name,float nominal, float min, float max){
    NormFactor *norm = GetNormFactor(name);
    if(norm==0x0){
        fNormFactors.push_back(new NormFactor(name,nominal,min,max));
        fNNorm ++;
    }
    else{
        norm = new NormFactor(name,nominal,min,max);;
    }
    return norm;
}

//_____________________________________________________________________________
//
SystematicHist* SampleHist::GetSystematic(string systName){
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(systName == fSyst[i_syst]->fName) return fSyst[i_syst];
    }
    return 0x0;
}

//_____________________________________________________________________________
//
NormFactor* SampleHist::GetNormFactor(string name){
    for(int i_syst=0;i_syst<fNNorm;i_syst++){
        if(name == fNormFactors[i_syst]->fName) return fNormFactors[i_syst];
    }
    return 0x0;
}

//_____________________________________________________________________________
//
bool SampleHist::HasSyst(string name){
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSyst[i_syst]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
bool SampleHist::HasNorm(string name){
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        if(fNormFactors[i_norm]->fName == name) return true;
    }
    return false;
}

//_____________________________________________________________________________
//
void SampleHist::WriteToFile(){
    WriteHistToFile(fHist,fFileName);
    WriteHistToFile(HistoTools::TranformHistogramBinning(fHist),fFileName);
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        fSyst[i_syst]->WriteToFile();
    }
}

//_____________________________________________________________________________
//
void SampleHist::ReadFromFile(){
    fHist = HistFromFile(fFileName,fHistoName);
}

//_____________________________________________________________________________
//
void SampleHist::FixEmptyBins(){
    for(int i_bin=1;i_bin<=fHist->GetNbinsX();i_bin++){
        if(fHist->GetBinContent(i_bin)<=0) fHist->SetBinContent(i_bin,1e-6);
    }
}

//_____________________________________________________________________________
//
void SampleHist::Print(){
    cout << "      Sample: " << fName << "\t" << fHist->GetName() << endl;
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        fSyst[i_syst]->Print();
    }
    for(int i_norm=0;i_norm<fNNorm;i_norm++){
        fNormFactors[i_norm]->Print();
    }
}

//_____________________________________________________________________________
//
void SampleHist::Rebin(int ngroup, const Double_t* xbins){
    fHist->Rebin(ngroup,"",xbins);
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(fSyst[i_syst]->fHistUp!=0x0) fSyst[i_syst]->fHistUp->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistDown!=0x0) fSyst[i_syst]->fHistDown->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistShapeUp!=0x0) fSyst[i_syst]->fHistShapeUp->Rebin(ngroup,"",xbins);
        if(fSyst[i_syst]->fHistShapeDown!=0x0) fSyst[i_syst]->fHistShapeDown->Rebin(ngroup,"",xbins);
    }
}

//_____________________________________________________________________________
// this draws the control plots (for each systematic) with the syst variations for this region & sample
void SampleHist::DrawSystPlot(string syst, const bool dumpSystPlots){
    
    //Perform the treatment of the systematics
    SmoothSyst(syst);
    
    if(!dumpSystPlots) return;
    
    float yield_syst_up, yield_syst_down, yield_nominal;
    TCanvas *c = new TCanvas("c","c",800,600);
    TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
    h_nominal->SetLineColor(kBlack);
    h_nominal->SetLineWidth(2);
    h_nominal->SetLineStyle(2);
    h_nominal->SetFillStyle(0);
    TH1* h_1 = (TH1*)h_nominal->Clone("h_1");
    TH1* h_syst_up;
    TH1* h_syst_down;
    TH1* h_syst_up_orig;
    TH1* h_syst_down_orig;
    
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;
        h_nominal->SetMinimum(0);
        h_nominal->SetMaximum(h_nominal->GetMaximum()*2);
        h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
        h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
        h_syst_up_orig = (TH1*)fSyst[i_syst]->fHistUp_original->Clone();
        h_syst_down_orig = (TH1*)fSyst[i_syst]->fHistDown_original->Clone();
        h_syst_up->SetLineColor(kRed);
        h_syst_up->SetLineWidth(2);
        h_syst_up->SetLineStyle(1);
        h_syst_up->SetFillStyle(0);
        h_syst_down->SetLineColor(kBlue);
        h_syst_down->SetLineWidth(2);
        h_syst_down->SetLineStyle(1);
        h_syst_down->SetFillStyle(0);
        h_syst_up_orig->SetLineColor(kRed);
        h_syst_up_orig->SetLineWidth(2);
        h_syst_up_orig->SetLineStyle(2);
        h_syst_up_orig->SetFillStyle(0);
        h_syst_down_orig->SetLineColor(kBlue);
        h_syst_down_orig->SetLineWidth(2);
        h_syst_down_orig->SetLineStyle(2);
        h_syst_down_orig->SetFillStyle(0);

        yield_nominal = h_nominal->Integral();
        yield_syst_up = h_syst_up->Integral();
        yield_syst_down = h_syst_down->Integral();
        // draw Relative difference
        h_1->Scale(0);
        h_syst_up->Add(h_nominal,-1);
        h_syst_down->Add(h_nominal,-1);
        h_syst_up->Divide(h_nominal);
        h_syst_down->Divide(h_nominal);
        h_syst_up->Scale(100);
        h_syst_down->Scale(100);
        h_syst_up_orig->Add(h_nominal,-1);
        h_syst_down_orig->Add(h_nominal,-1);
        h_syst_up_orig->Divide(h_nominal);
        h_syst_down_orig->Divide(h_nominal);
        h_syst_up_orig->Scale(100);
        h_syst_down_orig->Scale(100);
        h_1->Draw("HIST");
        h_syst_down_orig->Draw("same HIST");
        h_syst_up_orig->Draw("same HIST");
        h_syst_down->Draw("same HIST");
        h_syst_up->Draw("same HIST");
        
        double ymax = 0;
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_up->GetMaximum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_down->GetMaximum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_up->GetMinimum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_down->GetMinimum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_up_orig->GetMaximum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_down_orig->GetMaximum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_up_orig->GetMinimum()));
        ymax = TMath::Max( ymax,TMath::Abs(h_syst_down_orig->GetMinimum()));
        h_1->SetMinimum(-ymax*2.1);
        h_1->SetMaximum( ymax*2.1);
        h_1->GetYaxis()->SetTitle("Relative difference [%]");
        h_1->GetXaxis()->SetTitle(fVariableTitle.c_str());

        // Creates a legend for the plot
        TLatex *tex = new TLatex();
        tex->SetNDC();
        if(fSyst[i_syst]->fSystematic!=0x0) tex->DrawLatex(0.2,0.89,Form("%s, %s",fSyst[i_syst]->fSystematic->fTitle.c_str(),fSample->fTitle.c_str()));
        else                                tex->DrawLatex(0.2,0.89,Form("%s, %s",fSyst[i_syst]->fName.c_str(),fSample->fTitle.c_str()));
//         tex->DrawLatex(0.2,0.89,Form("%s, %s",fSample->fSystematics[i_syst]->fTitle.c_str(),fSample->fTitle.c_str()));
        tex->DrawLatex(0.2,0.84,fRegionLabel.c_str());
        
        //Legend of the histograms
        TLegend *leg = new TLegend(0.6,0.83,0.9,0.93);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetTextFont(gStyle->GetTextFont());
        
        float acc_up = (yield_syst_up-yield_nominal)/yield_nominal;
        string sign_up =  "+";
        if(acc_up<0) sign_up = "-";
        float acc_down = (yield_syst_down-yield_nominal)/yield_nominal;
        string sign_down =  "+";
        if(acc_down<0) sign_down = "-";
        leg->AddEntry(h_syst_up,  Form("+1#sigma (%s%.1f%)",sign_up.c_str(),  TMath::Abs(acc_up  *100)),"l");
        leg->AddEntry(h_syst_down,Form("-1#sigma (%s%.1f%)",sign_down.c_str(),TMath::Abs(acc_down*100)),"l");
        leg->Draw();
        
        //Legend to define the line style
        TLegend *leg2 = new TLegend(0.6,0.77,0.9,0.82);
        leg2->SetFillStyle(0);
        leg2->SetBorderSize(0);
        leg2->SetNColumns(2);
        leg2->SetTextSize(gStyle->GetTextSize());
        leg2->SetTextFont(gStyle->GetTextFont());
        TH1F* h_syst_up_black = (TH1F*)h_syst_up -> Clone();
        h_syst_up_black -> SetLineColor(kBlack);
        TH1F* h_syst_up_origin_black = (TH1F*)h_syst_up_orig -> Clone();
        h_syst_up_origin_black -> SetLineColor(kBlack);
        leg2->AddEntry(h_syst_up_origin_black,"Original","l");
        leg2->AddEntry(h_syst_up_black,"Modified","l");
        leg2 -> Draw();
        
        gSystem->mkdir(fFitName.c_str());
        gSystem->mkdir((fFitName+"/systPlots/").c_str());
        
        const char* saveName = Form("%s/systPlots/%s_%s.png",fFitName.c_str(),fHist->GetName(),fSyst[i_syst]->fName.c_str());
        c->SaveAs(saveName);
        
        delete h_syst_up_black;
        delete h_syst_up_origin_black;
    }
    delete c;
}

//_____________________________________________________________________________
//
void SampleHist::SmoothSyst(string syst,bool force){
    if(fSystSmoothed && !force) return;
    TH1* h_nominal = (TH1*)fHist->Clone("h_nominal");
    TH1* h_syst_up;
    TH1* h_syst_down;
    
    for(int i_syst=0;i_syst<fNSyst;i_syst++){
        
        if(syst!="all" && fSyst[i_syst]->fName.find(syst)==string::npos) continue;
        
        h_syst_up = (TH1*)fSyst[i_syst]->fHistUp->Clone();
        h_syst_down = (TH1*)fSyst[i_syst]->fHistDown->Clone();
        
        if(fSyst[i_syst]->fIsShape){
            HistoTools::ManageHistograms( fSyst[i_syst]->fSmoothType + fSyst[i_syst]->fSymmetrisationType,  h_nominal, fSyst[i_syst]->fHistUp, fSyst[i_syst]->fHistDown, h_syst_up, h_syst_down);
        }
        
        //
        // save stuff
        //
        fSyst[i_syst]->fHistUp_original = (TH1*)fSyst[i_syst]->fHistUp->Clone();
        fSyst[i_syst]->fHistUp = h_syst_up;
        fSyst[i_syst]->fHistDown_original = (TH1*)fSyst[i_syst]->fHistDown->Clone();
        fSyst[i_syst]->fHistDown = h_syst_down;
        
        //normalisation component first
        if(h_nominal->Integral()!=0){
            fSyst[i_syst]->fNormUp = fSyst[i_syst]->fHistUp->Integral()/h_nominal->Integral() - 1.;
            fSyst[i_syst]->fNormDown = fSyst[i_syst]->fHistDown->Integral()/h_nominal->Integral() - 1.;
        } else {
            std::cerr << "<!> In SampleHist::SmoothSyst(): A nominal histogram with 0 intergral has been found. Please check ! " << std::endl;
            std::cerr << "            -> Sample: " << fName << std::endl;
        }
        
        if(fSyst[i_syst]->fIsShape){
            // update shape hists as well
            fSyst[i_syst]->fHistShapeUp = (TH1*)h_syst_up->Clone(fSyst[i_syst]->fHistShapeUp->GetName());
            fSyst[i_syst]->fHistShapeDown = (TH1*)h_syst_down->Clone(fSyst[i_syst]->fHistShapeDown->GetName());
            fSyst[i_syst]->fHistShapeUp->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeUp->Integral());
            fSyst[i_syst]->fHistShapeDown->Scale(fHist->Integral() / fSyst[i_syst]->fHistShapeDown->Integral());
        }
    }
    fSystSmoothed = true;
}
