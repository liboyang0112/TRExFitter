// Class include
#include "TRExFitter/FitResults.h"

// framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/StatusLogbook.h"

// ROOT includes
#include "TBox.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"

//c++ includes
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//__________________________________________________________________________________
//
FitResults::FitResults(){
    fCorrMatrix = nullptr;
    fNLL = 0;
}

//__________________________________________________________________________________
//
FitResults::~FitResults(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
    delete fCorrMatrix;
    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        delete fNuisPar[i];
    }
    fNuisPar.clear();
}

//__________________________________________________________________________________
//
void FitResults::AddNuisPar(NuisParameter *par){
    fNuisPar.push_back(par);
    string p = par->fName;
    fNuisParIdx[p] = (int)fNuisParNames.size();
    fNuisParNames.push_back(p);
    fNuisParIsThere[p] = true;
}

//__________________________________________________________________________________
//
float FitResults::GetNuisParValue(const string& p){
    int idx = -1;
    if(fNuisParIsThere[p]){
        idx = fNuisParIdx[p];
    }
    else{
        return 0.;
    }
    return fNuisPar[idx]->fFitValue;
}

//__________________________________________________________________________________
//
float FitResults::GetNuisParErrUp(const std::string& p){
    int idx = -1;
    if(fNuisParIsThere[p]){
        idx = fNuisParIdx[p];
    }
    else{
        return 1.;
    }
    return fNuisPar[idx]->fPostFitUp;
}

//__________________________________________________________________________________
//
float FitResults::GetNuisParErrDown(const std::string& p){
    int idx = -1;
    if(fNuisParIsThere[p]){
        idx = fNuisParIdx[p];
    }
    else{
        return -1.;
    }
    return fNuisPar[idx]->fPostFitDown;
}

//__________________________________________________________________________________
//
void FitResults::ReadFromTXT(const std::string& fileName, const std::vector<std::string>& blinded){
    bool includeCorrelations = true;
    bool invertedCorrMatrix = true;
    bool print = true;
    //
    CorrelationMatrix *matrix = new CorrelationMatrix();
    //
    // get fitted NP's
    std::ifstream in;
    in.open(fileName.c_str());

    if (!in.is_open())	{
        delete matrix;
        WriteErrorStatus("FitResults::ReadFromTXT","Could not open the file \"" + fileName + "\"");
        return;
    }

    string input;
    string line;
    bool readingNP = false;
    bool readingCM = false;
    bool readingNLL = false;
    int i = 0;
    int j = 0;
    int Nsyst_corr = 0;
    //
    string name;
    float value, up, down;
    float corr;
    //
    // read file line by line
    while(std::getline(in, line)){
        if(line=="") continue;
        if(line=="NUISANCE_PARAMETERS"){
            WriteDebugStatus("FitResults::ReadFromTXT", "--------------------");
            WriteDebugStatus("FitResults::ReadFromTXT", "Reading Nuisance Parameters...");
            WriteDebugStatus("FitResults::ReadFromTXT", "--------------------");
            readingNP = true;
            continue;
        }
        else if(line=="CORRELATION_MATRIX"){
            WriteDebugStatus("FitResults::ReadFromTXT", "--------------------");
            WriteDebugStatus("FitResults::ReadFromTXT", "Reading Correlation Matrix...");
            WriteDebugStatus("FitResults::ReadFromTXT", "--------------------");
            readingNP = false;
            readingCM = true;
            std::getline(in, line); // skip 1 line
            Nsyst_corr = atof(line.substr(0,line.find(" ")).c_str());
            // add all NPs to corr matrix here (to keep decent order)
            for(auto npName : fNuisParNames) matrix->AddNuisPar(npName);
            continue;
        }
        else if(line=="NLL"){
            WriteDebugStatus("FitResults::ReadFromTXT", "--------------------");
            WriteDebugStatus("FitResults::ReadFromTXT", "Reading Negative Log-Likelihood (NLL) value...");
            WriteDebugStatus("FitResults::ReadFromTXT", "--------------------");
            readingNP = false;
            readingCM = false;
            readingNLL = true;
            continue;
        }
        std::istringstream iss(line);
        if(readingNP){
            iss >> input;
            if(input=="" || input=="CORRELATION_MATRIX"){
                readingNP = false;
            }
            while(input.find("\\")!=string::npos) input = input.replace(input.find("\\"),1,"");
            name = input;
            // clean the syst name...
            name = ReplaceString(name,"alpha_","");
            name = ReplaceString(name,"gamma_","");
            AddNuisPar(new NuisParameter(name));
            if (std::find(blinded.begin(), blinded.end(), name) == blinded.end()){
                iss >> value >> up >> down;
                NuisParameter *np = fNuisPar[fNuisParIdx[name]];
                np->fFitValue = value;
                np->fPostFitUp = up;
                np->fPostFitDown = down;
                if(print) WriteVerboseStatus("FitResults::ReadFromTXT", name + ": " + std::to_string(value) + " +" + std::to_string(up) + " " + std::to_string(down));
            } else {
                std::string hex;
                iss >> hex >> up >> down;
                NuisParameter *np = fNuisPar[fNuisParIdx[name]];
                np->fFitValue = HexToFloat(hex);
                np->fPostFitUp = up;
                np->fPostFitDown = down;
            }
            i++;
        }
        if(readingCM){
            if(!includeCorrelations) break;
            for(int i_sys=0;i_sys<Nsyst_corr;i_sys++){
                iss >> corr;
                if(invertedCorrMatrix){
                    matrix->SetCorrelation(fNuisParNames[Nsyst_corr-i_sys-1],fNuisParNames[j],corr);
                }
                else matrix->SetCorrelation(fNuisParNames[i_sys],fNuisParNames[j],corr);
            }
            j++;
        }
        if(readingNLL){
            iss >> fNLL;
        }
    }
    if(includeCorrelations){
        if(print){
            std::string temp_string = "";
            for(int j_sys=0;j_sys<Nsyst_corr;j_sys++){
                temp_string+= "\t " + fNuisParNames[j_sys];
            }
            WriteVerboseStatus("FitResults::ReadFromTXT",temp_string);
            temp_string = "";
            for(int i_sys=0;i_sys<Nsyst_corr;i_sys++){
                temp_string +=  fNuisParNames[i_sys];
                for(int j_sys=0;j_sys<Nsyst_corr;j_sys++){
                    temp_string += Form("\t%.4f",matrix->GetCorrelation(fNuisParNames[i_sys],fNuisParNames[j_sys]));
                }
                WriteVerboseStatus("FitResults::ReadFromTXT",temp_string);
            }
        }
    }
    fCorrMatrix = matrix;
    //
    int TOTsyst = fNuisParNames.size();
    WriteDebugStatus("FitResults::ReadFromTXT", "Found " + std::to_string(TOTsyst) + " systematics.");
    if (TOTsyst<=0) WriteDebugStatus("FitResults::ReadFromTXT", "No systematics found in fit result file. Stat-only fit-results?");
    WriteDebugStatus("FitResults::ReadFromTXT", "Negative Log-Likelihood value NLL = " + std::to_string(fNLL));
}

//__________________________________________________________________________________
//
void FitResults::DrawNormFactors( const string &path,
                                  const std::vector < NormFactor* > &normFactors, const std::vector<string>& blinded ) const {
    float xmin = 1000;
    float xmax = -1000;
    float max = 0;

    TGraphAsymmErrors g{};

    NuisParameter *par;
    std::vector< NuisParameter* > selected_norm_factors;

    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        par = fNuisPar[i];

        //skip the blinded NPs
        if (std::find(blinded.begin(), blinded.end(), par->fName) != blinded.end()) continue;

        bool isNorm = false;
        for( const auto *norm : normFactors ){
            if(norm->fName==par->fName){
                isNorm = true;
                break;
            }
        }
        if ( !isNorm ) continue;
        g.SetPoint(selected_norm_factors.size(),par->fFitValue,2*selected_norm_factors.size()+1);
        g.SetPointEXhigh(selected_norm_factors.size(), par->fPostFitUp);
        g.SetPointEXlow( selected_norm_factors.size(),-par->fPostFitDown);

        if( par->fFitValue+par->fPostFitUp > xmax ) xmax = par->fFitValue+par->fPostFitUp;
        if( par->fFitValue+par->fPostFitDown < xmin ) xmin = par->fFitValue+par->fPostFitDown;

        NuisParameter *nuis = new NuisParameter(par->fName);
        nuis->fFitValue =    par -> fFitValue;
        nuis->fPostFitUp =   par -> fPostFitUp;
        nuis->fPostFitDown = par -> fPostFitDown;
        nuis->fTitle =       par -> fTitle;
        selected_norm_factors.push_back(nuis);
        if(2*selected_norm_factors.size() > max)  max = 2*selected_norm_factors.size();
    }
    xmax *= (xmax<0 ? 0.5 : 1.5);
    xmin *= (xmin>0 ? 0.5 : 1.5);
    if(xmin>0) xmin = 0.;
    xmax += (xmax-xmin)*0.25;

    int lineHeight = 40;
    int offsetUp = 40;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy( "h_dummy_norm","h_dummy_norm",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0;
    TBox b1, b2;
    l0 = TLine(1,0,1,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    g.Draw("psame");

    TLatex systs{};
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(unsigned int i=0;i<selected_norm_factors.size();i++){
      systs.DrawLatex(xmax*1.05,2*i+0.75,(selected_norm_factors[i]->fTitle).c_str());
      systs.DrawLatex(xmax*0.7,2*i+0.75,
        Form(("%."+std::to_string(fPOIPrecision)+"f ^{%."+std::to_string(fPOIPrecision)+"f}_{%."+std::to_string(fPOIPrecision)+"f}").c_str(),selected_norm_factors[i]->fFitValue, selected_norm_factors[i]->fPostFitUp, selected_norm_factors[i]->fPostFitDown ) );
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    gPad->RedrawAxis();

    c.SaveAs(path.c_str());
}

//__________________________________________________________________________________
//
void FitResults::DrawGammaPulls( const string &path, const std::vector<std::string>& blinded ) const {
    float xmin = 10;
    float xmax = -10;
    float max = 0;

    TGraphAsymmErrors g{};

    NuisParameter *par;
    int idx = 0;
    std::vector< string > names;

    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        par = fNuisPar[i];

        std::string name = par->fName;
        name = ReplaceString(name,"gamma_","");
        if (std::find(blinded.begin(), blinded.end(), name) != blinded.end()) continue;

        if ( par->fName.find("stat_") == std::string::npos && par->fName.find("shape_") == std::string::npos ) continue;
        g.SetPoint(idx,par->fFitValue,idx+0.5);
        g.SetPointEXhigh(idx, par->fPostFitUp);
        g.SetPointEXlow( idx,-par->fPostFitDown);
        if( par->fFitValue+par->fPostFitUp > xmax ) xmax = par->fFitValue+par->fPostFitUp;
        if( par->fFitValue+par->fPostFitDown < xmin ) xmin = par->fFitValue+par->fPostFitDown;

        std::string clean_name = par->fTitle;
        clean_name = ReplaceString( clean_name, "stat_", "#gamma " );
        clean_name = ReplaceString( clean_name, "shape_", "#gamma " );
        clean_name = ReplaceString( clean_name, "#gamma #gamma ", "#gamma " );
        clean_name = ReplaceString( clean_name, "_", " " );
        names.push_back(clean_name);
        idx ++;
        if(idx > max)  max = idx;
    }
    xmax *= (1.2-(xmax<0));
    xmin *= (0.8+(xmin<0));

    int lineHeight = 20;
    int offsetUp = 40;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy( "h_dummy_gamma","h_dummy_gamma",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0;
    TBox b1, b2;
    l0 = TLine(1,0,1,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    g.Draw("psame");

    TLatex systs{};
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs.DrawLatex(xmax*1.05,i+0.25,names[i].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    gPad->RedrawAxis();

    c.SaveAs(path.c_str());
}

//__________________________________________________________________________________
//
void FitResults::DrawNPPulls( const string &path, const string &category, const std::vector < NormFactor* > &normFactors, const std::vector<std::string>& blinded ) const {
    float xmin = -2.9;
    float xmax = 2.9;
    float max = 0;
    std::vector<std::string> npToExclude = {"gamma_","stat_","shape_"};
    bool brazilian = true;

    TGraphAsymmErrors g{};

    NuisParameter *par = nullptr;
    int idx = 0;
    std::vector< string > names;

    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        par = fNuisPar[i];

        std::string name = par->fName;
        name = ReplaceString(name,"alpha_","");

        if (std::find(blinded.begin(), blinded.end(), name) != blinded.end()) continue;

        if( category != "all" && category != par->fCategory ) continue;
        if( FindInStringVector(fNuisParToHide,par->fName)>=0 ) continue;

        bool skip = false;
        for(const std::string& ii : npToExclude){
            if(par->fName.find(ii)!=string::npos){
                skip = true;
                break;
            }
        }
        for( const auto *norm : normFactors ){
          if(norm->fName==par->fName){
            skip = true;
            break;
          }
        }
        if(skip) continue;

        g.SetPoint(idx,par->fFitValue,idx+0.5);
        g.SetPointEXhigh(idx, par->fPostFitUp);
        g.SetPointEXlow( idx,-par->fPostFitDown);

        names.push_back(par->fTitle);

        idx ++;
        if(idx > max)  max = idx;
    }

    int lineHeight = 20;
    int offsetUp = 40;
    int offsetDown = 60;
    if (max < 10){
        offsetDown = 65;
    }
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy( ("h_dummy"+category).c_str(),("h_dummy"+category).c_str(),10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0;
    TBox b1, b2;
    if(brazilian){
        l0 = TLine(0,0,0,max);
        l0.SetLineStyle(7);
        l0.SetLineColor(kBlack);
        b1 = TBox(-1,0,1,max);
        b2 = TBox(-2,0,2,max);
        b1.SetFillColor(kGreen);
        b2.SetFillColor(kYellow);
        b2.Draw("same");
        b1.Draw("same");
        l0.Draw("same");
    }

    g.Draw("psame");

    TLatex systs{};
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs.DrawLatex(3.,i+0.25,names[i].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    h_dummy.GetXaxis()->CenterTitle();
    h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    if (max < 10){
        h_dummy.GetXaxis()->SetTitleOffset(0.9);
    } else {
        h_dummy.GetXaxis()->SetTitleOffset(1.15);
    }

    gPad->RedrawAxis();

    if(category!="all"){
        TLatex cat_legend{};
        cat_legend.DrawLatexNDC(0.5,1-0.8*offsetUp/newHeight,category.c_str());
    }

    c.SaveAs(path.c_str());
}

//__________________________________________________________________________________
//
void FitResults::DrawCorrelationMatrix(const std::string& path, const bool& useGammas, const double corrMin){
    if(fCorrMatrix){
        fCorrMatrix->fNuisParToHide = fNuisParToHide;
        fCorrMatrix->Draw(path, useGammas, corrMin);
    }
}
