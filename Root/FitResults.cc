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

//__________________________________________________________________________________
//
FitResults::FitResults() :
    fCorrMatrix(nullptr),
    fOutFolder(""),
    fPOIPrecision(2),
    fNLL() {
}

//__________________________________________________________________________________
//
FitResults::~FitResults(){
}

//__________________________________________________________________________________
//
void FitResults::AddNuisPar(NuisParameter *par){
    fNuisPar.emplace_back(par);
    const std::string p = par->fName;
    fNuisParIdx[p] = fNuisParNames.size();
    fNuisParNames.push_back(p);
    fNuisParIsThere[p] = true;
}

//__________________________________________________________________________________
//
double FitResults::GetNuisParValue(const std::string& p){
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
double FitResults::GetNuisParErrUp(const std::string& p){
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
double FitResults::GetNuisParErrDown(const std::string& p){
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
    std::unique_ptr<CorrelationMatrix> matrix(new CorrelationMatrix());
    //
    // get fitted NP's
    std::ifstream in;
    in.open(fileName.c_str());

    if (!in.is_open()) {
        WriteErrorStatus("FitResults::ReadFromTXT","Could not open the file \"" + fileName + "\"");
        return;
    }

    std::string input;
    std::string line;
    bool readingNP = false;
    bool readingCM = false;
    bool readingNLL = false;
    int i = 0;
    int j = 0;
    int Nsyst_corr = 0;
    //
    std::string name;
    double value, up, down;
    double corr;
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
            while(input.find("\\")!=std::string::npos) input = input.replace(input.find("\\"),1,"");
            name = input;
            // clean the syst name...
            name = Common::ReplaceString(name,"alpha_","");
            name = Common::ReplaceString(name,"gamma_","");
            AddNuisPar(new NuisParameter(name));
            if (std::find(blinded.begin(), blinded.end(), name) == blinded.end()){
                iss >> value >> up >> down;
                NuisParameter *np = fNuisPar[fNuisParIdx[name]].get();
                np->fFitValue = value;
                np->fPostFitUp = up;
                np->fPostFitDown = down;
                WriteVerboseStatus("FitResults::ReadFromTXT", name + ": " + std::to_string(value) + " +" + std::to_string(up) + " " + std::to_string(down));
            } else {
                std::string hex;
                iss >> hex >> up >> down;
                NuisParameter *np = fNuisPar[fNuisParIdx[name]].get();
                np->fFitValue = Common::HexToDouble(hex);
                np->fPostFitUp = up;
                np->fPostFitDown = down;
            }
            i++;
        }
        if(readingCM){
            matrix->Resize(Nsyst_corr);
            for(int i_sys=0;i_sys<Nsyst_corr;i_sys++){
                iss >> corr;
                matrix->SetCorrelation(fNuisParNames[Nsyst_corr-i_sys-1],fNuisParNames[j],corr);
            }
            j++;
        }
        if(readingNLL){
            iss >> fNLL;
        }
    }
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
    fCorrMatrix = std::unique_ptr<CorrelationMatrix>(matrix.release());
    //
    int TOTsyst = fNuisParNames.size();
    WriteDebugStatus("FitResults::ReadFromTXT", "Found " + std::to_string(TOTsyst) + " systematics.");
    if (TOTsyst<=0) WriteDebugStatus("FitResults::ReadFromTXT", "No systematics found in fit result file. Stat-only fit-results?");
    WriteDebugStatus("FitResults::ReadFromTXT", "Negative Log-Likelihood value NLL = " + std::to_string(fNLL));
}

//__________________________________________________________________________________
//
void FitResults::DrawNormFactors( const std::string &path,
                                  const std::vector < std::shared_ptr<NormFactor> > &normFactors, const std::vector<std::string>& blinded ) const {
    double xmin = 1000.;
    double xmax = -1000.;
    double max = 0.;

    TGraphAsymmErrors g{};

    std::vector< std::unique_ptr<NuisParameter> > selected_norm_factors;

    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        NuisParameter* par = fNuisPar[i].get();

        //skip the blinded NPs
        if (std::find(blinded.begin(), blinded.end(), par->fName) != blinded.end()) continue;

        bool isNorm = false;
        for( const auto& norm : normFactors ){
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

        std::unique_ptr<NuisParameter> nuis(new NuisParameter(par->fName));
        nuis->fFitValue =    par -> fFitValue;
        nuis->fPostFitUp =   par -> fPostFitUp;
        nuis->fPostFitDown = par -> fPostFitDown;
        nuis->fTitle =       par -> fTitle;
        selected_norm_factors.emplace_back(nuis.release());
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
void FitResults::DrawGammaPulls( const std::string & path, const std::vector<std::string> & blinded ) const {

    double xmin =  10.;
    double xmax = -10.;
    double max  =   0.;

    TGraphAsymmErrors g{};
      
    // get a copy of the NPs (I want an actual copy of the NP objects that I can
    // manipulate without changing the originals.  
    std::vector<NuisParameter> myNPs;
    
    // make the copies, dropping non-gamma NPs and blinded NPs
    for( unsigned i = 0; i < fNuisPar.size(); ++i )
    {
        std::string name = fNuisPar[i]->fName;
        name = Common::ReplaceString(name,"gamma_","");
        if (std::find(blinded.begin(), blinded.end(), name) != blinded.end()) continue;
        if ( fNuisPar[i]->fName.find("stat_") == std::string::npos && fNuisPar[i]->fName.find("shape_") == std::string::npos ) continue;
        
        // make a copy, clean the name and save it
        NuisParameter theNP( (*fNuisPar[i]) );
        
        std::string clean_name = theNP.fTitle;
        clean_name = Common::ReplaceString( clean_name, "stat_", "#gamma " );
        clean_name = Common::ReplaceString( clean_name, "shape_", "#gamma " );
        clean_name = Common::ReplaceString( clean_name, "#gamma #gamma ", "#gamma " );
        clean_name = Common::ReplaceString( clean_name, "_", " " );
        
        
        theNP.fTitle = Common::pad_trail( clean_name );
        
        myNPs.push_back( theNP );
  
    }
    
    // now sort myNPs by the cleaned and padded names
    // do a to_lower in the sort_func
    // sort using a custom struct
    struct {
        bool operator()(const NuisParameter & a, const NuisParameter & b) const
        {   
            return a.fTitle < b.fTitle;
        }   
    } my_sort_func;
        
    std::sort( myNPs.begin(), myNPs.end(), my_sort_func );
    
    // start making the plot
    for(unsigned int i = 0; i < myNPs.size(); ++i )
    {
        g.SetPoint(i, myNPs[i].fFitValue, i+0.5);
        g.SetPointEXhigh( i,  myNPs[i].fPostFitUp );
        g.SetPointEXlow ( i, -myNPs[i].fPostFitDown );
        if( (myNPs[i].fFitValue + myNPs[i].fPostFitUp) > xmax ) xmax = myNPs[i].fFitValue + myNPs[i].fPostFitUp;
        if( (myNPs[i].fFitValue + myNPs[i].fPostFitDown) < xmin ) xmin = myNPs[i].fFitValue + myNPs[i].fPostFitDown;
    }
    
    max   = myNPs.size();
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
    
    for(unsigned int i=0; i< myNPs.size(); ++i )
    {
        systs.DrawLatex( xmax*1.05, i+0.25, myNPs[i].fTitle.c_str() );
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    gPad->RedrawAxis();

    c.SaveAs(path.c_str());
}

//__________________________________________________________________________________
//
void FitResults::DrawNPPulls( const std::string &path, const std::string &category, const std::vector < std::shared_ptr<NormFactor> > &normFactors, const std::vector<std::string>& blinded ) const {
    double xmin = -2.9;
    double xmax = 2.9;
    double max = 0.;
    std::vector<std::string> npToExclude = {"gamma_","stat_","shape_"};
    
    // reorder the NPs
    std::vector < std::shared_ptr<NuisParameter> > nuisPar;
    if(fNuisParList.size()>0){
        for(auto& npName : fNuisParList){
            for(auto& np : fNuisPar){
                if(np->fName==npName) nuisPar.emplace_back(np);
            }
        }
    }
    else{
        for(auto& np : fNuisPar){
            nuisPar.emplace_back(np);
        }
    }

    TGraphAsymmErrors g{};

    NuisParameter *par = nullptr;
    int idx = 0;
    std::vector< std::string > names;

    for(unsigned int i = 0; i<nuisPar.size(); ++i){
        par = nuisPar[i].get();

        std::string name = par->fName;
        name = Common::ReplaceString(name,"alpha_","");

        if (std::find(blinded.begin(), blinded.end(), name) != blinded.end()) continue;

        if( category != "all" && category != par->fCategory ) continue;
        if( Common::FindInStringVector(fNuisParToHide,par->fName)>=0 ) continue;

        bool skip = false;
        for(const std::string& ii : npToExclude){
            if(par->fName.find(ii)!=std::string::npos){
                skip = true;
                break;
            }
        }
        for( const auto& norm : normFactors ){
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

    TLine l0(0,0,0,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    TBox b1(-1,0,1,max);
    TBox b2(-2,0,2,max);
    b1.SetFillColor(kGreen);
    b2.SetFillColor(kYellow);
    b2.Draw("same");
    b1.Draw("same");
    l0.Draw("same");

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
void FitResults::DrawCorrelationMatrix(const std::vector<std::string>& path, const bool& useGammas, const bool useHEPDataFormat, const double corrMin) {
    if(fCorrMatrix){
        fCorrMatrix->fOutFolder = fOutFolder;
        fCorrMatrix->fNuisParToHide = fNuisParToHide;
        fCorrMatrix->fNuisParList = fNuisParList;
        fCorrMatrix->Draw(path, useGammas, useHEPDataFormat, corrMin);
    }
}
