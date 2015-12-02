#include "TtHFitter/FitResults.h"

//__________________________________________________________________________________
//
FitResults::FitResults(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
    fCorrMatrix = 0;
    fNuisPar.clear();
}

//__________________________________________________________________________________
//
FitResults::~FitResults(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
    if(fCorrMatrix) delete fCorrMatrix;
    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        if(fNuisPar[i]) delete fNuisPar[i];
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
float FitResults::GetNuisParValue(string p){
    int idx = -1;
    if(fNuisParIsThere[p]){
        idx = fNuisParIdx[p];
    }
    else{
        if(TtHFitter::DEBUGLEVEL>0) cout << "  WARNING: NP " << p << " not found... Returning 0." << endl;
        return 0.;
    }
    return fNuisPar[idx]->fFitValue;
}

//__________________________________________________________________________________
//
float FitResults::GetNuisParErrUp(string p){
    int idx = -1;
    if(fNuisParIsThere[p]){
        idx = fNuisParIdx[p];
    }
    else{
        if(TtHFitter::DEBUGLEVEL>0) cout << "  WARNING: NP " << p << " not found... Returning error = 1." << endl;
        return 1.;
    }
    return fNuisPar[idx]->fPostFitUp;
}

//__________________________________________________________________________________
//
float FitResults::GetNuisParErrDown(string p){
    int idx = -1;
    if(fNuisParIsThere[p]){
        idx = fNuisParIdx[p];
    }
    else{
        if(TtHFitter::DEBUGLEVEL>0) cout << "  WARNING: NP " << p << " not found... Returning error = 1." << endl;
        return 1.;
    }
    return fNuisPar[idx]->fPostFitDown;
}

//__________________________________________________________________________________
//
void FitResults::ReadFromTXT(string fileName){
    bool includeCorrelations = true;
    bool invertedCorrMatrix = true;
    bool print = true;
    //
    CorrelationMatrix* matrix = new CorrelationMatrix();
    NuisParameter *np;
    //
    // get fitted NP's
    std::ifstream in;
    in.open(fileName.c_str());
    string input;
    string line;
    bool readingNP = false;
    bool readingCM = false;
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
            if(TtHFitter::DEBUGLEVEL>0){
                cout << "--------------------" << endl;
                cout << "Reading Nuisance Parameters..." << endl;
                cout << "--------------------" << endl;
            }
            readingNP = true;
            continue;
        }
        else if(line=="CORRELATION_MATRIX"){
            if(TtHFitter::DEBUGLEVEL>0){
                cout << "--------------------" << endl;
                cout << "Reading Correlation Matrix..." << endl;
                cout << "--------------------" << endl;
            }
            readingNP = false;
            readingCM = true;
            std::getline(in, line); // skip 1 line
            Nsyst_corr = atof(line.substr(0,line.find(" ")).c_str());
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
            iss >> value >> up >> down;
            np = fNuisPar[fNuisParIdx[name]];
            // set the title of this NP if there in the stored map
            //if(TtHFitter::SYSTMAP[name]!="") np->fTitle = TtHFitter::SYSTMAP[name];
            //
            np->fFitValue = value;
            np->fPostFitUp = up;
            np->fPostFitDown = down;
            if(TtHFitter::DEBUGLEVEL>0){
                if(print) cout << name << ": " << value << " +" << up << " " << down << endl;
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
    }
    if(includeCorrelations){
        if(TtHFitter::DEBUGLEVEL>0){
          if(print){
              for(int j_sys=0;j_sys<Nsyst_corr;j_sys++){
                  cout << "\t " << fNuisParNames[j_sys];
              }
              cout << endl;
              for(int i_sys=0;i_sys<Nsyst_corr;i_sys++){
                  cout << fNuisParNames[i_sys];
                  for(int j_sys=0;j_sys<Nsyst_corr;j_sys++){
                      cout << Form("\t%.4f",matrix->GetCorrelation(fNuisParNames[i_sys],fNuisParNames[j_sys]));
                  }
                  cout << endl;
              }
          }
        }
    }
    fCorrMatrix = matrix;
    //
    int TOTsyst = fNuisParNames.size();
    cout << "Found " << TOTsyst << " systematics." << endl;
    if(TOTsyst<=0) cout << "WARNING: No systematics found in fit result file..." << endl;
}

//__________________________________________________________________________________
//
void FitResults::DrawPulls(string path, string category){
    float xmin = -2.9;
    float xmax = 2.9;
    float max = 0;
    string npToExclude[] = {"SigXsecOverSM","gamma_","stat_"};
    bool brazilian = true;
    bool grayLines = false;
    
    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    
    NuisParameter *par;
    int idx = 0;
    std::vector< string > Names;
    Names.clear();
    
    for(unsigned int i = 0; i<fNuisPar.size(); ++i){
        par = fNuisPar[i];
        
        if( category != "all" && category != par->fCategory ) continue;
        
        bool skip = false;
        for(int ii=0; ii<sizeof(npToExclude)/sizeof(string); ii++){
            if(par->fName.find(npToExclude[ii])!=string::npos){
                skip = true;
                continue;
            }
        }
        if(skip) continue;
        
        g->SetPoint(idx,par->fFitValue,idx+0.5);
        g->SetPointEXhigh(idx, par->fPostFitUp);
        g->SetPointEXlow( idx,-par->fPostFitDown);
        
        Names.push_back(par->fTitle);
    
        idx ++;
        if(idx > max)  max = idx;      
    }

    int lineHeight = 20;
    int offsetUp = 40;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
//     TCanvas *c = new TCanvas("c","c",600,newHeight);
    TCanvas *c = new TCanvas("c","c",800,newHeight);
    c->SetTicks(1,0);
//     gPad->SetLeftMargin(0.05);
    gPad->SetLeftMargin(0.05/(8./6.));
//     gPad->SetRightMargin(0.33);
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);
    
    TH1F *h_dummy = new TH1F( ("h_dummy"+category).c_str(),("h_dummy"+category).c_str(),10,xmin,xmax);
    h_dummy->SetMaximum(max);
    h_dummy->SetLineWidth(0);
    h_dummy->SetFillStyle(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->SetFillColor(kWhite);
    h_dummy->SetMinimum(0.);
    h_dummy->GetYaxis()->SetLabelSize(0);
    h_dummy->Draw();
    h_dummy->GetYaxis()->SetNdivisions(0);

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
    
    g->Draw("psame");
    
    TLatex *systs = new TLatex();
    systs->SetTextSize( systs->GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs->DrawLatex(3.,i+0.25,Names[i].c_str());
    }
    h_dummy->GetXaxis()->SetLabelSize( h_dummy->GetXaxis()->GetLabelSize()*0.9 );
    
    gPad->RedrawAxis();
    
    if(category!="all"){
        TLatex *cat_legend = new TLatex();
//         cat_legend -> DrawLatexNDC(0.5,0.8,category.c_str());
        cat_legend -> DrawLatexNDC(0.5,1-0.8*offsetUp/newHeight,category.c_str());
    }
    
    c->SaveAs(path.c_str());
    delete c;
}

//__________________________________________________________________________________
//
void FitResults::DrawCorrelationMatrix(string path, const double corrMin){
    if(fCorrMatrix){
        fCorrMatrix->Draw(path, corrMin);
    }
}
