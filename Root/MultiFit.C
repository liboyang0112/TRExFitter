#include "TtHFitter/MultiFit.h"

#include "TtHFitter/ConfigParser.h"

// -------------------------------------------------------------------------------------------------
// class MultiFit

//__________________________________________________________________________________
//
MultiFit::MultiFit(string name){
    fFitList.clear();
    fName = name;
    fLabel = name;
    fShowObserved = false;
    fLimitTitle = "95% CL limit on #sigma/#sigma_{SM}(t#bar{t}H) at m_{H}=125 GeV";
    fPOITitle = "best fit #mu=#sigma^{t#bar{t}H}/#sigma^{t#bar{t}H}_{SM} for m_{H}=125 GeV";
    fConfig = new ConfigParser();
    fSaveSuf = "";
    fFitShowObserved.clear();
}

//__________________________________________________________________________________
//
MultiFit::~MultiFit(){
    fFitList.clear();
}

//__________________________________________________________________________________
//
void MultiFit::ReadConfigFile(string configFile,string options){
    fConfig->ReadFile(configFile);
    ConfigSet *cs; // to store stuff later
    string param;
    //
    // Read options (to skip stuff, or include only some regions, samples, systs...)
    // ...
    
    //
    // set multi-fit
    cs = fConfig->GetConfigSet("MultiFit");
    fName = cs->GetValue();
    param = cs->Get("Label");
    if(param!="") fLabel = param;
    else          fLabel = fName;
    param = cs->Get("LumiLabel"); if( param != "")  fLumiLabel = param;
    param = cs->Get("CmeLabel");  if( param != "")  fCmeLabel  = param;
    param = cs->Get("SaveSuf");   if( param != "")  fSaveSuf   = param;
    param = cs->Get("ShowObserved");   if( param != "" && param != "FALSE") fShowObserved = true;
    param = cs->Get("LimitTitle"); if( param != "") fLimitTitle = param;
    if(fLimitTitle.find("95CL")!=string::npos) fLimitTitle.replace(fLimitTitle.find("95CL"),4,"95% CL");
    
    //
    // fits
    int nFit = 0;
    while(true){
        cs = fConfig->GetConfigSet("Fit",nFit);
        if(cs==0x0) break;
        nFit++;
        // options
        string fullOptions;
        param = cs->Get("Options");
        if(param!="" && options!="") fullOptions = options+";"+param;
        else if(param!="") fullOptions = param;
        else fullOptions = options;
        // label
        param = cs->Get("Label");
        string label = cs->GetValue();
        if(param!="") label = param;
        // load suf
        param = cs->Get("LoadSuf");
        string loadSuf = "";
        if(param!="") loadSuf = param;
        // config file
        string confFile = "";
        param = cs->Get("ConfigFile");
        if(param!="") confFile = param;
        // show obs
        param = cs->Get("ShowObserved");
        if(param=="FALSE") fFitShowObserved.push_back(false);
        else fFitShowObserved.push_back(true);
        //
        AddFitFromConfig(confFile,fullOptions,label,loadSuf);
    }
}

//__________________________________________________________________________________
//
void MultiFit::AddFitFromConfig(string configFile,string options,string label,string loadSuf){
    fFitList.push_back(new TtHFit());
    fFitList[fFitList.size()-1]->ReadConfigFile(configFile,options);
    fFitLabels.push_back(label);
    fFitSuffs.push_back(loadSuf);
}

//__________________________________________________________________________________
//
void MultiFit::ComparePOI(string POI){
    float xmax = 2;
    string process = fLabel;
    
    // Fit titles
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        std::cout << "Adding Fit: " << fFitList[i_fit]->fName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << std::endl;
        names.push_back( fFitList[i_fit]->fName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
    }
  
    int N = names.size();
  
    float ymin = -0.5;
    float ymax = N+1-0.5;
  
    TCanvas *c = new TCanvas("c","c",700,500);
    gStyle->SetEndErrorSize(6.);
    
    TGraph *g_central    = new TGraph(N);
    TGraphErrors *g_stat = new TGraphErrors(N);
    TGraphErrors *g_tot  = new TGraphErrors(N);
    
    int Ndiv = N+1;
  
    NuisParameter *par;
    bool found = false;
    
    // get values
    for(int i=0;i<N;i++){
//         fFitList[i]->ReadFitResults(names[i]+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0"+suffs[i]+".txt");
        fFitList[i]->ReadFitResults(names[i]+"/Fits/"+names[i]+suffs[i]+".txt");
        found = false;
        for(unsigned int j = 0; j<fFitList[i]->fFitResults->fNuisPar.size(); ++j){
            par = fFitList[i]->fFitResults->fNuisPar[j];
            if(par->fName == POI){
                g_central->SetPoint(i,par->fFitValue,i);
                g_stat->SetPoint(i,par->fFitValue,i);
                g_tot->SetPoint(i,par->fFitValue,i);
                g_stat->SetPointError(i,0,0);
                g_tot->SetPointError(i,par->fPostFitUp,0);
                if(par->fFitValue+par->fPostFitUp > xmax) xmax = par->fFitValue+par->fPostFitUp;
                found = true;
                break;
            }
        }
        if(!found){
            g_central->SetPoint(i,-10,i);
            g_stat->SetPoint(i,-10,i);
            g_tot->SetPoint(i,-10,i);
            g_stat->SetPointError(i,0,0);
            g_tot->SetPointError(i,0,0);
        }
    }
    
    g_stat->SetLineWidth(2);
    g_tot->SetLineWidth(3);
    g_stat->SetLineColor(kGreen);
    g_tot->SetLineColor(kBlack);
    g_central->SetMarkerStyle(kFullCircle);
    g_central->SetMarkerColor(kRed);
    g_central->SetMarkerSize(1.5);
    g_tot->SetMarkerSize(0);
    g_stat->SetMarkerSize(0);    
  
    xmax *= 2.5;
    
    TH1F* h_dummy = new TH1F("h_dummy","h_dummy",1,0,xmax);
    h_dummy->Draw();
    h_dummy->SetMinimum(ymin);
    h_dummy->SetMaximum(ymax);
    h_dummy->SetLineColor(kWhite);
    h_dummy->GetYaxis()->Set(N+1,ymin,ymax);
    h_dummy->GetYaxis()->SetNdivisions(Ndiv);
    
    TLatex *tex = new TLatex();
//     tex->SetNDC();

    for(int i=0;i<N;i++){
        h_dummy->GetYaxis()->SetBinLabel(i+1,titles[i].c_str());
//         myText(0.5,(1.*i)/(1.*N),kBlack,Form("#mu= %.1f",g_central->GetY()[i]));
//                 tex->DrawLatex(0.5,(1.*i)/(1.*N),Form("#mu= %.1f",g_central->GetY()[i]));
                tex->DrawLatex(0.5*xmax,i,Form("#mu= %.1f",g_central->GetX()[i]));
                tex->DrawLatex(0.7*xmax,i,Form("^{+%.1f}",g_tot->GetErrorXhigh(i)));
                tex->DrawLatex(0.7*xmax,i,Form("_{-%.1f}",g_tot->GetErrorXlow(i)));
                tex->DrawLatex(0.85*xmax,i,Form("^{+%.1f}",g_stat->GetErrorXhigh(i)));
                tex->DrawLatex(0.85*xmax,i,Form("_{-%.1f}",g_stat->GetErrorXlow(i)));
    }
    
    g_tot->Draw("E same");
    g_stat->Draw("E same");
    g_central->Draw("P same");

    TLine *l_SM = new TLine(1,-0.5,1,N-0.5);
    l_SM->SetLineWidth(2);
    l_SM->SetLineColor(kGray);
    l_SM->Draw("same");
    
    c->RedrawAxis();

    gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
    gPad->SetBottomMargin( 1.15*gPad->GetBottomMargin() );
    gPad->SetTopMargin( 1.8*gPad->GetTopMargin() );
    h_dummy->GetXaxis()->SetTitle(fPOITitle.c_str());

    ATLASLabel(0.02,0.93,"    Internal",kBlack);
    myText(0.35,0.93,kBlack,process.c_str());
    myText(0.65,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    
    TLegend *leg;
    leg = new TLegend(0.35,0.775,0.7,0.9);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(g_tot,"tot.","l");
    leg->AddEntry(g_stat,"stat.","l");
    leg->Draw();
    
    tex->DrawLatex((0.7-0.02)*xmax,N,"( tot )");
    tex->DrawLatex((0.85-0.02)*xmax,N,"( stat )");
    
//     myText(0.75,0.4,kBlack,"Stat. only");
    
    c->SaveAs( (fName+"/POI.png").c_str() ); 
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::CompareLimit(){
    float xmax = 2;
    string process = fLabel;
    gStyle->SetEndErrorSize(0.);
    
    // ---
    
    // Fit titles
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        std::cout << "Adding Fit: " << fFitList[i_fit]->fName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << std::endl;
        names.push_back( fFitList[i_fit]->fName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
    }

    // ---
    
    bool showObs = fShowObserved;
  
    int N = names.size();
  
    float ymin = -0.5;
    float ymax = N-0.5;
  
    TCanvas *c = new TCanvas("c","c",700,500);
  
    TGraphErrors *g_obs = new TGraphErrors(N);
    TGraphErrors *g_exp = new TGraphErrors(N);
    TGraphAsymmErrors *g_1s = new TGraphAsymmErrors(N);
    TGraphAsymmErrors *g_2s = new TGraphAsymmErrors(N);
  
    int Ndiv = N+1;
  
    TFile *f;
    TH1* h;
  
    // get values
    for(int i=0;i<N;i++){        
        f = new TFile(Form("%s/Limits/%s.root",names[i].c_str(),(names[i]+suffs[i]).c_str()) );
        std::cout << "Reading file " << Form("%s/Limits/%s.root",names[i].c_str(),(names[i]+suffs[i]).c_str()) << std::endl;
        h = (TH1*)f->Get("limit");
        
        std::cout << " " << h->GetBinContent(1) << std::endl;
        
        if(fFitShowObserved[i]) g_obs->SetPoint(i,h->GetBinContent(1),i);
        else g_obs->SetPoint(i,-1,i);
        g_exp->SetPoint(i,h->GetBinContent(2),i);
        g_1s->SetPoint(i,h->GetBinContent(2),i);
        g_2s->SetPoint(i,h->GetBinContent(2),i);
        g_obs->SetPointError(i,0,0.5);
        g_exp->SetPointError(i,0,0.5);
        g_1s->SetPointError(i,h->GetBinContent(2)-h->GetBinContent(5),h->GetBinContent(4)-h->GetBinContent(2),0.5,0.5);
        g_2s->SetPointError(i,h->GetBinContent(2)-h->GetBinContent(6),h->GetBinContent(3)-h->GetBinContent(2),0.5,0.5);
        
        if(h->GetBinContent(1)>xmax) xmax = h->GetBinContent(1); 
        if(h->GetBinContent(2)>xmax) xmax = h->GetBinContent(2); 
        if(h->GetBinContent(3)>xmax) xmax = h->GetBinContent(3); 
        if(h->GetBinContent(4)>xmax) xmax = h->GetBinContent(4); 
        if(h->GetBinContent(5)>xmax) xmax = h->GetBinContent(5); 
        if(h->GetBinContent(6)>xmax) xmax = h->GetBinContent(6);
    }
    
    g_obs->SetLineWidth(3);
    g_exp->SetLineWidth(3);
    g_exp->SetLineStyle(2);
    g_1s->SetFillColor(kGreen);
    g_1s->SetLineWidth(3);
    g_1s->SetLineStyle(2);
    g_2s->SetFillColor(kYellow);
//     g_2s->SetLineColor(kYellow);
    g_2s->SetLineWidth(3);
    g_2s->SetLineStyle(2);
    
    g_2s->SetMarkerSize(0);
    g_1s->SetMarkerSize(0);
    g_exp->SetMarkerSize(0);
    g_obs->SetMarkerSize(0);
  
//     xmax *= 2;
    
    TH1F* h_dummy = new TH1F("h_dummy","h_dummy",1,0,xmax);
    h_dummy->Draw();
    h_dummy->SetMinimum(ymin);
    h_dummy->SetMaximum(ymax);
    h_dummy->SetLineColor(kWhite);
    h_dummy->GetYaxis()->Set(N,ymin,ymax);
    h_dummy->GetYaxis()->SetNdivisions(Ndiv);
    for(int i=0;i<N;i++){
        h_dummy->GetYaxis()->SetBinLabel(i+1,titles[i].c_str());
    }
    
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
    h_dummy->GetXaxis()->SetTitle(fLimitTitle.c_str());

    ATLASLabel(0.02,0.93,"    Internal",kBlack);
    myText(0.35,0.93,kBlack,process.c_str());
    myText(0.65,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    
    TLegend *leg;
    if(showObs) leg = new TLegend(0.65,0.2,0.95,0.40);
    else        leg = new TLegend(0.65,0.2,0.95,0.35);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(g_1s,"Expected #pm 1#sigma","lf");
    leg->AddEntry(g_2s,"Expected #pm 2#sigma","lf");
    if(showObs) leg->AddEntry(g_obs,"Observed","l");
    leg->Draw();
    
//     myText(0.75,0.4,kBlack,"Stat. only");
    
    c->SaveAs( (fName+"/Limits.png").c_str() );
    delete c;
}

//__________________________________________________________________________________
//
void MultiFit::ComparePulls(){
    float ydist = 0.2;
    
    // Fit titles
    vector<string> names;  names.clear();
    vector<string> suffs;  suffs.clear();
    vector<string> titles; titles.clear();
    vector<float>  yshift; yshift.clear();
//     vector<int>    color;  color.clear();
//     vector<int>    style;  style.clear();
    
    int color[] = {kBlack,kRed,kBlue,kViolet};
    int style[] = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown};
    
    unsigned int N = fFitList.size();
    
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        names.push_back( fFitList[i_fit]->fName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
        yshift.push_back( 0. - ydist*N/2. + ydist*i_fit );
    }

    float xmin = -2.9;
    float xmax = 2.9;
    float max = 0;
    string npToExclude[] = {"SigXsecOverSM","gamma_","stat_"};
    bool brazilian = true;
    bool grayLines = false;
    
    // create a list of Systematics
    std::vector< string > Names;  Names.clear();
    std::vector< string > Titles; Titles.clear();
    string systName;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        for(unsigned int i_syst=0;i_syst<fFitList[i_fit]->fNSyst;i_syst++){
            systName = fFitList[i_fit]->fSystematics[i_syst]->fName;
            if(FindInStringVector(Names,systName)<0){
                Names.push_back(systName);
                Titles.push_back(fFitList[i_fit]->fSystematics[i_syst]->fTitle);
            }
        }
    }
    unsigned int Nsyst = Names.size();
    
    // read fit resutls
    NuisParameter *par;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
//         fFitList[i_fit]->ReadFitResults(names[i_fit]+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0"+suffs[i_fit]+".txt");      
        fFitList[i_fit]->ReadFitResults(names[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".txt");      
    }
    
    // exclude unused systematics
    std::vector<string> NamesNew; NamesNew.clear();
    std::vector<string> TitlesNew; TitlesNew.clear();
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            for(unsigned int j = 0; j<fFitList[i_fit]->fFitResults->fNuisPar.size(); ++j){
                par = fFitList[i_fit]->fFitResults->fNuisPar[j];
                systName = par->fName;
                if(systName==Names[i_syst]){
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
        if(found){
            NamesNew.push_back(Names[i_syst]);
            TitlesNew.push_back(Titles[i_syst]);
        }
    }
    Nsyst = NamesNew.size();
    Names.clear();
    Titles.clear();
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        Names.push_back(NamesNew[i_syst]);
        Titles.push_back(TitlesNew[i_syst]);
    }
    
    // fill stuff
    std::vector< TGraphAsymmErrors* > g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,float> centralMap; centralMap.clear();
        std::map<string,float> errUpMap;   errUpMap.clear();
        std::map<string,float> errDownMap; errDownMap.clear();
//         fFitList[i_fit]->ReadFitResults(names[i_fit]+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0"+suffs[i_fit]+".txt");
        for(unsigned int j = 0; j<fFitList[i_fit]->fFitResults->fNuisPar.size(); ++j){
            par = fFitList[i_fit]->fFitResults->fNuisPar[j];
            systName = par->fName;
            centralMap[systName] = par->fFitValue;
            errUpMap[systName]   = par->fPostFitUp;
            errDownMap[systName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.push_back( new TGraphAsymmErrors(Nsyst) );
        for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
            systName = Names[i_syst];
            if(centralMap[systName]!=0 || (errUpMap[systName]!=0 || errDownMap[systName]!=0)){
                g[i_fit]->SetPoint(i_syst,centralMap[systName],i_syst+0.5+yshift[i_fit]);
                g[i_fit]->SetPointEXhigh(i_syst,  errUpMap[systName]  <1 ?  errUpMap[systName]   : 1);
                g[i_fit]->SetPointEXlow( i_syst, -errDownMap[systName]<1 ? -errDownMap[systName] : 1);
            }
            else{
                g[i_fit]->SetPoint(i_syst,-10,-10);
                g[i_fit]->SetPointEXhigh(i_syst, 0);
                g[i_fit]->SetPointEXlow( i_syst, 0);
            }
        }
    }
    
    max = Nsyst;
    
    int lineHeight = 20;
//     int offsetUp = 10;
    int offsetUp = 50;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas *c = new TCanvas("c","c",600,newHeight);
    c->SetTicks(1,0);
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.33);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);
    
    TH1F *h_dummy = new TH1F("h_dummy","h_dummy",10,xmin,xmax);
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
    
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        g[i_fit]->SetLineColor(color[i_fit]);
        g[i_fit]->SetMarkerColor(color[i_fit]);
        g[i_fit]->SetMarkerStyle(style[i_fit]);  
        g[i_fit]->Draw("P same");
    }
    
    TLatex *systs = new TLatex();
    systs->SetTextSize( systs->GetTextSize()*0.8 );
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        systs->DrawLatex(3.,i_syst+0.25,Titles[i_syst].c_str());
    }
    h_dummy->GetXaxis()->SetLabelSize( h_dummy->GetXaxis()->GetLabelSize()*0.9 );

    TLegend *leg;
    leg = new TLegend(0.01,0.97,0.99,0.99);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetNColumns(4);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        leg->AddEntry(g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg->Draw();
    
    gPad->RedrawAxis();
    
    c->SaveAs((fName+"/NuisPar_comp"+fSaveSuf+".png").c_str());
    delete c;
}
