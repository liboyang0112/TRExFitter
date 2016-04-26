#include "TtHFitter/MultiFit.h"

#include "TtHFitter/ConfigParser.h"
#include "TtHFitter/FittingTool.h"

//Roofit headers
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"

//HistFactory headers
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"

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
    fPOI = "SigXsecOverSM";
    fPOIMin = 0;
    fPOIMax = 10;
    //
    fCombine       = false;
    fCompare       = false;
    fCompareLimits = true;
    fComparePOI    = true;
    fComparePulls  = true;
    //
    fDataName      = "obsData";
    fFitType       = 1; // 1: S+B, 2: B-only
    //
    fCombineChByCh = false;
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
    param = cs->Get("ShowObserved");   if( param != "" && param != "FALSE" ) fShowObserved = true;
    param = cs->Get("LimitTitle"); if( param != "") fLimitTitle = param;
    if(fLimitTitle.find("95CL")!=string::npos) fLimitTitle.replace(fLimitTitle.find("95CL"),4,"95% CL");
    param = cs->Get("CompareLimits"); if( param != "" && param != "TRUE" )  fCompareLimits = false;
    param = cs->Get("ComparePOI");    if( param != "" && param != "TRUE" )  fComparePOI    = false;
    param = cs->Get("ComparePulls");  if( param != "" && param != "TRUE" )  fComparePulls  = false;
    //
    param = cs->Get("Combine"); if( param != "" && param != "FALSE" )  fCombine = true;
    param = cs->Get("Compare"); if( param != "" && param != "FALSE" )  fCompare = true;
    //
    param = cs->Get("POIName");  if( param != "" ) fPOI = param;
    param = cs->Get("POIRange"); if( param != "" && Vectorize(param,',').size()==2 ) {
        fPOIMin = atof( Vectorize(param,',')[0].c_str() );
        fPOIMax = atof( Vectorize(param,',')[1].c_str() );
    }
    param = cs->Get("DataName"); if( param != "" ) fDataName = param;
    param = cs->Get("FitType");  if( param != "" ){
        if(param=="SPLUSB") fFitType = 1;
        if(param=="BONLY")  fFitType = 2;
    }
    //
    param = cs->Get("CombineChByCh"); if( param != "" && param != "FALSE" )  fCombineChByCh = true;
    
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
        // name
        fFitNames.push_back(cs->GetValue());
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
        // workspace
        string wsFile = "";
        param = cs->Get("Workspace");
        if(param!="") wsFile = param;
        // show obs
        param = cs->Get("ShowObserved");
        if(param=="FALSE") fFitShowObserved.push_back(false);
        else fFitShowObserved.push_back(true);
        //
        AddFitFromConfig(confFile,fullOptions,label,loadSuf,wsFile);
    }
    
    // make directory
    gSystem->mkdir(fName.c_str());
}

//__________________________________________________________________________________
//
void MultiFit::AddFitFromConfig(string configFile,string options,string label,string loadSuf,string wsFile){
    fFitList.push_back(new TtHFit());
    fFitList[fFitList.size()-1]->ReadConfigFile(configFile,options);
    fFitLabels.push_back(label);
    fFitSuffs.push_back(loadSuf);
    fWsFiles.push_back(wsFile);
}

//__________________________________________________________________________________
//
RooWorkspace* MultiFit::CombineWS(){
    cout << "...................................." << endl;
    std::cout << "Combining workspaces..." << std::endl;
    
    std::vector < RooWorkspace* > vec_ws;
    std::vector < std::string > vec_chName;
    RooStats::HistFactory::Measurement *measurement = 0x0;
    TFile *rootFileCombined = 0x0;
    
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        std::string fitName = fFitList[i_fit]->fName;
        std::cout << "Adding Fit: " << fitName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << std::endl;
        
        RooStats::HistFactory::Measurement *meas;
        std::string fileName = fitName + "/RooStats/" + fitName + "_combined_" + fitName + fFitSuffs[i_fit] + "_model.root";
        if(fWsFiles[i_fit]!="") fileName = fWsFiles[i_fit];
        std::cout << "Opening file " << fileName << std::endl;
        TFile *rootFile = new TFile(fileName.c_str(),"read");
        RooWorkspace* m_ws = (RooWorkspace*) rootFile->Get("combined");
        std::cout << "Getting " << fitName+fFitSuffs[i_fit] << std::endl;
        meas = (RooStats::HistFactory::Measurement*) rootFile -> Get( (fitName+fFitSuffs[i_fit]).c_str());
        //
        // import measurement if not there yet
        if(!measurement){
            measurement = meas;
        }
        
        if(!fCombineChByCh){
            //
            // Combine combined workspaces directly
            std::vector<RooStats::HistFactory::Channel> chVec = meas->GetChannels();
            for(unsigned int i_ch=0;i_ch<chVec.size();i_ch++){
                vec_ws.push_back(m_ws);
                vec_chName.push_back(chVec[i_ch].GetName());
            }
        }
        
        // 
        // Alternative way: combine the individual workspaces for the different chanenels
        // Loop on all the regions in each fit
        if(fCombineChByCh){
            for(unsigned int i_reg=0;i_reg<fFitList[i_fit]->fRegions.size();i_reg++){
                Region *reg = fFitList[i_fit]->fRegions[i_reg];
                std::string fileName = fitName + "/RooStats/" + fitName + "_" + reg->fName + "_" + fitName + fFitSuffs[i_fit] + "_model.root";
                std::cout << "  Opening file " << fileName << std::endl;
                TFile *rootFile = new TFile(fileName.c_str(),"read");
                RooWorkspace* m_ws = (RooWorkspace*) rootFile->Get(reg->fName.c_str());
                std::cout << "  Getting " << reg->fName << std::endl;
                vec_ws.push_back(m_ws);
                vec_chName.push_back(reg->fName);
            }
        }
        //
    }
    
    //
    // Create the HistoToWorkspaceFactoryFast object to perform safely the combination
    //
    if(!measurement){
        std::cout << "<!> Error in MultiFit::CombineWS() : The measurement object has not been retrieved ! Please check." << std::endl;
        return 0;
    }
    RooStats::HistFactory::HistoToWorkspaceFactoryFast factory(*measurement);
    
    // Creating the combined model
    RooWorkspace* ws = factory.MakeCombinedModel( vec_chName, vec_ws );
    
    cout << "...................................." << endl;
    
    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws, *measurement );
    
    return ws;
}

//__________________________________________________________________________________
//
void MultiFit::SaveCombinedWS(){
    //
    // Creating the rootfile
    //
    TFile *f = new TFile( (fName+"/ws_combined.root").c_str() , "recreate" );
    //
    // Creating the workspace
    //
    RooWorkspace *ws = CombineWS();
    //
    // Save the workspace
    //
    f->cd();
    ws->Write("combWS");
    f->Close();
}

//__________________________________________________________________________________
//
std::map < std::string, double > MultiFit::FitCombinedWS(int fitType, string inputData){
    TFile *f = new TFile((fName+"/ws_combined.root").c_str() );
    RooWorkspace *ws = (RooWorkspace*)f->Get("combWS");
    
    std::map < std::string, double > result;
    
    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////
    
    //
    // Fit configuration (1: SPLUSB or 2: BONLY)
    //
    FittingTool *fitTool = new FittingTool();
    if(fitType==2){
        fitTool -> ValPOI(0.);
        fitTool -> ConstPOI(true);
    } else if(fitType==1){
        fitTool -> ValPOI(1.);
        fitTool -> ConstPOI(false);
    }

//     if(fVarNameMinos.size()>0){
//       std::cout << "Setting the variables to use MINOS with" << std::endl;
//       fitTool -> UseMinos(fVarNameMinos);
//     }
    
    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)ws->obj("ModelConfig");
    RooSimultaneous *simPdf = (RooSimultaneous*)(mc->GetPdf());
    
    //
    // Creates the data object
    //
    RooDataSet* data = 0;
    if(inputData!=""){
        data = (RooDataSet*)ws->data( inputData.c_str() );
    } else {
        std::cout << "In MultiFit::FitCombinedWS() function: you didn't specify inputData => will try with observed data !" << std::endl;
        data = (RooDataSet*)ws->data("obsData");
        if(!data){
            std::cout << "In MultiFit::FitCombinedWS() function: observed data not present => will use with asimov data !" << std::endl;
            data = (RooDataSet*)ws->data("asimovData");
        }
    }
    
    // Performs the fit
    gSystem -> mkdir((fName+"/Fits/").c_str(),true);
    fitTool -> MinimType("Minuit2");
    fitTool -> FitPDF( mc, simPdf, data );
    fitTool -> ExportFitResultInTextFile(fName+"/Fits/"+fName+".txt");
    result = fitTool -> ExportFitResultInMap();
    
    return result;
}
//__________________________________________________________________________________
//
void MultiFit::GetCombinedLimit(string inputData){ // or asimovData
    string wsFileName = fName+"/ws_combined.root";
//     gSystem->Exec( ("mkdir ") );
    string cmd;
    cmd = "root -l -b -q 'runAsymptoticsCLs.C+(\""+wsFileName+"\",\"combWS\",\"ModelConfig\",\""+inputData+"\",\"asimovData_0\",\"./"+fName+"/Limits/\",\""+fName+"\",0.95)'";
    
    //
    // Finally computing the limit
    //
    gSystem->Exec(cmd.c_str());
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
    if(fCombine){
        std::cout << "Adding Combined Fit" << std::endl;
        names.push_back( fName );
        titles.push_back( "Combined" );
        suffs.push_back( "" );
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
    
    bool isComb = false;
    
    // get values
    TtHFit *fit = 0x0;
    for(int i=0;i<N;i++){
        if(fCombine && i==N-1) isComb = true;
        //
        if(!isComb) fit = fFitList[i];
//         fFitList[i]->ReadFitResults(names[i]+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0"+suffs[i]+".txt");
        if(!isComb)       fit->ReadFitResults(names[i]+"/Fits/"+names[i]+suffs[i]+".txt");
        else              fit->ReadFitResults(fName+"/Fits/"+fName+".txt");
        found = false;
        for(unsigned int j = 0; j<fit->fFitResults->fNuisPar.size(); ++j){
            par = fit->fFitResults->fNuisPar[j];
            if(par->fName == POI){
                g_central->SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_stat->SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_tot->SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_stat->SetPointError(N-i-1,0,0);
                g_tot->SetPointError(N-i-1,par->fPostFitUp,0);
                if(par->fFitValue+par->fPostFitUp > xmax) xmax = par->fFitValue+par->fPostFitUp;
                found = true;
                break;
            }
        }
        if(!found){
            g_central->SetPoint(N-i-1,-10,N-i-1);
            g_stat->SetPoint(N-i-1,-10,N-i-1);
            g_tot->SetPoint(N-i-1,-10,N-i-1);
            g_stat->SetPointError(N-i-1,0,0);
            g_tot->SetPointError(N-i-1,0,0);
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
    
    TH1F* h_dummy = new TH1F("h_dummy","h_dummy",1,fPOIMin,fPOIMax + (fPOIMax-fPOIMin));
    h_dummy->Draw();
    h_dummy->SetMinimum(ymin);
    h_dummy->SetMaximum(ymax);
    h_dummy->SetLineColor(kWhite);
    h_dummy->GetYaxis()->Set(N+1,ymin,ymax);
    h_dummy->GetYaxis()->SetNdivisions(Ndiv);
    
    TLatex *tex = new TLatex();
//     tex->SetNDC();

    for(int i=0;i<N;i++){
        h_dummy->GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
//         myText(0.5,(1.*i)/(1.*N),kBlack,Form("#mu= %.1f",g_central->GetY()[i]));
//                 tex->DrawLatex(0.5,(1.*i)/(1.*N),Form("#mu= %.1f",g_central->GetY()[i]));
                tex->DrawLatex(0.5*xmax,N-i-1,Form("#mu= %.1f",g_central->GetX()[N-i-1]));
                tex->DrawLatex(0.7*xmax,N-i-1,Form("^{+%.1f}",g_tot->GetErrorXhigh(N-i-1)));
                tex->DrawLatex(0.7*xmax,N-i-1,Form("_{-%.1f}",g_tot->GetErrorXlow(N-i-1)));
                tex->DrawLatex(0.85*xmax,N-i-1,Form("^{+%.1f}",g_stat->GetErrorXhigh(N-i-1)));
                tex->DrawLatex(0.85*xmax,N-i-1,Form("_{-%.1f}",g_stat->GetErrorXlow(N-i-1)));
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
    if(fCombine){
        std::cout << "Adding combined limit" << std::endl;
        names.push_back( fName );
        titles.push_back( "Combined" );
        suffs.push_back( "" );
        if(fShowObserved) fFitShowObserved.push_back(true);
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
        
//         if(fFitShowObserved[i]) g_obs->SetPoint(i,h->GetBinContent(1),i);
//         else g_obs->SetPoint(i,-1,i);
//         g_exp->SetPoint(i,h->GetBinContent(2),i);
//         g_1s->SetPoint(i,h->GetBinContent(2),i);
//         g_2s->SetPoint(i,h->GetBinContent(2),i);
//         g_obs->SetPointError(i,0,0.5);
//         g_exp->SetPointError(i,0,0.5);
//         g_1s->SetPointError(i,h->GetBinContent(2)-h->GetBinContent(5),h->GetBinContent(4)-h->GetBinContent(2),0.5,0.5);
//         g_2s->SetPointError(i,h->GetBinContent(2)-h->GetBinContent(6),h->GetBinContent(3)-h->GetBinContent(2),0.5,0.5);
        if(fFitShowObserved[i]) g_obs->SetPoint(N-i-1,h->GetBinContent(1),N-i-1);
        else g_obs->SetPoint(N-i-1,-1,N-i-1);
        g_exp->SetPoint(N-i-1,h->GetBinContent(2),N-i-1);
        g_1s->SetPoint(N-i-1,h->GetBinContent(2),N-i-1);
        g_2s->SetPoint(N-i-1,h->GetBinContent(2),N-i-1);
        g_obs->SetPointError(N-i-1,0,0.5);
        g_exp->SetPointError(N-i-1,0,0.5);
        g_1s->SetPointError(N-i-1,h->GetBinContent(2)-h->GetBinContent(5),h->GetBinContent(4)-h->GetBinContent(2),0.5,0.5);
        g_2s->SetPointError(N-i-1,h->GetBinContent(2)-h->GetBinContent(6),h->GetBinContent(3)-h->GetBinContent(2),0.5,0.5);
        
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
//         h_dummy->GetYaxis()->SetBinLabel(i+1,titles[i].c_str());
        h_dummy->GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
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
    if(fCombine) N++;
    
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            std::cout << "Adding Combined Fit" << std::endl;
            names.push_back( fName );
            titles.push_back( "Combined" );
            suffs.push_back( "" );
        }
        else{
            names.push_back( fFitList[i_fit]->fName );
            titles.push_back( fFitLabels[i_fit] );
            suffs.push_back( fFitSuffs[i_fit] );
        }
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
        if(fCombine && i_fit==N-1) break;
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
        if(fCombine && i_fit==N-1) break;
//         fFitList[i_fit]->ReadFitResults(names[i_fit]+"/FitResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0"+suffs[i_fit]+".txt");      
        fFitList[i_fit]->ReadFitResults(names[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".txt");      
    }
    
    // exclude unused systematics
    std::vector<string> NamesNew; NamesNew.clear();
    std::vector<string> TitlesNew; TitlesNew.clear();
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        FitResults *fitRes;
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            fitRes = fFitList[i_fit]->fFitResults;
            for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
                par = fitRes->fNuisPar[j];
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
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            fitRes->ReadFromTXT(fName+"/Fits/"+fName+".txt");
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults;
        }
        for(unsigned int j = 0; j<fitRes->fNuisPar.size(); ++j){
            par = fitRes->fNuisPar[j];
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
//                 g[i_fit]->SetPointEXhigh(i_syst,  errUpMap[systName]  <1 ?  errUpMap[systName]   : 1);
//                 g[i_fit]->SetPointEXlow( i_syst, -errDownMap[systName]<1 ? -errDownMap[systName] : 1);
                g[i_fit]->SetPointEXhigh(i_syst,  errUpMap[systName]);
                g[i_fit]->SetPointEXlow( i_syst, -errDownMap[systName]);
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
//     TCanvas *c = new TCanvas("c","c",600,newHeight);
    TCanvas *c = new TCanvas("c","c",800,newHeight);
    c->SetTicks(1,0);
//     gPad->SetLeftMargin(0.05);
    gPad->SetLeftMargin(0.05/(8./6.));
//     gPad->SetRightMargin(0.33);
    gPad->SetRightMargin(0.5);
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
//     leg = new TLegend(0.01,1.-0.03*(30./max),0.99,0.99);
    leg = new TLegend(0.01,1.-0.03*(30./max),0.75,0.99);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
//     leg->SetNColumns(4);
    leg->SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        leg->AddEntry(g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg->Draw();
    
    gPad->RedrawAxis();
    
    c->SaveAs((fName+"/NuisPar_comp"+fSaveSuf+".png").c_str());
    delete c;
}
