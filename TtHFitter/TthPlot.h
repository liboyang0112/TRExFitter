#include "TFile.h"
#include "TObject.h"
#include "THStack.h"
#include "TH1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include <vector>
#include <string>
#include <iostream>

#ifndef __TthPlot__
#define __TthPlot__

const int MAXbins = 1000;

using namespace std;

class TthPlot {
  public:
//     TthPlot();
    TthPlot(string name="c",int canvasWidth=600,int canvasHeight=700);
    ~TthPlot(){};

//     void Init(int canvasWidth=600,int canvasHeight=700);
    
    void SetChannel(string name);
    void SetLumi(string name);
    void SetXaxis(string name,bool isNjet=false);
    void SetYaxis(string name);
    void SetYmaxScale(float scale);
    void SetBinLabel(int bin,string name);
    
    void SetData(TH1* h,string name="Data");
    void AddSignal(TH1* h,string name="Signal");
    void AddNormSignal(TH1* h,string name="Signal");
    void AddBackground(TH1* h,string name="MC");
    void SetTotBkgAsym(TGraphAsymmErrors* g);
    void SetTotBkg(TH1* h);

    void SetChi2KS(float chi2,float ks);
    
    void Draw(string options="");
    void SaveAs(string name);
    void WriteToFile(string name);
    
    TCanvas* GetCanvas();

//   private:
    string fName;
    TH1* h_data;
    TGraphAsymmErrors* g_data;
    TH1* h_mc;
    TH1* h_signal;
    TH1* h_normsig;
    TH1* h_bkg[100];
    THStack* h_stack;
    TH1* h_tot;
    TGraphAsymmErrors* g_tot;

    TCanvas* c;
    TLegend* leg;
    TLegend* leg1;
    TPad* pad0;
    TPad* pad1;
    
    string xtitle;
    string ytitle;
    string data_name;
    vector<string> sample_name;
    string leg_title;
    string fLumi;
    string fCME;
    string fATLASlabel;
    float yMaxScale;
    float Chi2prob;
    float KSprob;
    
    bool fIsNjet;
    
    string fBinLabel[MAXbins];
};

// function to get asymmetric error bars for hists
double GC_up(double data);
double GC_down(double data);
TGraphAsymmErrors* poissonize(TH1 *h);
TGraphAsymmErrors* histToGraph(TH1* h);

#endif
