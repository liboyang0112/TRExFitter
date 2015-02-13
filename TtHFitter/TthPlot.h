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

using namespace std;

class TthPlot {
  public:
    TthPlot();
    TthPlot(string name);
    ~TthPlot(){};

    void Init();
    
    void SetChannel(string name);
    void SetLumi(string name);
    void SetXaxis(string name,bool isNjet=false);
    void SetYaxis(string name);
    void SetYmaxScale(float scale);

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
};
 
#endif
