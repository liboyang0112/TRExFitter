#include "TtHFitter/Common.h"

#ifndef __TthPlot__
#define __TthPlot__

const int MAXbins = 1000;
const int MAXSAMPLES = 100;

using namespace std;

class TthPlot {
  public:
    TthPlot(string name="c",int canvasWidth=600,int canvasHeight=700);
    ~TthPlot(){};
    
    void SetChannel(string name);
    void AddLabel(string name);
    void SetLumi(string name);
    void SetLumiScale(float scale);
    void SetCME(string name);
    void SetXaxis(string name,bool isNjet=false);
    void SetYaxis(string name);
    void SetYmaxScale(float scale);
    void SetBinLabel(int bin,string name);
    void SetBinWidth(float width);
    
    void SetData(TH1* h,string name="Data");
    void AddSignal(TH1* h,string name="Signal");
    void AddNormSignal(TH1* h,string name="Signal");
    void AddOverSignal(TH1* h,string name="Signal");
    void AddBackground(TH1* h,string name="MC");
    void SetTotBkgAsym(TGraphAsymmErrors* g);
    void SetTotBkg(TH1* h);

    void SetChi2KS(float chi2prob,float ksprob=-1,float chi2val=-1,int ndf=-1);
    void BlindData();
    
    void Draw(string options="");
    void SaveAs(string name);
    void WriteToFile(string name);
    
    TCanvas* GetCanvas();
    
    void SetBinBlinding(bool on,float threshold=0.02);
    void SetBinBlinding(bool on,TH1F*h_blind);
  
//   private:
    string fName;
    TH1* h_data;
    TGraphAsymmErrors* g_data;
    TH1* h_mc;
    TH1* h_bkg[MAXSAMPLES];
    TH1* h_signal[MAXSAMPLES];
    TH1* h_normsig[MAXSAMPLES];
    TH1* h_oversig[MAXSAMPLES];
    THStack* h_stack;
    TH1* h_tot;
    TGraphAsymmErrors* g_tot;
    TH1F* h_blinding;
    TH1* h_tot_bkg_prefit;

    TCanvas* c;
    TLegend* leg;
    TLegend* leg1;
    TPad* pad0;
    TPad* pad1;
    
    string xtitle;
    string ytitle;
    string data_name;
    std::vector< string > fBkgNames;
    std::vector< string > fSigNames;
    std::vector< string > fNormSigNames;
    std::vector< string > fOverSigNames;
    std::vector< string > fLabels;
    string fLumi;
    string fCME;
    string fATLASlabel;
    float yMaxScale;
    int NDF;
    float Chi2val;
    float Chi2prob;
    float KSprob;
    
    float fYmax;
    float fYmin;
    float fRatioYmax;
    float fRatioYmin;
    float fBinWidth;
    bool fIsNjet;
    bool fShowYields;
    string fBinLabel[MAXbins];
    float fLumiScale;
    float fBlindingThreshold;
    int fLegendNColumns;

public:
  const TH1* GetTotal() const { return h_tot; };
  TH1* GetTotBkg() const;
  
};

// function to get asymmetric error bars for hists
double GC_up(double data);
double GC_down(double data);
TGraphAsymmErrors* poissonize(TH1 *h);
TGraphAsymmErrors* histToGraph(TH1* h);
void SetHistBinWidth(TH1* h,float width);
void SetGraphBinWidth(TGraphAsymmErrors* g,float width);

#endif
