#ifndef TRExPLOT_H
#define TRExPLOT_H

/// c++ includes
#include <string>
#include <vector>

/// Forwards class declaration
class TCanvas;
class TGraphAsymmErrors;
class TH1;
class TH1D;
class THStack;
class TLegend;
class TPad;

const int MAXbins = 1000;
// const int MAXSAMPLES = 100;

class TRExPlot {
  public:
    TRExPlot(std::string name="c",int canvasWidth=600,int canvasHeight=700);
    ~TRExPlot();

    void SetChannel(const std::string& name);
    void AddLabel(const std::string& name);
    void SetLumi(const std::string& name);
    void SetLumiScale(float scale);
    void SetCME(const std::string& name);
    void SetXaxis(const std::string& name,bool isNjet=false);
    void SetYaxis(const std::string& name);
    void SetYmaxScale(float scale);
    void SetBinLabel(int bin, const std::string& name);
    void SetBinWidth(float width);

    void SetData(TH1* h,std::string name="Data");
    void AddSignal(TH1* h,std::string name="Signal");
    void AddNormSignal(TH1* h,std::string name="Signal");
    void AddOverSignal(TH1* h,std::string name="Signal");
    void AddBackground(TH1* h,std::string name="MC");
    void SetTotBkgAsym(TGraphAsymmErrors* g);
    void SetTotBkg(TH1* h);

    void SetChi2KS(float chi2prob,float ksprob=-1,float chi2val=-1,int ndf=-1);
    void BlindData();

    void Draw(std::string options="");
    void SaveAs(const std::string& name) const;
    void WriteToFile(const std::string& name) const;

    TCanvas* GetCanvas() const;

    void SetBinBlinding(bool on,float threshold=0.02);
    void SetBinBlinding(bool on,TH1D*h_blind);

    std::string fName;
    TH1* h_data;
    TGraphAsymmErrors* g_data;
    TH1* h_mc;
    std::vector<TH1*> h_bkg;
    std::vector<TH1*> h_signal;
    std::vector<TH1*> h_normsig;
    std::vector<TH1*> h_oversig;
    THStack* h_stack;
    TH1* h_tot;
    TGraphAsymmErrors* g_tot;
    TH1D* h_blinding;
    TH1* h_tot_bkg_prefit;
    TH1* h_dummy;

    TCanvas* c;
    TLegend* leg;
    TLegend* leg1;
    TPad* pad0;
    TPad* pad1;

    std::string xtitle;
    std::string ytitle;
    std::string fDataName;
    std::vector< std::string > fBkgNames;
    std::vector< std::string > fSigNames;
    std::vector< std::string > fNormSigNames;
    std::vector< std::string > fOverSigNames;
    std::vector< std::string > fLabels;
    std::string fLumi;
    std::string fCME;
    std::string fATLASlabel;
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
    std::string fBinLabel[MAXbins];
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
