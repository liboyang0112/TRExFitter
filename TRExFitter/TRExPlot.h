#ifndef TRExPLOT_H
#define TRExPLOT_H

/// Framework includes
#include "TRExFitter/TRExFit.h"

/// c++ includes
#include <memory>
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

class TRExPlot {
  public:
    explicit TRExPlot(std::string name="c",int canvasWidth=600,int canvasHeight=700,bool hideRatioPad=false);

    ~TRExPlot();

    TRExPlot(const TRExPlot& t) = delete;
    TRExPlot(TRExPlot&& t) = delete;
    TRExPlot& operator=(const TRExPlot& t) = delete;
    TRExPlot& operator=(TRExPlot&& t) = delete;

    void SetChannel(const std::string& name);
    void AddLabel(const std::string& name);
    void SetLumi(const std::string& name);
    void SetLumiScale(double scale);
    void SetCME(const std::string& name);
    void SetXaxis(const std::string& name,bool isNjet=false);
    void SetYaxis(const std::string& name);
    void SetYmaxScale(double scale);
    void ResizeBinLabel(const int n);
    void SetBinLabel(int bin, const std::string& name);
    void SetBinWidth(double width);

    void SetData(TH1* h,std::string name="Data");
    void AddSignal(TH1* h,std::string name="Signal");
    void AddNormSignal(TH1* h,std::string name="Signal");
    void AddOverSignal(TH1* h,std::string name="Signal");
    void AddBackground(TH1* h,std::string name="MC");
    void SetTotBkgAsym(TGraphAsymmErrors* g);
    void SetTotBkg(TH1* h);

    void SetChi2KS(double chi2prob,double ksprob=-1.,double chi2val=-1.,int ndf=-1);
    void BlindData();

    void SetXaxisRange(const std::vector<double>& vec) {fXaxisRange = vec;}

    void Draw(std::string options="");
    void SaveAs(const std::string& name) const;
    void WriteToFile(const std::string& name) const;

    TCanvas* GetCanvas() const;


    std::string fName;
    TH1* h_data;
    TGraphAsymmErrors* g_data;
    std::vector<TH1*> h_bkg;
    std::vector<TH1*> h_signal;
    std::vector<TH1*> h_normsig;
    std::vector<TH1*> h_oversig;
    THStack* h_stack;
    TH1* h_tot;
    TGraphAsymmErrors* g_tot;
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
    double yMaxScale;
    int NDF;
    double Chi2val;
    double Chi2prob;
    double KSprob;

    double fYmax;
    double fYmin;
    double fRatioYmax;
    double fRatioYmin;
    double fBinWidth;
    bool fIsNjet;
    bool fShowYields;
    std::vector<std::string> fBinLabel;
    double fLumiScale;
    int fLegendNColumns;
    std::vector<double> fXaxisRange;
    std::string fRatioYtitle;

    std::string fRatioType;

    double fLabelX;
    double fLabelY;
    double fLegendX1;
    double fLegendX2;
    double fLegendY;

    std::vector<int> fBlindedBins;

public:
    TH1* GetTotBkg() const;

    void SetBinBlinding(const std::vector<int>& bins);
    
    const std::vector<int>& GetBlindedBins() const;

    const TH1D* GetBlindingHisto() const;
    
    void SetBlindingHisto(const TH1* h);

private:
    std::unique_ptr<TH1D> fBlinding;

};

// function to get asymmetric error bars for hists
double GC_up(double data);
double GC_down(double data);
TGraphAsymmErrors* poissonize(TH1 *h);
TGraphAsymmErrors* histToGraph(TH1* h);
void SetHistBinWidth(TH1* h,double width);
void SetGraphBinWidth(TGraphAsymmErrors* g,double width);

#endif
