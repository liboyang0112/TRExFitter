// Class include
#include "TRExFitter/CorrelationMatrix.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/YamlConverter.h"

// ROOT includes
#include "TCanvas.h"
#include "TH2D.h"
#include "TPad.h"
#include "TStyle.h"

using std::string;

//__________________________________________________________________________________
//
CorrelationMatrix::CorrelationMatrix() :
    fOutFolder("") {
}

//__________________________________________________________________________________
//
void CorrelationMatrix::AddNuisPar(const string& p){
    fNuisParIdx[p] = (int)fNuisParNames.size();
    fNuisParNames.push_back(p);
    fNuisParIsThere[p] = true;
}

//__________________________________________________________________________________
//
void CorrelationMatrix::Resize(const int size) {
    fMatrix.resize(size);
    for (auto& i : fMatrix) {
        i.resize(size);
    }
}

//__________________________________________________________________________________
//
void CorrelationMatrix::SetCorrelation(const string& p0, const string& p1,double corr){
    const std::size_t idx0 = fNuisParIdx[p0];
    const std::size_t idx1 = fNuisParIdx[p1];
    fMatrix[idx0][idx1] = corr;
}

//__________________________________________________________________________________
//
double CorrelationMatrix::GetCorrelation(const string& p0, const string& p1){
    bool isMorph_p0 = false;
    bool isMorph_p1 = false;
    if (p0.find("morph_") != std::string::npos) isMorph_p0 = true;
    if (p1.find("morph_") != std::string::npos) isMorph_p1 = true;
    // if one of the two is missing, return 1 or 0 (if name1==name2 ==> 1, not zero!)
    if(!fNuisParIsThere[p0]){
        if(!isMorph_p0) WriteVerboseStatus("CorrelationMatrix::GetCorrelation", "NP " + p0 + " not found in correlation matrix. Returning correlation = " + std::to_string(1.*(p0==p1)));
        else WriteVerboseStatus("CorrelationMatrix::GetCorrelation", "NP " + p0 + " not found in correlation matrix. The NP is for morphing. Returning correlation = " + std::to_string(1.*(p0==p1)));
        return 1.*(p0==p1);
    }
    if(!fNuisParIsThere[p1]){
        if(!isMorph_p1) WriteVerboseStatus("CorrelationMatrix::GetCorrelation", "NP " + p1 + " not found in correlation matrix. Returning correlation = " + std::to_string(1.*(p0==p1)));
        else WriteVerboseStatus("CorrelationMatrix::GetCorrelation", "NP " + p0 + " not found in correlation matrix. The NP is for morphing. Returning correlation = " + std::to_string(1.*(p0==p1)));
        return 1.*(p0==p1);
    }
    int idx0 = fNuisParIdx[p0];
    int idx1 = fNuisParIdx[p1];
    return fMatrix[idx0][idx1];
}


//__________________________________________________________________________________
//
void CorrelationMatrix::Draw(const std::vector<std::string>& path, const bool& useGammas, const bool useHEPDataFormat, const double minCorr){
    //
    // 0) Determines the number of lines/columns
    //
    std::vector < string > vec_NP = fNuisParNames;
    std::vector<std::vector<double> > correlations(fNuisParNames.size(), std::vector<double>(fNuisParNames.size()));
    for(unsigned int iNP = 0; iNP < fNuisParNames.size(); ++iNP){
        const string iSystName = fNuisParNames[iNP];
        for(unsigned int jNP = 0; jNP < fNuisParNames.size(); ++jNP){
            const string jSystName = fNuisParNames[jNP];
            const double corr = GetCorrelation(iSystName, jSystName);
            correlations.at(iNP).at(jNP) = corr;
        }
    }
    if(minCorr>-1){
        vec_NP.clear();
        for(unsigned int iNP = 0; iNP < fNuisParNames.size(); ++iNP){
            const string iSystName = fNuisParNames[iNP];
            for(unsigned int jNP = 0; jNP < fNuisParNames.size(); ++jNP){
                const string jSystName = fNuisParNames[jNP];
                if(jNP == iNP) continue;
                const double corr = GetCorrelation(iSystName, jSystName);
                if(std::fabs(corr)>=minCorr){
                    WriteVerboseStatus("CorrelationMatrix::Draw", iSystName + " " + std::to_string(minCorr) + "    " + std::to_string(corr) + " (" + jSystName + ")");
                    vec_NP.push_back(iSystName);
                    break;
                }
            }
        }
    }
    int N = vec_NP.size();

    // pass everythign to yaml converter
    YamlConverter converter{};
    converter.WriteCorrelation(fNuisParNames, correlations, fOutFolder);
    if (useHEPDataFormat) {
        converter.WriteCorrelationHEPData(fNuisParNames, correlations, fOutFolder);
    }

    //
    // 0.5) Skip some NPs
    //
    static const std::vector<string> npToExclude = {"gamma_","stat_"};
    std::vector < std::string > vec_NP_old = vec_NP;
    vec_NP.clear();
    for(unsigned int iNP = 0; iNP < vec_NP_old.size(); ++iNP){
        const std::string iSystName = vec_NP_old[iNP];
        bool skip = false;
        if (!useGammas){
            for(const std::string& ii : npToExclude){
                if(iSystName.find(ii)!=std::string::npos){
                    skip = true;
                    continue;
                }
            }
        }
        if(skip) continue;
        if( Common::FindInStringVector(fNuisParToHide,iSystName)>=0 ) continue;
        vec_NP.push_back(iSystName);
    }
    N = vec_NP.size();
    
    //
    // 0.75) Reorder NPs
    if(fNuisParList.size()>0){
        vec_NP_old = vec_NP;
        vec_NP.clear();
        for(auto& npName : fNuisParList){
            for(unsigned int iNP = 0; iNP < vec_NP_old.size(); ++iNP){
                if(vec_NP_old[iNP]==npName) vec_NP.emplace_back(vec_NP_old[iNP]);
            }
        }
        N = vec_NP.size();
    }

    //
    // 1) Performs the plot
    //
    TH2D h_corr("h_corr","",N,0,N,N,0,N);
    h_corr.SetDirectory(0);

    for(unsigned int iNP = 0; iNP < vec_NP.size(); ++iNP){//line number
        const string iSystName = vec_NP[iNP];

        if(TRExFitter::SYSTMAP[iSystName]!=""){
            h_corr.GetXaxis()->SetBinLabel(iNP+1,TRExFitter::SYSTMAP[iSystName].c_str());
            h_corr.GetYaxis()->SetBinLabel(N-iNP,TRExFitter::SYSTMAP[iSystName].c_str());
        }
        else{
            h_corr.GetXaxis()->SetBinLabel(iNP+1,iSystName.c_str());
            h_corr.GetYaxis()->SetBinLabel(N-iNP,iSystName.c_str());
        }

        for(unsigned int jNP = 0; jNP < vec_NP.size(); ++jNP){//column number
            const string jSystName = vec_NP[jNP];

            h_corr.SetBinContent(iNP+1,N-jNP,100.*GetCorrelation(iSystName, jSystName));

        }
    }
    h_corr.SetMinimum(-100.);
    h_corr.SetMaximum(100.);

    int size = 500;
    if(vec_NP.size()>10){
        size = vec_NP.size()*50;
    }

    //
    // 2) Style settings
    //
    TCanvas c1("","",0.,0.,size+300,size+300);
    gStyle->SetPalette(87);
    h_corr.SetMarkerSize(0.75*1000);
    gStyle->SetPaintTextFormat(".1f");
    gPad->SetLeftMargin(0.5*600/(size+300));
    gPad->SetBottomMargin(0.5*600/(size+300));
    gPad->SetRightMargin(0.);
    gPad->SetTopMargin(0.);

    h_corr.GetXaxis()->LabelsOption("v");
    h_corr.GetXaxis()->SetLabelSize( h_corr.GetXaxis()->GetLabelSize()*0.75 );
    h_corr.GetYaxis()->SetLabelSize( h_corr.GetYaxis()->GetLabelSize()*0.75 );
    c1.SetTickx(0);
    c1.SetTicky(0);
    h_corr.GetYaxis()->SetTickLength(0);
    h_corr.GetXaxis()->SetTickLength(0);
    c1.SetGrid();
    h_corr.Draw("col TEXT");
    c1.RedrawAxis("g");

    for (const auto& ipath : path) {
        c1.SaveAs(ipath.c_str());
    }
}
