#include "TtHFitter/CorrelationMatrix.h"
#include "TH2F.h"

//__________________________________________________________________________________
//
CorrelationMatrix::CorrelationMatrix(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
}

//__________________________________________________________________________________
//
CorrelationMatrix::~CorrelationMatrix(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
}

//__________________________________________________________________________________
//
void CorrelationMatrix::AddNuisPar(string p){
    fNuisParIdx[p] = (int)fNuisParNames.size();
    fNuisParNames.push_back(p);
    fNuisParIsThere[p] = true;
}

//__________________________________________________________________________________
//
void CorrelationMatrix::SetCorrelation(string p0,string p1,float corr){
    if(!fNuisParIsThere[p0]) AddNuisPar(p0);
    if(!fNuisParIsThere[p1]) AddNuisPar(p1);
    int idx0 = fNuisParIdx[p0];
    int idx1 = fNuisParIdx[p1];
    fMatrix[idx0][idx1] = corr;
}

//__________________________________________________________________________________
//
float CorrelationMatrix::GetCorrelation(string p0,string p1){
    if(!fNuisParIsThere[p0]){
        cout << "  WARNING: NP " << p0 << " not found in correlation matrix. Returning correlation = 0." << endl;
        return 0.;
    }
    if(!fNuisParIsThere[p1]){
        cout << "  WARNING: NP " << p1 << " not found in correlation matrix. Returning correlation = 0." << endl;
        return 0.;
    }
    int idx0 = fNuisParIdx[p0];
    int idx1 = fNuisParIdx[p1];
    return fMatrix[idx0][idx1];
}


//__________________________________________________________________________________
//
void CorrelationMatrix::Draw(const string &folder, const double minCorr){
    
    gStyle -> Reset();
    
    //
    // 0) Determines the number of lines/columns
    //
    std::vector < string > vec_NP = fNuisParNames;
    
    if(minCorr>-1){
        
        vec_NP.clear();
        for(unsigned int iNP = 0; iNP < fNuisParNames.size()-1; ++iNP){
            const string iSystName = fNuisParNames[iNP];
            
            for(unsigned int jNP = iNP+1; jNP < fNuisParNames.size(); ++jNP){
                const string jSystName = fNuisParNames[jNP];
                
                double corr = GetCorrelation(iSystName, jSystName);
                if(abs(corr)>=minCorr){
                    std::cout << minCorr << "    " << corr << std::endl;
                    vec_NP.push_back(iSystName);
                    break;
                }
                
            }
        }
    }
    const int N = vec_NP.size();
    
    //
    // 1) Performs the plot
    //
    TH2F *h_corr = new TH2F("h_corr","",N,0,N,N,0,N);
    h_corr->SetDirectory(0);
    
    for(unsigned int iNP = 0; iNP < vec_NP.size(); ++iNP){//line number
        const string iSystName = vec_NP[iNP];
        
        h_corr->GetXaxis()->SetBinLabel(iNP+1,(TString)iSystName);
        h_corr->GetYaxis()->SetBinLabel(N-iNP,(TString)iSystName);
        
        for(unsigned int jNP = 0; jNP < vec_NP.size(); ++jNP){//column number
            const string jSystName = vec_NP[jNP];
    
            h_corr -> SetBinContent(N-jNP,iNP+1,100.*GetCorrelation(iSystName, jSystName));
            
        }
    }
    h_corr->SetMinimum(-100.);
    h_corr->SetMaximum(100.);
    
    TCanvas *c1 = new TCanvas("","",0.,0.,600,600);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasBorderMode(0);
    c1->SetTickx();
    c1->SetTicky();
    c1->SetFillColor(0);
    
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);
    
    gStyle->SetPaintTextFormat(".1f");
    h_corr->Draw("colTEXT");
    h_corr->GetXaxis()->SetLabelSize(10);
    h_corr->GetYaxis()->SetLabelSize(10);
    h_corr->SetMarkerColor(1);
    h_corr->SetMarkerSize(1);
    c1 -> Print((TString)folder + "/CorrMarix.eps");
    c1 -> Print((TString)folder + "/CorrMarix.pdf");
    c1 -> Print((TString)folder + "/CorrMarix.png");
    
    delete c1;
    
}








