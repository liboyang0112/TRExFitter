#include "TRExFitter/RankingManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/YamlConverter.h"

#include "AtlasUtils/AtlasLabels.h"
#include "AtlasUtils/AtlasUtils.h"

#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"

#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooSimultaneous.h"

#include <algorithm>
#include <fstream>

//__________________________________________________________________________________
//
RankingManager::RankingManager() :
    fOutputPath(""),
    fInjectGlobalObservables(false),
    fFitStrategy(1),
    fCPU(1),
    fRndRange(0.1),
    fUseRnd(false),
    fRndSeed(-999),
    fStatOnly(false),
    fAtlasLabel("Internal"),
    fLumiLabel("139 fb^{-1}"),
    fCmeLabel("13 TeV"),
    fHEPDataFormat(false),
    fName("MyFit"),
    fSuffix(""),
    fRankingMaxNP(9999),
    fRankingPOIName(fName),
    fUsePOISinRanking(false)
{
}

//__________________________________________________________________________________
//
RankingManager::~RankingManager() {
}

//__________________________________________________________________________________
//
void RankingManager::AddNuisPar(const std::string& name, const bool isNF) {

    // check if the NP already exists by checking the names
    auto it = std::find_if(fNuisPars.begin(), fNuisPars.end(),
                [&name](const std::pair<std::string, bool>& element){return element.first == name;});

    if (it != fNuisPars.end()) {
        WriteWarningStatus("RankingManager::AddNuisPar", "NP " + name + " already exists in the list, not adding it");
        return;
    }
    if(name.find("_bin_") != std::string::npos) {
        fNuisPars.emplace_back("gamma_" + name, isNF);
    } else {
        fNuisPars.emplace_back(name, isNF);
    }
}

//__________________________________________________________________________________
//
void RankingManager::RunRanking(FitResults* fitResults,
                                RooWorkspace* ws,
                                RooDataSet* data,
                                const std::vector<std::shared_ptr<NormFactor> >& nfs) const {

    if (fOutputPath == "") {
        WriteErrorStatus("RankingManager::RunRanking", "OutputPath not set, plese set it via SetOutputPath()");
        exit(EXIT_FAILURE);
    }

    if (!ws) {
        WriteErrorStatus("RankingManager::RunRanking", "Workspace is nullptr");
        exit(EXIT_FAILURE);
    }

    std::vector<std::ofstream> outFiles;
    for (const auto& poi : fPOINames) { 
        outFiles.emplace_back((fOutputPath+"_"+poi+".txt").c_str());

        if (!outFiles.back().good() || !outFiles.back().is_open()) {
            WriteErrorStatus("RankingManager::RunRanking", "Cannot open file at " + fOutputPath +"_" + poi + ".txt");
            exit(EXIT_FAILURE);
        }
    }

    RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());
    if (!mc || !simPdf || !data){
        WriteErrorStatus("RankingManager::RunRanking","At least one of the objects that is needed to run ranking is not present");
        exit(EXIT_FAILURE);
    }
    
    if (fInjectGlobalObservables && !fFitValues.empty()) {
        FitUtils::InjectGlobalObservables(ws, fFitValues);
    }

    ws->saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));
   
    double poiInitial(0.);
    for(const auto& inf : nfs) {
        if (inf->fName == fPOINames.at(0)) {
            poiInitial = inf->fNominal;
            break;
        }
    }
 
    FittingTool fitTool{};
    fitTool.SetUseHesse(false);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetDebug(TRExFitter::DEBUGLEVEL);
    fitTool.ValPOI(poiInitial);
    fitTool.SetNCPU(fCPU);
    fitTool.ConstPOI(false);
    if(fStatOnly){
        fitTool.NoGammas();
        fitTool.NoSystematics();
    }
    fitTool.SetRandomNP(fRndRange, fUseRnd, fRndSeed);

    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : nfs) {
        if (!fUsePOISinRanking && Common::FindInStringVector(fPOINames, inf->fName) >= 0) continue;
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->fNominal);
    }
    fitTool.SetNPs(npNames,npValues);
    
    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, nfs);
    
    std::vector<double> muhats;
    for (const auto& poi : fPOINames) {
        muhats.emplace_back(fitResults -> GetNuisParValue(poi));
    }
    for(const auto& iNP : fNuisPars){

        std::string npName = iNP.first;
        if (npName.find("_bin_") != std::string::npos) {
            npName = Common::ReplaceString(npName, "gamma_", "");
        }

        //
        // Getting the postfit values of the nuisance parameter
        const double central = fitResults -> GetNuisParValue(  npName);
        const double up      = fitResults -> GetNuisParErrUp(  npName);
        const double down    = fitResults -> GetNuisParErrDown(npName);
        
        for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
            if (iNP.first == fPOINames.at(ipoi)) continue;
            
            outFiles.at(ipoi) << iNP.first << "   " << central << " +" << fabs(up) << " -" << fabs(down)<< "  ";
        }

        RankingManager::RankingValues values;
        values.central = central;
        values.up = up;
        values.down = down;

        //
        std::vector<double> dMuUp   = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, true, false, values, muhats);
        std::vector<double> dMuDown = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, false, false, values, muhats);

        for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
            if (iNP.first == fPOINames.at(ipoi)) continue;

            outFiles.at(ipoi) << dMuUp.at(ipoi) << "   " << dMuDown.at(ipoi) << "  ";
        }

        dMuUp   = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, true, true, values, muhats);
        dMuDown = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, false, true, values, muhats);
            
        for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
            if (iNP.first == fPOINames.at(ipoi)) continue;
            
            outFiles.at(ipoi) << dMuUp.at(ipoi) << "   " << dMuDown.at(ipoi) << " " << std::endl;
        }

    }
 
    ws->loadSnapshot("tmp_snapshot");
    for (auto& ifile : outFiles) {
        ifile.close();
    }
}

//__________________________________________________________________________________
//
std::vector<double> RankingManager::RunSingleFit(FittingTool* fitTool,
                                                 RooWorkspace* ws,       
                                                 RooStats::ModelConfig *mc,
                                                 RooSimultaneous *simPdf,
                                                 RooDataSet* data,
                                                 const std::pair<std::string, bool>& np,
                                                 const bool isUp,
                                                 const bool isPrefit,
                                                 const RankingManager::RankingValues& values,
                                                 const std::vector<double>& muhats) const {

    if (isPrefit && np.second) {
        std::vector<double> tmp(muhats.size(), 0);
        return tmp;
    }    

    ws->loadSnapshot("tmp_snapshot");
    fitTool->ResetFixedNP();
    if(fFitFixedNPs.size()>0){
        for(const auto& nuisParToFix : fFitFixedNPs){
            fitTool->FixNP(nuisParToFix.first,nuisParToFix.second);
        }
    }

    double shift = 1.0;
    if (isPrefit) {
       shift = isUp ? 1. : -1.;
    } else {
       shift = isUp ? std::fabs(values.up) : -std::fabs(values.down);
    }

    // Experimental: reduce the range of ranking
    if(TRExFitter::OPTION["ReduceRanking"]!=0){
        shift *= TRExFitter::OPTION["ReduceRanking"];
    }

    fitTool->FixNP( np.first, values.central + shift);
    fitTool->FitPDF( mc, simPdf, data );

    std::vector<double> result(fPOINames.size());
    for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
        result.at(ipoi) = fitTool->ExportFitResultInMap()[fPOINames.at(ipoi)] - muhats.at(ipoi);
    }

    if(TRExFitter::OPTION["ReduceRanking"]!=0){
        for (auto& i : result) {
            i /= TRExFitter::OPTION["ReduceRanking"];
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void RankingManager::PlotRanking(const std::vector<Region*>& regions,
                                 const bool flagSysts,
                                 const bool flagGammas) const {
    
    unsigned int maxNP = fRankingMaxNP;
    //
    std::string paramname;
    double nuiphat;
    double nuiperrhi;
    double nuiperrlo;
    double PoiUp;
    double PoiDown;
    double PoiNomUp;
    double PoiNomDown;
    std::vector<std::string> parname;
    std::vector<double> nuhat;
    std::vector<double> nuerrhi;
    std::vector<double> nuerrlo;
    std::vector<double> poiup;
    std::vector<double> poidown;
    std::vector<double> poinomup;
    std::vector<double> poinomdown;
    std::vector<double> number;

    std::vector<YamlConverter::RankingContainer> containerVec;

    std::ifstream fin(fOutputPath.c_str());
    fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    std::string temp_string = "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"";
    if (paramname=="Luminosity"){
        WriteErrorStatus("RankingManager::PlotRanking", temp_string);
        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
    }
    while (!fin.eof()){
        if(paramname.find("gamma")!=std::string::npos && !flagGammas){
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            if (paramname=="Luminosity"){
                WriteErrorStatus("RankingManager::PlotRanking", temp_string);
                fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            }
            continue;
        }
        if(paramname.find("gamma")==std::string::npos && !flagSysts){
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            if (paramname=="Luminosity"){
                WriteErrorStatus("RankingManager::PlotRanking", temp_string);
                fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
            }
            continue;
        }
        parname.push_back(paramname);
        nuhat.push_back(nuiphat);
        nuerrhi.push_back(nuiperrhi);
        nuerrlo.push_back(nuiperrlo);
        poiup.push_back(PoiUp);
        poidown.push_back(PoiDown);
        poinomup.push_back(PoiNomUp);
        poinomdown.push_back(PoiNomDown);

        YamlConverter::RankingContainer container;
        container.name = paramname;        
        container.nphat = nuiphat;        
        container.nperrhi = nuiperrhi;        
        container.nperrlo = nuiperrlo; 
        container.poihi = PoiUp;       
        container.poilo = PoiDown;       
        container.poiprehi = PoiNomUp;       
        container.poiprelo = PoiNomDown;       

        containerVec.emplace_back(std::move(container));

        fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        if (paramname=="Luminosity"){
            WriteErrorStatus("RankingManager::PlotRanking", temp_string);
            fin >> paramname >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown;
        }
    }

    {
        YamlConverter converter{};
        converter.WriteRanking(containerVec, fName+"/Ranking"+fSuffix+".yaml");
        if (fHEPDataFormat) {
            converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
            converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
            converter.WriteRankingHEPData(containerVec, fName, fSuffix);
        }
    }

    unsigned int SIZE = parname.size();
    WriteDebugStatus("RankingManager::PlotRanking", "NP ordering...");
    number.push_back(0.5);
    for (unsigned int i=1;i<SIZE;i++){
        number.push_back(i+0.5);
        double sumi = 0.0;
        int index=-1;
        sumi += std::max( std::abs(poiup[i]),std::abs(poidown[i]) );
        for (unsigned int j=1;j<=i;j++){
            double sumii = 0.0;
            sumii += std::max(std::abs(poiup[i-j]),std::abs(poidown[i-j]) );
            if (sumi<sumii){
                if (index==-1){
                    std::swap(poiup[i],poiup[i-j]);
                    std::swap(poidown[i],poidown[i-j]);
                    std::swap(poinomup[i],poinomup[i-j]);
                    std::swap(poinomdown[i],poinomdown[i-j]);
                    std::swap(nuhat[i],nuhat[i-j]);
                    std::swap(nuerrhi[i],nuerrhi[i-j]);
                    std::swap(nuerrlo[i],nuerrlo[i-j]);
                    std::swap(parname[i],parname[i-j]);
                    index=i-j;
                }
                else{
                    std::swap(poiup[index],poiup[i-j]);
                    std::swap(poidown[index],poidown[i-j]);
                    std::swap(poinomup[index],poinomup[i-j]);
                    std::swap(poinomdown[index],poinomdown[i-j]);
                    std::swap(nuhat[index],nuhat[i-j]);
                    std::swap(nuerrhi[index],nuerrhi[i-j]);
                    std::swap(nuerrlo[index],nuerrlo[i-j]);
                    std::swap(parname[index],parname[i-j]);
                    index=i-j;
                }
            }
            else{
                break;
            }
        }
    }
    number.push_back(parname.size()-0.5);

    double poimax = 0;
    for (unsigned int i=0;i<SIZE;i++) {
        poimax = std::max(poimax,std::max(std::abs(poiup[i]),std::abs(poidown[i]) ));
        poimax = std::max(poimax,std::max( std::abs(poinomup[i]),std::abs(poinomdown[i]) ));
        nuerrlo[i] = std::abs(nuerrlo[i]);
    }
    poimax *= 1.2;

    for (unsigned int i=0;i<SIZE;i++) {
        poiup[i]     *= (2./poimax);
        poidown[i]   *= (2./poimax);
        poinomup[i]  *= (2./poimax);
        poinomdown[i]*= (2./poimax);
    }

    // Restrict to the first N
    if(SIZE>maxNP) SIZE = maxNP;

    // Graphical part - rewritten taking DrawPulls in TRExFitter
    double lineHeight  =  30.;
    double offsetUp    =  60.; // external
    double offsetDown  =  60.;
    double offsetUp1   = 100.; // internal
    double offsetDown1 =  15.;
    int offset = offsetUp + offsetDown + offsetUp1 + offsetDown1;
    int newHeight = offset + SIZE*lineHeight;

    double xmin = -2.;
    double xmax =  2.;
    double max  =  0.;

    TGraphAsymmErrors g{};
    TGraphAsymmErrors g1{};
    TGraphAsymmErrors g2{};
    TGraphAsymmErrors g1a{};
    TGraphAsymmErrors g2a{};

    int idx = 0;
    std::vector< std::string > Names;
    std::string parTitle;

    for(unsigned int i = parname.size()-SIZE; i<parname.size(); ++i){
        g.SetPoint(idx, nuhat[i],  idx+0.5);
        g.SetPointEXhigh(      idx, nuerrhi[i]);
        g.SetPointEXlow(       idx, nuerrlo[i]);

        g1.SetPoint(      idx, 0.,idx+0.5);
        g1.SetPointEXhigh(idx, poiup[i]);
        g1.SetPointEXlow( idx, 0.);
        g1.SetPointEYhigh(idx, 0.4);
        g1.SetPointEYlow( idx, 0.4);

        g2.SetPoint(      idx, 0.,idx+0.5);
        g2.SetPointEXhigh(idx, poidown[i]);
        g2.SetPointEXlow( idx, 0.);
        g2.SetPointEYhigh(idx, 0.4);
        g2.SetPointEYlow( idx, 0.4);

        g1a.SetPoint(      idx, 0.,idx+0.5);
        g1a.SetPointEXhigh(idx, poinomup[i]);
        g1a.SetPointEXlow( idx, 0.);
        g1a.SetPointEYhigh(idx, 0.4);
        g1a.SetPointEYlow( idx, 0.4);

        g2a.SetPoint(      idx, 0.,idx+0.5);
        g2a.SetPointEXhigh(idx, poinomdown[i]);
        g2a.SetPointEXlow( idx, 0.);
        g2a.SetPointEYhigh(idx, 0.4);
        g2a.SetPointEYlow( idx, 0.4);

        if(parname[i].find("gamma")!=std::string::npos){
            // get name of the region
            std::vector<std::string> tmpVec = Vectorize(parname[i],'_');
            int nWords = tmpVec.size();
            std::string regName = tmpVec[2];
            for(int i_word=3;i_word<nWords-2;i_word++){
                regName += tmpVec[i_word];
            }
            // find the short label of this region
            std::string regTitle = regName;
            for(const auto& ireg: regions) {
                if(ireg->fName==regName){
                    regTitle = ireg->fShortLabel;
                    break;
                }
            }
            // build the title of the nuis par
            parTitle = "#gamma (" + regTitle + " bin " + tmpVec[nWords-1] + ")";
        }
        else parTitle = TRExFitter::SYSTMAP[ parname[i] ];

        Names.push_back(parTitle);

        idx ++;
        if(idx > max)  max = idx;
    }
    int newWidth = 600;
    if (fNPRankingCanvasSize.size() != 0){
        newWidth = fNPRankingCanvasSize.at(0);
        newHeight = fNPRankingCanvasSize.at(1);
    }
    TCanvas c("c","c",newWidth,newHeight);
    c.SetTicks(0,0);
    gPad->SetLeftMargin(0.4);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum( SIZE + offsetUp1/lineHeight   );
    h_dummy.SetMinimum(      - offsetDown1/lineHeight );
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);
    for(int i_bin=0;i_bin<h_dummy.GetNbinsX()+1;i_bin++){
        h_dummy.SetBinContent(i_bin,-10);
    }

    g1.SetFillColor(kAzure-4);
    g2.SetFillColor(kCyan);
    g1.SetLineColor(g1.GetFillColor());
    g2.SetLineColor(g2.GetFillColor());

    g1a.SetFillColor(kWhite);
    g2a.SetFillColor(kWhite);
    g1a.SetLineColor(kAzure-4);
    g2a.SetLineColor(kCyan);
    g1a.SetFillStyle(0);
    g2a.SetFillStyle(0);
    g1a.SetLineWidth(1);
    g2a.SetLineWidth(1);

    g.SetLineWidth(2);

    g1a.Draw("5 same");
    g2a.Draw("5 same");
    g1.Draw("2 same");
    g2.Draw("2 same");
    g.Draw("p same");

    TLatex systs{};
    systs.SetTextAlign(32);
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs.DrawLatex(xmin-0.1,i+0.5,Names[i].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    h_dummy.GetXaxis()->CenterTitle();
    h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy.GetXaxis()->SetTitleOffset(1.2);

    TGaxis axis_up( -2, SIZE + (offsetUp1)/lineHeight, 2, SIZE + (offsetUp1)/lineHeight, -poimax,poimax, 510, "-" );
    axis_up.SetLabelOffset( 0.01 );
    axis_up.SetLabelSize(   h_dummy.GetXaxis()->GetLabelSize() );
    axis_up.SetLabelFont(   gStyle->GetTextFont() );
    axis_up.Draw();
    axis_up.CenterTitle();
    axis_up.SetTitle(("#Delta"+fRankingPOIName).c_str());
    if(SIZE==20) axis_up.SetTitleOffset(1.5);
    axis_up.SetTitleSize(   h_dummy.GetXaxis()->GetLabelSize() );
    axis_up.SetTitleFont(   gStyle->GetTextFont() );

    TPad pad1("p1","Pad High",0,(newHeight-offsetUp-offsetUp1)/newHeight,0.4,1);
    pad1.Draw();

    pad1.cd();
    TLegend leg1(0.02,0.7,1,1.0,("Pre-fit impact on "+fRankingPOIName+":").c_str());
    leg1.SetFillStyle(0);
    leg1.SetBorderSize(0);
    leg1.SetMargin(0.25);
    leg1.SetNColumns(2);
    leg1.SetTextFont(gStyle->GetTextFont());
    leg1.SetTextSize(gStyle->GetTextSize());
    leg1.AddEntry(&g1a,"#theta = #hat{#theta}+#Delta#theta","f");
    leg1.AddEntry(&g2a,"#theta = #hat{#theta}-#Delta#theta","f");
    leg1.Draw();

    TLegend leg2(0.02,0.32,1,0.62,("Post-fit impact on "+fRankingPOIName+":").c_str());
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetMargin(0.25);
    leg2.SetNColumns(2);
    leg2.SetTextFont(gStyle->GetTextFont());
    leg2.SetTextSize(gStyle->GetTextSize());
    leg2.AddEntry(&g1,"#theta = #hat{#theta}+#Delta#hat{#theta}","f");
    leg2.AddEntry(&g2,"#theta = #hat{#theta}-#Delta#hat{#theta}","f");
    leg2.Draw();

    TLegend leg0(0.02,0.1,1,0.25);
    leg0.SetFillStyle(0);
    leg0.SetBorderSize(0);
    leg0.SetMargin(0.2);
    leg0.SetTextFont(gStyle->GetTextFont());
    leg0.SetTextSize(gStyle->GetTextSize());
    leg0.AddEntry(&g,"Nuis. Param. Pull","lp");
    leg0.Draw();

    c.cd();

    TLine l0(0,- offsetDown1/lineHeight,0,SIZE+0.5);// + offsetUp1/lineHeight);
    l0.SetLineStyle(kDashed);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    TLine l1 (-1,- offsetDown1/lineHeight,-1,SIZE+0.5);// + offsetUp1/lineHeight);
    l1.SetLineStyle(kDashed);
    l1.SetLineColor(kBlack);
    l1.Draw("same");
    TLine l2(1,- offsetDown1/lineHeight,1,SIZE+0.5);// + offsetUp1/lineHeight);
    l2.SetLineStyle(kDashed);
    l2.SetLineColor(kBlack);
    l2.Draw("same");

    if (fAtlasLabel!= "none") ATLASLabelNew(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.6*offsetUp1)/newHeight), fAtlasLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(       0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.3*offsetUp1)/newHeight), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    gPad->RedrawAxis();

    if(flagGammas && flagSysts){
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            c.SaveAs((fName+"/Ranking"+fSuffix+"."+format).c_str());
        }
    } else if(flagGammas){
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            c.SaveAs((fName+"/RankingGammas"+fSuffix+"."+format).c_str());
        }
    } else if(flagSysts){
        for(const auto& format: TRExFitter::IMAGEFORMAT) {
            c.SaveAs((fName+"/RankingSysts"+fSuffix+"."+format).c_str() );
        }
    } else{
        WriteWarningStatus("TRExFit::PlotNPRanking", "Your ranking plot felt in unknown category :s");
        for(const auto& format : TRExFitter::IMAGEFORMAT) {
            c.SaveAs((fName+"/RankingUnknown"+fSuffix+"."+format).c_str() );
        }
    }
}
