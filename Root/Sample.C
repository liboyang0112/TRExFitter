#include "TtHFitter/Sample.h"

// -------------------------------------------------------------------------------------------------
// Sample

Sample::Sample(string name,int type){
  fName = name;
  fTitle = name;
  fType = type;
  fFillColor = kWhite;
  fLineColor = kBlack;
  fNSyst = 0;
  fNNorm = 0;
  fNormalizedByTheory = true;
  //
  // ntuples
  fSelection = "1";
  fMCweight = "1";
  fNtuplePaths.clear();
  fNtupleFiles.clear();
  fNtupleNames.clear();
  fNtuplePathSuffs.clear();
  fNtupleFileSuffs.clear();
  fNtupleNameSuffs.clear();
  //
  // histograms
  fHistoPaths.clear();
  fHistoFiles.clear();
  fHistoNames.clear();
  fHistoPathSuffs.clear();
  fHistoFileSuffs.clear();
  fHistoNameSuffs.clear();
    
  fNormFactors.clear();
  fSystematics.clear();
    
}
Sample::~Sample(){
    for(unsigned int i=0;i<fNormFactors.size();++i){
        if(fNormFactors[i]){
            delete fNormFactors[i];
        }
    }
    for(unsigned int i=0;i<fSystematics.size();++i){
        if(fSystematics[i]){
            delete fSystematics[i];
        }
    }
}

// cosmetics
void Sample::SetTitle(string title){
  fTitle = title;
}
void Sample::SetFillColor(int color){
  fFillColor = color;
}
void Sample::SetLineColor(int color){
  fLineColor = color;
}
void Sample::NormalizedByTheory(const bool norm ){
  fNormalizedByTheory = norm;
}


// read from ntuples
void Sample::SetMCweight(string weight){
  fMCweight = weight;
}
void Sample::SetSelection(string selection){
  fSelection = selection;
}
void Sample::AddNtuplePath(string path){
  fNtuplePaths.push_back(path);
}
void Sample::AddNtupleFile(string file){
  fNtupleFiles.push_back(file);
}
void Sample::AddNtupleName(string name){
  fNtupleNames.push_back(name);
}

// read from histograms
void Sample::AddHistoPath(string path){
  fHistoPaths.push_back(path);
}
void Sample::AddHistoFile(string file){
  fHistoFiles.push_back(file);
}
void Sample::AddHistoName(string name){
  fHistoNames.push_back(name);
}


// norm factors and systs
void Sample::AddNormFactor(NormFactor* normFactor){
  fNormFactors.push_back(normFactor);
  fNNorm ++;
}
void Sample::AddSystematic(Systematic* syst){
  fSystematics.push_back(syst);
  fNSyst++;
}
NormFactor* Sample::AddNormFactor(string name,float nominal,float min,float max){
  fNormFactors.push_back(new NormFactor(name,nominal,min,max));
  fNNorm ++;
  return fNormFactors[fNNorm-1];
}
Systematic* Sample::AddSystematic(string name,int type,float up,float down){
  fSystematics.push_back(new Systematic(name,type,up,down));
  fNSyst++;
  return fSystematics[fNSyst-1];
}
