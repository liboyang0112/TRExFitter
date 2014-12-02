#include "TtHFitter/Sample.h"

// -------------------------------------------------------------------------------------------------
// Sample

Sample::Sample(string name){
  fName = name;
  fTitle = name;
  fFillColor = kWhite;
  fLineColor = kBlack;
  fIsData = false;
  fIsSignal = false;
  fNSyst = 0;
  fNNorm = 0;
//   fNOverallSysts = 0;
//   fNWeightSysts = 0;
//   fNNtupleSysts = 0;
  fNtuplePaths.clear();
  fNtupleNames.clear();
  fTreeName = "";
}
Sample::~Sample(){}

void Sample::SetIsData(bool isData){
  fIsData = isData;
}

void Sample::SetIsSignal(bool isSignal){
  fIsSignal = isSignal;
}

void Sample::SetTitle(string title){
  fTitle = title;
}

void Sample::SetFillColor(int color){
  fFillColor = color;
}

void Sample::SetLineColor(int color){
  fLineColor = color;
}

void Sample::AddNormFactor(string name,float nominal,float min,float max){
  fNormFactors[fNNorm] = new NormFactor(name,nominal,min,max);
  fNNorm ++;
}

void Sample::AddNormFactor(NormFactor* normFactor){
  fNormFactors[fNNorm] = normFactor;
  fNNorm ++;
}

void Sample::AddNtupleName(string name){
  fNtupleNames.push_back(name);
}

void Sample::AddNtuplePath(string path){
  fNtuplePaths.push_back(path);
}

void Sample::AddSystematic(Systematic* syst){
  fSystematics[fNSyst] = syst;
  fNSyst++;
}

void Sample::SetMCweight(string weight){
  fMCweight = weight;
}
