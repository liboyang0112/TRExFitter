#include "TtHFitter/Common.h"

#ifndef __NormFactor__
#define __NormFactor__

class NormFactor{
public:
  NormFactor();
  NormFactor(string name, float nominal=1, float min=0, float max=10, bool isConst=false);
  ~NormFactor();
  void Set(string name, float nominal=1, float min=0, float max=10, bool isConst=false);

  void Print();
  
  string fName;
  string fNuisanceParameter;
  string fTitle;
  
  float fNominal;
  float fMin;
  float fMax;
  bool fConst;
  
  std::vector<string> fRegions;
  std::vector<string> fExclude;  
};
 
#endif
