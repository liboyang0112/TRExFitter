#include "TtHFitter/Common.h"

#ifndef __ShapeFactor__
#define __ShapeFactor__

class ShapeFactor{
public:
  ShapeFactor();
  ShapeFactor(string name, float nominal=1, float min=0, float max=10, bool isConst=false);
  ~ShapeFactor();
  void Set(string name, float nominal=1, float min=0, float max=10, bool isConst=false);

  void Print();
  
  string fName;
  string fNuisanceParameter;
  string fTitle;
  string fCategory;
  
  float fNominal;
  float fMin;
  float fMax;
  bool fConst;
  
  std::vector<string> fRegions;
  std::vector<string> fExclude;  
};
 
#endif
