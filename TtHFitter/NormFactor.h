#include "TtHFitter/Common.h"

#ifndef __NormFactor__
#define __NormFactor__

class NormFactor{
public:
  NormFactor();
  NormFactor(string name, float nominal, float min, float max);
  ~NormFactor();
  void Set(string name, float nominal, float min, float max);

  void Print();
  
  string fName;
  float fNominal;
  float fMin;
  float fMax;
};
 
#endif
