#include "TtHFitter/Common.h"

#ifndef __CorrelationMatrix__
#define __CorrelationMatrix__

class CorrelationMatrix {
public:
  CorrelationMatrix();
  ~CorrelationMatrix();
  
  void AddNuisPar(string p);
  void SetCorrelation(string p0,string p1,float corr);
  float GetCorrelation(string p0,string p1);
  
  vector<string> fNuisParNames;
  map<string,int> fNuisParIdx;
  map<string,bool> fNuisParIsThere;
  float fMatrix[MAXsyst][MAXsyst];
};

#endif
