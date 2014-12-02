#ifndef __NuisParameter__
#define __NuisParameter__

class NuisParameter {
public:
  NuisParameter(string name);
  ~NuisParameter();

  string fName;
  float fStartValue;
  float fFitValue;
  float fPostFitUp; // this should be like +0.8... So alpha+deltaAlpha = fFitValue + fPostFitUp
  float fPostFitDown; // this like -0.7...
  int fConstrainType;
  int fInterpCode;
};

#endif
