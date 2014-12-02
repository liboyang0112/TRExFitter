#ifndef __Systematic__
#define __Systematic__

class Systematic {
public:
  Systematic(string name);
  ~Systematic();

  bool IsOverallOnly();
  
  string fName;
  bool fIsOverall;
  bool fIsShape;
  float fOverallUp;
  float fOverallDown;
  
  string fWeightUp;
  string fWeightSufUp;  
  vector<string> fNtuplePathsUp;
  string fNtuplePathSufUp;
  vector<string> fNtupleFilesUp;
  string fNtupleFileSufUp;
  string fTreeNameUp;
  string fTreeNameSufUp;

  string fWeightDown;
  string fWeightSufDown;  
  vector<string> fNtuplePathsDown;
  string fNtuplePathSufDown;
  vector<string> fNtupleFilesDown;
  string fNtupleFileSufDown;
  string fTreeNameDown;
  string fTreeNameSufDown;
};

#endif
