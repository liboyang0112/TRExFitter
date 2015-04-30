#include "TtHFitter/Common.h"

#ifndef __Systematic__
#define __Systematic__

class Systematic {
public:

    enum SystType{
        OVERALL, // 0
        SHAPE, // 1
        HISTO // 2
    };
  
    Systematic(string name,int type=0,float up=0,float down=0);
    ~Systematic();

    // -------
    // Members
    // -------
  
    string fName;
    string fTitle;
    int fType;
    int fSmoothType;
    int fSymmetrisationType;
      
    float fOverallUp;
    float fOverallDown;
    
    // from ntuples - up
    string fWeightUp;
    string fWeightSufUp;  
    std::vector<string> fNtuplePathsUp;
    string fNtuplePathSufUp;
    std::vector<string> fNtupleFilesUp;
    string fNtupleFileSufUp;
    std::vector<string> fNtupleNamesUp;
    string fNtupleNameSufUp;

    // from ntuples - down
    string fWeightDown;
    string fWeightSufDown;  
    std::vector<string> fNtuplePathsDown;
    string fNtuplePathSufDown;
    std::vector<string> fNtupleFilesDown;
    string fNtupleFileSufDown;
    std::vector<string> fNtupleNamesDown;
    string fNtupleNameSufDown;

    // from histos - up
    std::vector<string> fHistoPathsUp;
    string fHistoPathSufUp;
    std::vector<string> fHistoFilesUp;
    string fHistoFileSufUp;
    std::vector<string> fHistoNamesUp;
    string fHistoNameSufUp;
    
    // from histos - down
    std::vector<string> fHistoPathsDown;
    string fHistoPathSufDown;
    std::vector<string> fHistoFilesDown;
    string fHistoFileSufDown;
    std::vector<string> fHistoNamesDown;
    string fHistoNameSufDown;
};

#endif
