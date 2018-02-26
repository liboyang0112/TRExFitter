#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H

#include "TtHFitter/Common.h"

class CorrelationMatrix {

public:
    CorrelationMatrix();
    ~CorrelationMatrix();

    //
    // Functions
    //
    void AddNuisPar(std::string p);
    void SetCorrelation(std::string p0,std::string p1,float corr);
    float GetCorrelation(std::string p0,std::string p1);
    void Draw(std::string path, const double corrMin = -1.);

    //
    // Data members
    //
    std::vector<std::string> fNuisParNames;
    std::map<std::string,int> fNuisParIdx;
    std::map<std::string,bool> fNuisParIsThere;
    std::vector<std::string> fNuisParToHide;
    float fMatrix[MAXsyst][MAXsyst];
};

#endif
