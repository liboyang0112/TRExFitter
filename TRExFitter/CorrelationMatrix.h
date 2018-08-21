#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H

/// Framework includes
#include "TRExFitter/Common.h"

/// c++ includes
#include <string>
#include <vector>
#include <map>

class CorrelationMatrix {

public:
    CorrelationMatrix();
    ~CorrelationMatrix();

    //
    // Functions
    //
    void AddNuisPar(const std::string& p);
    void SetCorrelation(const std::string& p0, const std::string& p1,float corr);
    float GetCorrelation(const std::string& p0, const std::string& p1);
    void Draw(const std::string& path, const double corrMin = -1.);

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
