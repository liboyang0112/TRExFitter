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
    void SetCorrelation(const std::string& p0, const std::string& p1,double corr);
    double GetCorrelation(const std::string& p0, const std::string& p1);

    /**
      * Function to draw correlation matrix
      * @param Path to the output file
      * @param Flag to include gammas on the matrix
      * @param Minimum correlation considered for plotting
      */ 
    void Draw(const std::string& path, const bool& useGammas, const double corrMin = -1.);

    //
    // Data members
    //
    std::vector<std::string> fNuisParNames;
    std::map<std::string,int> fNuisParIdx;
    std::map<std::string,bool> fNuisParIsThere;
    std::vector<std::string> fNuisParToHide;
    double fMatrix[MAXsyst][MAXsyst];
};

#endif
