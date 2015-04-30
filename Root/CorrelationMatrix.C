#include "TtHFitter/CorrelationMatrix.h"

//__________________________________________________________________________________
//
CorrelationMatrix::CorrelationMatrix(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
}

//__________________________________________________________________________________
//
CorrelationMatrix::~CorrelationMatrix(){
    fNuisParNames.clear();
    fNuisParIdx.clear();
    fNuisParIsThere.clear();
}

//__________________________________________________________________________________
//
void CorrelationMatrix::AddNuisPar(string p){
    fNuisParIdx[p] = (int)fNuisParNames.size();
    fNuisParNames.push_back(p);
    fNuisParIsThere[p] = true;
}

//__________________________________________________________________________________
//
void CorrelationMatrix::SetCorrelation(string p0,string p1,float corr){
    if(!fNuisParIsThere[p0]) AddNuisPar(p0);
    if(!fNuisParIsThere[p1]) AddNuisPar(p1);
    int idx0 = fNuisParIdx[p0];
    int idx1 = fNuisParIdx[p1];
    fMatrix[idx0][idx1] = corr;
}

//__________________________________________________________________________________
//
float CorrelationMatrix::GetCorrelation(string p0,string p1){
    if(!fNuisParIsThere[p0]){
        cout << "  WARNING: NP " << p0 << " not found in correlation matrix. Returning correlation = 0." << endl;
        return 0.;
    }
    if(!fNuisParIsThere[p1]){
        cout << "  WARNING: NP " << p1 << " not found in correlation matrix. Returning correlation = 0." << endl;
        return 0.;
    }
    int idx0 = fNuisParIdx[p0];
    int idx1 = fNuisParIdx[p1];
    return fMatrix[idx0][idx1];
}
