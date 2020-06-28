#ifndef LIKELIHOOHSVANMANAGER_H
#define LIKELIHOOHSVANMANAGER_H

#include <memory>
#include <string>
#include <vector>

class RooDataSet;
class RooWorkspace;
class NormFactor;

class LikelihoodScanManager {
public:

    explicit LikelihoodScanManager();
    ~LikelihoodScanManager();

    LikelihoodScanManager(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager(LikelihoodScanManager&& m) = delete;
    LikelihoodScanManager& operator=(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager& operator=(LikelihoodScanManager&& m) = delete;

    inline void SetScanParamsX(const double min, const double max, const int steps) {
        fScanMinX = min;
        fScanMaxX = max;
        fStepsX = steps;
    }
    
    inline void SetScanParamsY(const double min, const double max, const int steps) {
        fScanMinY = min;
        fScanMaxY = max;
        fStepsY = steps;
    }

    inline void SetNCPU(const int n){fCPU = n;}
    inline void SetOffSet(const bool flag){fUseOffset = flag;}

    typedef std::pair<std::vector<double>, std::vector<double> > scanResult1D;

    scanResult1D Run1DScan(const RooWorkspace* ws,
                           const std::string& varName,
                           RooDataSet* data,
                           const std::vector<std::shared_ptr<NormFactor> >& nfs) const;

private:
    double fScanMinX;
    double fScanMinY;
    int fStepsX;
    double fScanMaxX;
    double fScanMaxY;
    int fStepsY;
    bool fUseOffset;
    int fCPU;
};

#endif
