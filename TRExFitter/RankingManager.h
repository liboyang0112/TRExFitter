#ifndef RANKINGMANAGER_H
#define RANKINGMANAGER_H

#include <string>
#include <map>
#include <memory>
#include <vector>

class FitResults;
class FittingTool;
class RooDataSet;
class RooWorkspace;
class RooSimultaneous;
class NormFactor;

namespace RooStats {
    class ModelConfig;
}

class RankingManager {
public:

    struct RankingValues {
        double central;
        double up;
        double down;
    };

    explicit RankingManager();
    ~RankingManager();
    RankingManager(const RankingManager& r) = delete;
    RankingManager(RankingManager&& r) = delete;
    RankingManager& operator=(const RankingManager& r) = delete;
    RankingManager& operator=(RankingManager&& r) = delete;

    inline void SetOutputPath(const std::string& path){fOutputPath = path;}
    inline void SetInjectGlobalObservables(const bool flag){fInjectGlobalObservables = flag;}
    inline void SetNPValues(const std::map<std::string, double>& m){fFitValues = m;}
    inline void SetFixedNPs(const std::map<std::string, double>& m){fFitFixedNPs = m;}
    inline void SetFitStrategy(const int s){fFitStrategy = s;}
    inline void SetPOIName(const std::string& name){fPOIName = name;}
    inline void SetNCPU(const int n){fCPU = n;}
    inline void SetRng(const double range, const bool use, const int seed) {
        fRndRange = range;
        fUseRnd = use;
        fRndSeed = seed;
    }
    inline void SetStatOnly(const bool flag){fStatOnly = flag;}
    void AddNuisPar(const std::string& name, const bool isNF);

    void RunRanking(FitResults* fitResults,
                    RooWorkspace* ws,
                    RooDataSet* data,
                    const std::vector<std::shared_ptr<NormFactor> >& nfs) const;

    double RunSingleFit(FittingTool* fitTool,
                        RooWorkspace* ws,       
                        RooStats::ModelConfig *mc,
                        RooSimultaneous *simPdf,
                        RooDataSet* data,
                        const std::pair<std::string, bool>& np,
                        const bool isUp,
                        const bool isPrefit,
                        const RankingValues& values,
                        const double muhat) const;
private:

    std::string fOutputPath;
    std::vector<std::pair<std::string,bool> > fNuisPars;
    bool fInjectGlobalObservables;
    std::map<std::string, double> fFitValues;
    std::map<std::string, double> fFitFixedNPs;
    int fFitStrategy;
    std::string fPOIName;
    int fCPU;
    double fRndRange;
    bool fUseRnd;
    int fRndSeed;
    bool fStatOnly;
};

#endif
