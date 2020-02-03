#ifndef PruningUtil_H
#define PruningUtil_H

#include "TH1F.h"
#include <iostream>

class PruningUtil {
public:
    explicit PruningUtil();
    ~PruningUtil() = default;

    PruningUtil(const PruningUtil& p) = delete;
    PruningUtil(PruningUtil&& p) = delete;
    PruningUtil& operator=(const PruningUtil& p) = delete;
    PruningUtil& operator=(PruningUtil&& p) = delete;

    // 0 = sample-by-sample, 1 = relative to tot background, 2 = relative to tot S+B
    void SetStrategy(const int strat);
    inline int GetStrategy() const {return fStrategy;}
    void SetThresholdNorm(const double thres);
    void SetThresholdShape(const double thres);
    void SetThresholdIsLarge(const double thres);

    // if hTot is not set, just prun w.r.t. only this nominal,
    // otherwise use hTot and do relative pruning
    // return scheme:
    // 0 : nothing pruned
    // 1 : shape pruned
    // 2 : norm pruned
    // 3 : all pruned
    // -2 : bad norm
    // -3 : bad shape
    // -4 : all bad
    int CheckSystPruning(const TH1* const hUp,
                         const TH1* const hDown,
                         const TH1* const hNom,
                         const TH1* hTot = nullptr);

    bool HasShapeRelative(const TH1* const hNom,
                          const TH1* const hUp,
                          const TH1* const hDown,
                          const TH1* const combined,
                          const double threshold) const;

 private:
    int fStrategy;
    double fThresholdNorm;
    double fThresholdShape;
    double fThresholdIsLarge;
};

#endif
