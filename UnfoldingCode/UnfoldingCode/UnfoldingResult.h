#ifndef UNFOLDINGRESULT_H_
#define UNFOLDINGRESULT_H_

#include <memory>
#include <vector>

class TGraphAsymmErrors;
class TH1;
class TH1D;

class UnfoldingResult {

public:

    struct FitValue {
        double nominal;
        double up;
        double down;
    };

    explicit UnfoldingResult();
    ~UnfoldingResult();
    UnfoldingResult(const UnfoldingResult& r) = delete;
    UnfoldingResult& operator=(const UnfoldingResult& r) = delete;
    UnfoldingResult(UnfoldingResult&& r) = delete;
    UnfoldingResult& operator=(UnfoldingResult&& r) = delete;

    inline void SetFitValues(const std::vector<FitValue>& v) {fFitValues = v;}
    inline const std::vector<FitValue>& GetFitValues() const {return fFitValues;}
    inline void AddFitValue(const FitValue& v) {fFitValues.emplace_back(v);}
    void AddFitValue(const double nominal, const double up, const double down);
    inline void ResetFitValues() {fFitValues.clear();}

    void SetTruthDistribution(const TH1* truth);
    inline const TH1D* GetTruthDistribution() const {return fTruthDistribution.get();}

    std::unique_ptr<TGraphAsymmErrors> GetUnfoldedResultErrorBand() const;
    std::unique_ptr<TH1D> GetUnfoldedResult() const;

private:

    std::vector<FitValue> fFitValues;
    std::unique_ptr<TH1D> fTruthDistribution;

};

#endif
