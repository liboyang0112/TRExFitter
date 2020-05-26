#ifndef YAMLCOVERTER_H_
#define YAMLCOVERTER_H_

#include <fstream>
#include <string>
#include <vector>

namespace YAML {
    class Emitter;
}

class TGraphAsymmErrors;

class YamlConverter {

public:
    struct RankingContainer {
        std::string name;
        double nphat;
        double nperrhi;
        double nperrlo;
        double poihi;
        double poilo;
        double poiprehi;
        double poiprelo;
    };

    struct TableContainer {
        std::vector<std::string> regionNames;
        std::vector<std::string> sampleNames;
        std::vector<std::vector<double> > mcYields;
        std::vector<std::vector<double> > mcErrors;
        std::vector<std::string> dataNames;
        std::vector<std::vector<double> > dataYields;
    };

    struct PlotContainer {
        std::vector<std::string> samples;
        std::vector<std::vector<double> > signalYields;
        std::vector<std::vector<double> > backgroundYields;
        TGraphAsymmErrors* errors;
        std::vector<double> data;
        std::vector<int> blindedBins;
        std::string xAxis;
        std::string yAxis;
        std::string region;
    };

    struct SubmissionContainer {
        bool useTables;
        bool isUnfolding;
        std::string folder;
        std::vector<std::string> regionNames;
    };

    explicit YamlConverter();
    ~YamlConverter() = default;

    YamlConverter(const YamlConverter& c) = delete;
    YamlConverter& operator=(const YamlConverter& c) = delete;
    YamlConverter(const YamlConverter&& c) = delete;
    YamlConverter& operator=(const YamlConverter&& c) = delete;

    void SetLumi(const std::string& lumi){ m_lumi = lumi; }
    void SetCME(const std::string& cme){ m_cme = cme; }

    void WriteRanking(const std::vector<RankingContainer>& ranking,
                      const std::string& path) const;
    
    void WriteRankingHEPData(const std::vector<RankingContainer>& ranking,
                             const std::string& folder) const;
    
    void WriteCorrelation(const std::vector<std::string>& np,
                          const std::vector<std::vector<double> >& corr,
                          const std::string& path) const;
    
    void WriteCorrelationHEPData(const std::vector<std::string>& np,
                                 const std::vector<std::vector<double> >& corr,
                                 const std::string& path) const;
    
    void WriteTables(const TableContainer& container,
                     const std::string& directory,
                     const bool isPostFit) const;

    void WriteTablesHEPData(const TableContainer& container,
                            const std::string& directory,
                            const bool isPostFit) const;

    void WriteUnfolding(const TGraphAsymmErrors* const graph,
                        const std::string& directory) const;
    
    void WriteUnfoldingHEPData(const TGraphAsymmErrors* const graph,
                               const std::string& xAxis,
                               const std::string& directory) const;

    void WritePlot(const PlotContainer& container,
                   const std::string& directory,
                   const bool isPostFit) const;

    void WritePlotHEPData(const PlotContainer& container,
                          const std::string& directory,
                          const bool isPostFit) const;

    void WriteHEPDataSubmission(const SubmissionContainer& container) const;

private:
    std::string m_lumi;
    std::string m_cme;

    void AddQualifiers(YAML::Emitter& out) const;

    void AddValueErrors(YAML::Emitter& out,
                       const double mean,
                       const double up,
                       const double down) const;

    void Write(const YAML::Emitter& out, const std::string& type, const std::string& path) const;

    bool TableContainerIsOK(const TableContainer& container) const;
    
    bool PlotContainerIsOK(const PlotContainer& container) const;
};

#endif
