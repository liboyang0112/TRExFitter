#ifndef YAMLCOVERTER_H_
#define YAMLCOVERTER_H_

#include <fstream>
#include <string>
#include <vector>

namespace YAML {
    class Emitter;
}

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
    

private:
    std::string m_lumi;
    std::string m_cme;

    void AddQualifiers(YAML::Emitter& out) const;

    void AddValueErrors(YAML::Emitter& out,
                       const double mean,
                       const double up,
                       const double down) const;
};

#endif
