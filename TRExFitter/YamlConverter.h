#ifndef YAMLCOVERTER_H_
#define YAMLCOVERTER_H_

#include <string>
#include <vector>

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

    explicit YamlConverter() = default;
    ~YamlConverter() = default;

    YamlConverter(const YamlConverter& c) = delete;
    YamlConverter& operator=(const YamlConverter& c) = delete;
    YamlConverter(const YamlConverter&& c) = delete;
    YamlConverter& operator=(const YamlConverter&& c) = delete;

    void WriteRanking(const std::vector<RankingContainer>& ranking,
                      const std::string& path) const;
    

private:

};

#endif
