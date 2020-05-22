#include "TRExFitter/YamlConverter.h"

#include "TRExFitter/StatusLogbook.h"

#include "yaml-cpp/include/yaml-cpp/yaml.h"

#include <fstream>

void YamlConverter::WriteRanking(const std::vector<YamlConverter::RankingContainer>& ranking,
                                 const std::string& path) const {

    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (const auto& irank : ranking) {
        out << YAML::BeginMap;
            out << YAML::Key << "Name";
            out << YAML::Value << irank.name;
            out << YAML::Key << "NPhat";
            out << YAML::Value << irank.nphat;
            out << YAML::Key << "NPerrHi";
            out << YAML::Value << irank.nperrhi;
            out << YAML::Key << "NPerrLo";
            out << YAML::Value << irank.nperrlo;
            out << YAML::Key << "POIup";
            out << YAML::Value << irank.poihi;
            out << YAML::Key << "POIdown";
            out << YAML::Value << irank.poilo;
            out << YAML::Key << "POIupPreFit";
            out << YAML::Value << irank.poiprehi;
            out << YAML::Key << "POIdownPreFit";
            out << YAML::Value << irank.poiprelo;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    // Write to the file
    WriteInfoStatus("YamlConverter::WriteRanking", "Writing ranking yaml file to: " + path);
    std::ofstream file;
    file.open(path.c_str());
    if (!file.is_open()) {
        WriteWarningStatus("YamlConverter::WriteRanking", "Cannot open yaml file at: " + path);
        return;
    }

    file << out.c_str();
    file.close();
}
