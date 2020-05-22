#include "TRExFitter/YamlConverter.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"

#include "yaml-cpp/include/yaml-cpp/yaml.h"

#include "TSystem.h"

#include <fstream>

YamlConverter::YamlConverter() :
    m_lumi("139"),
    m_cme("13000") {
}

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
    
void YamlConverter::WriteRankingHEPData(const std::vector<RankingContainer>& ranking,
                                        const std::string& folder) const {

    gSystem->mkdir((folder+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value <<  "parameter" << YAML::EndMap; 
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << irank.name;
                    out << YAML::EndMap;
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
            // NP value
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "NP value, error" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    AddValueErrors(out, irank.nphat, irank.nperrhi, irank.nperrlo);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
            
            // NP PosFit impact up
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI high" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << Common::KeepSignificantDigits(irank.poihi,2);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
            // NP Postfit imapct down
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI low" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << Common::KeepSignificantDigits(irank.poilo,2);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
            // NP PreFit impact up
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI prefit high" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << Common::KeepSignificantDigits(irank.poiprehi, 2);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
            // NP Prefit imapct down
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI prefit high" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << Common::KeepSignificantDigits(irank.poiprelo,2);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;
    out << YAML::EndMap;
    
    // Write to the file
    WriteInfoStatus("YamlConverter::WriteRankingHEPData", "Writing HEPData ranking yaml file to: " + folder + "/HEPData/Ranking.yaml");
    std::ofstream file;
    file.open((folder + "/HEPData/Ranking.yaml").c_str());
    if (!file.is_open()) {
        WriteWarningStatus("YamlConverter::WriteRankingHEPData", "Cannot open yaml file at: " + folder + "/HEPData/Ranking.yaml");
        return;
    }

    file << out.c_str();
    file.close();
}

void YamlConverter::AddQualifiers(YAML::Emitter& out) const {
    out << YAML::Key << "qualifiers";
    out << YAML::Value << YAML::BeginSeq;
    out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "SQRT(s)";
        out << YAML::Key << "units";
        out << YAML::Value << "GeV";
        out << YAML::Key << "value";
        out << YAML::Value << std::stoi(m_cme);
    out << YAML::EndMap;
    out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "LUMINOSITY";
        out << YAML::Key << "units";
        out << YAML::Value << "fb$^{-1}$";
        out << YAML::Key << "value";
        out << YAML::Value << std::stoi(m_lumi);
    out << YAML::EndMap;
    out << YAML::EndSeq;
}

void YamlConverter::AddValueErrors(YAML::Emitter& out,
                                   const double mean,
                                   const double up,
                                   const double down) const {
    
    double value = mean;
    if ((std::fabs(up) - std::fabs(down)) < 0.1) {
        // are symmetric
        double error = up;
        const int n = Common::ApplyATLASrounding(value, error);
        out << YAML::Key << "value";
        out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),value);
        out << YAML::Key << "errors";
        out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "symerror";
        if (n >= 0) {
            out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),error);
        } else {
            out << YAML::Value << Form("%.f",error);
        }
        out << YAML::EndMap;
        out << YAML::EndSeq;
    } else {
        double error = 0.5*(up-down);
        const int n = Common::ApplyATLASrounding(value, error);
        out << YAML::Key << "value";
        out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),value);
        out << YAML::Key << "errors";
        out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "asymerror";
        out << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "plus";
            out << YAML::Key << Form(("%."+std::to_string(n)+"f").c_str(),up);
            out << YAML::Key << "minus";
            out << YAML::Key << Form(("%."+std::to_string(n)+"f").c_str(),down);
        out << YAML::EndMap;
        out << YAML::EndMap;
        out << YAML::EndSeq;
    }
}
