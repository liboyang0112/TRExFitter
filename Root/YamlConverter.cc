#include "TRExFitter/YamlConverter.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/StatusLogbook.h"

#include "yaml-cpp/include/yaml-cpp/yaml.h"

#include "TGraphAsymmErrors.h"
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
    Write(out, "ranking", path);
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
    Write(out, "HEPData ranking", folder + "/HEPData/Ranking.yaml");
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
    if ((std::fabs(down) > 1e-6) && (std::fabs(up/down) > 0.9) && (std::fabs(up/down) < 1.1)) {
        double error = 0.5*(std::fabs(up) + std::fabs(down)) ;
        const int n = Common::ApplyATLASrounding(value, error);
        // are symmetric
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
        double error = std::min(std::fabs(up),std::fabs(down));
        const int n = Common::ApplyATLASrounding(value, error);
        out << YAML::Key << "value";
        out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),value);
        out << YAML::Key << "errors";
        out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "asymerror";
        out << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "plus";
            if (n >= 0) {
                out << YAML::Key << Form(("%."+std::to_string(n)+"f").c_str(),up);
            } else {
                out << YAML::Key << Form("%.f",up);
            }
            out << YAML::Key << "minus";
            if (n >= 0) {
                out << YAML::Key << Form(("%."+std::to_string(n)+"f").c_str(),down);
            } else {
                out << YAML::Key << Form("%.f",down);
            }
        out << YAML::EndMap;
        out << YAML::EndMap;
        out << YAML::EndSeq;
    }
}


void YamlConverter::Write(const YAML::Emitter& out, const std::string& type, const std::string& path) const {
    WriteInfoStatus("YamlConverter::Write", "Writing " + type + " yaml file to: " + path);
    std::ofstream file;
    file.open(path.c_str());
    if (!file.is_open() || !file.good()) {
        WriteWarningStatus("YamlConverter::Write", "Cannot open yaml file at: " + path);
        return;
    }

    file << out.c_str();
    file.close();
}

void YamlConverter::WriteCorrelation(const std::vector<std::string>& np,
                                     const std::vector<std::vector<double> >& corr,
                                     const std::string& path) const {

    const std::size_t n = np.size();
    if (corr.size() != n) {
        WriteWarningStatus("YamlConverter::WriteCorrelation", "Inconsistent inputs!");
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginSeq;
    out << YAML::BeginMap;
    out << YAML::Key << "parameters";
    out << YAML::Value << YAML::BeginSeq;
    for (const auto& iname : np) {
        out << iname;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    out << YAML::BeginMap;
    out << YAML::Key << "correlation_rows";
    out << YAML::Value << YAML::BeginSeq;
    for (const auto& icorr : corr) {
        if (icorr.size() != n) {
            WriteWarningStatus("YamlConverter::WriteCorrelation", "Inconsistent inputs for correlation!");
            return;
        }
        out << YAML::Flow << YAML::BeginSeq;
        for (const auto& i : icorr) {
            out << Form("%.4f",i);
        }
        out << YAML::EndSeq;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    out << YAML::EndSeq;

    Write(out, "correlation", path+"/CorrelationMatrix.yaml");
}
    
void YamlConverter::WriteCorrelationHEPData(const std::vector<std::string>& np,
                                            const std::vector<std::vector<double> >& corr,
                                            const std::string& folder) const {

    gSystem->mkdir((folder+"/HEPData").c_str());

    const std::size_t n = np.size();
    if (corr.size() != n) {
        WriteWarningStatus("YamlConverter::WriteCorrelationHEPData", "Inconsistent inputs!");
        return;
    }

    if (n > 100) {
        WriteInfoStatus("YamlConverter::WriteCorrelationHEPData", "Processing matrix of more than 100x100 elements, this may take some time...");
    }
    
    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "dependent_variables";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "NP correlations" << YAML::EndMap;
        AddQualifiers(out);
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (const auto& icorr : corr) {
            if (icorr.size() != n) {
                WriteWarningStatus("YamlConverter::WriteCorrelationHEPData", "Inconsistent inputs for correlation!");
                return;
            }
            for (const auto& i : icorr) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << YAML::Value << Form("%.2f", i);
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::Key << "independent_variables";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "NPs" << YAML::EndMap;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t inp = 0; inp < n; ++inp) {
            for (std::size_t jnp = 0; jnp < n; ++jnp) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << np.at(inp);
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "NPs" << YAML::EndMap;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t inp = 0; inp < n; ++inp) {
            for (std::size_t jnp = 0; jnp < n; ++jnp) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << np.at(jnp);
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::EndMap;

    Write(out, "HEPData correlation", folder+"/HEPData/Correlation.yaml");

}
    
void YamlConverter::WriteTables(const YamlConverter::TableContainer& container,
                                const std::string& directory,
                                const bool isPostFit) const {

    if (!YamlConverter::TableContainerIsOK(container)) {
        WriteWarningStatus("YamlConverter::WriteTables", "Inconsistent inputs for tables!");
        return;
    } 
    gSystem->mkdir((directory+"/HEPData").c_str());
    
    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (std::size_t ireg = 0; ireg < container.regionNames.size(); ++ireg) {
        out << YAML::BeginMap;
        out << YAML::Key << "Region";
        out << YAML::Value << container.regionNames.at(ireg);
        out << YAML::Key << "Samples";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t isample = 0; isample < container.sampleNames.size(); ++isample) {
            out << YAML::BeginMap;
                out << YAML::Key << "Sample";
                out << YAML::Value << container.sampleNames.at(isample);
                out << YAML::Key << "Yield";
                out << YAML::Value << container.mcYields.at(isample).at(ireg);
                out << YAML::Key << "Error";
                out << YAML::Value << container.mcErrors.at(isample).at(ireg);
            out << YAML::EndMap;
        }
        for (std::size_t idata = 0; idata < container.dataNames.size(); ++idata) {
            out << YAML::BeginMap;
                out << YAML::Key << "Data";
                out << YAML::Value << container.dataNames.at(idata);
                out << YAML::Key << "Yield";
                out << YAML::Value << container.dataYields.at(idata).at(ireg);
            out << YAML::EndMap;
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;
    // Write to the file
    if (isPostFit) {
        Write(out, "postfit yield tables", directory + "/Tables/Table_postfit.yaml");
    } else {
        Write(out, "prefit yield tables", directory + "/Tables/Table_prefit.yaml");
    }
}

    
void YamlConverter::WriteTablesHEPData(const YamlConverter::TableContainer& container,
                                       const std::string& directory,
                                       const bool isPostFit) const {

    if (!YamlConverter::TableContainerIsOK(container)) {
        WriteWarningStatus("YamlConverter::WriteTablesHEPData", "Inconsistent inputs for tables!");
        return;
    } 
    
    gSystem->mkdir((directory+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value <<  "process" << YAML::EndMap; 
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& isample : container.sampleNames) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << isample;
                    out << YAML::EndMap;
                }
                for (const auto& idata : container.dataNames) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << idata;
                    out << YAML::EndMap;
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;
        
        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
        // loop over regions
        for (std::size_t ireg = 0; ireg < container.regionNames.size(); ++ireg) {
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << container.regionNames.at(ireg) << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (std::size_t ivalue = 0; ivalue < container.mcYields.size(); ++ivalue) {
                    out << YAML::BeginMap;
                    AddValueErrors(out, container.mcYields.at(ivalue).at(ireg), container.mcErrors.at(ivalue).at(ireg), container.mcErrors.at(ivalue).at(ireg));
                    out << YAML::EndMap;
                }
                for (std::size_t ivalue = 0; ivalue < container.dataYields.size(); ++ivalue) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << Form("%.f",container.dataYields.at(ivalue).at(ireg));
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        }
        out << YAML::EndSeq;
    out << YAML::EndMap;
    
    // Write to the file
    if (isPostFit) {
        Write(out, "HEPData postfit yield tables", directory + "/HEPData/Table_postfit.yaml");
    } else {
        Write(out, "HEPData prefit yield tables", directory + "/HEPData/Table_prefit.yaml");
    }
}

bool YamlConverter::TableContainerIsOK(const YamlConverter::TableContainer& container) const {

    const std::size_t nRegions = container.regionNames.size();
    const std::size_t nSamples = container.sampleNames.size();

    if (nRegions == 0 || nSamples == 0) return false;
    if (container.mcYields.size() != nSamples) return false;
    if (container.mcErrors.size() != nSamples) return false;
    for (const auto& ivec : container.mcYields) {
        if (ivec.size() != nRegions) return false;
    }
    for (const auto& ivec : container.mcErrors) {
        if (ivec.size() != nRegions) return false;
    }
    for (const auto& ivec : container.dataYields) {
        if (ivec.size() != nRegions) return false;
    }
    
    return true;
}
    
void YamlConverter::WriteUnfolding(const TGraphAsymmErrors* const graph,
                                   const std::string& directory) const {

    const int n = graph->GetN();
    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (int i = 0; i < n; ++i) {
        double x;
        double y;
        graph->GetPoint(i, x, y);

        const double x_min = graph->GetErrorXlow(i);
        const double x_max = graph->GetErrorXhigh(i);
        const double y_min = graph->GetErrorYlow(i);
        const double y_max = graph->GetErrorYhigh(i);
        out << YAML::BeginMap;
            out << YAML::Key << "range";
            out << YAML::Value << YAML::Flow << YAML::BeginSeq << x-x_min << x+x_max << YAML::EndSeq;
            out << YAML::Key << "mean";
            out << YAML::Value << y;
            out << YAML::Key << "uncertaintyUp";
            out << YAML::Value << y_max;
            out << YAML::Key << "uncertaintyDown";
            out << YAML::Value << -y_min;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    Write(out, "unfolding result", directory + "/UnfoldingData.yaml");
}
    
void YamlConverter::WriteUnfoldingHEPData(const TGraphAsymmErrors* const graph,
                                          const std::string& directory) const {
}

