// includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// headers
#include "parser.h"
#include "json.h"

// types
using json = nlohmann::json;

SpecData::SpecData(int Z_, std::string symbol_): Z(Z_), symbol(symbol_) {}

std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first) return str;
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

Config parseConfig(const std::string& filename) {
    std::ifstream file(filename);
    Config config;
    std::string line;
    std::string section;

    // return an empty map on error
    if (!file) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return config;
    }

    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line[0] == '[' && line[line.length() - 1] == ']') {
            section = line.substr(1, line.length() - 2);
            continue;
        }
        // parse key-value pairs and store them in config
        size_t pos_delimit = line.find('=');
        if (pos_delimit == std::string::npos) continue;
        std::string key = trim(line.substr(0, pos_delimit));
        std::string value = trim(line.substr(pos_delimit + 1));
        config[section][key] = value;
    }

    file.close();
    return config;
};

int parseEEDL(const std::string& filename, EEDLData& eedl_data) {
    std::ifstream file(filename);

    // return 1 on error
    if (!file) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return 1;
    }
    
    // read json file
    json eedl_json;
    file >> eedl_json;

    Vector1d1d sig_ion_data_1ion;
    Vector1d2d2d spec_ion_data_1ion;

    std::vector<int> Z_list = eedl_json["Z"].get<std::vector<int>>();
    std::vector<std::string> symbol_list = eedl_json["symbol"].get<std::vector<std::string>>();
    for ( size_t i = 0; i < Z_list.size(); i++ ) {
        auto spec_json = eedl_json["data"][i];
        SpecData spec_data(Z_list[i], symbol_list[i]);
        spec_data.sig_scat_data.first = spec_json["sig_scat"]["ener"].get<Vector1d>();
        spec_data.sig_scat_data.second = spec_json["sig_scat"]["sig"].get<Vector1d>();
        spec_data.sig_scat_la_data.first = spec_json["sig_scat_la"]["ener"].get<Vector1d>();
        spec_data.sig_scat_la_data.second = spec_json["sig_scat_la"]["sig"].get<Vector1d>();
        spec_data.sig_brem_data.first = spec_json["sig_brem"]["ener"].get<Vector1d>();
        spec_data.sig_brem_data.second = spec_json["sig_brem"]["sig"].get<Vector1d>();
        spec_data.sig_exc_data.first = spec_json["sig_exc"]["ener"].get<Vector1d>();
        spec_data.sig_exc_data.second = spec_json["sig_exc"]["sig"].get<Vector1d>();
        spec_data.th_scat_data.first = spec_json["th_scat"]["ener"].get<Vector1d>();
        spec_data.th_scat_data.second = spec_json["th_scat"]["cos_th"].get<Vector2d>();
        spec_data.th_scat_data.third = spec_json["th_scat"]["cos_th_dist"].get<Vector2d>();
        spec_data.spec_brem_data.first = spec_json["spec_brem"]["ener"].get<Vector1d>();
        spec_data.spec_brem_data.second = spec_json["spec_brem"]["ener_loss"].get<Vector1d>();
        spec_data.spec_exc_data.first = spec_json["spec_exc"]["ener"].get<Vector1d>();
        spec_data.spec_exc_data.second = spec_json["spec_exc"]["ener_loss"].get<Vector1d>();
        spec_data.ion_list = spec_json["ion_list"].get<std::vector<std::string>>();
        for ( size_t j = 0; j < spec_data.ion_list.size(); j++ ) {
            sig_ion_data_1ion.first = spec_json["sig_ion_list"][j]["ener"].get<Vector1d>();
            sig_ion_data_1ion.second = spec_json["sig_ion_list"][j]["sig"].get<Vector1d>();
            spec_data.sig_ion_data.push_back(sig_ion_data_1ion);
            spec_ion_data_1ion.first = spec_json["spec_ion_list"][j]["ener"].get<Vector1d>();
            spec_ion_data_1ion.second = spec_json["spec_ion_list"][j]["ener_loss"].get<Vector2d>();
            spec_ion_data_1ion.third = spec_json["spec_ion_list"][j]["ener_loss_dist"].get<Vector2d>();
            spec_data.spec_ion_data.push_back(spec_ion_data_1ion);
        }
        eedl_data.push_back(spec_data);
    }

    return 0;

}
