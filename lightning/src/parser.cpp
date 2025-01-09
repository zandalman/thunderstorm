// includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <mpi.h>

// headers
#include "json.h"
#include "parser.h"
#include "functions.h"

// types
using json = nlohmann::json;

/// @brief A constructor for the SpecData structure.
SpecData::SpecData(int Zelem_, std::string symbol_) 
  : Zelem(Zelem_)   // The proton number.
  , symbol(symbol_) // The symbol.
  {}

/**
 * @brief Remove spaces at the beginning and end of a string.
 * 
 * @param str The string.
 * @return The trimmed string.
*/
std::string trim(const std::string &str) {
  size_t first = str.find_first_not_of(' ');
  if (std::string::npos == first) return str;
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

/**
 * @brief Parse the configuration file.
 * 
 * @param filename The configuration file name.
 * @param config   The structure to store the parsed data.
*/
void parseConfig(const std::string &filename, Config &config) {
  std::ifstream file(filename);
  std::string line;
  std::string section;

  // return 1 on error
  if (!file) {
    std::cerr << "Error: Failed to open file " << filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  while (std::getline(file, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;
    if (line[0] == '[' && line[line.length() - 1] == ']') {
      section = line.substr(1, line.length() - 2);
      continue;
    }
    // parse key-value pairs and store them in config
    size_t pos_delimit = line.find('=');
    if (pos_delimit == std::string::npos)
      continue;
    std::string key = trim(line.substr(0, pos_delimit));
    std::string value = trim(line.substr(pos_delimit + 1));
    config[section][key] = value;
  }

  file.close();
};

/**
 * @brief Parse the abundance data.
 * 
 * @param filename The abundance data filename.
 * @param time     The time at which to sample the abundance data.
 * @param ab       The structure to store the parsed data.
 * @return 0 if success, 1 if failure.   
*/
void parseAb(const std::string &filename, double time, Vector1d &ab) {
  std::ifstream file(filename);

  if (!file) {
    std::cerr << "Error: Failed to open file " << filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // read json file
  json ab_json;
  file >> ab_json;

  // find the index corresponding to the desired time
  Vector1d time_list = ab_json["time"].get<Vector1d>();
  int idx = findIdx(time, time_list);

  // retrieve the abundance data
  ab = ab_json["ab"][idx].get<Vector1d>();

  file.close();
}

/**
 * @brief Parse the EEDL data.
 * 
 * @param filename The EEDL data filename.
 * @param eedl     The structure to store the parsed data.
*/
void parseEEDL(const std::string &filename, EEDLData &eedl) {
  std::ifstream file(filename);

  if ( !file ) {
    std::cerr << "Error: Failed to open file " << filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // read json file
  json eedl_json;
  file >> eedl_json;

  Vector1d1d sig_ion_data_1ion;
  Vector1d2d2d spec_ion_data_1ion;
  
  std::vector<int> Z_list = eedl_json["Z"].get<std::vector<int>>();
  std::vector<std::string> symbol_list =
      eedl_json["symbol"].get<std::vector<std::string>>();
  for (size_t i = 0; i < Z_list.size(); i++) {
    auto spec_json = eedl_json["data"][i];
    SpecData spec_data(Z_list[i], symbol_list[i]);
    spec_data.sig_tot_data.first      = spec_json["sig_tot"]["ener"].get<Vector1d>();
    spec_data.sig_tot_data.second     = spec_json["sig_tot"]["sig"].get<Vector1d>();
    spec_data.sig_scat_data.first     = spec_json["sig_scat_la"]["ener"].get<Vector1d>();
    spec_data.sig_scat_data.second    = spec_json["sig_scat_la"]["sig"].get<Vector1d>();
    spec_data.sig_brem_data.first     = spec_json["sig_brem"]["ener"].get<Vector1d>();
    spec_data.sig_brem_data.second    = spec_json["sig_brem"]["sig"].get<Vector1d>();
    spec_data.sig_exc_data.first      = spec_json["sig_exc"]["ener"].get<Vector1d>();
    spec_data.sig_exc_data.second     = spec_json["sig_exc"]["sig"].get<Vector1d>();
    spec_data.sig_ion_tot_data.first  = spec_json["sig_ion"]["ener"].get<Vector1d>();
    spec_data.sig_ion_tot_data.second = spec_json["sig_ion"]["sig"].get<Vector1d>();
    spec_data.th_scat_data.first      = spec_json["th_scat"]["ener"].get<Vector1d>();
    spec_data.th_scat_data.second     = spec_json["th_scat"]["cos_th"].get<Vector2d>();
    spec_data.th_scat_data.third      = spec_json["th_scat"]["cos_th_dist"].get<Vector2d>();
    spec_data.spec_brem_data.first    = spec_json["spec_brem"]["ener"].get<Vector1d>();
    spec_data.spec_brem_data.second   = spec_json["spec_brem"]["ener_loss"].get<Vector1d>();
    spec_data.spec_exc_data.first     = spec_json["spec_exc"]["ener"].get<Vector1d>();
    spec_data.spec_exc_data.second    = spec_json["spec_exc"]["ener_loss"].get<Vector1d>();
    spec_data.ion_list                = spec_json["ion_list"].get<std::vector<std::string>>();
    spec_data.ener_bind_list          = spec_json["ener_bind_list"].get<Vector1d>();
    for (size_t j = 0; j < spec_data.ion_list.size(); j++) {
      sig_ion_data_1ion.first   = spec_json["sig_ion_list"][j]["ener"].get<Vector1d>();
      sig_ion_data_1ion.second  = spec_json["sig_ion_list"][j]["sig"].get<Vector1d>();
      spec_ion_data_1ion.first  = spec_json["spec_ion_list"][j]["ener"].get<Vector1d>();
      spec_ion_data_1ion.second = spec_json["spec_ion_list"][j]["ener_loss"].get<Vector2d>();
      spec_ion_data_1ion.third  = spec_json["spec_ion_list"][j]["ener_loss_dist"].get<Vector2d>();
      spec_data.sig_ion_data.push_back(sig_ion_data_1ion);
      spec_data.spec_ion_data.push_back(spec_ion_data_1ion);
    }
    eedl.push_back(spec_data);
  }
  file.close();
}
