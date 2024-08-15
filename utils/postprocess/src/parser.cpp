// includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <mpi.h>

// headers
#include "parser.h"
#include "functions.h"

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
 * @param configfile_name The configuration file name.
 * @param config          The structure to store the parsed data.
*/
void parseConfig(const std::string &configfile_name, Config &config) {
  std::ifstream configfile(configfile_name);
  std::string line;
  std::string section;

  if ( !configfile ) {
    std::cerr << "Error: Failed to open file " << configfile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  while ( std::getline(configfile, line) ) {
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

  configfile.close();
};

/**
 * @brief Make a list using data from the configuration file.
 * 
 * @param dict The dictionary containing the bin information.
 * @param list The vector to store the list.
 * @param num  The number of elements in the list.
 * @param unit The data unit.
 */
void makeList(Dict &dict, std::vector<double> &list, size_t &num, double unit) {
  double min = std::stod(dict["min"]) * unit;
  double max = std::stod(dict["max"]) * unit;
  bool log = dict["log"] == "true";
  num = std::stoul(dict["num"]);
  linspace(min, max, num, log, list);
} 
