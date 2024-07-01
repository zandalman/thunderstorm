// includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// headers
#include "parser.h"

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
int parseConfig(const std::string &filename, Config &config) {
  std::ifstream file(filename);
  std::string line;
  std::string section;

  // return 1 on error
  if (!file) {
    std::cerr << "Error: Failed to open file " << filename << std::endl;
    // MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
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
  return 0;
};