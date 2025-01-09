#ifndef PARSER_H
#define PARSER_H

// includes
#include <string>
#include <vector>
#include <unordered_map>

typedef std::unordered_map<std::string, std::string> Dict;
typedef std::unordered_map<std::string, Dict> Config;

std::string trim(const std::string &str);
void parseConfig(const std::string &configfile_name, Config &config);
void makeList(Dict &dict, std::vector<double> &list, size_t &num, double unit = 1.0);

#endif
