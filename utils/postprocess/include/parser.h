#ifndef PARSER_H
#define PARSER_H

// includes
#include <string>
#include <unordered_map>
#include <vector>

typedef std::unordered_map<std::string, std::string> Dict;
typedef std::unordered_map<std::string, Dict> Config;

std::string trim(const std::string &str);
void parseConfig(const std::string &configfile_name, Config &config);
std::vector<std::string> splitString(const std::string& str, char delimiter);
void makeList(Dict &dict, double unit, size_t &num, std::vector<double> &list);

#endif
