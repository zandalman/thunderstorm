#ifndef PARSER_H
#define PARSER_H

// includes
#include <string>
#include <unordered_map>
#include <vector>

typedef std::unordered_map<std::string, std::string> Dict;
typedef std::unordered_map<std::string, Dict> Config;

std::string trim(const std::string &str);
int parseConfig(const std::string &filename, Config &config);

#endif
