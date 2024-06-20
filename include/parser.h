#ifndef PARSER_H
#define PARSER_H

// includes
#include <string>
#include <vector>
#include <unordered_map>

// types
typedef std::vector<double> Vector1d;
typedef std::pair<Vector1d, Vector1d> Vector1d1d;
typedef std::vector<Vector1d> Vector2d;
typedef std::pair<Vector2d, Vector2d> Vector2d2d;
typedef std::vector<Vector1d1d> Vector1d1dVector;
typedef std::unordered_map<std::string, std::unordered_map<std::string, std::string>> Config;

struct Vector1d2d2d {
    Vector1d first;
    Vector2d second, third;
};
typedef std::vector<Vector1d2d2d> Vector1d2d2dVector;

struct SpecData {
    int Z;
    std::string symbol;
    std::vector<std::string> ion_list;
    Vector1d1d sig_scat_data;
    Vector1d1d sig_scat_la_data;
    Vector1d1d sig_brem_data;
    Vector1d1d sig_exc_data;
    Vector1d1dVector sig_ion_data;
    Vector1d2d2d th_scat_data;
    Vector1d1d spec_brem_data;
    Vector1d1d spec_exc_data;
    Vector1d2d2dVector spec_ion_data;
    SpecData() = default;
    SpecData(int Z_, std::string symbol_);
};
typedef std::vector<SpecData> EEDLData;

std::string trim(const std::string& str);
Config parseConfig(const std::string& filename);
int parseEEDL(const std::string& filename, EEDLData& eedl_data);

#endif
