#ifndef PARSER_H
#define PARSER_H

// includes
#include <string>
#include <unordered_map>
#include <vector>

// types
typedef std::vector<double> Vector1d;
typedef std::pair<Vector1d, Vector1d> Vector1d1d;
typedef std::vector<Vector1d> Vector2d;
typedef std::pair<Vector2d, Vector2d> Vector2d2d;
typedef std::vector<Vector1d1d> Vector1d1dVector;
typedef std::unordered_map<std::string, std::string> Dict;
typedef std::unordered_map<std::string, Dict> Config;

/// @brief A structure to hold one 1D vector and two 2D vectors
struct Vector1d2d2d {
  Vector1d first;
  Vector2d second, third;
};
typedef std::vector<Vector1d2d2d> Vector1d2d2dVector;

/// @brief A structure to represent the EEDL data for a single species.
struct SpecData {
  int Zelem;                         // The proton number.
  std::string symbol;                // The symbol.
  std::vector<std::string> ion_list; // The ionization name data.
  Vector1d ener_bind_list;           // The ionization binding energy data.
  Vector1d1d sig_tot_data;           // The total cross section data.
  Vector1d1d sig_scat_data;          // The scattering cross section data.
  Vector1d1d sig_brem_data;          // The Bremsstrahlung cross section data.
  Vector1d1d sig_exc_data;           // The excitation cross section data.
  Vector1d1d sig_ion_tot_data;       // The total ionization cross section data.
  Vector1d1dVector sig_ion_data;     // The ionization cross section data.
  Vector1d2d2d th_scat_data;         // The scattering angle data.
  Vector1d1d spec_brem_data;         // The Bremsstrahlung energy loss data.
  Vector1d1d spec_exc_data;          // The excitation energy loss data.
  Vector1d2d2dVector spec_ion_data;  // The ionization energy loss data.
  
  SpecData() = default;
  SpecData(int Zelem_, std::string symbol_);
};
typedef std::vector<SpecData> EEDLData;

std::string trim(const std::string &str);
int parseConfig(const std::string &filename, Config &config);
int parseAb(const std::string &filename, double time, Vector1d &ab);
int parseEEDL(const std::string &filename, EEDLData &eedl_data);

#endif
