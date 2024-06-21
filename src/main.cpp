// defines
#define _USE_MATH_DEFINES

// includes
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// headers
#include "main.h"
#include "parser.h"
#include "const.h"
#include "random.h"
#include "sim.h"

int main() {

  std::string config_dir = "../config/config.ini";
  Config config;
  parseConfig(config_dir, config);

  Vector1d ab;
  double time = std::stod(config["Simulation"]["ab_time"]) * constants::day;
  parseAb(config["IO"]["ab"], time, ab);

  EEDLData eedl;
  parseEEDL(config["IO"]["EEDL"], eedl);

  int id = 0;
  double ener = std::stod(config["Particle"]["ener"]);
  Part part = Part(id, constants::m_e, constants::e, ener);

  std::string outfile = config["IO"]["outfile"];
  int fsave = std::stoi(config["IO"]["fsave"]);
  double tmax = std::stod(config["Simulation"]["tmax"]);
  double rho = std::stod(config["Simulation"]["rho"]);
  double Bmag_turb = std::stod(config["Bfield"]["Bmag_turb"]);
  double Bmag_co = std::stod(config["Bfield"]["Bmag_co"]);
  Sim sim = Sim(part, eedl, ab, outfile, rho, Bmag_co, Bmag_turb);

  clearOutfile(outfile);
  while ( sim.time < tmax ) {
    sim.multistep(fsave);
  }
}