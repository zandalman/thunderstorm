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
#include <chrono>
#include <mpi.h>

// headers
#include "main.h"
#include "parser.h"
#include "const.h"
#include "random.h"
#include "sim.h"

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto start = std::chrono::steady_clock::now();

  std::string config_dir = "../config/config.ini";
  Config config;
  parseConfig(config_dir, config);

  Vector1d ab;
  double time = std::stod(config["Simulation"]["ab_time"]) * constants::day;
  parseAb(config["IO"]["ab"], time, ab);

  EEDLData eedl;
  parseEEDL(config["IO"]["EEDL"], eedl);

  double ener = std::stod(config["Particle"]["ener"]);
  Part part = Part(0, constants::m_e, constants::e, ener);

  std::string outfile = config["IO"]["outfile"] + "." + std::to_string(rank);
  double tpart_max = std::stod(config["Particle"]["tpart_max"]);
  double enerpart_min = std::stod(config["Particle"]["enerpart_min"]);
  double rho = std::stod(config["Simulation"]["rho"]);
  double Bmag_turb = std::stod(config["Bfield"]["Bmag_turb"]);
  double Bmag_co = std::stod(config["Bfield"]["Bmag_co"]);
  double q = std::stod(config["Bfield"]["q"]);
  double Lmax = std::stod(config["Bfield"]["Lmax"]);
  Sim sim = Sim(part, eedl, ab, outfile, rho, Bmag_co, Bmag_turb, q, Lmax);

  int npac = std::stoi(config["Simulation"]["npac"]);
  int tsim_max = std::stoi(config["Simulation"]["tsim_max"]) * 3600;
  int count_loc = 0, count_glob = 0;

  clearOutfile(outfile);
  MPI_Barrier(MPI_COMM_WORLD);
  
  while ( true ) {
    
    // reset the simulation with a new particle
    int id = 1000 * count_loc + rank;
    Part part = Part(id, constants::m_e, constants::e, ener);
    sim.reset(part);

    // run the simulation
    if ( tpart_max > 0 ) {
      while ( sim.time < tpart_max ) sim.step();
    } else if ( enerpart_min > 0 ) {
      while ( sim.part.ener > enerpart_min ) sim.step();
    }
    writeEvent(sim.outfile, sim.event_list);
    sim.event_list.clear();

    // compute the global packet count
    count_loc++;
    MPI_Reduce(&count_loc, &count_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&count_glob, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ( npac > 0 && count_glob >= npac ) break;

    // compute the elapsed time
    auto now = std::chrono::steady_clock::now();
    auto tsim = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    if ( tsim_max > 0 && tsim >= tsim_max) break;

  }

  MPI_Finalize();
  return 0;
}