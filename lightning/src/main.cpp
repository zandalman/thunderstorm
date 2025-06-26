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
#include "parser.h"
#include "const.h"
#include "random.h"
#include "sim.h"
#include "functions.h"

// types
typedef std::vector<double> Vector1d;

int main(int argc, char** argv) {

  // initialize MPI
  MPI_Init(&argc, &argv);

  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // start timer
  auto start = std::chrono::steady_clock::now();

  // read config file
  std::string config_dir = "../config/config.ini";
  Config config;
  parseConfig(config_dir, config);
  if ( rank == 0 ) std::cout << "Read config file." << std::endl;
  
  // read abundance data
  Vector1d ab;
  double mmw;
  double time = std::stod(config["Background"]["ab_time"]) * constants::day;
  parseAb(config["IO"]["ab"], time, ab, mmw);
  if ( rank == 0 ) std::cout << "Read abundance data." << std::endl;
  
  // read EEDL data
  EEDLData eedl;
  parseEEDL(config["IO"]["EEDL"], eedl);
  if ( rank == 0 ) std::cout << "Read EEDL data." << std::endl;

  // read the config file
  int tsim = std::stoi(config["Simulation"]["tsim"]);
  int count_max = std::stoi(config["Simulation"]["count_max"]);
  double ener_min = std::stod(config["Simulation"]["ener_min"]);
  bool cont = config["Simulation"]["continue"] == "true";
  double ener = std::stod(config["Particle"]["ener"]);
  double tmax = std::stod(config["Particle"]["tmax"]);
  double rho = std::stod(config["Background"]["rho"]);
  double temp = std::stod(config["Background"]["temp"]);
  double B0 = std::stod(config["Background"]["B0"]);
  double cos_th_cut = std::stod(config["Simulation"]["cos_th_cut"]);

  // ionization state
  Vector1d ion_state;
  double q_avg = 0.0;
  double qsq_avg = 0.0;
  bool neutral = false;
  ion_state.push_back(std::stod(config["Ionization"]["I"]));
  ion_state.push_back(std::stod(config["Ionization"]["II"]));
  ion_state.push_back(std::stod(config["Ionization"]["III"]));
  ion_state.push_back(std::stod(config["Ionization"]["IV"]));
  normIonState(ion_state, q_avg, qsq_avg);
  if ( q_avg == 0.0 ) neutral = true;

  // physics
  bool do_cerenkov = config["Physics"]["cerenkov"] == "true";
  bool do_sync = config["Physics"]["sync"] == "true";

  // compute useful quantities
  double n_i, n_e_free, lam_deb;
  calcLamDeb(ab, rho, temp, q_avg, qsq_avg, n_i, n_e_free, lam_deb);
  
  // clear files and write info file
  std::string infofile = config["IO"]["outpath"] + "/info.txt";
  std::string outfile = config["IO"]["outpath"] + "/data.bin." + std::to_string(rank);
  if ( !cont ) {
    clearOutfile(outfile);
    if ( rank == 0 ) { 
      clearInfo(infofile); 
      writeInfo(infofile, size, config, ab, n_i, n_e_free, lam_deb, mmw, q_avg, qsq_avg);
    }
  }

  // initialize the simulation
  int id;
  double cos_alpha;
  Part part = Part(0, constants::m_e, constants::e, ener, 1.0);
  Sim sim = Sim(
    part, 
    eedl, 
    ab, 
    outfile, 
    rho, 
    temp, 
    ion_state,
    B0, 
    cos_th_cut, 
    neutral, 
    do_cerenkov, 
    do_sync,
    mmw
  );
  int count_loc = 0, count_dead_loc = 0;
  if ( rank == 0 ) {
    std::cout << "Starting simulation." << std::endl << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "rank " << rank << std::flush;
  
  while ( true ) {
    // start a timer for the particle
    auto start_part = std::chrono::steady_clock::now();
    
    // reset the simulation with a new particle
    id = size * count_loc + rank;
    cos_alpha = 2.0 * xi() - 1.0;
    Part part = Part(id, constants::m_e, constants::e, ener, cos_alpha);
    sim.reset(part);

    // run the simulation
    while ( sim.part.alive ) {
      sim.step();
      std::cout << "\rrank " << rank << ": " << count_loc << " (" << sim.part.ener << "eV)" << std::flush;
      if ( sim.part.ener < ener_min ) { 
        sim.kill(); count_dead_loc++; 
      }
      if ( tmax > 0. && sim.time > tmax ) break;
    }
    writeEvent(sim.outfile, sim.event_list);
    sim.event_list.clear();

    // check if the simulation is over
    count_loc++;
    auto now = std::chrono::steady_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    auto tpart = std::chrono::duration_cast<std::chrono::seconds>(now - start_part).count();  
    if ( runtime >= (tsim - 2.0 * tpart) ) break;
    if ( count_max > 0 && count_loc >= count_max ) break;
  }
  std::cout << std::endl << std::endl;

  // compute the packet counts
  int count_glob, count_dead_glob;
  std::vector<int> count_loc_list(size);
  MPI_Gather(&count_loc, 1, MPI_INT, count_loc_list.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_loc, &count_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_dead_loc, &count_dead_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if ( rank == 0 ) {
    auto now = std::chrono::steady_clock::now();
    auto tsim = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    std::cout << "Simulation complete." << std::endl;
    std::cout << "Total packet count: " << count_glob << std::endl;
    std::cout << "Total dead count:   " << count_dead_glob << std::endl;
    std::cout << "Runtime [s]:        " << tsim << std::endl;
    std::cout << std::endl;
    std::cout << "Local packet counts" << std::endl;
    std::cout << "rank,count" << std::endl;
    for ( size_t i = 0; i < count_loc_list.size(); i++ ) {
      std::cout << i << "," << count_loc_list[i] << std::endl;
    }
  }

  MPI_Finalize();
  return 0;
}
