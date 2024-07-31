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
  double time = std::stod(config["Background"]["ab_time"]) * constants::day;
  parseAb(config["IO"]["ab"], time, ab);
  if ( rank == 0 ) std::cout << "Read abundance data." << std::endl;

  // read EEDL data
  EEDLData eedl;
  parseEEDL(config["IO"]["EEDL"], eedl);
  if ( rank == 0 ) std::cout << "Read EEDL data." << std::endl;

  // read the config file
  int tmax = std::stoi(config["Simulation"]["tmax"]);
  double ener_min = std::stod(config["Simulation"]["ener_min"]);
  bool cont = config["Simulation"]["continue"] == "true";
  double ener = std::stod(config["Particle"]["ener"]);
  double tsim_part = std::stod(config["Particle"]["tpart"]);
  double rho = std::stod(config["Background"]["rho"]);
  double temp = std::stod(config["Background"]["temp"]);
  double ion_state_avg = std::stod(config["Background"]["ion_state_avg"]);
  double L = std::stod(config["Bfield"]["L"]);
  double beta = std::stod(config["Bfield"]["beta"]);
  double mach_A = std::stod(config["Bfield"]["mach_A"]);
  double sig_turb_frac = std::stod(config["Bfield"]["sig_turb_frac"]);
  double alpha = std::stod(config["Bfield"]["alpha"]);
  double cos_th_cut = std::stod(config["Simulation"]["cos_th_cut"]);

  // compute useful quantities
  double n_i, n_e_free, lam_deb, B0;
  calcLamDeb(ab, rho, temp, ion_state_avg, n_i, n_e_free, lam_deb);
  B0 = calcB0(n_i + n_e_free, temp, beta, mach_A);
  
  // clear files and write info file
  std::string infofile = config["IO"]["outpath"] + "/info.txt";
  std::string outfile = config["IO"]["outpath"] + "/data.bin." + std::to_string(rank);
  if ( !cont ) {
    clearOutfile(outfile);
    if ( rank == 0 ) { 
      clearInfo(infofile); 
      writeInfo(infofile, size, config, ab, eedl, n_i, n_e_free, lam_deb, B0); 
    }
  }

  // initialize the simulation
  Part part = Part(0, constants::m_e, constants::e, ener, B0);
  Sim sim = Sim(part, eedl, ab, outfile, rho, temp, ion_state_avg, L, beta, mach_A, sig_turb_frac, alpha, cos_th_cut);
  int count_loc = 0, count_dead_loc = 0;
  if ( rank == 0 ) {
    std::cout << "Starting simulation." << std::endl << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  while ( true ) {
    // start a timer for the particle
    auto start_part = std::chrono::steady_clock::now();
    
    // reset the simulation with a new particle
    int id = size * count_loc + rank;
    Part part = Part(id, constants::m_e, constants::e, ener, B0);
    sim.reset(part);

    // run the simulation
    while ( sim.part.alive ) {
      sim.step();
      if ( sim.part.ener < ener_min ) { 
        sim.kill(); count_dead_loc++; 
      }
      if ( tsim_part > 0. && sim.time > tsim_part ) break;
    }
    writeEvent(sim.outfile, sim.event_list);
    sim.event_list.clear();

    // check if the simulation is over
    count_loc++;
    auto now = std::chrono::steady_clock::now();
    auto tsim = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    auto tpart = std::chrono::duration_cast<std::chrono::seconds>(now - start_part).count();  
    if ( tsim >= (tmax - 1.5*tpart) ) break;
  }

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
