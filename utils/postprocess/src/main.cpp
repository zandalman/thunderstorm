// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <mpi.h>

// headers
#include "io.h"
#include "functions.h"
#include "parser.h"
#include "const.h"

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
  Config config;
  std::string configfile_name = "../config/config.ini";
  parseConfig(configfile_name, config);
  if ( rank == 0 ) std::cout << "Read config file." << std::endl;
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);

  // get statistics list
  size_t num_stat = 0;
  bool do_stat;
  std::vector<bool> do_stat_list, do_stat_cat;
  std::vector<std::string> stat_list = {
    "ener_loss_mech",
    "num_ion_elem",
    "num_ener_sec",
    "dis_par_time",
    "dis_perp_time",
    "ener_time",
    "ener_loss_time",
    "ener_loss_dis2d"
  };
  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    do_stat = config["Statistics"][stat_list[i]] == "true";
    do_stat_list.push_back(do_stat);
    if ( do_stat ) num_stat += 1;
  }
  do_stat_cat.push_back(do_stat_list[2]);
  do_stat_cat.push_back(do_stat_list[3] || do_stat_list[4] || do_stat_list[5] || do_stat_list[6]);
  do_stat_cat.push_back(do_stat_list[7]);
  do_stat_cat.push_back(do_stat_list[7]);

  // make lists
  size_t num_ener, num_ener_sec, num_time, num_dis_par, num_dis_perp;
  std::vector<double> ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list;
  makeList(config["Energy"], 1., num_ener, ener_list);
  makeList(config["EnergySec"], 1., num_ener_sec, ener_sec_list);
  makeList(config["Time"], constants::hr, num_time, time_list);
  makeList(config["DistancePar"], constants::AU, num_dis_par, dis_par_list);
  makeList(config["DistancePerp"], constants::AU, num_dis_perp, dis_perp_list);

  // write info file
  if ( rank == 0 ) {
    std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, config, ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list, num_stat, do_stat_list, do_stat_cat);
  }

  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(num_file / size));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file);

  MPI_Barrier(MPI_COMM_WORLD);
  int count_loc = 0;
  for ( int i = idx_file_min; i < idx_file_max; i++ ) {
    std::string datafile_name = data_path + "/data.bin." + std::to_string(i);
    std::string outfile_name = config["IO"]["outpath"] + "/data.txt." + std::to_string(i);
    clearFile(outfile_name);
    count_loc += processFile(datafile_name, outfile_name, num_event_per_chunk, ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list, do_stat_list, do_stat_cat);
  }

  // compute the particle number
  int count_glob;
  std::vector<int> count_loc_list(size);
  MPI_Gather(&count_loc, 1, MPI_INT, count_loc_list.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_loc, &count_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if ( rank == 0 ) {
    concatenateFiles(config["IO"]["outpath"], num_file);
    auto now = std::chrono::steady_clock::now();
    auto tpp = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    std::cout << "Post processing complete." << std::endl;
    std::cout << "Packet count: " << count_glob << std::endl;
    std::cout << "Runtime [s]:  " << tpp << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
