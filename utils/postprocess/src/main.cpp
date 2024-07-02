// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

// headers
#include "io.h"
#include "functions.h"
#include "parser.h"
#include "const.h"

int main(int argc, char** argv) {

  // initialize MPI
  // MPI_Init(&argc, &argv);

  // get MPI rank and size
  // int rank, size;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &size);

  // read command line arguments
  if ( argc < 2 ) {
    // if ( rank == 0 ) 
    std::cerr << "Usage: " << argv[0] << " <data_path>" << std::endl;
    // MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  std::string data_path = std::string(argv[1]);

  // read config file
  std::string config_dir = "../config/config.ini";
  Config config;
  parseConfig(config_dir, config);
  std::cout << "Read config file." << std::endl;

  std::string filename = data_path + "/data.bin.0";

  // make energy list
  std::vector<double> ener_list;
  double ener_min = std::stod(config["Energy"]["min"]);
  double ener_max = std::stod(config["Energy"]["max"]);
  size_t num_ener = std::stoul(config["Energy"]["num"]);
  bool log_ener = config["Energy"]["log"] == "true";
  linspace(ener_min, ener_max, num_ener, log_ener, ener_list);

  // make time list
  std::vector<double> time_list;
  double time_min = std::stod(config["Time"]["min"]) * constants::hr;
  double time_max = std::stod(config["Time"]["max"]) * constants::hr;
  size_t num_time = std::stoul(config["Time"]["num"]);
  bool log_time = config["Time"]["log"] == "true";
  linspace(time_min, time_max, num_time, log_time, time_list);

  // make distance list
  std::vector<double> dis_list;
  double dis_min = std::stod(config["Distance"]["min"]) * constants::AU;
  double dis_max = std::stod(config["Distance"]["max"]) * constants::AU;
  size_t num_dis = std::stoul(config["Distance"]["num"]);
  bool log_dis = config["Distance"]["log"] == "true";
  linspace(dis_min, dis_max, num_dis, log_dis, dis_list);

  std::string infofile = config["IO"]["outpath"] + "/info.txt";
  std::string outfile = config["IO"]["outpath"] + "/data.txt";
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);
  writeInfo(infofile, ener_list, time_list, dis_list);
  clearFile(outfile);
  processFile(filename, outfile, num_event_per_chunk, ener_list, time_list, dis_list);

  // MPI_Finalize();
  return 0;
}
