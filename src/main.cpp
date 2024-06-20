// defines
#define _USE_MATH_DEFINES

// includes
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <random>

// header files
#include "parser.h"
#include "main.h"

int main() {
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::cout << dis(gen) << std::endl;

    std::string config_dir = "../config/config.ini";
    Config config = parseConfig(config_dir);
    std::cout << config["Simulation"]["npac"] << std::endl;

    EEDLData eedl_data;
    std::cout << parseEEDL(config["Directories"]["EEDL"], eedl_data) << std::endl;
}