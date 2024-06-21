// includes
#include <iostream>
#include <fstream>
#include <vector>

// headers
#include "io.h"

/// @brief A constructor to initialize the Event structure.
Event::Event(int id_, int nstep_) 
  : id(id_)                // The particle ID.
  , nstep(nstep_)          // The step number.
  , Zelem(-1)              // The proton number of the element.
  , interaction(-1)        // The interaction flag.
  , ion(-1)                // The ion index.
  , time(0.0)              // The event time.
  , x(0.0), y(0.0), z(0.0) // The event coordinates.
  , ener(0.0)              // The particle energy.
  , cos_th(1.0)            // The cosine of the scattering angle.
  , ener_loss(0.0)         // The energy lost.
  , ener_sec(0.0)          // The energy of the secondary.
  , ener_loss_sync(0.0)    // The energy lost due to synchrotron.
 {}

/**
 * @brief Clear the outfile.
 * 
 * @param outfile The outfile name.
*/
void clearOutfile(const std::string& outfile) {
    std::ofstream file(outfile, std::ios::binary | std::ios::trunc);
    file.close();
}

/**
 * @brief Write a list of events to the outfile.
 * 
 * @param outfile    The outfile name.
 * @param event_list A vector of events.
 * 
 * @return 0 if success, 1 if failure.
*/
int writeEvent(std::string outfile, std::vector<Event> event_list) {
  std::ofstream file(outfile, std::ios::out | std::ios::binary | std::ios::app);

  if ( !file ) {
    std::cerr << "Failed to open file for writing." << std::endl;
    return 1;
  }

  file.write(reinterpret_cast<const char*>(event_list.data()), event_list.size() * sizeof(Event));

  if ( !file.good() ) {
    std::cerr << "Error writing to file." << std::endl;
    return 1;
  }
  
  file.close();
  return 0;
}
