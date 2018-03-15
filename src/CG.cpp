/*!
 * @file CG.cpp
 * @author Maxim Shugaev
 * @date 27 Jan 2018
 * @brief The file implements Crystal Generator */

#include <string>
#include "Simulation.h"
using namespace MSE6270_MD;

int main() {
  static const std::string output_file = "Ar.data";  // output file
  static const int nx = 5;                           // number of unit cells in x,y, and z direction
  static const int ny = 5;
  static const int nz = 5;
  static const double a_lat = 5.405;                    // lattice parameter
  static const auto lattice = Reader_CG::Lattice::FCC;  // lattice type; check Reader_CG in IO.cpp for more details

  System system;                               // an object for storing MD system
  Reader_CG cg({nx, ny, nz}, a_lat, lattice);  // crystal generator
  Writer_Data writer(output_file);             // a writer capable for saving .data files
  cg.load(system);  // generate system (in the MD code, you also can generate system instead of loading a file)
  // for (auto& particle : system.particles) {					//defining the rigid layers and "isotopes" for diffusion
  // simulation (Homework #5) 	if (particle.pos.z > system.Zcenter) particle.khist = 1; 	if ((particle.pos.z <
  //(a_lat + 0.01)) || 		(particle.pos.z > (system.XL - a_lat - 0.01))) particle.khist = 3;
  //}
  writer.save(system);  // save system to a file
                        // write .d snapshot to versify if system is generated correctly
  Writer_Snapshot snapshot("system.d");  // a writer capable for saving .d files
  snapshot.save(system);

  return EXIT_SUCCESS;
}