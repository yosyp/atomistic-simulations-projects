/*!
 * @file IO.cpp
 * @author Maxim Shugaev
 * @date 22 Jan 2018
 * @brief The file contains implementation of methods of IO classes declared in IO.h */

#include "IO.h"
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace MSE6270_MD;

// Reader_Data==================================================================================================

Reader_Data::Reader_Data(const std::string& filename_) : filename(filename_) {}

int Reader_Data::load(System& system) {  // load system form a restarting (.data) file
  std::ifstream file(filename);          // open input file
  if (!file.is_open()) {
    std::cout << "Cannot open " << filename << std::endl;
    return EXIT_FAILURE;
  };
  std::cout << "Reading " << filename << std::endl;
  int n;
  file >> n >> system.time;                            // read the number of particles and time
  for (int i = 0; i < 3; i++) file >> system.size[i];  // read system dimensions and the position of the system center
  for (int i = 0; i < 3; i++) file >> system.center[i];

  auto& particles = system.particles;
  particles.resize(n, Particle());  // resize array to match the number of particles
                                    // read type, khist (group), position, and velocities of particles
  for (int i = 0; i < n; i++) file >> particles[i].type;
  for (int i = 0; i < n; i++)
    file >> particles[i].khist >> particles[i].pos.x >> particles[i].pos.y >> particles[i].pos.z;
  for (int i = 0; i < n; i++) file >> particles[i].vel.x >> particles[i].vel.y >> particles[i].vel.z;

  return EXIT_SUCCESS;
}

// Reader_CG==================================================================================================
// initialize Reader_CG, *"ncells" is the number of unit cells{ nx, ny, nz }, "a_lat" is the lattice constant, "lattice"
// is the crystal structure of the generated system
Reader_CG::Reader_CG(const std::array<int, 3>& ncells_, double a, Lattice lattice) {
  std::copy(ncells_.begin(), ncells_.end(), ncells.begin());  // set value of ncells
  a_lat = a;
  switch (lattice) {  // assign parameters corresponding to a particular lattice
    case Lattice::SC:
      cell_n = 1;
      types = {1};
      basis = {{0.0, 0.0, 0.0}};
      shift.x = 0.5;
      shift.y = 0.5;
      shift.z = 0.5;
      lattice_type = "Simple Cubic";
      break;
    case Lattice::BCC:
      cell_n = 2;
      types = {1, 1};
      basis = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
      shift.x = 0.25;
      shift.y = 0.25;
      shift.z = 0.25;
      lattice_type = "Body Centered Cubic";
      break;
    case Lattice::BCC_CsCl:
      cell_n = 2;
      types = {1, 2};
      basis = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
      shift.x = 0.25;
      shift.y = 0.25;
      shift.z = 0.25;
      lattice_type = "Body Centered Cubic (CsCl)";
      break;
    case Lattice::FCC:
      cell_n = 4;
      types = {1, 1, 1, 1};
      basis = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
      shift.x = 0.25;
      shift.y = 0.25;
      shift.z = 0.25;
      lattice_type = "Face Centered Cubic";
      break;
    case Lattice::FCC_Ni3Al:
      cell_n = 4;
      types = {2, 1, 1, 1};
      basis = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
      shift.x = 0.25;
      shift.y = 0.25;
      shift.z = 0.25;
      lattice_type = "Face Centered Cubic (Ni3Al)";
      break;
    case Lattice::Diamond:
      cell_n = 8;
      types = {1, 1, 1, 1, 1, 1, 1, 1};
      basis = {{0.0, 0.0, 0.0},    {0.0, 0.5, 0.5},    {0.5, 0.0, 0.5},    {0.5, 0.5, 0.0},
               {0.25, 0.25, 0.25}, {0.25, 0.75, 0.75}, {0.75, 0.25, 0.75}, {0.75, 0.75, 0.25}};
      shift.x = 0.125;
      shift.y = 0.125;
      shift.z = 0.125;
      lattice_type = "Diamond lattice";
      break;
    case Lattice::Diamond_GaAs:
      cell_n = 8;
      types = {1, 1, 1, 1, 2, 2, 2, 2};
      basis = {{0.0, 0.0, 0.0},    {0.0, 0.5, 0.5},    {0.5, 0.0, 0.5},    {0.5, 0.5, 0.0},
               {0.25, 0.25, 0.25}, {0.25, 0.75, 0.75}, {0.75, 0.25, 0.75}, {0.75, 0.75, 0.25}};
      shift.x = 0.125;
      shift.y = 0.125;
      shift.z = 0.125;
      lattice_type = "Diamond lattice (GaAs)";
      break;
    default:
      cell_n = 0;
      shift.x = 0.0;
      shift.y = 0.0;
      shift.z = 0.0;
      lattice_type = "Unknown";
  }
}

int Reader_CG::load(System& system) {                  // load system form a Crystal Generator pseudo file
  int n = cell_n * ncells[0] * ncells[1] * ncells[2];  // the total number of particles

  std::cout << "Building " << lattice_type << " crystal with nx = " << ncells[0] << ", ny = " << ncells[1]
            << ", nz = " << ncells[2] << std::endl;
  std::cout << "lattice constant = " << a_lat << ", number of atoms = " << n << std::endl;

  system.time = 0.0;  // set time, system dimensions, and the position of the system center
  for (int i = 0; i < 3; i++) system.size[i] = ncells[i] * a_lat;
  for (int i = 0; i < 3; i++) system.center[i] = 0.5 * ncells[i] * a_lat;

  auto& particles = system.particles;
  particles.resize(n, Particle());  // resize array to match the number of particles

  int id = 0;
  for (int iz = 0; iz < ncells[2]; iz++) {
    for (int iy = 0; iy < ncells[1]; iy++) {
      for (int ix = 0; ix < ncells[0]; ix++) {
        for (int ic = 0; ic < cell_n; ic++) {
          // generate a particle based on the cell indexes and the index of atom within the cell
          particles[id].pos = (basis[ic] + shift +
                               DPVector3(static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz))) *
                              a_lat;
          particles[id].vel = DPVector3();
          particles[id].type = types[ic];
          particles[id].khist = 0;
          id++;
        }
      }
    }
  }
  return EXIT_SUCCESS;
}

// Writer_Data==================================================================================================

Writer_Data::Writer_Data(const std::string& filename_) : filename(filename_) {}

int Writer_Data::save(const System& system) {  // save restarting (.data) file
  std::cout << "Writing " << filename << std::endl;
  std::ofstream file(filename);  // open a file
  auto& particles = system.particles;
  int n = particles.size();  // number of particles

  file << n << " " << system.time << std::endl;               // write number of particles and time
  for (int i = 0; i < 3; i++) file << system.size[i] << " ";  // write system dimensions
  file << std::endl;
  for (int i = 0; i < 3; i++) file << system.center[i] << " ";  // write positions of the center of the system
  file << std::endl;

  for (int i = 0; i < n; i++) {  // write particle types
    file << particles[i].type << " ";
    if (i % 100 == 99) file << std::endl;  // limit the maximum number of numbers in one line
  }
  file << std::endl;
  for (int i = 0; i < n; i++) {  // write khist and positions of particles
    file << particles[i].khist << " " << particles[i].pos.x << " " << particles[i].pos.y << " " << particles[i].pos.z
         << " ";
    if (i % 100 == 99) file << std::endl;
  }
  file << std::endl;
  for (int i = 0; i < n; i++) {  // write velocities of particles
    file << particles[i].vel.x << " " << particles[i].vel.y << " " << particles[i].vel.z << " ";
    if (i % 100 == 99) file << std::endl;
  }
  return EXIT_SUCCESS;
}

// Writer_Snapshot==================================================================================================

Writer_Snapshot::Writer_Snapshot(const std::string& filename_) : filename(filename_) {}

int Writer_Snapshot::save(const System& system) {  // Save atomic snapshot (.d) file
  std::cout << "Writing " << filename << std::endl;
  std::ofstream file(filename);  // open  file
  auto& particles = system.particles;
  int n = particles.size();  // number of particles
  for (int i = 0; i < n;
       i++) {  // write: step, khist, pos[0-2] [A], vel[0-2] [A/ps], Ep [eV], Ek [eV], T [K], Etot [eV], type
    double Etot = particles[i].Ep + particles[i].Ek;
    double T = 2.0 * particles[i].Ek * kb / 3.0;
    file << i << " " << particles[i].khist << " " << particles[i].pos.x << " " << particles[i].pos.y << " "
         << particles[i].pos.z << " " << particles[i].vel.x << " " << particles[i].vel.y << " " << particles[i].vel.z
         << " " << particles[i].Ep << " " << particles[i].Ek << " " << T << " " << Etot << " " << particles[i].type
         << std::endl;
  }
  return EXIT_SUCCESS;
}