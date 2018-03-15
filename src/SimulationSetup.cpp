/*!
 * @file SimulationSetup.cpp
 * @author Maxim Shugaev
 * @date 22 Jan 2018
 * @brief The file contains implementation of methods of SimulationSetup class declared in SimulationSetup.h */

#include "SimulationSetup.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
//#include <experimental/filesystem>

using namespace MSE6270_MD;
// initializing the setup based on the files listed in "filelist" file
SimulationSetup::SimulationSetup(const std::string& filelist) : open_flag(false) {
  std::cout << "Loading a list of files from " << filelist << std::endl;
  std::ifstream in(filelist);
  if (!in.is_open()) {
    std::cout << "Cannot open " << filelist << std::endl;
    return;
  }
  while (1) {  // load a list of files
    std::string tmp_str, filename;
    int file_id;
    in >> tmp_str >> file_id >> filename >> tmp_str;
    if (in.eof()) break;
    set_names(file_id, filename);  // set file names according to the id
  }
  // std::experimental::filesystem::create_directory(output_dir);//create output directory (atm supported only by VS)
  if (load_info() != EXIT_SUCCESS) return;  // load information about the simulation (md.input)

  open_flag = true;  // information is loaded correctly
}

Boundary get_boundary(int val) {  // convert integer value into boundary type
  return (val == 0) ? Boundary::free : Boundary::periodic;
}

int SimulationSetup::set(Simulation& simulation) const {  // build "simulation" based on the loaded information
  auto& system = simulation.get_system();
  system.dt = delta;                                                                 // set time step
  system.dim = ndim;                                                                 // set the number of dimensions
  system.boundaries = {get_boundary(lidx), get_boundary(lidy), get_boundary(lidz)};  // set boundary conditions

  auto integrator = std::make_shared<Integrator_Nord5>();  // create a smart pointer to Nord5 integrator
  simulation.set_integrator(integrator);                   // set integrator used in the simulation

  auto nb_list =
      std::make_shared<NB_list>();        // create a smart pointer to NB_list object used for nearest neighbour search
  nb_list->set_Rskin(Rskin);              // set the thickness of a skin layer
  nb_list->set_update_frequency(newtab);  // set the update frequency of the neighbour list
  simulation.set_nsearch(nb_list);        // set a nearest neighbour search algorithm used in the simulation

  switch (kflag) {  // add handlers to the simulation to reproduce behaviour corresponding to a particular value of
                    // kflag
    case 1: {  // Quench
      auto handler = std::make_shared<Handler_Quench>();
      simulation.add_handler(handler);
      break;
    }
    case 2: {  // Velocity distribution
      auto handler = std::make_shared<Handler_Velocity>(qtem);
      simulation.add_handler(handler);
      break;
    }
    case 3: {  // Heating
      auto handler = std::make_shared<Handler_Heating>(qtem);
      simulation.add_handler(handler);
      break;
    }
    case 5: {  // Berendsen thermostat
      auto handler = std::make_shared<Handler_BerendsenT>(qtem);
      simulation.add_handler(handler);
      break;
    }
  }

  if (lflag == 1) {  // Berendsen barostat
    using Mode = Handler_BerendsenP::Mode;
    Mode mode = Mode::XYZ;  // All components together
    switch (ipcon) {        // select a way how pressure is controlled
      case 2:               // Together along XY, and Z independent
        mode = Mode::XY_Zindependent;
        break;
      case 3:  // Independent along XYZ
        mode = Mode::XYZindependent;
        break;
      case 4:  // Pressure control only along Z
        mode = Mode::Z;
        break;
    }
    auto handler = std::make_shared<Handler_BerendsenP>(1e9 * qpress, mode);
    simulation.add_handler(handler);
  }

  if (kbound == 1) {  // Rigid boundary
    auto& system = simulation.get_system();
    system.groups_rigid.push_back(3);  // Atoms with khist = 3 are considered as "rigid"
  }

  switch (keybs) {  // Select interatomic potential
    case 1: {       // Stillinger-Weber potential
      system.set_ntypes(2);
      using Material = Potential_SW::Material;
      auto potential_Si_Si = std::make_shared<Potential_SW>(Material::Si, Material::Si);  // create potentials
      auto potential_Ge_Ge = std::make_shared<Potential_SW>(Material::Ge, Material::Ge);
      auto potential_Si_Ge = std::make_shared<Potential_SW>(Material::Si, Material::Ge);

      system.potentials[0] = potential_Si_Si;          // type1 - type1
      system.potentials[1] = potential_Ge_Ge;          // type2 - type2
      system.potentials[2] = potential_Si_Ge;          // type1 - type2
      system.masses[0] = potential_Si_Si->get_mass();  // assign masses
      system.masses[1] = potential_Ge_Ge->get_mass();
      break;
    }
    case 2: {  // EAM potential
      std::string path = "potentials/CuAg.eam.alloy";
      system.set_ntypes(2);  // An example of CuAg alloy with using Mishin potential (DOI: 10.1088/0965-0393/14/5/002)
      auto potential_Cu_Cu = std::make_shared<Potential_EAM_alloy>(path, "Cu", "Cu");  // create potentials
      auto potential_Ag_Ag = std::make_shared<Potential_EAM_alloy>(path, "Ag", "Ag");
      auto potential_Cu_Ag = std::make_shared<Potential_EAM_alloy>(path, "Cu", "Ag");

      // assign potentials, you may use system.get_interaction_id(i, j) to get id in the array
      system.potentials[0] = potential_Cu_Cu;          // type1 - type1
      system.potentials[1] = potential_Ag_Ag;          // type2 - type2
      system.potentials[2] = potential_Cu_Ag;          // type1 - type2
      system.masses[0] = potential_Cu_Cu->get_mass();  // assign masses
      system.masses[1] = potential_Ag_Ag->get_mass();
      break;
    }
    default: {  // Lennard-Jones Ar potential
      using Material = Potential_LJ::Material;
      auto potential = std::make_shared<Potential_LJ>(Material::Ar);  // create potential for Ar
      system.set_ntypes(1);                      // set the maximum number of particle types in the simulation
      system.potentials[0] = potential;          // assign potential
      system.masses[0] = potential->get_mass();  // assign mass or Ar

      break;
    }
  }

  return EXIT_SUCCESS;
}

bool SimulationSetup::is_open() const { return open_flag; }

std::string SimulationSetup::get_input() const { return input_snapshot; }

std::string SimulationSetup::get_output() const { return output_snapshot; }

std::string SimulationSetup::get_output_dir() const { return output_dir; }

std::string SimulationSetup::get_statistics() const { return simulation_statistics; }

int SimulationSetup::get_nsteps() const { return nstep; }

int SimulationSetup::get_nprint() const { return neprt; }

int SimulationSetup::get_nwrite() const { return nwrite; }
// Set file names according to the convention used in the Fortran version of MD code, id - file id, name - file name
void SimulationSetup::set_names(int id, const std::string& name) {
  switch (id) {
    case 14:
      info = name;
      break;
    case 15:
      input_snapshot = name;
      break;
    case 16:
      output_snapshot = name;
      break;
    case 17:
      simulation_statistics = name;
      break;
    case 1:
      output_dir = name;
      break;
  }
}

template <typename T>  // read a string from "in" stream and return the first number in the string as "val"
void read_line(std::istream& in,
               T& val) {  // template is used to allow the compiler to generate a code for both int and double "val"
  std::string buf;
  std::getline(in, buf);      // read line
  std::stringstream ss(buf);  // create a string stream based on the read line
  ss >> val;                  // get the first value from the string stream
};

int SimulationSetup::load_info() {  // read the simulation setup (md.input)
  std::ifstream in(info);
  if (!in.is_open()) {
    std::cout << "Cannot open " << info << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Loading a simulation setup from " << info << std::endl;
  read_line(in, nstep);
  read_line(in, newtab);
  read_line(in, neprt);
  read_line(in, nwrite);
  read_line(in, nper);
  read_line(in, kflag);
  read_line(in, lflag);
  read_line(in, ipcon);
  read_line(in, keybs);
  read_line(in, lidx);
  read_line(in, lidy);
  read_line(in, lidz);
  read_line(in, kbound);
  read_line(in, ndim);
  read_line(in, qtem);
  read_line(in, qpress);
  read_line(in, delta);
  read_line(in, Rskin);

  return EXIT_SUCCESS;
}

std::string step_to_str(int val) {  // Convert "val" into a string with adding leading zeros to match the name
                                    // convention from the Fortran version of MD code
  int length = 0;
  if (val < 1000)
    length = 3;  // select the length of the returned string
  else if (val < 10000)
    length = 4;
  else if (val < 100000)
    length = 5;
  else
    return std::string();

  std::stringstream ss;
  ss << std::setw(length) << std::setfill('0')
     << val;  // convert "val" into string of length "length" with adding leading zeroes if needed
  return ss.str();
};
