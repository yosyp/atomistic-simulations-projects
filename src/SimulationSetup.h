/*!
* @file SimulationSetup.h
* @author Maxim Shugaev
* @date 20 Jan 2018
* @brief The file contains SimulationSetup class used for compatibility with Fortran version of MD code */

#pragma once
#include "Simulation.h"
#include <string>

/*! A class for loading input files used in Fortran version of MD code
	and building a corresponding simulation setup */
class SimulationSetup {
public:
	SimulationSetup(const std::string& filelist =  "md.rc");	//!< Initialize the setup based on the files listed in "filelist" file
	int set(MSE6270_MD::Simulation& simulation) const;			//!< Build "simulation" based on the loaded information
	bool is_open() const;										//!< Check if the information is loaded
	std::string get_input() const;								//!< Return the name of a file with the initial atomic configuration
	std::string get_output() const;								//!< Return the name of a file used for saving restarting files
	std::string get_output_dir() const;							//!< Return the name of a directory used for saving atomic snapshots (should be created manually)
	std::string get_statistics() const;							//!< Return the name of a file used for writing simulation statistics: Step, Time, Energy, Kinetic, Potential, Temperature, Pressure
	int get_nsteps() const;										//!< Return the maximum number of steps in the simulation
	int get_nprint() const;										//!< Return the frequency of printing simulations statistics
	int get_nwrite() const;										//!< Return the frequency of writing snapshots and restarting files
private:
	void set_names(int id, const std::string & name);			// Set file names according to the convention used in the Fortran version of MD code, id - file id, name - file name
	int load_info();											// Read the simulation setup (md.input)
	
	std::string input_snapshot;									/*!< The name of a file with the initial atomic configuration */
	std::string output_snapshot;								/*!< The name of a file used for saving restarting files */
	std::string output_dir;										/*!< The name of a directory used for saving atomic snapshots (should be created manually) */
	std::string info;											/*!< The name of a file with the simulation setup (md.input) */
	std::string simulation_statistics;							/*!< The name of a file used for writing simulation statistics */				
	bool open_flag;												/*!< A flag showing that the input is loaded correctly */

																//following variables are described in "md.input"
	int nstep, newtab, neprt, nwrite, nper, kflag, lflag, ipcon, keybs, lidx, lidy, lidz, kbound, ndim;
	double qtem, qpress, delta, Rskin;
};

																/*! Convert "val" into a string with adding leading zeros to match the name convention from the Fortran version of MD code */
std::string step_to_str(int val);