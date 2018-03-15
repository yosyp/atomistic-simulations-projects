/*!
* @file main_direct.cpp
* @author Maxim Shugaev
* @date 27 Jan 2018
* @brief An example of code with manual building the simulation setup (using MSE6270_MD as a library), without loading configuration from "md.rc".
*		 The following example generates an identical setup to one obtained from md.input file, provided for Homework #3.
*		 After modifying this file, recompile your code. This method provides more flexibility and may be useful for your projects. */

#include "Simulation.h"
#include "SimulationSetup.h"
#include <cstdlib>
#include <memory>
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <string>
#include <chrono>

using namespace MSE6270_MD;

int main() {
	static const std::string input_file = "Ar.data";			//input file with system configuration (.data)
	static const std::string output_file = "Art.data";			//restarting file
	static const std::string info_file = "Art.out";				//a file with simulation statistics
	static const std::string output_dir = "data";				//a directory for writing simulation snapshots (.d)

	static const double dt = 0.004;								//time step
	static const int nsteps = 5000;								//number of steps in the simulation
	static const double T_init = 300.0;							//initial temperature
	static const double Rskin = 2.0;							//the thickness of skin layer used in neighbour list algorithm
	static const int nblist_frequency = 50;						//the frequency of neighbour list update
	static const int write_frequency = 500;						//the frequency of saving snapshots and restarting files
	static const int print_frequency = 50;						//the frequency of printing simulation statistics


	Simulation simulation;										//create "simulation" object
	auto& system = simulation.get_system();						//an abbreviated name for the simulated system
	system.dt = dt;
	system.boundaries = { Boundary::periodic, Boundary::periodic, Boundary::periodic };

	auto integrator = std::make_shared<Integrator_Nord5>();		//create a smart pointer to Nord5 integrator
	simulation.set_integrator(integrator);						//set integrator used in the simulation

	auto nb_list = std::make_shared<NB_list>();					//create a smart pointer to NB_list object used for nearest neighbour search
	nb_list->set_Rskin(Rskin);									//set the thickness of a skin layer
	nb_list->set_update_frequency(nblist_frequency);			//set the update frequency of the neighbour list
	simulation.set_nsearch(nb_list);							//set a nearest neighbour search algorithm used in the simulation

	auto handler_vel = std::make_shared<Handler_Velocity>(T_init); //create a handler that assigned initial velocities
	simulation.add_handler(handler_vel);						//add handler_vel handler to the simulation

	using Material = Potential_LJ::Material;
	auto potential = std::make_shared<Potential_LJ>(Material::Ar); //create potential for Ar
	system.set_ntypes(1);										//set the maximum number of particle types in the simulation
	system.potentials[0] = potential;							//assign potential
	system.masses[0] = potential->get_mass();					//assign mass or Ar
	//system.potentials[0]->print();								//print (type1-type1) interaction potential; check "print" in "Potentials.cpp" for more details

	Reader_Data loader(input_file);
	if (loader.load(system) != EXIT_SUCCESS) {					//load initial atomic configuration
		std::cout << "Cannot load " << input_file << std::endl;
		return EXIT_FAILURE;
	}																				

	auto handler_E = std::make_shared<Handler_Energy>();		//create an energy handler, which computes Etot, Ep, and Ek
	auto handler_T = std::make_shared<Handler_Temperature>();	//create a temperature handler, which computes T
	auto handler_P = std::make_shared<Handler_Pressure>();		//create a pressure handler, which computes P
	handler_E->set_update_frequency(print_frequency);			//set update frequency
	handler_T->set_update_frequency(print_frequency);
	handler_P->set_update_frequency(print_frequency);
	simulation.add_handler(handler_E);							//add the handlers to the simulation
	simulation.add_handler(handler_T);
	simulation.add_handler(handler_P);

	std::ofstream info(info_file);								//open a file for writing simulation statistics
	info << "Loading the initial atomic configuration from \"" << input_file << "\"" << std::endl;
	info << "Number of steps: " << nsteps << std::endl;
	simulation.print(info);										//write information about the simulation setup
	info << "Step   Time   Energy   Kinetic   Potential   Temperature   Pressure" << std::endl;
	info << std::setprecision(10);								//increase the number of digits written to the file

	std::cout << "Starting the simulation" << std::endl;
	std::cout << "Number steps: " << nsteps << std::endl << std::endl;
	
	using clock = std::chrono::high_resolution_clock;
	auto t_start = clock::now();								//a starting point for the measurement of the runtime

	for (int i = 0; i <= nsteps; i++) {
		simulation.run();										//perform 1 MD step

		if (i % print_frequency == 0) {							//write simulation statistics		
			double time = system.time;
			double Etot = handler_E->get_Etot();
			double Ek = handler_E->get_Ek();
			double Ep = handler_E->get_Ep();
			double T = handler_T->get_T();
			double P = handler_P->get_P();
			info << i << " " << time << " " << Etot << " " << Ek << " " << Ep << " " << T << " " << P << std::endl;
			std::cout << "time = " << time << ", Etot = " << Etot << ", T = " << T << std::endl;
		}

		if (i % write_frequency == 0) {
			Writer_Data writer(output_file);
			writer.save(system);								//write a file for restarting

			int time = static_cast<int>(system.time);
			std::string filename = output_dir + "/time" + step_to_str(time) + ".d";
			Writer_Snapshot snapshot(filename);
			snapshot.save(system);								//write a snapshot
		}
	}
	
	auto t_end = clock::now();									//an ending point for the measurement of the runtime
	std::cout << "Runtime: " << std::chrono::duration<float>(t_end - t_start).count() << " s" << std::endl;

	return EXIT_SUCCESS;
}
