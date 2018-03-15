/*!
* @file main.cpp
* @author Maxim Shugaev
* @date 24 Jan 2018
* @brief main file */

#include "Simulation.h"
#include "SimulationSetup.h"
#include <cstdlib>
#include <memory>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace MSE6270_MD;

int main() {
	SimulationSetup setup("md.rc");								//initialize "setup" objects by data from "md.rc"
	if (!setup.is_open()) {										//check if the input information is loaded
		std::cout << "Cannot load input files" << std::endl;
		return EXIT_FAILURE;
	}
	Simulation simulation;										//create "simulation" object
	auto& system = simulation.get_system();						//an abbreviated name for the simulated system
	setup.set(simulation);										//setup the simulation based on the parameters loaded by "setup"
	//system.potentials[0]->print();							//print (type1-type1) interaction potential; check "print" in "Potentials.cpp" for more details

	Reader_Data loader(setup.get_input());
	if (loader.load(system) != EXIT_SUCCESS) {					//load initial atomic configuration
		std::cout << "Cannot load " << setup.get_input() << std::endl;
		return EXIT_FAILURE;
	}

	auto handler_E = std::make_shared<Handler_Energy>();		//create an energy handler, which computes Etot, Ep, and Ek
	auto handler_T = std::make_shared<Handler_Temperature>();	//create a temperature handler, which computes T
	auto handler_P = std::make_shared<Handler_Pressure>();		//create a pressure handler, which computes P
	handler_E->set_update_frequency(setup.get_nprint());		//set update frequency
	handler_T->set_update_frequency(setup.get_nprint());
	handler_P->set_update_frequency(setup.get_nprint());
	simulation.add_handler(handler_E);							//add the handlers to the simulation
	simulation.add_handler(handler_T);
	simulation.add_handler(handler_P);

	std::ofstream info(setup.get_statistics());					//open a file for writing simulation statistics
	info << "Loading the initial atomic configuration from \"" << setup.get_input() << "\"" << std::endl;
	info << "Number of steps: " << setup.get_nsteps() << std::endl;
	simulation.print(info);										//write information about the simulation setup
	info << "Step   Time   Energy   Kinetic   Potential   Temperature   Pressure" << std::endl;
	info << std::setprecision(10);								//increase the number of digits written to the file

	std::cout << "Starting the simulation" << std::endl;
	std::cout << "Number steps: " << setup.get_nsteps() << std::endl << std::endl;

	using clock = std::chrono::high_resolution_clock;
	auto t_start = clock::now();								//a starting point for the measurement of the runtime
	
	for (int i = 0; i <= setup.get_nsteps(); i++) {
		simulation.run();										//perform 1 MD step

		if (i % setup.get_nprint() == 0) {						//write simulation statistics		
			double time = system.time;
			double Etot = handler_E->get_Etot();
			double Ek = handler_E->get_Ek();
			double Ep = handler_E->get_Ep();
			double T = handler_T->get_T();
			double P = handler_P->get_P();
			info << i << " " << time << " " << Etot << " " << Ek << " " << Ep << " " << T << " " << P << std::endl;
			std::cout << "time = " << time << ", Etot = " << Etot << ", T = " << T << std::endl;
		}

		if (i % setup.get_nwrite() == 0) {
			Writer_Data writer(setup.get_output());
			writer.save(system);								//write a file for restarting

			int time = static_cast<int>(system.time);
			std::string filename = setup.get_output_dir() + "/time" + step_to_str(time) + ".d";
			Writer_Snapshot snapshot(filename);
			snapshot.save(system);								//write a snapshot
		}
	}
	
	auto t_end = clock::now();									//an ending point for the measurement of the runtime
	std::cout << "Runtime: " << std::chrono::duration<float>(t_end - t_start).count() << " s" << std::endl;
	
	return EXIT_SUCCESS;
}
