/*!
* @file Simulation.h
* @author Maxim Shugaev
* @date 20 Jan 2018
* @brief The file contains Simulation class combining different parts of MSE6270_MD code to perform MD simulations */

#pragma once
#include "Handlers.h"
#include "System.h"
#include "Integrators.h"
#include "IO.h"
#include "NSearch.h"
#include <vector>
#include <memory>
#include <iostream>

namespace MSE6270_MD {
	/*! This class is the heart piece of MSE6270_MD that ties together different parts of the code */
	class Simulation {
	public:
		using Handler_ptr = std::shared_ptr<Handler_Base>;				//An abbreviated name for a smart pointer to Handler_Base
		using Integrator_ptr = std::shared_ptr<Integrator_Base>;		//An abbreviated name for a smart pointer to Integrator_Base
		using NSearch_ptr = std::shared_ptr<NSearch_Base>;				//An abbreviated name for a smart pointer to NN_Search_Base

		Simulation();													//!< Initialize an object of the Simulation class
		int run(int nsteps = 1);										//!< Perform "nsteps" steps of MD
		int add_handler(const Handler_ptr &handler);					//!< Attach a handler to the simulation, which will be called after each MD step. Check Handlers.h for more details.
		int set_integrator(const Integrator_ptr &integrator);			//!< Set an integrator used in the simulation
		int set_nsearch(const NSearch_ptr &n_search);					//!< Set an nearest neighbour search algorithm
		System& get_system();											//!< Return a reference to a system used in the simulation. Check System.h for more details.
		void print(std::ostream& os) const;								//!< Print information about the simulation setup into "os" stream
	private:		
		Simulation(const Simulation &simulation) = default;				//forbid creation of copies
		void gather();													//!< Put particles back into the computational cell when periodic boundary conditions are used
		void stress();													//!< Calculate Sxx, Syy, and Szz components of stress tensor for each atom 

		std::vector<Handler_ptr> handlers;								/*!< An array of smart pointers to handlers called after each MD step. Check Handlers.h for more details.*/
		Integrator_ptr integrator;										/*!< A smart pointer to an integrator used in the simulation.*/
		NSearch_ptr n_search;											/*!< A smart pointer to a neighbour search object.*/
		System system;													/*!< System used in the simulation. Check System.h for more details.*/
	};
}