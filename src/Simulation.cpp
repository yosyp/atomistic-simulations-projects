/*!
* @file Simulation.cpp
* @author Maxim Shugaev
* @date 22 Jan 2018
* @brief The file contains implementation of methods of Simulation class declared in Simulation.h */

#include "Simulation.h"
#include "Force.h"
#include <cstdlib>
#include <algorithm>
#include <iostream>

//Simulation==================================================================================================

using namespace MSE6270_MD;

Simulation::Simulation() : system(), integrator(nullptr), n_search(nullptr), handlers() {}

int Simulation::run(int nsteps) {												//perform "nsteps" steps of MD
	if ((n_search == nullptr) || (integrator == nullptr)) return EXIT_FAILURE;	//check if the simulation is set up correctly
	if(!system.verify()) return EXIT_FAILURE;									//check if the system is set up correctly
	for (int i = 0; i < nsteps; i++) {
		n_search->update();														//nearest neighbour search
		integrator->pre_step();													//pre step of integrator, used in such algorithms like predictor-corrector
		
		for (auto& particle : system.particles) particle.reset();				//set force, stress, Ep, Ek, ro equal to zero before the next force calculation
		Force::eam_density(system, n_search->get_neighbours());					//calculate the electron density of EAM potentials
		Force::force(system, n_search->get_neighbours());						//calculate forces and potential energies
		Force::force_3b(system, n_search->get_neighbours());					//calculate forces and potential energies for 3 body interaction part if one is present

		stress();																//calculate Sxx, Syy, and Szz components of stress tensor for each atom 
		integrator->step();														//integration of equations of motion
		gather();																//put particles back into the computational cell
		for (auto& handler : handlers) handler->update(system);					//update attached handlers
		system.step++;
	}
	return EXIT_SUCCESS;
}

int Simulation::add_handler(const Simulation::Handler_ptr& handler) {			//attach a handler to the simulation, which will be called after each MD step
	if (handler == nullptr) return EXIT_FAILURE;								//exit if null pointer is passed
	handlers.push_back(handler);												//add the handler to the array
	return EXIT_SUCCESS;
}

int Simulation::set_integrator(const Simulation::Integrator_ptr& integrator_) {	//set an integrator used in the simulation
	if (integrator_ == nullptr) return EXIT_FAILURE;							//exit if null pointer is passed
	integrator = integrator_;
	integrator->set_system(&system);
	return EXIT_SUCCESS;
}

int Simulation::set_nsearch(const Simulation::NSearch_ptr& n_search_) {			//set an nearest neighbour search algorithm
	if (n_search_ == nullptr) return EXIT_FAILURE;								//exit if null pointer is passed
	n_search = n_search_;
	n_search->set_system(&system);
	return EXIT_SUCCESS;
}

System & Simulation::get_system() {
	return system;
}

void Simulation::print(std::ostream & os) const {								//print information about the simulation setup into os stream
	os << system;
	if (system.get_ntypes() > 0) {
		os << "Atomic masses: m[type1] = " << system.masses[0] << " Da";		//print atomic masses
		for (int i = 1; i < system.get_ntypes(); i++)
			os << ", m[type" << i + 1 << "] = " << system.masses[i] << " Da";
		os << std::endl;
	}
	if(n_search != nullptr) os << n_search->info() << std::endl;				//print information about the neighbour search algorithm
	else os << "neighbour search algorithm is not set" << std::endl;
	if (integrator != nullptr) os << integrator->info() << std::endl;			//print information about the integration algorithm
	else os << "Integration algorithm is not set" << std::endl;
	if (!handlers.empty()) {
		os << "Handlers:" << std::endl;											//print information about handlers attached to the simulation
		for (auto& handler : handlers) {
			auto info = handler->info();
			if(info.size() > 0) os << info << std::endl;
		}
	}
	os << std::endl;

	os << "Potentials:" << std::endl;											//print information about interatomic interactions
	int ntypes = system.get_ntypes();
	for (int i = 1; i <= ntypes; i++) {
		for (int j = i; j <= ntypes; j++) {
			int id = system.get_interaction_id(i, j);
			os << "type" << i << "-type" << j << " : ";
			if (system.potentials[id] != nullptr) os << system.potentials[id]->info();
			else os << "null";
			os << std::endl;
		}
	}
	os << std::endl;
}

void Simulation::gather(){														//put particles back into the computational cell when periodic boundary conditions are used
	auto& particles = system.particles;											//an abbreviated name
	int n = particles.size();													//number of particles
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < DPVector3::dim; j++) {								//dim = 3			
			if (system.boundaries[j] != Boundary::periodic) continue;
			auto& pos = particles[i].pos;
			if (pos.data[j] < 0.0) pos.data[j] += system.size[j];				//put particles back into the computational cell
			else if(pos.data[j] > system.size[j]) pos.data[j] -= system.size[j];
		}
	}
}

void Simulation::stress() {														//calculate Sxx, Syy, and Szz components of stress tensor for each atom
	double volume;
	if (system.dim == 2) volume = system.XL*system.YL*1e-20;
	else volume = system.XL*system.YL*system.ZL*1e-30;							//volume of the computational domain
	for (auto& particle : system.particles) {									//loop over all particles
		double m = system.masses[particle.type - 1];							//mass; particle type ids start from 1 wile indexes in array start form 0, therefore "-1"
		auto& velocity = particle.vel;											//velocity
																				//calculation of stress based on the Virial theorem
																				//s_xx = sum_i(-m_i*V_i^2 + 0.5*sum_j((x_j - x_i)*F_ij))/V
																				//before call of this subroutine particle.stress contains dr*F
		particle.stress = (-m * enunit*(velocity | velocity) + 0.5*particle.stress)*(e / volume);	//[J/m3] = [Pa]
	}
}
