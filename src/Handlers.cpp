/*!
* @file Handlers.cpp
* @author Maxim Shugaev
* @date 22 Jan 2018
* @brief The file contains implementation of methods of Handler classes declared in Handlers.h */

#include "Handlers.h"
#include <cstdlib>
#include <random>
#include <vector>
#include <iostream>
#include <cmath>

using namespace MSE6270_MD;

//Handler_Base==================================================================================================

Handler_Base::Handler_Base() : active(true), update_n(1) {}

void Handler_Base::update(System & system) {
	if (active && (system.step % update_n == 0)) update_impl(system);							//call update_impl(system) every "update_n" steps if handler is active
}

void Handler_Base::activate() {
	active = true;
}

void Handler_Base::deactivate() {
	active = false;
}

void Handler_Base::set_update_frequency(int n) {
	update_n = n;
}

std::string Handler_Base::info() const {														//return a string with information about a handler							
	return std::string();																		//return an empty string by default
}

//Handler_Velocity==================================================================================================

Handler_Velocity::Handler_Velocity(double T0_) : Handler_Base(), T0(T0_) {}

std::string Handler_Velocity::info() const {													//return a string with information about a handler
	return "Assigning the initial velocity corresponding to T = " + std::to_string(T0) + " K";
}

//this class is responsible for assignment of the velocity to particles according to the Maxwell distribution.
void Handler_Velocity::update_impl(System & system){
	if (system.step == 0) {																		//assign velocities only at the 0th step
		std::cout << info() << std::endl;

		std::mt19937 random_generator;															//create a random generator
		int ntypes = system.get_ntypes();														//number of particle types
		using Distribution = std::normal_distribution<double>;									//an abbreviated name for a normal distribution
		std::vector<Distribution> distributions(ntypes);										//an array of normal distributions
		for (int i = 0; i < ntypes; i++) distributions[i] = Distribution(0.0, sqrt(2.0*kb*T0 / (system.masses[i] * enunit))); //initialize velocity distribution functions
																								//2 is used since half of the energy will go to the potential energy

		auto& particles = system.particles;														//an abbreviated name
		int n = particles.size();																//number of particles
		for (int i = 0; i < n; i++) {															
			auto& distribution = distributions[particles[i].type - 1];							//an abbreviated name; particle type ids start from 1 wile indexes in array start form 0, therefore "-1"
			for (int j = 0; j < DPVector3::dim; j++) particles[i].vel.data[j] = distribution(random_generator); //assign the initial velocity according to the normal distributions
		}
																				
		DPVector3 p_avr = DPVector3();
		for (int i = 0; i < n; i++) p_avr += particles[i].vel*system.masses[particles[i].type - 1]; //the total momentum
		if(n > 0) p_avr /= n;
		for (int i = 0; i < n; i++) particles[i].vel -= p_avr/system.masses[particles[i].type - 1];	//correct the momentum to make sure that the total momentum is zero

		for (int i = 0; i < n; i++) {
			if (system.rigid(i)) particles[i].vel = DPVector3();								//"rigid" particles have zero velocity
			if (system.dim == 2) particles[i].vel.z = 0.0;										//zero velocity in z direction for 2d simulation
		}

		double Ek_avr = 0.0;
		for (int i = 0; i < n; i++) {															
			double m = system.masses[particles[i].type - 1];									//mass; particle type ids start from 1 wile indexes in array start form 0, therefore "-1"
			particles[i].Ek = 0.5*m*enunit*particles[i].vel.length_s();							//kinetic energy
			Ek_avr += particles[i].Ek;
		}
		if (n > 0) Ek_avr /= n;
		double T = Ek_avr / (kb*system.dim);													//half of the energy goes to Ep, therefore "3.0" instead of "1.5"
		double correction = sqrt(T0/T);
		for (int i = 0; i < n; i++) {
			particles[i].vel *= correction;														//correct the kinetic energy to make sure that energy corresponding to T0 is assigned
			particles[i].Ek *= pow(correction,2);
		}
	}
}

//Handler_Energy==================================================================================================

Handler_Energy::Handler_Energy() : Ep(0.0), Ek(0.0), Etot(0.0) {}

void Handler_Energy::update_impl(System & system) {												//calculate Ep, Ek, and Etot
	Ep = 0.0; Ek = 0.0;																			
	for (auto& particle : system.particles) {													//loop over all particles
		Ep += particle.Ep;																		//calculate total Ep and Ek
		Ek += particle.Ek;
	}
	Etot = Ep + Ek;
}

double Handler_Energy::get_Ep() const {
	return Ep;
}

double Handler_Energy::get_Ek() const {
	return Ek;
}

double Handler_Energy::get_Etot() const {
	return Etot;
}

//Handler_Temperature==================================================================================================

Handler_Temperature::Handler_Temperature() : Handler_Base(), T(0.0) {}

void Handler_Temperature::update_impl(System & system) {										//Calculate the Temperature [K], estimated based on the average atomic kinetic energy
	double Ek = 0; 	T = 0.0;
	auto &particles = system.particles;
	int n = particles.size();														
	int n_count = 0;
	for (int i = 0; i < n; i++) {
		if (system.rigid(i) || exclude(particles[i].khist)) continue;							//exclude "rigid" particles and particles from the exclusion list
		Ek += particles[i].Ek;																	//kinetic energy
		n_count++;																				//number of considered atoms
	}
	if(n_count > 0) T = Ek / (0.5*system.dim*kb*n_count);										//<Ek> = 1.5*kb*T
}

bool Handler_Temperature::exclude(int group) const {											//check if the group is in the exclusion list
	return std::find(exclude_list.begin(), exclude_list.end(), group) != exclude_list.end();
}

double Handler_Temperature::get_T() const {
	return T;
}

void Handler_Temperature::exclude_groups(const std::vector<int>& list) {
	exclude_list = list;
}

//Handler_Pressure==================================================================================================

Handler_Pressure::Handler_Pressure() : Handler_Base(), P(0.0) {}

void Handler_Pressure::update_impl(System & system) {											//calculate the Pressure based on the Virial theorem
	P = 0.0; Pcomponents = DPVector3();
	for (auto& particle : system.particles) {													//loop over all particles
		auto& stress = particle.stress;															//stress corresponding to each atom is already calculated in Simulation.cpp with applying Virial theorem
		if (system.dim == 3) P -= (stress.x + stress.y + stress.z) / 3.0;						//the only thing that has to be done is summing everything together
		else P -= 0.5*(stress.x + stress.y);
		Pcomponents -= stress;																	//scale velocities
	}
}

double Handler_Pressure::get_P() const {
	return P;
}

DPVector3 Handler_Pressure::get_Pcomponents() const {
	return Pcomponents;
}

//Handler_BerendsenT==================================================================================================

Handler_BerendsenT::Handler_BerendsenT(double T, double tau_) : Handler_Temperature(), T0(T), tau(tau_) {
	std::cout << info() << std::endl;
}

void Handler_BerendsenT::update_impl(System & system) {											//update particle velocities to maintain T0
	Handler_Temperature::update_impl(system);													//call update_impl(system) from Handler_Temperature to get the value of the current temperature
	double scale = sqrt(1.0 + (T0 / T - 1.0)*system.dt / tau);									//a scaling coefficient used in Berendsen algorithm
	for (auto& particle : system.particles) {													//loop over all particles
		if (exclude(particle.khist)) continue;													//skip if a particle in the exclusion list
		particle.vel *= scale;																	//scale velocities
	}
}

void Handler_BerendsenT::set_T(double T) {
	T0 = T;
}

std::string Handler_BerendsenT::info() const {													//return a string with information about a handler
	return "Berendsen thermostat: T = " + std::to_string(T0) + " K, tau = " + std::to_string(tau) + " ps";
}

//Handler_BerendsenT==================================================================================================

Handler_BerendsenP::Handler_BerendsenP(double P, Mode mode_, double beta_) : Handler_Pressure(), P0(P), beta(beta_), mode(mode_) {
	std::cout << info() << std::endl;
}

void Handler_BerendsenP::update_impl(System & system) {											//update system dimensions to maintain P0
	Handler_Pressure::update_impl(system);														//call update_impl(system) from Handler_Pressure to get the value of the current pressure
	DPVector3 scale = DPVector3();
	double d;

	switch (mode) {																				//select the method of pressure control
	case Mode::XYZ:																				//X, Y, and Z sizes are scaled together
		d = 1.0 - beta * system.dt*(P0 - P);
		scale.x = d; scale.y = d; scale.z = d;
		break;
	case Mode::XY_Zindependent:																	//X and Y sizes are controlled together based on 0.5*(Pxx + Pyy), and Z size is controlled based on Pzz
		d = 1.0 - beta * system.dt*(P0 - 0.5*(Pcomponents.x + Pcomponents.y));
		scale.x = d; scale.y = d;
		scale.z = 1.0 - beta * system.dt*(P0 - Pcomponents.z);
		break;
	case Mode::XYZindependent:																	//X, Y, and Z sizes are scaled independently based on Pxx, Pyy, and Pzz
		for (int i = 0; i < 3; i++) scale.data[i] = 1.0 - beta * system.dt*(P0 - Pcomponents.data[i]);
	case Mode::Z:																				//control only Z size based on Pzz
		scale.z = 1.0 - beta * system.dt*(P0 - Pcomponents.z);
	}

	int n = system.particles.size();															//number of particles
	for (int j = 0; j < 3; j++) {
		if (system.boundaries[j] != Boundary::periodic) continue;								//perform pressure control only if a periodic boundary condition is applied in a particular direction
		system.size[j] *= scale.data[j];														//scale system dimensions
		for (int i = 0; i < n; i++) {
			auto& pos = system.particles[i].pos;												
			pos.data[j] = (pos.data[j] - system.center[j])*scale.data[j] + system.center[j];	//scale particle positions
		}
	}
}

std::string Handler_BerendsenP::to_string(Mode mode) const {
	std::string result;
	switch (mode) {																				//select the method of pressure control
	case Mode::XYZ:
		result = "XYZ together";
		break;
	case Mode::XY_Zindependent:
		result = "XY together, Z independent";
		break;
	case Mode::XYZindependent:
		result = "XYZ independent";
		break;
	case Mode::Z:
		result = "only Z direction";
		break;
	}
	return result;
}

void Handler_BerendsenP::set_P(double P) {
	P0 = P;
}

std::string Handler_BerendsenP::info() const {
	return "Berendsen barostat: P = " + std::to_string(1e-9*P0) + " GPa; pressure control " + to_string(mode) +
		"; beta = " + std::to_string(beta*1e9) + " GPa^-1";
}

//Handler_Heating==================================================================================================

Handler_Heating::Handler_Heating(double T, double tau_) : Handler_Temperature(), T0(T), tau(tau_), done(false) {
	std::cout << info() << std::endl;
}

void Handler_Heating::update_impl(System & system) {											//perform heating of the system
	if (done) return;																			//do nothing if the system is already heated to T0
	Handler_Temperature::update_impl(system);													//call update_impl(system) from Handler_Temperature to get the value of the current temperature
	if (T > T0) {																				//check if the system is heated to T0
		done = true;																			
		std::cout << "The system is heated to " << T0 << " K" << std::endl;
	}
	else {
		double scale = sqrt(1.0 + (1.0 - T / T0)*system.dt / tau);								//scaling coefficient to perform slow heating
		for (auto& particle : system.particles) particle.vel *= scale;							//velocity scaling
	}
}

std::string Handler_Heating::info() const {														//return a string with information about a handler
	return "Slow heating to T = " + std::to_string(T0) + " K, tau = " + std::to_string(tau) + " ps";
}

//Handler_Quench==================================================================================================

Handler_Quench::Handler_Quench() {
	std::cout << info() << std::endl;
}

std::string MSE6270_MD::Handler_Quench::info() const {
	return "Quench algorithm";
}

void Handler_Quench::update_impl(System & system) {												//quench atomic motion
	auto& particles = system.particles;															//an abbreviated name
	int n = particles.size();																	//number of particles
	if (Ek_prev.size() != particles.size()) Ek_prev.resize(n);									//On the first iteration, set Ek_prev equal to the current Ek
	else {																						//Else set particle velocity to zero if Ek decreased in comparison with the previous step
		for (int i = 0; i < n; i++) {
			if (particles[i].Ek < Ek_prev[i]) {
				particles[i].vel = DPVector3();
				particles[i].Ek = 0.0;
			}
		}
	}
	for (int i = 0; i < n; i++) Ek_prev[i] = particles[i].Ek;									//Update Ek_prev
}

