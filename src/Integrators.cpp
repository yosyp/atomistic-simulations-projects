/*!
* @file Integrators.cpp
* @author Maxim Shugaev
* @date 24 Jan 2018
* @brief This file contains classes implementation of numerical algorithms declared in Integrators.h */

#include "Integrators.h"
#include <cstdlib>
#include <iostream>

using namespace MSE6270_MD;

//Integrator_Base==================================================================================================

Integrator_Base::Integrator_Base() : system(nullptr) {}

void Integrator_Base::init() {}										//initialize internal variables of the integrator when system is assign
																	//Declare an empty method since it can also be empty in children classes.
																	//So, no multiple overriding of a purely virtual function is needed.

int Integrator_Base::set_system(System *system_) {					//assign a system (MD setup) to the integrator
	if (system_ == nullptr) return EXIT_FAILURE;
	system = system_;
	init();															//initialize internal variables of integrator according to the assigned MD setup
	return EXIT_SUCCESS;
}

void Integrator_Base::pre_step() {}									//perform a pre step of integration used in such algorithms like predictor-corrector
																	//Declare an empty method since it can also be empty in children classes.
																	//So, no multiple overriding of a purely virtual function is needed.

//Integrator_Verlet==================================================================================================

void Integrator_Verlet::step() {									//perform a step of integration
	if (system == nullptr) {										//exit if system is not set
		std::cout << "MD system is not assigned to the integrator" << std::endl;
		return;
	}
	auto &particles = system->particles;							//an abbreviated name
	int n = particles.size();										//number of particles
	double dt = system->dt;											//time step

	for (int i = 0; i < n; i++) {
		if (system->rigid(i)) continue;								//"rigid" particles do not move
		double m = system->masses[particles[i].type - 1];			//particle mass; particle type ids start from 1 wile indexes in array start form 0, therefore "-1"
		DPVector3 a = particles[i].force/(m*enunit);				//acceleration
		if (system->dim == 2) {										//in 2D simulation particles cannot move along z
			a.z = 0.0;
			particles[i].vel.z = 0.0;
		}
		particles[i].vel += 0.5*dt*a;								//velocity at t step
		particles[i].Ek = 0.5*m*enunit*particles[i].vel.length_s(); //kinetic energy
		particles[i].pos += particles[i].vel*dt + 0.5*dt*dt*a;		//position at t+1 step
		particles[i].vel += 0.5*dt*a;								//velocity at t+1/2 step
																	//This implementation allows to avoid using additional temporal variable for storing values on the previous step
																	//and performing a prestep before force calculation. However, the velocity outside the integrator is shifted by a half step.
	}
	system->time += dt;												//increase time by one time step
}

std::string Integrator_Verlet::info() const {						//return a string with information about an integration algorithm
	return "Verlet integration algorithm: time step = " + std::to_string(system->dt) + " ps";
}

//Integrator_Nord5==================================================================================================

//Following constant are taken form p. 154 in C. William Gear, Numerical initial value problems
//in ordinary differential equations, Prentice-Hall Inc., Englewood Cliffs, NJ, 1971 
const double Integrator_Nord5::C[6] = {3.0/20.0, 251.0/360.0, 1.0, 11.0/18.0, 1.0/6.0, 1.0/60.0};

void Integrator_Nord5::init(){
	q.resize(0);													//set size to zero to make sure that after resizing all values will be equal to zero
	q.resize(system->particles.size(), { DPVector3(), DPVector3(), DPVector3(), DPVector3() }); //resize array high order derivatives according to the number of particles
}

void Integrator_Nord5::pre_step(){									//perform a pre step of integration(before force calculation): predictor step
	if (system == nullptr) return;									//exit if system is not set
	if (q.size() != system->particles.size()) init();				//check size of the array with high order derivatives
	auto& particles = system->particles;							//an abbreviated name
	int n = particles.size();										//number of particles
	double dt = system->dt;											//time step

	for (int i = 0; i < n; i++) {
		if (system->rigid(i)) continue;								//"rigid" particles do not move

		if (system->dim == 2) {										//in 2D simulation particles cannot move along z
			for (int j = 0; j < 4; j++) q[i][j].z = 0.0;
			particles[i].vel.z = 0.0;
		}
																	//q[particle][k] is a high order (m-th) derivative multiplied by dt^m and divided by m!
		particles[i].pos += particles[i].vel*dt + q[i][0] + q[i][1] + q[i][2] + q[i][3];	//Update high order derivatives. Since v has units [A/dt], it multiplied by dt.
		particles[i].vel += (2.0*q[i][0] + 3.0*q[i][1] + 4.0*q[i][2] + 5.0*q[i][3]) / dt;	//qi_next = sum(C_ij*qj). The sum goes from i to 4, and C_ij - is a binomial coefficient.
		q[i][0] += 3.0*q[i][1] + 6.0*q[i][2] + 10.0*q[i][3];								//Since position and velocity are set, index of q is shifted by 2.
		q[i][1] += 4.0*q[i][2] + 10.0*q[i][3];						//Check lecture notes for more details on the predictor-corrector algorithm
		q[i][2] += 5.0*q[i][3];
	}
}

void Integrator_Nord5::step() {
	if (system == nullptr) {										//exit if system is not set
		std::cout << "MD system is not assigned to the integrator" << std::endl;
		return;
	}
	if (q.size() != system->particles.size()) init();				//check size of the array with high order derivatives
	auto& particles = system->particles;							//an abbreviated name
	int n = particles.size();										//number of particles
	double dt = system->dt;											//time step

	for (int i = 0; i < n; i++) {
		if (system->rigid(i)) continue;								//"rigid" particles do not move

		double m = system->masses[particles[i].type - 1];			//particle mass; particle type ids start from 1 wile indexes in array start form 0, therefore "-1"
		DPVector3 correction = 0.5*dt*dt/(m*enunit)*particles[i].force - q[i][0]; //correction calculated based on the calculated value of force
		if (system->dim == 2) {										//in 2D simulation particles cannot move along z
			for (int j = 0; j < 4; j++) q[i][j].z = 0.0;
			particles[i].vel.z = 0.0;
			correction.z = 0.0;
		}

		particles[i].Ek = 0.5*m*enunit*particles[i].vel.length_s();	//kinetic energy		
		particles[i].pos += C[0] * correction;						//correction of the predicted values
		particles[i].vel += C[1] * correction / dt;
		for (int j = 0; j < 4; j++) q[i][j] += C[j + 2] * correction;
	}
	system->time += dt;												//increase time by one time step
}

std::string Integrator_Nord5::info() const {						//return a string with information about an integration algorithm
	return "Nord5 (5th order predictor-corrector integration algorithm): time step = " + std::to_string(system->dt) + " ps";
}
