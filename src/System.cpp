/*!
* @file System.cpp
* @author Maxim Shugaev
* @date 22 Jan 2018
* @brief The file contains implementation of methods of classes declared in System.h: DPVector3, Particle, System */

#include "System.h"
#include "Potentials.h"
#include <algorithm>
#include <iostream>

using namespace MSE6270_MD;

//DPVector3==================================================================================================

DPVector3::DPVector3() : x(0.0), y(0.0), z(0.0) {}

DPVector3::DPVector3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

//DPVector3==================================================================================================

Particle::Particle() : pos(), vel(), force(), type(0), khist(0), Ep(0.0), Ek(0.0), ro(0.0) {}

void MSE6270_MD::Particle::reset() {								//set force, stress, Ep, Ek, ro equal to zero
	force = DPVector3();
	stress = DPVector3();
	Ep = 0.0;
	Ek = 0.0;
	ro = 0.0;
}

//System==================================================================================================

System::System() : XL(0.0), YL(0.0), ZL(0.0), Xcenter(0.0), Ycenter(0.0), Zcenter(0.0), //set default values
	boundaries{ Boundary::periodic, Boundary::periodic, Boundary::periodic }, particles(),
	masses(), potentials(), nptypes(0), dim(3), dt(1e-3), step(0) {}

void System::set_ntypes(int n){										//set the number of particles types
	nptypes = n;
	masses.resize(n, 1.0);											//resize array of particle masses based on the number of types
	potentials.resize(get_n_interactions());						//resize array of potentials based on the number of possible interactions 
	std::fill(potentials.begin(), potentials.end(), Potential_ptr());	//fill array of potentials with null pointer; after "set_ntypes(n)" call, potentials must be reassigned
}

int System::get_ntypes() const {
	return nptypes;
}

int System::get_n_interactions() const {							//return the number of interaction types for the given number of particle types: (type1, type1), (type2, type2), (type1, type2), etc.
	return nptypes * (nptypes + 1) / 2;
}

bool System::rigid(int id) const {									//check if a particle with "id" group (khist) is considered as rigid
	return std::find(groups_rigid.begin(), groups_rigid.end(), particles[id].khist) != groups_rigid.end();
}

bool System::verify() {												//if system is set up correctly, returns true
	bool result = true;
	if (masses.size() < nptypes) {									//check the size of array of masses				
		std::cout << "Particle masses are not set correctly" << std::endl;
		result = false;
	}
	if (get_n_interactions() < potentials.size()) {					//check the size of array of interactions
		std::cout << "Unspecified interatomic interactions are detected" << std::endl;
		result = false;
	}
	for (int i = 0; i < nptypes; i++) {								
		for (int j = i + 1; j < nptypes; j++) {
			int id = get_interaction_id(i + 1, j + 1);				//type id starts form 1, therefore "+ 1"
																	//check if EAM interaction is set correctly
			if ((potentials[id] != nullptr) && potentials[id]->eam() && //if cross interaction ij is eam, interactions ii and jj must be also eam
				((potentials[i] == nullptr) || !potentials[i]->eam() ||
				(potentials[j] == nullptr) || !potentials[j]->eam())) {
				std::cout << "EAM cross interaction is set incorrectly" << std::endl;
				result = false;
			}
																	//check if 3 body interaction is set correctly
			if ((potentials[id] != nullptr) && potentials[id]->threebody() && //if cross interaction ij is 3 body, interactions ii and jj must be also 3 body
				((potentials[i] == nullptr) || !potentials[i]->threebody() ||
				(potentials[j] == nullptr) || !potentials[j]->threebody())) {
				std::cout << "3 Body cross interaction is set incorrectly" << std::endl;
				result = false;
			}
		}
	}

	if (step == 0) {
		bool wrong_type = false;									//check if the particle types in the file are correct
		for (auto& particle : particles) {							//loop over all particles
			if (particle.type > nptypes) {
				particle.type = 1;									//assign type = 1 if type is larger than the number of particle types 
				wrong_type = true;
			}
		}
		if (wrong_type) {
			std::cout << "The number of particle types in the loaded file is larger ";
			std::cout << "than number of types used in the simulation" << std::endl;
			std::cout << "Incorrect types are reassigned to the type 1" << std::endl;
		}
	}
	return result;
}

//print boundary into os stream
std::ostream & MSE6270_MD::operator << (std::ostream & os, Boundary boundaty) {
	if (boundaty == Boundary::periodic) os << "periodic";
	else os << "free";
	return os;
}

//print information about a system into "os" stream
std::ostream & MSE6270_MD::operator << (std::ostream & os, const System & system) {
	os << "Boundary: " << system.boundaries[0] << " " << system.boundaries[1] << " " << system.boundaries[2] << std::endl;
	os << "System dimensions: " << system.XL << " " << system.YL << " " << system.ZL << " A" << std::endl;
	os << "The total number of atoms: " << system.particles.size();
	if (system.get_ntypes() > 1) os << " (" << system.get_ntypes() << " atomic types)";
	os << std::endl;

	auto& rigid = system.groups_rigid;
	if (!rigid.empty()) {											//print a list of "rigid" atomic groups (khist)
		os << "Rigid: khist = " << rigid[0];
		for (int i = 1; i < rigid.size(); i++) os << ", " << rigid[i];
		os << std::endl;
	}

	if(system.dim == 2) os << "Two dimensional simulation" << std::endl;
	return os;
}

