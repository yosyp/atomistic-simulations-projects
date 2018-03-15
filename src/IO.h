/*!
* @file IO.h
* @author Maxim Shugaev
* @date 21 Jan 2018
* @brief This file contains a set of classes used for input/output.
*  It includes reading/writing restarting files, writing atomic snapshots, and reading from crystal generator pseudo file. */

#pragma once
#include "System.h"
#include <array>
#include <vector>
#include <string>

namespace MSE6270_MD {
	/*! A reader capable to load restarting (.data) files */
	class Reader_Data {
	public:
		Reader_Data(const std::string &filename);				//!< Initialize Reader_Data
		int load(System& system);								//!< Load system form a restarting (.data) file
	protected:
		std::string filename;
	};

	/*! A reader capable to load system from a crystal generator pseudo files */
	class Reader_CG {
	public:
		/*! This is an enum class describing a crystal structure */
		enum class Lattice { SC, BCC, BCC_CsCl, FCC, FCC_Ni3Al, Diamond, Diamond_GaAs };
		Reader_CG(const std::array<int, 3>& ncells, double a_lat, Lattice lattice); /**< Initialize Reader_CG,
																					* "ncells" is the number of unit cells {nx, ny, nz}, "a_lat" is the lattice constant,
																					* "lattice" is the crystal structure of the generated system */
		int load(System& system);								//!< Load system form a Crystal Generator pseudo file
	protected:
		std::array<int, 3> ncells;								/*!< The number of unit cells {nx, ny, nz} */
		std::vector<int> types;									/*!< An array describing the types of atoms within a unit cell a selected crystal structure */
		std::vector<DPVector3> basis;							/*!< An array describing the basis of a selected crystal structure */
		DPVector3 shift;										/*!< Shift atoms from the origin of the unit cell */
		double a_lat;											/*!< Lattice constant */
		int cell_n;												/*!< The number of atoms in the unit cell */
		std::string lattice_type;								/*!< The name of the generated crystal structure */
	};

	/*! A writer capable to save restarting (.data) files */
	class Writer_Data {	
	public:
		Writer_Data(const std::string &filename);				//!< Initialize Writer_Data, "filename" is the name of the saved file
		int save(const System& system);							//!< Save restarting (.data) file
	protected:
		std::string filename;									/*!< The name of the file */
	};

	/*! A writer capable to save atomic snapshots, (.d) files: step, khist, pos[0-2] [A], vel[0-2] [A/ps], Ep [eV], Ek [eV], T [K], Etot [eV], type */
	class Writer_Snapshot {
	public:
		Writer_Snapshot(const std::string &filename);			//!< Initialize Writer_Snapshot, "filename" is the name of the saved file
		int save(const System& system);							//!< Save atomic snapshot (.d) file
	protected:
		std::string filename;									/*!< The name of the file */
	};
}
