/*!
* @file Potentials.h
* @author Maxim Shugaev
* @date 23 Jan 2018
* @brief This file contains classes used for description of interatomic interactions */

#pragma once
#include "Table.h"
#include <string>
#include <vector>
#include <algorithm>

namespace MSE6270_MD {
	/*! A class that implements a general behaviour of interatomic interaction potential. Children classes are responsible for assignment of values to the tables  and setting corresponding flags. */
	class Potential_Base {
	public:
		virtual ~Potential_Base() = default;								//destructor must be virtual if inheritance is used
	
		const Table& U() const;												//!< Return a table with pair interaction function U
		const Table& f() const;												//!< Return a table with electron density function f, used for EAM potentials
		const Table& F() const;												//!< Return a radial dependent part f(r_ij) of 3 body potential
		const Table& tb_f() const;											//!< Return an angular dependent part g[cos(theta)] of 3 body potential
		const Table& tb_g() const;

		bool full_list() const;												//!< Return "true" if full list neighbour list is required for a particular potential
		bool eam() const;													//!< Return "true" for EAM potential
		bool threebody() const;												//!< Return "true" for 3 body potential
		double r_cut() const;												//!< Return cutoff distance 
		void print() const;													//!< Print tables into a file: r[A], U[eV], dU[eV/A]... check function implementation for more details
		bool is_build() const;												//!< Return "true" if potential is built
		virtual std::string info() const = 0;								//!< Return a string with information about the potential
		double get_mass() const;											//!< Return mass, which is set if the potential is created by "Potential_LJ(material)"
	protected:
		Potential_Base();													//!< Initialize Potential_Base, creation of Potential_Base objects is forbidden

		bool full_list_flag;												/*!< Full list flag */
		bool eam_flag;														/*!< EAM potential flag */
		bool tb_flag;														/*!< 3 body potential flag */
		double cut_off;														/*!< Cutoff distance */
		std::string name;													/*!< Name of a particular interatomic potential */
		bool build_flag;													/*!< Status flag, "true" if potential is built */
		double mass;														/*!< Mass of atoms */

		Table U_table, f_table, F_table;									//Tables for U, dU/dr, f, df/dr, F, dF/dro, linear interpolation is used
		Table tb_f_table, tb_g_table;										//Tables used in three body potentials
	};

	/*! A class generating Lennard-Jones potential */
	class Potential_LJ : public Potential_Base {							//expand functionality of Potential_Base
	public:
		/*! This is an enum class describing avalible materials */
		enum class Material {Ar, Al, Ca, Au, Pb, Ni, Pd, Pt, Ag, Cu, Cr, Fe, Li, Mo, W, Na, K, Si, Ge};
		Potential_LJ(Material material, double Rcut = 2.5);					//!< Initialize Potential_LJ for a given "material", Rcut is cutoff distance in units of sigma
		Potential_LJ(double epsilon, double sigma, const std::string& name, double Rcut = 2.5);	/**< Initialize Potential_LJ for given "epsilon" and "sigma", 
																			* "name" is the name of the interaction, Rcut is cutoff distance in units of sigma
																			* (2.5 is the default value) */
		virtual ~Potential_LJ() = default;									//destructor must be virtual if inheritance is used

		double get_epsilon() const;											//!< Return value of epsilon
		double get_sigma() const;											//!< Return value of sigma
		virtual std::string info() const override;							//!< Return a string with information about the potential
	protected:
		Potential_LJ() = default;
		void set_parameters(Material material);								//!< Assign "epsilon", "sigma", "name", and "mass" for a given material
		virtual void build();												//!< Assign values of U and dU tables for given LJ potential parameters

		double epsilon, sigma;												//parameters of LJ potential
		double Rcut0;														//cutoff distance in units of sigma
	};

	/*! A class generating Lennard-Jones potential with cut off function. 
		See S. D. Stoddard and J. Ford, Phys. Rev. A 8, 1504, 1973. for more details*/
	class Potential_LJ_cut : public Potential_LJ {
	public:
		using Material = Potential_LJ::Material;							//make Material available in Potential_LJ_cut namespace
		Potential_LJ_cut(Material material, double Rcut = 2.5);				//!< Initialize Potential_LJ_cut for a given "material", Rcut is cutoff distance in units of sigma
		Potential_LJ_cut(double epsilon, double sigma, const std::string& name, double Rcut = 2.5);	/**< Initialize Potential_LJ for given "epsilon" and "sigma",
																			* "name" is the name of the interaction, Rcut is cutoff distance in units of sigma
																			* (2.5 is the default value)*/
		virtual std::string info() const override;							//!< Return a string with information about the potential
	protected:
		virtual void build() override;										//!< Assign values of U and dU tables for given LJ potential parameters with using cutoff function
	};

	/*! A class generating a pair potential based on a *.tab table. Checke Table::load for more details. Do not forget to assign correct a mass to the system when set the potential */
	/*	The file format is the following:
		table name, number of points(np), arg_min, arg_max
		np lines: r, U, dU / dr */
	class Potential_Pair_tab : public Potential_Base {						//expand functionality of Potential_Base
	public:
		Potential_Pair_tab(const std::string& filename);					//!< Initialize Potential_Pair based on a table from a file with a name "filename"
		virtual ~Potential_Pair_tab() = default;							//destructor must be virtual if inheritance is used
		virtual std::string info() const override;							//!< Return a string with information about the potential
	};

	/*! A class generating EAM potential based on *.alloy table. It is LAMMPS format, and tables for different potentials are available online. 
		You may check https://www.ctcms.nist.gov/potentials/ to find the tables. But before using, try to plot a loaded potential to make sure that everything is Ok. 
		Also, you may compare corresponding cohesive energies and lattice parameters at zero temperature with ones given on the web site. */
	class Potential_EAM_alloy : public Potential_Base {						//expand functionality of Potential_Base
	public:
																			/*! Generate potential based on an *.alloy file with a name "filename". 
																				"El1" and "El2" are names of elements that should be loaded form the file.
																				Data corresponding to the selected elements must be present in the file. */
		Potential_EAM_alloy(const std::string& filename, const std::string& El1, const std::string& El2);
		virtual ~Potential_EAM_alloy() = default;							//destructor must be virtual if inheritance is used

		virtual std::string info() const override;							//!< Return a string with information about the potential
	};

	/*! A class generating EAM potential based on *.tab table (a table format used in the Fortran version of the code). */
	/*	The format of *.tab files :
		number of particle types(N)
		-----> N blocks for each type
			mass, atomic_volume, material_id
			(atomic_volume and material_id are not used in the current code and introduced for backward compatibility of the table format)
			U : table name, number of points(np), arg_min, arg_max
			U_data(np) : r, U, dU / dr
			f : table name, number of points(np), arg_min, arg_max
			f_data(np) : r, f, df / dr
			F : table name, number of points(np), arg_min, arg_max
			F_data(np) : ro, F, dF / dro
		------
		-----> N*(N - 1) / 2 blocks for cross interaction potential(t1 - t2, t1 - t3, ... t1 - tN, t2 - t3, ...)
			U : table name, number of points(np), arg_min, arg_max
			U_data(np) : r, U, dU / dr
		------ */
	class Potential_EAM_tab : public Potential_Base {						//expand functionality of Potential_Base
	public:
																			/*! Generate potential based on an *.tab file with a name "filename".
																			"id1" and "id2" are indexes of elements listed in the file. For example, CuAg.tab includes 2 elements: 
																			(1,1) -> Cu-Cu, (2,2) -> Ag-Ag, (1,2) -> Cu-Ag */
		Potential_EAM_tab(const std::string& filename, int id1, int id2);
		virtual ~Potential_EAM_tab() = default;								//destructor must be virtual if inheritance is used

		virtual std::string info() const override;							//!< Return a string with information about the potential
	protected:
	};

	/*! A class generating Stillinger-Weber potential */
	class Potential_SW : public Potential_Base {							//expand functionality of Potential_Base
	private :
		static const int n_materials = 2;									//number of materials is Material enum class
	public:
		/*! This is an enum class describing available materials */
		enum class Material { Si, Ge };
		Potential_SW(Material El1, Material El2);							//!< Build potential for a given pair of elements (El1, El2)
		virtual ~Potential_SW() = default;									//destructor must be virtual if inheritance is used

		virtual std::string info() const override;							//!< Return a string with information about the potential
	protected:
		void set_parameters(Material El1, Material El2);					//!< Assign parameters for a given pair of elements (El1, El2)
		void build();														//!< Build potential tables for given parameters
		Potential_SW();

		double sigma, epsilon, A, B, gamma, lamda, p, q, cos_c, lamda_p;	//parameters of SW potential
	};


	//The following functions are defined in .h file to allow inlining
	//================================================================================================

	inline double Table::min() const {
		return arg_min;
	}

	inline double Table::max() const {
		return arg_max;
	}

	inline const Table & Potential_Base::U() const {
		return U_table;
	}

	inline const Table & Potential_Base::f() const {
		return f_table;
	}

	inline const Table & Potential_Base::F() const 	{
		return F_table;
	}

	inline const Table & Potential_Base::tb_f() const {
		return tb_f_table;
	}

	inline const Table & Potential_Base::tb_g() const {
		return tb_g_table;
	}

	inline bool Potential_Base::full_list() const {
		return full_list_flag;
	}

	inline bool Potential_Base::eam() const {
		return eam_flag;
	}

	inline bool Potential_Base::threebody() const {
		return tb_flag;
	}

	inline double Potential_Base::r_cut() const {
		return cut_off;
	}
}
