/*!
* @file NSearch.h
* @author Maxim Shugaev
* @date 23 Jan 2018
* @brief This file contains a set of classes used for neighbour search */

#pragma once
#include "System.h"
#include "Potentials.h"
#include <vector>

namespace MSE6270_MD {
	/*! An interface that describes methods used for neighbour search */
	class NSearch_Base {
	public:
		using Narray = std::vector<std::vector<std::vector<int>>>;			//An abbreviated name for 3D std::vector [interaction_id][particle_id][neighbour]

		NSearch_Base();														//!< Initialize NN_Search_Base
		virtual ~NSearch_Base() = default;									//destructor must be virtual if inheritance is used
		int update();														//!< Update neighbour list
		const Narray& get_neighbours() const;										//!< Return neighbour list represented as a 3D std::vector [interaction_id][particle_id][neighbour]
		void set_update_frequency(int n);									//!< Assign the frequency of neighbour list update (every "n"-th steps)
		int set_system(System const* system);								//!< Assign system (MD setup) used for neighbour search
		void set_Rskin(double Rskin);										//!< Assign the vlaue of skin layer that is added to cutoff distance during neighbour search
		virtual std::string info() const = 0;								//!< Return a string with information about a neighbour search algorithm
	protected:
		virtual int update_impl() = 0;										//!< Implementation of the neighbour list update. Override this function to implement a particular algorithm

		Narray neighbours;													/*!< Neighbour list represented as a 3D std::vector [interaction_id][particle_id][neighbour].
																				 Since each std::vector expands automatically based on the number of neighbours of a particular atom,
																				 the consumed memory can be significantly reduced in comparison with a simple 3D array, even if numerous particle types are considered. */
		System const* system;												/*!< Pointer to a system (MD setup) */
		int update_n;														/*!< Frequency of neighbour list update */
		double Rskin;														/*!< Thickness of a skin layer that is added to cutoff distance during neighbour search */
	};

	/*! This class implements neighbour search algorithm. This algorithm checks all pairs in the system and adds a particular pair to the list if distance in the pair is less than Rcut + Rskin. 
		The complexity of the algorithm is O(N^2). For more details check Frenkel & Smit pp. 363-367. */
	class NB_list : public NSearch_Base {
	public:
		virtual ~NB_list() = default;										//destructor must be virtual if inheritance is used
		virtual std::string info() const override;							//!< Return a string with information about a neighbour search algorithm
	protected:
		virtual int update_impl() override;									//!< Implementation of the neighbour list update. Override this function to implement a particular algorithm
	};
}
