/*!
 * @file Force.h
 * @author Maxim Shugaev
 * @date 24 Jan 2018
 * @brief This file contains Force class implementing algorithms for calculation of interatomic interaction energies and
 * forces for pair, EAM, and 3 body potentials */

#pragma once
#include "NSearch.h"
#include "System.h"

namespace MSE6270_MD {
/*! A class implementing algorithms for calculation of interatomic interaction energies and forces for pair, EAM, and 3
 * body potentials */
class Force {
 public:
  using Narray = NSearch_Base::Narray;  // An abbreviated name for neighbour list (3D std::vector
                                        // [interaction_id][particle_id][neighbour])
  static int eam_density(System& system,
                         const Narray& neighbours);           /**< Update EAM density, should be called before force().
                                                               *	 "neighbours" is a neighbour list (3D vector
                                                               *neighbours[interaction_id][atom_id][neighbour_id]) */
  static int force(System& system, const Narray& neighbours); /**< Calculate forces and potential energies, "system" is
                                                               *MD simulation setup, "neighbours" is a neighbour list
                                                               *(3D vector
                                                               *neighbours[interaction_id][atom_id][neighbour_id]) */
  static int force_3b(System& system, const Narray& neighbours); /**< Calculate forces and potential energies for 3 body
                                                                  *interaction part, if one is present. "system" is MD
                                                                  *simulation setup, "neighbours" is a neighbour list
                                                                  *(3D vector
                                                                  *neighbours[interaction_id][atom_id][neighbour_id])*/
 private:
  static void apply_bounday(
      DPVector3& dr,
      const System& system);  //!< Correct the interatomic distance if periodic boundary conditions are applied
};

// The following functions are defined in .h file to allow inlining
//================================================================================================
inline void Force::apply_bounday(
    DPVector3& dr,
    const System& system) {  // correct the interatomic distance if periodic boundary conditions are applied
  for (int k = 0; k < DPVector3::dim; k++) {  // dim = 3
    if (system.boundaries[k] == Boundary::periodic) {
      if (dr.data[k] > 0.5 * system.size[k])
        dr.data[k] -= system.size[k];
      else if (dr.data[k] < -0.5 * system.size[k])
        dr.data[k] += system.size[k];
    }
  }
}
}  // namespace MSE6270_MD