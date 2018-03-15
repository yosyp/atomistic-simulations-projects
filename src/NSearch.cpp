/*!
 * @file NSearch.cpp
 * @author Maxim Shugaev
 * @date 23 Jan 2018
 * @brief The file contains implementation of classes used for neighbour search and declared in NSearch.h */

#include "NSearch.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>

using namespace MSE6270_MD;

// NSearch_Base==================================================================================================

NSearch_Base::NSearch_Base() : update_n(10), system(nullptr), Rskin(1.0) {}

int NSearch_Base::update() {
  if (system == nullptr) {
    std::cout << "System (MD setup) is not assigned to the neighbour search algorithm" << std::endl;
    return EXIT_FAILURE;
  }
  if (system->step % update_n == 0)
    return update_impl();  // update neighbour list every "update_n" steps
  else
    return EXIT_SUCCESS;
}

const NSearch_Base::Narray& NSearch_Base::get_neighbours() const { return neighbours; }

void NSearch_Base::set_update_frequency(int n) { update_n = n; }

int NSearch_Base::set_system(System const* system_) {  // assign system (MD setup) used for neighbour search
  if (system_ == nullptr) return EXIT_FAILURE;         // exit if null pointer is passed
  system = system_;
  int n_interactions = system->get_n_interactions();  // number of possible interactions in the system, like (type1,
                                                      // type1), (type2, type2), (type1, type2), etc.
  int n = system->particles.size();   // number of particles
  neighbours.resize(n_interactions);  // set the size of the array equal to the number of particles
  for (int i = 0; i < n_interactions; i++)
    neighbours[i].resize(n);  // set first 2 dimensions in 3D std::vector [interaction_id][particle_id][neighbour]
                              // the size in the third dimension will be selected individually for each particle based
                              // on the number of neighbours
  return EXIT_SUCCESS;
}

void NSearch_Base::set_Rskin(double Rskin_) { Rskin = Rskin_; }

// NB_list==================================================================================================

std::string NB_list::info() const {  // return a string with information about a neighbour search algorithm
  return "Neigbor List algorithm: O(N^2) complexity, Rskin = " + std::to_string(Rskin) + " A, update every " +
         std::to_string(update_n) + " steps";
}

int NB_list::update_impl() {                                          // implementation of the neighbour list update
  if (neighbours.size() != system->get_ntypes()) set_system(system);  // reset system if number of types is updated
  int n_interactions = neighbours.size();  // number of possible interactions in the system, like (type1, type1),
                                           // (type2, type2), (type1, type2), etc.
  int n = system->particles.size();  // number of particles
  for (int i = 0; i < n_interactions; i++) {
    for (int j = 0; j < n; j++)
      neighbours[i][j].resize(0);  // empty neighbour list without deallocation of reserved memory, which allow to avoid
                                   // time consuming memory allocation
  }
  auto& potentials = system->potentials;  // abbreviated names
  auto& particles = system->particles;
  std::vector<double> cut_off2(n_interactions);  // array of square cutoff distances
  for (int i = 0; i < n_interactions; i++) {
    if (potentials[i] != nullptr)
      cut_off2[i] = pow(potentials[i]->r_cut() + Rskin, 2);
    else
      cut_off2[i] = 0.0;  // if particular atoms do not interact, set zero
  }

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {  // loop over all pairs in the system
      int interaction =
          system->get_interaction_id(particles[i].type, particles[j].type);  // get interaction id for a pair (i,j)
      auto dr = particles[i].pos - particles[j].pos;                         // distance between atoms i and j
      for (int k = 0; k < DPVector3::dim; k++) {                             // dim = 3
        if (system->boundaries[k] ==
            Boundary::periodic) {  // correct distance between atoms i and j if periodic boundary conditions are applied
          if (dr.data[k] > 0.5 * system->size[k])
            dr.data[k] -= system->size[k];
          else if (dr.data[k] < -0.5 * system->size[k])
            dr.data[k] += system->size[k];
        }
      }
      // skip pair if dr^2 > cutoff^2
      if (dr.length_s() > cut_off2[interaction]) continue;
      neighbours[interaction][i].push_back(j);  // add particle j to the i-th particle neighbour list
                                                // add particle i to the j-th particle neighbour list if full list is
                                                // requested by (i,j) cross interaction potential
      if (potentials[interaction]->full_list()) neighbours[interaction][j].push_back(i);
    }
  }
  return EXIT_FAILURE;
}
