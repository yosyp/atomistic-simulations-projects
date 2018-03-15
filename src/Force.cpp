/*!
 * @file Force.cpp
 * @author Maxim Shugaev
 * @date 24 Jan 2018
 * @brief The file contains implementation of algorithms for calculation of interatomic interaction energies and forces
 * for pair, EAM, and 3 body potentials, declared in Force.h */

#include "Force.h"
#include <cmath>

using namespace MSE6270_MD;

int Force::eam_density(System& system, const Narray& neighbours) {  // update EAM density
  int n_interactions = system.get_n_interactions();  // number of different interatomic interactions: (type1, type1),
                                                     // (type2, type2), (type1, type2), etc.
  int n = system.particles.size();       // number of particles
  auto& particles = system.particles;    // an abbreviated name
  auto& potentials = system.potentials;  // an abbreviated name

  //	Calculation of electron density only, applied only for EAM potentials
  for (int interaction = 0; interaction < n_interactions;
       interaction++) {  // loop over all interactomic interactions present in the system, like (type1, type1), (type2,
                         // type2), (type1, type2), etc.
    auto& potential = potentials[interaction];  // an abbreviated name
    if ((potential == nullptr) || (!potential->eam()))
      continue;  // check if a particular interaction involves calculation of electron density; if potential is not set,
                 // particles do not interact
    double Rcut2 = pow(potential->r_cut(), 2);  // square cutoff distance
    for (int i = 0; i < n; i++) {
      int n_neighbours = neighbours[interaction][i].size();  // number of neighbours
      for (int j_id = 0; j_id < n_neighbours; j_id++) {
        int j = neighbours[interaction][i][j_id];
        auto& p1 = particles[i];
        auto& p2 = particles[j];
        auto dr = p2.pos - p1.pos;  // vector distance between particles i and j
        apply_bounday(dr, system);  // correct the interatomic distance if periodic boundary conditions are applied
        double r2 = dr.length_s();
        if (r2 > Rcut2) continue;  // skip the pair if the distance is larger than cutoff distance
        double r = sqrt(r2);

        p1.ro += potentials[p2.type - 1]->f().get(r);  // add electron density to particle1; particle type ids start
                                                       // from 1 wile indexes in array start form 0, therefore "-1"
        if (!potentials[interaction]->full_list()) {  // add electron density to particle2 if half list is used
          p2.ro += potentials[p1.type - 1]->f().get(r);
        }
      }
    }
  }

  //	Electron contribution to the total energy, applied only for EAM potentials
  for (int i = 0; i < n; i++) {
    auto& potential = system.potentials[particles[i].type - 1];  // interaction potential; particle type ids start from
                                                                 // 1 wile indexes in array start form 0, therefore "-1"
    if ((potential != nullptr) &&
        (potential->eam())) {  // skip if interaction of atoms i is not described by EAM potential
      particles[i].Ep += potential->F().get(particles[i].ro);   // add F(ro) contribution to the potential energy
      particles[i].ro = potential->F().get_d(particles[i].ro);  // replace variable ro by dF(ro), to avoid calculation
                                                                // of dF/dro for each atomic pair during force
                                                                // calculation
    }
  }
  return EXIT_SUCCESS;
}

int Force::force(System& system, const Narray& neighbours) {  // calculate forces and potential energies
  int n_interactions = system.get_n_interactions();  // number of different interatomic interactions: (type1, type1),
                                                     // (type2, type2), (type1, type2), etc.
  int n = system.particles.size();       // number of particles
  auto& particles = system.particles;    // an abbreviated name
  auto& potentials = system.potentials;  // an abbreviated name

  //	Pair force calculation
  for (int interaction = 0; interaction < n_interactions;
       interaction++) {  // loop over all interactomic interactions present in the system, like (type1, type1), (type2,
                         // type2), (type1, type2), etc.
    auto& potential = potentials[interaction];  // an abbreviated name
    if (potential == nullptr) continue;         // if potential is not set, particles do not interact
    double Rcut2 = pow(potential->r_cut(), 2);  // square cutoff distance
    for (int i = 0; i < n; i++) {
      int n_neighbours = neighbours[interaction][i].size();  // number of neighbours
      for (int j_id = 0; j_id < n_neighbours; j_id++) {
        int j = neighbours[interaction][i][j_id];
        auto& p1 = particles[i];
        auto& p2 = particles[j];
        auto dr = p2.pos - p1.pos;  // vector distance between particles i and j
        apply_bounday(dr, system);  // correct the interatomic distance if periodic boundary conditions are applied
        double r2 = dr.length_s();
        if (r2 > Rcut2) continue;  // skip the pair if the distance is larger than cutoff distance
        double r = sqrt(r2);

        double force, U;
        potential->U().get_tot(
            r, U, force);  // get values of U and dU/dr from the table
                           // if the considered interaction is EAM, add an electron density contribution to the force
        if (potential->eam())
          force += p1.ro * potentials[p2.type - 1]->f().get_d(r) +
                   p2.ro * potentials[p1.type - 1]->f().get_d(r);  // dF(ro)/dro reassigned to ro several lines above

        DPVector3 force_v = -(force / r) * dr;  // vector force
        p1.Ep += 0.5 * U;
        p1.force -= force_v;
        p1.stress -= force_v | dr;  // sress_xx += force_x*dx ...

        if (!potential->full_list()) {  // if half list is used, add corresponding contributions to the second particle
          p2.Ep += 0.5 * U;
          p2.force += force_v;
          p2.stress -= force_v | dr;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

int Force::force_3b(System& system, const Narray& neighbours) {  // calculate forces and potential energies for 3 body
                                                                 // interaction part, if one is present
  int n_interactions = system.get_n_interactions();  // number of different interatomic interactions: (type1, type1),
                                                     // (type2, type2), (type1, type2), etc.
  auto& particles = system.particles;  // abreviated names
  auto& potentials = system.potentials;
  int n = particles.size();  // number of particles

  for (int interaction_ij = 0; interaction_ij < n_interactions;
       interaction_ij++) {  // loop over all interactomic interactions present in the system, like (type1, type1),
                            // (type2, type2), (type1, type2), etc.
    auto& potential_ij = potentials[interaction_ij];
    if ((potential_ij == nullptr) || (!potential_ij->threebody()))
      continue;  // check if a particular interaction is defined and involves calculation of 3 body term
    double Rcut2_ij = pow(potential_ij->r_cut(), 2);  // square cutoff distance
    for (int interaction_ik = interaction_ij; interaction_ik < n_interactions;
         interaction_ik++) {  // loop over interactomic interactions present in the system
      auto& potential_ik = potentials[interaction_ik];
      if ((potential_ik == nullptr) || (!potential_ik->threebody()))
        continue;  // check if a particular interaction is defined and involves calculation of 3 body term
      double Rcut2_ik = pow(potential_ik->r_cut(), 2);  // square cutoff distance

      for (int i = 0; i < n; i++) {
        int n_neighbors_ij = neighbours[interaction_ij][i].size();  // number of neighbors for a particular interaction
        int n_neighbors_ik = neighbours[interaction_ik][i].size();
        if (n_neighbors_ik == 0) continue;
        int type_i = particles[i].type;  // type of particle i
        for (int j_id = 0; j_id < n_neighbors_ij; j_id++) {
          int j = neighbours[interaction_ij][i][j_id];
          auto dr_ij = particles[j].pos - particles[i].pos;  // vector distance between particles i and j
          apply_bounday(dr_ij, system);  // correct the interatomic distance if periodic boundary conditions are applied
          double r2_ij = dr_ij.length_s();  // square distance
          if (r2_ij > Rcut2_ij) continue;   // skip the pair if the distance is larger than cutoff distance
          double r_ij = sqrt(r2_ij);
          double f_ij, df_ij;
          potential_ij->tb_f().get_tot(r_ij, f_ij, df_ij);  // radial dependent parts of the interaction
          int k_id_start = (interaction_ij == interaction_ik)
                               ? j_id + 1
                               : 0;  // start "k" loop from "j_id + 1" for interaction between atoms of the same kind
          for (int k_id = k_id_start; k_id < n_neighbors_ik; k_id++) {  // otherwise start from 0
            int k = neighbours[interaction_ik][i][k_id];
            auto dr_ik = particles[k].pos - particles[i].pos;  // vector distance between particles i and k
            apply_bounday(dr_ik,
                          system);  // correct the interatomic distance if periodic boundary conditions are applied
            double r2_ik = dr_ik.length_s();  // square distance
            if (r2_ik > Rcut2_ik) continue;   // skip the pair if the distance is larger than cutoff distance
            double r_ik = sqrt(r2_ik);
            // jik triplet of atoms is found
            auto dr_jk = dr_ik - dr_ij;
            double i_divider = 1.0 / (r_ij * r_ik);  // a temporal variable used in the following calculations
            double cos_jik = dr_ij * dr_ik * i_divider;

            double g_cos_jik, dg_cos_jik, f_ik, df_ik;
            potentials[type_i - 1]->tb_g().get_tot(cos_jik, g_cos_jik,
                                                   dg_cos_jik);  // angle dependent part of the interaction
            potential_ik->tb_f().get_tot(
                r_ik, f_ik, df_ik);  // radial dependent parts of the interaction, types start from 1, therefore "-1"

            double U = f_ij * f_ik * g_cos_jik;  // potential energy
            double prefactor = f_ij * f_ik * dg_cos_jik * i_divider;
            double force_jk = -prefactor;  // forces between different pairs of atoms in the triplet
            double force_ij = (prefactor * (r_ij - cos_jik * r_ik) + df_ij * f_ik * g_cos_jik) * r_ik * i_divider;
            double force_ik = (prefactor * (r_ik - cos_jik * r_ij) + f_ij * df_ik * g_cos_jik) * r_ij * i_divider;

            particles[i].Ep += U;  // add calculated values of energy, force, and sterss
            particles[i].force += force_ij * dr_ij + force_ik * dr_ik;
            particles[j].force += force_jk * dr_jk - force_ij * dr_ij;
            particles[k].force += -force_jk * dr_jk - force_ik * dr_ik;
            particles[i].stress += force_ij * (dr_ij | dr_ij) + force_ik * (dr_ik | dr_ik);
            particles[j].stress += force_ij * (dr_ij | dr_ij) + force_jk * (dr_jk | dr_jk);
            particles[k].stress += force_ik * (dr_ik | dr_ik) + force_jk * (dr_jk | dr_jk);
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}