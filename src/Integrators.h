/*!
 * @file Integrators.h
 * @author Maxim Shugaev
 * @date 24 Jan 2018
 * @brief This file contains classes implementing algorithms for numerical integration of equations of motion in MD */

#pragma once
#include <array>
#include <vector>
#include "System.h"

namespace MSE6270_MD {
/*! An interface that declares methods used for numerical integration of equations of motion in MD */
class Integrator_Base {
 public:
  virtual ~Integrator_Base() = default;  // destructor must be virtual if inheritance is used

  int set_system(System* system);  //!< Assign a system (MD setup) to the integrator
  virtual void pre_step();  //!< Perform a pre step of integration (before force calculation) used in some integration
                            //!< algorithms like predictor-corrector
  virtual void step() = 0;  //!< Perform a step of integration
  virtual std::string info() const = 0;  //!< Return a string with information about an integration algorithm
 protected:
  Integrator_Base();    //!< Set default values
  virtual void init();  //!< Initialize internal variables of the integrator when system is assign

  System* system; /*!< Pointer to a system (MD setup) which the integrator is applied for */
};

/*! This class implements velocity Verlet algorithm  */
class Integrator_Verlet : public Integrator_Base {  // implement Integrator_Base interface
 public:
  virtual ~Integrator_Verlet() = default;  // destructor must be virtual if inheritance is used
  virtual void step() override;            //!< Perform a step of integration
  virtual std::string info() const;        //!< Return a string with information about an integration algorithm
};

/*! This class implements 5-th order Predictor-Corrector algorithm.
        See [Nordsieck A., Math.Comput.16, 22, 1962], p. 154 in [C. William Gear, Numerical initial value problems in
   ordinary differential equations, Prentice-Hall Inc., Englewood Cliffs, NJ, 1971], and [Allen & Tildesley, p.340] for
   more details */
class Integrator_Nord5 : public Integrator_Base {
 public:
  virtual ~Integrator_Nord5() = default;  // destructor must be virtual if inheritance is used
  virtual void pre_step() override;  //!< Perform a pre step of integration(before force calculation): predictor step
  virtual void step() override;      //!< Perform a step of integration: corrector step
  virtual std::string info() const;  //!< Return a string with information about an integration algorithm
 protected:
  static const double C[6];                /*!< Constants used during correction step */
  std::vector<std::array<DPVector3, 4>> q; /*!< High order derivatives of coordinates */

  virtual void init() override;  //!< Initialize internal variables of the integrator when system is assign
};
}  // namespace MSE6270_MD