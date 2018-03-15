/*!
* @file System.h
* @author Maxim Shugaev
* @date 22 Jan 2018
* @brief The file describes classes-containers used in the simulations, including DPVector3 (3D double precision
vector), Particle (information about particles), and System (information about a particular computational system
simulated with MD)*/

#pragma once
#include <array>
#include <iostream>
#include <memory>
#include <vector>
#include "Potentials.h"

namespace MSE6270_MD {
// MSE6270_MD uses the following units:
// x - [A], t - [ps], m - [Da], E - [eV], F - [eV/A], P - [Pa]
// define several constants for unit conversion
static const double e = 1.60219e-19;       /*!< Electron charge */
static const double kb = 8.617385e-5;      /*!< Boltzman constant [eV/K] */
static const double m0 = 1.6605402e-27;    /*!< Atomic mass unit [kg] */
static const double enunit = 1e4 * m0 / e; /*!< A constant used to convert [Da*A^2/ps] to [eV]  */

/*! A class representing a 3D double precision vector and correponding vector operations */
class DPVector3 {
 public:
  DPVector3();                              //!< Construct {0, 0, 0} vector
  DPVector3(double X, double Y, double Z);  //!< Construct {X, Y, Z} vector

  DPVector3& operator+=(const DPVector3& rhs);                   //!< Increase current vector by a rhs vector
  DPVector3& operator-=(const DPVector3& rhs);                   //!< Decrease current vector by a rhs vector
  const DPVector3 operator+(const DPVector3& rhs) const;         //!< Sum of two vectors
  const DPVector3 operator-(const DPVector3& rhs) const;         //!< Residual of two vectors
  DPVector3& operator*=(double rhs);                             //!< Multiply current vector by a rhs scalar
  DPVector3& operator/=(double rhs);                             //!< Divide current vector by a rhs scalar
  const DPVector3 operator*(double rhs) const;                   //!< Product of the current vector and a rhs scalar
  const DPVector3 operator/(double rhs) const;                   //!< Quotient of the current vector and a rhs scalar
  friend DPVector3 operator*(double lhs, const DPVector3& rhs);  //!< Product of a lhs scalar and a rhs vector
  double operator*(const DPVector3& rhs) const;                  //!< Scalar product of two vectors
  const DPVector3 operator|(const DPVector3& rhs) const;         //!< Element wise product of two vectors
  double length_s() const;                                       //!< Square length of the vector
  double length() const;                                         //!< Length of the vector

  static const int dim = 3; /*!< The number of dimensions, equal to 3  */

  /** A union that allows a dual representation of the elements of the vector */
  union {
    struct {
      double x; /*!< x elements of the vector */
      double y; /*!< y elements of the vector */
      double z; /*!< z elements of the vector */
    };
    double data[dim]; /*!< Elements of the vector accessed as data[0], data[1], data[2] */
  };
};

/*! A class container representing a Particle */
class Particle {
 public:
  Particle();    //!< Initialize particle
  void reset();  //!< Set force, stress, Ep, Ek, ro equal to zero

  DPVector3 pos;    /*!< Position of the particle */
  DPVector3 vel;    /*!< Velocity of the particle */
  DPVector3 force;  /*!< Force acting on the particle */
  DPVector3 stress; /*!< s_xx, s_yy, s_zz components of the stress tensor */
  int type;         /*!< Type of the particle */
  int khist;        /*!< Group that the particle belong to */
  double Ep;        /*!< Potential energy of the particle */
  double Ek;        /*!< Kinetic energy of the particle */
  double ro;        /*!< Electron density used in EAM potentials */
 private:
};

/*! This is an enum class to describe the method used for pressure control */
enum class Boundary { free, periodic };
/*! Print boundary into os stream */
std::ostream& operator<<(std::ostream& os, Boundary boundaty);

/*! A class container with information about a particular computational system */
class System {
 public:
  using Potential_ptr = std::shared_ptr<Potential_Base>;  // An abbreviated name for a smart pointer to Potential_Base

  System();                                       //!< System particle
  void set_ntypes(int n);                         //!< Set the number of particles types
  int get_ntypes() const;                         //!< Return the number of particles types
  int get_interaction_id(int id1, int id2) const; /*! Return the interaction id for a (id1, id2) pair.
                                                                                                          For example,
                                                     if nptypes = 2, (1,1) -> 0, (2,2) -> 1, (1,2) -> 2. If nptypes = 3,
                                                     (1,1) -> 0, (2,2) -> 1, (3,3) -> 2, (1,2) -> 3, (1,3) -> 4, (2,3)
                                                     -> 5. */
  int get_n_interactions() const;  //!< Return the number of interaction types for the given number of particle types:
                                   //!< (type1, type1), (type2, type2), (type1, type2), etc.
  bool rigid(int id) const;        //!< Check if a particle with "id" group (khist) is considered as rigid
  bool verify();                   //!< If system is set up correctly, return true

  /** A union that allows a dual representation of system dimensions */
  union {
    struct {
      double XL; /*!< System size in X direction */
      double YL; /*!< System size in Y direction */
      double ZL; /*!< System size in Z direction*/
    };
    double size[3]; /*!< System dimensions accessed as size[0], size[1], size[2] */
  };
  /** A union that allows a dual representation of the position of the system center */
  union {
    struct {
      double Xcenter; /*!< X coordinate of the system center */
      double Ycenter; /*!< Y coordinate of the system center */
      double Zcenter; /*!< Z coordinate of the system center */
    };
    double center[3]; /*!< The position of the system center accessed as center[0], center[1], center[2] */
  };
  std::vector<Particle> particles;       /*!< An array of particles composing the system */
  std::array<Boundary, 3> boundaries;    /*!< Boundary conditions {Boundary_x, Boundary_y, Boundary_z} */
  std::vector<double> masses;            /*!< An array of masses corresponding to particular particle types */
  std::vector<int> groups_rigid;         /*!< An array of groups of atoms (khist) that are considered as rigid */
  std::vector<Potential_ptr> potentials; /*!< An array of smart pointers to potentials used in the simulation. The
                                            potential must be assigned for each cross interection type. Different
                                            potential types can be used together. Null pointer corresponds to
                                            noninteracting particles */

  int dim;     /*!< The dimensionality of the simulation (2 - 2D, 3 - 3D) */
  double time; /*!< Current time */
  int step;    /*!< Current step */
  double dt;   /*!< Time step */
 private:
  int nptypes; /*!< The number of particle types */
};

/*! Print information about a system into "os" stream */
std::ostream& operator<<(std::ostream& os, const System& system);

// The following functions are defined in .h file to allow inlining
//================================================================================================

inline DPVector3& DPVector3::operator+=(const DPVector3& rhs) {
  for (int i = 0; i < dim; i++) data[i] += rhs.data[i];
  return *this;
}

inline DPVector3& DPVector3::operator-=(const DPVector3& rhs) {
  for (int i = 0; i < dim; i++) data[i] -= rhs.data[i];
  return *this;
}

inline const DPVector3 DPVector3::operator+(const DPVector3& rhs) const { return DPVector3(*this) += rhs; };

inline const DPVector3 DPVector3::operator-(const DPVector3& rhs) const { return DPVector3(*this) -= rhs; }

inline DPVector3& DPVector3::operator*=(double rhs) {
  for (int i = 0; i < dim; i++) data[i] *= rhs;
  return *this;
}

inline DPVector3& DPVector3::operator/=(double rhs) {
  for (int i = 0; i < dim; i++) data[i] /= rhs;
  return *this;
}

inline const DPVector3 DPVector3::operator*(double rhs) const { return DPVector3(*this) *= rhs; }

inline const DPVector3 DPVector3::operator/(double rhs) const { return DPVector3(*this) /= rhs; }

inline DPVector3 MSE6270_MD::operator*(double lhs, const DPVector3& rhs) { return rhs * lhs; }

inline double DPVector3::operator*(const DPVector3& rhs) const {
  double sum = 0.0;
  for (int i = 0; i < dim; i++) sum += data[i] * rhs.data[i];
  return sum;
}

inline const DPVector3 DPVector3::operator|(const DPVector3& rhs) const {
  DPVector3 result;
  for (int i = 0; i < 3; i++) result.data[i] = data[i] * rhs.data[i];
  return result;
}

inline double DPVector3::length_s() const {
  double sum = 0.0;
  for (int i = 0; i < dim; i++) sum += data[i] * data[i];
  return sum;
}

inline double DPVector3::length() const { return sqrt(length_s()); }

inline int System::get_interaction_id(int t1,
                                      int t2) const {  // return id of a given interatomic interaction, id1 and id2 are
                                                       // ids of particle1 and particle2 in the considered pair of atoms
  if (t1 == t2) return t1 - 1;
  int id1 = std::min(t1, t2);
  int id2 = std::max(t1, t2);
  return id1 * (2 * nptypes - id1 - 1) / 2 + id2 - 1;
}
}  // namespace MSE6270_MD