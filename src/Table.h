/*!
 * @file Table.h
 * @author Maxim Shugaev
 * @date 23 Jan 2018
 * @brief This file contains Table class used to tabulate function and its derivative */

#pragma once
#include <iostream>
#include <string>
#include <vector>

namespace MSE6270_MD {
/*! A class used to tabulate function and its derivative */
class Table {
 public:
  Table();  //!< Initialize Spline_linear
  int set(
      const std::vector<double>& f, const std::vector<double>& df, double arg_min,
      double arg_max); /**< Set arrays representing tabulated values of function "f(x)" and its derivative "df(x)/dx"
                        *	 Argument (x) is assumed to be uniformally distributed within arg_min and arg_max */
  double get(
      double arg) const; /**< Return a value of function corresponding to the argument "arg". For out of bound "arg",
                          *	 the result is extrapolated according to the function derivative at the boundary */
  double get_d(double arg) const;  //!< Return a value of derivative corresponding to the argument "arg". For out of
                                   //!< bound "arg", the result corresponds to the dirivative at the boundary
  void get_tot(double arg, double& f, double& df)
      const;  //!< Return both a value of function "f" and its derivative "df" corresponding to the argument "arg"
  double min() const;                      //!< Return the minimum value of argument used for building the spline
  double max() const;                      //!< Return the maximum value of argument used for building the spline
  const std::string& get_name() const;     //!< Return the name of the table
  void set_name(const std::string& name);  //!< Set the name of the table
  int size() const;                        //!< Return the number of points in the table

  void load(std::istream& in);         //!< Load table from a file stream "in"
  void save(std::ostream& out) const;  //!< Save table to a file stream "out"
  void print() const;                  //!< Print table to a file
 private:
  std::vector<double> data, d_data; /*!< Array of values */
  double arg_min;                   /*!< Minimum value argument */
  double arg_max;                   /*!< Maximum value argument */
  std::string name;                 /*!< Name of the table */
};

// The following functions are defined in .h file to alow inlining
//================================================================================================

// Return a value of function corresponding to the argument "arg". For out of bound "arg", the result is extrapolated
// according to the function derivative at the boundary
inline double Table::get(double arg) const {
  int n = data.size();
  if (n < 1) return 0.0;
  if (arg < arg_min)
    return data.front() +
           d_data.front() * (arg - arg_min);  // out of bound check, return a value based on the derivative value
  else if (arg >= arg_max)
    return data.back() + d_data.back() * (arg - arg_max);
  else {                                                             // linear interpolation of tabulated values
    double idx_d = (n - 1) * (arg - arg_min) / (arg_max - arg_min);  // idx_d is a float point index corresponding to
                                                                     // arg
    int idx = static_cast<int>(idx_d);  // idx_d lies in the [idx, idx + 1] interval
    double correction = idx_d - idx;    // correction represents the position of idx_d within [idx, idx + 1]
    return data[idx] * (1.0 - correction) +
           data[idx + 1] *
               correction;  // estimation of val[idx_d] based on the linear interpolation within [idx, idx + 1] interval
  }
}

// Return a value of derivative corresponding to the argument "arg". For out of bound "arg", the result corresponds to
// the derivative at the boundary
inline double Table::get_d(double arg) const {
  int n = data.size();
  if (n < 1) return 0.0;
  if (arg < arg_min)
    return d_data.front();  // out of bound check, return values of derivatives at the boundaries
  else if (arg >= arg_max)
    return d_data.back();
  else {  // linear interpolation of tabulated values
    double idx_d = (n - 1) * (arg - arg_min) / (arg_max - arg_min);
    int idx = static_cast<int>(idx_d);
    double correction = idx_d - idx;
    return d_data[idx] * (1.0 - correction) + d_data[idx + 1] * correction;
  }
}

// Return both a value of function "f" and its derivative "df" corresponding to the argument "arg"
inline void Table::get_tot(double arg, double& f, double& df) const {
  int n = data.size();
  if (n < 1) return;
  if (arg < arg_min) {  // out of bound check, return a value based on the derivative value
    f = data.front() + d_data.front() * (arg - arg_min);
    df = d_data.front();
  } else if (arg >= arg_max) {
    f = data.back() + d_data.back() * (arg - arg_max);
    df = d_data.back();
  } else {                                                           // linear interpolation of tabulated values
    double idx_d = (n - 1) * (arg - arg_min) / (arg_max - arg_min);  // idx_d is a float point index corresponding to
                                                                     // arg
    int idx = static_cast<int>(idx_d);  // idx_d lies in the [idx, idx + 1] interval
    double correction = idx_d - idx;    // correction represents the position of idx_d within [idx, idx + 1]
    f = data[idx] * (1.0 - correction) +
        data[idx + 1] *
            correction;  // estimation of val[idx_d] based on the linear interpolation within [idx, idx + 1] interval
    df = d_data[idx] * (1.0 - correction) + d_data[idx + 1] * correction;
  }
  return;
}
}  // namespace MSE6270_MD
