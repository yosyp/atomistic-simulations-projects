/*!
 * @file Table.cpp
 * @author Maxim Shugaev
 * @date 23 Jan 2018
 * @brief This file contains implementation of methods of Table class used to tabulate function and its derivative */

#include "Table.h"
#include <fstream>
#include <iostream>

using namespace MSE6270_MD;

// Table==================================================================================================

Table::Table() : arg_min(0.0), arg_max(0.0), data(), d_data(), name("unnamed") {}

// Set arrays representing tabulated values of function "f(x)" and its derivative "df(x)/dx". Argument(x) is assumed to
// be uniformally distributed within arg_min and arg_max
int Table::set(const std::vector<double>& f, const std::vector<double>& df, double arg_min_, double arg_max_) {
  int n = f.size();
  if ((n < 1) || (df.size() != n)) return EXIT_FAILURE;
  arg_min = arg_min_;
  arg_max = arg_max_;  // assign arg_min, and arg_max
  data.resize(n);
  d_data.resize(n);
  std::copy(f.begin(), f.end(), data.begin());      // set values of the function f
  std::copy(df.begin(), df.end(), d_data.begin());  // set values of the function derivative df/dx
  return EXIT_SUCCESS;
}

const std::string& Table::get_name() const { return name; }

void Table::set_name(const std::string& name_) { name = name_; }

int Table::size() const { return data.size(); }

void Table::load(std::istream& in) {  // load table from a file stream "in"
  int n;
  in >> name >> n >> arg_min >> arg_max;  // read header
  // std::cout << "loading " << name << " table: Np = " << n << ", arg_min = " << arg_min << ", arg_max = " << arg_max
  // << std::endl;
  data.resize(n);  // resize tables
  d_data.resize(n);
  for (int i = 0; i < n; i++) {  // read table
    double x;
    in >> x >> data[i] >> d_data[i];
  }
}

void Table::save(std::ostream& out) const {  // save table to a file stream "out"
  std::cout << "saving " << name << " table" << std::endl;
  int n = data.size();                                                       // number of points in the table
  out << name << " " << n << " " << arg_min << " " << arg_max << std::endl;  // write header
  for (int i = 0; i < n; i++) {                                              // write table
    double x = arg_min + (arg_max - arg_min) * i / (n - 1);
    out << x << " " << data[i] << " " << d_data[i] << std::endl;
  }
}

void Table::print() const {  // print table to a file
  std::cout << "printing " << name << " table" << std::endl;
  std::ofstream out(name + ".tab");
  int n = data.size();           // number of points in the table
  for (int i = 0; i < n; i++) {  // write table
    double x = arg_min + (arg_max - arg_min) * i / (n - 1);
    out << x << " " << data[i] << " " << d_data[i] << std::endl;
  }
}
