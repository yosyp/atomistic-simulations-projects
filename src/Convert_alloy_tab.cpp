/*!
* @file Convert_alloy_tab.cpp
* @author Maxim Shugaev
* @date 24 Jan 2018
* @brief A short program for conversion of *.alloy to *.tab file */

#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include "Potentials.h"

using namespace MSE6270_MD;

int main() {
	static const std::string in_file = "potentials/CuAg.eam.alloy";				//input file
	static const std::string out_file = "potentials/CuAg.tab";					//output file
	static const std::vector<std::string> elements = {"Cu", "Ag"};				//list of elements written in the output file
	using Potential_ptr = std::shared_ptr<Potential_Base>;
	int ntypes = elements.size();												//number of partcle types
	int n_interactions = ntypes * (ntypes + 1) / 2;
	std::vector<std::vector<Potential_ptr>> potentials(ntypes, std::vector<Potential_ptr>(ntypes)); //empty array
	for (int i = 0; i < ntypes; i++) {
		for (int j = i; j < ntypes; j++) {										//load potentials
			potentials[i][j] = std::make_shared<Potential_EAM_alloy>(in_file, elements[i], elements[j]);
			if (!potentials[i][j]->is_build()) {								//check if a potential is loaded
				std::cout << "cannot load " << elements[i] << "-" << elements[j] << " potential"<< std::endl;
				return EXIT_FAILURE;
			}
		}
	}

	std::ofstream out(out_file);												//write output
	out << ntypes << std::endl;
	for (int i = 0; i < ntypes; i++) {
		out << potentials[i][i]->get_mass() << " " << 0.0 << " " << 0 << std::endl;
		potentials[i][i]->U().save(out);
		potentials[i][i]->f().save(out);
		potentials[i][i]->F().save(out);
	}
	for (int i = 0; i < ntypes; i++) {
		for (int j = i + 1; j < ntypes; j++) {
			potentials[i][j]->U().save(out);
		}
	}

	return EXIT_SUCCESS;
}