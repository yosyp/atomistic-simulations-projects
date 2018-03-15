/*!
* @file Potentials.cpp
* @author Maxim Shugaev
* @date 24 Jan 2018
* @brief This file contains implementation of classes used for description of interatomic interactions, 
   which declared in Potentials.h */

#include "Potentials.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>

using namespace MSE6270_MD;

//Potential_Base==================================================================================================
																			//set default values
Potential_Base::Potential_Base() : full_list_flag(false), eam_flag(false), tb_flag(false), cut_off(0.0), name(), build_flag(false), mass(1.0) {} 

bool Potential_Base::is_build() const {
	return build_flag;
}

double Potential_Base::get_mass() const {
	return mass;
}

void Potential_Base::print() const {										//print tables into a file: r[A], U[eV], dU/dr[eV/A]...
	std::cout << "Printing " << name << " potential" << std::endl;
	std::cout << std::boolalpha;											//output bool as true/false
	std::cout << "eam = " << eam_flag << ", 3 body = " << tb_flag << ", full list = " << full_list_flag << std::endl;
	std::ofstream out(name + ".pot");										//open output file
	static const int n_points = 500;										//number of points printed in the file
	double r1 = U_table.min(), r2 = U_table.max();							//minimum and maximum values of r
	double ro1, ro2, costh1, costh2;
	if (eam_flag) ro1 = F_table.min(), ro2 = F_table.max();					//minimum and maximum values of electron density for EAM potentials
	if (tb_flag) costh1 = tb_g_table.min(), costh2 = tb_g_table.max();		//minimum and maximum values of electron density for EAM potentials

	for (int i = 0; i < n_points; i++) {
		double r = r1 + (r2 - r1) * i/ (n_points - 1);
		out << r << " " << U_table.get(r) << " " << U_table.get_d(r);		//print r[A], U[eV], dU/dr[eV/A]
		if (eam_flag) {
			double ro = ro1 + (ro2 - ro1) * i / (n_points - 1);				//if EAM, print f, df/dr, ro, F[eV], dF/r[eV/A]
			out << " " << f_table.get(r) << " " << f_table.get_d(r) << " " << ro << " " << F_table.get(ro) << " " << F_table.get_d(ro);
		}
		if (tb_flag) {
			double costh = costh1 + (costh2 - costh1) * i / (n_points - 1); //if 3 body potential, print f (radial dependent part of 3 body term), df/dr, 
																			//cos[theta], g(angular dependent part of the interaction), dg/d(cos[theta])
			out << " " << tb_f_table.get(r) << " " << tb_f_table.get_d(r) << " " << costh << " " << tb_g_table.get(costh) << " " << tb_g_table.get_d(costh);
		}
		out << std::endl;
	}
}

//Potential_Potential_LJ==================================================================================================

Potential_LJ::Potential_LJ(Material material, double Rcut) : Rcut0(Rcut){	//initialize Potential_LJ for a given "material", Rcut is the cutoff distance in units of sigma
	set_parameters(material);												//assign "epsilon", "sigma", "name", and "mass" for a given material
	build();																//build LJ potential for given parameters
}

Potential_LJ::Potential_LJ(double epsilon_, double sigma_, const std::string & name_, double Rcut)
	: Rcut0(Rcut), epsilon(epsilon_), sigma(sigma_){
	name = name_;
	build();																//build LJ potential for given parameters
}

double Potential_LJ::get_epsilon() const {
	return epsilon;
}

double Potential_LJ::get_sigma() const {
	return sigma;
}

std::string Potential_LJ::info() const {									//return a string with information about the potential
	return "Lennard-Jones potential for " + name + " interaction, Rcut = " + std::to_string(cut_off) + " A";
}

//assign "epsilon", "sigma", "name", and "mass" for a given material
void Potential_LJ::set_parameters(Material material) {
	switch (material) {														//select potential properties for a given material
	case Material::Ar:														//Ar : D.V.Matyushov and R. Schmid, J.Chem.Phys. 104, 8627 (1996)
		name = "Ar-Ar"; epsilon = 0.0103; sigma = 3.405; mass = 39.948;
		break;
	case Material::Al:														//Metals : T.Halicioglu and G.M.Pound, Phys.Stat.Sol.A 30, 619 (1975)
		name = "Al-Al"; epsilon = 0.3922; sigma = 2.62; mass = 27.0;
		break;
	case Material::Ca:
		name = "Ca-Ca"; epsilon = 0.2152; sigma = 3.60; mass = 40.1;
		break;
	case Material::Au:
		name = "Au-Au"; epsilon = 0.4415; sigma = 2.637; mass = 197.0;
		break;
	case Material::Pb:
		name = "Pb-Pb"; epsilon = 0.2364; sigma = 3.197; mass = 207.2;
		break;
	case Material::Ni:
		name = "Ni-Ni"; epsilon = 0.5197; sigma = 2.282; mass = 58.7;
		break;
	case Material::Pd:
		name = "Pd-Pd"; epsilon = 0.4267; sigma = 2.52; mass = 106.4;
		break;
	case Material::Pt:
		name = "Pt-Pt"; epsilon = 0.6817; sigma = 2.542; mass = 195.1;
		break;
	case Material::Ag:
		name = "Ag-Ag"; epsilon = 0.3448; sigma = 2.644; mass = 107.9;
		break;
	case Material::Cu:
		name = "Cu-Cu"; epsilon = 0.4094; sigma = 2.338; mass = 63.6;
		break;
	case Material::Cr:
		name = "Cr-Cr"; epsilon = 0.5020; sigma = 2.336; mass = 52.0;
		break;
	case Material::Fe:
		name = "Fe-Fe"; epsilon = 0.5264; sigma = 2.321; mass = 55.9;
		break;
	case Material::Li:
		name = "Li-Li"; epsilon = 0.2054; sigma = 2.839; mass = 6.94;
		break;
	case Material::Mo:
		name = "Mo-Mo"; epsilon = 0.8379; sigma = 2.551; mass = 95.9;
		break;
	case Material::W:
		name = "W-W"; epsilon = 1.0678; sigma = 2.562; mass = 183.85;
		break;
	case Material::Na:
		name = "Na-Na"; epsilon = 0.1379; sigma = 3.475; mass = 23.0;
		break;
	case Material::K:
		name = "K-K"; epsilon = 0.1144; sigma = 4.285; mass = 39.1;
		break;
	case Material::Si:
		name = "Si-Si"; epsilon = 3.3900; sigma = 2.0936; mass = 28.09;
		break;
	case Material::Ge:
		name = "Ge-Ge"; epsilon = 2.7400; sigma = 2.1827; mass = 72.59;
		break;
	default:
		name = "unknown"; epsilon = 1.0; sigma = 1.0; mass = 1.0;
	}
}

void Potential_LJ::build() {												//build LJ potential for given parameters
	std::cout << "Building LJ " << name << " potential" << std::endl;
	std::cout << "sigma = " << sigma << ", epsilon = " << epsilon << std::endl;

	cut_off = Rcut0*sigma;													//cutoff distance
	static const int n_points = 8192;										//number of points in the table
	std::vector<double> U(n_points), dU(n_points);							//temporal arrays
	double r1 = 0.5*sigma, r2 = cut_off;									//set minimum and maximum r
	for (int i = 0; i < n_points; i++) {
		double r = r1 + (r2 - r1)*i / (n_points - 1);
		U[i] = 4.0*epsilon*(pow(sigma/ r, 12) - pow(sigma / r, 6));			//U and dU for Lennard-Jones potential
		dU[i] = -24.0*epsilon*(2.0*pow(sigma / r, 12) - pow(sigma / r, 6))/r;
	}
	U_table.set(U, dU, r1, r2);												//build table
	U_table.set_name("U_" + name);
	build_flag = true;
}

//Potential_Potential_LJ==================================================================================================

Potential_LJ_cut::Potential_LJ_cut(Material material, double Rcut) {		//initialize Potential_LJ for a given "material", Rcut is the cutoff distance in units of sigma
	Rcut0 = Rcut;
	set_parameters(material);												//assign "epsilon", "sigma", "name", and "mass" for a given material
	build();																//build LJ potential for given parameters with using cutoff function
}

Potential_LJ_cut::Potential_LJ_cut(double epsilon_, double sigma_, const std::string & name_, double Rcut) {
	epsilon = epsilon_;
	sigma = sigma_;
	name = name_;
	Rcut0 = Rcut;
	build();																//build LJ potential for given parameters with using cutoff function
}

std::string Potential_LJ_cut::info() const {								//return a string with information about the potential
	return "Lennard-Jones potential with cutoff function for " + name + " interaction, Rcut = " + std::to_string(cut_off) + " A";
}

void Potential_LJ_cut::build() {											//assign values of U and dU tables for given LJ potential parameters with using cutoff function
	std::cout << "Building LJ " << name << " potential with cutoff function" << std::endl;
	std::cout << "sigma = " << sigma << ", epsilon = " << epsilon << std::endl;

	cut_off = Rcut0*sigma;													//2.5 sigma cutoff distance
	static const int n_points = 8192;										//number of points in the table
	std::vector<double> U(n_points), dU(n_points);							//temporal arrays
	double r1 = 0.5*sigma, r2 = cut_off;									//set minimum and maximum r
	for (int i = 0; i < n_points; i++) {
		double r = r1 + (r2 - r1)*i / (n_points - 1);
		double sr6 = pow(sigma / r, 6);
		double src6 = pow(sigma / cut_off, 6);
		double XTRC2 = pow(r / cut_off, 2);
																			//U and dU for Lennard-Jones potential with cutoff function
		U[i] = 4.0*epsilon*(sr6*(sr6 - 1.0) + 3.0*src6*(2.0*src6 - 1.0)*XTRC2 - src6 * (7.0*src6 - 4.0));
		dU[i] = -(epsilon * (48.0 / r)*sr6*(sr6 - 0.5) - epsilon * 24.0*src6*(2.0*src6 - 1.0)*r / (cut_off*cut_off));
	}
	U_table.set(U, dU, r1, r2);												//build U table
	U_table.set_name("U_" + name);
	std::cout << "Rcut = " << cut_off << std::endl;
	build_flag = true;
}

//Potential_Pair_tab==================================================================================================

Potential_Pair_tab::Potential_Pair_tab(const std::string & filename) {
	std::ifstream in(filename);												//open file
	if (!in.is_open()) {
		std::cout << "Cannot open " << filename << std::endl;
		return;
	}
	U_table.load(in);
	cut_off = U_table.max();
	std::cout << "Rcut = " << cut_off << std::endl;
	build_flag = true;
}

std::string Potential_Pair_tab::info() const {								//return a string with information about the potential
	return name + " interaction, Rcut = " + std::to_string(cut_off) + " A";
}

//Potential_EAM_alloy==================================================================================================
//genrate EAM potential based on an *.alloy file with a name "filename", "El1" and "El2" are names of elements that should be loaded form the file
Potential_EAM_alloy::Potential_EAM_alloy(const std::string & filename, const std::string & El1, const std::string & El2) {
	std::ifstream in(filename);												//open file
	if (!in.is_open()) {
		std::cout << "Cannot open " << filename << std::endl;
		return;
	}
	std::string buf;
	for (int i = 0; i < 4; i++) std::getline(in, buf);						//skip 3 lines with information about the file
	std::stringstream str_stream(buf);
	int ntypes;
	str_stream >> ntypes;													//number of types described in the files
	std::vector<std::string> materials(ntypes);
	for (int i = 0; i < ntypes; i++) str_stream >> materials[i];			//read a list of elements described in the file
	
	int n_points_r, n_points_ro;
	double dro, dr, r_max;
	in >> n_points_ro >> dro >> n_points_r >> dr >> r_max;					//the number of points in tables tabulated based on ro, dro,
																			//the number of points in tables tabulated based on r, dr, r_max
																			//[0,n_points_ro*dro] - range of ro, [r_max - n_points_r*dr,r_max] - range of r
	std::vector<double> masses(ntypes);										//define several temporal arrays used during file reading
	int n_interactions = ntypes * (ntypes + 1) / 2;							//the number of possible interactions for "ntypes" particle types
	std::vector<std::vector<double>> F(ntypes), dF(ntypes), f(ntypes), df(ntypes), U(n_interactions), dU(n_interactions);
	for (auto &el : F) el.resize(n_points_ro);
	for (auto &el : dF) el.resize(n_points_ro);
	for (auto &el : f) el.resize(n_points_r);
	for (auto &el : df) el.resize(n_points_r);
	for (auto &el : U) el.resize(n_points_r);
	for (auto &el : dU) el.resize(n_points_r);

																			//the remaining part of the file has the following structure:												
	for (int i = 0; i < ntypes; i++) {										//ntypes blocks:
		std::getline(in, buf);
		std::getline(in, buf); str_stream.clear(); str_stream.str(buf);
		int tmp_i; 
		str_stream >> tmp_i >> masses[i];									//element number, mass, lattice constant, lattice type

		for (int j = 0; j < n_points_ro; j++) in >> F[i][j];				//F table (n_points_ro value), f table (n_points_r value)
		for (int j = 0; j < n_points_r; j++) in >> f[i][j];
		F[i][0] = 0.0;														//energy of a single atom in vacuum is zero
	}
	for (int i = 0; i < n_interactions; i++) {								//n_interactions blocks:
		for (int j = 0; j < n_points_r; j++) {								//U*r table (n_points_r value)
			in >> U[i][j];
			double r = r_max - dr * (n_points_r - j);
			U[i][j] /= std::max(r, dr);										//file written for U*r
		}
	}

	for (int i = 0; i < ntypes; i++) {										//calculate derivatives for F and f
		for (int j = 1; j < n_points_ro - 1; j++) dF[i][j] = (F[i][j + 1] - F[i][j - 1]) / (2.0*dro);	//dF(i)/dr = (F(i+1) - F(i-1))/2dr
		dF[i][0] = dF[i][1]; dF[i][n_points_ro - 1] = dF[i][n_points_ro - 2];	//set derivatives at boundaries
		for (int j = 1; j < n_points_r - 1; j++) df[i][j] = (f[i][j + 1] - f[i][j - 1]) / (2.0*dr);
		df[i][0] = df[i][1]; df[i][n_points_r - 1] = 0.0;					//set derivatives at boundaries
	}
	for (int i = 0; i < n_interactions; i++) {								//calculate derivatives for U
		for (int j = 1; j < n_points_r - 1; j++) dU[i][j] = (U[i][j + 1] - U[i][j - 1]) / (2.0*dr);
		dU[i][0] = dU[i][1]; dU[i][n_points_r - 1] = 0.0;					//set derivatives at boundaries
	}

	auto idx1 = std::find(materials.begin(), materials.end(), El1);			//find "El1" and "El2" among loaded elements
	auto idx2 = std::find(materials.begin(), materials.end(), El2);
	if ((idx1 == materials.end()) || (idx2 == materials.end())) {			//cannot find "El1" or "El2"
		std::cout << "Cannot find an appropriate element in the file " << filename << std::endl;
		std::cout << "Element1 = " << El1 << ", Element2 = " << El2 << std::endl;
		std::cout << "Availible elements: ";
		for (int i = 0; i < ntypes; i++) std::cout << materials[i] << ", ";
		std::cout << std::endl;
		return;
	}
	int material1 = idx1 - materials.begin();								//get ids of "El1" and "El2" element in the read file
	int material2 = idx2 - materials.begin();

	name = materials[material1] + "-" + materials[material2];				//generate name "El1"-"El2"
	std::cout << filename << ": loading " << name << " potential" << std::endl;
	mass = masses[material1];												//set mass equal to the atomic mass of El1
	double r1 = r_max - dr*n_points_r, r2 = r_max;							//minimum and maximum values of r
	double ro1 = 0.0, ro2 = dro * n_points_ro;								//minimum and maximum values of ro

	int interaction = 0;													//calculate id for U table in .alloy file
	for (int i = 0; i <= std::min(material1, material2); i++) {				//for 2 material types: (type1, type1)->0, (type1, type2)->1, (type2, type2)->2
		for (int j = i; j < ntypes; j++) {									//this rule is different form one used in the current code that would be: (type1, type1)->0, (type2, type2)->1, (type1, type2)->2													
			if ((i == std::min(material1, material2)) && (j == std::max(material1, material2))) break;
			interaction++;
		}
	}

	U_table.set(U[interaction], dU[interaction], r1, r2);					//build U table
	U_table.set_name("U_" + El1 + "-" + El2);
																			//F, dF/dro, f, df/dr are needed only if "El1"="El2"; 
	if (material1 == material2) {											//cross interaction is described based on F, dF/dro, f, df/dr of pure elements
		f_table.set(f[material1], df[material1], r1, r2);					//build f table
		f_table.set_name("f_" + El1 + "-" + El2);
		F_table.set(F[material1], dF[material1], ro1, ro2);					//build F table
		F_table.set_name("F_" + El1 + "-" + El2);
	}
	cut_off = r2;
	std::cout << "Rcut = " << cut_off << std::endl;
	eam_flag = true;														//it is EAM potential
	build_flag = true;
}


std::string Potential_EAM_alloy::info() const {								//return a string with information about the potential
	return "EAM potential for " + name + " interaction, Rcut = " + std::to_string(cut_off) + " A";
}

//Potential_EAM_tab==================================================================================================

Potential_EAM_tab::Potential_EAM_tab(const std::string & filename, int id1_, int id2_) {
	std::ifstream in(filename);												//open file
	if (!in.is_open()) {
		std::cout << "Cannot open " << filename << std::endl;
		return;
	}
	int id1 = std::min(id1_, id2_);											//ordered pair of indexes
	int id2 = std::max(id1_, id2_);
	int ntypes;
	in >> ntypes;
	if (id2 > ntypes) {
		std::cout << "The number of materials specified in " << filename << " is less than passed index = " << id2 << std::endl;
		return;
	}
	name = std::to_string(id1) + "-" + std::to_string(id2) + "_in_" + filename;
	std::replace(name.begin(), name.end(), '/', '_');						//replace "/" by "_" in the name
	std::cout << filename << ": loading " << id1 << "-" << id2 << " potential" << std::endl;

	int n_interactions = ntypes * (ntypes + 1) / 2;
	std::vector<Table> U(n_interactions), f(ntypes), F(ntypes);
	std::vector<double> m(ntypes);

	for (int i = 0; i < ntypes; i++) {
		double atomic_volume;												//atomic_volume and material_id are not used in the current code
		int material_id;													//and introduced for backward compatibility of the table format
		in >> m[i] >> atomic_volume >> material_id;
		U[i].load(in);														//load tables for i-i interactions
		f[i].load(in);
		F[i].load(in);
	}
	int idx = ntypes;
	int cross_id;
	for (int i = 1; i <= ntypes; i++) {
		for (int j = i + 1; j <= ntypes; j++) {
			U[idx].load(in);												//load tables for cross interactions
			if ((i == id1) && (j == id2)) cross_id = idx;
			idx++;
		}
	}

	if (id1 == id2) {														//assign tables
		U_table = U[id1 - 1];
		f_table = f[id1 - 1];
		F_table = F[id1 - 1];
	}
	else {
		U_table = U[cross_id];
	}

	mass = m[id1 - 1];
	cut_off = U_table.max();
	std::cout << "Rcut = " << cut_off << std::endl;
	eam_flag = true;														//it is EAM potential
	build_flag = true;
}

std::string Potential_EAM_tab::info() const {								//return a string with information about the potential
	return "EAM potential for " + name + " interaction, Rcut = " + std::to_string(cut_off) + " A";
}

//Potential_SW==================================================================================================

Potential_SW::Potential_SW(Material El1, Material El2) {
	set_parameters(El1, El2);
	build();
}

std::string Potential_SW::info() const {									//return a string with information about the potential
	return "Stillinger-Weber potential for " + name + " interaction, Rcut = " + std::to_string(cut_off) + " A";
}

Potential_SW::Potential_SW() : sigma(1.0), epsilon(0.0), A(0.0), B(0.0), gamma(0.0), lamda(0.0), p(0.0), q(0.0) {} //set default values

void Potential_SW::set_parameters(Material El1, Material El2) {				//assign parameters for a given pair of elements (El1, El2)
	if ((El1 == Material::Si) && (El2 == Material::Si)) {					//F.H.Stillinger and T.A.Weber, Phys.Rev.B 31, 5262-5271 (1985)
		name = "Si-Si";
		mass = 28.0855;
		epsilon = 2.167222; sigma = 2.0951;	//A
		p = 4.0; q = 0.0; A = 15.2779; B = 0.60222;
		lamda = 21.0; gamma = 1.2; cut_off = 3.7712; 
		cos_c = -1.0 / 3.0; lamda_p = lamda;
	}
	else if ((El1 == Material::Ge) && (El2 == Material::Ge)) {				//K.Ding and H.Andersen.Phys.Rev.B 34, 6987 (1986)
		name = "Ge-Ge";
		mass = 72.64;
		epsilon = 1.93; sigma = 2.181;
		p = 4.0; q = 0.0; A = 13.6056; B = 0.60222;
		lamda = 31.0; gamma = 1.2; cut_off = 3.9258; 
		cos_c = -1.0 / 3.0; lamda_p = lamda;
	}
	else {
		name = "Si-Ge";														//M.Laradji, D.P.Landau and B.Dunweg.Phys.Rev.B 51 4894 (1995)
		epsilon = 2.0427; sigma = 2.13805;
		p = 4.0; q = 0.0; A = 14.40013; B = 0.60222;
		gamma = 1.2; cut_off = 3.8437; 
		cos_c = -1.0 / 3.0; lamda_p = sqrt(31.0*21.0); //sqrt(lamda_Si*lamda_Ge)
	}
}

void Potential_SW::build() {
	std::cout << "Building SW " << name << " potential" << std::endl;

	static const int n_points = 8192;										//number of points in the table
																			//build pair part
	std::vector<double> U(n_points), dU(n_points);							//temporal arrays
	double r1 = 0.5*sigma, r2 = cut_off;									//set minimum and maximum r
	for (int i = 0; i < n_points - 1; i++) {
		double r = r1 + (r2 - r1)*i / (n_points - 1);
		double rs = sigma / r;
		U[i] = A*(B*pow(rs, p) - pow(rs, q))*exp(sigma/(r - cut_off));
		dU[i] = -A*(B*p*pow(rs, p + 1.0) - q*pow(rs, q + 1.0))*exp(sigma / (r - cut_off)) / sigma - U[i] * sigma/pow(r - cut_off,2);
	}
	U[n_points - 1] = 0.0; dU[n_points - 1] = 0.0;							//set the last point manually, otherwise there will be dividing by zero
	U_table.set(U, dU, r1, r2);												//build U table
	U_table.set_name("U_" + name);
																			
																			//build 3 body part, introduce these tables to make the description general enough for other potentials, like Tersoff
	std::vector<double> tb_f(n_points), tb_df(n_points);					//radial dependent part
	for (int i = 0; i < n_points - 1; i++) {
		double r = r1 + (r2 - r1)*i / (n_points - 1);
		tb_f[i] = sqrt(epsilon*lamda_p)*exp(sigma*gamma / (r - cut_off));
		tb_df[i] = -tb_f[i] * sigma*gamma / pow(r - cut_off, 2);
	}
	tb_f[n_points - 1] = 0.0; tb_df[n_points - 1] = 0.0;					//set the last point manually, otherwise there will be dividing by zero
	tb_f_table.set(tb_f, tb_df, r1, r2);									//build a table for radial dependent part of 3 body interaction
	tb_f_table.set_name("tb_f_" + name);

	std::vector<double> tb_g(n_points), tb_dg(n_points);					//angular dependent part g[cos(theta)]
	double cos1 = -1.0, cos2 = 1.0;
	for (int i = 0; i < n_points; i++) {
		double cos_theta = cos1 + (cos2 - cos1)*i / (n_points - 1);
		tb_g[i] = pow(cos_theta - cos_c, 2);
		tb_dg[i] = 2.0*(cos_theta - cos_c);
	}
	tb_g_table.set(tb_g, tb_dg, cos1, cos2);								//build a table for angular dependent part of 3 body interaction
	tb_g_table.set_name("tb_g_" + name);

	tb_flag = true;
	full_list_flag = true;
	build_flag = true;
}

