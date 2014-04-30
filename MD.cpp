/*!
 @file MD.cpp
 @author Danielle Rothermel
 @date 4/27/14

 @brief MD contains all of the methods involved in an MD run.

 */

#include "MD.h"
using namespace std;

MD::MD(double box_length_angstroms, int num_molecules, int init_temp, int sim_length){
	m_N = num_molecules; //num particles
	m_L = box_length_angstroms; //box length angstroms
	m_T = init_temp; //initial temp
	m_simL = sim_length; //number of fs to run the sim

	m_step = 0; //current step number
	m_molecs = new Water[m_N](); //an array of N Waters
	m_xcoords = new Point[m_N*(m_simL/COORDTS)]; //a collapsed array of arrays of N points, access with array[i*m_N+j]
	m_vcoords = new Point[m_N*(m_simL/COORDTS)]; //a collapsed array of arrays of N points
	m_temps = new double[m_simL/ENERTS]; //an array of sim_length/ENERTS temperatures
	m_KE = new double[m_simL/ENERTS]; //an array of sim_length/ENERTS KEs
	m_PE = new double[m_simL/ENERTS]; //an array of sim_length/ENERTS PEs

	m_neigh = new Water*[m_N*m_N]; //an NxN array of neighbor pointers to Water molecules in m_molecs array
	m_num_neigh = new double[m_N]; //an N array of number of neighbors for each molecule
}

MD::~MD(){
	delete m_neigh;
	delete m_num_neigh;
	delete m_PE;
	delete m_KE;
	delete m_temps;
	delete m_vcoords;
	delete m_xcoords;
	delete m_molecs;	
}

void MD::init_position(){

}

void MD::init_velocity(){

}

void MD::evolve(){

}

void MD::update_neighbors(){

}

void MD::single_timestep(){

}

double MD::calc_sumforces(){

}

double MD::calc_sumtorques(){

}

double MD::calc_KE(){

}

double MD::calc_PE(){

}

double calc_T(){

}
