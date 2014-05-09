/*!
 @file MD.cpp
 @author Danielle Rothermel
 @date 4/27/14

 @brief MD contains all of the methods involved in an MD run.

 */

#include "MD.h"
using namespace std;

MD::MD(double box_length_angstroms, int num_molecules, double init_temp, int sim_length){

	//Error Check
	if( sqrt((double)num_molecules) !=  (double)((int) sqrt((double)num_molecules)_){
		printf("Number of Molecules not a perfect square ==> problems will occur, you should end sim");
	}

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

//Sets the initial pos in lattice, random angle
void MD::init_position(){
	int sqN = (int)sqrt((double)m_N);
	double dbet = (double)n_L/(double)sqN;

	double xloc = dbet/2;
	double yloc = dbet/2;
	int moln = 0;
	for(int r=0; r < sqN; r++){
		for (int c=0; c < sqN; c++){
			moln = r*sqN+c;
			Point loc = {
				xloc,
				yloc,
				((double)(rand() % 360))*M_PI/180.0
			};
			Point vel = {
				(double)(rand() % 100)-50,
				(double)(rand() % 100)-50,
				(double)(rand() % 100)-50
			}
			m_molecs[moln].update_pos(loc);
			m_molecs[moln].update_vel(vel);

			xloc += dbet;
		}
		xloc = dbet/2;
		yloc += dbet;
	}

	zero_total_momentum();
	scale_init_temps();
}

void MD::evolve(){

}

// Updates each molecule's neighbors list and number of neighbors
void MD::update_neighbors(){

	Water* wp = &m_molecs[0];
	Water* op = &m_molecs[0];
	double dx = 0;
	double dy = 0;
	double dist = 0;
	double halfL = m_L/2.0;
	int m;
	int o;

	//for all molecules, set m_num_neigh to zero
	for(m = 0; m < m_N; m++){
		m_num_neigh[m] = 0;
	}

	//for all molecules except last
	for(m = 0; m < m_N-1; m++){
		wp = &m_molecs[m];
		//for all molecules of higher number
		for(o = m+1; o < m_N; o++){
			op = &m_molecs[o];

			dx = op->m_x.x - wp->m_x.x;
			if (dx > halfL) dx = dx - m_L;
			if (dx < -halfL) dx = dx + m_L;

			dy = op->m_x.y - wp->m_x.y;
			if (dy > halfL) dy = dy - m_L;
			if (dy < -halfL) dy = dy + m_L;

			dist = sqrt(dx*dx+dy*dy);

			// If w and o are close enough to each other
			if(dist < RC){
				m_neigh[m*m_N + m_num_neigh[m]] = wo; //set m's next neighbor pointer to wo
				m_neigh[o*m_N + m_num_neigh[o]] = wp; //set o's next neighbor pointer to wp


				m_num_neigh[m] += 1; //update the number of neighbors
				m_num_neigh[o] += 1;
			}
		}
	}

}

void MD::single_timestep(){

}

double MD::calc_sumforces(){

}

double MD::calc_sumtorques(){

}

double MD::calc_PE(){

}

// Calculate KE by summing the KEx, KEy, KEw for each molecule
double MD::calc_KE(){
	double ret = 0;
	Water* wp = &m_molecs[0];
	for(int m = 0; m < m_N; m++){
		wp = &m_molecs[m];
		ret += .5*18*wp->m_v.x*wp->m_v.x; // (1/2)m vx^2
		ret += .5*18*wp->m_v.y*wp->m_v.y; // (1/2)m vy^2
		ret += .5*wp->m_I*wp->m_v.th*wp->m_v.th; // (1/2)I w^2
	}

	return ret;
}

// Calculate current temperature using Equipartition Theorem
double MD::calc_T(){
	double ke = calc_KE();

	return (ke/m_N)*(2.0/3.0);
}

// Calculates the mean of every form of velocity and subtracts it from each velocity
void MD::zero_total_momentum(){ // O(2N)
	double tot_vx = 0;
	double tot_vy = 0;
	double tot_w = 0;

	Water* wp = &m_molecs[0];
	for(int m = 0; m < m_N; m++){
		wp = &m_molecs[m];
		tot_vx += wp->m_v.x;
		tot_vy += wp->m_v.y;
		tot_w += wp->m_v.th;
	}

	double meanx = tot_vx/(double)m_N;
	double meany = tot_vy/(double)m_N;
	double meanth = tot_w/(double)m_N;
	for(int m = 0; m < m_N; m++){
		wp = &m_molecs[m];
		wp->m_v.x = wp->m_v.x - meanx;
		wp->m_v.y = wp->m_v.y - meany;
		wp->m_v.th = wp->m_v.th - meanth;
	}

	//while it might not be necessary to zero total th, the total th momentum should be constant and an easy way to
	//ensure that is to make it constant and zero.
}

// Rescales the velocities of all molecules to set temperature to initial temp
void MD::scale_init_temps(){
	double tc = calc_T();
	double scale = sqrt(m_T/tc);

	Water* wp = &m_molecs[0];
	for(int m = 0; m < m_N; m++){
		wp = &m_molecs[m];
		wp->m_v.x = wp->m_v.x*scale;
		wp->m_v.y = wp->m_v.y*scale;
		wp->m_v.th = wp->m_v.th*scale;
	}
}