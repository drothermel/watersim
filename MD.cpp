/*!
 @file MD.cpp
 @author Danielle Rothermel
 @date 4/27/14

 @brief MD contains all of the methods involved in an MD run.

 */

#include "MD.h"
using namespace std;

////////////// REMEMBER TO CALCULATE FORCES BEFORE FIRST TIME STEP!!!!!! //////////////////////////////

MD::MD(double box_length_angstroms, int num_molecules, double init_temp, int sim_length){

	//Error Check
	if( sqrt((double)num_molecules) !=  (double)((int) sqrt((double)num_molecules) )) {
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
	m_num_neigh = new int[m_N]; //an N array of number of neighbors for each molecule
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
	double dbet = (double)m_L/(double)sqN;

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
			};
			m_molecs[moln].m_moln = moln;
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

void MD::single_timestep(){
	Water* wp;
	int m;
	double ax, ay, alpha;
	double posx, posy, theta;

	//First half of v step
	for(m = 0; m < m_N; m++){
		wp = &m_molecs[m];

		// get (d^2 _ / d _^2)
		ax = (*wp).m_f.fx / 18.0; // sum(F_x)/m = a_x
		ay = (*wp).m_f.fy / 18.0; // sum(F_y)/m = a_y
		alpha = (*wp).m_f.t / (*wp).m_I; // sum(Torque)/I = alpha

		// update velcoities 1/2 step
		(*wp).m_v.x = (*wp).m_v.x + .5*TS*ax;
		(*wp).m_v.y = (*wp).m_v.y + .5*TS*ay;
		(*wp).m_v.th = (*wp).m_v.th + .5*TS*alpha;

		// update position
		posx = ( (*wp).m_x.x + TS * ((*wp).m_v.x) );
		posy = ( (*wp).m_x.y + TS * ((*wp).m_v.y) );
		theta = ( (*wp).m_x.th + TS * ((*wp).m_v.th) );

		// periodic boundary conditions
		while( posx > m_L) posx = posx - m_L;
		while( posx < 0) posx = posx + m_L;
		while( posy > m_L) posy = posy - m_L;
		while( posy < 0) posy = posy + m_L;
		while( theta < 0) theta = theta + 2*M_PI;
		while( theta > 2*M_PI) theta = theta - 2*M_PI;

		(*wp).m_x.x = posx;
		(*wp).m_x.y = posx;
		(*wp).m_x.th = theta;
	}

	//Update forces
	calc_sumforces();

	//Second half of v step
	for(m = 0; m < m_N; m++){
		wp = &m_molecs[m];

		// get (d^2 _ / d _^2)
		ax = (*wp).m_f.fx / 18.0; // sum(F_x)/m = a_x
		ay = (*wp).m_f.fy / 18.0; // sum(F_y)/m = a_y
		alpha = (*wp).m_f.t / (*wp).m_I; // sum(Torque)/I = alpha

		// update velcoities 1/2 step
		(*wp).m_v.x = (*wp).m_v.x + .5*TS*ax;
		(*wp).m_v.y = (*wp).m_v.y + .5*TS*ay;
		(*wp).m_v.th = (*wp).m_v.th + .5*TS*alpha;
	}

	//total momentum to zero
	zero_total_momentum();
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
				m_neigh[m*m_N + m_num_neigh[m]] = op; //set m's next neighbor pointer to wo
				m_neigh[o*m_N + m_num_neigh[o]] = wp; //set o's next neighbor pointer to wp


				m_num_neigh[m] += 1; //update the number of neighbors
				m_num_neigh[o] += 1;
			}
		}
	}

}

// Calculates the forces between every pair of molecules withing r_c of each other
// Puts the sum of fx, fy, t into the Force (m_f) of each molecule.
double MD::calc_sumforces(){
	int m, n;
	Water* wp;
	Water* np;
	Force* wpf_p;
	Force* npf_p;
	Point* cm1p_p;
	Point* cm2p_p;

	//first clear molecules' forces
	for( m = 0; m < m_N; m++){
		wp = &m_molecs[m];

		(*wp).m_f.fx = 0;
		(*wp).m_f.fy = 0;
		(*wp).m_f.t = 0;
	}

	//sum forces for all molecules
	for( m = 0; m < m_N; m++){
		wp = &m_molecs[m];

		// for all neighbors
		for ( n = 0; n < m_num_neigh[m]; n++){
			np = &(*m_neigh[m*m_N + n]); //the mth molecule's nth neighbor
			if ((*np).m_moln < m) continue; //only calculate forces once

			wpf_p = &((*wp).m_f); // pointer to force on particle 1
			npf_p = &((*np).m_f); // pointer to force on particle 2
			cm1p_p = &((*wp).m_x); // pointer to CM location of particle 1
			cm2p_p = &((*np).m_x); // pointer to CM location of particle 2

			lj_force( wpf_p, npf_p, &((*wp).m_ol), &((*np).m_ol), cm1p_p, cm2p_p );// LJ o1 and o2
			c_force( wpf_p, npf_p, QO, QH, &((*wp).m_ol), &((*np).m_al), cm1p_p, cm2p_p );// C o1 and h2a
			c_force( wpf_p, npf_p, QO, QH, &((*wp).m_ol), &((*np).m_bl), cm1p_p, cm2p_p );// C o1 and h2b
			c_force( wpf_p, npf_p, QO, QH, &((*wp).m_al), &((*np).m_ol), cm1p_p, cm2p_p );// C h1a and o2
			c_force( wpf_p, npf_p, QO, QH, &((*wp).m_bl), &((*np).m_ol), cm1p_p, cm2p_p );// C h1b and o2
			c_force( wpf_p, npf_p, QO, QO, &((*wp).m_ol), &((*np).m_ol), cm1p_p, cm2p_p );// C o1 and o2
			c_force( wpf_p, npf_p, QH, QH, &((*wp).m_al), &((*np).m_al), cm1p_p, cm2p_p );// C h1a and h2a
			c_force( wpf_p, npf_p, QH, QH, &((*wp).m_al), &((*np).m_bl), cm1p_p, cm2p_p );// C h1a and h2b
			c_force( wpf_p, npf_p, QO, QH, &((*wp).m_bl), &((*np).m_al), cm1p_p, cm2p_p );// C h1b and h2a
			c_force( wpf_p, npf_p, QO, QH, &((*wp).m_bl), &((*np).m_bl), cm1p_p, cm2p_p );// C h1b and h2b 
		}
	}
	// at this point the molecule's Force (m_f) contains a sum of (fx, fy, t) for all molecules
}

double MD::calc_PE(){
	int m, n;
	Water* wp;
	Water* np;
	double pe_tot = 0;

	for( m = 0; m < m_N; m++){ //for all molecules
		wp = &m_molecs[m];

		for( n = 0; n < m_num_neigh[m]; n++){ //for all neighbors
			np = &(*m_neigh[m*m_N + n]); //the mth molecule's nth neighbor
			if ((*np).m_moln < m) continue; //only calculate potential between two molecs once

			pe_tot += lj_pe( &((*wp).m_ol), &((*np).m_ol) );// LJ PE between o1 and o2
			pe_tot += c_pe( QO, QH, &((*wp).m_ol), &((*np).m_al) );// C PE between o1 and h2a
			pe_tot += c_pe( QO, QH, &((*wp).m_ol), &((*np).m_bl) );// C PE between o1 and h2b
			pe_tot += c_pe( QO, QH, &((*wp).m_al), &((*np).m_ol) );// C PE between h1a and o2
			pe_tot += c_pe( QO, QH, &((*wp).m_bl), &((*np).m_ol) );// C PE between h1b and o2
			pe_tot += c_pe( QO, QO, &((*wp).m_ol), &((*np).m_ol) );// C PE between o1 and o2
			pe_tot += c_pe( QH, QH, &((*wp).m_al), &((*np).m_al) );// C PE between h1a and h2a
			pe_tot += c_pe( QH, QH, &((*wp).m_al), &((*np).m_bl) );// C PE between h1a and h2b
			pe_tot += c_pe( QO, QH, &((*wp).m_bl), &((*np).m_al) );// C PE between h1b and h2a
			pe_tot += c_pe( QO, QH, &((*wp).m_bl), &((*np).m_bl) );// C PE between h1b and h2b 
		}
	}

	return pe_tot;
}

// Calcs LJ force between o1 and o2 and adds it to f1 and f2 force structs
void MD::lj_force(Force* f1_p, Force* f2_p, Point* o1_p, Point* o2_p, Point* cm1_p, Point* cm2_p){

	// Work out normalized direction vectors
	Point r12 = {};
	Point rc1 = {};
	Point rc2 = {};
	Point f12 = {};
	Point f21 = {};
	double dist12;
	double mag_f;

	(*o2_p).sub(*o1_p, &r12, m_L); // make r12 (not normalized)
	(*o1_p).sub(*cm1_p, &rc1, m_L); // make rc1
	(*o2_p).sub(*cm2_p, &rc2, m_L); // make rc2

	dist12 = r12.get_mag();
	r12.normalize(); //normalize 
	mag_f = 24*(EPS/SIGMA)*( 2*pow((SIGMA/dist12), 13) - pow((SIGMA/dist12), 7) ); //magnitude of force
	r12.mult(mag_f, &f12); // make f12
	f12.mult(-1.0, &f21); // make f21

	// force on o1 from o2
	(*f1_p).fx += f12.x; // fx in positive (o1 to o2) direction
	(*f1_p).fy += f12.y; // fy in positive (o1 to o2) direction
	(*f1_p).t += rc1.zcross(f12);

	// force on o2 from o1
	(*f2_p).fx += f21.x; // force in (o2 to o1) direction
	(*f2_p).fy += f21.y;
	(*f2_p).t += rc2.zcross(f21);

}

// Calculates the Coulomb force between particle 1 and particle 2 and adds it to f1 and f2 force structs
void MD::c_force(Force* f1_p, Force* f2_p, double q1, double q2, Point* loc1_p, Point* loc2_p, Point* cm1_p, Point* cm2_p){
	Point r12 = {};
	Point rc1 = {};
	Point rc2 = {};
	Point f12 = {};
	Point f21 = {};
	double dist12;
	double mag_f;

	(*loc2_p).sub(*loc1_p, &r12, m_L); //make r12
	(*loc1_p).sub(*cm1_p, &rc1, m_L); //make rc1
	(*loc2_p).sub(*cm2_p, &rc2, m_L); //make rc2

	dist12 = r12.get_mag();
	r12.normalize(); //normalize
	mag_f = q1*q2/(dist12*dist12); //magnitude of force
	r12.mult(mag_f, &f12); //make f12 = |f|*r_hat12 = force on 2
	f12.mult(-1.0, &f21); //make f21 = |f|*r_hat21 = force on 1

	// force on 1 from 2
	(*f1_p).fx += f21.x;
	(*f1_p).fy += f21.y;
	(*f1_p).t += rc1.zcross(f21);

	// force on 2 from 1
	(*f2_p).fx += f12.x;
	(*f2_p).fy += f12.y;
	(*f2_p).t += rc2.zcross(f12);
}

// calcs and returns the LJ Potential between o1 and o2
double MD::lj_pe(Point* o1_p, Point* o2_p){

	// Work out normalized direction vectors
	Point r12 = {};
	double dist12;

	(*o2_p).sub(*o1_p, &r12, m_L); // make r12
	dist12 = r12.get_mag();
	return 4*EPS*( pow((SIGMA/dist12), 12) - pow((SIGMA/dist12), 6) ); //magnitude of PE
}

// calcs and returns the Coloumb Potential between particle 1 and particle 2
double MD::c_pe(double q1, double q2, Point* loc1_p, Point* loc2_p){
	Point r12 = {};
	double dist12;

	(*loc2_p).sub(*loc1_p, &r12, m_L); //make r12
	dist12 = r12.get_mag();
	return q1*q2/dist12; //magnitude of force
}

// Calculate KE by summing the KEx, KEy, KEw for each molecule
double MD::calc_KE(){
	double ret = 0;
	Water* wp = &m_molecs[0];
	for(int m = 0; m < m_N; m++){
		wp = &m_molecs[m];
		ret += .5*18*(wp->m_v.x)*(wp->m_v.x); // (1/2)m vx^2
		ret += .5*18*(wp->m_v.y)*(wp->m_v.y); // (1/2)m vy^2
		ret += .5*(wp->m_I)*(wp->m_v.th)*(wp->m_v.th); // (1/2)I w^2
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