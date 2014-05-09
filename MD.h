#ifndef MD_H
#define MD_H

#include "Library.h"
#include "Water.h"

class MD{
	
	public:
		MD(double box_length_angstroms, int num_molecules, double init_temp, int sim_length);
		virtual ~MD();

		int m_N; //num particles
		double m_L; //box length angstroms
		double m_T; //initial temp
		int m_simL; //number of fs to run the sim, equivalent to number of timesteps since 1fs = timestep

		void evolve();

	private:
		int m_step; //current step number
		Water* m_molecs; //an array of N Waters
		Point* m_xcoords; //an array of (sim_length/COORDTS) arrays of N Points (collapsed into one array)
		Point* m_vcoords; //an array of (sim_length/COORDTS) arrays of N Points (collapsed into one array)
		Water** m_neigh; //an NxN array of neighbor pointers to Water molecules in m_molecs array
		int* m_num_neigh; //an N array of number of neighbors for each molecule
		double* m_temps; //an array of sim_length/ENERTS temperatures
		double* m_KE; //an array of sim_length/ENERTS KEs
		double* m_PE; //an array of sim_length/ENERTS PEs

		void init_position();
		void init_velocity();
		void update_neighbors();
		void single_timestep();
		double calc_sumforces();
		double calc_sumtorques();
		double calc_KE();
		double calc_PE();
		double calc_T();
		void zero_total_momentum();
		void scale_init_temps();

};

#endif //MD_H
