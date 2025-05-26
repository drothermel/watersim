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
		Water** m_neigh; //an NxN array of neighbor pointers to Water molecules in m_molecs array
		int* m_num_neigh; //an N array of number of neighbors for each molecule
		bool m_zeroTM;

		void init_position();
		void init_pos_vel();
		void single_timestep();
		void update_neighbors();
                void calc_sumforces();
		double calc_PE();
		void lj_force(Force* f1_p, Force* f2_p, Point* o1_p, Point* o2_p, Point* cm1_p, Point* cm2_p);
		void c_force(Force* f1_p, Force* f2_p, double q1, double q2, Point* loc1_p, Point* loc2_p, Point* cm1_p, Point* cm2_p);
		double lj_pe(Point* o1_p, Point* o2_p);
		double c_pe(double q1, double q2, Point* loc1_p, Point* loc2_p);
		double calc_KE();
		double calc_T();
		void zero_total_momentum();
		void scale_init_temps();

};

#endif //MD_H
