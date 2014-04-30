#ifndef WATER_H
#define WATER_H

#include "Library.h"

class Water{
	
	public:
		Water() {}
		Water(double xc, double yc, double theta_rad);
		
		Point m_x; //COM location
		Point m_v; //Velocity of molecule

		Point m_ol; //Oxygen location (th = 0)
		Point m_al; //Ha location (th = 0)
		Point m_bl; //Hb location (th = 0)

		void update_pos(Point newpos);
		void update_vel(Point newvel);

	private:
		double m_docm; //Distance from oxygen and CM
		double m_I; //Moment of inertia
		void update_oab();

};

#endif //WATER_H
