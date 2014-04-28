#ifndef WATER_H
#define WATER_H

const double HOH_A = 109.47;
const double ROH = 1.0;

struct Point{
	double x;
	double y;
};

class Water{
	
	public:
		Water(double xc, double yc, double theta_rad);
		
		Point m_cl; //COM location
		Point m_ol; //Oxygen location
		Point m_al; //Ha location
		Point m_bl; //Hb location
		double m_th; //H20 orientation

		void update_pos(double disp, double angle, double dtheta);

	private:
		double m_docm; //Distance from oxygen and CM
		void update_oab();

};

#endif //WATER_H