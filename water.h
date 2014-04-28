#ifndef WATER_H
#define WATER_H

const float HOH_A = 109.47;
const float ROH = 1.0;

struct Point{
	float x;
	float y;
};

class Water{
	
	public:
		Water(float xc, float yc, float theta_rad);
		virtual ~Water();
		
		Point m_cl; //COM location
		Point m_ol; //Oxygen location
		Point m_al; //Ha location
		Point m_bl; //Hb location
		float m_th; //H20 orientation

		void update_pos(float disp, float angle, float dtheta);

	private:
		float m_docm; //Distance from oxygen and CM
		void update_oab();

};

#endif //WATER_H