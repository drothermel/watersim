/*!
 @file Water.cpp
 @author Danielle Rothermel
 @date 4/27/14

 @brief Water contains all of the information about a single water molecule.

 Note that I'm generally using the z 

 */

#include "Water.h"
using namespace std;

Water::Water(double xc, double yc, double theta_rad){

	m_x.x = xc;
	m_x.y = yc;
	m_x.th = theta_rad;

	m_docm = (1.0/9.0) * cos(HOH_A / 2.0);
	m_I = 2 - (2/9)*cos(HOH_A/2.0)*cos(HOH_A/2.0);

	m_ol.th = 0.0;
	m_al.th = 0.0;
	m_bl.th = 0.0;

	update_oab();
}

void Water::init(){
	m_docm = (1.0/9.0) * cos(HOH_A / 2.0);
	m_I = 2 - (2/9)*cos(HOH_A/2.0)*cos(HOH_A/2.0);

	m_ol.th = 0.0;
	m_al.th = 0.0;
	m_bl.th = 0.0;
}

void Water::update_pos(Point nx){

	m_x.x = nx.x;
	m_x.y = nx.y;
	m_x.th = nx.th;

	update_oab();
}

void Water::update_vel(Point nv){
	m_v.x = nv.x;
	m_v.y = nv.y;
	m_v.th = nv.th;
}

void Water::update_oab(){
	m_ol.x = m_x.x + m_docm * cos(m_x.th);
	m_ol.y = m_x.y + m_docm * sin(m_x.th);

	m_al.x = m_ol.x - cos(HOH_A/2.0 - m_x.th);
	m_al.y = m_ol.y + sin(HOH_A/2.0 - m_x.th);

	m_bl.x = m_ol.x - cos(HOH_A/2.0 + m_x.th);
	m_bl.y = m_ol.y - sin(HOH_A/2.0 + m_x.th);
}

