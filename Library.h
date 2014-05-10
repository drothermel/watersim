#ifndef LIBRARY_H
#define LIBRARY_H
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

#define HOH_A (109.47*M_PI/180.0) // Degrees
#define ROH 1.0 // Angstrom
#define EPS 0.650 // KJ/mol
#define SIGMA 3.166 // Angstrom
#define QO -0.820// Charge of Oxygen
#define QH 0.410// Charge of Hydrogen
#define TS 1 // timestep fs
#define NEIGHTS 20 // time between updating neighbors fs
#define COORDTS 500 // time between storing coords fs
#define ENERTS 100 // time between storing energy fs
#define RC 9.0 // cutoff distance in angstroms

struct Force{
	double fx, fy, t;
};

struct Point{
	double x, y, th;

	void normalize(){
		double mag = sqrt(x*x+y*y);
		x = x/mag;
		y = y/mag;
	}

	double get_mag(){
		return sqrt(x*x+y*y);
	}

	void add(Point addthis, Point* put_p){
		(*put_p).x = x+addthis.x;
		(*put_p).y = y+addthis.y;
	}

	//Includes boundary conditions
	void sub(Point subthis, Point* put_p, double L){
		double dx = x - subthis.x;
		double dy = y - subthis.y;
		double halfL = L/2.0;
		if(dx > halfL){
			(*put_p).x = dx - L;
		} else if(dx < -halfL){ 
			(*put_p).x = dx + L;
		} else{
			(*put_p).x = dx;
		}

		if(dy > halfL){
			(*put_p).y = dy - L;
		} else if(dy < -halfL){
			(*put_p).y = dy + L;
		} else{
			(*put_p).y = dy;
		}
	}

	void mult(double scalar, Point* put_p){
		(*put_p).x = x*scalar;
		(*put_p).y = y*scalar;
	}

	double dot(Point dotwthis){
		return x*dotwthis.x + y*dotwthis.y;
	}

	double zcross(Point crosswthis){
		return x*crosswthis.y - y*crosswthis.x;
	}

	void print(){
		printf("Point: %.2f, %.2f, %.2f\n", x, y, th);
	}
};

#endif //LIBRARY_H