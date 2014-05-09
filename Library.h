#ifndef LIBRARY_H
#define LIBRARY_H
#include <cmath>
#include <cstdio>
#include <Cstdlib>

#define HOH_A (109.47*M_PI/180) // Degrees
#define ROH 1.0 // Angstrom
#define EPS 0.650 // KJ/mol
#define SIGMA 3.166 // Angstrom
#define TS 1 // timestep fs
#define NEIGHTS 20 // time between updating neighbors fs
#define COORDTS 500 // time between storing coords fs
#define ENERTS 100 // time between storing energy fs
#define RC 9 // cutoff distance in angstroms

struct Point{
	double x, y, th;

	Point get_dir(){
		double mag = sqrt(x*x+y*y);
		Point p = {
			x/mag,
			y/mag,
			0.0
		};
		return p;
	}

	double get_mag(){
		return sqrt(x*x+y*y);
	}

	Point add(Point addthis){
		Point p = {
			x+addthis.x,
			y+addthis.y,
			0.0
		};
		return p;
	}

	Point sub(Point subthis){
		Point p = {
			x-subthis.x,
			y-subthis.y,
			0.0
		};
		return p;
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