#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end
#include <limits>
#include <random>
#include "ann_data.hpp"
#include "ann_coordinate.hpp"
using namespace std;

void ann_coordinate::frac_to_cart(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->cx[i] = d->x[i]*d->lv[0][0] + d->y[i]*d->lv[1][0] + d->z[i]*d->lv[2][0];
		d->cy[i] = d->x[i]*d->lv[0][1] + d->y[i]*d->lv[1][1] + d->z[i]*d->lv[2][1];
		d->cz[i] = d->x[i]*d->lv[0][2] + d->y[i]*d->lv[1][2] + d->z[i]*d->lv[2][2];
		//printf("x y z: %15.7f%15.7f%15.7f\n", d->cx[i], d->cy[i], d->cz[i]);
	}
}

void ann_coordinate::cart_to_frac(ann_data *d) {
	double xx, yy, zz, xd, yd, zd;

	double th = 1e-5;
	double lx1 = -th / d->lt[0];
	double lx2 = 1. - th / d->lt[0];
	double ly1 = -th / d->lt[1];
	double ly2 = 1. - th / d->lt[1];
	double lz1 = -th / d->lt[2];
	double lz2 = 1. - th / d->lt[2];
	
	for (int i=0; i<d->natom; i++) {
		d->x[i] = d->cx[i]*d->ilv[0][0] + d->cy[i]*d->ilv[1][0] + d->cz[i]*d->ilv[2][0];
		d->y[i] = d->cx[i]*d->ilv[0][1] + d->cy[i]*d->ilv[1][1] + d->cz[i]*d->ilv[2][1];
		d->z[i] = d->cx[i]*d->ilv[0][2] + d->cy[i]*d->ilv[1][2] + d->cz[i]*d->ilv[2][2];
		
		xx = modf(d->x[i], &xd);
		yy = modf(d->y[i], &yd);
		zz = modf(d->z[i], &zd);
		
		d->x[i] = xx;
		d->y[i] = yy;
		d->z[i] = zz;

		/*
		if (d->x[i] < lx1) {
			d->x[i] += 1.;
		} else if (d->x[i] > lx2) {
			d->x[i] -= 1.;
		}
		
		if (d->y[i] < ly1) {
			d->y[i] += 1.;
		} else if (d->y[i] > ly2) {
			d->y[i] -= 1.;
		}
		
		if (d->z[i] < lz1) {
			d->z[i] += 1.;
		} else if (d->z[i] > lz2) {
			d->z[i] -= 1.;
		}

		if (d->x[i] <= lx1) d->x[i] = 0.;
		if (d->y[i] <= ly1) d->y[i] = 0.;
		if (d->z[i] <= lz1) d->z[i] = 0.;
		if (lx2 <= d->x[i]) d->x[i] = 0.;
		if (ly2 <= d->y[i]) d->y[i] = 0.;
		if (lz2 <= d->z[i]) d->z[i] = 0.;
		*/
		
		if (d->x[i] < 0.) {
			d->x[i] += 1.;
		} else if (d->x[i] >= 1.) {
			d->x[i] -= 1.;
		}
		
		if (d->y[i] < 0.) {
			d->y[i] += 1.;
		} else if (d->y[i] >= 1.) {
			d->y[i] -= 1.;
		}
		
		if (d->z[i] < 0.) {
			d->z[i] += 1.;
		} else if (d->z[i] >= 1.) {
			d->z[i] -= 1.;
		}

		if (d->x[i] < lx1) d->x[i] = 0.;
		if (d->y[i] < ly1) d->y[i] = 0.;
		if (d->z[i] < lz1) d->z[i] = 0.;

		if (lx2 < d->x[i]) d->x[i] = 0.;
		if (ly2 < d->y[i]) d->y[i] = 0.;
		if (lz2 < d->z[i]) d->z[i] = 0.;
	}
}

void ann_coordinate::calc_lattice_const(ann_data *d) {
	double lt02, lt12, lt22;
	double t1, t2, t3, a, b,   b1, b2, b3;

	d->lt[0] = sqrt(d->lv[0][0] * d->lv[0][0] + d->lv[0][1] * d->lv[0][1] + d->lv[0][2] * d->lv[0][2]);
	d->lt[1] = sqrt(d->lv[1][0] * d->lv[1][0] + d->lv[1][1] * d->lv[1][1] + d->lv[1][2] * d->lv[1][2]);
	d->lt[2] = sqrt(d->lv[2][0] * d->lv[2][0] + d->lv[2][1] * d->lv[2][1] + d->lv[2][2] * d->lv[2][2]);
	//cout << d->lt[0] <<" " << d->lt[1] << " " << d->lt[2] << endl;

	lt02 = d->lt[0] * d->lt[0];
	lt12 = d->lt[1] * d->lt[1];
	lt22 = d->lt[2] * d->lt[2];
	//cout << lt02 <<" " << lt12 << " " << lt22 << endl;

	t1 = d->lv[0][0] - d->lv[1][0];
	t2 = d->lv[0][1] - d->lv[1][1];
	t3 = d->lv[0][2] - d->lv[1][2];
	a = sqrt(t1*t1 + t2*t2 + t3*t3);
	b = (lt02 + lt12 - a*a) / (2 * d->lt[0] * d->lt[1]);
	b3 = (lt02 + lt12 - a*a) / (2 * d->lt[0] * d->lt[1]); 
	d->lt[5] = rad_to_deg(acos(b));

	t1 = d->lv[0][0] - d->lv[2][0];
	t2 = d->lv[0][1] - d->lv[2][1];
	t3 = d->lv[0][2] - d->lv[2][2];
	a = sqrt(t1*t1 + t2*t2 + t3*t3);
	b = (lt02 + lt22 - a*a) / (2 * d->lt[0] * d->lt[2]);
	b2 = (lt02 + lt22 - a*a) / (2 * d->lt[0] * d->lt[2]);   
	d->lt[4] = rad_to_deg(acos(b));

	t1 = -d->lv[1][0] + d->lv[2][0];
	t2 = -d->lv[1][1] + d->lv[2][1];
	t3 = -d->lv[1][2] + d->lv[2][2];
	a = sqrt(t1*t1 + t2*t2 + t3*t3);
	b = (lt12 + lt22 - a*a) / (2 * d->lt[1] * d->lt[2]);
	b1 = (lt12 + lt22 - a*a) / (2 * d->lt[1] * d->lt[2]);   
	d->lt[3] = rad_to_deg(acos(b));
	
	//d->volume = d->lt[0] * d->lt[1] * d->lt[2] * sqrt(1 - b1 * b1 -b2 * b2 -b3 * b3 +2.0 * b1 * b2 * b3); 
	//if (d->rank == 0 ) printf("primitive cell volume = %-15.7f\n", d->volume);

	d->volume = d->lv[0][0]*d->lv[1][1]*d->lv[2][2] + d->lv[0][2]*d->lv[1][0]*d->lv[2][1] + d->lv[0][1]*d->lv[1][2]*d->lv[2][1]
				- (d->lv[0][2]*d->lv[1][1]*d->lv[2][1] + d->lv[0][1]*d->lv[1][0]*d->lv[2][2] + d->lv[0][0]*d->lv[1][2]*d->lv[2][1]);
	//if (d->rank == 0 ) printf("cell volume = %-15.7f\n", d->volume);
}

void ann_coordinate::calc_ilv(ann_data *d) {
	double buf;
	double a[3][3];
	int i, j, k;

	for (i=0; i<3; i++){
		for (j=0; j<3; j++) {
			a[i][j] = d->lv[i][j];
		}
	}

	for (i=0; i<3; i++){
		for (j=0; j<3; j++) {
			d->ilv[i][j] = (i==j)?1.:0.;
		}
	}

	for (i=0; i<3; i++) {
		buf = 1./a[i][i];
		for (j=0; j<3; j++) {
			a[i][j] *= buf;
			d->ilv[i][j] *= buf;
		}

		for (j=0; j<3; j++) {
			if (i!=j) {
				buf = a[j][i];
				for (k=0; k<3; k++) {
					a[j][k] -= a[i][k] * buf;
					d->ilv[j][k] -= d->ilv[i][k] * buf;
				}
			}
		}
	}
}