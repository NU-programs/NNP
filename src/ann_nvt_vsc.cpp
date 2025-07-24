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
#include "ann_nvt_vsc.hpp"
using namespace std;

void ann_nvt_vsc::copy_prev_force(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->pfcx[i] = d->fcx[i];
		d->pfcy[i] = d->fcy[i];
		d->pfcz[i] = d->fcz[i];
	}
}

void ann_nvt_vsc::update_position(ann_data *d) {
	//cout << "mass const =" << d->mass_const1[0] << endl; 
	int i;
	for (i=0; i<d->natom; i++) {
		//printf("x y z: %15.7f%15.7f%15.7f\n", d->vx[i], d->vy[i], d->vz[i]);
		d->cx[i] = d->cx[i] + d->vx[i] * d->md_time + d->mass_const1[i] * d->fcx[i];
		d->cy[i] = d->cy[i] + d->vy[i] * d->md_time + d->mass_const1[i] * d->fcy[i];
		d->cz[i] = d->cz[i] + d->vz[i] * d->md_time + d->mass_const1[i] * d->fcz[i];
	}
}

void ann_nvt_vsc::update_velocity(ann_data *d) {
	//cout << "mass const =" << d->mass_const2[0] << endl; 
	int i;
	for (i=0; i<d->natom; i++) {
		d->vx[i] = d->vx[i] + d->mass_const2[i] * (d->pfcx[i] + d->fcx[i]);
		d->vy[i] = d->vy[i] + d->mass_const2[i] * (d->pfcy[i] + d->fcy[i]);
		d->vz[i] = d->vz[i] + d->mass_const2[i] * (d->pfcz[i] + d->fcz[i]);
	}
}

void ann_nvt_vsc::velocity_scaling(ann_data *d) {
	int i;
	double cvx, cvy, cvz;
	
	d->ke = 0;
	for (i=0; i<d->natom; i++) {
		d->ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->ke *= 0.5 * d->c_unit2;
	d->c_temp = d->ke / (1.5 * (double)d->natom * kB);
	//printf("Before ke [eV], c_temp [K] = %20.10f%20.10f\n", d->ke, d->c_temp);
	
	if (d->md_temp != 0.) {
		if (d->c_temp > 1e-10) { 
			//d->md_vscale = sqrt((d->md_temp + (d->c_temp - d->md_temp)*1e-4) / d->c_temp);
			d->md_vscale = sqrt(d->md_temp / d->c_temp);
			//printf("md_vscale = %20.10f\n", d->md_vscale);
			
			for (i=0; i<d->natom; i++) {
				d->vx[i] *= d->md_vscale;
				d->vy[i] *= d->md_vscale;
				d->vz[i] *= d->md_vscale;
			}
	
			d->ke = 0;
			for (i=0; i<d->natom; i++) {
				d->ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
			}
			d->ke *= 0.5 * d->c_unit2;
			d->c_temp = d->ke / (1.5 * (double)d->natom * kB);
			//printf("After ke [eV], c_temp [K] = %20.10f%20.10f\n", d->ke, d->c_temp);
		}
	} else {
		for (i=0; i<d->natom; i++) {
			d->vx[i] = 0.;
			d->vy[i] = 0.;
			d->vz[i] = 0.;
		}
	}
}