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
#include "ann_nvt_lv.hpp"
using namespace std;

void ann_nvt_lv::initialize_lv(ann_data *d) {
	d->gzi_x = new double[d->natom];
	d->gzi_y = new double[d->natom];
	d->gzi_z = new double[d->natom];
	
	d->eta_x = new double[d->natom];
	d->eta_y = new double[d->natom];
	d->eta_z = new double[d->natom];
	
	d->eg_x = new double[d->natom];
	d->eg_y = new double[d->natom];
	d->eg_z = new double[d->natom];
	
	for (int i=0; i<d->natom; i++) {
		d->gzi_x[i] = 0.;
		d->gzi_y[i] = 0.;
		d->gzi_z[i] = 0.;
		d->eta_x[i] = 0.;
		d->eta_y[i] = 0.;
		d->eta_z[i] = 0.;
		d->eg_x[i] = 0.;
		d->eg_y[i] = 0.;
		d->eg_z[i] = 0.;
	}
}

void ann_nvt_lv::generate_gausian_variable(ann_data *d){
	double num1, num2, x1, x2;
	double randm2 = (double)RAND_MAX + 2.;
	
	for (int i=0; i<d->natom; i++) {
		for (int j=0; j<3; j++) {
			num1 = ((double)rand()+1.) / randm2;
			num2 = ((double)rand()+1.) / randm2;
			
			//x1 = sqrt(-2.0 * log(num1)) * cos(2.*PI*num2);
			//x2 = sqrt(-2.0 * log(num1)) * sin(2.*PI*num2);       
			//d->box_m1[j] = 0.0 + 1.0 * x1;     //average = 0.0, SD = 1.0
			//d->box_m2[j] = 0.0 + 1.0 * x2;
			
			d->box_m1[j] = sqrt(-2.0 * log(num1)) * cos(PI2 * num2);
			d->box_m2[j] = sqrt(-2.0 * log(num1)) * sin(PI2 * num2);
		}
        
		d->gzi_x[i] = d->box_m1[0];
		d->gzi_y[i] = d->box_m2[0];
		d->gzi_z[i] = d->box_m1[1];
		d->eta_x[i] = d->box_m2[1];
		d->eta_y[i] = d->box_m1[2];
		d->eta_z[i] = d->box_m2[2];
                
		d->eg_x[i] = d->mass_const6[i] * (0.5*d->gzi_x[i] + d->eta_x[i]/s3);
		d->eg_y[i] = d->mass_const6[i] * (0.5*d->gzi_y[i] + d->eta_y[i]/s3);
		d->eg_z[i] = d->mass_const6[i] * (0.5*d->gzi_z[i] + d->eta_z[i]/s3);       
                //printf("gzi: %20.10f%20.10f%20.10f\n", d->gzi_x[i], d->gzi_y[i], d->gzi_z[i]);
	}
}

void ann_nvt_lv::update_velocity(ann_data *d){
	for (int i=0; i<d->natom; i++) {
		d->vx[i] += d->mass_const2[i] * d->fcx[i] - d->time_h_fric * d->vx[i] + d->mass_const3[i] * d->gzi_x[i] - 
					d->mass_const4[i] * d->fcx[i] + d->md_const4_fric * d->vx[i] - d->eg_x[i];
		
		d->vy[i] += d->mass_const2[i] * d->fcy[i] - d->time_h_fric * d->vy[i] + d->mass_const3[i] * d->gzi_y[i] - 
					d->mass_const4[i] * d->fcy[i] + d->md_const4_fric * d->vy[i] - d->eg_y[i];
		
		d->vz[i] += d->mass_const2[i] * d->fcz[i] - d->time_h_fric * d->vz[i] + d->mass_const3[i] * d->gzi_z[i] - 
					d->mass_const4[i] * d->fcz[i] + d->md_const4_fric * d->vz[i] - d->eg_z[i];
	}
}

void ann_nvt_lv::update_position(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		//c1 = d->mass_const11[i] * d->eta[i];
		d->cx[i] += d->md_time * d->vx[i] + d->mass_const5[i] * d->eta_x[i];
		d->cy[i] += d->md_time * d->vy[i] + d->mass_const5[i] * d->eta_y[i];
		d->cz[i] += d->md_time * d->vz[i] + d->mass_const5[i] * d->eta_z[i];
	}
}

void ann_nvt_lv::copy_prev_force(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->pfcx[i] = d->fcx[i];
		d->pfcy[i] = d->fcy[i];
		d->pfcz[i] = d->fcz[i];
	}
}

void ann_nvt_lv::calc_current_temperature(ann_data *d) {
	int i;
	
	d->c_temp = 0;
	for (i=0; i<d->natom; i++) {
		d->c_temp += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->c_temp = d->c_temp * 0.5 * d->c_unit2 /(1.5  * (double)d->natom * kB);
	
	if (d->md_dtemp != 0.) {
		d->md_temp += d->md_dtemp;
		
		d->lv_sigma = sqrt(2. * kB * d->md_temp * d->md_friction * d->c_unit1);
		d->md_const5 = pow(d->md_time, 1.5) * d->lv_sigma / (2.*sqrt(3.));
		d->md_const6 = 1./4. * pow(d->md_time, 1.5) * d->md_friction * d->lv_sigma;
		
		for (i=0; i<d->natom; i++) {
			d->mass_const3[i] = 0.5 * sqrt(d->md_time) * d->lv_sigma / sqrt(d->mass[i]);
			d->mass_const5[i] = d->md_const5 / sqrt(d->mass[i]);
			d->mass_const6[i] = d->md_const6 / sqrt(d->mass[i]);
		}
	}
}