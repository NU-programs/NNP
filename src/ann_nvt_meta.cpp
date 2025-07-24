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
#include "ann_nvt_meta.hpp"
using namespace std;

void ann_nvt_meta::copy_prev_position(ann_data *d){
	for (int i=0; i<d->natom; i++) {
		d->px[i] = d->cx[i];
		d->py[i] = d->cy[i];
		d->pz[i] = d->cz[i];
	}
}

void ann_nvt_meta::cutoff_particle_list(ann_data *d){
	double dx1, dx2, jx, jy, jz, jj;
	int n;

	dx1 = 0.11627 - 1e-3;
	dx2 = 0.44258 + 1e-3;
	//printf("dz = %15.7f\n", dz);
	
	n = 1;
	for (int i=0; i<d->natom; i++){
		if (dx1 <= d->x[i] && d->x[i] <= dx2) {
			d->cutoff_list[i] = n;
			for (int j=0; j<d->natom; j++){
				jx = d->x[j] - (1. - d->x[i]);
				jy = d->y[j] - d->y[i];
				jz = d->z[j] - d->z[i];
				jj = sqrt(jx*jx + jy*jy + jz*jz);
				if (jj < 1e-3){
					d->cutoff_list[j] = n;
				}
			}
			d->mt_atom = n;
			n += 1;
		}
		//printf("cutoff_list: %d  , %d\n", i, d->cutoff_list[i]);
	}
	//printf("metadynamics_atom: %d\n", d->mt_atom);
}

void ann_nvt_meta::initialize_lv(ann_data *d) {
	d->gzi_x = new double[d->natom];
	d->gzi_y = new double[d->natom];
	d->gzi_z = new double[d->natom];
	
	d->eta_x = new double[d->natom];
	d->eta_y = new double[d->natom];
	d->eta_z = new double[d->natom];
	
	d->eg_x = new double[d->natom];
	d->eg_y = new double[d->natom];
	d->eg_z = new double[d->natom];
	
	d->ds_dx = new double[d->natom];
	d->ds_dy = new double[d->natom];
	d->ds_dz = new double[d->natom];
	
	d->x_s_of_t = new double[d->md_step];
	d->y_s_of_t = new double[d->md_step];
	d->z_s_of_t = new double[d->md_step];
	
	d->cutoff_list = new int[d->natom];
	d->wh = new double[d->md_step];
	
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

		d->ds_dx[i] = 0.;
		d->ds_dy[i] = 0.;
		d->ds_dz[i] = 0.;
		
		d->cutoff_list[i] = 0;
	}
	
	for (int i=0; i<d->md_step; i++) {
		d->x_s_of_t[i] = 0.;
		d->y_s_of_t[i] = 0.;
		d->z_s_of_t[i] = 0.;
		d->wh[i] = 0.;
	}
	
	d->NG = 0;
	d->tmpV = 0;
	d->plus_Vg = 0.;
}

void ann_nvt_meta::generate_gaussian_variable(ann_data *d){
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

void ann_nvt_meta::update_velocity(ann_data *d){
	for (int i=0; i<d->natom; i++) {
		d->vx[i] += d->mass_const2[i] * d->fcx[i] - d->time_h_fric * d->vx[i] + d->mass_const3[i] * d->gzi_x[i] - 
					d->mass_const4[i] * d->fcx[i] + d->md_const4_fric * d->vx[i] - d->eg_x[i];
		
		d->vy[i] += d->mass_const2[i] * d->fcy[i] - d->time_h_fric * d->vy[i] + d->mass_const3[i] * d->gzi_y[i] - 
					d->mass_const4[i] * d->fcy[i] + d->md_const4_fric * d->vy[i] - d->eg_y[i];
		
		d->vz[i] += d->mass_const2[i] * d->fcz[i] - d->time_h_fric * d->vz[i] + d->mass_const3[i] * d->gzi_z[i] - 
					d->mass_const4[i] * d->fcz[i] + d->md_const4_fric * d->vz[i] - d->eg_z[i];
	}
}

void ann_nvt_meta::update_position(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->cx[i] += d->md_time * d->vx[i] + d->mass_const5[i] * d->eta_x[i];
		d->cy[i] += d->md_time * d->vy[i] + d->mass_const5[i] * d->eta_y[i];
		d->cz[i] += d->md_time * d->vz[i] + d->mass_const5[i] * d->eta_z[i];
		
		d->px[i] += d->md_time * d->vx[i] + d->mass_const5[i] * d->eta_x[i];
		d->py[i] += d->md_time * d->vy[i] + d->mass_const5[i] * d->eta_y[i];
		d->pz[i] += d->md_time * d->vz[i] + d->mass_const5[i] * d->eta_z[i];
	}
}

void ann_nvt_meta::copy_prev_force(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->pfcx[i] = d->fcx[i];
		d->pfcy[i] = d->fcy[i];
		d->pfcz[i] = d->fcz[i];
	}
}

void ann_nvt_meta::calc_current_temperature(ann_data *d) {
	d->c_temp = 0;
	for (int i=0; i<d->natom; i++) {
			d->c_temp += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->c_temp = d->c_temp * 0.5 * d->c_unit2 /(1.5  * (double)d->natom * kB);

	if (d->md_dtemp != 0.) {
		d->md_temp += d->md_dtemp;
		
		d->lv_sigma = sqrt(2. * kB * d->md_temp * d->md_friction * d->c_unit1);
		d->md_const5 = pow(d->md_time, 1.5) * d->lv_sigma / (2.*sqrt(3.));
		d->md_const6 = 1./4. * pow(d->md_time, 1.5) * d->md_friction * d->lv_sigma;
		
		for (int i=0; i<d->natom; i++) {
				d->mass_const3[i] = 0.5 * sqrt(d->md_time) * d->lv_sigma / sqrt(d->mass[i]);
				d->mass_const5[i] = d->md_const5 / sqrt(d->mass[i]);
				d->mass_const6[i] = d->md_const6 / sqrt(d->mass[i]);
		}
	}
}

void ann_nvt_meta::set_parameter(ann_data *d){
	d->sy = new double[d->mt_atom+1];
	d->sz = new double[d->mt_atom+1];
	
	for (int i=1; i<d->mt_atom+1; i++){
		d->sy[i] = 0.;
		d->sz[i] = 0.;
	}
}

void ann_nvt_meta::make_CV(ann_data *d){
	double temp_dy[2], temp_dz[2];
	int n;
	double ty, tz;	
	
	d->yg = 0.;
	d->zg = 0.;
	
	for (int i=1; i<d->mt_atom+1; i++){
		n = 0;
		for (int j=0; j<d->natom; j++){
			if (d->cutoff_list[j] == i){
				temp_dy[n] = d->py[j];
				temp_dz[n] = d->pz[j];
				n += 1;
			}
		}
		ty = temp_dy[0] - temp_dy[1];
		tz = temp_dz[0] - temp_dz[1];
		
		if (ty > 0.5*d->lv[1][1]) {
			ty -= d->lv[1][1];
		} else if (ty < -0.5*d->lv[1][1]) {
			ty += d->lv[1][1];
		}
		
		if (tz > 0.5*d->lv[2][2]) {
			tz -= d->lv[2][2];
		} else if (tz < -0.5*d->lv[2][2]) {
			tz += d->lv[2][2];
		}
		
		d->sy[i] = ty;
		d->sz[i] = tz;
		d->yg += ty*ty;
		d->zg += tz*tz;
	}
	
	d->yg /= d->mt_atom;
	d->zg /= d->mt_atom;
	//cout << "d->zg = " << d->yg << endl;
}

void ann_nvt_meta::partial_diff(ann_data *d){
	double half_x;
	half_x = 0.5 * d->lv[0][0];
	
	for (int i=1; i<d->mt_atom+1; i++) {
		for (int j=0; j<d->natom; j++) {
			if (d->cutoff_list[j] == i) {
				if (d->cx[j] < half_x) {
					d->ds_dy[j] = (2.* d->sy[i]) / d->mt_atom;
					d->ds_dz[j] = (2.* d->sz[i]) / d->mt_atom;
				} else {
					d->ds_dy[j] = (-2.* d->sy[i]) / d->mt_atom;
					d->ds_dz[j] = (-2.* d->sz[i]) / d->mt_atom;
				}
				//cout << "ds_dz=  " << d->ds_dz[j] << endl;
			}
		}
	}
}

void ann_nvt_meta::save_s_of_t(ann_data *d){
	d->NG += 1;
	//d->x_s_of_t[d->NG] = d->xg;
	d->y_s_of_t[d->NG] = d->yg;
	d->z_s_of_t[d->NG] = d->zg;
	d->tmpV = d->plus_Vg;
	
	if (d->NG == 1) {
		d->wh[d->NG] = d->gau_h;
	} else {
		d->wh[d->NG] = d->gau_h * exp(-d->tmpV/d->delT);
	}
	
	//printf("s_of_t  = %20.10f%20.10f%20.10f\n", d->x_s_of_t[d->NG], d->y_s_of_t[d->NG], d->z_s_of_t[d->NG]);
	printf("WTmeta_w = %20.10f\n", d->wh[d->NG]);
}

void ann_nvt_meta::calc_vias_potential(ann_data *d){
	double dx, dy, dz, dd;
	d->dVg_ds_y = 0.;
	d->dVg_ds_z = 0.;
	d->plus_Vg = 0.;
	dd = d->gau_w * d->gau_w;
	
	for (int i=1; i<d->NG+1; i++){
		dy = d->yg - d->y_s_of_t[i];
		dz = d->zg - d->z_s_of_t[i];
		
		d->gauss = d->wh[d->NG] * exp(-(dy*dy + dz*dz)/(2. *dd));
		d->plus_Vg += d->gauss;
		
		d->dVg_ds_y += -d->gauss * (dy / dd);
		d->dVg_ds_z += -d->gauss * (dz / dd);
		//printf("gauss, dVg_x, dVg_y, dVg_z  = %20.10f%20.10f%20.10f%20.10f\n", d->gauss, d->dVg_ds_x, d->dVg_ds_y, d->dVg_ds_z);
	}

}

void ann_nvt_meta::force_plus_vias(ann_data *d){
	for (int i=0; i<d->natom; i++){
		if (d->cutoff_list[i] >= 1 ) {
			//d->fcx[i] += -d->dVg_ds_x * d->ds_dx[i];
			d->fcy[i] += -d->dVg_ds_y * d->ds_dy[i];
			d->fcz[i] += -d->dVg_ds_z * d->ds_dz[i];
			//printf("fcx, fcy, fcz  = %20.10f%20.10f%20.10f\n", d->fcx[i], d->fcy[i], d->fcz[i]);
		}
	}
}