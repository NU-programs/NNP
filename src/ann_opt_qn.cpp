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
#include "ann_opt.hpp"
using namespace std;

void ann_opt::opt_init_qn(ann_data *d){
	int i, j;
	//set matrix
	d->H = new double*[d->natom3];
	d->ykdxT_T = new double*[d->natom3];
	d->dxdxT = new double*[d->natom3];
	d->H1 = new double*[d->natom3];

	for (i=0; i<d->natom3; i++){
		d->H[i] = new double[d->natom3];
		d->ykdxT_T[i] = new double[d->natom3];
		d->dxdxT[i] = new double[d->natom3];
		d->H1[i] = new double[d->natom3];
	}

	for (i=0; i<d->natom3; i++){
		for (j=0; j<d->natom3; j++){
			d->ykdxT_T[i][j] = 0.0;
			d->dxdxT[i][j] = 0.0;
			d->H1[i][j] = 0.0;
			//initialize Hessian(inverse) = identity matrix
			if (i == j){
				d->H[i][j] = d->H_init;
			}else{
				d->H[i][j] = 0.0;
			}
		}
	}

	//set array
	d->dx = new double[d->natom3];
	d->yk = new double[d->natom3];

	for (i=0; i<d->natom3; i++){
		d->dx[i] = 0.0;
		d->yk[i] = 0.0;
	}

}

void ann_opt::update_position_qn(ann_data *d){
	int i, j;
	//calculaiton delta
	for (i=0; i<d->natom3; i++){
		d->dx[i] = 0.;
		for (j=0; j<d->natom; j++){
			d->dx[i] += (d->H[i][3*j] * d->fcx[j] + d->H[i][3*j+1] * d->fcy[j] + d->H[i][3*j+2] * d->fcz[j]);
		}
		d->dx[i] *= d->step_length;
	}
	//update position
	for (int i=0; i<d->natom; i++) {
		d->cx[i] += d->dx[3*i];
		d->cy[i] += d->dx[3*i+1];
		d->cz[i] += d->dx[3*i+2];
	}
}

void ann_opt::approximate_hessian_inverse(ann_data *d){
	int i, j, k;
	double ykTdx = 0.;
	double r_ykTdx;
	double term;

	for (i=0; i<d->natom; i++){
		d->yk[3*i] = -(d->fcx[i]-d->pfcx[i]);
		d->yk[3*i+1] = -(d->fcy[i]-d->pfcy[i]);
		d->yk[3*i+2] = -(d->fcz[i]- d->pfcz[i]);
	}

	for (i=0; i<d->natom3; i++){
		ykTdx += d->yk[i] * d->dx[i];
	}
	r_ykTdx = 1. / ykTdx;

	for (i=0; i<d->natom3; i++){
		for (j=0; j<d->natom3; j++){
			if (i==j){
				d->ykdxT_T[i][j] = 1.0 - d->yk[j] * d->dx[i] * r_ykTdx;
			}else {
				d->ykdxT_T[i][j] = -d->yk[j] * d->dx[i] * r_ykTdx;
			}
		}
	}

	for (i=0; i<d->natom3; i++){
		for (j=0; j<=i; j++){
			d->dxdxT[i][j] = d->dx[i] * d->dx[j] * r_ykTdx;
			d->dxdxT[j][i] = d->dxdxT[i][j];
		}
	}

	for (i=0; i<d->natom3; i++){
		for (j=0; j<d->natom3; j++) {
			term = 0;
			for (k=0; k<d->natom3; k++) {
				term += d->ykdxT_T[i][k]*d->H[j][k];
			}
			d->H1[i][j] = term;
		}
	}

	for (i=0; i<d->natom3; i++) {
		for (j=0; j<d->natom3; j++) {
			term = 0;
			for (k=0; k<d->natom3; k++) {
				term += d->H1[i][k]*d->ykdxT_T[j][k];
			}
			d->H[i][j] = term;
		}
	}

	for (i=0; i<d->natom3; i++){
		for(j=0; j<=i; j++){
			d->H[i][j] += d->dxdxT[i][j];
			d->H[j][i] = d->H[i][j];
		}
	}
}

void ann_opt::reset_hessian(ann_data *d){
	d->H_init *= 0.1;
	d->step_length = d->step_length_init;

	for (int i=0; i<d->natom3; i++) {
		d->H[i][i] = d->H_init;

		for (int j=0; j<i; j++) {
			d->H[i][j] = 0.0;
			d->H[j][i] = 0.0;
		}
	}
}
