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

void ann_opt::update_cg_direction(ann_data *d) {
	int i;
	double beta, b1, b2;
	
	beta = 0.;
	b1 = 0.;
	b2 = 0.;
	
	// Polak-Ribiere
	for (i=0; i<d->natom; i++) {
		b1 = b1 + d->fcx[i]*(d->fcx[i]-d->pfcx[i]) + d->fcy[i]*(d->fcy[i]-d->pfcy[i]) + d->fcz[i]*(d->fcz[i]-d->pfcz[i]);
		b2 = b2 + d->pfcx[i]*d->pfcx[i] + d->pfcy[i]*d->pfcy[i] + d->pfcz[i]*d->pfcz[i];
	}
	//printf("b1, b2 = %15.7f%15.7f\n", b1, b2);
	
	if (fabs(b2) > 1e-12) {
		beta = b1 / b2;
	
		for (i=0; i<d->natom; i++) {
			d->cdtx[i] = d->fcx[i] + beta * d->cdtx[i];
			d->cdty[i] = d->fcy[i] + beta * d->cdty[i];
			d->cdtz[i] = d->fcz[i] + beta * d->cdtz[i];
			
			//d->cx[i] += d->step_length * d->cdtx[i];
			//d->cy[i] += d->step_length * d->cdty[i];
			//d->cz[i] += d->step_length * d->cdtz[i];
		}
	} else {
		for (i=0; i<d->natom; i++) {
			d->cdtx[i] = d->fcx[i];
			d->cdty[i] = d->fcy[i];
			d->cdtz[i] = d->fcz[i];
			
			//d->cx[i] += d->step_length * d->cdtx[i];
			//d->cy[i] += d->step_length * d->cdty[i];
			//d->cz[i] += d->step_length * d->cdtz[i];
		}
	}
}

void ann_opt::update_position_cg(ann_data *d) {
	/*
	if (d->selective_opt == "yes") {
		for (int i=0; i<d->natom; i++) {
			if (fabs(d->fcx[i]) >= d->fc_diff) d->cx[i] += d->step_length * d->cdtx[i];
			if (fabs(d->fcy[i]) >= d->fc_diff) d->cy[i] += d->step_length * d->cdty[i];
			if (fabs(d->fcz[i]) >= d->fc_diff) d->cz[i] += d->step_length * d->cdtz[i];
		}
	} else {
		for (int i=0; i<d->natom; i++) {
			d->cx[i] += d->step_length * d->cdtx[i];
			d->cy[i] += d->step_length * d->cdty[i];
			d->cz[i] += d->step_length * d->cdtz[i];
		}
	}
	*/

	for (int i=0; i<d->natom; i++) {
		d->cx[i] += d->step_length * d->cdtx[i];
		d->cy[i] += d->step_length * d->cdty[i];
		d->cz[i] += d->step_length * d->cdtz[i];
	}
}

void ann_opt::update_cg_direction_h(ann_data *d) {   ///2023.1.31
	int i,j;
	double beta, b1, b2;
	
	beta = 0.;
	b1 = 0.;
	b2 = 0.;
	
	// Polak-Ribiere
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++){
			//b1 = b1 + d->dEdh[i][j]*(d->dEdh[i][j] - d->pdEdh[i][j]);
			b1 = b1 - d->dEdh[i][j]*(-d->dEdh[i][j] + d->pdEdh[i][j]);
			b2 = b2 + d->pdEdh[i][j]*d->pdEdh[i][j];
		}
	}
	
	if (fabs(b2) > 1e-12) {
		beta = b1 / b2;
		//beta = max(0.0,beta);  //2023.2.16
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++){
				d->cdlv[i][j] = -d->dEdh[i][j] + beta * d->cdlv[i][j];
			}
		}
	} else {
		beta = 0.0;  //2023.2.16
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++){			
				d->cdlv[i][j] = -d->dEdh[i][j];
			}	
		}
	}

	//printf("b1, b2, beta: %15.7f%15.7f%15.7f\n", b1, b2, beta);
}

void ann_opt::update_lv_cg(ann_data *d) {   //2023.1.31
	if (d->lv_update == "all") {
		for (int i=0; i<3; i++) {
			for(int j=0; j<3; j++){
				d->lv[i][j] += d->step_length_lv * d->cdlv[i][j];
				//d->lv[i][j] -= d->step_length_lv * d->cdlv[i][j];
			}
		}
	} else if (d->lv_update == "xx") {
		d->lv[0][0] += d->step_length_lv * d->cdlv[0][0];
	} else if (d->lv_update == "yy") {
		d->lv[1][1] += d->step_length_lv * d->cdlv[1][1];
	} else if (d->lv_update == "zz") {
		d->lv[2][2] += d->step_length_lv * d->cdlv[2][2];
	} else if (d->lv_update == "xyz") {
		d->lv[0][0] += d->step_length_lv * d->cdlv[0][0];
		d->lv[1][1] += d->step_length_lv * d->cdlv[1][1];
		d->lv[2][2] += d->step_length_lv * d->cdlv[2][2];
	} else if (d->lv_update == "wo_rot") {
		d->lv[0][0] += d->step_length_lv * d->cdlv[0][0];
		d->lv[1][1] += d->step_length_lv * d->cdlv[1][1];
		d->lv[2][2] += d->step_length_lv * d->cdlv[2][2];
		d->lv[1][0] += d->step_length_lv * d->cdlv[1][0];
		d->lv[2][0] += d->step_length_lv * d->cdlv[2][0];
		d->lv[2][1] += d->step_length_lv * d->cdlv[2][1];
	}

	/*
	printf("lv : ");
	for (int i=0; i<3; i++) {
		for(int j=0; j<3; j++){
			printf("%15.7f", d->lv[i][j]);
		}
	}
	printf("\n");
	*/
}