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
#include "ann_ti.hpp"
using namespace std;

vector<string> ann_ti::split (string line) {
	istringstream iss(line);
	string s;
	vector<string> ss;
	while (iss >> s) ss.push_back(s);
	return ss;
}

void ann_ti::alloc(ann_data *d) {
	int i, j;
	
	d->x0 = new double[d->natom];
	d->y0 = new double[d->natom];
	d->z0 = new double[d->natom];
	for (i=0; i<d->natom; i++) {
		d->x0[i] = 0.;
		d->y0[i] = 0.;
		d->z0[i] = 0.;
	}
	
	d->disx = new double[d->natom];
	d->disy = new double[d->natom];
	d->disz = new double[d->natom];
	
	d->fcx_h = new double[d->natom];
	d->fcy_h = new double[d->natom];
	d->fcz_h = new double[d->natom];
	
	//d->fcx_ha = new double[d->natom];
	//d->fcy_ha = new double[d->natom];
	//d->fcz_ha = new double[d->natom];
	
	for (i=0; i<d->natom; i++) {
		d->disx[i] = 0.;
		d->disy[i] = 0.;
		d->disz[i] = 0.;
		
		d->fcx_h[i] = 0.;
		d->fcy_h[i] = 0.;
		d->fcz_h[i] = 0.;
		
		//d->fcx_ha[i] = 0.;
		//d->fcy_ha[i] = 0.;
		//d->fcz_ha[i] = 0.;
	}
	
	d->fc_const = new double*[d->natom3];
	for (i=0; i<d->natom3; i++) {
		d->fc_const[i] = new double[d->natom3];
		for (j=0; j<d->natom3; j++) {
			d->fc_const[i][j] = 0.;
		}
	}
}

void ann_ti::input_ti_pos(ann_data *d) {
	int i;
	string line;
	vector<string> s;
	ifstream fin;

	fin.open(d->ti_infile.c_str(), ios::in);
	if (!fin) {
		cout << d->ti_infile << " cannot be opened." << endl;
		exit(1);
	}

	for (i=0; i<7; i++) getline(fin, line);

	getline(fin, line);
	s = split(line);
	if (s[0] != "Direct") {
		cout << "not direct." << endl;
		exit(1);
	}

	for (i=0; i<d->natom; i++) {
		getline(fin, line); s = split(line);
		d->x0[i] = atof(s[0].c_str());
		d->y0[i] = atof(s[1].c_str());
		d->z0[i] = atof(s[2].c_str());
		//printf("%15.7f%15.7f%15.7f\n", d->x0[i], d->y0[i], d->z0[i]);
	}
	
	fin.close();
}

void ann_ti::init(ann_data *d) {
	d->ti_energy = 0.;
	d->ha_energy = 0.;
	d->ti_diff   = 0.;
}

void ann_ti::input_fc(ann_data *d) {
	int i, j;
	int m1, m2, m3, n1, n2, n3;

	string line;
	vector<string> s;
	vector<int> s_na;
	ifstream fin;
	fin.open(d->fc_infile.c_str(), ios::in);
	
	if (!fin) {
		cout << d->fc_infile << " cannot be opened." << endl;
		exit(1);
	}
	
	getline(fin, line);
	
	for (i=0; i<d->natom; i++) {
		m1 = 3*i;
		m2 = 3*i + 1;
		m3 = 3*i + 2;
		
		for (j=0; j<d->natom; j++) {
			getline(fin, line);
			
			n1 = 3*j;
			n2 = 3*j + 1;
			n3 = 3*j + 2;
			
			getline(fin, line); s = split(line);
			d->fc_const[m1][n1] = atof(s[0].c_str());
			d->fc_const[m1][n2] = atof(s[1].c_str());
			d->fc_const[m1][n3] = atof(s[2].c_str());
			
			getline(fin, line); s = split(line);
			d->fc_const[m2][n1] = atof(s[0].c_str());
			d->fc_const[m2][n2] = atof(s[1].c_str());
			d->fc_const[m2][n3] = atof(s[2].c_str());
			
			getline(fin, line); s = split(line);
			d->fc_const[m3][n1] = atof(s[0].c_str());
			d->fc_const[m3][n2] = atof(s[1].c_str());
			d->fc_const[m3][n3] = atof(s[2].c_str());
		}
	}
		
	fin.close();
	
	//m1 = 3*(d->natom-1);
	//printf("%15.7f%15.7f%15.7f\n", d->fc_const[m1][m1], d->fc_const[m1][m1+1], d->fc_const[m1][m1+2]);
}

void ann_ti::disp(ann_data *d) {
	double dx, dy, dz;
	
	for (int i=0; i<d->natom; i++) {
		dx = d->x[i] - d->x0[i];
		dy = d->y[i] - d->y0[i];
		dz = d->z[i] - d->z0[i];
		
		if (dx < -0.5) { dx = (d->x[i] + 1.) - d->x0[i]; }
		else if (dx > 0.5) { dx = (d->x[i] - 1.) - d->x0[i]; }
		if (dy < -0.5) { dy = (d->y[i] + 1.) - d->y0[i]; }
		else if (dy > 0.5) { dy = (d->y[i] - 1.) - d->y0[i]; }
		if (dz < -0.5) { dz = (d->z[i] + 1.) - d->z0[i]; }
		else if (dz > 0.5) { dz = (d->z[i] - 1.) - d->z0[i]; }
		
		d->disx[i] = dx*d->lv[0][0] + dy*d->lv[1][0] + dz*d->lv[2][0];
		d->disy[i] = dx*d->lv[0][1] + dy*d->lv[1][1] + dz*d->lv[2][1];
		d->disz[i] = dx*d->lv[0][2] + dy*d->lv[1][2] + dz*d->lv[2][2];
		//printf("%10d%20.15f%20.15f%20.15f\n", i, d->disx[i], d->disy[i], d->disz[i]);
	}
}
	
void ann_ti::calc_fc(ann_data *d) {
	int i, j;
	int m1, m2, m3, n1, n2, n3;
	
	//for (i=0; i<d->natom; i++) {
	//	d->fcx_ha[i] = d->fcx[i];
	//	d->fcy_ha[i] = d->fcy[i];
	//	d->fcz_ha[i] = d->fcz[i];
	//}
	
	// calc force	
	for (i=0; i<d->natom; i++) {
		m1 = 3*i;
		m2 = 3*i + 1;
		m3 = 3*i + 2;
		
		d->fcx_h[i] = 0.;
		d->fcy_h[i] = 0.;
		d->fcz_h[i] = 0.;
		
		for (j=0; j<d->natom; j++) {
			n1 = 3*j;
			n2 = 3*j + 1;
			n3 = 3*j + 2;
			
			d->fcx_h[i] += (d->fc_const[m1][n1]*d->disx[j] + d->fc_const[m1][n2]*d->disy[j] + d->fc_const[m1][n3]*d->disz[j]);
			d->fcy_h[i] += (d->fc_const[m2][n1]*d->disx[j] + d->fc_const[m2][n2]*d->disy[j] + d->fc_const[m2][n3]*d->disz[j]);
			d->fcz_h[i] += (d->fc_const[m3][n1]*d->disx[j] + d->fc_const[m3][n2]*d->disy[j] + d->fc_const[m3][n3]*d->disz[j]);
		}
		
		d->fcx[i] = -d->ti_lambda_h * d->fcx_h[i] + d->ti_lambda * d->fcx[i];
		d->fcy[i] = -d->ti_lambda_h * d->fcy_h[i] + d->ti_lambda * d->fcy[i];
		d->fcz[i] = -d->ti_lambda_h * d->fcz_h[i] + d->ti_lambda * d->fcz[i];
	}
	
	// calc_energy
	d->ha_energy = 0;
	for (i=0; i<d->natom; i++) d->ha_energy += (d->fcx_h[i]*d->disx[i] + d->fcy_h[i]*d->disy[i] + d->fcz_h[i]*d->disz[i]);
	d->ha_energy = 0.5*d->ha_energy + d->ha_energy0;
	
	d->ti_energy = d->ti_lambda_h * d->ha_energy + d->ti_lambda * d->energy;
	d->ti_diff = d->energy - d->ha_energy;
}
