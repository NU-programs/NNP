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
#include "ann_npt_nh.hpp"
using namespace std;

void ann_npt_nh::calc_npt_const(ann_data *d) {
	d->mass_sr = new double[d->natom];
	for (int i=0; i<d->natom; i++) d->mass_sr[i] = sqrt(d->mass[i]);

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			d->ba_fric[i][j] = 0.;
		}
	}

	d->ba_mass = 3.*d->natom * (kB * d->md_temp) * d->ba_mass*d->ba_mass;
	d->md_pressure_cx = d->md_pressure_x * 6.24151 * 1e-3;
	d->md_pressure_cy = d->md_pressure_y * 6.24151 * 1e-3;
	d->md_pressure_cz = d->md_pressure_z * 6.24151 * 1e-3;
	//printf("th_mass = %20.10f\n", d->th_mass);

	d->md_time_4 = 0.25 * d->md_time;
	d->md_time_8 = 0.125 * d->md_time;
	//d->md_time_4c = d->md_time_4 / d->ba_mass;
	//d->md_time_8c = d->md_time_8 / d->th_mass;
	//printf("md_time_4c: %15.7f %15.7f %15.7f\n", d->md_time_4, d->ba_mass, d->md_time_4/d->ba_mass);
	d->t_ke = 0.5 * (3.*d->natom - 3.) * (kB * d->md_temp);
	d->th_fric_c = 2.*d->t_ke + 9.*kB*d->md_temp;
}

void ann_npt_nh::calc_ba_fric_tr(ann_data *d) {
    d->ba_fric_tr = d->ba_fric[0][0]*d->ba_fric[0][0] + d->ba_fric[1][0]*d->ba_fric[1][0] + d->ba_fric[2][0]*d->ba_fric[2][0] +
                    d->ba_fric[0][1]*d->ba_fric[0][1] + d->ba_fric[1][1]*d->ba_fric[1][1] + d->ba_fric[2][1]*d->ba_fric[2][1] +
                    d->ba_fric[0][2]*d->ba_fric[0][2] + d->ba_fric[1][2]*d->ba_fric[1][2] + d->ba_fric[2][2]*d->ba_fric[2][2];
}

void ann_npt_nh::update_th_friction(ann_data *d) {
	d->c_ke = 0.;
	for (int i=0; i<d->natom; i++) {
		d->c_ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	//d->c_ke *= (0.5 * d->c_unit2);
	d->c_ke *= d->c_unit2;
	
	//d->fric += d->md_time_8c * (2.*d->c_ke + d->ba_mass*d->ba_fric_tr - 2.*d->t_ke - 9.*kB*d->md_temp);
	d->fric += d->md_time_8 * (d->c_ke + d->ba_mass*d->ba_fric_tr - d->th_fric_c) / d->th_mass;
}

void ann_npt_nh::update_th_velocity(ann_data *d) {
    double c = exp(-d->fric*d->md_time_4);

    for (int i=0; i<d->natom; i++) {
        d->vx[i] *= c;
        d->vy[i] *= c;
        d->vz[i] *= c;
    }
}

void ann_npt_nh::calc_total_stress(ann_data *d) {
    int i, j;
    double dx, dy, dz;
    double kvxx, kvyy, kvzz, kvxy, kvxz, kvyz;

    d->stress_t[0][0] = d->vstress[0] * d->volume;    // Sxx 
    d->stress_t[1][1] = d->vstress[1] * d->volume;    // Syy
    d->stress_t[2][2] = d->vstress[2] * d->volume;    // Szz
    d->stress_t[1][2] = d->vstress[3] * d->volume;    // Syz
    d->stress_t[0][2] = d->vstress[4] * d->volume;    // Sxz
    d->stress_t[0][1] = d->vstress[5] * d->volume;    // Sxy
    //d->stress_t[1][0] = d->vstress[5] * d->volume;    // Syx   Syx = Sxy 
    //d->stress_t[2][0] = d->vstress[4] * d->volume;    // Szx   Szx = Sxz
    //d->stress_t[2][1] = d->vstress[3] * d->volume;    // Szy   Szy = Syz

    /*
    d->vxa = 0.; d->vya = 0.; d->vza = 0.;
    for (i=0; i<d->natom; i++) {
        d->vxa += d->vx[i];
        d->vya += d->vy[i];
        d->vza += d->vz[i];
    }
    d->vxa /= d->natom;
    d->vya /= d->natom;
    d->vza /= d->natom;
    */

    kvxx = 0.; kvyy = 0.; kvzz = 0.;
    kvxy = 0.; kvxz = 0.; kvyz = 0.;
    for (i=0; i<d->natom; i++) {
    	//dx = d->mass_sr[i] * (d->vx[i] - d->vxa);
    	//dy = d->mass_sr[i] * (d->vy[i] - d->vya);
    	//dz = d->mass_sr[i] * (d->vz[i] - d->vza);
    	//kvxx -= dx * dx; kvyy -= dy * dy; kvzz -= dz * dz;
    	//kvxy -= dx * dy; kvxz -= dx * dz; kvyz -= dy * dz;
    	
    	dx = d->mass_sr[i] * d->vx[i];
    	dy = d->mass_sr[i] * d->vy[i];
    	dz = d->mass_sr[i] * d->vz[i];
    	kvxx += dx * dx; kvyy += dy * dy; kvzz += dz * dz;
    	kvxy += dx * dy; kvxz += dx * dz; kvyz += dy * dz;
    }
    kvxx *= d->c_unit2;  kvyy *= d->c_unit2;  kvzz *= d->c_unit2;
    kvxy *= d->c_unit2;  kvxz *= d->c_unit2;  kvyz *= d->c_unit2;

    d->stress_t[0][0] += kvxx;
    d->stress_t[1][1] += kvyy;
    d->stress_t[2][2] += kvzz;
    d->stress_t[0][1] += kvxy;
    d->stress_t[1][0] = d->stress_t[0][1];
    d->stress_t[0][2] += kvxz;
    d->stress_t[2][0] = d->stress_t[0][2];
    d->stress_t[1][2] += kvyz;
    d->stress_t[2][1] = d->stress_t[1][2];
    //printf("stress_t: %15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n",
    //	d->stress_t[0][0]/d->volume*evang3_to_kbar, d->stress_t[1][1]/d->volume*evang3_to_kbar, d->stress_t[2][2]/d->volume*evang3_to_kbar,
    //	d->stress_t[0][1]/d->volume*evang3_to_kbar, d->stress_t[0][2]/d->volume*evang3_to_kbar, d->stress_t[1][2]/d->volume*evang3_to_kbar);
}

void ann_npt_nh::update_ba_friction(ann_data *d) {
    int i, j;
    double c1 = exp(-d->fric * d->md_time_8);
    double c2x = d->md_pressure_cx * d->volume;
    double c2y = d->md_pressure_cy * d->volume;
    double c2z = d->md_pressure_cz * d->volume;

    //printf("fric, md_time_8, c1 = %15.7f%15.7f%15.7f\n", d->fric, d->md_time_8, c1);
    //printf("c2 = %15.7f\n", c2);

    //printf("ba_fric1: %15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n",
    // 	d->ba_fric[0][0], d->ba_fric[1][1], d->ba_fric[2][2], d->ba_fric[0][1], d->ba_fric[0][2], d->ba_fric[1][2]);

    d->ba_fric[0][0] *= c1;
    d->ba_fric[1][1] *= c1;
    d->ba_fric[2][2] *= c1;
    d->ba_fric[0][1] *= c1;
    //d->ba_fric[1][0] = d->ba_fric[0][1];
    d->ba_fric[0][2] *= c1;
    //d->ba_fric[2][0] = d->ba_fric[0][2];
    d->ba_fric[1][2] *= c1;
    //d->ba_fric[2][1] = d->ba_fric[1][2];

    //printf("ba_fric2: %15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n",
    //	d->ba_fric[0][0], d->ba_fric[1][1], d->ba_fric[2][2], d->ba_fric[0][1], d->ba_fric[0][2], d->ba_fric[1][2]);

    //printf("md_time_4c, stress_t[0][0]: %15.7f%15.7f\n", d->md_time_4c, d->stress_t[0][0]);

    d->ba_fric[0][0] += d->md_time_4 * (d->stress_t[0][0] - c2x) / d->ba_mass;
    d->ba_fric[1][1] += d->md_time_4 * (d->stress_t[1][1] - c2y) / d->ba_mass;
    d->ba_fric[2][2] += d->md_time_4 * (d->stress_t[2][2] - c2z) / d->ba_mass;
    d->ba_fric[0][1] += d->md_time_4 * d->stress_t[0][1] / d->ba_mass;
    //d->ba_fric[1][0] = d->ba_fric[0][1];
    d->ba_fric[0][2] += d->md_time_4 * d->stress_t[0][2] / d->ba_mass;
    //d->ba_fric[2][0] = d->ba_fric[0][2];
    d->ba_fric[1][2] += d->md_time_4 * d->stress_t[1][2] / d->ba_mass;
    //d->ba_fric[2][1] = d->ba_fric[2][1];

    //printf("ba_fric3: %15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n",
    //	d->ba_fric[0][0], d->ba_fric[1][1], d->ba_fric[2][2], d->ba_fric[0][1], d->ba_fric[0][2], d->ba_fric[1][2]);

    d->ba_fric[0][0] *= c1;
    d->ba_fric[1][1] *= c1;
    d->ba_fric[2][2] *= c1;
    d->ba_fric[0][1] *= c1;
    d->ba_fric[1][0] = d->ba_fric[0][1];
    d->ba_fric[0][2] *= c1;
    d->ba_fric[2][0] = d->ba_fric[0][2];
    d->ba_fric[1][2] *= c1;
    d->ba_fric[2][1] = d->ba_fric[1][2];

    //printf("ba_fric4: %15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n",
    //	d->ba_fric[0][0], d->ba_fric[1][1], d->ba_fric[2][2], d->ba_fric[0][1], d->ba_fric[0][2], d->ba_fric[1][2]);
}

void ann_npt_nh::update_ba_velocity(ann_data *d) {
	int i, j, k;
	double vxx, vyy, vzz;
	double aaa[3][3], bbb[3][3];

	//printf("aaa: ");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			aaa[i][j] = -d->ba_fric[i][j] * d->md_time_h;
			//printf("%15.5f", aaa[i][j]);
		}
	}
	//printf("\n");

	//printf("bbb: ");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			bbb[i][j] = 0.;
			for (k=0; k<3; k++) {
				bbb[i][j] += aaa[i][k]*aaa[k][j];
			}
			//printf("%15.5f", bbb[i][j]);
		}
	}
	//printf("\n");

	//printf("ba_exp1: ");
	for (i=0; i<3; i++) {
		for (j=0; j<=i; j++) {
			d->ba_exp1[i][j] = aaa[i][j] + 0.5*bbb[i][j];
			d->ba_exp1[j][i] = d->ba_exp1[i][j];
			if (i == j) {
				d->ba_exp1[i][i] += 1.;
			}
			//printf("%15.5f", d->ba_exp1[i][j]);
		}
	}
	//printf("\n");

    for (i=0; i<d->natom; i++) {
    	vxx = d->ba_exp1[0][0]*d->vx[i] + d->ba_exp1[0][1]*d->vy[i] + d->ba_exp1[0][2]*d->vz[i];
    	vyy = d->ba_exp1[1][0]*d->vx[i] + d->ba_exp1[1][1]*d->vy[i] + d->ba_exp1[1][2]*d->vz[i];
    	vzz = d->ba_exp1[2][0]*d->vx[i] + d->ba_exp1[2][1]*d->vy[i] + d->ba_exp1[2][2]*d->vz[i];
        d->vx[i] = vxx;
        d->vy[i] = vyy;
        d->vz[i] = vzz;
    }
}

void ann_npt_nh::calc_exp2(ann_data *d) {
	int i, j, k;
	double aaa[3][3], bbb[3][3];

	//printf("aaa: ");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			aaa[i][j] = d->ba_fric[i][j] * d->md_time;
			//printf("%15.5f", aaa[i][j]);
		}
	}
	//printf("\n");

	//printf("bbb: ");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			bbb[i][j] = 0.;
			for (k=0; k<3; k++) {
				bbb[i][j] += aaa[i][k]*aaa[k][j];
			}
			//printf("%15.5f", bbb[i][j]);
		}
	}
	//printf("\n");

	//printf("ba_exp2: ");
	for (i=0; i<3; i++) {
		for (j=0; j<=i; j++) {
			d->ba_exp2[i][j] = aaa[i][j] + 0.5*bbb[i][j];
			d->ba_exp2[j][i] = d->ba_exp2[i][j];
			if (i == j) {
				d->ba_exp2[i][i] += 1.;
			}
			//printf("%15.5f", d->ba_exp2[i][j]);
		}
	}
	//printf("\n");
}

void ann_npt_nh::update_lv(ann_data *d) {
	int i, j, k;
	double aaa[3][3];

	//printf("aaa: ");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			aaa[i][j] = 0.;
			for (k=0; k<3; k++) {
				aaa[i][j] += d->lv[i][k]*d->ba_exp2[k][j];
			}
			//printf("%15.5f", aaa[i][j]);
		}
	}
	//printf("\n");

	/*
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			d->lv[i][j] = aaa[i][j];
		}
	}
	*/

	if (d->npt_shear == "no") {
		d->lv[0][0] = aaa[0][0];
		d->lv[1][1] = aaa[1][1];
		d->lv[2][2] = aaa[2][2];
	} else {
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				d->lv[i][j] = aaa[i][j];
			}
		}
	}
}

void ann_npt_nh::update_position(ann_data *d) {
	int i;
	double cxx, cyy, czz;
	double pxx, pyy, pzz;

	d->mcx = 0.;
	d->mcy = 0.;
	d->mcz = 0.;
	for (i=0; i<d->natom; i++) {
		d->mcx += d->mass[i]*d->cx[i];
		d->mcy += d->mass[i]*d->cy[i];
		d->mcz += d->mass[i]*d->cz[i];
	}
	d->mcx /= d->total_mass;
	d->mcy /= d->total_mass;
	d->mcz /= d->total_mass;

	for (i=0; i<d->natom; i++) {
		cxx = (d->cx[i] - d->mcx);
		cyy = (d->cy[i] - d->mcy);
		czz = (d->cz[i] - d->mcz);

		pxx = d->ba_exp2[0][0]*cxx + d->ba_exp2[0][1]*cyy + d->ba_exp2[0][2]*czz;
		pyy = d->ba_exp2[1][0]*cxx + d->ba_exp2[1][1]*cyy + d->ba_exp2[1][2]*czz;
		pzz = d->ba_exp2[2][0]*cxx + d->ba_exp2[2][1]*cyy + d->ba_exp2[2][2]*czz;

		d->cx[i] = pxx + d->md_time * d->vx[i] + d->mcx;
		d->cy[i] = pyy + d->md_time * d->vy[i] + d->mcy;
		d->cz[i] = pzz + d->md_time * d->vz[i] + d->mcz; 
	}
}

void ann_npt_nh::update_velocity(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		//c = d->md_const2 / d->mass[i];
		d->vx[i] += d->mass_const2[i] * d->fcx[i];
		d->vy[i] += d->mass_const2[i] * d->fcy[i];
		d->vz[i] += d->mass_const2[i] * d->fcz[i];
	}
}

void ann_npt_nh::calc_current_temperature(ann_data *d) {
	d->c_temp = 0;
	for (int i=0; i<d->natom; i++) {
		d->c_temp += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->c_temp = d->c_temp * 0.5 * d->c_unit2 / (1.5  * (double)d->natom * kB);
	
	if (d->md_dtemp != 0.) {
		d->md_temp += d->md_dtemp;
		d->t_ke = 0.5 * (3.*d->natom - 3.) * (kB * d->md_temp);
	}
}
