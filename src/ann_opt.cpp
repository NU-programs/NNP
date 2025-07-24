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

void ann_opt::opt_init(ann_data *d) {
	d->cdtx = new double[d->natom];
	d->cdty = new double[d->natom];
	d->cdtz = new double[d->natom];

	for (int i=0; i<d->natom; i++) {
		d->cdtx[i] = 0.;
		d->cdty[i] = 0.;
		d->cdtz[i] = 0.;
	}

	for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            d->cdlv[i][j] = 0.;
        }
    }
}

void ann_opt::copy_current_position(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->px[i] = d->x[i];
		d->py[i] = d->y[i];
		d->pz[i] = d->z[i];
	}
}

void ann_opt::copy_previous_position(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->x[i] = d->px[i];
		d->y[i] = d->py[i];
		d->z[i] = d->pz[i];
	}
}

void ann_opt::copy_current_force(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->pfcx[i] = d->fcx[i];
		d->pfcy[i] = d->fcy[i];
		d->pfcz[i] = d->fcz[i];
	}
}

void ann_opt::copy_previous_force(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->fcx[i] = d->pfcx[i];
		d->fcy[i] = d->pfcy[i];
		d->fcz[i] = d->pfcz[i];
	}
}

void ann_opt::copy_current_lv(ann_data *d) {    // 2023.1.30  Calc. for beta for conjugate direction
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            d->plv[i][j] = d->lv[i][j];
        }
    }
}

void ann_opt::copy_previous_lv(ann_data *d) {    // 2023.1.30  Calc. for beta for conjugate direction
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            d->lv[i][j] = d->plv[i][j];
        }
    }
}

void ann_opt::copy_current_dEdh(ann_data *d) {    // 2023.1.30  Calc. for beta for conjugate direction
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            d->pdEdh[i][j] = d->dEdh[i][j];
        }
    }
}

void ann_opt::copy_previous_dEdh(ann_data *d) {    // 2023.1.30  Calc. for beta for conjugate direction
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++){
            d->dEdh[i][j] = d->pdEdh[i][j];
        }
    }
}

int ann_opt::check_force_conv(ann_data *d) {
	int i, flag;
	
	flag = 0;
	for (i=0; i<d->natom; i++) {
		if (d->fc_diff < fabs(d->fcx[i])) {
			flag = 1;
			break;
		}
		
		if (d->fc_diff < fabs(d->fcy[i])) {
			flag = 1;
			break;
		}
		
		if (d->fc_diff < fabs(d->fcz[i])) {
			flag = 1;
			break;
		}
	}
	
	return flag;
}

int ann_opt::check_stress_conv(ann_data *d) {
    int i, j, flag;
    
    flag = 0;
    if (d->lv_update == "xx") {
        if (d->stress_diff < fabs(d->vstress[0])*1602.18) {
            flag = 1;
        }
    } else if (d->lv_update == "yy") {
        if (d->stress_diff < fabs(d->vstress[1])*1602.18) {
            flag = 1;
        }
    } else if (d->lv_update == "zz") {
        if (d->stress_diff < fabs(d->vstress[2])*1602.18) {
            flag = 1;
        }
    } else if (d->lv_update == "xyz") {
        for (i=0; i<3; i++) {
            if (d->stress_diff < fabs(d->vstress[i])*1602.18) {
                flag = 1;
                break;
            }
        }
    } else {
        for (i=0; i<6; i++) {
            if (d->stress_diff < fabs(d->vstress[i])*1602.18) {
                flag = 1;
                break;
            }
        }
    }
    
    return flag;
}

void ann_opt::calc_dEdh(ann_data *d) {
	int i,j,k;
	double inv_h_T[3][3],b[3][3];
	double stress[3][3];
	double buf;
	int n = 3;

	stress[0][0] = d->vstress[0];    // Sxx 
	stress[1][1] = d->vstress[1];    // Syy
	stress[2][2] = d->vstress[2];    // Szz
	stress[1][2] = d->vstress[3];    // Syz
	stress[0][2] = d->vstress[4];    // Sxz
	stress[0][1] = d->vstress[5];	 // Sxy
	stress[1][0] = d->vstress[5];	 // Syx   Syx = Sxy 
	stress[2][0] = d->vstress[4];	 // Szx   Szx = Sxz
	stress[2][1] = d->vstress[3];	 // Szy   Szy = Syz

	//printf("lattice vector lv[][]    transposed matrix originaly\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            //printf(" %10.5f",d->lv[i][j]);
            b[i][j] = d->lv[i][j];
        }
        //printf("\n");   
    }

    //making invert matrix
    //this->inv_matrix(d->lv,inv_h_T,n)  ... can not excute,because error   2023.1.26

    //1.making unit matrix
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            inv_h_T[i][j]=(i==j)?1.0:0.0;
        }
    }
    //2.sweeping method
    for(i=0;i<n;i++){
        buf=1/b[i][i];
        for(j=0;j<n;j++){
            b[i][j]*=buf;
            inv_h_T[i][j]*=buf;
        }
        for(j=0;j<n;j++){
            if(i!=j){
                buf=b[j][i];
                for(k=0;k<n;k++){
                    b[j][k]-=b[i][k]*buf;
                    inv_h_T[j][k]-=inv_h_T[i][k]*buf;
                }
            }
        }
    }

    //printf("inv_h_T =\n");
    //for (i=0; i<3; i++) {
    //    printf("%15.7f%15.7f%15.7f\n", inv_h_T[i][0], inv_h_T[i][1], inv_h_T[i][2]);
    //}

    //Calc. of dEdh[][] = -vol. * (inner product of stress[][] & inv_h_T[][])    2023.1.30
    for(i=0;i<3;i++){
    	for(j=0;j<3;j++){
            d->dEdh[i][j] = 0.0;
            for(k=0;k<3;k++){
                d->dEdh[i][j] -= d->volume * stress[i][k] * inv_h_T[k][j];
                //d->dEdh[i][j] -= d->volume * inv_h_T[i][k] * stress[k][j];
            }
        }
    }

    //printf("dEdh1 = \n");
    //for (i=0; i<3; i++) {
    //    printf("%15.7f%15.7f%15.7f\n", d->dEdh[i][0], d->dEdh[i][1], d->dEdh[i][2]);
    //}

    //Calc. of dEdh[][]= results of above calc. + (direct product of force[] and fractional coodinate[][])   2023.1.30
    for(i=0; i<d->natom; i++) {
        d->dEdh[0][0] -= d->fcx[i] * d->x[i];
        d->dEdh[0][1] -= d->fcx[i] * d->y[i];        
        d->dEdh[0][2] -= d->fcx[i] * d->z[i];        
        d->dEdh[1][0] -= d->fcy[i] * d->x[i];        
        d->dEdh[1][1] -= d->fcy[i] * d->y[i];
        d->dEdh[1][2] -= d->fcy[i] * d->z[i];
        d->dEdh[2][0] -= d->fcz[i] * d->x[i];
        d->dEdh[2][1] -= d->fcz[i] * d->y[i];
        d->dEdh[2][2] -= d->fcz[i] * d->z[i];
    }

    //printf("dEdh2 = \n");
    //for (i=0; i<3; i++) {
    //    printf("%15.7f%15.7f%15.7f\n", d->dEdh[i][0], d->dEdh[i][1], d->dEdh[i][2]);
    //}
}