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
#include "ann_sf.hpp"
using namespace std;

void ann_sf::init_array(ann_data *d) {
	int i, j, k;

	for (i=0; i<d->n_atproc; i++) {
		for (j=0; j<d->n_symmetry; j++) {
			d->sfv[i][j]  = 0.;
		}
			
		for (j=0; j<d->cf_natom2[i]; j++) {
			for (k=0; k<d->n_symmetry; k++) {
				d->sfdx[i][j][k] = 0.;
				d->sfdy[i][j][k] = 0.;
				d->sfdz[i][j][k] = 0.;
			}
		}
		
		//for (j=0; j<d->n_symmetry; j++) {
		//	d->sfdx[i][d->max_natom][j] = 0.;
		//	d->sfdy[i][d->max_natom][j] = 0.;
		//	d->sfdz[i][d->max_natom][j] = 0.;
		//}

		//for (j=0; j<d->cf_natom2[i]; j++) {
		//	for (k=0; k<d->n_symmetry; k++) {
		//		d->sfdx_i[i][j][k] = 0.;
		//		d->sfdy_i[i][j][k] = 0.;
		//		d->sfdz_i[i][j][k] = 0.;
		//	}
		//}
	}

	for (i=0; i<d->n_atproc; i++) {
		d->cf_natom2[i] = 0;
		
		for (j=0; j<d->cf_natom2[i]; j++) {
			d->cf_pair1[i][j] = 0;
		}
	}
}

void ann_sf::calc_distance(ann_data *d){
	int w, c0, n, mc[3], ii, jj;
	int i, j, k, c, q, l, m, o, c1, p, s;
	double xi, yi, zi;
	double tdx, tdy, tdz;
	double delx, dely, delz;
	double delxx, delyy, delzz;
	double dis;
	
	for (i=0; i<d->lcxyz; i++) d->cell_natom[i] = 0;
	
	for (w=0; w<d->natom; w++) {
		mc[0] = int(floor(d->x[w] * d->lc[0]));
		mc[1] = int(floor(d->y[w] * d->lc[1]));
		mc[2] = int(floor(d->z[w] * d->lc[2]));
		//if (d->rank == 0) {
		//	printf("%15.7f%15.7f%15.7f%10d%10d%10d\n", d->x[w], d->y[w], d->z[w], mc[0], mc[1], mc[2]);
		//}
		c0 = mc[0]*d->lcyz + mc[1]*d->lc[2] + mc[2];
		
		d->cell_atom_list[c0][d->cell_natom[c0]] = w;
		d->cell_num[w][0] = mc[0];
		d->cell_num[w][1] = mc[1];
		d->cell_num[w][2] = mc[2];
		
		d->cell_natom[c0]++;

		if (d->cell_natom[c0] == d->max_lc_natom) {
			if (d->rank == 0) printf("cell_natom[c0] is larger than max_lc_natom!!!\n");
			MPI_Abort(MPI_COMM_WORLD, 22);
		}
		//cout << c0 << " " << d->cell_natom[c0] << endl;
	}

	//if (d->rank == 0) {
	//	for (i=0; i<d->lcxyz; i++) printf("cell: %10d%10d\n", i, d->cell_natom[i]);
	//}
	
	for (i=0; i<d->n_atproc; i++) {
		n = i + d->n_start;
		
		q = d->cell_num[n][0];
		j = d->cell_num[n][1];
		k = d->cell_num[n][2];
		//c = q*d->lcyz + j*d->lc[2] + k;
			
		for (l = q-d->neighbor[0]; l <= q+d->neighbor[0]; l++) {
			if (d->flagx == 0){
				if (l < 0) { xi = -1.; } 
				else if (l >= d->lc[0]) { xi = 1.; } 
				else { xi = 0.; }
			} else {
				xi = float(l);
			}
			
			for (m = j-d->neighbor[1]; m <= j+d->neighbor[1]; m++) {
				if (d->flagy == 0){
					if (m < 0) { yi = -1.; } 
					else if (m >= d->lc[1]) { yi = 1.; } 
					else { yi = 0.; }
				} else {
					yi = float(m);
				}
				
				for (o = k-d->neighbor[2]; o <= k+d->neighbor[2]; o++) {
					if (d->flagz == 0){
						if (o < 0) { zi = -1.; } 
						else if (o >= d->lc[2]) { zi = 1.; } 
						else { zi = 0.; }
					} else {
						zi = float(o);
					}
					
					c1 = ((l + d->lc[0]) % d->lc[0]) * d->lcyz + ((m + d->lc[1]) % d->lc[1]) * d->lc[2] + (o + d->lc[2]) % d->lc[2];
					//cout << c << " " << c1 << endl;
					
					for (p=0; p<d->cell_natom[c1]; p++){
						s = d->cell_atom_list[c1][p];
						//cout << c << " " << c1 << " " << n << " " << s << endl;
						
						tdx = d->x[s] + xi;
						tdy = d->y[s] + yi;
						tdz = d->z[s] + zi;
						//cout << tdx << endl;							
						
						delx = tdx - d->x[n];
						dely = tdy - d->y[n];
						delz = tdz - d->z[n];
						//cout << delx << endl;
								
						delxx = d->lv[0][0]*delx + d->lv[1][0]*dely + d->lv[2][0]*delz;
						//cout << delxx <<endl;					
							
						if (fabs(delxx) <= d->cutoff) {
							delyy = d->lv[0][1]*delx + d->lv[1][1]*dely + d->lv[2][1]*delz;
								
							if (fabs(delyy) <= d->cutoff) {
								delzz = d->lv[0][2]*delx + d->lv[1][2]*dely + d->lv[2][2]*delz;
										
								if (fabs(delzz) <= d->cutoff) {
									dis = delxx*delxx + delyy*delyy + delzz*delzz;
								
									if (dis <= d->cutoff2) {
										if (1e-3 <= dis) {
											if (d->cf_natom2[i] == d->max_natom) {
												cout << "The number of cf_natom2 is larger than max_natom" << endl;
												MPI_Abort(MPI_COMM_WORLD, 26);
											}
	
											ii = d->cf_natom2[i];
											//cout << n << " " << s << " " << ii << " " << d->cf_natom2[i] << endl;
											d->cf_pair1[i][ii] = s;
													
											d->dx1[i][ii] = delxx;
											d->dy1[i][ii] = delyy;
											d->dz1[i][ii] = delzz;
											//d->dis12[i][ii] = dis;
											d->dis12[i][ii] = dis;
											d->dis1[i][ii] = sqrt(dis);
												
											//d->x1[i][ii] = tdx;
											//d->y1[i][ii] = tdy;
											//d->z1[i][ii] = tdz;
							
											d->cf_natom2[i]++;
										}
									}					
								}
							}
						}
					}
				}
			}
		}
	}
	//printf("cf_natom[0] = %10d\n", d->cf_natom2[0]);
	//for (i=0; i<d->cf_natom2[0]; i++) {
	//	printf("%10d%10d\n", i, d->cf_pair1[0][i]);
	//}
}

///////////////////////////////////////////////////////////////////////////////////
void ann_sf::calc_g2(ann_data *d) {
	double vv, vx, vy, vz;
	int i, j, ii, jj, n, nij0, nij;
	int e_i, e_j;
	double c1, c2;
	double chv[3], chvdev[3];
	
	for (i=0; i<d->n_atproc; i++) {
		ii = i + d->n_start;

		for (j=0; j<d->cf_natom2[i]; j++) {
			e_j = d->e[d->cf_pair1[i][j]];
			jj = d->cf_pair1[i][j];
			nij0 = d->ng2 * e_j;

			c1 = tanh(1. + d->tanh_const1 * d->dis1[i][j]);
			c2 = c1 * c1;

			d->fc1[i][j] = d->tanh_const0 * c1 * c2;
			//d->fc1d[i][j] = d->tanh_const2 * c2 * (1. - c2);
			d->fc1d[i][j] = d->tanh_const2 * c2 * (1. - c2) / d->dis1[i][j];
			d->gv2[i][j] = d->ch_const1 * d->dis1[i][j] - 1.;

			chv[0] = 1.;
			chv[1] = d->gv2[i][j];

			chvdev[0] = 0.;
			chvdev[1] = 1.;
			
			for (n=0; n<d->ng2; n++) {
				nij = n + nij0;

				if (n == 0) {
					chv[2] = 1.;
					chvdev[2] = 0.;
				} else if (n == 1) {
					chv[2] = d->gv2[i][j];
					chvdev[2] = 1.;
				} else {
					chv[2] = 2. * d->gv2[i][j] * chv[1] - chv[0];
					chvdev[2] = 2.*chv[1] + 2.*d->gv2[i][j]*chvdev[1] - chvdev[0];

					chv[0] = chv[1];
					chv[1] = chv[2];
					chvdev[0] = chvdev[1];
					chvdev[1] = chvdev[2];
				}

				d->sfv[i][nij] += (chv[2] * d->fc1[i][j]);

				vv = (chvdev[2] * d->ch_const1 * d->fc1[i][j] / d->dis1[i][j] + chv[2] * d->fc1d[i][j]);
				vx = vv * (-d->dx1[i][j]);
				vy = vv * (-d->dy1[i][j]);
				vz = vv * (-d->dz1[i][j]);

				//d->sfdx[i][d->max_natom][nij] += vx;
				//d->sfdy[i][d->max_natom][nij] += vy;
				//d->sfdz[i][d->max_natom][nij] += vz;

				//d->sfdx_i[i][j][nij] += vx;
				//d->sfdy_i[i][j][nij] += vy;
				//d->sfdz_i[i][j][nij] += vz;

				d->sfdx[i][j][nij] -= vx;
				d->sfdy[i][j][nij] -= vy;
				d->sfdz[i][j][nij] -= vz;
			}
		}
	}
}

void ann_sf::calc_g3_term(ann_data *d) {
	double d_jk, c1, c2;
	int i, j, k, ii, jj, kk;
	
	for (i=0; i<d->n_atproc; i++) {
		ii = i + d->n_start;

		for (j=1; j<d->cf_natom2[i]; j++) {
			jj = d->cf_pair1[i][j];
			
			//for (k=j+1; k<d->cf_natom2[i]; k++) {
			for (k=0; k<j; k++) {
				kk = d->cf_pair1[i][k];
				
				d_jk = d->dis1[i][j] * d->dis1[i][k];
				d->cosine[i][j][k] = (d->dx1[i][j]*d->dx1[i][k] + d->dy1[i][j]*d->dy1[i][k] + d->dz1[i][j]*d->dz1[i][k]) / d_jk;
				d->gv3[i][j][k] = 0.5*(d->cosine[i][j][k] + 1.);

				// derivative with respect to x_i
				c1 = d->cosine[i][j][k] / d->dis12[i][j];
				c2 = d->cosine[i][j][k] / d->dis12[i][k];
				//d->dx_cosine[i][j][k] = -(d->dx1[i][j] + d->dx1[i][k]) / d_jk + c1 * d->dx1[i][j] + c2 * d->dx1[i][k];
				//d->dy_cosine[i][j][k] = -(d->dy1[i][j] + d->dy1[i][k]) / d_jk + c1 * d->dy1[i][j] + c2 * d->dy1[i][k];
				//d->dz_cosine[i][j][k] = -(d->dz1[i][j] + d->dz1[i][k]) / d_jk + c1 * d->dz1[i][j] + c2 * d->dz1[i][k];

				// derivative with respect to x_j
				//c1 = d->cosine[i][j][k] / d->dis12[i][j];
				d->pdx_cosine[i][j][k] = d->dx1[i][k] / d_jk - c1 * d->dx1[i][j];
				d->pdy_cosine[i][j][k] = d->dy1[i][k] / d_jk - c1 * d->dy1[i][j];
				d->pdz_cosine[i][j][k] = d->dz1[i][k] / d_jk - c1 * d->dz1[i][j];

				// derivative with respect to x_k
				//c1 = d->cosine[i][j][k] / d->dis12[i][k];
				d->pdx_cosine[i][k][j] = d->dx1[i][j] / d_jk - c2 * d->dx1[i][k];
				d->pdy_cosine[i][k][j] = d->dy1[i][j] / d_jk - c2 * d->dy1[i][k];
				d->pdz_cosine[i][k][j] = d->dz1[i][j] / d_jk - c2 * d->dz1[i][k];
			}
		}
	}
}

void ann_sf::calc_g3(ann_data *d) {
	double chv[3], chvdev[3];
	double c1, c2, c3;
	double fcij_fcik, h_fcij_fcik, dfcij_fcik, fcij_dfcik;
	double vx, vy, vz;

	int i, j, k, ii, jj, kk, n, njk, njk0;
	int e_j, e_k;
	
	for (i=0; i<d->n_atproc; i++) {
		ii = i + d->n_start;
		//printf("%10d%10d\n", ii, d->cf_natom2[i]);

		for (j=1; j<d->cf_natom2[i]; j++) {
			e_j = d->e[d->cf_pair1[i][j]];
			jj = d->cf_pair1[i][j];

			for (k=0; k<j; k++) {
				e_k = d->e[d->cf_pair1[i][k]];
				kk = d->cf_pair1[i][k];
				
				njk0 = d->t_ng2 + d->ng3 * d->g3x[e_j][e_k];

				chv[0] = 1.;
				chv[1] = d->gv3[i][j][k];

				chvdev[0] = 0.;
				chvdev[1] = 1.;

				fcij_fcik = d->fc1[i][j] * d->fc1[i][k];
				h_fcij_fcik = 0.5 * fcij_fcik;
				dfcij_fcik = d->fc1d[i][j] * d->fc1[i][k];
				fcij_dfcik = d->fc1[i][j] * d->fc1d[i][k];

				for (n=0; n<d->ng3; n++) {
					njk = n + njk0;

					if (n == 0) {
						chv[2] = 1.;
						chvdev[2] = 0.;
					} else if (n == 1) {
						chv[2] = d->gv3[i][j][k];
						chvdev[2] = 1.;
					} else {
						chv[2] = 2. * d->gv3[i][j][k] * chv[1] - chv[0];
						chvdev[2] = 2.*chv[1] + 2.*d->gv3[i][j][k]*chvdev[1] - chvdev[0];
							
						chv[0] = chv[1];
						chv[1] = chv[2];
						chvdev[0] = chvdev[1];
						chvdev[1] = chvdev[2];
					}

					d->sfv[i][njk] += chv[2] * fcij_fcik;
					//d->sfv[i][njk] += (chv[2] * d->fc1[i][j] * d->fc1[i][k]);

					c1 = chvdev[2] * h_fcij_fcik;
					//c1 = chvdev[2] * 0.5 * fcij_fcik;
					//c1 = chvdev[2] * 0.5 * d->fc1[i][j] * d->fc1[i][k];
						
					// derivative with respect to atom_i
					c2 = chv[2] * dfcij_fcik;
					c3 = chv[2] * fcij_dfcik;
					//c2 = chv[2] * d->fc1d[i][j] * d->fc1[i][k];
					//c3 = chv[2] * d->fc1[i][j] * d->fc1d[i][k];

					// derivative with respect to atom_i
					//d->sfdx[i][d->max_natom][njk] += (c1 * d->dx_cosine[i][j][k] - c2 * d->dx1[i][j] - c3 * d->dx1[i][k]);
					//d->sfdy[i][d->max_natom][njk] += (c1 * d->dy_cosine[i][j][k] - c2 * d->dy1[i][j] - c3 * d->dy1[i][k]);
					//d->sfdz[i][d->max_natom][njk] += (c1 * d->dz_cosine[i][j][k] - c2 * d->dz1[i][j] - c3 * d->dz1[i][k]);

					// derivative with respect to atom_j
					//c2 = chv[2] * d->fc1d[i][j] * d->fc1[i][k];
					vx = (c1 * d->pdx_cosine[i][j][k] + c2 * d->dx1[i][j]);
					vy = (c1 * d->pdy_cosine[i][j][k] + c2 * d->dy1[i][j]);
					vz = (c1 * d->pdz_cosine[i][j][k] + c2 * d->dz1[i][j]);
					d->sfdx[i][j][njk] += vx;
					d->sfdy[i][j][njk] += vy;
					d->sfdz[i][j][njk] += vz;
					//d->sfdx_i[i][j][njk] -= vx;
					//d->sfdy_i[i][j][njk] -= vy;
					//d->sfdz_i[i][j][njk] -= vz;
					//d->sfdx[i][j][njk] += (c1 * d->pdx_cosine[i][j][k] + c2 * d->dx1[i][j]);
					//d->sfdy[i][j][njk] += (c1 * d->pdy_cosine[i][j][k] + c2 * d->dy1[i][j]);
					//d->sfdz[i][j][njk] += (c1 * d->pdz_cosine[i][j][k] + c2 * d->dz1[i][j]);

					// derivative with respect to atom_k
					//c2 = chv[2] * d->fc1[i][j] * d->fc1d[i][k];
					vx = (c1 * d->pdx_cosine[i][k][j] + c3 * d->dx1[i][k]);
					vy = (c1 * d->pdy_cosine[i][k][j] + c3 * d->dy1[i][k]);
					vz = (c1 * d->pdz_cosine[i][k][j] + c3 * d->dz1[i][k]);
					d->sfdx[i][k][njk] += vx;
					d->sfdy[i][k][njk] += vy;
					d->sfdz[i][k][njk] += vz;
					//d->sfdx_i[i][k][njk] -= vx;
					//d->sfdy_i[i][k][njk] -= vy;
					//d->sfdz_i[i][k][njk] -= vz;
					//d->sfdx[i][k][njk] += (c1 * d->pdx_cosine[i][k][j] + c3 * d->dx1[i][k]);
					//d->sfdy[i][k][njk] += (c1 * d->pdy_cosine[i][k][j] + c3 * d->dy1[i][k]);
					//d->sfdz[i][k][njk] += (c1 * d->pdz_cosine[i][k][j] + c3 * d->dz1[i][k]);
				}
			}
		}
	}

	//printf("sfdx = \n");
	//for (i=0; i<d->n_symmetry; i++) {
	//	for (j=0; j<d->max_natom_p; j++) {
	//		printf("%10d%10d%15.7f%15.7f%15.7f\n", i, j, d->sfdx[i][0][j], d->sfdy[i][0][j], d->sfdz[i][0][j]);
	//	}
	//}
}