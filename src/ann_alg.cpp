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
#include "ann_alg.hpp"
using namespace std;

void ann_alg::zscore_normalization(ann_data *d) {
	for (int i=0; i<d->n_atproc; i++) {
		for (int j=0; j<d->n_symmetry; j++) {
			d->sfv[i][j] = (d->sfv[i][j] - d->sf_mean) / d->sf_std;
		}
	}
}

void ann_alg::forward_propagation(ann_data *d) {
	int i, ii, j, jj, k, l, m, ei;
	double sumf, v;
	
	d->t_energy = 0.;
	
	// forward propagation
	for (i=0; i<d->n_atproc; i++) {
		ei = d->e[i + d->n_start];

		for (j=0; j<d->n_hidden1; j++) {
			sumf = 0.;
			
			for (k=0; k<d->n_symmetry; k++) {
				sumf += d->sfv[i][k] * d->weight0[ei][k][j];
			}
			
			d->out1[i][j] = tanh(sumf + d->bias0[ei][j]);
			d->d_out1[i][j] = 1. - d->out1[i][j] * d->out1[i][j];
		}
		
		for (j=0; j<d->n_hidden2; j++) {
			sumf = 0.;
			
			for (k=0; k<d->n_hidden1; k++) {
				sumf += d->out1[i][k] * d->weight1[ei][k][j];
			}
			
			d->out2[i][j] = tanh(sumf + d->bias1[ei][j]);
			d->d_out2[i][j] = 1. - d->out2[i][j] * d->out2[i][j];
		}
		
		for (j=0; j<d->n_hidden2; j++) {
			d->t_energy += d->out2[i][j] * d->weight2[ei][j];
		}
		d->t_energy += d->bias2[ei];
	}

	////////////////////////////////////////////	
	for (i=0; i<d->natom; i++) {
		d->t_fcx[i] = 0.;
		d->t_fcy[i] = 0.;
		d->t_fcz[i] = 0.;
	}

	for (i=0; i<d->n_atproc; i++) {
		ei = d->e[i + d->n_start];

		for (k=0; k<d->n_hidden1; k++) {
			d->w1w2_dout2[k] = 0.;
			for (l=0; l<d->n_hidden2; l++) {
				d->w1w2_dout2[k] += d->w1w2[ei][k][l] * d->d_out2[i][l];
			}
			d->w1w2_dout2[k] *= d->d_out1[i][k];
		}

		for (k=0; k<d->n_symmetry; k++) {
			d->dedg[i][k] = 0.;
					
			for (l=0; l<d->n_hidden1; l++) {
				d->dedg[i][k] += d->w1w2_dout2[l] * d->weight0[ei][k][l];
				//d->dedg[i][k] += d->w1w2_dout2[l] * d->weight0[ei][k][l] * d->d_out1[i][l];
				//v = 0.;
				//for (m=0; m<d->n_hidden2; m++) {
				//	v += d->w1w2[ei][l][m] * d->d_out2[i][m];
				//}
				//d->dedg[i][k] += v * d->weight0[ei][k][l] * d->d_out1[i][l];
			}
		}
	}
	
	for (i=0; i<d->n_atproc; i++) {
		ii = i + d->n_start;
		//ei = d->e[i + d->n_start];
		
		/*
		for (k=0; k<d->n_symmetry; k++) {
			d->t_fcx[ii] -= d->dedg[i][k] * d->sfdx[i][d->max_natom][k];
			d->t_fcy[ii] -= d->dedg[i][k] * d->sfdy[i][d->max_natom][k];
			d->t_fcz[ii] -= d->dedg[i][k] * d->sfdz[i][d->max_natom][k];
		}
		*/
		
		for (j=0; j<d->cf_natom2[i]; j++) {
			jj = d->cf_pair1[i][j];
			
			for (k=0; k<d->n_symmetry; k++) {
				d->t_fcx[ii] += d->dedg[i][k] * d->sfdx[i][j][k];
				d->t_fcy[ii] += d->dedg[i][k] * d->sfdy[i][j][k];
				d->t_fcz[ii] += d->dedg[i][k] * d->sfdz[i][j][k];
				
				d->t_fcx[jj] -= d->dedg[i][k] * d->sfdx[i][j][k];
				d->t_fcy[jj] -= d->dedg[i][k] * d->sfdy[i][j][k];
				d->t_fcz[jj] -= d->dedg[i][k] * d->sfdz[i][j][k];
			}
		}
	}
	
	for (i=0; i<d->natom; i++) {
		d->t_fcx[i] /= d->sf_std;
		d->t_fcy[i] /= d->sf_std;
		d->t_fcz[i] /= d->sf_std;
	}
	//for (i=0; i<d->n_atproc; i++) {	
	//	if (d->rank == 0) printf("%10d%20.10f%20.10f%20.10f\n", i, d->t_fcx[i], d->t_fcy[i], d->t_fcz[i]);
	//}
}

void ann_alg::forward_propagation_stress(ann_data *d) {
	int i, ii, j, jj, k, l, ei;
	double sumf;
	double fx2b, fy2b, fz2b;
	double ftx, fty, ftz;

	//double sf_std_c = 0.5 / d->sf_std;

	d->t_energy = 0.;
	for (i=0; i<6; i++) d->t_vstress[i] = 0.;

	// forward propagation
	for (i=0; i<d->n_atproc; i++) {
		ei = d->e[i + d->n_start];

		for (j=0; j<d->n_hidden1; j++) {
			sumf = 0.;
			
			for (k=0; k<d->n_symmetry; k++) {
				sumf += d->sfv[i][k] * d->weight0[ei][k][j];
			}
			
			d->out1[i][j] = tanh(sumf + d->bias0[ei][j]);
			d->d_out1[i][j] = 1. - d->out1[i][j] * d->out1[i][j];
		}
		
		for (j=0; j<d->n_hidden2; j++) {
			sumf = 0.;
			
			for (k=0; k<d->n_hidden1; k++) {
				sumf += d->out1[i][k] * d->weight1[ei][k][j];
			}
			
			d->out2[i][j] = tanh(sumf + d->bias1[ei][j]);
			d->d_out2[i][j] = 1. - d->out2[i][j] * d->out2[i][j];
		}
		
		sumf = d->bias2[ei];
		for (j=0; j<d->n_hidden2; j++) {
			//d->t_energy += d->out2[i][j] * d->weight2[ei][j];
			sumf += d->out2[i][j] * d->weight2[ei][j];
		}

		d->t_energy += sumf;
	}

	////////////////////////////////////////////	
	for (i=0; i<d->natom; i++) {
		d->t_fcx[i] = 0.;
		d->t_fcy[i] = 0.;
		d->t_fcz[i] = 0.;

		//for (j=0; j<d->natom; j++) {
		//	d->t_fcx_ij[i][j] = 0.;
		//	d->t_fcy_ij[i][j] = 0.;
		//	d->t_fcz_ij[i][j] = 0.;
		//}
	}

	for (i=0; i<d->n_atproc; i++) {
		ei = d->e[i + d->n_start];

		for (k=0; k<d->n_hidden1; k++) {
			d->w1w2_dout2[k] = 0.;
			for (l=0; l<d->n_hidden2; l++) {
				d->w1w2_dout2[k] += d->w1w2[ei][k][l] * d->d_out2[i][l];
			}
			d->w1w2_dout2[k] *= d->d_out1[i][k];
		}

		for (k=0; k<d->n_symmetry; k++) {
			d->dedg[i][k] = 0.;
					
			for (l=0; l<d->n_hidden1; l++) {
				d->dedg[i][k] += d->w1w2_dout2[l] * d->weight0[ei][k][l];
				//d->dedg[i][k] += d->w1w2_dout2[l] * d->weight0[ei][k][l] * d->d_out1[i][l];
				//v = 0.;
				//for (m=0; m<d->n_hidden2; m++) {
				//	v += d->w1w2[ei][l][m] * d->d_out2[i][m];
				//}
				//d->dedg[i][k] += v * d->weight0[ei][k][l] * d->d_out1[i][l];
			}
		}
	}
	
	for (i=0; i<d->n_atproc; i++) {
		ii = i + d->n_start;
		//ei = d->e[i + d->n_start];
		
		//for (k=0; k<d->n_symmetry; k++) {
		//	d->t_fcx[ii] -= d->dedg[i][k] * d->sfdx[i][d->max_natom][k];
		//	d->t_fcy[ii] -= d->dedg[i][k] * d->sfdy[i][d->max_natom][k];
		//	d->t_fcz[ii] -= d->dedg[i][k] * d->sfdz[i][d->max_natom][k];
		//}

		for (j=0; j<d->cf_natom2[i]; j++) {
			jj = d->cf_pair1[i][j];
			fx2b = 0.;
			fy2b = 0.;
			fz2b = 0.;

			for (k=0; k<d->n_symmetry; k++) {
				/*
				d->t_fcx[ii] -= d->dedg[i][k] * d->sfdx_i[i][j][k];
				d->t_fcy[ii] -= d->dedg[i][k] * d->sfdy_i[i][j][k];
				d->t_fcz[ii] -= d->dedg[i][k] * d->sfdz_i[i][j][k];

				fx2b -= d->dedg[i][k] * d->sfdx_i[i][j][k];
				fy2b -= d->dedg[i][k] * d->sfdy_i[i][j][k];
				fz2b -= d->dedg[i][k] * d->sfdz_i[i][j][k];
				*/

				ftx = d->dedg[i][k] * d->sfdx[i][j][k];
				fty = d->dedg[i][k] * d->sfdy[i][j][k];
				ftz = d->dedg[i][k] * d->sfdz[i][j][k];

				d->t_fcx[ii] += ftx;
				d->t_fcy[ii] += fty;
				d->t_fcz[ii] += ftz;

				d->t_fcx[jj] -= ftx;
				d->t_fcy[jj] -= fty;
				d->t_fcz[jj] -= ftz;

				fx2b += ftx;
				fy2b += fty;
				fz2b += ftz;

				//d->t_fcx[ii] += d->dedg[i][k] * d->sfdx[i][j][k];
				//d->t_fcy[ii] += d->dedg[i][k] * d->sfdy[i][j][k];
				//d->t_fcz[ii] += d->dedg[i][k] * d->sfdz[i][j][k];

				//d->t_fcx[jj] -= d->dedg[i][k] * d->sfdx[i][j][k];
				//d->t_fcy[jj] -= d->dedg[i][k] * d->sfdy[i][j][k];
				//d->t_fcz[jj] -= d->dedg[i][k] * d->sfdz[i][j][k];

				//fx2b += d->dedg[i][k] * d->sfdx[i][j][k];
				//fy2b += d->dedg[i][k] * d->sfdy[i][j][k];
				//fz2b += d->dedg[i][k] * d->sfdz[i][j][k];
			}
			fx2b /= d->sf_std;
			fy2b /= d->sf_std;
			fz2b /= d->sf_std;

			// two-body force
			//d->t_fcx_ij[ii][jj] += fx2b;
			//d->t_fcy_ij[ii][jj] += fy2b;
			//d->t_fcz_ij[ii][jj] += fz2b;
			//d->t_fcx_ij[jj][ii] -= fx2b;
			//d->t_fcy_ij[jj][ii] -= fy2b;
			//d->t_fcz_ij[jj][ii] -= fz2b;

			//fx2b *= 0.5;
			//fy2b *= 0.5;
			//fz2b *= 0.5;

			d->t_vstress[0] -= fx2b * d->dx1[i][j];     // Sxx
			d->t_vstress[1] -= fy2b * d->dy1[i][j];     // Syy
			d->t_vstress[2] -= fz2b * d->dz1[i][j];     // Szz 
			d->t_vstress[3] -= 0.5*(fz2b * d->dy1[i][j] + fy2b * d->dz1[i][j]);     // Syz
			d->t_vstress[4] -= 0.5*(fz2b * d->dx1[i][j] + fx2b * d->dz1[i][j]);     // Sxz
			d->t_vstress[5] -= 0.5*(fy2b * d->dx1[i][j] + fx2b * d->dy1[i][j]);     // Sxy
		}

		/*
		for (j=0; j<d->cf_natom2[i]; j++) {
			jj = d->cf_pair1[i][j];
			fx2b = 0.;
			fy2b = 0.;
			fz2b = 0.;
			
			for (k=0; k<d->n_symmetry; k++) {
				d->t_fcx[jj] -= d->dedg[i][k] * d->sfdx[i][j][k];
				d->t_fcy[jj] -= d->dedg[i][k] * d->sfdy[i][j][k];
				d->t_fcz[jj] -= d->dedg[i][k] * d->sfdz[i][j][k];

				fx2b -= d->dedg[i][k] * d->sfdx[i][j][k];
				fy2b -= d->dedg[i][k] * d->sfdy[i][j][k];
				fz2b -= d->dedg[i][k] * d->sfdz[i][j][k];
			}
			fx2b /= d->sf_std;
			fy2b /= d->sf_std;
			fz2b /= d->sf_std;

			d->t_fcx_ij[jj][ii] += fx2b;
			d->t_fcy_ij[jj][ii] += fy2b;
			d->t_fcz_ij[jj][ii] += fz2b;

			fx2b *= 0.5;
			fy2b *= 0.5;
			fz2b *= 0.5;

			d->t_vstress[0] += fx2b * d->dx1[i][j];     // Sxx
			d->t_vstress[1] += fy2b * d->dy1[i][j];     // Syy
			d->t_vstress[2] += fz2b * d->dz1[i][j];     // Szz 
			d->t_vstress[3] += 0.5*(fz2b * d->dy1[i][j] + fy2b * d->dz1[i][j]);     // Syz
			d->t_vstress[4] += 0.5*(fz2b * d->dx1[i][j] + fx2b * d->dz1[i][j]);     // Sxz
			d->t_vstress[5] += 0.5*(fy2b * d->dx1[i][j] + fx2b * d->dy1[i][j]);     // Sxy
		}
		*/
	}

	for (i=0; i<d->natom; i++) {
		d->t_fcx[i] /= d->sf_std;
		d->t_fcy[i] /= d->sf_std;
		d->t_fcz[i] /= d->sf_std;
	}

	//for (i=0; i<6; i++) printf("t_vstress[%d]=%15.7f\n", i, d->t_vstress[i]);
}

void ann_alg::forward_propagation_pmd(ann_data *d) {
	int i, ii, j, jj, k, l, ei;
	double sumf;
	double vxx, vyy, vzz, vyz, vxz, vxy;
	double fx2b, fy2b, fz2b;
	double ftx, fty, ftz;

	//double sf_std_c = 0.5 / d->sf_std;

	d->t_energy = 0.;
	for (i=0; i<6; i++) d->t_vstress[i] = 0.;

	// forward propagation
	for (i=0; i<d->n_atproc; i++) {
		ei = d->e[i + d->n_start];

		for (j=0; j<d->n_hidden1; j++) {
			sumf = 0.;
			
			for (k=0; k<d->n_symmetry; k++) {
				sumf += d->sfv[i][k] * d->weight0[ei][k][j];
			}
			
			d->out1[i][j] = tanh(sumf + d->bias0[ei][j]);
			d->d_out1[i][j] = 1. - d->out1[i][j] * d->out1[i][j];
		}
		
		for (j=0; j<d->n_hidden2; j++) {
			sumf = 0.;
			
			for (k=0; k<d->n_hidden1; k++) {
				sumf += d->out1[i][k] * d->weight1[ei][k][j];
			}
			
			d->out2[i][j] = tanh(sumf + d->bias1[ei][j]);
			d->d_out2[i][j] = 1. - d->out2[i][j] * d->out2[i][j];
		}
		
		sumf = d->bias2[ei];
		for (j=0; j<d->n_hidden2; j++) {
			//d->t_energy += d->out2[i][j] * d->weight2[ei][j];
			sumf += d->out2[i][j] * d->weight2[ei][j];
		}

		d->p_at_en[i + d->n_start] = sumf;
		d->t_energy += sumf;
		//d->t_energy += d->bias2[ei];
	}

	////////////////////////////////////////////	
	for (i=0; i<d->natom; i++) {
		d->t_fcx[i] = 0.;
		d->t_fcy[i] = 0.;
		d->t_fcz[i] = 0.;
		
        d->p_at_stress[i][0] = 0.;
        d->p_at_stress[i][1] = 0.;
        d->p_at_stress[i][2] = 0.;
        d->p_at_stress[i][3] = 0.;
        d->p_at_stress[i][4] = 0.;
        d->p_at_stress[i][5] = 0.;
	}

	for (i=0; i<d->n_atproc; i++) {
		ei = d->e[i + d->n_start];

		for (k=0; k<d->n_hidden1; k++) {
			d->w1w2_dout2[k] = 0.;
			for (l=0; l<d->n_hidden2; l++) {
				d->w1w2_dout2[k] += d->w1w2[ei][k][l] * d->d_out2[i][l];
			}
			d->w1w2_dout2[k] *= d->d_out1[i][k];
		}

		for (k=0; k<d->n_symmetry; k++) {
			d->dedg[i][k] = 0.;
					
			for (l=0; l<d->n_hidden1; l++) {
				d->dedg[i][k] += d->w1w2_dout2[l] * d->weight0[ei][k][l];
				//d->dedg[i][k] += d->w1w2_dout2[l] * d->weight0[ei][k][l] * d->d_out1[i][l];
				//v = 0.;
				//for (m=0; m<d->n_hidden2; m++) {
				//	v += d->w1w2[ei][l][m] * d->d_out2[i][m];
				//}
				//d->dedg[i][k] += v * d->weight0[ei][k][l] * d->d_out1[i][l];
			}
		}
	}
	
	for (i=0; i<d->n_atproc; i++) {
		ii = i + d->n_start;
		
		//for (k=0; k<d->n_symmetry; k++) {
		//	d->t_fcx[ii] -= d->dedg[i][k] * d->sfdx[i][d->max_natom][k];
		//	d->t_fcy[ii] -= d->dedg[i][k] * d->sfdy[i][d->max_natom][k];
		//	d->t_fcz[ii] -= d->dedg[i][k] * d->sfdz[i][d->max_natom][k];
		//}
		for (j=0; j<d->cf_natom2[i]; j++) {
			jj = d->cf_pair1[i][j];
			fx2b = 0.;
			fy2b = 0.;
			fz2b = 0.;

			for (k=0; k<d->n_symmetry; k++) {
				//d->t_fcx[ii] -= d->dedg[i][k] * d->sfdx_i[i][j][k];
				//d->t_fcy[ii] -= d->dedg[i][k] * d->sfdy_i[i][j][k];
				//d->t_fcz[ii] -= d->dedg[i][k] * d->sfdz_i[i][j][k];

				//fx2b -= d->dedg[i][k] * d->sfdx_i[i][j][k];
				//fy2b -= d->dedg[i][k] * d->sfdy_i[i][j][k];
				//fz2b -= d->dedg[i][k] * d->sfdz_i[i][j][k];

				ftx = d->dedg[i][k] * d->sfdx[i][j][k];
				fty = d->dedg[i][k] * d->sfdy[i][j][k];
				ftz = d->dedg[i][k] * d->sfdz[i][j][k];

				d->t_fcx[ii] += ftx;
				d->t_fcy[ii] += fty;
				d->t_fcz[ii] += ftz;

				d->t_fcx[jj] -= ftx;
				d->t_fcy[jj] -= fty;
				d->t_fcz[jj] -= ftz;

				fx2b += ftx;
				fy2b += fty;
				fz2b += ftz;
			}
			fx2b /= d->sf_std;
			fy2b /= d->sf_std;
			fz2b /= d->sf_std;

			vxx = fx2b * d->dx1[i][j];
			vyy = fy2b * d->dy1[i][j];
			vzz = fz2b * d->dz1[i][j];
			vyz = 0.5*(fz2b * d->dy1[i][j] + fy2b * d->dz1[i][j]);
			vxz = 0.5*(fz2b * d->dx1[i][j] + fx2b * d->dz1[i][j]);
			vxy = 0.5*(fy2b * d->dx1[i][j] + fx2b * d->dy1[i][j]);

			//d->p_at_stress[ii][0] -= 0.5*vxx;     // Sxx
			//d->p_at_stress[ii][1] -= 0.5*vyy;     // Syy
			//d->p_at_stress[ii][2] -= 0.5*vzz;     // Szz 
			//d->p_at_stress[ii][3] -= 0.5*vyz;     // Syz
			//d->p_at_stress[ii][4] -= 0.5*vxz;     // Sxz
			//d->p_at_stress[ii][5] -= 0.5*vxy;     // Sxy

			d->p_at_stress[jj][0] -= vxx;     // Sxx
			d->p_at_stress[jj][1] -= vyy;     // Syy
			d->p_at_stress[jj][2] -= vzz;     // Szz 
			d->p_at_stress[jj][3] -= vyz;     // Syz
			d->p_at_stress[jj][4] -= vxz;     // Sxz
			d->p_at_stress[jj][5] -= vxy;     // Sxy

			d->t_vstress[0] -= vxx;     // Sxx
			d->t_vstress[1] -= vyy;     // Syy
			d->t_vstress[2] -= vzz;     // Szz 
			d->t_vstress[3] -= vyz;     // Syz
			d->t_vstress[4] -= vxz;     // Sxz
			d->t_vstress[5] -= vxy;     // Sxy
		}

		/*
		for (j=0; j<d->cf_natom2[i]; j++) {
			jj = d->cf_pair1[i][j];
			fx2b = 0.;
			fy2b = 0.;
			fz2b = 0.;
			
			for (k=0; k<d->n_symmetry; k++) {
				d->t_fcx[jj] -= d->dedg[i][k] * d->sfdx[i][j][k];
				d->t_fcy[jj] -= d->dedg[i][k] * d->sfdy[i][j][k];
				d->t_fcz[jj] -= d->dedg[i][k] * d->sfdz[i][j][k];

				fx2b -= d->dedg[i][k] * d->sfdx[i][j][k];
				fy2b -= d->dedg[i][k] * d->sfdy[i][j][k];
				fz2b -= d->dedg[i][k] * d->sfdz[i][j][k];
			}
			fx2b /= d->sf_std;
			fy2b /= d->sf_std;
			fz2b /= d->sf_std;

			vxx = fx2b * d->dx1[i][j];
			vyy = fy2b * d->dy1[i][j];
			vzz = fz2b * d->dz1[i][j];
			vyz = 0.5*(fz2b * d->dy1[i][j] + fy2b * d->dz1[i][j]);
			vxz = 0.5*(fz2b * d->dx1[i][j] + fx2b * d->dz1[i][j]);
			vxy = 0.5*(fy2b * d->dx1[i][j] + fx2b * d->dy1[i][j]);
			
			d->p_at_stress[jj][0] += 0.5*vxx;     // Sxx
			d->p_at_stress[jj][1] += 0.5*vyy;     // Syy
			d->p_at_stress[jj][2] += 0.5*vzz;     // Szz 
			d->p_at_stress[jj][3] += 0.5*vyz;     // Syz
			d->p_at_stress[jj][4] += 0.5*vxz;     // Sxz
			d->p_at_stress[jj][5] += 0.5*vxy;     // Sxy

			d->t_vstress[0] += 0.5*vxx;     // Sxx
			d->t_vstress[1] += 0.5*vyy;     // Syy
			d->t_vstress[2] += 0.5*vzz;     // Szz 
			d->t_vstress[3] += 0.5*vyz;     // Syz
			d->t_vstress[4] += 0.5*vxz;     // Sxz
			d->t_vstress[5] += 0.5*vxy;     // Sxy
		}
		*/
	}

	for (i=0; i<d->natom; i++) {
		d->t_fcx[i] /= d->sf_std;
		d->t_fcy[i] /= d->sf_std;
		d->t_fcz[i] /= d->sf_std;
	}
	//for (i=0; i<6; i++) printf("t_vstress[%d]=%15.7f\n", i, d->t_vstress[i]);
}
