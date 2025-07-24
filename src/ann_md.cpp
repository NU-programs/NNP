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
#include "ann_md.hpp"
using namespace std;

vector<string> ann_md::split (string line) {
    istringstream iss(line);
    string s;
    vector<string> ss;
    while (iss >> s) ss.push_back(s);
    return ss;
}

void ann_md::md_init(ann_data *d) {
	int i;
	double cvx, cvy, cvz;
	double total_mass;

	d->natom_moved = 0;
	d->natom_fixed = 0;
	d->total_mass = 0.;
	
	if (d->sel_flag == "yes") {
		for (i=0; i<d->natom; i++) {
			if (d->mflx[i] == 1) {
				d->natom_moved += 1;
				d->total_mass += d->mass[i];
			} else {
				d->natom_fixed += 1;
			}
		}
	} else {
		d->natom_moved = d->natom;
		for (i=0; i<d->natom; i++) d->total_mass += d->mass[i];
	}
	
	if (d->init_vel == "init") {
		if (d->sel_flag == "yes") {
			for (i=0; i<d->natom; i++) {
				if (d->mflx[i] == 1) {
					d->vx[i] = 0.1 * rand() / RAND_MAX - 0.1;
					d->vy[i] = 0.1 * rand() / RAND_MAX - 0.1;
					d->vz[i] = 0.1 * rand() / RAND_MAX - 0.1;
				} else {
					d->vx[i] = 0.;
					d->vy[i] = 0.;
					d->vz[i] = 0.;
				}
			}
		} else {
			for (i=0; i<d->natom; i++) {
				d->vx[i] = 0.1 * rand() / RAND_MAX - 0.1;
				d->vy[i] = 0.1 * rand() / RAND_MAX - 0.1;
				d->vz[i] = 0.1 * rand() / RAND_MAX - 0.1;
			}
		}
	
		d->ke = 0;
		for (i=0; i<d->natom; i++) {
			d->ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
		}
		d->ke *= 0.5 * d->c_unit2;
		d->c_temp = d->ke / (1.5 * (double)d->natom_moved * kB);
		//printf("ke [eV], c_temp [K] = %20.10f%20.10f\n", d->ke, d->c_temp);
			
		//d->md_vscale = sqrt((d->md_temp + (d->c_temp - d->md_temp)*1e-4) / d->c_temp);
		d->md_vscale = sqrt(d->md_temp / d->c_temp);
		//printf("md_vscale = %20.10f\n", d->md_vscale);
			
		for (i=0; i<d->natom; i++) {
			d->vx[i] *= d->md_vscale;
			d->vy[i] *= d->md_vscale;
			d->vz[i] *= d->md_vscale;
		}

	} else if (d->init_vel == "read") {
		d->ke = 0;
		for (i=0; i<d->natom; i++) {
			d->ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
		}
		d->ke *= 0.5 * d->c_unit2;
		d->c_temp = d->ke / (1.5 * (double)d->natom_moved * kB);
		//printf("ke [eV], c_temp [K] = %20.10f%20.10f\n", d->ke, d->c_temp);
	}
	
	/*
	d->ke = 0;
	for (i=0; i<d->natom; i++) {
		d->ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->ke *= 0.5 * d->c_unit2;
	d->c_temp = d->ke / (1.5 * (double)d->natom * kB);
	*/
	//printf("ke [eV], c_temp [K] = %20.10f%20.10f\n", d->ke, d->c_temp);
}

void ann_md::allocation_pmd(ann_data *d) {
	int i, j;

	//cout << "allocation_pmd" << endl;

	// allocation of atomic energy
	d->p_at_en = new double[d->natom];
	d->at_en = new double[d->natom];
	for (i=0; i<d->natom; i++) {
		d->p_at_en[i] = 0.;
		d->at_en[i] = 0.;
	}

	d->at_stress = new double*[d->natom];
	d->p_at_stress = new double*[d->natom];
	d->at_stress[0] = new double[d->natom * 6];
	d->p_at_stress[0] = new double[d->natom * 6];
	for (i=0; i<d->natom; i++) {
		d->at_stress[i] = d->at_stress[0] + i*6;
		d->p_at_stress[i] = d->p_at_stress[0] + i*6;
	}

	for (i=0; i<d->natom; i++) {
		for (j=0; j<6; j++) {
			d->at_stress[i][j] = 0.;
			d->p_at_stress[i][j] = 0.;
			//if (d->rank == 0) {
			//	printf("%10d%10d%15.7f%15.7f\n", i, j, d->p_at_stress[i][j], d->at_stress[i][j]);
			//}
		}
	}

    d->d_tensor = new double*[d->natom];
    d->d_tensor[0] = new double[d->natom*6];
    for (i=1; i < d->natom; i++) {
        d->d_tensor[i] = d->d_tensor[0] + i * 6;
    }

	if (d->rank == 0) {
		for (i=0; i<d->natom; i++) {
			for (j=0; j<6; j++) d->d_tensor[i][j] = 0.;
		}

		for (i=0; i<6; i++) d->d_average[i] = 0.;

		d->at_hflux_avg = new double*[d->natom];
		for (i=0; i<d->natom; i++) {
			d->at_hflux_avg[i] = new double[9];
			for (j=0; j<9; j++) d->at_hflux_avg[i][j] = 0.;
		}
    }

    if (d->mode_decomp == "yes") {

        d->nmode = 3 * d->natom;
        d->local_eigenvectors  = new double[d->n_atproc*3*d->natom*3];
        d->local_at_hflux_mode_avg = new double[d->n_atproc*3*d->natom*4];

        double* eigenvectors;
        if (d->rank == 0) {
            d->frequencies = new double[d->nmode];
            eigenvectors = new double[d->nmode*d->natom*3];
            d->at_hflux_mode_avg = new double[d->nmode*d->natom*4];

            ifstream fin;
            string line;
            vector<string> s;
    
            fin.open(d->mode_infile.c_str(), ios::in);
            if (!fin) {
                cout << d->mode_infile << " cannot be opened." << endl;
                exit(1);
            }

            cout << "Now reading eigenvectors and frequencies..." << flush;
            for (int m=0; m<d->nmode; m++){
                getline(fin, line); s = split(line);
                d->frequencies[m] = atof(s[0].c_str());
                for (i=0; i<d->natom; i++){
                    getline(fin, line); s = split(line);
                    eigenvectors[(m*d->natom+i)*3+0] = atof(s[0].c_str()) / sqrt(d->mass[i]);
                    eigenvectors[(m*d->natom+i)*3+1] = atof(s[1].c_str()) / sqrt(d->mass[i]);
                    eigenvectors[(m*d->natom+i)*3+2] = atof(s[2].c_str()) / sqrt(d->mass[i]);
                }
            }
            fin.close();
            cout << "finished" << endl;
        }

        MPI_Bcast(&d->mass[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // for eigenvectors
        d->all_n_mode_proc  = new int[d->nprocs];
        d->all_n_mode_start = new int[d->nprocs];
        if (d->rank == 0) {
            for (int i=0; i<d->nprocs; i++) {
                d->all_n_mode_proc[i] = d->all_n_atproc[i] * 3 * d->natom * 3;
            }
            d->all_n_mode_start[0] = 0;
            for (int i=1; i<d->nprocs; i++) {
                d->all_n_mode_start[i] = d->all_n_mode_start[i-1] + d->all_n_mode_proc[i-1];
            }
        }

        MPI_Scatterv(
            eigenvectors,
            d->all_n_mode_proc,
            d->all_n_mode_start,
            MPI_DOUBLE,
            d->local_eigenvectors,
            d->n_atproc * 3 * d->natom * 3,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        if (d->rank == 0) delete[] eigenvectors;

        // for local_at_hflux_avg
        if (d->rank == 0) {
            for (int i=0; i<d->nprocs; i++) {
                d->all_n_mode_proc[i] = d->all_n_atproc[i] * 3 * d->natom * 4;
            }
            d->all_n_mode_start[0] = 0;
            for (int i=1; i<d->nprocs; i++) {
                d->all_n_mode_start[i] = d->all_n_mode_start[i-1] + d->all_n_mode_proc[i-1];
            }
        }
    }
}

void ann_md::set_md_const(ann_data *d) {
	d->c_unit1 = 1.602177 * 1e-2 / 1.660539; // [eV->J] * [m^2->ang^2] * [s^2->fs^2] / [atomic_mass->Kg]
	d->c_unit2 = 1.660539 * 6.24151 * 10.; // [atomic_mass->Kg] * [fs^2->s^2] * [ang^2->m^2] * [J->eV]

	d->md_dtemp = (d->md_fin_temp - d->md_init_temp) / d->md_step;
	d->md_temp = d->md_init_temp;

	d->fric = 0.;
	d->th_mass = (3.*d->natom - 3.) * (kB * d->md_temp) * d->th_mass * d->th_mass;

	d->md_time_h = 0.5 * d->md_time;
	d->md_time2_h = 0.5 * d->md_time * d->md_time;
	d->t_ke = 0.5 * (3.*d->natom + 1.) * (kB * d->md_temp);

	d->lv_sigma = sqrt(2. * kB * d->md_temp * d->md_friction * d->c_unit1);
	d->time_h_fric = d->md_time_h * d->md_friction;

	d->md_const1 = d->md_time2_h * d->c_unit1;
	d->md_const2 = d->md_time_h * d->c_unit1;
	//d->md_const3 = 1. / (3. * d->natom * kB) * d->c_unit2;
	
	d->md_const4 = 1./8. * d->md_time * d->md_time * d->md_friction;
	d->md_const4_fric = d->md_const4 * d->md_friction;
	d->md_const5 = pow(d->md_time, 1.5) * d->lv_sigma / (2.*sqrt(3.));
	d->md_const6 = 1./4. * pow(d->md_time, 1.5) * d->md_friction * d->lv_sigma;

	//d->mass = new double[d->natom];
	//d->mass_unit1 = new double[d->natom];
	d->mass_const1 = new double[d->natom];
	d->mass_const2 = new double[d->natom];
	d->mass_const3 = new double[d->natom];
	d->mass_const4 = new double[d->natom];
	d->mass_const5 = new double[d->natom];
	d->mass_const6 = new double[d->natom];
	
	for (int i=0; i<d->natom; i++) {
		d->mass[i] = d->md_mass[d->e[i]];
		d->mass_const1[i] = d->md_const1 / d->mass[i];
		d->mass_const2[i] = d->md_const2 / d->mass[i];
		d->mass_const3[i] = 0.5 * sqrt(d->md_time) * d->lv_sigma / sqrt(d->mass[i]);
		d->mass_const4[i] = d->md_const4 * d->c_unit1 / d->mass[i];
		d->mass_const5[i] = d->md_const5 / sqrt(d->mass[i]);
		d->mass_const6[i] = d->md_const6 / sqrt(d->mass[i]);
		//printf("mass = %20.10f%20.10f%20.10f\n", d->mass[i], d->mass_const1[i], d->mass_const2[i]);
	}
}

void ann_md::zero_velocity_center_of_mass(ann_data *d) {
	int i;
	double cvx, cvy, cvz;
	
	// set center of velocity to zero
	cvx = 0.;
	cvy = 0.;
	cvz = 0.;
	for (i=0; i<d->natom; i++) {
		cvx += d->mass[i] * d->vx[i];
		cvy += d->mass[i] * d->vy[i];
		cvz += d->mass[i] * d->vz[i];
	}
		
	cvx /= d->total_mass;
	cvy /= d->total_mass;
	cvz /= d->total_mass;
		
	//printf("center_of_velocity = %15.7f%15.7f%15.7f\n", cvx, cvy, cvz);
	
	if (d->sel_flag == "yes") {
		for (i=0; i<d->natom; i++) {
			if (d->mflx[i] == 1) {
				d->vx[i] -= cvx;
				d->vy[i] -= cvy;
				d->vz[i] -= cvz;
			}
		}
	} else {
		for (i=0; i<d->natom; i++) {
			d->vx[i] -= cvx;
			d->vy[i] -= cvy;
			d->vz[i] -= cvz;
		}
	}
}

void ann_md::allocation_velauto(ann_data *d) {
	int i, j, k;
	
	d->tva = new double[d->velauto_step];
	for (i=0; i<d->velauto_step; i++) d->tva[i] = 0.;

	d->va0 = new double[d->velauto_n];
	for (i=0; i<d->velauto_n; i++) d->va0[i] = 0.;

	if (d->velauto_at_flag == "yes") {
		d->at_tva = new double*[d->velauto_step];
		for (i=0; i<d->velauto_step; i++) {
			d->at_tva[i] = new double[d->natom];
			for (j=0; j<d->natom; j++) {
				d->at_tva[i][j] = 0.;
			}
		}

		d->at_va0 = new double*[d->velauto_n];
		for (i=0; i<d->velauto_n; i++) {
			d->at_va0[i] = new double[d->natom];
			for (j=0; j<d->natom; j++) {
				d->at_va0[i][j] = 0.;
			}
		}
	}

	d->ivx = new double*[d->velauto_n];
	d->ivy = new double*[d->velauto_n];
	d->ivz = new double*[d->velauto_n];
	for (i=0; i<d->velauto_n; i++) {
		d->ivx[i] = new double[d->natom];
		d->ivy[i] = new double[d->natom];
		d->ivz[i] = new double[d->natom];
		for (j=0; j<d->natom; j++) {
			d->ivx[i][j] = 0.;
			d->ivy[i][j] = 0.;
			d->ivz[i][j] = 0.;
		}
	}
}

void ann_md::calc_velauto(ann_data *d, int step) {
	int i, j, n, ns;
	double v, vv;

	(step-1) / d->velauto_inv < (d->velauto_n-1) ? ns = (step-1) / d->velauto_inv : ns = d->velauto_n-1;

	if ( ((step-1) % d->velauto_inv == 0) && ((step-1) / d->velauto_inv <= (d->velauto_n-1)) ) {
		if (d->velauto_at_flag == "yes") {
			for (i=0; i<d->natom; i++) {
				d->ivx[ns][i] = d->vx[i];
				d->ivy[ns][i] = d->vy[i];
				d->ivz[ns][i] = d->vz[i];
				d->at_va0[ns][i] = d->ivx[ns][i]*d->ivx[ns][i] + d->ivy[ns][i]*d->ivy[ns][i] + d->ivz[ns][i]*d->ivz[ns][i];
				d->va0[ns] += d->at_va0[ns][i]; 
			}
		} else {
			for (i=0; i<d->natom; i++) {
				d->ivx[ns][i] = d->vx[i];
				d->ivy[ns][i] = d->vy[i];
				d->ivz[ns][i] = d->vz[i];
				//d->at_va0[ns][i] = d->ivx[ns][i]*d->ivx[ns][i] + d->ivy[ns][i]*d->ivy[ns][i] + d->ivz[ns][i]*d->ivz[ns][i];
				d->va0[ns] += d->ivx[ns][i]*d->ivx[ns][i] + d->ivy[ns][i]*d->ivy[ns][i] + d->ivz[ns][i]*d->ivz[ns][i];
			}
		}
	}

	for (i=0; i<=ns; i++) {
		n = (step-1) - i*d->velauto_inv;

		if (n < d->velauto_step) {
			vv = 0.;

			if (d->velauto_at_flag == "yes") {
				for (j=0; j<d->natom; j++) {
					v = d->ivx[i][j]*d->vx[j] + d->ivy[i][j]*d->vy[j] + d->ivz[i][j]*d->vz[j];
					d->at_tva[n][j] += (v / d->at_va0[i][j]);
					vv += v;
					//d->tva[n] += (v / d->va0[i]); 
				}
				d->tva[n] += (vv / d->va0[i]);
			} else {
				for (j=0; j<d->natom; j++) {
					//v = d->ivx[i][j]*d->vx[j] + d->ivy[i][j]*d->vy[j] + d->ivz[i][j]*d->vz[j];
					vv += d->ivx[i][j]*d->vx[j] + d->ivy[i][j]*d->vy[j] + d->ivz[i][j]*d->vz[j];
				}
				d->tva[n] += (vv / d->va0[i]);
			}
		}
	}
}

void ann_md::average_velauto(ann_data *d) {
	int i, j, k;

	for (i=0; i<d->velauto_step; i++) {
		d->tva[i] /= (double)d->velauto_n;
	}

	if (d->velauto_at_flag == "yes") {
		for (i=0; i<d->velauto_step; i++) {
			for (j=0; j<d->natom; j++) {
				d->at_tva[i][j] /= (double)d->velauto_n;
			}
		}
	}
}

void ann_md::output_velauto(ann_data *d) {
	int i, j;
	FILE *fout;
	char cm[256], fn[256]; 

	sprintf(cm, "mkdir %s", d->velauto_outdir.c_str());
	system(cm);

	/// output total va
	sprintf(fn, "%s/total.dat", d->velauto_outdir.c_str());

	if ((fout = fopen(fn, "w")) != NULL) {
	} else {
		printf("fail to open %s.\n", fn);
		exit(1);
	}

	fprintf(fout, "# velauto_n    = %d\n", d->velauto_n);
	fprintf(fout, "# velauto_step = %d\n", d->velauto_step);

	fprintf(fout, "# time total_va\n");
	for (i=0; i<d->velauto_step; i++) {
		fprintf(fout, "%15.3f%15.7f\n", d->md_time*i, d->tva[i]);
	}
	fprintf(fout, "\n");

	fclose(fout);

	if (d->velauto_at_flag == "yes") {
		/// output individual va
		for (i=0; i<d->natom; i++) {
			sprintf(fn, "%s/%05d.dat", d->velauto_outdir.c_str(), i);
		
			if ((fout = fopen(fn, "w")) != NULL) {
			} else {
				printf("fail to open %s.\n", fn);
				exit(1);
			}

			for (j=0; j<d->velauto_step; j++) {
				fprintf(fout, "%15.7f", d->at_tva[j][i]);
			}
		
			fclose(fout);
		}
	}
}

void ann_md::allocation_msd(ann_data *d) {
	int i, j, k;

	d->msdt = new double*[d->msd_step];
	d->msdx = new double*[d->msd_step];
	d->msdy = new double*[d->msd_step];
	d->msdz = new double*[d->msd_step];

	for (i=0; i<d->msd_step; i++) {
		d->msdt[i] = new double[d->natom];
		d->msdx[i] = new double[d->natom];
		d->msdy[i] = new double[d->natom];
		d->msdz[i] = new double[d->natom];
		for (j=0; j<d->natom; j++) {
			d->msdt[i][j] = 0.;
			d->msdx[i][j] = 0.;
			d->msdy[i][j] = 0.;
			d->msdz[i][j] = 0.;
		}
	}

	d->ix = new double*[d->msd_n];
	d->iy = new double*[d->msd_n];
	d->iz = new double*[d->msd_n];
	for (i=0; i<d->msd_n; i++) {
		d->ix[i] = new double[d->natom];
		d->iy[i] = new double[d->natom];
		d->iz[i] = new double[d->natom];
		for (j=0; j<d->natom; j++) {
			d->ix[i][j] = 0.;
			d->iy[i][j] = 0.;
			d->iz[i][j] = 0.;
		}
	}
}

void ann_md::calc_msd(ann_data *d, int step) {
	int i, j, n, ns;
	double dx, dy, dz, dxx, dyy, dzz;

	(step-1) / d->msd_inv < (d->msd_n-1) ? ns = (step-1) / d->msd_inv : ns = d->msd_n-1;

	if ( ((step-1) % d->msd_inv == 0) && ((step-1) / d->msd_inv <= (d->msd_n-1)) ) { 
		for (i=0; i<d->natom; i++) {
			d->ix[ns][i] = d->x[i];
			d->iy[ns][i] = d->y[i];
			d->iz[ns][i] = d->z[i];
		}
	}

	for (i=0; i<=ns; i++) {
		n = (step-1) - i*d->msd_inv;

		if (n < d->msd_step) {
			for (j=0; j<d->natom; j++) {
				dx = d->x[j] - d->ix[i][j];
				if (dx > 0.5) {
					dx -= 1.;
				} else if (dx < -0.5) {
					dx += 1.;
				}

				dy = d->y[j] - d->iy[i][j];
				if (dy > 0.5) {
					dy -= 1.;
				} else if (dy < -0.5) {
					dy += 1.;
				}

				dz = d->z[j] - d->iz[i][j];
				if (dz > 0.5) {
					dz -= 1.;
				} else if (dz < -0.5) {
					dz += 1.;
				}

				dxx = dx*d->lv[0][0] + dy*d->lv[1][0] + dz*d->lv[2][0];
				dyy = dx*d->lv[0][1] + dy*d->lv[1][1] + dz*d->lv[2][1];
				dzz = dx*d->lv[0][2] + dy*d->lv[1][2] + dz*d->lv[2][2];
				dxx = dxx * dxx;
				dyy = dyy * dyy;
				dzz = dzz * dzz;

				//dx = (d->x[j] - d->ix[i][j]) * (d->x[j] - d->ix[i][j]);
				//dy = (d->y[j] - d->iy[i][j]) * (d->y[j] - d->iy[i][j]);
				//dz = (d->z[j] - d->iz[i][j]) * (d->z[j] - d->iz[i][j]);
				d->msdx[n][j] += dxx;
				d->msdy[n][j] += dyy;
				d->msdz[n][j] += dzz;
				d->msdt[n][j] += (dxx + dyy + dzz);
			}
		}
	}
}

void ann_md::average_msd(ann_data *d) {
	int i, j;
	for (i=0; i<d->msd_step; i++) {
		for (j=0; j<d->natom; j++) {
			d->msdx[i][j] = sqrt(d->msdx[i][j] / (double)d->msd_n);
			d->msdy[i][j] = sqrt(d->msdy[i][j] / (double)d->msd_n);
			d->msdz[i][j] = sqrt(d->msdz[i][j] / (double)d->msd_n);
			d->msdt[i][j] = sqrt(d->msdt[i][j] / (double)d->msd_n);
		}
	}
}

void ann_md::output_msd(ann_data *d) {
	int i, j;
	FILE *fout;
	char cm[256], fn[256]; 

	sprintf(cm, "mkdir %s", d->msd_outdir.c_str());
	system(cm);

	for (i=0; i<d->natom; i++) {
		sprintf(fn, "%s/%05d.dat", d->msd_outdir.c_str(), i);
		
		if ((fout = fopen(fn, "w")) != NULL) {
		} else {
			printf("fail to open %s.\n", fn);
			exit(1);
		}

		fprintf(fout, "# md_time, msd_step = %15.7f%10d\n", d->md_time, d->msd_step);
		for (j=0; j<d->msd_step; j++) {
			fprintf(fout, "%15.7f%15.7f%15.7f%15.7f\n", d->msdx[j][i], d->msdy[j][i], d->msdz[j][i], d->msdt[j][i]);
		}
		
		fclose(fout);
	}
}

void ann_md::allocation_average_position(ann_data *d) {
	d->x0 = new double[d->natom];
	d->y0 = new double[d->natom];
	d->z0 = new double[d->natom];
	d->xa = new double[d->natom];
	d->ya = new double[d->natom];
	d->za = new double[d->natom];

	for (int i=0; i<d->natom; i++) {
		d->x0[i] = d->x[i];
		d->y0[i] = d->y[i];
		d->z0[i] = d->z[i];
		d->xa[i] = d->x[i];
		d->ya[i] = d->y[i];
		d->za[i] = d->z[i];
	}
}

void ann_md::calc_average_position(ann_data *d) {
	double dx, dy, dz;

	for (int i=0; i<d->natom; i++) {
		dx = d->x[i] - d->x0[i];
		dy = d->y[i] - d->y0[i];
		dz = d->z[i] - d->z0[i];

		if (dx > 0.5) {
			d->xa[i] += (d->x[i] - 1.0);
		} else if (dx < -0.5) {
			d->xa[i] += (d->x[i] + 1.0);
		} else {
			d->xa[i] += d->x[i];
		}

		if (dy > 0.5) {
			d->ya[i] += (d->y[i] - 1.0);
		} else if (dy < -0.5) {
			d->ya[i] += (d->y[i] + 1.0);
		} else {
			d->ya[i] += d->y[i];
		}

		if (dz > 0.5) {
			d->za[i] += (d->z[i] - 1.0);
		} else if (dz < -0.5) {
			d->za[i] += (d->z[i] + 1.0);
		} else {
			d->za[i] += d->z[i];
		}
	}
}

void ann_md::output_average_position(ann_data *d) {
	int i, j;
	FILE *fout;

	for (i=0; i<d->natom; i++) {
		d->xa[i] /= d->md_step;
		d->ya[i] /= d->md_step;
		d->za[i] /= d->md_step;

		if (d->xa[i] < 0.) {
			d->xa[i] += 1.;
		} else if (d->xa[i] >= 1.) {
			d->xa[i] -= 1.;
		}
		
		if (d->ya[i] < 0.) {
			d->ya[i] += 1.;
		} else if (d->ya[i] >= 1.) {
			d->ya[i] -= 1.;
		}
		
		if (d->za[i] < 0.) {
			d->za[i] += 1.;
		} else if (d->za[i] >= 1.) {
			d->za[i] -= 1.;
		}
	}

	if ((fout = fopen(d->average_outfile.c_str(), "w")) != NULL) {
	} else {
		printf("fail to open %s.\n", d->average_outfile.c_str());
		exit(1);
	}

	fprintf(fout, "average_position\n");
	fprintf(fout, "1.\n");
	for (i=0; i<3; i++) {
		fprintf(fout, "%25.15f%25.15f%25.15f\n", d->lv[i][0], d->lv[i][1], d->lv[i][2]);
	}
	
	for (i=0; i<d->nspecies; i++) fprintf(fout, "%10s", d->species[i].c_str());
	fprintf(fout, "\n");
	for (i=0; i<d->nspecies; i++) fprintf(fout, "%10d", d->e_natom[i]);
	fprintf(fout, "\n");

	fprintf(fout, "Direct\n");
	for (i=0; i<d->nspecies; i++) {
		for (j=0; j<d->natom; j++) {
			if (d->e[j] == d->species_id[i]) {
				fprintf(fout, "%25.15f%25.15f%25.15f\n", d->xa[j], d->ya[j], d->za[j]);
			}
		}
	}

	fclose(fout);
}

void ann_md::input_dist_pos(ann_data *d) {
	int i;
	string line;
	vector<string> s;
	ifstream fin;

	d->dist_id = new int[d->natom];
	d->x0 = new double[d->natom];
	d->y0 = new double[d->natom];
	d->z0 = new double[d->natom];
	for (i=0; i<d->natom; i++) {
		d->dist_id[i] = 0;
		d->x0[i] = 0.;
		d->y0[i] = 0.;
		d->z0[i] = 0.;
	}

	fin.open(d->dist_infile.c_str(), ios::in);
	if (!fin) {
		cout << d->dist_infile << " cannot be opened." << endl;
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

	if (d->dist_mode == "region") {
		d->dist_natom = 0;
		for (i=0; i<d->natom; i++) {
			if ((d->dist_range[0] <= d->x0[i]) && (d->x0[i] <= d->dist_range[1])) {
				if ((d->dist_range[2] <= d->y0[i]) && (d->y0[i] <= d->dist_range[3])) {
					if ((d->dist_range[4] <= d->z0[i]) && (d->z0[i] <= d->dist_range[5])) {
						d->dist_id[d->dist_natom] = i;
						d->dist_natom++;
					}
				}
			}
		}
	} else if (d->dist_mode == "atom") {
		d->dist_natom = 0;

		fin.open(d->dist_id_infile.c_str(), ios::in);
		if (!fin) {
			cout << d->dist_id_infile << " cannot be opened." << endl;
			exit(1);
		}

		getline(fin, line);
		s = split(line);
		for (i=0; i<s.size(); i++) {
			d->dist_id[i] = atoi(s[i].c_str());
		}
		d->dist_natom = s.size();

		fin.close();
	} else {
		for (i=0; i<d->natom; i++) d->dist_id[i] = i;
		d->dist_natom = d->natom;
	}

	printf("# dist_natom = %d\n", d->dist_natom);
}

void ann_md::allocation_dist(ann_data *d) {
	int i, j, k, l;

	d->dist_nstep = d->md_step / d->dist_interv;
	d->dist_ndata = int((d->dist_lim[1] - d->dist_lim[0]) / d->dist_dl) + 1;
	d->dist_c1 = 1. / (2.*sqrt(2) * pow(PI, 3./2.) * d->dist_sigma*d->dist_sigma*d->dist_sigma);
	d->dist_c2 = -1. / (2. * d->dist_sigma*d->dist_sigma);
	printf("# dist_c1 = %15.7f\n", d->dist_c1);
	printf("# dist_c2 = %15.7f\n", d->dist_c2);
	
	d->dist_xyz = new double***[d->dist_natom];
	for (i=0; i<d->dist_natom; i++) {
		d->dist_xyz[i] = new double**[d->dist_ndata];
		for (j=0; j<d->dist_ndata; j++) {
			d->dist_xyz[i][j] = new double*[d->dist_ndata];
			for (k=0; k<d->dist_ndata; k++) {
				d->dist_xyz[i][j][k] = new double[d->dist_ndata];
				for (l=0; l<d->dist_ndata; l++) {
					d->dist_xyz[i][j][k][l] = 0.;
				}
			}
		}
	}
}

void ann_md::calc_dist(ann_data *d) {
	int i, j, k, l, ii;
	double dx, dy, dz, dxx, dyy, dzz, px, py, pz, dxpx, dypy, dzpz;

	for (i=0; i<d->dist_natom; i++) {
		ii = d->dist_id[i];

		dx = d->x[ii] - d->x0[ii];
		if (dx > 0.5) {
			dx -= 1.;
		} else if (dx < -0.5) {
			dx += 1.;
		}

		dy = d->y[ii] - d->y0[ii];
		if (dy > 0.5) {
			dy -= 1.;
		} else if (dy < -0.5) {
			dy += 1.;
		}

		dz = d->z[ii] - d->z0[ii];
		if (dz > 0.5) {
			dz -= 1.;
		} else if (dz < -0.5) {
			dz += 1.;
		}

		dxx = dx*d->lv[0][0] + dy*d->lv[1][0] + dz*d->lv[2][0];
		dyy = dx*d->lv[0][1] + dy*d->lv[1][1] + dz*d->lv[2][1];
		dzz = dx*d->lv[0][2] + dy*d->lv[1][2] + dz*d->lv[2][2];

		for (j=0; j<d->dist_ndata; j++) {
			px = j * d->dist_dl + d->dist_lim[0];
			dxpx = (dxx - px) * (dxx - px);

			for (k=0; k<d->dist_ndata; k++) {
				py = k * d->dist_dl + d->dist_lim[0];
				dypy = (dyy - py) * (dyy - py);

				for (l=0; l<d->dist_ndata; l++) {
					pz = l * d->dist_dl + d->dist_lim[0];
					dzpz = (dzz - pz) * (dzz - pz);

					d->dist_xyz[i][j][k][l] += d->dist_c1 * exp(d->dist_c2 * (dxpx + dypy + dzpz));
				}
			}
		}
	}
}

void ann_md::output_dist(ann_data *d) {
	int i, j, k, l, ii;
	FILE *fout;
	char cm[256], fn[256]; 

	sprintf(cm, "mkdir %s", d->dist_outdir.c_str());
	system(cm);

	for (i=0; i<d->dist_natom; i++) {
		ii = d->dist_id[i];
		printf("ii = %d\n", ii);
		sprintf(fn, "%s/%05d.dat", d->dist_outdir.c_str(), ii);
		
		if ((fout = fopen(fn, "w")) != NULL) {
		} else {
			printf("fail to open %s.\n", fn);
			exit(1);
		}

		fprintf(fout, "# %10s%15.10f%15.10f%15.10f\n", d->species[d->e[ii]].c_str(), d->x0[ii], d->y0[ii], d->z0[ii]);
		fprintf(fout, "# dist_lim   = %15.5f%15.5f\n", d->dist_lim[0], d->dist_lim[1]);
		fprintf(fout, "# dist_dl    = %15.5f\n", d->dist_dl);
		fprintf(fout, "# dist_ndata = %10d\n", d->dist_ndata);

		for (j=0; j<d->dist_ndata; j++) {
			for (k=0; k<d->dist_ndata; k++) {
				for (l=0; l<d->dist_ndata; l++) {
					fprintf(fout, "%15.7f", d->dist_xyz[i][j][k][l] / (double)d->dist_nstep);
				}
			}
		}
		
		fclose(fout);
	}
}
