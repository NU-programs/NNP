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
#include "ann_nvt_nh.hpp"
#include <omp.h>
using namespace std;

void ann_nvt_nh::update_position(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->cx[i] += d->md_time * d->vx[i];
		d->cy[i] += d->md_time * d->vy[i];
		d->cz[i] += d->md_time * d->vz[i];
	}
}

void ann_nvt_nh::update_friction(ann_data *d) {
	int i;
	double c;
	
	d->c_ke = 0.;
	for (i=0; i<d->natom; i++) {
		d->c_ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->c_ke *= (0.5 * d->c_unit2);
	
	d->fric = d->fric + d->md_time_h * (d->c_ke - d->t_ke) / d->th_mass;
	
	c = exp(-d->fric*d->md_time_h);
	for (i=0; i<d->natom; i++) {
		d->vx[i] *= c;
		d->vy[i] *= c;
		d->vz[i] *= c;
	}
	
	d->c_ke = 0.;
	for (i=0; i<d->natom; i++) {
		d->c_ke += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->c_ke *= (0.5 * d->c_unit2);
	
	d->fric = d->fric + d->md_time_h * (d->c_ke - d->t_ke) / d->th_mass;
}

void ann_nvt_nh::update_velocity(ann_data *d) {
	double cvx, cvy, cvz;	
	int i;
	
	for (i=0; i<d->natom; i++) {
		//c = d->md_const2 / d->mass[i];
		d->vx[i] += d->mass_const2[i] * d->fcx[i];
		d->vy[i] += d->mass_const2[i] * d->fcy[i];
		d->vz[i] += d->mass_const2[i] * d->fcz[i];
	}
}

void ann_nvt_nh::copy_prev_force(ann_data *d) {
	for (int i=0; i<d->natom; i++) {
		d->pfcx[i] = d->fcx[i];
		d->pfcy[i] = d->fcy[i];
		d->pfcz[i] = d->fcz[i];
	}
}

void ann_nvt_nh::calc_current_temperature(ann_data *d) {
	d->c_temp = 0;
	for (int i=0; i<d->natom; i++) {
		d->c_temp += d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
	}
	d->c_temp = d->c_temp * 0.5 * d->c_unit2 /(1.5  * (double)d->natom * kB);
	
	if (d->md_dtemp != 0.) {
		d->md_temp += d->md_dtemp;
		d->t_ke = 0.5 * (3.*d->natom - 3.) * (kB * d->md_temp);
	}
}

void ann_nvt_nh::calc_d_tensor(ann_data *d) {
    // init D_tensor
    for (int j=0; j<6; j++) {
        for (int i=0; i<d->natom; i++) {
            d->d_tensor[i][j] = 0.0;
        }
    }
    // add E terms (ke + pe)
    double ke;
    for (int i=0; i<d->natom; i++) {
        ke = 0.5 * d->c_unit2 * d->mass[i] * (d->vx[i]*d->vx[i] + d->vy[i]*d->vy[i] + d->vz[i]*d->vz[i]);
        for (int j=0; j<3; j++) {
            d->d_tensor[i][j] += ke + d->at_en[i];
       }
    }
    // add Q term
    for (int i=0; i<d->natom; i++) {
        d->d_tensor[i][0] -= d->at_stress[i][0];
        d->d_tensor[i][1] -= d->at_stress[i][1];
        d->d_tensor[i][2] -= d->at_stress[i][2];
        d->d_tensor[i][3] -= d->at_stress[i][5];
        d->d_tensor[i][4] -= d->at_stress[i][4];
        d->d_tensor[i][5] -= d->at_stress[i][3];
    }
    // substract D_average
    int e_natom;
    for (int i=0; i<d->nspecies; i++) {
        e_natom = 0;
        for (int j=0; j<6; j++) d->d_average[j] = 0.0;
        for (int j=0; j<d->natom; j++) {
            if (d->e[j] == i) {
                for (int k=0; k<6; k++) d->d_average[k] += d->d_tensor[j][k];
                e_natom += 1;
            }
        }
        for (int j=0; j<6; j++) d->d_average[j] /= e_natom;
        for (int j=0; j<d->natom; j++) {
            if (d->e[j] == i) {
                for (int k=0; k<6; k++) d->d_tensor[j][k] -= d->d_average[k];
            }
        }
    }
    // add perturbations to atoms
    // d_tensor x: Jxy, Jxz, Jxz = 0, 3, 4
    // d_tensor y: Jyx, Jyy, Jyz = 3, 1, 5
    // d_tensor z: Jzx, Jzy, Jzz = 4, 5, 2
    // d->cjx, cjy, and cjz are set in ann_init

    for (int i=0; i<d->natom; i++) {
        d->fcx[i] += d->d_tensor[i][d->cjx] * d->fext;
        d->fcy[i] += d->d_tensor[i][d->cjy] * d->fext;
        d->fcz[i] += d->d_tensor[i][d->cjz] * d->fext;
    }
}

void ann_nvt_nh::calc_heatflux(ann_data *d) {
    for (int i=0; i<d->natom; i++){
        d->at_hflux_avg[i][0] += d->vx[i] * d->d_tensor[i][0];
        d->at_hflux_avg[i][1] += d->vy[i] * d->d_tensor[i][3];
        d->at_hflux_avg[i][2] += d->vz[i] * d->d_tensor[i][4];
        d->at_hflux_avg[i][3] += d->vx[i] * d->d_tensor[i][3];
        d->at_hflux_avg[i][4] += d->vy[i] * d->d_tensor[i][1];
        d->at_hflux_avg[i][5] += d->vz[i] * d->d_tensor[i][5];
        d->at_hflux_avg[i][6] += d->vx[i] * d->d_tensor[i][4];
        d->at_hflux_avg[i][7] += d->vy[i] * d->d_tensor[i][5];
        d->at_hflux_avg[i][8] += d->vz[i] * d->d_tensor[i][2];
    }
}

void ann_nvt_nh::calc_heatflux_mode(ann_data *d) {

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&d->vx[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d->vy[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d->vz[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(d->d_tensor[0], d->natom * 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int m, i;
    double tmp, amplitude;
    double at_hflux[3];
    // eigenvectors are already devided by sqrt(mass)

    for (m=0; m<d->n_atproc*3; m++) {
        amplitude = 0.0;
        for (i=0; i<d->natom; i++) {
            tmp =  d->local_eigenvectors[(m*d->natom+i)*3+0] * d->vx[i];
            tmp += d->local_eigenvectors[(m*d->natom+i)*3+1] * d->vy[i];
            tmp += d->local_eigenvectors[(m*d->natom+i)*3+2] * d->vz[i];
            amplitude += tmp * d->mass[i];
        }
        for (i=0; i<d->natom; i++) {
            at_hflux[0] = d->local_eigenvectors[(m*d->natom+i)*3+0] * amplitude * d->d_tensor[i][d->cjx];
            at_hflux[1] = d->local_eigenvectors[(m*d->natom+i)*3+1] * amplitude * d->d_tensor[i][d->cjy];
            at_hflux[2] = d->local_eigenvectors[(m*d->natom+i)*3+2] * amplitude * d->d_tensor[i][d->cjz];
            d->local_at_hflux_mode_avg[(m*d->natom+i)*4+0] += at_hflux[0] + at_hflux[1] + at_hflux[2];
            d->local_at_hflux_mode_avg[(m*d->natom+i)*4+1] += at_hflux[0];
            d->local_at_hflux_mode_avg[(m*d->natom+i)*4+2] += at_hflux[1];
            d->local_at_hflux_mode_avg[(m*d->natom+i)*4+3] += at_hflux[2];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void ann_nvt_nh::output_kappa(ann_data *d, int step) {

    if (step % d->n_write_pmd == 0) {
        FILE *fout;
        char filename[256];
        if (step == d->n_write_pmd) {
            if ( abs(d->fext) > 1e-10 ) {
                fout = fopen("heatflux.dat", "w");
                fprintf(fout, "# Time-averaged data for fix output %d\n", d->n_write_pmd);
                if      (d->ptb_direc == "x") fprintf(fout, "# TimeStep jx jxx jxy jxz\n");
                else if (d->ptb_direc == "y") fprintf(fout, "# TimeStep jy jyx jyy jyz\n");
                else if (d->ptb_direc == "z") fprintf(fout, "# TimeStep jz jzx jzy jzz\n");
                fclose(fout);
                fout = fopen("kappa.dat", "w");
                fprintf(fout, "# Time-averaged data for fix output %d\n", d->n_write_pmd);
                if      (d->ptb_direc == "x") fprintf(fout, "# TimeStep kx kxx kxy kxz\n");
                else if (d->ptb_direc == "y") fprintf(fout, "# TimeStep ky kyx kyy kyz\n");
                else if (d->ptb_direc == "z") fprintf(fout, "# TimeStep kz kzx kzy kzz\n");
                fclose(fout);
                if ( d->n_write_pmd > 1 ) { 
                    fout = fopen("hflux_atom.dmp", "w");
                    fclose(fout);
                }
            }
            else {
                fout = fopen("heatflux_x.dat", "w");
                fprintf(fout, "# Time-averaged data for fix output %d\n", d->n_write_pmd);
                fprintf(fout, "# TimeStep jx jxx jxy jxz\n");
                fclose(fout);
                fout = fopen("heatflux_y.dat", "w");
                fprintf(fout, "# Time-averaged data for fix output %d\n", d->n_write_pmd);
                fprintf(fout, "# TimeStep jy jyx jyy jyz\n");
                fclose(fout);
                fout = fopen("heatflux_z.dat", "w");
                fprintf(fout, "# Time-averaged data for fix output %d\n", d->n_write_pmd);
                fprintf(fout, "# TimeStep jz jzx jzy jzz\n");
                fclose(fout);
                if ( d->n_write_pmd > 1 ) { 
                    fout = fopen("hflux_atom_x.dmp", "w");
                    fclose(fout);
                    fout = fopen("hflux_atom_y.dmp", "w");
                    fclose(fout);
                    fout = fopen("hflux_atom_z.dmp", "w");
                    fclose(fout);
                }
            }
        }

        double eV2J = 1.6021766208e-19;
        double AAfs2mms = 1.0e-10 * 1.0e-10 * 1.0e-15;
        double coeff_hflux = eV2J / AAfs2mms / d->volume / (double)d->n_write_pmd;
        double coeff_kappa = 1.0 / (d->fext * 1.0e+10) / d->md_temp;

        double total_hflux[9];
        for (int i=0; i<9; i++) total_hflux[i] = 0.0;
        for (int i=0; i<d->natom; i++) {
            for (int j=0; j<9; j++) {
                d->at_hflux_avg[i][j] *= coeff_hflux;
                total_hflux[j] += d->at_hflux_avg[i][j];
            }
        }

        int cj;
        if ( abs(d->fext) > 1e-10 ) {

            if      (d->ptb_direc == "x") cj = 0;
            else if (d->ptb_direc == "y") cj = 3;
            else if (d->ptb_direc == "z") cj = 6;

            fout = fopen("heatflux.dat", "a+");
            fprintf(fout, "%08d %15.8e  %15.8e %15.8e %15.8e\n", step, \
                                                                 total_hflux[cj+0] + total_hflux[cj+1] + total_hflux[cj+2], \
                                                                 total_hflux[cj+0], total_hflux[cj+1], total_hflux[cj+2]);
            fclose(fout);
            fout = fopen("kappa.dat", "a+");
            fprintf(fout, "%08d %15.8f  %15.8f %15.8f %15.8f\n", step, \
                                                                 (total_hflux[cj+0] + total_hflux[cj+1] + total_hflux[cj+2]) * coeff_kappa, \
                                                                 total_hflux[cj+0] * coeff_kappa, \
                                                                 total_hflux[cj+1] * coeff_kappa, \
                                                                 total_hflux[cj+2] * coeff_kappa);
            fclose(fout);

            if ( d->n_write_pmd > 1) {
                double total_at_hflux_avg;
                fout = fopen("hflux_atom.dmp", "a+");
                fprintf(fout, "ITEM: TIMESTEP\n");
                fprintf(fout, "%d\n", step);
                fprintf(fout, "ITEM: NUMBER OF ATOMS\n");
                fprintf(fout, "%d\n", d->natom);
                fprintf(fout, "ITEM: BOX BOUNDS\n");
                for (int i=0; i<3; i++) fprintf(fout, "dammy\n");
                if      (d->ptb_direc == "x") fprintf(fout, "ITEM: ATOMS id type kx kxx kxy kxz\n");
                else if (d->ptb_direc == "y") fprintf(fout, "ITEM: ATOMS id type ky kyx kyy kyz\n");
                else if (d->ptb_direc == "z") fprintf(fout, "ITEM: ATOMS id type kz kzx kzy kzz\n");
                for (int i=0; i<d->natom; i++) {
                    total_at_hflux_avg = d->at_hflux_avg[i][cj+0] + d->at_hflux_avg[i][cj+1] + d->at_hflux_avg[i][cj+2];
                    fprintf(fout, "%d %d %15.8e  %15.8e %15.8e %15.8e\n", i+1, d->e[i]+1, \
                                                                               total_at_hflux_avg, \
                                                                               d->at_hflux_avg[i][cj+0], \
                                                                               d->at_hflux_avg[i][cj+1], \
                                                                               d->at_hflux_avg[i][cj+2]);
                }
                fclose(fout);
            }
        }
        else {

            cj = 0;
            fout = fopen("heatflux_x.dat", "a+");
            fprintf(fout, "%08d %15.8e  %15.8e %15.8e %15.8e\n", step, \
                                                                 total_hflux[cj+0] + total_hflux[cj+1] + total_hflux[cj+2], \
                                                                 total_hflux[cj+0], total_hflux[cj+1], total_hflux[cj+2]);
            fclose(fout);
            cj = 3;
            fout = fopen("heatflux_y.dat", "a+");
            fprintf(fout, "%08d %15.8e  %15.8e %15.8e %15.8e\n", step, \
                                                                 total_hflux[cj+0] + total_hflux[cj+1] + total_hflux[cj+2], \
                                                                 total_hflux[cj+0], total_hflux[cj+1], total_hflux[cj+2]);
            fclose(fout);
            cj = 6;
            fout = fopen("heatflux_z.dat", "a+");
            fprintf(fout, "%08d %15.8e  %15.8e %15.8e %15.8e\n", step, \
                                                                 total_hflux[cj+0] + total_hflux[cj+1] + total_hflux[cj+2], \
                                                                 total_hflux[cj+0], total_hflux[cj+1], total_hflux[cj+2]);
            fclose(fout);

            if ( d->n_write_pmd > 1) {

                double total_at_hflux_avg;
                fout = fopen("hflux_atom_x.dmp", "a+");
                fprintf(fout, "ITEM: TIMESTEP\n");
                fprintf(fout, "%d\n", step);
                fprintf(fout, "ITEM: NUMBER OF ATOMS\n");
                fprintf(fout, "%d\n", d->natom);
                fprintf(fout, "ITEM: BOX BOUNDS\n");
                for (int i=0; i<3; i++) fprintf(fout, "dammy\n");
                fprintf(fout, "ITEM: ATOMS id type kx kxx kxy kxz\n");
                for (int i=0; i<d->natom; i++) {
                    total_at_hflux_avg = d->at_hflux_avg[i][0] + d->at_hflux_avg[i][1] + d->at_hflux_avg[i][2];
                    fprintf(fout, "%d %d %15.8e  %15.8e %15.8e %15.8e\n", i+1, d->e[i]+1, total_at_hflux_avg, \
                                                                          d->at_hflux_avg[i][0], d->at_hflux_avg[i][1], d->at_hflux_avg[i][2]);
                }
                fclose(fout);

                fout = fopen("hflux_atom_y.dmp", "a+");
                fprintf(fout, "ITEM: TIMESTEP\n");
                fprintf(fout, "%d\n", step);
                fprintf(fout, "ITEM: NUMBER OF ATOMS\n");
                fprintf(fout, "%d\n", d->natom);
                fprintf(fout, "ITEM: BOX BOUNDS\n");
                for (int i=0; i<3; i++) fprintf(fout, "dammy\n");
                fprintf(fout, "ITEM: ATOMS id type ky kyx kyy kyz\n");
                for (int i=0; i<d->natom; i++) {
                    total_at_hflux_avg = d->at_hflux_avg[i][3] + d->at_hflux_avg[i][4] + d->at_hflux_avg[i][5];
                    fprintf(fout, "%d %d %15.8e  %15.8e %15.8e %15.8e\n", i+1, d->e[i]+1, total_at_hflux_avg, \
                                                                          d->at_hflux_avg[i][3], d->at_hflux_avg[i][4], d->at_hflux_avg[i][5]);
                }
                fclose(fout);

                fout = fopen("hflux_atom_z.dmp", "a+");
                fprintf(fout, "ITEM: TIMESTEP\n");
                fprintf(fout, "%d\n", step);
                fprintf(fout, "ITEM: NUMBER OF ATOMS\n");
                fprintf(fout, "%d\n", d->natom);
                fprintf(fout, "ITEM: BOX BOUNDS\n");
                for (int i=0; i<3; i++) fprintf(fout, "dammy\n");
                fprintf(fout, "ITEM: ATOMS id type kz kzx kzy kzz\n");
                for (int i=0; i<d->natom; i++) {
                    total_at_hflux_avg = d->at_hflux_avg[i][6] + d->at_hflux_avg[i][7] + d->at_hflux_avg[i][8];
                    fprintf(fout, "%d %d %15.8e  %15.8e %15.8e %15.8e\n", i+1, d->e[i]+1, total_at_hflux_avg, \
                                                                          d->at_hflux_avg[i][6], d->at_hflux_avg[i][7], d->at_hflux_avg[i][8]);
                }
                fclose(fout);
            }
        }

        for (int i=0; i<d->natom; i++) {
            for (int j=0; j<9; j++) d->at_hflux_avg[i][j] = 0.0;
        }
    }
}  

void ann_nvt_nh::output_kappa_mode(ann_data *d, int step) {

    MPI_Barrier(MPI_COMM_WORLD);

    if (step % d->n_write_pmd == 0) {

        double eV2J = 1.6021766208e-19;
        double AAfs2mms = 1.0e-10 * 1.0e-10 * 1.0e-15;
        double coeff_hflux = eV2J / AAfs2mms / d->volume / (double)d->n_write_pmd;

        for (int i=0; i<d->n_atproc * 3 * d->natom * 4; i++) {
            d->local_at_hflux_mode_avg[i] *= coeff_hflux;
        }

        MPI_Gatherv(
            d->local_at_hflux_mode_avg,
            d->n_atproc * 3 * d->natom * 4,
            MPI_DOUBLE,
            d->at_hflux_mode_avg,
            d->all_n_mode_proc,
            d->all_n_mode_start,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        if (d->rank == 0) {
            FILE *fout;
            char filename[256];
            if (step == d->n_write_pmd) {
                for (int i=0; i<d->natom; i++){
                    sprintf(filename, "hflux_mode_atom_%06d.dmp", i+1);
                    fout = fopen(filename, "w");
                    fclose(fout);
                }
            }


            for (int i=0; i<d->natom; i++){
                sprintf(filename, "hflux_mode_atom_%06d.dmp", i+1);
                fout = fopen(filename, "a+");
                fprintf(fout, "ITEM: TIMESTEP\n");
                fprintf(fout, "%d\n", step);
                fprintf(fout, "ITEM: NUMBER OF MODES\n");
                fprintf(fout, "%d\n", d->nmode);
                fprintf(fout, "ITEM: BOX BOUNDS\n");
                for (int j=0; j<3; j++) fprintf(fout, "dammy\n");
                if      (d->ptb_direc == "x") fprintf(fout, "ITEM: MODES id kx kxx kxy kxz\n");
                else if (d->ptb_direc == "y") fprintf(fout, "ITEM: MODES id ky kyx kyy kyz\n");
                else if (d->ptb_direc == "z") fprintf(fout, "ITEM: MODES id kz kzx kzy kzz\n");
                for (int m=0; m<d->nmode; m++) {
                    fprintf(fout, "%d %15.8e  %15.8e %15.8e %15.8e\n", m+1, \
                                                                       d->at_hflux_mode_avg[(m*d->natom+i)*4+0], \
                                                                       d->at_hflux_mode_avg[(m*d->natom+i)*4+1], \
                                                                       d->at_hflux_mode_avg[(m*d->natom+i)*4+2], \
                                                                       d->at_hflux_mode_avg[(m*d->natom+i)*4+3]);
                }
                fclose(fout);
            }

            for (int i=0; i< d->nmode * d->natom * 4; i++) {
                d->at_hflux_mode_avg[i] = 0.0;
            }
        }

        for (int i=0; i< d->n_atproc * 3 * d->natom * 4; i++) {
            d->local_at_hflux_mode_avg[i] = 0.0;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
