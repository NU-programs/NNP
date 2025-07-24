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
#include "ann_data.hpp"
#include "ann_mole.hpp"
#include "ann_init.hpp"
#include "ann_coordinate.hpp"
#include "ann_sf.hpp"
#include "ann_alg.hpp"
#include "ann_nvt_vsc.hpp"
#include "ann_nvt_nh.hpp"
#include "ann_nvt_lv.hpp"
#include "ann_nvt_meta.hpp"     //metadynamics
#include "ann_opt.hpp"
#include "ann_md.hpp"
#include "ann_output.hpp"
#include "ann_ti.hpp"
using namespace std;

void ann_mole::set_initial_condition(string infile) {
	init->input_parameters(infile, d);
	
	init->calc_const(d);
	if (d->rank == 0) init->check_parameters(d);
	init->input_at_pos(d);
	coord->calc_lattice_const(d);
	coord->calc_ilv(d);
	init->create_cell_list(d);

	MPI_Barrier(MPI_COMM_WORLD);

	init->calc_n_atproc(d);
	init->alloc_nn(d);
	init->alloc_descriptor(d);
	init->input_weight_bias(d);

	if (d->rank == 0) init->print_parameters(d);
}

void ann_mole::calc_dis_sf() {
	//if (d->rank == 0) cout << "01" << endl;
	MPI_Bcast(&d->x[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d->y[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d->z[0], d->natom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d->lv[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//if (d->rank == 0) cout << "02" << endl;
	sf->init_array(d);
	MPI_Barrier(MPI_COMM_WORLD);

	//if (d->rank == 0) cout << "03" << endl;
	sf->calc_distance(d);
	MPI_Barrier(MPI_COMM_WORLD);

	//if (d->rank == 0) cout << "04" << endl;
	sf->calc_g2(d);
	sf->calc_g3_term(d);
	sf->calc_g3(d);
}

void ann_mole::forward_propagation() {
	d->energy = 0.;

	alg->zscore_normalization(d);
	alg->forward_propagation(d);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(&d->t_energy, &d->energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcx[0], &d->fcx[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcy[0], &d->fcy[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcz[0], &d->fcz[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		if (d->sel_flag == "yes") {
			for (int i=0; i<d->natom; i++) {
				if (d->mflx[i] == 0) d->fcx[i] = 0.;
				if (d->mfly[i] == 0) d->fcy[i] = 0.;
				if (d->mflz[i] == 0) d->fcz[i] = 0.;
			}
		}
	}
}

void ann_mole::forward_propagation_stress() {
	d->energy = 0.;

	alg->zscore_normalization(d);
	alg->forward_propagation_stress(d);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(&d->t_energy, &d->energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);;
	MPI_Reduce(&d->t_fcx[0], &d->fcx[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcy[0], &d->fcy[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcz[0], &d->fcz[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_vstress[0], &d->vstress[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&d->t_fcx_ij[0][0], &d->fcx_ij[0][0], d->natom*d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&d->t_fcy_ij[0][0], &d->fcy_ij[0][0], d->natom*d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&d->t_fcz_ij[0][0], &d->fcz_ij[0][0], d->natom*d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		d->vstress[0] /= d->volume;
		d->vstress[1] /= d->volume;
		d->vstress[2] /= d->volume;
		d->vstress[3] /= d->volume;
		d->vstress[4] /= d->volume;
		d->vstress[5] /= d->volume;

		if (d->sel_flag == "yes") {
			for (int i=0; i<d->natom; i++) {
				if (d->mflx[i] == 0) d->fcx[i] = 0.;
				if (d->mfly[i] == 0) d->fcy[i] = 0.;
				if (d->mflz[i] == 0) d->fcz[i] = 0.;
			}
		}

		/*
		double v;

		printf("# pair force\n");
		for (int i=0; i<d->natom; i++) {
			for (int j=0; j<d->natom; j++) {
				printf("%10d%10d%15.7f\n", i, j, d->fcx_ij[i][j]);
			}
		}

		for (int i=0; i<d->natom; i++) {
			v = 0.;
			for (int j=0; j<d->natom; j++) {
				v += d->fcx_ij[i][j];
				//printf("%10d%10d%15.7f\n", i, j, d->fcx_ij[i][j]);
			}
			printf("# Total: %10d%15.7f%15.7f\n", i, d->fcx[i], v);
		}
		*/
	}
}

void ann_mole::forward_propagation_pmd() {
	double t;
	double tv[6];
	d->energy = 0.;

	alg->zscore_normalization(d);
	alg->forward_propagation_pmd(d);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(&d->t_energy, &d->energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);;
	MPI_Reduce(&d->t_fcx[0], &d->fcx[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcy[0], &d->fcy[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_fcz[0], &d->fcz[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->t_vstress[0], &d->vstress[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&d->p_at_en[0], &d->at_en[0], d->natom, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&d->p_at_stress[0][0], &d->at_stress[0][0], d->natom*6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		//t = 0.;
		//for (int i=0; i<d->natom; i++) t += d->at_en[i];
		//printf("total energy obtained from at_en = %15.10f%15.10f\n", t, t-d->energy);

		/*
		printf("# Total stress\n");
		for (int i=0; i<6; i++) {
			tv[i] = 0.;
			for (int j=0; j<d->natom; j++) {
				tv[i] += d->at_stress[j][i];
			}
			printf("%10d%15.7f%15.7f%15.7f\n", i, tv[i], d->vstress[i], d->vstress[i] / d->volume);

			d->vstress[i] = d->vstress[i] / d->volume;
			//d->vstress[i] = d->vstress[i] / d->volume / d->sf_std;
			//printf("%20.10f", d->vstress[i] * 160.2176634);
			//printf("%20.10f", d->vstress[i] * 1602.18);   // eV/ang3 to kbar
		}

		printf("# Atomic stress\n");
		for (int i=0; i<20; i++) {
			printf("%10d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n",
				i, d->at_stress[i][0], d->at_stress[i][1], d->at_stress[i][2], d->at_stress[i][3], d->at_stress[i][4], d->at_stress[i][5]);
		}
		*/
	}
}

void ann_mole::opt_cg() {
	int i, force_flag, step_length_flag;
	double prev_en, en_diff, length_lim;

	FILE *fout;

	if (d->rank == 0) {
		opt->opt_init(d);
		prev_en = 100.;

		fout = fopen(d->energy_outfile.c_str(), "w");
	}

	this->calc_dis_sf();
	this->forward_propagation();

	if (d->rank == 0) {
		opt->copy_current_position(d);
		opt->copy_current_force(d);

		prev_en = d->energy;
	}

	step_length_flag = 1;
	length_lim = 1e-6;

	for (i=0; i<d->n_iteration; ++i) {
		if (d->rank == 0) {
			if (i % d->n_write == 0) output->output_at_pos(d, i);

			if (i % d->n_write_en == 0) {
				printf("%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
				fprintf(fout, "%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
			}

			force_flag = opt->check_force_conv(d);
		}

		MPI_Bcast(&force_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (force_flag == 0) {
			break;
		} else {
			if (d->rank == 0) {
				// conjugate gradient method
				coord->frac_to_cart(d);
				opt->update_cg_direction(d);
				opt->update_position_cg(d);
				coord->cart_to_frac(d);
			}
		}

		this->calc_dis_sf();
		this->forward_propagation();

		if (d->rank == 0) {
			en_diff = prev_en - d->energy;

			if (en_diff < 0.) {
				cout << "!!! step_length is decreased !!!" << endl;
				opt->copy_previous_position(d);
				opt->copy_previous_force(d);
				d->step_length *= d->decay_rate;
			} else {
				//cout << "copy_current_position" << endl;
				opt->copy_current_position(d);
				opt->copy_current_force(d);
				prev_en = d->energy;
			}

			if (d->step_length <= length_lim) step_length_flag = 0;
		}

		MPI_Bcast(&step_length_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (step_length_flag == 0) {
			if (d->rank == 0) printf("Step length is smaller than a length min of %e.\n", length_lim);
			break;
		}
	}

	if (d->rank == 0) {
		output->output_at_pos(d, i);

		printf("\n");
		printf("##### Final potential energy #####\n");
		printf("%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);

		fprintf(fout, "\n");
		fprintf(fout, "##### Final potential energy #####\n");
		fprintf(fout, "%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);

		fclose(fout);
	}
}

void ann_mole::opt_cg_conp() {
	int i, force_flag, stress_flag, step_length_flag;     //2023.01.23 counter n added
	double prev_en, en_diff, length_lim;

	FILE *fout;

	step_length_flag = 1;
	length_lim       = 1e-6;
	force_flag       = 1;
	
	if (d->rank == 0) {
		opt->opt_init(d);
		prev_en = 100.;

		fout = fopen(d->energy_outfile.c_str(), "w");
	}

	this->calc_dis_sf();
	this->forward_propagation_stress();

	if (d->rank == 0) {
		opt->copy_current_position(d);
		opt->copy_current_force(d);

		opt->calc_dEdh(d);
		opt->copy_current_lv(d);
		opt->copy_current_dEdh(d);

		prev_en = d->energy;
	}

	for (i=0; i<d->n_iteration; ++i) {
		if (d->rank == 0) {
			if (i % d->n_write == 0) output->output_at_pos(d, i);

			if (i % d->n_write_en == 0) {
				printf("Energy: %10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
				printf("Stress: %15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", 
					d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar); // eV/ang3 to kbar
				
				fprintf(fout, "%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
			}

			force_flag = opt->check_force_conv(d);
			stress_flag = opt->check_stress_conv(d);
		}

		MPI_Bcast(&force_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&stress_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if ((force_flag == 0) && (stress_flag == 0)) {
			break;
		} else {
			if (d->rank == 0) {
				coord->frac_to_cart(d);
				opt->update_cg_direction(d);
				opt->update_position_cg(d);
				coord->cart_to_frac(d);

				opt->update_cg_direction_h(d);
				opt->update_lv_cg(d);

				coord->calc_lattice_const(d);        // 2023.03.13 d->volume update
				coord->calc_ilv(d);
			}
		}

		this->calc_dis_sf();
		this->forward_propagation_stress();
		if (d->rank == 0) opt->calc_dEdh(d);

		if (d->rank == 0) {
			en_diff = prev_en - d->energy;
			
			if (en_diff < 0.) {
				cout << "!!! step_length is decreased !!!" << endl;
				opt->copy_previous_position(d);
				opt->copy_previous_force(d);
				opt->copy_previous_lv(d);
				opt->copy_previous_dEdh(d);
				d->step_length *= d->decay_rate;
				d->step_length_lv *= d->decay_rate;
			} else {
				//cout << "copy_current_position" << endl;
				opt->copy_current_position(d);
				opt->copy_current_force(d);
				opt->copy_current_lv(d);
				opt->copy_current_dEdh(d);
				prev_en = d->energy;
			}

			if ((d->step_length <= length_lim) && (d->step_length_lv <= length_lim)) step_length_flag = 0;
		}

		MPI_Bcast(&step_length_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (step_length_flag == 0) {
			if (d->rank == 0) printf("Step length is smaller than a length min of %e.\n", length_lim);
			break;
		}
	}

	if (d->rank == 0) {
		output->output_at_pos(d, i);

		printf("\n");
		printf("##### Final potential energy #####\n");
		printf("Energy: %10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
		printf("Stress: %15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", 
			d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar); // eV/ang3 to kbar
		
		fprintf(fout, "\n");
		fprintf(fout, "##### Final potential energy #####\n");
		fprintf(fout, "%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
		fprintf(fout, "Stress: %15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", 
			d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar); // eV/ang3 to kbar

		fclose(fout);
	}
}

void ann_mole::opt_qn() {
	int i, force_flag, step_length_flag;
	double prev_en, en_diff, length_lim;

	FILE *fout;

	if (d->rank == 0) {
		opt->opt_init_qn(d);
		prev_en = 100.;

		fout = fopen(d->energy_outfile.c_str(), "w");
	}

	this->calc_dis_sf();
	this->forward_propagation();

	if (d->rank == 0) {
		opt->copy_current_position(d);
		opt->copy_current_force(d);

		prev_en = d->energy;
	}

	step_length_flag = 1;
	length_lim = 1e-4;
	d->step_length_init = d->step_length;

	for (i=0; i<d->n_iteration; ++i) {
		if (d->rank == 0) {
			if (i % d->n_write == 0) output->output_at_pos(d, i);

			printf("%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);
			fprintf(fout, "%10d%20.10f%20.10f%20e%20e\n", i, d->energy, d->energy/d->natom, en_diff, d->step_length);

			force_flag = opt->check_force_conv(d);
		}

		MPI_Bcast(&force_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (force_flag == 0) {
			break;
		} else {
			if (d->rank == 0) {
				// conjugate gradient method
				coord->frac_to_cart(d);
				opt->update_position_qn(d);
				coord->cart_to_frac(d);
			}
		}

		this->calc_dis_sf();
		this->forward_propagation();

		if (d->rank == 0) {
			en_diff = prev_en - d->energy;

			if (en_diff < 0.) {
				//cout << "copy_previous_position" << endl;
				opt->copy_previous_position(d);
				opt->copy_previous_force(d);
				d->step_length *= d->decay_rate;
			} else {
				//cout << "copy_current_position" << endl;
				opt->approximate_hessian_inverse(d);
				opt->copy_current_position(d);
				opt->copy_current_force(d);
				prev_en = d->energy;
			}

			//if (d->step_length <= length_lim) step_length_flag = 0;
			if (d->step_length <= length_lim) {
				opt->reset_hessian(d);
			}
		}

		//MPI_Bcast(&step_length_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//if (step_length_flag == 0) {
		//	if (d->rank == 0) printf("Step length is smaller than a length min of %e.\n", length_lim);
		//	break;
		//}
	}

	if (d->rank == 0) {
		output->output_at_pos(d, i);
		fclose(fout);
	}
}

void ann_mole::run_md() {
	int i;

	if (d->rank == 0) {
		md->set_md_const(d);
		md->md_init(d);
		//md->zero_velocity_center_of_mass(d);

		if (d->simulation_type == "ti") {
			ti->alloc(d);
			ti->input_ti_pos(d);
			ti->input_fc(d);
		}

		if (d->velauto_flag == "yes") {
			md->allocation_velauto(d);
		}

		if (d->msd_flag == "yes") {
			md->allocation_msd(d);
		}

		if (d->average_flag == "yes") {
			md->allocation_average_position(d);
		}

		if (d->dist_flag == "yes") {
			md->input_dist_pos(d);
			md->allocation_dist(d);
		}
	}

	this->calc_dis_sf();
	this->forward_propagation_stress();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		if (d->simulation_type == "ti") ti->init(d);
		output->output_at_pos(d, 0);
	}

	// velocity scaling
	if (d->vsc_step != 0) {
		if (d->rank == 0) printf("Velocity scaling starts...\n");

		for (i=1; i<=d->vsc_step; i++) {
			this->md_nvt_vsc();
			this->output_md_result(i);
		}
	}

	// Nose-Hoover or Langevin
	if (d->md_condition == "nvt_nh") {
		if (d->rank == 0) printf("Nose-hoover nvt starts...\n");

		for (i=1; i<=d->md_step; i++) {
			this->md_nvt_nh(i);
			this->output_md_result(i);
		}

	} else if (d->md_condition == "npt_nh") {
		if (d->rank == 0) {
			printf("Nose-hoover npt starts...\n");
			npt_nh->calc_npt_const(d);
		}

		for (i=1; i<=d->md_step; i++) {
			this->md_npt_nh(i);
			this->output_md_result(i);
		}

	} else if (d->md_condition == "nvt_lv") {
		if (d->rank == 0) nvt_lv->initialize_lv(d);

		for (i=1; i<=d->md_step; i++) {
			this->md_nvt_lv();
			this->output_md_result(i);
		}

	} else if (d->md_condition == "nve") {
		for (i=1; i<=d->md_step; i++) {
			this->md_nve(i);
			this->output_md_result(i);
		}

	} else if (d->md_condition == "nvt_pmd") {
		if (d->rank == 0) printf("Nose-hoover with pmd starts...\n");
		
    	md->allocation_pmd(d);

		for (i=1; i<=d->md_step; i++) {
			this->md_nvt_pmd(i);
			this->output_md_result(i);
		}

	} else {
		cout << "Undefine parameter is set to md_condition." << endl;
		MPI_Abort(MPI_COMM_WORLD, 32);
	}

	if (d->rank == 0) {
		if ((i-1) % d->n_write != 0) {
			output->output_at_pos(d, i-1);
		}

		if (d->velauto_flag == "yes") {
			md->average_velauto(d);
			md->output_velauto(d);
		}

		if (d->msd_flag == "yes") {
			md->average_msd(d);
			md->output_msd(d);
		}

		if (d->average_flag == "yes") {
			md->output_average_position(d);
		}

		if (d->dist_flag == "yes") {
			md->output_dist(d);
		}
	}
}

void ann_mole::output_md_result(int i) {
	if (d->rank == 0) {
		FILE *fout;

		fout = fopen(d->energy_outfile.c_str(), "a+");

		if (d->simulation_type == "md") {
			if (i % d->n_write_en == 0) {
				if (d->md_condition == "npt_nh") {
					printf("%10.2f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.7f\n",
						i*d->md_time, d->c_temp, d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar,
						d->lt[0], d->lt[1], d->lt[2], d->lt[3], d->lt[4], d->lt[5], d->volume, d->energy);

					fprintf(fout, "%10.2f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.7f\n",
						i*d->md_time, d->c_temp, d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar,
						d->lt[0], d->lt[1], d->lt[2], d->lt[3], d->lt[4], d->lt[5], d->volume, d->energy);
                } else if (d->nvt_stress == "yes") {
                    printf("%10d%15.2f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.7f\n", i, i*d->md_time, d->md_temp, d->c_temp, d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar, d->energy);
                    fprintf(fout, "%10d%15.2f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.7f\n", i, i*d->md_time, d->md_temp, d->c_temp, d->vstress[0]*evang3_to_kbar, d->vstress[1]*evang3_to_kbar, d->vstress[2]*evang3_to_kbar, d->vstress[3]*evang3_to_kbar, d->vstress[4]*evang3_to_kbar, d->vstress[5]*evang3_to_kbar, d->energy);
				} else {
					printf("%10d%15.2f%20.7f%20.7f%20.7f\n", i, i*d->md_time, d->md_temp, d->c_temp, d->energy);
					fprintf(fout, "%10d%15.2f%20.7f%20.7f%20.7f\n", i, i*d->md_time, d->md_temp, d->c_temp, d->energy);
				}
			}
		} else if (d->simulation_type == "ti") {
			if (i % d->n_write_en == 0) {
				printf("%10d%15.2f%20.7f%20.7f%20.7f%20.7f%20.7f%20.7f\n", i, i*d->md_time, d->md_temp, d->c_temp, d->ti_energy, d->energy, d->ha_energy, d->ti_diff);
				fprintf(fout, "%10d%15.2f%20.7f%20.7f%20.7f%20.7f%20.7f%20.7f\n", i, i*d->md_time, d->md_temp, d->c_temp, d->ti_energy, d->energy, d->ha_energy, d->ti_diff);
			}
		}

		fclose(fout);

		if (i % d->n_write == 0) {
			output->output_at_pos(d, i);
			//if (d->simulation_type == "ti") {
			//	output->output_force(d, i);
			//}
		}
	}
}

void ann_mole::md_nvt_vsc() {
	if (d->rank == 0) {
		vsc->velocity_scaling(d);
		coord->frac_to_cart(d);
		vsc->update_position(d);
		coord->cart_to_frac(d);
		vsc->copy_prev_force(d);
	}

	this->calc_dis_sf();
	//this->forward_propagation();
	this->forward_propagation_stress();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->simulation_type == "ti") this->calc_ti();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		vsc->update_velocity(d);
		md->zero_velocity_center_of_mass(d);
	}
}

void ann_mole::md_nvt_nh(int step) {
	if (d->rank == 0) {
		nvt_nh->update_friction(d);
		nvt_nh->update_velocity(d);
		md->zero_velocity_center_of_mass(d);

		coord->frac_to_cart(d);
		nvt_nh->update_position(d);
		coord->cart_to_frac(d);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	this->calc_dis_sf();
	if (d->nvt_stress == "yes") {
		this->forward_propagation_stress();
	} else {
		this->forward_propagation();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->simulation_type == "ti") this->calc_ti();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		nvt_nh->update_velocity(d);
		md->zero_velocity_center_of_mass(d);
		nvt_nh->update_friction(d);
		nvt_nh->calc_current_temperature(d);

		if (d->velauto_flag == "yes") md->calc_velauto(d, step);
		if (d->msd_flag == "yes") md->calc_msd(d, step);
		if (d->average_flag == "yes") md->calc_average_position(d);
		
		if (d->dist_flag == "yes") {
			if (step % d->dist_interv == 0) md->calc_dist(d);
		}
	}
}

void ann_mole::md_npt_nh(int step) {
	if (d->rank == 0) {
		// thermostat
		//cout << "thermostat01" << endl;
		npt_nh->calc_ba_fric_tr(d);
		npt_nh->update_th_friction(d);
		npt_nh->update_th_velocity(d);
		//md->zero_velocity_center_of_mass(d);
		npt_nh->update_th_friction(d);

		// barostat
		//cout << "barostat01" << endl;
		npt_nh->calc_total_stress(d);
		npt_nh->update_ba_friction(d);
		npt_nh->update_ba_velocity(d);
		//md->zero_velocity_center_of_mass(d);
		npt_nh->calc_total_stress(d);
		npt_nh->update_ba_friction(d);

		// thermostat
		//cout << "thermostat02" << endl;
		npt_nh->calc_ba_fric_tr(d);
		npt_nh->update_th_friction(d);
		npt_nh->update_th_velocity(d);
		//md->zero_velocity_center_of_mass(d);
		npt_nh->update_th_friction(d);

		// update velocity, lattice vector and position
		//cout << "update01" << endl;
		npt_nh->update_velocity(d);
		//md->zero_velocity_center_of_mass(d);

		coord->frac_to_cart(d);
		npt_nh->calc_exp2(d);
		npt_nh->update_lv(d);
		coord->calc_lattice_const(d);
		coord->calc_ilv(d);

		npt_nh->update_position(d);
		coord->cart_to_frac(d);

		//output->output_at_pos(d, step);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	this->calc_dis_sf();
	this->forward_propagation_stress();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		// update velocity
		npt_nh->update_velocity(d);
		//md->zero_velocity_center_of_mass(d);

		// themostat
		//npt_nh->calc_ba_fric_tr(d);
		npt_nh->update_th_friction(d);
		npt_nh->update_th_velocity(d);
		//md->zero_velocity_center_of_mass(d);
		npt_nh->update_th_friction(d);

		// barostat
		npt_nh->calc_total_stress(d);
		npt_nh->update_ba_friction(d);
		npt_nh->update_ba_velocity(d);
		npt_nh->calc_total_stress(d);
		//md->zero_velocity_center_of_mass(d);
		npt_nh->update_ba_friction(d);

		// thermostat
		npt_nh->calc_ba_fric_tr(d);
		npt_nh->update_th_friction(d);
		npt_nh->update_th_velocity(d);
		//md->zero_velocity_center_of_mass(d);
		npt_nh->update_th_friction(d);

		md->zero_velocity_center_of_mass(d);

		nvt_nh->calc_current_temperature(d);
	}
}

void ann_mole::md_nve(int step) {
	if (d->rank == 0) {
		//nvt_nh->update_friction(d);
		nvt_nh->update_velocity(d);
		md->zero_velocity_center_of_mass(d);

		coord->frac_to_cart(d);
		nvt_nh->update_position(d);
		coord->cart_to_frac(d);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	this->calc_dis_sf();
	if (d->nvt_stress == "yes") {
		this->forward_propagation_stress();
	} else {
		this->forward_propagation();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->simulation_type == "ti") this->calc_ti();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		nvt_nh->update_velocity(d);
		md->zero_velocity_center_of_mass(d);
		//nvt_nh->update_friction(d);
		nvt_nh->calc_current_temperature(d);

		if (d->velauto_flag == "yes") md->calc_velauto(d, step);
	}
}

void ann_mole::md_nvt_pmd(int step) {
	if (d->rank == 0) {
		nvt_nh->update_friction(d);
		nvt_nh->update_velocity(d);
		//md->zero_velocity_center_of_mass(d);

		coord->frac_to_cart(d);
		nvt_nh->update_position(d);
		coord->cart_to_frac(d);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	this->calc_dis_sf();
	this->forward_propagation_pmd();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
        nvt_nh->calc_d_tensor(d);

		nvt_nh->update_velocity(d);
		//md->zero_velocity_center_of_mass(d);
		nvt_nh->update_friction(d);
		nvt_nh->calc_current_temperature(d);

        nvt_nh->calc_heatflux(d);
        nvt_nh->output_kappa(d, step);
    }

	MPI_Barrier(MPI_COMM_WORLD);
    if (d->mode_decomp == "yes") {
        nvt_nh->calc_heatflux_mode(d);
        nvt_nh->output_kappa_mode(d, step);
    }
}

void ann_mole::md_nvt_lv() {
	if (d->rank == 0) {
		nvt_lv->generate_gausian_variable(d);
		nvt_lv->update_velocity(d);
		md->zero_velocity_center_of_mass(d);

		coord->frac_to_cart(d);
		nvt_lv->update_position(d);
		coord->cart_to_frac(d);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	this->calc_dis_sf();
	if (d->nvt_stress == "yes") {
		this->forward_propagation_stress();
	} else {
		this->forward_propagation();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->simulation_type == "ti") this->calc_ti();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		nvt_lv->update_velocity(d);
		nvt_lv->calc_current_temperature(d);
	}
}

void ann_mole::calc_ti() {
	if (d->rank == 0) {
		ti->disp(d);
		ti->calc_fc(d);
	}
}

/*
void ann_mole::run_meta() {
	int i;
	FILE *fout;

	if (d->rank == 0) {
		md->set_md_const(d);
		md->md_init(d);
		fout = fopen(d->energy_outfile.c_str(), "w");
	}

	this->calc_dis_sf();
	this->forward_propagation();
	MPI_Barrier(MPI_COMM_WORLD);

	if(d->rank == 0) {
		nvt_meta->initialize_lv(d);
		coord->frac_to_cart(d);
		nvt_meta->copy_prev_position(d);
		nvt_meta->cutoff_particle_list(d);
		nvt_meta->set_parameter(d);
		coord->cart_to_frac(d);
		nvt_meta->make_CV(d);
		nvt_meta->partial_diff(d);
		nvt_meta->save_s_of_t(d);
		nvt_meta->calc_vias_potential(d);
		nvt_meta->force_plus_vias(d);

		printf("%10d%15.2f%15.7f%15.7f%15.7f\n", 0, 0., d->md_temp, d->c_temp, d->energy);
		//fprintf(fout, "%10d%15.2f%15.7f%15.7f%15.7f\n", 0, 0., d->md_temp, d->c_temp, d->energy);

		output->output_at_pos(d, 0);
	}

	for (i=1; i<=d->md_step; i++) {
		d->md_time_number = i*d->md_time;
		this->md_nvt_meta();
		if (d->rank == 0) {
			printf("%10d%15.2f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n",
			       i, i*d->md_time, d->md_temp, d->c_temp, d->energy, d->xg, d->yg, d->zg, d->plus_Vg, d->wh[d->NG]);
			fprintf(fout, "%10d%15.2f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n",
			        i, i*d->md_time, d->md_temp, d->c_temp, d->energy, d->xg, d->yg, d->zg, d->plus_Vg, d->wh[d->NG]);

			if (i % d->n_write == 0) output->output_at_pos(d, i);
			// if (i % d->n_write == 0) output->output_at_xyz(d, i);
		}
	}

	if (d->rank == 0) {
		output->output_at_pos(d, i);
		fclose(fout);
	}
}

void ann_mole::md_nvt_meta() {
	if (d->rank == 0){
		nvt_meta->generate_gaussian_variable(d);
		nvt_meta->update_velocity(d);
		coord->frac_to_cart(d);
		nvt_meta->update_position(d);
		coord->cart_to_frac(d);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	this->calc_dis_sf();
	this->forward_propagation();
	MPI_Barrier(MPI_COMM_WORLD);

	if (d->rank == 0) {
		coord->frac_to_cart(d);
		nvt_meta->make_CV(d);
		coord->cart_to_frac(d);
		nvt_meta->partial_diff(d);
		if ((d->md_time_number % d->tau_G) == 0) {
			nvt_meta->save_s_of_t(d);
		}
		nvt_meta->calc_vias_potential(d);
		nvt_meta->force_plus_vias(d);

		nvt_meta->update_velocity(d);
		nvt_meta->calc_current_temperature(d);
	}
}
*/
