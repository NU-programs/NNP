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
#include "ann_init.hpp"
using namespace std;

vector<string> ann_init::split (string line) {
	istringstream iss(line);
	string s;
	vector<string> ss;
	while (iss >> s) ss.push_back(s);
	return ss;
}

void ann_init::input_parameters(string fn, ann_data *d) {
	int i, j;
	string line;
	vector<string> s;

	d->sel_flag      = "no";
	d->fc_diff       = 1e-3;
	d->decay_rate    = 0.9;

	d->max_natom     = 80;
	d->max_lc_natom  = 150;

	d->n_write       = 100;
	d->n_write_en    = 1;

	d->simulation_type = "";
	d->at_infile       = "";
	d->wb_infile       = "";
	
	d->fc_infile       = "";
	d->ti_infile       = "";

	d->cutoff    = 0.;
	d->cutoff2   = 0.;
	d->ng2       = 0;
	d->ng3       = 0;
	d->n_hidden1 = 0;
	d->n_hidden2 = 0;

	// md
	d->nvt_stress = "no";
    d->md_pressure_x = 0.;
    d->md_pressure_y = 0.;
    d->md_pressure_z = 0.;
	d->npt_shear = "yes";

	// velocity autocorrelation
	d->velauto_flag    = "no";
	d->velauto_at_flag = "no";
	d->velauto_outdir  = "velauto";
	d->velauto_n       = 10;
	d->velauto_step    = 1000;
	d->velauto_inv     = 100;

	// msd
	d->msd_flag    = "no";
	d->msd_outdir  = "msd";
	d->msd_n       = 10;
	d->msd_step    = 1000;
	d->msd_inv     = 100;

	// average position
	d->average_flag    = "no";
	d->average_outfile = "average.vasp";

	// distribution
	d->dist_flag      = "no";
	d->dist_outdir    = "dist";
	d->dist_mode      = "all";
	d->dist_infile    = "POSCAR.vasp";
	d->dist_id_infile = "dist_id.dat";
	d->dist_lim[0]    = -1.;
	d->dist_lim[1]    =  1.;
	d->dist_dl        = 0.02;
	d->dist_sigma     = 0.02;
	d->dist_interv    = 1;
	d->dist_range[0]  = 0.;
	d->dist_range[1]  = 1.;
	d->dist_range[2]  = 0.;
	d->dist_range[3]  = 1.;
	d->dist_range[4]  = 0.;
	d->dist_range[5]  = 1.;

	// stress
	d->lv_update = "all";

    // pmd
    d->mode_decomp = "no";
    d->mode_infile = "";

	ifstream fin;
	fin.open(fn.c_str(), ios::in);
	if(!fin) {
		printf("%s cannot be opened.\n", fn.c_str());
		exit(1);
	}

	while(getline(fin, line)) {
		s = split(line);

		if (s.size() >= 2) {
			if (s[0] == "simulation_type")    { d->simulation_type = s[1]; }

			else if (s[0] == "at_infile")     { d->at_infile = s[1]; }
			else if (s[0] == "wb_infile")     { d->wb_infile = s[1]; }

			else if (s[0] == "species")     {
				d->nspecies = s.size()-1;
				d->species = new string[s.size()-1];
				d->species_id = new int[s.size()-1];
				d->e_natom = new int[s.size()-1];

				for (i=0; i<s.size()-1; i++) {
					d->species[i] = s[i+1];
					d->species_id[i] = i;
				}
			}

			else if (s[0] == "ng2")              { d->ng2            = atoi(s[1].c_str()); }
			else if (s[0] == "ng3")              { d->ng3            = atoi(s[1].c_str()); }
			else if (s[0] == "n_hidden1")        { d->n_hidden1      = atoi(s[1].c_str()); }
			else if (s[0] == "n_hidden2")        { d->n_hidden2      = atoi(s[1].c_str()); }
			else if (s[0] == "cutoff")           { d->cutoff         = atof(s[1].c_str()); }
			else if (s[0] == "max_natom")        { d->max_natom      = atoi(s[1].c_str()); }
			else if (s[0] == "max_lc_natom")     { d->max_lc_natom   = atoi(s[1].c_str());}
			else if (s[0] == "n_write")          { d->n_write        = atoi(s[1].c_str()); }
			else if (s[0] == "n_write_en")       { d->n_write_en     = atoi(s[1].c_str()); }
			else if (s[0] == "sel_flag")         { d->sel_flag       = s[1]; }
			else if (s[0] == "opt_alg")          { d->opt_alg        = s[1]; }
			else if (s[0] == "step_length")      { d->step_length    = atof(s[1].c_str()); }
			else if (s[0] == "step_length_lv")   { d->step_length_lv = atof(s[1].c_str()); }
			else if (s[0] == "lv_update")        { d->lv_update      = s[1]; }
			else if (s[0] == "H_init")           { d->H_init         = atof(s[1].c_str()); }
			else if (s[0] == "decay_rate")       { d->decay_rate     = atof(s[1].c_str()); }
			else if (s[0] == "en_diff")          { d->en_diff        = atof(s[1].c_str()); }
			else if (s[0] == "fc_diff")          { d->fc_diff        = atof(s[1].c_str()); }
			else if (s[0] == "stress_diff")      { d->stress_diff    = atof(s[1].c_str()); }
			else if (s[0] == "n_iteration")      { d->n_iteration    = atoi(s[1].c_str()); }
			else if (s[0] == "energy_outfile")   { d->energy_outfile = s[1]; }
			else if (s[0] == "vsc_step")         { d->vsc_step       = atoi(s[1].c_str()); }
			else if (s[0] == "md_step")          { d->md_step        = atoi(s[1].c_str()); }
			else if (s[0] == "init_vel")         { d->init_vel       = s[1]; }
			else if (s[0] == "md_condition")     { d->md_condition   = s[1]; }
			else if (s[0] == "nvt_stress")       { d->nvt_stress     = s[1]; }
			else if (s[0] == "md_time")          { d->md_time        = atof(s[1].c_str()); }
			else if (s[0] == "md_init_temp")     { d->md_init_temp   = atof(s[1].c_str()); }
			else if (s[0] == "md_fin_temp")      { d->md_fin_temp    = atof(s[1].c_str()); }
			else if (s[0] == "th_mass")          { d->th_mass        = atof(s[1].c_str()); }
			else if (s[0] == "md_sigma")         { d->md_sigma       = atof(s[1].c_str()); }
			else if (s[0] == "md_friction")      { d->md_friction    = atof(s[1].c_str()); }
			else if (s[0] == "tau_G")            { d->tau_G          = atoi(s[1].c_str()); }
			else if (s[0] == "gau_h")            { d->gau_h          = atof(s[1].c_str()); }
			else if (s[0] == "gau_w")            { d->gau_w          = atof(s[1].c_str()); }
			else if (s[0] == "meta_cutoff")      { d->meta_cutoff    = atof(s[1].c_str()); }
			else if (s[0] == "delT")             { d->delT           = atof(s[1].c_str()); }
			else if (s[0] == "ha_energy0")       { d->ha_energy0     = atof(s[1].c_str()); }
			else if (s[0] == "ti_lambda")        { d->ti_lambda      = atof(s[1].c_str()); d->ti_lambda_h = 1. - d->ti_lambda; }
			else if (s[0] == "fc_infile")        { d->fc_infile      = s[1].c_str(); }
			else if (s[0] == "ti_infile")        { d->ti_infile      = s[1].c_str(); }

			// npt
			else if (s[0] == "ba_mass")          { d->ba_mass        = atof(s[1].c_str()); }
			else if (s[0] == "md_pressure")      { d->md_pressure_x  = atof(s[1].c_str());
			                                       d->md_pressure_y  = atof(s[2].c_str());
			                                       d->md_pressure_z  = atof(s[3].c_str()); }
			else if (s[0] == "npt_shear")        { d->npt_shear      = s[1].c_str(); }

			// velocity autocorrelation
			else if (s[0] == "velauto_flag")     { d->velauto_flag    = s[1].c_str(); }
			else if (s[0] == "velauto_at_flag")  { d->velauto_at_flag = s[1].c_str(); }
			else if (s[0] == "velauto_outdir")   { d->velauto_outdir  = s[1].c_str(); }
			else if (s[0] == "velauto_n")        { d->velauto_n       = atoi(s[1].c_str()); }
			else if (s[0] == "velauto_step")     { d->velauto_step    = atoi(s[1].c_str()); }
			else if (s[0] == "velauto_inv")      { d->velauto_inv     = atoi(s[1].c_str()); }

			// msd
			else if (s[0] == "msd_flag")         { d->msd_flag    = s[1].c_str(); }
			else if (s[0] == "msd_outdir")       { d->msd_outdir  = s[1].c_str(); }
			else if (s[0] == "msd_n")            { d->msd_n       = atoi(s[1].c_str()); }
			else if (s[0] == "msd_step")         { d->msd_step    = atoi(s[1].c_str()); }
			else if (s[0] == "msd_inv")          { d->msd_inv     = atoi(s[1].c_str()); }

			// average position
			else if (s[0] == "average_flag")     { d->average_flag    = s[1].c_str(); }
			else if (s[0] == "average_outfile")  { d->average_outfile = s[1].c_str(); }

			// distribution
			else if (s[0] == "dist_flag")        { d->dist_flag      = s[1].c_str(); }
			else if (s[0] == "dist_outdir")      { d->dist_outdir    = s[1].c_str(); }
			else if (s[0] == "dist_infile")      { d->dist_infile    = s[1].c_str(); }
			else if (s[0] == "dist_mode")        { d->dist_mode      = s[1].c_str(); }
			else if (s[0] == "dist_id_infile")   { d->dist_id_infile = s[1].c_str(); }
			else if (s[0] == "dist_dl")          { d->dist_dl        = atof(s[1].c_str()); }
		    else if (s[0] == "dist_sigma")       { d->dist_sigma     = atof(s[1].c_str()); }
		    else if (s[0] == "dist_interv")      { d->dist_interv    = atoi(s[1].c_str()); }
			else if (s[0] == "dist_lim")         { d->dist_lim[0]    = atof(s[1].c_str());
		                                           d->dist_lim[1]    = atof(s[2].c_str()); }
		    else if (s[0] == "dist_range")       { d->dist_range[0]  = atof(s[1].c_str());
		                                           d->dist_range[1]  = atof(s[2].c_str());
		                                           d->dist_range[2]  = atof(s[3].c_str());
		                                           d->dist_range[3]  = atof(s[4].c_str());
		                                           d->dist_range[4]  = atof(s[5].c_str());
		                                           d->dist_range[5]  = atof(s[6].c_str()); }

			else if (s[0] == "md_mass") {
				d->md_mass = new double[s.size()-1];
				for (i=0; i<s.size()-1; i++) d->md_mass[i] = atof(s[i+1].c_str());
			}

			// pmd
            else if (s[0] == "fext") {
                d->fext =  atof(s[1].c_str());
                if (abs(d->fext) < 1e-10) {
                    if (d->rank == 0) cout << "fext = 0 defected. In PMD, heat flux is output in all directions, and kappa is not calculated" << endl;
                }
            }
			else if (s[0] == "ptb_direc"){
                d->ptb_direc = s[1].c_str();
                if (d->ptb_direc == "x") {
                    d->cjx = 0; d->cjy = 3; d->cjz = 4;
                }
                else if (d->ptb_direc == "y") {
                    d->cjx = 3; d->cjy = 1; d->cjz = 5;
                }
                else if (d->ptb_direc == "z") {
                    d->cjx = 4; d->cjy = 5; d->cjz = 2;
                }
            }
			else if (s[0] == "n_write_pmd")      {
                d->n_write_pmd = atoi(s[1].c_str());
                if (d->n_write_pmd == 1) {
                    if (d->rank == 0) cout << "n_write_pmd = 1 defected. In PMD, katom is not output as hflux_atom.dmp will be huge." << endl;
                }
            }
            else if (s[0] == "mode_decomp")      { d->mode_decomp = s[1].c_str(); }
            else if (s[0] == "mode_infile")      { d->mode_infile = s[1].c_str(); }
		}
	}

	fin.close();

	if ((d->velauto_flag == "yes") && ((d->msd_flag == "no"))) {
		d->md_step = d->velauto_step + d->velauto_n * d->velauto_inv;

	} else if ((d->velauto_flag == "no") && ((d->msd_flag == "yes"))) {
		d->md_step = d->msd_step + d->msd_n * d->msd_inv;

	} else if ((d->velauto_flag == "yes") && ((d->msd_flag == "yes"))) {
		i = d->velauto_step + d->velauto_n * d->velauto_inv;
		j = d->msd_step + d->msd_n * d->msd_inv;
		i >= j ? (d->md_step = i) : (d->md_step = j);
	}
}

void ann_init::calc_const(ann_data *d) {
	int i, j;
	int c;

	d->cutoff2 = d->cutoff * d->cutoff;

	d->tanh_const0 = 2.;
	d->tanh_const1 = -1. / d->cutoff;
	d->tanh_const2 = -3. / d->cutoff * d->tanh_const0;
	d->ch_const1 = 2. / d->cutoff;

	d->t_ng2 = d->ng2 * d->nspecies;
	d->t_ng3 = d->ng3 * d->nspecies * (d->nspecies + 1) / 2;
	d->n_symmetry = d->t_ng2 + d->t_ng3;

	d->g3x = new int*[d->nspecies];
	for (i=0; i<d->nspecies; i++) {
		d->g3x[i] = new int[d->nspecies];

		for (j=0; j<d->nspecies; j++) {
			d->g3x[i][j] = 0;
		}
	}

	c = 0;
	for (i=0; i<d->nspecies; i++) {
		for (j=i; j<d->nspecies; j++) {
			d->g3x[i][j] = c;
			d->g3x[j][i] = c;
			c++;
		}
	}
}

void ann_init::check_parameters(ann_data *d) {
	if (d->simulation_type == "") {
		printf("!!!!! simulation type is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 51);
	}

	if ((d->simulation_type != "opt") &&
		(d->simulation_type != "md") &&
		(d->simulation_type != "ti") &&
		(d->simulation_type != "meta")) {
		printf("!!!!! Undefined simulation type is set : %s !!!!!\n", d->simulation_type.c_str());
		MPI_Abort(MPI_COMM_WORLD, 56);
	}

	if (d->at_infile == "") {
		printf("!!!!! at_infile is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 52);
	}

	if (d->wb_infile == "") {
		printf("!!!!! wb_infile is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 53);
	}

	if (d->cutoff == 0) {
		printf("!!!!! cutoff is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 54);
	}

	if (d->ng2 == 0) {
		printf("!!!!! ng2 is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 54);
	}

	if (d->ng3 == 0) {
		printf("!!!!! ng3 is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 54);
	}

	if (d->n_hidden1 == 0) {
		printf("!!!!! n_hidden1 is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 54);
	}

	if (d->n_hidden2 == 0) {
		printf("!!!!! n_hidden2 is not set !!!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 55);
	}
}

void ann_init::print_parameters(ann_data *d) {
	printf("\n############# Default parameters ######################\n");

	printf("# simulation_type  = %s\n", d->simulation_type.c_str());
	printf("# at_infile        = %s\n", d->at_infile.c_str());
	printf("# wb_infile        = %s\n", d->wb_infile.c_str());

	printf("# nspecies         = %d\n", d->nspecies);
	printf("# species          =");
	for (int i=0; i<d->nspecies; i++) printf(" %s", d->species[i].c_str());
	printf("\n");

	printf("# ng2              = %d\n", d->ng2);
	printf("# t_ng2            = %d\n", d->t_ng2);
	printf("# ng3              = %d\n", d->ng3);
	printf("# t_ng3            = %d\n", d->t_ng3);
	printf("# n_symmetry       = %d\n", d->n_symmetry);
	printf("# n_hidden1        = %d\n", d->n_hidden1);
	printf("# n_hidden2        = %d\n", d->n_hidden2);
	printf("# sf_mean          = %f\n", d->sf_mean);
	printf("# sf_std           = %f\n", d->sf_std);
	printf("# cutoff           = %f\n", d->cutoff);

	printf("# tanh_const1      = %f\n", d->tanh_const1);
	printf("# tanh_const2      = %f\n", d->tanh_const2);
	printf("# ch_const1        = %f\n", d->ch_const1);

	printf("# max_natom        = %d\n", d->max_natom);
	printf("# max_lc_natom     = %d\n", d->max_lc_natom);
	printf("# natom            = %d\n", d->natom);
	printf("# natom3           = %d\n", d->natom3);
	printf("# sel_flag         = %s\n", d->sel_flag.c_str());
	printf("# energy_outfile   = %s\n", d->energy_outfile.c_str());
	printf("# n_write          = %d\n", d->n_write);

	if (d->simulation_type == "opt") {
		printf("# step_length     = %e\n", d->step_length);
		printf("# step_length_lv  = %e\n", d->step_length_lv);
		printf("# decay_rate      = %e\n", d->decay_rate);
		printf("# en_diff         = %e\n", d->en_diff);
		printf("# fc_diff         = %e\n", d->fc_diff);
		printf("# stress_diff     = %e\n", d->stress_diff);
		printf("# n_iteration     = %d\n", d->n_iteration);
	}
	if ((d->simulation_type == "md") || (d->simulation_type == "ti")) {
		printf("# vsc_step         = %d\n", d->vsc_step);
		printf("# md_condition     = %s\n", d->md_condition.c_str());
		printf("# md_step          = %d\n", d->md_step);
		printf("# md_time          = %f\n", d->md_time);
		printf("# md_init_temp     = %f\n", d->md_init_temp);
		printf("# md_fin_temp      = %f\n", d->md_fin_temp);
		printf("# md_pressure      = %f %f %f\n", d->md_pressure_x, d->md_pressure_y, d->md_pressure_z);
		printf("# md_mass          =");
		for (int i=0; i<d->nspecies; i++) {
			printf(" %.7f", d->md_mass[i]);
		}
		printf("\n");

		if (d->md_condition == "nvt_nh") {
			printf("# th_mass          = %f\n", d->th_mass);
		} else if (d->md_condition == "npt_nh") {
			printf("# ba_mass          = %f\n", d->ba_mass);
		} else if (d->md_condition == "nvt_lv") {
			printf("# md_sigma         = %f\n", d->md_sigma);
			printf("# md_friction      = %f\n", d->md_friction);
		}
	}
	if (d->simulation_type == "ti") {
		printf("# ha_energy0       = %-15.7f\n", d->ha_energy0);
		printf("# ti_lambda        = %e\n", d->ti_lambda);
		printf("# ti_lambda_h      = %e\n", d->ti_lambda_h);
		printf("# fc_infile        = %s\n", d->fc_infile.c_str());
		printf("# ti_infile        = %s\n", d->ti_infile.c_str());
	}

	if (d->velauto_flag == "yes") {
		printf("# velauto_flag     = %s\n", d->velauto_flag.c_str());
		printf("# velauto_at_flag  = %s\n", d->velauto_at_flag.c_str());
		printf("# velauto_outdir   = %s\n", d->velauto_outdir.c_str());
		printf("# velauto_n        = %d\n", d->velauto_n);
		printf("# velauto_step     = %d\n", d->velauto_step);
		printf("# velauto_inv      = %d\n", d->velauto_inv);
	}

	if (d->msd_flag == "yes") {
		printf("# msd_flag         = %s\n", d->msd_flag.c_str());
		printf("# msd_outdir       = %s\n", d->msd_outdir.c_str());
		printf("# msd_n            = %d\n", d->msd_n);
		printf("# msd_step         = %d\n", d->msd_step);
		printf("# msd_inv          = %d\n", d->msd_inv);
	}

	printf("\n");
}

void ann_init::input_at_pos(ann_data *d) {
	int i, j, k, c;
	int flag;

	string line;
	vector<string> s, s_el;
	vector<int> s_na;
	ifstream fin;

	fin.open(d->at_infile.c_str(), ios::in);
	if (!fin) {
		cout << d->at_infile << " cannot be opened." << endl;
		exit(1);
	}

	getline(fin, line);
	getline(fin, line);
	for (i=0; i<3; i++) {
		getline(fin, line);
		s = split(line);
		d->lv[i][0] = atof(s[0].c_str());
		d->lv[i][1] = atof(s[1].c_str());
		d->lv[i][2] = atof(s[2].c_str());
	}

	getline(fin, line);
	s_el = split(line);

	for (i=0; i<s_el.size(); i++) {
		flag = d->nspecies;

		for (j=0; j<d->nspecies; j++) {
			if (s_el[i] != d->species[j]) flag--;
		}

		if (flag == 0) {
			printf("Undefined species is cotained in at_infile.\n");
			MPI_Abort(MPI_COMM_WORLD, 99);
		}
	}

	getline(fin, line);
	s = split(line);
	for (i=0; i<s.size(); i++) s_na.push_back(atoi(s[i].c_str()));

	d->natom = 0;
	for (i=0; i<s_na.size(); i++) {
		d->natom += s_na[i];
	}
	d->natom3 = 3 * d->natom;

	d->el = new string[d->natom];
	d->e = new int[d->natom];
	d->x = new double[d->natom];
	d->y = new double[d->natom];
	d->z = new double[d->natom];
    d->mass = new double[d->natom];

	c = 0;
	for (i=0; i<s_na.size(); i++) {
		for (j=0; j<d->nspecies; j++) {
			if (s_el[i] == d->species[j]) {
				k = d->species_id[j];
				break;
			}
		}

		for (j=0; j<s_na[i]; j++) {
			d->e[c] = k;
			d->el[c] = s_el[i];
			c++;
		}
	}

	for (i=0; i<d->nspecies; i++) d->e_natom[i] = 0;
	for (i=0; i<d->natom; i++) d->e_natom[d->e[i]]++;

	//if (d->rank == 0) {
	//	for (i=0; i<d->nspecies; i++) printf("e_natom = %d\n", d->e_natom[i]);
	//}

	//if (d->rank == 0) {
	//	for (i=0; i<d->natom; i++) {
	//		printf("e[%d] = %s %d\n", i, d->el[i].c_str(), d->e[i]);
	//	}
	//}

	if (d->rank == 0) {
		if (d->sel_flag == "yes") {
			d->mflx = new int[d->natom];
			d->mfly = new int[d->natom];
			d->mflz = new int[d->natom];
		}
	}

	getline(fin, line);
	s = split(line);
	if (s[0] == "Direct") {}
    else if (s[0] == "direct") {}
    else{
		cout << "not direct." << endl;
		exit(1);
	}

	if (d->sel_flag == "no") {
		for (i=0; i<d->natom; i++) {
			getline(fin, line); s = split(line);
			d->x[i] = atof(s[0].c_str());
			d->y[i] = atof(s[1].c_str());
			d->z[i] = atof(s[2].c_str());
		}
	} else if (d->sel_flag == "yes") {
		for (i=0; i<d->natom; i++) {
			getline(fin, line); s = split(line);

			if (s.size() < 6) {
				cout << "sel_flag is yes, but the vector size is less than 6." << endl;
				MPI_Abort(MPI_COMM_WORLD, 81);
			}

			d->x[i] = atof(s[0].c_str());
			d->y[i] = atof(s[1].c_str());
			d->z[i] = atof(s[2].c_str());

			if (d->rank == 0) {
				d->mflx[i] = atoi(s[3].c_str());
				d->mfly[i] = atoi(s[4].c_str());
				d->mflz[i] = atoi(s[5].c_str());
			}
		}
	}

	//if (d->rank == 0) {
	//	for (i=0; i<3; i++) printf("%-15.7f%-15.7f%-15.7f\n", d->lv[i][0], d->lv[i][1], d->lv[i][2]);
	//	for (i=0; i<d->natom; i++) printf("%10d%15.7f%15.7f%15.7f\n", d->e[i], d->x[i], d->y[i], d->z[i]);
	//}

	// wrap atoms
	for (i=0; i<d->natom; i++) {
		if (d->x[i] < 0.) { d->x[i] += 1.; }
		else if (d->x[i] >= 1.) { d->x[i] -= 1.; }

		if (d->y[i] < 0.) { d->y[i] += 1.; }
		else if (d->y[i] >= 1.) { d->y[i] -= 1.; }

		if (d->z[i] < 0.) { d->z[i] += 1.; }
		else if (d->z[i] >= 1.) { d->z[i] -= 1.; }
	}

	d->vx = new double[d->natom];
	d->vy = new double[d->natom];
	d->vz = new double[d->natom];

	if (d->rank == 0) {
		d->cx = new double[d->natom];
		d->cy = new double[d->natom];
		d->cz = new double[d->natom];
/*
		d->vx = new double[d->natom];
		d->vy = new double[d->natom];
		d->vz = new double[d->natom];
*/
		d->px = new double[d->natom];
		d->py = new double[d->natom];
		d->pz = new double[d->natom];

		for (i=0; i<d->natom; i++) {
			d->vx[i] = 0.;
			d->vy[i] = 0.;
			d->vz[i] = 0.;

			d->px[i] = 0.;
			d->py[i] = 0.;
			d->pz[i] = 0.;
		}

		if ((d->simulation_type == "md") || (d->simulation_type == "ti")) {
			if (d->init_vel == "read") {
				cout << "read velocity" << endl;
				getline(fin, line);

				for (i=0; i<d->natom; i++) {
					getline(fin, line); s = split(line);

					if (s.size() != 3) {
						printf("The number of lines for velocity is incorrect, or velocities are not contained.\n");
						MPI_Abort(MPI_COMM_WORLD, 86);
					}

					d->vx[i] = atof(s[0].c_str());
					d->vy[i] = atof(s[1].c_str());
					d->vz[i] = atof(s[2].c_str());
					//printf("%25.15f%25.15f%25.15f\n", d->vx[i], d->vy[i], d->vz[i]);
				}
			}
		}
	}

	fin.close();
}

void ann_init::create_cell_list(ann_data *d) {
	int i, j, s;
	int cell_n_size;
	double csina, asinb, bsing;
	//double ec = 1.2;
	double ec = 1.2;

	d->flagx = 0;
	d->flagy = 0;
	d->flagz = 0;

	if (d->rank == 0) {
		printf("# lt = %-15.7f%-15.7f%-15.7f%-15.7f%-15.7f%-15.7f\n", 
			d->lt[0], d->lt[1], d->lt[2], d->lt[3], d->lt[4], d->lt[5]);
	}

	s = (int) floor(d->lt[0] * sin(d->lt[5]*PI/180.) / (ec*d->cutoff));
	if (s == 0) { d->lc[0] = 1; }
	else { d->lc[0] = s; }

	s = (int) floor(d->lt[1] * sin(d->lt[3]*PI/180.) / (ec*d->cutoff));
	if (s == 0) { d->lc[1] = 1; }
	else { d->lc[1] = s; }

	s = (int) floor(d->lt[2] * sin(d->lt[4]*PI/180.) / (ec*d->cutoff));
	if (s == 0) { d->lc[2] = 1; }
	else { d->lc[2] = s; }

	if (d->rank == 0) {
		printf("# lc = %d %d %d\n", d->lc[0], d->lc[1], d->lc[2]);
	}

	d->lcyz = d->lc[1]*d->lc[2];
	d->lcxyz = d->lcyz*d->lc[0];

	d->cell_atom_list = new int*[d->lcxyz];
	d->cell_num = new int*[d->natom];
	d->cell_natom = new int[d->lcxyz];

	cell_n_size = (int) ceil((double)d->natom/(double)d->lcxyz);
	//cout << "# cell_n_size = " << cell_n_size << endl;

	for (i=0; i<d->lcxyz; i++) {
		d->cell_atom_list[i] = new int[d->max_lc_natom];
		//d->cell_atom_list[i] = new int[cell_n_size];

		for (j=0; j<d->max_lc_natom; j++) {
		//for (j=0; j<cell_n_size; j++) {
			d->cell_atom_list[i][j] = 0;
			//cout << d->cell_atom_list[i][j] << endl;
		}
	}

	for (i=0; i<d->natom; i++) {
		d->cell_num[i] = new int[3];
		for (j=0; j<3; j++) {
			d->cell_num[i][j] = 0;
		}
	}

	//for (i=0; i<3; i++) d->neighbor[i] = 1;
	d->neighbor[0] = (int) ceil(ec*d->cutoff / (d->lt[0]*sin(d->lt[5]*PI/180.)));
	if (d->neighbor[0] != 1) d->flagx = 1;

	d->neighbor[1] = (int) ceil(ec*d->cutoff / (d->lt[1]*sin(d->lt[3]*PI/180.)));
	if (d->neighbor[1] != 1) d->flagy = 1;

	d->neighbor[2] = (int) ceil(ec*d->cutoff / (d->lt[2]*sin(d->lt[4]*PI/180.)));
	if (d->neighbor[2] != 1) d->flagz = 1;

	//if (d->rank == 0) printf("# neighbor = %d %d %d\n", d->neighbor[0], d->neighbor[1], d->neighbor[2]);
}

void ann_init::calc_n_atproc(ann_data *d) {
	double remainder;
	d->all_n_atproc = new int[d->nprocs];

	remainder = d->natom % d->nprocs;
	if (remainder - (d->rank+1) < 0) {
		d->n_atproc = d->natom / d->nprocs;
	} else {
		d->n_atproc = d->natom / d->nprocs + 1;
	}

	MPI_Allgather(&d->n_atproc, 1, MPI_INT, &d->all_n_atproc[0], 1, MPI_INT, MPI_COMM_WORLD);

	d->n_start = 0;

	for (int i=0; i<d->rank; i++) {
		d->n_start += d->all_n_atproc[i];
	}

	MPI_Barrier(MPI_COMM_WORLD);

	printf("# %d: n_atproc=%d n_start=%d\n", d->rank, d->n_atproc, d->n_start);
}

void ann_init::alloc_nn(ann_data *d) {
	int i, j, k;

	d->out1 = new double*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->out1[i] = new double[d->n_hidden1];
		for (j=0; j<d->n_hidden1; j++) {
			d->out1[i][j] = 0.;
		}
	}

	d->out2 = new double*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->out2[i] = new double[d->n_hidden2];
		for (j=0; j<d->n_hidden2; j++) {
			d->out2[i][j] = 0.;
		}
	}

	d->d_out1 = new double*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->d_out1[i] = new double[d->n_hidden1];
		for (j=0; j<d->n_hidden1; j++) {
			d->d_out1[i][j] = 0.;
		}
	}

	d->d_out2 = new double*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->d_out2[i] = new double[d->n_hidden2];
		for (j=0; j<d->n_hidden2; j++) {
			d->d_out2[i][j] = 0.;
		}
	}

	d->dedg = new double*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->dedg[i] = new double[d->n_symmetry];
		for (j=0; j<d->n_symmetry; j++) {
			d->dedg[i][j] = 0.;
		}
	}

	d->weight0 = new double**[d->nspecies];
	d->weight1 = new double**[d->nspecies];
	d->weight2 = new double*[d->nspecies];

	d->bias0 = new double*[d->nspecies];
	d->bias1 = new double*[d->nspecies];
	d->bias2 = new double[d->nspecies];

	d->w1w2 = new double**[d->nspecies];

	d->w1w2_dout2 = new double[d->n_hidden1];

	for (i=0; i<d->nspecies; i++) {
		d->weight0[i] = new double*[d->n_symmetry];
		for (j=0; j<d->n_symmetry; j++) {
			d->weight0[i][j] = new double[d->n_hidden1];
			for (k=0; k<d->n_hidden1; k++) {
				d->weight0[i][j][k] = 0.;
			}
		}

		d->weight1[i] = new double*[d->n_hidden1];
		for (j=0; j<d->n_hidden1; j++) {
			d->weight1[i][j] = new double[d->n_hidden2];
			for (k=0; k<d->n_hidden2; k++) {
				d->weight1[i][j][k] = 0.;
			}
		}

		d->weight2[i] = new double[d->n_hidden2];
		for (j=0; j<d->n_hidden2; j++) {
			d->weight2[i][j] = 0.;
		}

		d->bias0[i] = new double[d->n_hidden1];
		for (j=0; j<d->n_hidden1; j++) {
			d->bias0[i][j] = 0.;
		}

		d->bias1[i] = new double[d->n_hidden2];
		for (j=0; j<d->n_hidden2; j++) {
			d->bias1[i][j] = 0.;
		}

		d->bias2[i] = 0.;

		d->w1w2[i] = new double*[d->n_hidden1];
		for (j=0; j<d->n_hidden1; j++) {
			d->w1w2[i][j] = new double[d->n_hidden2];
			for (k=0; k<d->n_hidden2; k++) {
				d->w1w2[i][j][k] = 0.;
			}
		}
	}
}

void ann_init::input_weight_bias(ann_data *d) {
	int i, j, k;
	string line;
	vector<string> s;

	ifstream fin;
	fin.open(d->wb_infile.c_str(), ios::in);

	if(!fin) {
		printf("%s cannot be opened.\n", d->wb_infile.c_str());
		exit(1);
	}

	for (i=0; i<6; i++) getline(fin, line);

	getline(fin, line);
	s = split(line);
	d->sf_mean = atof(s[s.size()-2].c_str());
	d->sf_std = atof(s[s.size()-1].c_str());

	for (i=0; i<d->nspecies; i++) {
		getline(fin, line); // species
		
		getline(fin, line);
		for (j=0; j<d->n_symmetry; j++) {
			for (k=0; k<d->n_hidden1; k++) {
				getline(fin, line);
				s = split(line);
				d->weight0[i][j][k] = atof(s[0].c_str());
			}
		}
		getline(fin, line);

		getline(fin, line);
		for (j=0; j<d->n_hidden1; j++) {
			for (k=0; k<d->n_hidden2; k++) {
				getline(fin, line);
				s = split(line);
				d->weight1[i][j][k] = atof(s[0].c_str());
			}
		}
		getline(fin, line);

		getline(fin, line);
		for (j=0; j<d->n_hidden2; j++) {
			getline(fin, line);
			s = split(line);
			d->weight2[i][j] = atof(s[0].c_str());
		}
		getline(fin, line);

		getline(fin, line);
		for (j=0; j<d->n_hidden1; j++) {
			getline(fin, line);
			s = split(line);
			d->bias0[i][j] = atof(s[0].c_str());
		}
		getline(fin, line);

		getline(fin, line);
		for (j=0; j<d->n_hidden2; j++) {
			getline(fin, line);
			s = split(line);
			d->bias1[i][j] = atof(s[0].c_str());
		}
		getline(fin, line);
		
		getline(fin, line);
		getline(fin, line);
		s = split(line);
		d->bias2[i] = atof(s[0].c_str());
		getline(fin, line);
	}

	fin.close();

	for (i=0; i<d->nspecies; i++) {
		for (j=0; j<d->n_hidden1; j++) {
			for (k=0; k<d->n_hidden2; k++) {
				d->w1w2[i][j][k] = d->weight1[i][j][k] * d->weight2[i][k];
			}
		}
	}
}

void ann_init::alloc_descriptor(ann_data *d) {
	int i, j, jj, k;

	//if (d->rank == 0) cout << "test5_1" << endl;
	d->cf_natom2 = new int[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->cf_natom2[i] = 0;
	}

	//if (d->rank == 0) cout << "test5_2" << endl;
	d->cf_pair1 = new int*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->cf_pair1[i] = new int[d->max_natom];
		for (j=0; j<d->max_natom; j++) {
			d->cf_pair1[i][j] = 0;
		}
	}

	//if (d->rank == 0) cout << "test5_3" << endl;
	/*
	d->x1 = new double*[d->n_atproc];
	d->y1 = new double*[d->n_atproc];
	d->z1 = new double*[d->n_atproc];

	for (i=0; i<d->n_atproc; i++) {
		d->x1[i] = new double[d->max_natom];
		d->y1[i] = new double[d->max_natom];
		d->z1[i] = new double[d->max_natom];

		for (j=0; j<d->max_natom; j++) {
			d->x1[i][j] = -1000.;
			d->y1[i][j] = -1000.;
			d->z1[i][j] = -1000.;
		}
	}
	*/

	d->dx1   = new double*[d->n_atproc];
	d->dy1   = new double*[d->n_atproc];
	d->dz1   = new double*[d->n_atproc];
	d->dis1  = new double*[d->n_atproc];
	d->dis12 = new double*[d->n_atproc];
	d->fc1   = new double*[d->n_atproc];
	d->fc1d  = new double*[d->n_atproc];
	d->gv2   = new double*[d->n_atproc];

	for (i=0; i<d->n_atproc; i++) {
		d->dx1[i]   = new double[d->max_natom];
		d->dy1[i]   = new double[d->max_natom];
		d->dz1[i]   = new double[d->max_natom];
		d->dis1[i]  = new double[d->max_natom];
		d->dis12[i] = new double[d->max_natom];
		d->fc1[i]   = new double[d->max_natom];
		d->fc1d[i]  = new double[d->max_natom];
		d->gv2[i]   = new double[d->max_natom];

		for (j=0; j<d->max_natom; j++) {
			d->dx1[i][j]   = 0.;
			d->dy1[i][j]   = 0.;
			d->dz1[i][j]   = 0.;
			d->dis1[i][j]  = 0.;
			d->dis12[i][j] = 0.;
			d->fc1[i][j]   = 0.;
			d->fc1d[i][j]  = 0.;
			d->gv2[i][j]   = 0.;
		}
	}

	d->gv3       = new double**[d->n_atproc];
	d->cosine    = new double**[d->n_atproc];
	//d->dx_cosine = new double**[d->n_atproc];
	//d->dy_cosine = new double**[d->n_atproc];
	//d->dz_cosine = new double**[d->n_atproc];

	for (i=0; i<d->n_atproc; i++) {
		d->gv3[i]       = new double*[d->max_natom];
		d->cosine[i]    = new double*[d->max_natom];
		//d->dx_cosine[i] = new double*[d->max_natom];
		//d->dy_cosine[i] = new double*[d->max_natom];
		//d->dz_cosine[i] = new double*[d->max_natom];

		for (j=0; j<d->max_natom; j++) {
			jj = j + 1;
			d->gv3[i][j]       = new double[jj];
			d->cosine[i][j]    = new double[jj];
			//d->dx_cosine[i][j] = new double[jj];
			//d->dy_cosine[i][j] = new double[jj];
			//d->dz_cosine[i][j] = new double[jj];

			for (k=0; k<jj; k++) {
				d->gv3[i][j][k]       = 0.;
				d->cosine[i][j][k]    = 0.;
				//d->dx_cosine[i][j][k] = 0.;
				//d->dy_cosine[i][j][k] = 0.;
				//d->dz_cosine[i][j][k] = 0.;
			}
		}
	}

	d->pdx_cosine = new double**[d->n_atproc];
	d->pdy_cosine = new double**[d->n_atproc];
	d->pdz_cosine = new double**[d->n_atproc];

	for (i=0; i<d->n_atproc; i++) {
		d->pdx_cosine[i] = new double*[d->max_natom];
		d->pdy_cosine[i] = new double*[d->max_natom];
		d->pdz_cosine[i] = new double*[d->max_natom];

		for (j=0; j<d->max_natom; j++) {
			d->pdx_cosine[i][j] = new double[d->max_natom];
			d->pdy_cosine[i][j] = new double[d->max_natom];
			d->pdz_cosine[i][j] = new double[d->max_natom];

			for (k=0; k<d->max_natom; k++) {
				d->pdx_cosine[i][j][k] = 0.;
				d->pdy_cosine[i][j][k] = 0.;
				d->pdz_cosine[i][j][k] = 0.;
			}
		}
	}

	d->sfv  = new double*[d->n_atproc];
	for (i=0; i<d->n_atproc; i++) {
		d->sfv[i] = new double[d->n_symmetry];
		for (j=0; j<d->n_symmetry; j++) {
			d->sfv[i][j] = 0.;
		}
	}

	d->sfdx = new double**[d->n_atproc];
	d->sfdy = new double**[d->n_atproc];
	d->sfdz = new double**[d->n_atproc];

	for (i=0; i<d->n_atproc; i++) {
		d->sfdx[i] = new double*[d->max_natom];
		d->sfdy[i] = new double*[d->max_natom];
		d->sfdz[i] = new double*[d->max_natom];

		for (j=0; j<d->max_natom; j++) {
			d->sfdx[i][j] = new double[d->n_symmetry];
			d->sfdy[i][j] = new double[d->n_symmetry];
			d->sfdz[i][j] = new double[d->n_symmetry];

			for (k=0; k<d->n_symmetry; k++) {
				d->sfdx[i][j][k] = 0.;
				d->sfdy[i][j][k] = 0.;
				d->sfdz[i][j][k] = 0.;
			}
		}
	}

	/*
	d->sfdx_i = new double**[d->n_atproc];
	d->sfdy_i = new double**[d->n_atproc];
	d->sfdz_i = new double**[d->n_atproc];

	for (i=0; i<d->n_atproc; i++) {
		d->sfdx_i[i] = new double*[d->max_natom];
		d->sfdy_i[i] = new double*[d->max_natom];
		d->sfdz_i[i] = new double*[d->max_natom];

		for (j=0; j<d->max_natom; j++) {
			d->sfdx_i[i][j] = new double[d->n_symmetry];
			d->sfdy_i[i][j] = new double[d->n_symmetry];
			d->sfdz_i[i][j] = new double[d->n_symmetry];

			for (k=0; k<d->n_symmetry; k++) {
				d->sfdx_i[i][j][k] = 0.;
				d->sfdy_i[i][j][k] = 0.;
				d->sfdz_i[i][j][k] = 0.;
			}
		}
	}
	*/

	d->fcx = new double[d->natom];
	d->fcy = new double[d->natom];
	d->fcz = new double[d->natom];
	d->t_fcx = new double[d->natom];
	d->t_fcy = new double[d->natom];
	d->t_fcz = new double[d->natom];

	for (i=0; i<d->natom; i++) {
		d->fcx[i] = 0.;
		d->fcy[i] = 0.;
		d->fcz[i] = 0.;
		d->t_fcx[i] = 0.;
		d->t_fcy[i] = 0.;
		d->t_fcz[i] = 0.;
	}

	if (d->rank == 0) {
		d->pfcx = new double[d->natom];
		d->pfcy = new double[d->natom];
		d->pfcz = new double[d->natom];

		for (i=0; i<d->natom; i++) {
			d->pfcx[i] = 0.;
			d->pfcy[i] = 0.;
			d->pfcz[i] = 0.;
		}
	}

	/*
	d->fcx_ij = new double*[d->natom];
	d->fcx_ij[0] = new double[d->natom * d->natom];
	d->fcy_ij = new double*[d->natom];
	d->fcy_ij[0] = new double[d->natom * d->natom];
	d->fcz_ij = new double*[d->natom];
	d->fcz_ij[0] = new double[d->natom * d->natom];
	d->t_fcx_ij = new double*[d->natom];
	d->t_fcx_ij[0] = new double[d->natom * d->natom];
	d->t_fcy_ij = new double*[d->natom];
	d->t_fcy_ij[0] = new double[d->natom * d->natom];
	d->t_fcz_ij = new double*[d->natom];
	d->t_fcz_ij[0] = new double[d->natom * d->natom];

	for (i=0; i<d->natom; i++) {
		d->fcx_ij[i] = d->fcx_ij[0] + i*d->natom;
		d->fcy_ij[i] = d->fcy_ij[0] + i*d->natom;
		d->fcz_ij[i] = d->fcz_ij[0] + i*d->natom;
		d->t_fcx_ij[i] = d->t_fcx_ij[0] + i*d->natom;
		d->t_fcy_ij[i] = d->t_fcy_ij[0] + i*d->natom;
		d->t_fcz_ij[i] = d->t_fcz_ij[0] + i*d->natom;
	}
	*/

	/*
	d->sfv  = new double*[d->n_symmetry];
	d->sfv[0]  = new double[d->n_symmetry * d->n_atproc];

	for (i=0; i<d->n_symmetry; i++) {
		d->sfv[i] = d->sfv[0] + i * d->n_atproc;
	}

	for (i=0; i<d->n_symmetry; i++) {
		for (j=0; j<d->n_atproc; j++) {
			d->sfv[i][j] = 0.;
		}
	}

	d->sfdx = new double**[d->n_symmetry];
	d->sfdy = new double**[d->n_symmetry];
	d->sfdz = new double**[d->n_symmetry];
	d->sfdx[0] = new double*[d->n_symmetry * d->n_atproc];
	d->sfdy[0] = new double*[d->n_symmetry * d->n_atproc];
	d->sfdz[0] = new double*[d->n_symmetry * d->n_atproc];
	d->sfdx[0][0] = new double[d->n_symmetry * d->n_atproc * d->natom];
	d->sfdy[0][0] = new double[d->n_symmetry * d->n_atproc * d->natom];
	d->sfdz[0][0] = new double[d->n_symmetry * d->n_atproc * d->natom];

	for (i=0; i<d->n_symmetry; i++) {
		d->sfdx[i] = d->sfdx[0] + i * d->n_atproc;
		d->sfdy[i] = d->sfdy[0] + i * d->n_atproc;
		d->sfdz[i] = d->sfdz[0] + i * d->n_atproc;

		for (j=0; j<d->n_atproc; j++) {
			d->sfdx[i][j] = d->sfdx[0][0] + i * d->n_atproc * d->natom + j * d->natom;
			d->sfdy[i][j] = d->sfdy[0][0] + i * d->n_atproc * d->natom + j * d->natom;
			d->sfdz[i][j] = d->sfdz[0][0] + i * d->n_atproc * d->natom + j * d->natom;
		}
	}
	*/
}
