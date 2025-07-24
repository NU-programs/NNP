#ifndef _ANN
#define _ANN

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
#include "ann_init.hpp"
#include "ann_coordinate.hpp"
#include "ann_sf.hpp"
#include "ann_alg.hpp"
#include "ann_nvt_vsc.hpp"
#include "ann_nvt_nh.hpp"
#include "ann_npt_nh.hpp"
#include "ann_nvt_lv.hpp"
//#include "ann_nvt_meta.hpp"     //metadynamics
#include "ann_opt.hpp"
#include "ann_md.hpp"
#include "ann_output.hpp"
#include "ann_ti.hpp"
using namespace std;

class ann_mole {
public:
	ann_data       *d;
	ann_init       *init;
	ann_coordinate *coord;
	ann_sf         *sf;
	ann_alg        *alg;
	ann_opt        *opt;
	ann_md         *md;
	ann_nvt_vsc    *vsc;
	ann_nvt_nh     *nvt_nh;
	ann_npt_nh     *npt_nh;
	ann_nvt_lv     *nvt_lv;
	ann_ti         *ti;
	ann_output     *output;
	
	/*
	ann_coordinate *coord;
	ann_sf         *sf;
	ann_alg        *alg;
	ann_opt        *opt;
	ann_nvt_meta   *nvt_meta;        //metadynamics
	ann_output     *output;
	*/

	void set_initial_condition(string infile);
	
	void calc_dis_sf();
	void forward_propagation();
	void forward_propagation_stress();
	void forward_propagation_pmd();
	void calc_ti();

	void run_md();
	void output_md_result(int i);
	void md_nvt_vsc();
	void md_nvt_nh(int step);
	void md_npt_nh(int step);
	void md_nvt_lv();
	void md_nve(int step);
	void md_nvt_pmd(int i);
	void opt_cg();
	void opt_cg_conp();
	void opt_qn();
	
	//void run_meta();       //metadynamics
	//void md_nvt_meta();

	void output_at_pos(int n) { output->output_at_pos(d, n); }

	ann_mole(int rank, int nprocs)  {
		d = new ann_data;
		d->rank = rank;
		d->nprocs = nprocs;
		//printf("rank=%d, nprocs=%d\n", d->rank, d->nprocs);
	}

	~ann_mole() { delete d; }
};

#endif
