#ifndef _ANN_DATA
#define _ANN_DATA

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
//#include "/usr/local/openmpi-1.8.8/include/mpi.h"
//#include "/home/yokoi/intel/impi/2019.2.187/intel64/include/mpi.h"
#include "mpi.h"
using namespace std;

#define PI     3.14159265359
#define kB     8.61733e-5 // [eV/K]
#define rad_to_deg(rad)   (rad * 180./PI)
#define deg_to_rad(deg)   (deg * PI/180.)
#define evang3_to_kbar    1602.18

class ann_data {
public:
	string simulation_type;

	int rank, nprocs;
	int n_atproc, n_init, n_start;
    int *all_n_atproc;

	string at_infile;
	string wb_infile;
	//int plx, ply, plz, nlx, nly, nlz;

	int *cf_natom2;
	int **cf_pair1;

	double **dx1, **dy1, **dz1, **dis1, **dis12, **fc1, **fc1d, **gv2;
	double ***gv3;
	double ***cosine;
	//double ***dx_cosine, ***dy_cosine, ***dz_cosine;
	double **sfv;
	double ***sfdx, ***sfdy, ***sfdz;
	//double ***sfdx_i, ***sfdy_i, ***sfdz_i;
	double ***pdx_cosine, ***pdy_cosine, ***pdz_cosine;

	int **g3x;
	int ng2, ng3;
	int e_ng2, e_ng3;
	int t_ng2, t_ng3, t_ng2g3;

	double cutoff, cutoff2;
	double tanh_const0, tanh_const1, tanh_const2;
	double ch_const1;
	int max_natom, max_natom2, max_natom_p;
	int max_lc_natom;

	double *x, *y, *z;
	double *px, *py, *pz;
	double *cx, *cy, *cz;
	double *cdtx, *cdty, *cdtz;
	double *fcx, *fcy, *fcz;
	double *t_fcx, *t_fcy, *t_fcz;
	double *pfcx, *pfcy, *pfcz;
	int *mflx, *mfly, *mflz;
	string *el;
	int *e;
	string *species;
	int *species_id;
	int *e_natom;
	int natom, natom3;
	int natom_moved, natom_fixed;
	int nspecies;
	double lt[6], lv[3][3], ilv[3][3], plv[3][3];

	// two-body force
	//double **fcx_ij, **fcy_ij, **fcz_ij;
	//double **t_fcx_ij, **t_fcy_ij, **t_fcz_ij;

	// cell link list
	int lc[3], flagx, flagy, flagz;
	int lcyz, lcxyz;
	int **cell_atom_list;
	int neighbor[3];
	int **cell_num, *cell_natom;

	// neural network
	int n_symmetry, n_hidden1, n_hidden2;
	double ***weight0, ***weight1, **weight2, ***w1w2;
	double **bias0, **bias1, *bias2;
	double **out1, **out2;
	double **d_out1, **d_out2;
	double **dedg;
	double *w1w2_dout2;
	double sf_mean, sf_std;

	double energy, p_energy, t_energy;

	//double slength_min, slength_max;
	double step_length;
	double step_length_lv; //optimum step in inner coodinate  2023.1.30
	double step_length_init;
	double decay_rate;
	double en_diff;
	double fc_diff;
	double stress_diff;
	int n_iteration;
	int n_write;
	int n_write_en;

	string energy_outfile;

	// optimization
	string opt_alg;

	double volume;                      
	double t_vstress[6], vstress[6];
	double dEdh[3][3];
	double pdEdh[3][3];
	double cdlv[3][3];
	string lv_update;
	
	//quasi newton;
	double **H, *dx;
	double H_init;
	double *yk, **ykdxT_T, **dxdxT, **H1; //calculation procass of Hessian inverse

	// molecular dynamics
	int md_step, vsc_step;
	double md_temp, c_temp;
	double md_init_temp, md_fin_temp;
	double md_dtemp;
	double md_time, md_time_h, md_time2_h;
	double *md_mass;
	double *mass;
	double *vx, *vy, *vz;
	double *gzi_x, *gzi_y, *gzi_z, *eta_x, *eta_y, *eta_z, *eg_x, *eg_y, *eg_z;
	double md_vscale;
	double md_sigma, md_friction;
	double fric, p_fric;
	double th_mass;
	double t_ke, c_ke, ke;
	double total_mass;

	// npt
	double ba_mass;
	double ba_fric[3][3], ba_exp1[3][3], ba_exp2[3][3];
	double ba_fric_tr;
	double md_pressure_x, md_pressure_cx;
	double md_pressure_y, md_pressure_cy;
	double md_pressure_z, md_pressure_cz;
	double vxa, vya, vza;
	double stress_t[3][3];
	double *mass_sr;
	double md_time_4, md_time_8;
	double md_time_4c, md_time_8c;
	double th_fric_c;
	double mcx, mcy, mcz;
	string npt_shear;

	double c_unit1, c_unit2;
	double md_const1, md_const2, md_const4, md_const5, md_const6, md_const4_fric;
	double *mass_const1, *mass_const2, *mass_const3, *mass_const4, *mass_const5, *mass_const6;
	double lv_sigma, time_h_fric;
	double box_m1[3], box_m2[3];

	string sel_flag;
	string init_vel;
	string md_condition;
	string nvt_stress;

	// velocity autocorrelation
	string velauto_flag;
	string velauto_at_flag;
	string velauto_outdir;
	int velauto_n;
	int velauto_step;
	int velauto_inv;
	double *tva, **at_tva;
	double *va0, **at_va0;
	double **ivx, **ivy, **ivz;

	// msd
	string msd_flag, msd_outdir;
	int msd_step, msd_n, msd_inv;
	double **msdt, **msdx, **msdy, **msdz;
	double **ix, **iy, **iz;

	// average position
	string average_flag, average_outfile;
	double *xa, *ya, *za;

	// distribution
	string dist_flag;
	string dist_outdir;
	string dist_infile;
	string dist_mode;
	string dist_id_infile;
	double dist_range[6];
	double dist_lim[2];
	double dist_dl;
	double dist_sigma;
	int dist_interv;
	int dist_natom;
	int *dist_id;
	int dist_ndata;
	int dist_nstep;
	double ****dist_xyz;
	double dist_c1, dist_c2;

	// metadynamics
	double xg, yg, zg;
	double *ds_dx, *ds_dy, *ds_dz;
	double *x_s_of_t, *y_s_of_t, *z_s_of_t;
	double dVg_ds_x, dVg_ds_y, dVg_ds_z;
	double gauss, plus_Vg;
	int md_time_number;
	int tau_G, NG;
	double gau_h, gau_w, meta_cutoff;
	int *cutoff_list;
	double *wh, tmpV, delT;
	int mt_atom;
	double *sy, *sz;

	// thermodynamic integration
	double *x0, *y0, *z0;
	double *disx, *disy, *disz;
	double *fcx_h, *fcy_h, *fcz_h;
	//double *fcx_ha, *fcy_ha, *fcz_ha;
	double **fc_const;
	double ti_lambda, ti_lambda_h;
	double ti_energy, ha_energy, ha_energy0, ti_diff;
	string ti_infile;
	string fc_infile;

	// thermal conductivity
	double **d_tensor;
	double d_average[6];
	double *at_en, *p_at_en;
	double **at_stress, **p_at_stress;
	double fext;
	string ptb_direc;
    int cjx, cjy, cjz; // component of x, y, z waves
	int n_write_pmd;
	double **at_hflux_avg;

    // mode-decomposed thermal conductivity
    string mode_decomp;
    string mode_infile;
    int nmode;
    double *frequencies;
    double *local_eigenvectors;
    double *at_hflux_mode_avg;
    double *local_at_hflux_mode_avg;
    int *all_n_mode_proc;
    int *all_n_mode_start;
};

#endif
