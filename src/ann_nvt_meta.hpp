#ifndef _ANN_NVT_META
#define _ANN_NVT_META

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
using namespace std;

#define s3 sqrt(3.)
#define PI2 2.*3.1415926535

class ann_nvt_meta {
public:
	void copy_prev_position(ann_data *d);
	void cutoff_particle_list(ann_data *d);
	void initialize_lv(ann_data *d);
	void generate_gaussian_variable(ann_data *d);
	void update_velocity(ann_data *d);
	void update_position(ann_data *d);
	void copy_prev_force(ann_data *d);
	void calc_current_temperature(ann_data *d);
	void set_parameter(ann_data *d);
	void make_CV(ann_data *d);
	void partial_diff(ann_data *d);
	void save_s_of_t(ann_data *d);
	void calc_vias_potential(ann_data *d);
	void force_plus_vias(ann_data *d);
};

#endif