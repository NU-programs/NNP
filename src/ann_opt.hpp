#ifndef _ANN_OPT
#define _ANN_OPT

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

class ann_opt {
public:
	void opt_init(ann_data *d);
	void copy_current_position(ann_data *d);
	void copy_previous_position(ann_data *d);
	void copy_current_force(ann_data *d);
	void copy_previous_force(ann_data *d);
	void update_cg_direction(ann_data *d);
	void update_position_cg(ann_data *d);
	void opt_init_qn(ann_data *d);
	void update_position_qn(ann_data *d);
	void approximate_hessian_inverse(ann_data *d);
	void reset_hessian(ann_data *d);
	int check_force_conv(ann_data *d);
	int check_stress_conv(ann_data *d);

	void calc_dEdh(ann_data *d);
	void copy_current_lv(ann_data *d);
	void copy_current_dEdh(ann_data *d);
	void copy_previous_lv(ann_data *d);
	void copy_previous_dEdh(ann_data *d);
	void update_cg_direction_h(ann_data *d);
	void update_lv_cg(ann_data *d);
};

#endif
