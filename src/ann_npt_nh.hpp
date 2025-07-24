#ifndef _ANN_NPT_NH
#define _ANN_NPT_NH

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

class ann_npt_nh {
public:
	void calc_npt_const(ann_data *d);
	void calc_ba_fric_tr(ann_data *d);
	void calc_total_stress(ann_data *d);
	void calc_exp2(ann_data *d);

	void update_position(ann_data *d);
	void update_velocity(ann_data *d);
	void update_lv(ann_data *d);

	void update_th_friction(ann_data *d);
	void update_ba_friction(ann_data *d);
	void update_th_velocity(ann_data *d);
	void update_ba_velocity(ann_data *d);
	
	void calc_current_temperature(ann_data *d);
};

#endif
