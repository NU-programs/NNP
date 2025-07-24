#ifndef _ANN_NVT_NH
#define _ANN_NVT_NH

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

class ann_nvt_nh {
public:
	void update_position(ann_data *d);
	void update_friction(ann_data *d);
	void update_velocity(ann_data *d);
	void copy_prev_force(ann_data *d);
	void calc_current_temperature(ann_data *d);

	void calc_d_tensor(ann_data *d);
	void calc_heatflux(ann_data *d);
	void output_kappa(ann_data *d, int step);
	void calc_heatflux_mode(ann_data *d);
	void output_kappa_mode(ann_data *d, int step);
};

#endif
