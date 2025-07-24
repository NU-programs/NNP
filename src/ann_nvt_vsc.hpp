#ifndef _ANN_NVT_VSC
#define _ANN_NVT_VSC

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

class ann_nvt_vsc {
public:
	void copy_prev_force(ann_data *d);
	void update_position(ann_data *d);
	void update_velocity(ann_data *d);
	void velocity_scaling(ann_data *d);
};

#endif