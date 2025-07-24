#ifndef _ANN_MD
#define _ANN_MD

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

class ann_md {
public:
    vector<string> split(string line);
	void md_init(ann_data *d);
	void allocation_pmd(ann_data *d);
	void set_md_const(ann_data *d);
	void zero_velocity_center_of_mass(ann_data *d);

	void allocation_velauto(ann_data *d);
	void calc_velauto(ann_data *d, int step);
	void average_velauto(ann_data *d);
	void output_velauto(ann_data *d);

	void allocation_msd(ann_data *d);
	void calc_msd(ann_data *d, int step);
	void average_msd(ann_data *d);
	void output_msd(ann_data *d);

	void allocation_average_position(ann_data *d);
	void calc_average_position(ann_data *d);
	void output_average_position(ann_data *d);

	void input_dist_pos(ann_data *d);
	void allocation_dist(ann_data *d);
	void calc_dist(ann_data *d);
	void output_dist(ann_data *d);
};

#endif
