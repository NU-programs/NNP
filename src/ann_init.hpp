#ifndef _ANN_INIT
#define _ANN_INIT

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

class ann_init {
public:
	vector<string> split (string line);

	void input_parameters(string fn, ann_data *d);
	void calc_const(ann_data *d);
	void check_parameters(ann_data *d);
	void print_parameters(ann_data *d);

	void input_at_pos(ann_data *d);
	void calc_lattice_const(ann_data *d);
	void calc_ilv(ann_data *d);
	void calc_n_atproc(ann_data *d);
	void alloc_nn(ann_data *d);
	void input_weight_bias(ann_data *d);
	void alloc_descriptor(ann_data *d);
	void create_cell_list(ann_data *d);
};

#endif
