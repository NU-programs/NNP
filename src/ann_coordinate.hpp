#ifndef _ANN_OPT_MD
#define _ANN_OPT_MD

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

class ann_coordinate {
public:
	void frac_to_cart(ann_data *d);
	void cart_to_frac(ann_data *d);
	void calc_lattice_const(ann_data *d);
	void calc_ilv(ann_data *d);
};

#endif