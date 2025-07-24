#ifndef _ANN_SF
#define _ANN_SF

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

class ann_sf {
public:
	void init_array(ann_data *d);
	void calc_distance(ann_data *d);
	void calc_g2(ann_data *d);
	void calc_g3_term(ann_data *d);
	void calc_g3(ann_data *d);
};

#endif