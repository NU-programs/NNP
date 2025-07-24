#ifndef _ANN_OUTPUT
#define _ANN_OUTPUT

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
#include "ann_data.hpp"
using namespace std;

class ann_output {
public:
	void output_at_pos(ann_data *d, int n);
	void output_weight_bias(ann_data *d);
};
#endif
