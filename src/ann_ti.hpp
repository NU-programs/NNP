#ifndef _ANN_TI
#define _ANN_TI

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

class ann_ti {
public:
	vector<string> split (string line);
	void alloc(ann_data *d);
	void input_ti_pos(ann_data *d);
	void init(ann_data *d);
	void input_fc(ann_data *d);
	void disp(ann_data *d);
	void calc_fc(ann_data *d);
};

#endif
