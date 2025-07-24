#ifndef _ANN_ALG
#define _ANN_ALG

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

class ann_alg {
public:
	void zscore_normalization(ann_data *d);
	void forward_propagation(ann_data *d);
	void forward_propagation_stress(ann_data *d);
	void forward_propagation_pmd(ann_data *d);
};

#endif