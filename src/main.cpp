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
#include "unistd.h"
#include "ann_mole.hpp"
//#include "/usr/local/openmpi-1.8.8/include/mpi.h"
//#include "/home/yokoi/intel/impi/2019.2.187/intel64/include/mpi.h"
#include "mpi.h"
using namespace std;

int main(int argc, char *argv[]) {
	int rank, nprocs;
	string input_file = argv[1];

	srand((unsigned) time(NULL));
	//srand(1);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	ann_mole *ob;
	ob = new ann_mole(rank, nprocs);

	ob->set_initial_condition(input_file);

	if (ob->d->simulation_type == "md") {
		ob->run_md();
	} else if (ob->d->simulation_type == "opt")  {
		if (ob->d->opt_alg == "qn") {
			ob->opt_qn();
		} else if ((ob->d->opt_alg == "cg")||(ob->d->opt_alg == "cg_conv")) {
			ob->opt_cg();
		} else if (ob->d->opt_alg == "cg_conp") {
			ob->opt_cg_conp();
		}
	} else if (ob->d->simulation_type == "ti") {
		ob->run_md();
	}

	/*
	if (ob->d->simulation_type == "md")        { ob->run_md(); }
	else if (ob->d->simulation_type == "ti")   { ob->run_md(); }
	else if (ob->d->simulation_type == "opt")  {
		if (ob->d->opt_alg == "qn") {
			ob->opt_qn();
		} else if (ob->d->opt_alg == "cg") {
			ob->opt_cg();
		}
	}
	else if (ob->d->simulation_type == "meta") { ob->run_meta(); }
	*/

	if (rank == 0) cout << "### correctly finish ###" << endl;

	MPI_Finalize();

	return 0;
}
