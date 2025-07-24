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
#include "ann_output.hpp"
using namespace std;

void ann_output::output_at_pos(ann_data *d, int n) {
	int i, j, k;
	FILE *fout;
	char fn[256];

	//double t = 0.;
	//for (i=0; i<d->natom; i++) t+= d->t_atomic_en[i];
	//printf("total_en = %15.7f\n", t);

	sprintf(fn, "NN_POSCAR.%08d.vasp", n);
	
	if ((fout = fopen(fn, "w")) != NULL) {
	} else {
		printf("fail to open %s.\n", fn);
		exit(1);
	}
	
	fprintf(fout, "test_at_pos\n");
	fprintf(fout, "1.\n");
	for (i=0; i<3; i++) {
		fprintf(fout, "%25.15f%25.15f%25.15f\n", d->lv[i][0], d->lv[i][1], d->lv[i][2]);
	}
	
	for (i=0; i<d->nspecies; i++) fprintf(fout, "%10s", d->species[i].c_str());
	fprintf(fout, "\n");
	for (i=0; i<d->nspecies; i++) fprintf(fout, "%10d", d->e_natom[i]);
	fprintf(fout, "\n");
	
	fprintf(fout, "Direct\n");
	if (d->sel_flag == "no") {
		for (i=0; i<d->nspecies; i++) {
			for (j=0; j<d->natom; j++) {
				if (d->e[j] == d->species_id[i]) {
					fprintf(fout, "%25.15f%25.15f%25.15f\n", d->x[j], d->y[j], d->z[j]);
				}
			}
		}
	
		fprintf(fout, "\n");
		for (i=0; i<d->nspecies; i++) {
			for (j=0; j<d->natom; j++) {
				if (d->e[j] == d->species_id[i]) {
					fprintf(fout, "%25.15f%25.15f%25.15f\n", d->vx[j], d->vy[j], d->vz[j]);
				}
			}
		}
	
		fprintf(fout, "\n");
		for (i=0; i<d->nspecies; i++) {
			for (j=0; j<d->natom; j++) {
				if (d->e[j] == d->species_id[i]) {
					fprintf(fout, "%25.15f%25.15f%25.15f\n", d->fcx[j], d->fcy[j], d->fcz[j]);
				}
			}
		}
	} else if (d->sel_flag == "yes") {
		for (i=0; i<d->nspecies; i++) {
			for (j=0; j<d->natom; j++) {
				if (d->e[j] == d->species_id[i]) {
					fprintf(fout, "%25.15f%25.15f%25.15f", d->x[j], d->y[j], d->z[j]);
					fprintf(fout, "%5d%5d%5d\n", d->mflx[j], d->mfly[j], d->mflz[j]);
				}
			}
		}
	
		fprintf(fout, "\n");
		for (i=0; i<d->nspecies; i++) {
			for (j=0; j<d->natom; j++) {
				if (d->e[j] == d->species_id[i]) {
					fprintf(fout, "%25.15f%25.15f%25.15f\n", d->vx[j], d->vy[j], d->vz[j]);
				}
			}
		}
	
		fprintf(fout, "\n");
		for (i=0; i<d->nspecies; i++) {
			for (j=0; j<d->natom; j++) {
				if (d->e[j] == d->species_id[i]) {
					fprintf(fout, "%25.15f%25.15f%25.15f\n", d->fcx[j], d->fcy[j], d->fcz[j]);
				}
			}
		}
	}

	fprintf(fout, "\n");
	fprintf(fout, "%15.7f%15.7f%15.7f\n", d->vstress[0], d->vstress[5], d->vstress[4]);
	fprintf(fout, "%15.7f%15.7f%15.7f\n", d->vstress[5], d->vstress[1], d->vstress[3]);
	fprintf(fout, "%15.7f%15.7f%15.7f\n", d->vstress[4], d->vstress[3], d->vstress[2]);
	
	fclose(fout);
}

void ann_output::output_weight_bias(ann_data *d) {
	int i, j, k;
	char fn[256] = "wb.test.dat";
	FILE *fout;
	
	fout = fopen(fn, "w");
	
	fprintf(fout, "# n_str     = %-10d\n", 1);
	fprintf(fout, "# n_h1      = %-10d\n", d->n_hidden1);
	fprintf(fout, "# n_h2      = %-10d\n", d->n_hidden2);
	fprintf(fout, "# n_sym     = %-10d\n", d->n_symmetry);
	fprintf(fout, "# nspecies = %-10d\n", d->nspecies);
	fprintf(fout, "# species   = ");
	for (i=0; i<d->nspecies; i++) {
		fprintf(fout, "%s ", d->species[i].c_str());
	}
	fprintf(fout, "\n");
	fprintf(fout, "energy_diff =\n");

	for (i=0; i<d->nspecies; i++) {
		fprintf(fout, "%s %d\n", d->species[i].c_str(), d->species_id[i]);

		fprintf(fout, "w0[n_sym][n_h1]\n");
		for (j=0; j<d->n_symmetry; j++) {
			for (k=0; k<d->n_hidden1; k++) {
				fprintf(fout, "%15.10f\n", d->weight0[i][j][k]);
			}
		}
		fprintf(fout, "\n");
	
		fprintf(fout, "w1[n_h1][n_h2]\n");
		for (j=0; j<d->n_hidden1; j++) {
			for (k=0; k<d->n_hidden2; k++) {
				fprintf(fout, "%15.10f\n", d->weight1[i][j][k]);
			}
		}
		fprintf(fout, "\n");
		
		fprintf(fout, "w2[n_h2]\n");
		for (j=0; j<d->n_hidden2; j++) {
			fprintf(fout, "%15.10f\n", d->weight2[i][j]);
		}
		fprintf(fout, "\n");
	
		fprintf(fout, "b0[n_h1]\n");
		for (j=0; j<d->n_hidden1; j++) {
			fprintf(fout, "%15.10f\n", d->bias0[i][j]);
		}
		fprintf(fout, "\n");
	
		fprintf(fout, "b1[n_h2]\n");
		for (j=0; j<d->n_hidden2; j++) {
			fprintf(fout, "%15.10f\n", d->bias1[i][j]);
		}
		fprintf(fout, "\n");
	
		fprintf(fout, "b2\n");
		fprintf(fout, "%15.10f\n", d->bias2[i]);
		fprintf(fout, "\n");
	}
	
	fclose(fout);
}