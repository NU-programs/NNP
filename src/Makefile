#!/bin/bash

CC       = mpiicpx
CXXFLAGS = -O3 -std=c++11
INCLUDE  =

DEPS     = ann_data.hpp ann_coordinate.hpp ann_mole.hpp ann_init.hpp ann_sf.hpp ann_alg.hpp ann_opt.hpp \
           ann_md.hpp ann_nvt_nh.hpp ann_npt_nh.hpp ann_nvt_vsc.hpp ann_nvt_lv.hpp ann_ti.hpp ann_output.hpp

OBJ      = main.o ann_mole.o ann_init.o ann_sf.o ann_alg.o ann_coordinate.o ann_opt.o ann_opt_cg.o ann_opt_qn.o \
           ann_md.o ann_nvt_lv.o ann_nvt_nh.o ann_npt_nh.o ann_nvt_vsc.o ann_ti.o ann_output.o

%.o: %.cpp $(DEPS)
	$(CC) $(CXXFLAGS) -c -o $@ $<

ann_mole: $(OBJ)
	$(CC) $(CXXFLAGS) $(INCLUDE) -o $@ $^

clean:
	rm -f $(OBJ)
