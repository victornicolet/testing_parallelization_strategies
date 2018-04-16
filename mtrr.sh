#!/bin/sh
source ~/python35/bin/activate

cmake CMakeLists.txt
make
./parallel_strategies_testing
python plot_max_top_right_rectangle.py parallel_strategies_test_mtrr.csv

# Second set of experiments: flag off
cmake CMakeLists.txt -DDEFINE_MTRR_LWD_AUX=ON
make
./parallel_strategies_testing
python plot_max_top_right_rectangle.py parallel_strategies_test_mtrr_aux.csv aux_rw
