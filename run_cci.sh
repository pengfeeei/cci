#!/bin/bash

./bin/s1_dis $1
./bin/s2_dis_bf $1
./bin/s3_cor $1
./bin/s4_cor_bf $1
echo "dimension reduction ... "
Rscript script/s5_dr_pca.r $1 > temp_out/r_log.txt
echo "done."
./bin/s6_embedding_dis $1
./bin/s7_embedding_dis_bf $1
./bin/s8_combine_bf $1

# an example for visualization 
./bin/s9_readscore_for_plot $1
Rscript script/s10_plot_skeleton_ref.r > temp_out/r_log.txt
Rscript script/s11_plot_on_skeleton.r $1 > temp_out/r_log.txt
echo "plot cell types and abundance with top CCI score on lineage skeleton."
echo "CCI calculation finished. results could be found in result/ and temp_out/."

