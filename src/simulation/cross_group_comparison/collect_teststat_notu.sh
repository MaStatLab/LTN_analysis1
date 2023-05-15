#!/usr/bin/bash 

notu=50
export WORK_DIR=/work/zw122/LTN_analysis1
INPUT_DIR_ltn=/work/zw122/LTN_analysis1/cache/cross_group_comparison/otu${notu}/LTNoutput/
INPUT_DIR_dirfactor=/work/zw122/LTN_analysis1/cache/cross_group_comparison/otu${notu}/dirfactor/
INPUT_DIR_maaslin=/work/zw122/LTN_analysis1/cache/cross_group_comparison/otu${notu}/MaAsLin2/
J_grid=(33)
nj_grid=(10 20)
signal_grid=(0.5 0.75 1 2 4)
lambda_grid=(10)
ntop=20
niter_ltn=10000
niter_dirfactor=100000
nseed=500
for J in ${J_grid[@]};do
for nj in ${nj_grid[@]};do
for lambda in ${lambda_grid[@]};do
H="H0"
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_ltn_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_ltn --J $J --nj $nj --niter $niter_ltn --lambda $lambda --H $H --nseed $nseed
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_dirfactor_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_dirfactor --J $J --nj $nj --niter $niter_dirfactor --r 5 --H $H --nseed $nseed
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_maaslin_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_maaslin --J $J --nj $nj --H $H --nseed $nseed
for signal in ${signal_grid[@]};do
H="H1_single"
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_ltn_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_ltn --J $J --nj $nj --signal $signal --ntop $ntop --niter $niter_ltn --lambda $lambda --H $H --nseed $nseed
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_dirfactor_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_dirfactor --J $J --nj $nj --niter $niter_dirfactor --r 5 --H $H --nseed $nseed --signal $signal --ntop $ntop 
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_maaslin_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_maaslin --J $J --nj $nj --H $H --nseed $nseed --signal $signal --ntop $ntop 
H="H1_multi"
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_ltn_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_ltn --J $J --nj $nj --signal $signal --ntop $ntop --niter $niter_ltn --lambda $lambda --H $H --nseed $nseed
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_dirfactor_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_dirfactor --J $J --nj $nj --niter $niter_dirfactor --r 5 --H $H --nseed $nseed --signal $signal --ntop $ntop 
/work/zw122/LTN_analysis1/src/simulation/cross_group_comparison/collect_teststat_maaslin_notu.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_maaslin --J $J --nj $nj --H $H --nseed $nseed --signal $signal --ntop $ntop 
echo ${J},${nj},${signal},${H}
done 
done
done
done
