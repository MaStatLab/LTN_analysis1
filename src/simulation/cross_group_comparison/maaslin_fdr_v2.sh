#!/usr/bin/bash 

notu=100
export WORK_DIR=...
INPUT_DIR_maaslin=$WORK_DIR/cache/cross_group_comparison/otu${notu}/MaAsLin2/
J_grid=(33)
nj_grid=(10 20)
signal_grid=(0.5 0.75 1 2 4)
lambda_grid=(10)
ntop=20
nseed=500
for J in ${J_grid[@]};do
for nj in ${nj_grid[@]};do
for lambda in ${lambda_grid[@]};do
H="H0"
$WORK_DIR/src/simulation/cross_group_comparison/maaslin_fdr_v2.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_maaslin --J $J --nj $nj --H $H --nseed $nseed
for signal in ${signal_grid[@]};do
H="H1_single"
$WORK_DIR/src/simulation/cross_group_comparison/maaslin_fdr_v2.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_maaslin --J $J --nj $nj --H $H --nseed $nseed --signal $signal --ntop $ntop 
H="H1_multi"
$WORK_DIR/src/simulation/cross_group_comparison/maaslin_fdr_v2.R --notu $notu --WORK_DIR $WORK_DIR --INPUT_DIR $INPUT_DIR_maaslin --J $J --nj $nj --H $H --nseed $nseed --signal $signal --ntop $ntop 
echo ${J},${nj},${signal},${H}
done 
done
done
done
