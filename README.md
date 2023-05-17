# Description

This repository contains all datasets and source code needed to reproduce the simulation examples and case study in "Differential abundance analysis for microbiome compositions with logistic-tree normal models".

The folder "res" contains the raw data.

The folder "src" contains the source code. 
- src/data: code for data processing.
- src/simulation: code for numerical examples in the paper (section 3).
- src/application: code for case study in the paper (section 4).

The simulation study and case study were run on a computing cluster with SLURM job scheduler.

The data used in case study is available at https://diabimmune.broadinstitute.org/diabimmune/t1d-cohort. 

# Workflow

Set up the work directory first. 
```
export WORK_DIR=...
mkdir $WORK_DIR/cache
mkdir $WORK_DIR/results
```

## Data processing

Run the following for data preprocessing.

```
$WORK_DIR/src/data/process_data_full.R $WORK_DIR 100
$WORK_DIR/src/data/process_data_simulation.R $WORK_DIR 100
```

## Simulation study

For the cross-group comparison simulation:
1. Generate synthetic datasets with "src/simulation/cross_group_comparison/sim_v3.R" and specify simulation parameters as described in the file.
2. Fit LTN to the simulated datasets:
- Generate a configuration csv file, containing a dataframe with four columns: raw (the file name of dataset simulated in step 1), input (file name for intermdediate data), output_dir (output directory), output_pjap_dir (output directory for PJAP only). Each row represents a simulated dataset. See "src/simulation/cross_group_comparison/config_ltn.R" for an example.
- For each simulated dataset, run "src/simulation/cross_group_comparison/LTNinput_slurm_v2.R" then "src/simulation/cross_group_comparison/fit_ltn_v3.R" to fit LTN on each of the datasets. (This step were run on a computing cluster with SLURM job scheduler.)
3. Fit DirFactor to the simulated datasets:
- Generate a configuration csv file, containing a dataframe with two columns: input (the file name of dataset simulated in step 1), output (output directory). See "src/simulation/cross_group_comparison/dirfactor/config_dirfactor.R" for an example. 
- For each simulated dataset, run "src/simulation/cross_group_comparison/dirfactor/fit_dirfactor_v2.R" to fit DirFactor on each of the datasets. (This step were run on a computing cluster with SLURM job scheduler.)
4. Fit MaAsLin to the simulated datasets:
- Generate a configuration csv file, containing a dataframe with two columns: input (the file name of dataset simulated in step 1), output (output directory). See "src/simulation/cross_group_comparison/MaAsLin2/config_maaslin.R" for an example. 
- For each simulated dataset, run "src/simulation/cross_group_comparison/MaAsLin2/fit_v2.R" to fit DirFactor on each of the datasets. 
5. Collect and visualize results with the rest of the code files in the "src/simulation/cross_group_comparison/" folder.

Please replace the string "REPLACE_WITH_WORKDIR" in the code with the work directory. 

## Case study: DIABIMMUNE 

- "src/application/application_005.R": Fit LTN mixed effects model to the DIABIMMUNE data.
- "src/application/dirfactor_perm.R": Fit DirFactor to the DIABIMMUNE data.
- "src/application/dirfactor_results.R": Collect results from DirFactor.
- "src/application/maaslin_application_default.R": Fit MaAsLin to the DIABIMMUNE data and visualize the comparison between LTN and MaAsLin. 

Please replace the string "REPLACE_WITH_WORKDIR" in the code with the work directory. 






