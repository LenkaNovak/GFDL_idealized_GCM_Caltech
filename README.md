# Caltech's Idealized GCM

- Caltech's idealized GCM is based on GFDL's [FMS](https://www.gfdl.noaa.gov/fms/)
- this is atruncation including only the spectral model with vertical finite-difference scheme (hybrid vertical coordinates) and the full description can be found in `idealized/src/atmos_spectral/documentation/spectral_core.pdf`

## General structure:
```
                                   fms_default_public
                                          |
                      ---------------------------------------------------------
                      |                   |                    |              |
                 idealized               exp                analysis         sub
                      |                   |                    |        ,_____|_____,
             (DO NOT MODIFY src)          |                       run_exp.sh     list_exp        
                                          |
                   ----------------------------------------------
                   |                      |                     | 
                 exp1                   exp2                  <your-exp-name>...
                                          |
           ------------------------------------------------------------------
           |                              |                                 |
         input                           run                             srcmods
           |                              |                                 |
   namelists, diag table             run_scripts         put here any modification to the source code;
                                                       it will automatically replace any sources in src
```
Experiments (`exp1`, `exp2`, `<your-exp-name>` ...) may differ by source code, and all modifications to the source code are located in the `srcmods` directory within the specific experiment. Once in that directory they overwrite any source code in the model directory (e.g. `~/fms_default_public/idealized/src/`). Every experiment can have several runs which share the same source code. Experiment sweeps are facilitated using the `run_exp.sh` script. Runs typically differ in namelist parameters.

Once an experiment is run, two more directories will be created on the level of fms_default_public:

1. `fms_tmp`: A temporary directory where the experiment will run and computations will be executed. All model output will be automatically exported to `fms_output` once the experiment is done. Note that for speed some of the computation uses a local `/scratch` on the compute nodes. Unless a job is resubmitted immediately, this directory can be deleted once the run is over. For future restarts, a restart file archive (.cpio) is copied to `fms_output`. 

2. `fms_output`: Contains three directories:
    - 2.1. `history`: `4xday` and `1x10days` output of the model (output frequency and diagnostics can be adjusted at exp?/input/diag_table).
    - 2.2. `logfiles`: logfiles from each output stage.
    - 2.3. `restart`: restart files for each output stage.

## General instructions for setting up an experiment

(adjustments for specific computer configurations will be needed).

1. `cd exp`
2. Create a duplicate directory of the default experiments. e.g.: `cp -rf default_idealized test_moist_1`.
3. Copy any modifications to the source code into the `srcmods` directory of the new experiment.
4. Adjust if needed any `namelists` (see source code for param ...`_nml`s) or `diag_tables` (for model variables search src code for `register_diag_field`) in the run or input directories within the new experiment tree. These diagnostics are computed during the run time, using the `diag_manager` module.
5. Customize the run file (e.g. `exp/<your-exp-name>/run/run_test`) for your directory names (e.g., `workdir`, `outdir`) and experiment parameters 
6. Make sure the correct `mkmf` file for your specific compiler is used (e.g. `idealized/bin/mkmf.template.ifc_hpc_mpi`). For `mppnccombine`, which combines output files from separate processors, a GNU C++ (gcc) compiler is used.
7. Submit your job with the script in the `<your-exp-name>/run/` directory (e.g. `sbatch run_exp`), OR from the `sub/` directory (useful for multiple-run experiments sweeping over many parameters). In the latter, submit `./run_exp.sh` with parameters user-modified defined in the `list_exp` file.
8. Output diagnostics are in `<your-exp-name>/output/combine/<your-time-segment>/<filename>.nc` (cange the path in the `run_exp` file if needed)

## Analysis program
The directory analysis contains code for an offline computation of various kinds of averages (isentropic, sigma coordinate, etc.) of flow fields simulated in an experiment, after the model has run. This postprocessing step can be run using the module run file (e.g. `run_test`, which generates a `posprocessing_info` file with the relevant config info) or run separately after model is finished. The steps to run the offline version are:

1. ensure the model run has properly copied the postprocessing info, or alter the analysis runfile manually (e.g. from `analysis/analysis_3d/run/run_analysis_dry_3d_hpc` to `<your-scratch-path>/workdir/<date-stamp>/run_analysis_dry_3d_hpc` and `<your-scratch-path>/workdir/<date-stamp>/postprocessing_info`)
2. `offlinediag` in `OfflineDiag.f90` is the analysis program, and the executable analysis is built from it. Modify it to add/remove variables.
3. `cd <your-analysis-working-directory>` and `./run_analysis_dry_3d_hpc`
4. To run for all time segments, include this at the end of the model run script:
```
cd $run_analysis
echo "#\!/bin/bash" >> run_all
echo ' ' >> run_all
echo 'module load intel/18.1' >> run_all
echo 'for d in day*; do' >> run_all
echo '    test -d "${d}" && cd "${d}" && sbatch run_analysis_* && cd "../"'  >> run_all
echo 'done' >> run_all
```

- then do `chmod 775 run_all`
- then `./run_all`

Farid Ait Chaalal, Ian Eisenman, Yohai Kaspi, Xavier Levine, Tim Merlis, April 2011--October 2012
Lenka Novak 2020

## Troubleshooting
- modify `VERY_LARGE_FILE_FREQ` parameter for very long runs, otherwise it errors due to storage limits.
- right now there seems to be a bug when INPUT not generated upon first compile

## References:
- https://www.gfdl.noaa.gov/fms/
- https://github.com/NOAA-GFDL/FMS


