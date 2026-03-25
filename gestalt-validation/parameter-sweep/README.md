# `sperm-data-analysis/gestalt_validation`

This directory contains files needed to conduct parameter sweeps on the GESTALT pipeline. The program 'run_gestalt_parameter_sweep.py' runs the GESTALT pipeline for a set of parameters defined by the user. This code was used to generate the data used for the supplemental heatmap figure. 'run_gestalt_parameter_sweep.py' is run instead of 'run_gestalt_pipeline.sh' in the GESTALT workflow. This workflow requires the creation of tearsheet templates described below.

## Tearsheet Templates
*  The output directory of a tearsheet must be replaced with \*\*\*OUTPUT_FOLDER***, as seen in 'Sample6tearsheet_allReps_template.txt'
*  Tearsheets templates must point to fastq data as if it were a normal GESTALT run