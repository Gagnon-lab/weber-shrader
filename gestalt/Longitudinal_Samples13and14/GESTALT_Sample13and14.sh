#!/bin/bash

SING_DIR=<Path to singularity container>
DATA_DIR=<Path to containing directory>/Longitudinal_Sample13and14

mkdir -p Longitudinal_Sample13and14-OUTPUT

# Load singularity if needed
# module load singularity

cd $DATA_DIR
chmod u+x run_gestalt_pipeline.sh

# my_data is a "mounting point" i.e. a way to map host file systems into the container
singularity exec --bind $DATA_DIR:/my_data $SING_DIR/GESTALT.sif /my_data/run_gestalt_pipeline.sh