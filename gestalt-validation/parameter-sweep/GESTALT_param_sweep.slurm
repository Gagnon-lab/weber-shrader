#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH -o slurmjob-%j.out-%N
#SBATCH -e slurmjob-%j.err-%N
#SBATCH --account=gagnon
#SBATCH --partition=lonepeak
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=u0771757@umail.utah.edu

GROUPDATA=/uufs/chpc.utah.edu/common/home/gagnon-group2
# Change this to your data directory
DATA_DIR=/uufs/chpc.utah.edu/common/home/gagnon-group2/data/drew/Small_Sample_Test

cd $DATA_DIR

module load singularity
chmod u+x run_gestalt_parameter_sweep.py

python3 run_gestalt_parameter_sweep.py $GROUPDATA $DATA_DIR