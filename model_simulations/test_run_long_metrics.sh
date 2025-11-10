#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --account=gagnon
#SBATCH --partition=lonepeak
#SBATCH --output=logs/test.out
#SBATCH --error=logs/test.err

source activate tm
conda info
cd ~/Documents/sperm-data-analysis/model_simulations
python long_metrics.py --overwrite --nsims 2 --scenarios 0 3
