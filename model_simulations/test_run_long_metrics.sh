#!/usr/bin/env bash

#SBATCH --time=0:30:00
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --account=gagnon
#SBATCH --partition=lonepeak
#SBATCH --output=logs/long_metrics.out
#SBATCH --error=logs/long_metrics.err

source activate tm_test2
conda info
cd ~/Documents/sperm-data-analysis/model_simulations
python long_metrics.py --nsims 5 --scenarios 0 3 --overwrite
