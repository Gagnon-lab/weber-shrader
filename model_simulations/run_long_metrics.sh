#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --ntasks=16
#SBATCH --mem=16G
#SBATCH --account=gagnon
#SBATCH --partition=lonepeak
#SBATCH --output=logs/long_metrics.out
#SBATCH --error=logs/long_metrics.err

source activate tm_test
conda info
cd ~/Documents/sperm-data-analysis/model_simulations
python long_metrics.py --scenarios 0 14
