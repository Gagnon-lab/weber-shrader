# `sperm-data-analysis/model_simulations`

This directory contains files needed to conduct a large, computationally intensive
simulation study. The general workflow is as follows:
* A large number of model simulations are done in `generate_simulated_data.ipynb`,
with results saved in the `scenarios/` sub-directory.
* Experimental data is converted into the same format as model simulations in
`convert_experimental_data.ipynb`, with results saved in the `scenarios/`
sub-directory.
* Drift rate parameters are estimated on all model and simulated data in
`long_metrics.py` ("long metrics" meaning metrics that take a long time to compute,
i.e. the max likelihood algorithm). This was run on University of Utah Center
for High-Performance Computing resources. Each individual scenario takes a few
hours to run. Results are saved in `scenarios/`.
* All of the individual results are collected in `integrate_results.ipynb` and
output to a single file, `all_clone_results.csv`. It is this file that is read
for most of the analysis in the paper, see the `sperm-data-analysis/notebooks` directory.

## Sub-directories
* `logs/`: Contains output logs from running the model inference.
* `scenarios/`: Contains all of the simulation scenarios that were run.

## Files
* `all_clone_results.csv`: All the results after doing maximum likelihood inference
on simulated/experimental data. This file is output from `integrate_results.ipynb`
* `convert_experimental_data.ipynb`: This file converts the processed experimental
data in `sperm-data-analysis/experiments` into a format suitable for model inference.
This converted data is saved in the `scenarios/` sub-directory.
* `generate_simulated_data.ipynb`: This file runs all the model simulations used
to compare to experimental data. Outputs are saved in `scenarios/`.
* `integrate_results.ipynb`: This file collects all the maximum likelihood
inference results (stored as one file per scenario in `scenarios/`) into one
file and outputs this to `all_clone_results.csv`.
* `long_metrics.py`: This file is run to perform maximum likelihood inference on
simulated or experimental scenarios. The general process is to run `python long_metrics.py -scenarios 0 100`.
This will perform maximum likelihood estimation on the first 100 sub-directories in `scenarios/`, which should be
all of them. More information and optional settings are available by running `python long_metrics.py -h`