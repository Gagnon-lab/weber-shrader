# weber-shrader
This repository contains all of the data cleaning, mathematical modeling, and
statistical analysis done by Weber et. al. (2026?).

## Sub-directories

The main directories are as follows:
* `experiments/`: Stores the output of the `pipeline.ipynb` data cleaning file.
* `gestalt/`: Code to run the GESTALT pipeline.
* `gestalt_output/`: Contains the serial sperm data that is output from the GESALT pipeline.
* `misc_data/`: Miscellaneous data files used for some results and figures.
* `model_simulations/`: Contains data from running simulations of the mathematical model.
* `notebooks/`: Contains code for all later analysis, including statistical tests
and figure generation.
* `testis_model/`: Contains the python package to simulate the model and perform inference.

Each of these directories has its own `README.md` file with more details.

The remaining files in this directory are:
* `.gitignore`: Describes files and directories to be ignored by `git`.
* `pipeline.ipynb`: The main computational pipeline to transform GESTALT output into final data.
* `pipeline_utility.py`: A collection of helper functions to streamline `pipeline.ipynb`.
* `requirements.yml`: The file read by `conda` to recreate the Python environment.
* `setup.py`: This file defines the `testis_model/` directory as a Python package so that the user can install it using `pip install .`

## Installation

To install the Python environment, we highly recommend using [Anaconda](https://anaconda.org/) and a Unix command line (if using Windows, consider using [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install)). Open the terminal, navigate to the repository directory and run:
```
conda env create --name weber --file environment.yml
conda activate weber
```
Note that the local `testis_model` package is installed in editable mode (meaning that changes you make to the files will affect the installation).

To open the files with Jupyter notebooks, you can then run, for example:
```
python -m ipykernel install --user --name weber --display-name "weber"
jupyter lab
```

## Reproducing our analysis

Our analysis code consists of 4 main steps:
* Running GESTALT on raw sequencing data (as in McKenna et al., 2016; Raj et al., 2018).
* A secondary data cleaning pipeline that was developed for this project.
* Computing clonal metrics from both simulated and experimental data (e.g. first/last time of appearance, drift strength).
* Statistical testing and figure generation.

Each step produces intermediate result files, so readers interested in simply re-generating the figures can skip ahead to that section.

## Downloading the raw sequencing data
A GEO link to download the raw `fastq` data will be added in the near future.

### Reproducing GESTALT barcode analysis pipeline
Download the `fastq` files using the instructions above, and place them into the `gestalt/` directory according to the sample number. For example, the data folder for sample 9 can be placed in the `gestalt/Longitudinal_Sample9/` directory.

Next, download the Docker container. The Center for High Performance Computing (CHPC) at the University of Utah offers Singularity as a container alternative. You can download the container as a Singularity image using:
```
singularity build GESTALT.sif docker://aaronmck/genomics:sc_GESTALT 
```
The other necessary files to run the GESTALT barcode analysis pipeline are included in the sample directories.

A shell script is included in each sample directory ('GESTALT_sample{number}.sh') that runs the GESTALT pipeline using the singularity container. Edit the following shell variables:
* GROUPDATA: should hold the directory holding GESTALT.sif
* DATA_DIR: should hold the enclosing directory, i.e. the absolute path to gestalt/sample#/
If running the GESTALT barcode analysis pipeline on an HPC cluster that requires job scheduling (e.g. Slurm, OpenPBS, etc.), further edit this file so that it is compatible with the scheduler. We found that 2 cpu cores were sufficient on CHPC at the University of Utah.

For each sample, run the script. This produces an `OUTPUT/` directory. After running on all samples, copy the `.allReadCounts` and `merged.info` files from HPC cluster to your local `gestalt_output/mp0.8_tq10`.

To run on other parameter values, edit the shell script 'run_gestalt_pipeline.sh' accordingly, changing the settings to `--matchProportion 0.95` or `--trimQual 14`.

### Secondary quality control pipeline
1. Make sure that the `gestalt_output/` directory contains the subdirectories `mp0.8_tq10`, `mp0.95_tq10`, and `mp0.8_tq14`, each containing `.allReadCounts` and `merged.info` files. These can either be downloaded from this repository or re-generated using the instructions above.
2. Install the Python environment using the steps above.
3. Run the pipeline. The `pipeline.ipynb` notebook contains a `preset` parameter that will configure different settings for the run:
    - `'full_pipeline'`: the default data cleaning pipeline.
    - `'full_pipeline_TQ14'`: the pipeline using the Trim Quality 14 GESTALT pipeline parameter.
    - `'full_pipeline_MP0.95'`: the pipeline using the Match Proportion 0.95 GESTALT pipeline parameter.
    - `'no_cluster'`: skip the clustering step.
    - `'keep_questionable'`: does not remove barcodes appearing in multiple fish.
    - `'no_filter'`: ignores both the used entry proportion from the GESTALT barcode analysis pipeline and the kept entry cutoff at the end of the secondary pipeline.
    - `'used_filter'`: ignores the kept entry cutoff at the end of the pipeline.
Note that if the results already exist in `experiments/` (for example, if you just cloned this repository), then you should either delete these results or change the `name` parameter above. The quality control code will not overwrite this data automatically.

Running this code will (re-)produce tables of final clonal reads in the `experiments/` subfolder.

### Clonal metric analysis
1. Move to the `model_simulations/` directory.
2. Run the notebooks `convert_experimental_data.ipynb` and `generate_simulated_data.ipynb`. This generates model data and formats the pipeline data to match the model output.
3. Run `model_simulations/python long_metrics.py --overwrite --scenarios 0 100`. This computes all the metrics about the experimental & simulated clone data, including drift strength estimates. This is computationally intensive and is best done on high performance computing resources. Using 16 parallel nodes, expect it to take 24-48 hours to run. This can be run in smaller batches by changing the `--scenarios [min] [max]` parameter.
4. Run `model_simulations/integrate_results.ipynb`. This generates `model_simulations/all_clone_results.csv` which contains all the results generated by `long_metrics.py` in a single table.

### Statistics and figure generation
1. Check that `model_simulations/all_clone_results.csv` exists and that you have installed the environment as described above.
2. Move to `notebooks/` and run any of the files there.

## References
* McKenna, A., Findlay, G. M., Gagnon, J. A., Horwitz, M. S., Schier, A. F., & Shendure, J. (2016). Whole-organism lineage tracing by combinatorial and cumulative genome editing. Science, 353(6298), aaf7907.
* Raj, B., Wagner, D. E., McKenna, A., Pandey, S., Klein, A. M., Shendure, J., ... & Schier, A. F. (2018). Simultaneous single-cell profiling of lineages and cell types in the vertebrate brain. Nature biotechnology, 36(5), 442-450.
