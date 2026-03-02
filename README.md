# sperm-data-analysis
This repository contains all of the data cleaning, mathematical modeling, and
statistical analysis done by Weber et. al. (2026).

## Sub-directories

The main directories are as follows:
* `gestalt/`: Code to run the GESTALT pipeline.
* `gestalt_output/`: Contains the serial sperm data that is output from the GESALT pipeline.
* `experiments/`: Stores the output of the `pipeline.ipynb` data cleaning file.
* `model_simulations/`: Contains data from running simulations of the mathematical model.
* `testis_model/`: Contains the python package to simulate the model and perform inference.
* `notebooks/`: Contains code for all later analysis, including statistical tests
and figure generation.

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
Note that the local `testis_model` package is installed in editable mode (meaning that changes you make will affect the installation).

To open the files with Jupyter notebooks, you can then run, for example:
```
python -m ipykernel install --user --name weber --display-name "weber"
jupyter lab
```

## Reproducing our analysis

Our analysis code consists of 4 main steps:
* Running GESTALT on raw sequencing data (McKenna et. al. 2016).
* A secondary data cleaning pipeline for this paper.
* Computing clonal metrics from both simulated and experimental data (e.g. first/last time of appearance, drift strength).
* Statistical testing and figure generation.

Each step produces intermediate result files, so readers interested in simply re-generating the figures can skip ahead to that section.

### Reproducing GESTALT runs
Download the fastq files from (TODO: link) and place them into the `gestalt/` directory according to the sample number. For example, the data folder for sample 9 can be placed in the `gestalt/sample9/` directory.

Next, download the Docker container. The Center for High Performance Computing (CHPC) at the University of Utah offers Singularity as a container alternative. You can download the container as a Singularity image using:
```
singularity build GESTALT.sif docker://aaronmck/genomics:sc_GESTALT 
```
The other necessary files to run GESTALT are included in the sample directories.

For each sample, you can run the bash (`.sh`) script file (or edit as necessary). This produces an `OUTPUT/` directory. After running on all samples, copy the `.allReadCounts` and `merged.info` files from HPC cluster to your local `gestalt_output/mp0.8_tq10`.

To run on other parameter values, edit the shell scripts accordingly, changing the settings to `--matchProportion 0.95` or `--trimQual 14`.


<!-- Install the GESTALT pipeline from [github.com/aaronmck/Cas9FateMapping](https://github.com/aaronmck/Cas9FateMapping) and the `pipeline/gestalt_pipeline_no_trees.scala` file from [github.com/Gagnon-lab/takasugi-genetics](https://github.com/Gagnon-lab/takasugi-genetics) -->

### Secondary quality control pipeline
1. Make sure that the `gestalt_output/` directory contains the subdirectories `mp0.8_tq10`, `mp0.95_tq10`, and `mp0.8_tq14`, each containing `.allReadCounts` and `merged.info` files. These can either be downloaded from this repository or re-generated using the instructions above.
2. Install the Python environment using the steps above.
3. Run the pipeline using the following settings:
    a. `name = 'full_pipeline'`, `dataset = 'mp0.8_tq10`, `cluster_barcodes = True`
    b. `name = 'full_pipeline_TQ14'`, `dataset = 'mp0.8_tq14`, `cluster_barcodes = True`
    c. `name = 'full_pipeline_MP0.95'`, `dataset = 'mp0.95_tq10`, `cluster_barcodes = True`
    d. `name = 'no_cluster'`, `dataset = 'mp0.8_tq10`, `cluster_barcodes = False`
Note that if the results already exist in `experiments/` (for example, if you just cloned this repository), then you should either delete these results or change the `name` parameter above. The quality control code will not overwrite this data automatically.

Running this code will (re-)produce tables of final clonal reads in the `experiments/` subfolder.

### Clonal metric analysis
1. Move to the `model_simulations/` directory.
2. Run the notebooks `convert_experimental_data.ipynb` and `generate_simulated_data.ipynb`. This generates model data and formats the pipeline data to match the model output.
3. Run `model_simulations/python long_metrics.py --overwrite --scenarios 0 100`. This computes all the metrics about the experimental & simulated clone data, including drift strength estimates. This is computationally intensive and is best done on high performance computing resources. Using 16 parallel nodes, expect it to take 24-48 hours to run. This can be run in smaller batches by changing the `--scenarios [min] [max]` parameter.
4. Run `model_simulations/integrate_results.ipynb`. This generates `model_simulations/all_clone_results.csv` which contains all the results generated by `long_metrics.py` in a single table.

### Statistics and figure generation
1. Check that `model_simulations/all_clone_results.csv` exists and that you have installed the environement as described above.
2. Move to `notebooks/` and run any of the files there.

## References
* McKenna, A., Findlay, G. M., Gagnon, J. A., Horwitz, M. S., Schier, A. F., & Shendure, J. (2016). Whole-organism lineage tracing by combinatorial and cumulative genome editing. Science, 353(6298), aaf7907.
