# `sperm-data-analysis/testis_model`

This directory contains the package to simulate and perform inference on the 
neutral drift mathematical model.

## Installation

Move into the parent directory (`sperm-data-analysis/`) and run
```pip install .```
Alternatively, the package can be installed in editable mode using
```pip install -e .```

## Sub-packages

This package contains five subpackages:
* `data_processing.py`: Contains functions used to manipulate both simulated and
experimental data.
* `gillespie.py`: Contains the `AdvancedTestisModel` class that is used to run
Gillespie simulations of the model.
* `likelihood.py`: Contains the `LikelihoodModel` class that is used to estimate
rates of genetic drift.
* `metrics.py`: Various metrics that can be computed from the simulations. 
* `visualize.py`: Contains code to create Muller plots of model simulations.