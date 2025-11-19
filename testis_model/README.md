# testis_model

This directory contains the package to simulate and perform inference on the 
neutral drift mathematical model.

## Installation

Move into the parent directory (`sperm-data-analysis/`) and run
```pip install .```
Alternatively, the package can be installed in editable mode using
```pip install -e .```

## Subpackages

This package contains four subpackages:
* `data_processing`: Contains functions used to manipulate both simulated and
experimental data.
* `gillespie`: Contains the `AdvancedTestisModel` class that is used to run
Gillespie simulations of the model.
* `likelihood`: Contains the `LikelihoodModel` class that is used to estimate
rates of genetic drift.
* `visualize`: Contains code to create Muller plots of model simulations.

An detailed example of the model code being used is given in
`sperm-data-analysis/notebooks/tour`.