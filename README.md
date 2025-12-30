# sperm-data-analysis
This repository contains all of the data cleaning, mathematical modeling, and
statistical analysis done by Weber et. al. (2026).

To install the Python environment, open the terminal, navigate to this directory and run:
```
conda env create --name <name> --file python_requirements.yml
conda activate <name>
pip install -e .
python -m ipykernel install --user --name neutral_model --display-name Python (neutral_model)
```

Helpful commands if this fails:
```
conda clean --packages --tarballs --yes
```

The main directories are as follows:
* `data/`: Contains the serial sperm data that is output from the GESALT pipeline.
* `experiments/`: Stores the output of the `pipeline.ipynb` data cleaning file.
* `model_simulations/`: Contains data from running simulations of the mathematical model.
* `testis_model/`: Contains the python package to simulate the model and perform inference.
* `notebooks/`: Contains code for all later analysis, including statistical tests
and figure generation.

Each of these directories has its own `README.md` file with more details.