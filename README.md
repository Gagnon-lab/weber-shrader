# sperm-data-analysis
A common space to analyze GESTALT sperm data. 

The purpose of this repository is to provide a common space to share code to process, analyze, and visualize serial sperm sampling data.

To install the Python environment, open the terminal, navigate to this directory and run:
```
conda create --name <name> --file python_requirements.yml
conda activate <name>
pip install -e .
python -m ipykernel install --user --name neutral_model --display-name Python (neutral_model)
```

The main directories are as follows:
* `data`: The data output from the GESALT pipeline.
* `experiments`: Contains processed data.
* `model_simulations`: Contains data from running simulations of the mathematical model.
* `neutral_model`: Contains the python package to simulate the model and perform inference.
* `notebooks`: Contains code to generate paper figures.

```
sperm-data-analysis
├── data/
    ├── viz-Output/
        ├── fish1_s1_rep1/
            ├── fish1_s1_rep1.allReadCounts (input)
            ├── fish1_s1_rep1.passReadCounts (intermediate output from step A)
            ├── fish1_s1_rep1.filteredReadCounts (intermediate output from step B)
            ├── fish1_s1_rep1.filteredQualityReadCounts (output from step D)
            ├── fish1_s1_rep1_stepB_removed_barcodes.txt (related to .filteredReadCounts)
        ├── fish1_s1_rep2/
	    ├── fish1_s1_rep3/
	    ├── fish1_s2_rep1/
            .
	    .
	    .
├── experiments/
	├── 11.11.2024_01.56_PM/
        ├── matrix/
            ├── fish1_barcodeMatrix_3filters.tsv
            ├── fish5_barcodeMatrix_3filters.tsv
            ├── fish10_barcodeMatrix_3filters.tsv
            ├── fishctrl3_barcodeMatrix_3filters.tsv
        ├── full_job_report.txt
        ├── stepA_summary.txt
        ├── stepB_summary.txt
        ├── stepD_summary.txt
        ├── stepE_summary.txt
├── pipeline/
    ├── process_sperm_barcodes.ipynb
```
