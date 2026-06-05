# `sperm-data-analysis/experiments`

This directory contains the outputs from running `pipeline.ipynb` in the main
directory. There should be seven sub-directories:

* `full_pipeline`: The standard data processing pipeline (used for the majority of results and figures).
* `full_pipeline_MP0.95`: Processed data with increased match proportion GESTALT parameter (pertains to Fig. S7, S10).
* `full_pipeline_TQ14`: Processed data with increased match trim quality GESTALT parameter (pertains to Fig. S7, S10).
* `keep_questionable`: Processed data that keeps all barcodes that appear in multiple fish (used for Fig. 1C, 1D, S2)
* `no_cluster`: Processed data that skips step C (clustering) (Fig. S7)
* `no_filter`: The standard data without the used/kept read filter (Fig. S8).
* `used_filter`: The standard data without the kept read filter (Fig. S8).