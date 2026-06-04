# `sperm-data-analysis/notebooks`

This directory contains several notebooks to conduct the main analysis for the paper.

## Sub-directories
* `images/`: Images generated are outputted to this directory.
* `OUTPUT_MP_TQ_extracted/`: Data for `singleton_proportion.ipynb`.

Figure-generating files:
* `anomaly_detection.ipynb`: Anomaly detection in Figure S8.
* `barcodes_per_fish.ipynb`: Swarm plots for Figures 1E and 2G.
* `barcodes_s1.ipynb`: Figure 1F.
* `caterpillar.ipynb`: Caterpillar plots for Figure 4C and S7F.
* `cosine_s1.ipynb`: Cosine similarity plot for Figure 1C, before removing barcodes that appear in multiple fish.
* `cosines.iypnb`: Cosine similarity plots for Figures 1H, 1I, and S5.
* `correlation.ipynb`: Correlation plots for Figure 4B and S7E.
* `drift_vs_samples.ipynb`: Plotting drift confidence vs num. samples for Figures 3E and S7D.
* `gestalt_params.ipynb`: Gestalt parameter sensitivity in Figure S10.
* `inference_strategy.iypnb`: Inferred drift rates for Figure 3C.
* `model_bias.ipynb`: The bias of model estimated drift rates, for Supplemental Table S2.
* `muller.ipynb`: Muller plots in Figures 2B-F and S4.
* `pie_charts.ipynb`: Creating pie charts for Figures 1D and S1.
* `pseudosperm.ipynb`: Pseudosperm parameter fit in Figure S11A.
* `replicate_validation.ipynb`: Validating the robustness of technical replicates for Figure S2.
* `sample_simulations.ipynb`: Sample simulations for Figure 3B.
* `shannon_diversity.ipynb`: Shannon diversity plots for Figures 2J and S6.
* `significant_drift.ipynb`: Counting clones with significant drift for Figures 3D and S7B.
* `singleton_proportion.ipynb`: Singleton proportion plot in Figure S9.
* `sperm_decay.ipynb`: Sperm decay parameter fit in Figure S11B.
* `survivorship.ipynb.ipynb`: Survivorship plots for Figure 4D and S7G.
* `violin.ipynb`: Violin plots for Figures 4A and S7C.

Other files:
* `analysis.py`: Helper functions and data loading that is used by most other notebooks in this directory.
* `anomaly_detection.py`: Helper code for `anomaly_detection.ipynb`.
* `cattelan_data.csv`: Data from Cattelan & Gasparini (2021) used in `sperm_decay.ipynb`.

## References
* Cattelan, S., & Gasparini, C. (2021). Male sperm storage impairs sperm quality in the zebrafish. Scientific reports, 11(1), 16689.
