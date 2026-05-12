# `sperm-data-analysis/notebooks`

This directory contains several notebooks to conduct the main analysis for the paper.

## Sub-directories
* `images/`: Images generated are outputted to this directory.
* `OUTPUT_MP_TQ_extracted/`: Data for `singleton_proportion.ipynb`.

Figure-generating files:
* `fig_4a.ipynb`: Figure 4A and S7C.
* `fig_4b.ipynb`: TODO
* `fig_4c.ipynb`: TODO
* `fig_4d.ipynb`: Figure 4D and S7G.
* `fig_s7a.ipynb`: Figure S7A
* `fig_s11a.ipynb`: Figure S11A
* `fig_s11b.ipynb`: Figure S11B
* `tour.ipynb`: Not used in the paper, this file just gives an overview of the
model interface.
* `replicate_validation.ipynb`: Validating the robustness of technical replicates for Figure S2.
* `pie_charts.ipynb`: Creating pie charts for Figures [*TODO*] and S1.
* `cosines.iypnb`: Cosine similarity plots for Figures [*TODO*] and S5.
* `cosine_s1.ipynb`: Cosine similarity plot for Figure [*TODO*], before removing barcodes that appear in multiple fish.
* `sample_simulations.ipynb`: Sample simulations for Figure 3B.
* `inference_strategy.iypnb`: Inferred drift rates for Figure 3C.
* `significant_drift.ipynb`: Counting clones with significant drift for Figures 3D and S7B.
* `drift_vs_samples.ipynb`: Plotting drift confidence vs num. samples for Figures 3E and S7D.
* `violin.ipynb`: Violin plots for Figures 4A and S7C.
* `correlation.ipynb`: Correlation plots for Figure 4B and S7E.
* `caterpillar.ipynb`: Caterpillar plots for Figure 4C and S7F.
* `survivorship.ipynb.ipynb`: Survivorship plots for Figure 4D and S7G.
* `anomaly_detection.ipynb`: Anomaly detection in Figure S8.
* `singleton_proportion.ipynb`: Singleton proportion plot in Figure S9.
* `gestalt_params.ipynb`: Gestalt parameter sensitivity in Figure S10.
* `pseudosperm.ipynb`: Pseudosperm parameter fit in FIgure S11A.
* `sperm_decay.ipynb`: Sperm decay parameter fit in Figure S11B.

Other files:
* `analysis.py`: Helper functions and data loading that is used by most other notebooks in this directory.
* `cattelan_data.csv`: Data from Cattelan & Gasparini (2021) used in `sperm_decay.ipynb`.

## References
* Cattelan, S., & Gasparini, C. (2021). Male sperm storage impairs sperm quality in the zebrafish. Scientific reports, 11(1), 16689.
