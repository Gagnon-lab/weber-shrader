# `sperm-data-analysis/gestalt`

This directory contains files needed to run the GESTALT barcode analysis pipeline. Before running the GESTALT pipeline, the `fastq` files should be downloaded from GEO (using the instructions in the main README).

This directory should contain folders:
* `Longitudinal_Samples1and2`
* `Longitudinal_Sample2`
* `Longitudinal_Sample3`
* `Longitudinal_Sample4`
* `Longitudinal_Sample5`
* `Longitudinal_Sample6`
* `Longitudinal_Sample7`
* `Longitudinal_Sample8`
* `Longitudinal_Sample9`
* `Longitudinal_Sample10`
* `Longitudinal_Sample11and12`
* `Longitudinal_Sample13and14`
* `Longitudinal_Sample15and16`
* `Longitudinal_Sample17and18`

Each of these folders should then contain the files:
* `gestalt_pipeline_no_trees.scala`: the GESTALT script.
* `GESTALT_Sample{number}.sh`: An example bash script to run the GESTALT pipeline.
* `run_gestalt_pipeline.sh`: Another script that directly runs the GESTALT pipeline.
* `Sample{number}tearsheet.txt`: metadata linking each `fastq` file to a fish, sample, replicate.
* `sperm.referenceseq.fa`: The reference sequence.
* `sperm.referenceseq.fa.cutSites`: The sites and cut positions.
* `sperm.referenceseq.fa.primers`: The primers.
