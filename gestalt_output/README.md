# `sperm-data-analysis/gestalt`

This directory contains output from the GESTALT pipeline. To re-generate this data, follow the instructions in the
main README file.

There should be three folders:
* `mp0.8_tq10`: The `.allReadCounts` and `merged.info` files from the standard GESTALT run (match proportion 0.8, trim quality 10).
* `mp0.8_tq10`: The `.allReadCounts` and `merged.info` files from the GESTALT run with increased trim quality (match proportion 0.8, trim quality 14).
* `mp0.8_tq10`: The `.allReadCounts` and `merged.info` files from the GESTALT run with increased match proportion (match proportion 0.95, trim quality 10).

Note that GESTALT outputs more result files, but only the `.allReadCounts` and `merged.info` files are needed for downstream analysis.