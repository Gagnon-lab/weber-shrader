# Pipeline configurations

Each YAML file is a complete, reproducible settings record for `pipeline.ipynb`.
In the notebook, set `config_name` to a file stem (for example,
`"full_pipeline"`) to load `configs/full_pipeline.yaml`.

To make a new run configuration, copy an existing YAML file, give it a
descriptive filename, update its `name` (the experiment output directory), and
change the desired settings. All settings used by the pipeline must remain at
the YAML document's top level.
