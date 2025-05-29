# Extracts the files necessary for the data visualization notebooks.

import sys
import os
import shutil
import pandas as pd

ENCLOSING_DIR = sys.argv[1]

script_dir = os.path.dirname(__file__)

# The csv should contain the name each output directory.
parameters_df = pd.read_csv(os.path.join(ENCLOSING_DIR, 'parameters.csv')).drop_duplicates().set_index('Directory name')

output_dir = os.path.join(script_dir, str.split(ENCLOSING_DIR, os.sep)[-1] + '_extracted')
os.makedirs(output_dir)

# Copy parameters csv
parameters_sourcepath = os.path.join(ENCLOSING_DIR, 'parameters.csv')
parameters_destpath = os.path.join(output_dir, 'parameters.csv')
shutil.copyfile(parameters_sourcepath, parameters_destpath)

for dir in parameters_df.index:
    # Create a directory in the new output directory
    extracted_info_dir = os.path.join(output_dir, dir)
    os.makedirs(extracted_info_dir)
    
    # Copy the merged.info file into the new output directory
    dir_inside_html = os.listdir(os.path.join(ENCLOSING_DIR, dir, 'html'))[0]
    merged_info_sourcepath = os.path.join(ENCLOSING_DIR, dir, 'html', dir_inside_html, 'merged.info')
    merged_info_destpath = os.path.join(output_dir, dir, 'merged.info')
    shutil.copyfile(merged_info_sourcepath, merged_info_destpath)
    
    # Create a list of subdirectories
    dirs = os.listdir(os.path.join(ENCLOSING_DIR, dir))
    
    # Remove all directories that don't hold the outputs we need
    dirs.remove('html')
    try:
        dirs.remove('.DS_Store')
    except:
        pass
        
    # Copy all allreadcounts files to the new output directory
    for subdir in dirs:
        sample_output_dir = os.path.join(output_dir, dir, subdir)
        os.makedirs(sample_output_dir)
        source_path = os.path.join(ENCLOSING_DIR, dir, subdir, f'{subdir}.allReadCounts')
        dest_path = os.path.join(sample_output_dir, f'{subdir}.allReadCounts')
        shutil.copyfile(source_path, dest_path)
