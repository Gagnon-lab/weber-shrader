import os
import sys

GROUPDATA = sys.argv[1]
DATA_DIR = sys.argv[2]

# List of arguments to run
match_proportion_values = [0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
match_count_values = [5, 15, 25, 35, 45]
trim_qual_values = [12, 13, 14, 15, 16, 17]
trim_window_values = [4]

# Use arguments. If true, will run gestalt for all above values. If false, uses only default the value.
use_MP = True
use_MC = False
use_TQ = True
use_TW = True

# Default arguments
default_MP = 0.8
default_MC = 25
default_TQ = 10
default_TW = 5

# A class holding the arguments for one GESTALT run
class Gestalt_parameter_set:
    def __init__(self, match_proportion, match_count, trim_quality, trim_window):
        self.match_proportion = match_proportion
        self.match_count = match_count
        self.trim_quality = trim_quality
        self.trim_window = trim_window

# Create the array of all combinations of parameters
if not use_MP:
    match_proportion_values = [default_MP]
if not use_MC:
    match_count_values = [default_MC]
if not use_TQ:
    trim_qual_values = [default_TQ]
if not use_TW:
    trim_window_values = [default_TW]
parameter_set_array = []
for mp in match_proportion_values:
    for mc in match_count_values:
        for tq in trim_qual_values:
            for tw in trim_window_values:
                set = Gestalt_parameter_set(mp, mc, tq, tw)
                parameter_set_array.append(set)

# Read in the template shell file for running the gestalt pipeline
source_file = "Templates/run_gestalt_pipeline_template.sh"
try:
    with open(source_file, 'r') as file:
        shell_template_file_content = file.read()
except FileNotFoundError:
    print(f"Error: The template file '{source_file}' was not found.")
    exit(1)

# Read in the tearsheet template file
source_file = "Templates/Sample6tearsheet_rep1_template.txt"
try:
    with open(source_file, 'r') as file:
        tearsheet_template_file_content = file.read()
except FileNotFoundError:
    print(f"Error: The template file '{source_file}' was not found.")
    exit(1)

batch_name = 'Longitudinal_Sample_6'
containing_directory = f'mkdir {batch_name}-OUTPUT'

os.system(f'mkdir {containing_directory}')
parameter_info_str = "Directory name,Match proportion,Match count,Trim quality,Trim window\n"
info_destination_file = "parameters.csv"

for parameter_set in parameter_set_array:
    MATCH_PROPORTION = str(parameter_set.match_proportion)
    MATCH_COUNT = str(parameter_set.match_count)
    TRIM_QUAL = str(parameter_set.trim_quality)
    TRIM_WINDOW = str(parameter_set.trim_window)
    EXP_NAME = f"{batch_name}_MP{MATCH_PROPORTION}_MC{MATCH_COUNT}_TQ{TRIM_QUAL}_TW{TRIM_WINDOW}"
    OUTPUT_FOLDER = EXP_NAME + '-OUTPUT'
    TEAR_SHEET = 'Sample6tearsheet_Fish38.txt'
    
    # Add parameter values to the CSV
    parameter_info_str += f"{OUTPUT_FOLDER},{MATCH_PROPORTION},{MATCH_COUNT},{TRIM_QUAL},{TRIM_WINDOW}\n"

    # Destination file paths
    shell_destination_file = "run_gestalt_pipeline.sh"
    tearsheet_destination_file = TEAR_SHEET

    # Store a temporary variable with the content from the template shell file
    fc = shell_template_file_content

    # Perform the replacements
    fc = fc.replace("***EXP_NAME***", EXP_NAME)
    fc = fc.replace("***TEAR_SHEET***", TEAR_SHEET)
    fc = fc.replace("***MATCH_PROPORTION***", MATCH_PROPORTION)
    fc = fc.replace("***OUTPUT_FOLDER***", OUTPUT_FOLDER)
    fc = fc.replace("***MATCH_COUNT***", MATCH_COUNT)
    fc = fc.replace("***TRIM_QUAL***", TRIM_QUAL)
    fc = fc.replace("***TRIM_WINDOW***", TRIM_WINDOW)
    
    # Write the new file (this will overwrite any existing file with the same name)
    try:
        with open(shell_destination_file, 'w') as file:
            file.write(fc)
        print(f"Successfully created '{shell_destination_file}' with custom arguments.")
    except IOError:
        print(f"Error: Could not write to file '{shell_destination_file}'.")
        exit(1)
        
    # Store a temporary variable with the content from the template tearsheet file
    fc = tearsheet_template_file_content

    # Perform the replacements
    fc = fc.replace("***OUTPUT_FOLDER***", OUTPUT_FOLDER)
    
    # Write the new file (this will overwrite any existing file with the same name)
    try:
        with open(tearsheet_destination_file, 'w') as file:
            file.write(fc)
        print(f"Successfully created '{tearsheet_destination_file}' with custom arguments.")
    except IOError:
        print(f"Error: Could not write to file '{tearsheet_destination_file}'.")
        exit(1)

    # my_data is a "mounting point" i.e. a way to map host file systems into the container
    singularity_command = f"singularity exec --bind {DATA_DIR}:/my_data {GROUPDATA}/LabResources/sc_GESTALT.sif /my_data/run_gestalt_pipeline.sh"

    move_folder_command = f"mv {DATA_DIR}/{OUTPUT_FOLDER} {DATA_DIR}/{containing_directory}/"

    os.system('chmod u+x run_gestalt_pipeline.sh')
    os.system(f"mkdir -p {OUTPUT_FOLDER}")
    os.system(singularity_command)
    os.system(move_folder_command)
os.system('rm run_gestalt_pipeline.sh')

# Write to the info file
try:
    with open(info_destination_file, 'w') as file:
        file.write(parameter_info_str)
    print(f"Successfully created '{tearsheet_destination_file}' with custom arguments.")
except IOError:
    print(f"Error: Could not write to file '{tearsheet_destination_file}'.")
    exit(1)

os.system(f'mv {info_destination_file} {DATA_DIR}/{containing_directory}/')
