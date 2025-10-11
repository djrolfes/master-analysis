import os
import sys
import yaml

def analyze_directory(directory):
    """
    Analyzes the output of a klft simulation in a given directory.
    """
    print(f"Analyzing directory: {directory}")

    # --- Read Input YAML ---
    input_yaml_path = os.path.join(directory, 'input.yaml')
    if not os.path.exists(input_yaml_path):
        print(f"Error: 'input.yaml' not found in {directory}")
        return

    with open(input_yaml_path, 'r') as f:
        config = yaml.safe_load(f)

    print("Successfully loaded input.yaml:")
    print(yaml.dump(config, indent=2))

    # --- Check for Output Files ---
    output_files = config.get('output_files', {})
    if not output_files:
        print("No output files specified in input.yaml")
        return

    for key, filename in output_files.items():
        file_path = os.path.join(directory, filename)
        if os.path.exists(file_path):
            print(f"Found '{key}': {filename}")
            # --- Trigger Analysis ---
            # You can add your analysis logic here based on the file's existence
            # For example, if 'data_filename' exists, run a plotting script.
            if key == 'data_filename':
                print("  -> Triggering analysis for data file...")
                # plot_data(file_path)
        else:
            print(f"Warning: Did not find '{key}': {filename}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python analysis.py <output_directory>")
        sys.exit(1)

    output_dir = sys.argv[1]
    analyze_directory(output_dir)