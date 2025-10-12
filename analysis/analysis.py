import os
import sys
import yaml
import subprocess
from data_io import read_data_file

def read_yaml_config(yaml_path):
    """Reads a YAML config file and returns the config dict."""
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)

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

    config = read_yaml_config(input_yaml_path)

    print("Successfully loaded input.yaml:")
    print(yaml.dump(config, indent=2))

    # --- Check for Output Files ---
    output_files = config.get('GaugeObservableParams', {})
    if not output_files:
        print("No output files specified in input.yaml")
        return

    for key, filename in output_files.items():
        if key.endswith('_filename'):
            file_path = os.path.join(directory, filename)
            if os.path.exists(file_path):
                print(f"Found '{key}': {filename}")
                # --- Trigger Analysis ---
                if key == 'plaquette_filename':
                    print("  -> Reading plaquette data...")
                    df = read_data_file(file_path)
                    print(df.head())
            else:
                print(f"Warning: Did not find '{key}': {filename}")
    
    skip_steps = "250"
    # --- Dispatch R Analyses ---
    r_scripts = [
        ["Rscript", "analysis_SimulationLoggingParams_log_file.R", directory, skip_steps],
        ["Rscript", "analysis_Plaquette.R", directory, "125"],
        ["Rscript", "analysis_wilsonflow_tests.R", directory, "250"],
        ["Rscript", "analysis_W_temp.R", directory]
    ]

    # Set the R_LIBS_USER environment variable to match the library path in check_and_install_hadron.sh
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    r_libs_user = os.path.join(project_root, ".R/library")
    env = os.environ.copy()
    env["R_LIBS_USER"] = r_libs_user

    for r_script in r_scripts:
        try:
            print(f"Dispatching R analysis with {r_script}...")
            subprocess.run(
                r_script,
                check=True,
                env=env
            )
            print(f"R analysis with {r_script[0]} completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error: R analysis with {r_script[0]} failed with exit code {e.returncode}")
        except FileNotFoundError:
            print("Error: Rscript not found. Ensure R is installed and available in PATH.")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python analysis.py <output_directory>")
        sys.exit(1)

    output_dir = sys.argv[1]
    analyze_directory(output_dir)