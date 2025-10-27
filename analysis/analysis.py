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

def analyze_directory(directory, skip_steps=None):
    """
    Analyzes the output of a klft simulation in a given directory.
    If skip_steps is None, it will be determined from plaquette analysis.
    """
    print(f"Analyzing directory: {directory}")
    if skip_steps is not None:
        print(f"Using manually specified skip_steps: {skip_steps}")
    else:
        print("Skip steps will be determined from plaquette thermalization analysis")

    # --- Read Input YAML ---
    yaml_files = [f for f in os.listdir(directory) if f.endswith('.yaml')]
    if not yaml_files:
        print(f"Error: No .yaml file found in {directory}")
        return

    input_yaml_path = os.path.join(directory, yaml_files[0])
    print(f"Using YAML file: {yaml_files[0]}")

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
    
    # Set the R_LIBS_USER environment variable to match the library path in check_and_install_hadron.sh
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    r_libs_user = os.path.join(project_root, ".R/library")
    env = os.environ.copy()
    env["R_LIBS_USER"] = r_libs_user
    
    # --- Step 1: Run Plaquette Analysis First to Determine Thermalization ---
    print("\n" + "="*60)
    print("STEP 1: Running plaquette analysis for thermalization estimate")
    print("="*60)
    
    initial_skip = skip_steps if skip_steps is not None else 0
    plaquette_cmd = ["Rscript", "analysis_Plaquette.R", directory, str(initial_skip)]
    
    try:
        print(f"Running: {' '.join(plaquette_cmd)}")
        subprocess.run(plaquette_cmd, check=True, env=env)
        print("Plaquette analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: Plaquette analysis failed with exit code {e.returncode}")
        print("Continuing with other analyses using initial skip value...")
    except FileNotFoundError:
        print("Error: Rscript not found. Ensure R is installed and available in PATH.")
        return
    
    # --- Step 2: Read Recommended Skip Steps (if not manually specified) ---
    if skip_steps is None:
        recommended_skip_file = os.path.join(directory, "recommended_skip.txt")
        if os.path.exists(recommended_skip_file):
            with open(recommended_skip_file, 'r') as f:
                skip_steps = int(f.read().strip())
            print(f"\n*** Using recommended skip_steps from plaquette analysis: {skip_steps} ***\n")
        else:
            skip_steps = initial_skip
            print(f"\n*** Warning: Could not read recommended skip steps. Using default: {skip_steps} ***\n")
    else:
        print(f"\n*** Using manually specified skip_steps: {skip_steps} ***\n")
    
    # --- Step 3: Run Acceptance Analysis (may update skip_steps) ---
    print("\n" + "="*60)
    print(f"STEP 2: Running acceptance analysis with skip_steps={skip_steps}")
    print("="*60 + "\n")
    
    # Determine if skip_steps was manually provided
    was_manual_skip = skip_steps is not None and len(sys.argv) == 3
    acceptance_cmd = ["Rscript", "analysis_acceptance.R", directory, str(skip_steps), 
                     "TRUE" if was_manual_skip else "FALSE"]
    
    try:
        print(f"Running: {' '.join(acceptance_cmd)}")
        subprocess.run(acceptance_cmd, check=True, env=env)
        print("Acceptance analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: Acceptance analysis failed with exit code {e.returncode}")
        print("Continuing with other analyses...")
    except FileNotFoundError:
        print("Error: Rscript not found. Ensure R is installed and available in PATH.")
    
    # Re-read recommended_skip.txt in case acceptance updated it
    if not was_manual_skip:
        recommended_skip_file = os.path.join(directory, "recommended_skip.txt")
        if os.path.exists(recommended_skip_file):
            with open(recommended_skip_file, 'r') as f:
                updated_skip = int(f.read().strip())
            if updated_skip != skip_steps:
                print(f"\n*** Acceptance analysis updated skip_steps: {skip_steps} -> {updated_skip} ***\n")
                skip_steps = updated_skip
    
    # --- Step 4: Dispatch Remaining R Analyses ---
    print("\n" + "="*60)
    print(f"STEP 3: Running remaining analyses with skip_steps={skip_steps}")
    print("="*60 + "\n")
    r_scripts = [
        ["Rscript", "analysis_SimulationLoggingParams_log_file.R", directory, str(skip_steps)],
        ["Rscript", "analysis_wilsonflow_tests.R", directory, str(skip_steps)],
        # ["Rscript", "analysis_wilsonflow_improv.R", directory, str(skip_steps)],
        ["Rscript", "analysis_W_temp.R", directory],
        ["Rscript", "analysis_topological_charge.R", directory, str(skip_steps)],
        ["Rscript", "analysis_wilsonflow_details.R", directory, str(skip_steps)],
        ["Rscript", "analysis_ptbc_log.R", directory, str(skip_steps)],
    ]

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
    if len(sys.argv) not in (2, 3):
        print("Usage: python analysis.py <output_directory> [skip_steps]")
        print("  If skip_steps is not provided, it will be determined from plaquette analysis")
        sys.exit(1)

    output_dir = sys.argv[1]
    skip_steps = int(sys.argv[2]) if len(sys.argv) == 3 else None
    analyze_directory(output_dir, skip_steps)
