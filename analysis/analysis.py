import os
import sys
import yaml
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from data_io import read_data_file

def run_r_script(args):
    """Helper function to run a single R script (must be at module level for pickling)"""
    r_script, env = args
    script_name = r_script[1]  # Extract the script filename
    try:
        print(f"Starting: {script_name}")
        result = subprocess.run(
            r_script,
            check=True,
            env=env,
            capture_output=True,
            text=True
        )
        print(f"Completed: {script_name}")
        return (script_name, True, None)
    except subprocess.CalledProcessError as e:
        error_msg = f"Exit code {e.returncode}\nStderr: {e.stderr}"
        print(f"Failed: {script_name} - {error_msg}")
        return (script_name, False, error_msg)
    except FileNotFoundError:
        error_msg = "Rscript not found. Ensure R is installed and available in PATH."
        print(f"Failed: {script_name} - {error_msg}")
        return (script_name, False, error_msg)

def read_yaml_config(yaml_path):
    """Reads a YAML config file and returns the config dict."""
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)

def analyze_directory(directory, skip_steps="1250"):
    """
    Analyzes the output of a klft simulation in a given directory.
    """
    print(f"Analyzing directory: {directory} (skip={skip_steps})")

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
    
    # Create logs directory for individual analysis logs
    logs_dir = os.path.join(directory, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    print(f"Created logs directory: {logs_dir}")
    
    # --- Dispatch R Analyses ---
    r_scripts = [
        ["Rscript", "analysis_SimulationLoggingParams_log_file.R", directory, str(skip_steps)],
        ["Rscript", "analysis_Plaquette.R", directory, str(skip_steps)],
        ["Rscript", "analysis_action_density.R", directory, str(skip_steps)],
        ["Rscript", "analysis_wilsonflow_tests.R", directory, str(skip_steps)],
        ["Rscript", "analysis_wilson_temp.R", directory, str(skip_steps)],
        # ["Rscript", "analysis_wilsonflow_improv.R", directory, str(skip_steps)],
        ["Rscript", "analysis_W_temp.R", directory],
        ["Rscript", "analysis_topological_charge.R", directory, str(skip_steps)],
        ["Rscript", "analysis_wilsonflow_details.R", directory, str(skip_steps)],
        ["Rscript", "analysis_ptbc_log.R", directory, str(skip_steps)],
    ]

    # Set the R_LIBS_USER environment variable to match the library path in check_and_install_hadron.sh
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    r_libs_user = os.path.join(project_root, ".R/library")
    env = os.environ.copy()
    env["R_LIBS_USER"] = r_libs_user

    # Determine number of cores to use:
    # 1. Check SLURM_CPUS_PER_TASK environment variable (set by SLURM)
    # 2. Fall back to len(os.sched_getaffinity(0)) for cgroup-limited environments
    # 3. Fall back to multiprocessing.cpu_count() as last resort
    if "SLURM_CPUS_PER_TASK" in os.environ:
        num_cores = int(os.environ["SLURM_CPUS_PER_TASK"])
        print(f"Detected SLURM allocation: using {num_cores} cores")
    elif hasattr(os, "sched_getaffinity"):
        num_cores = len(os.sched_getaffinity(0))
        print(f"Detected available cores (sched_getaffinity): using {num_cores} cores")
    else:
        num_cores = multiprocessing.cpu_count()
        print(f"Using all available cores: {num_cores} cores")
    
    print(f"Running R analyses in parallel using {num_cores} workers...")

    # Prepare arguments for each R script (including env for pickling)
    script_args = [(r_script, env) for r_script in r_scripts]

    # Execute R scripts in parallel
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = {executor.submit(run_r_script, args): args[0] for args in script_args}
        
        results = []
        for future in as_completed(futures):
            script_name, success, error = future.result()
            results.append((script_name, success, error))
    
    # Print summary
    print("\n" + "="*60)
    print("Analysis Summary:")
    print("="*60)
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    print(f"Successful: {len(successful)}/{len(results)}")
    for script_name, _, _ in successful:
        print(f"  ✓ {script_name}")
    
    if failed:
        print(f"\nFailed: {len(failed)}/{len(results)}")
        for script_name, _, error in failed:
            print(f"  ✗ {script_name}")
            if error:
                print(f"    Error: {error}")
    print("="*60)

if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        print("Usage: python analysis.py <output_directory> [skip_steps]")
        sys.exit(1)

    output_dir = sys.argv[1]
    skip_steps = sys.argv[2] if len(sys.argv) == 3 else "1000"
    analyze_directory(output_dir, skip_steps)
