import pandas as pd
import os

def read_data_file(filepath):
    """Reads a whitespace-separated data file into a pandas DataFrame, ignoring '#' in the header line."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Data file not found: {filepath}")

    # Read the first line, strip leading '#' if present, then read the rest
    with open(filepath, 'r') as f:
        header = f.readline().lstrip('#').strip()
        columns = header.split()
        df = pd.read_csv(f, delim_whitespace=True, names=columns)

    # Ensure the first column is named 'hmc_step'
    if df.columns[0] != 'hmc_step':
        df.rename(columns={df.columns[0]: 'hmc_step'}, inplace=True)

    return df

# Example: Specific IO functions using YAML variable names
def read_data_plaquette_filename(config, base_dir):
    filename = get_output_file_map(config).get('plaquette_filename')
    if filename:
        return read_data_file(os.path.join(base_dir, filename))
    else:
        raise ValueError("plaquette_filename not found in config")

def read_data_action_density_filename(config, base_dir):
    filename = get_output_file_map(config).get('action_density_filename')
    if filename:
        return read_data_file(os.path.join(base_dir, filename))
    else:
        raise ValueError("action_density_filename not found in config")

# For files not present in example_data, e.g. W_temp_filename:
def read_data_W_temp_filename(config, base_dir):
    """
    Example data for W_temp_filename should be provided.
    This function may need to be amended for the actual data format.
    """
    filename = get_output_file_map(config).get('W_temp_filename')
    if filename:
        # Placeholder: amend as needed for actual data format
        return read_data_file(os.path.join(base_dir, filename))
    else:
        raise ValueError("W_temp_filename not found in config")