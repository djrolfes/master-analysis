#!/usr/bin/env python3
"""
Merge two simulation data directories into a new directory.

This script merges text data files from two simulation directories, ensuring
continuous step numbering across the merged dataset. It handles:
- Regular data files with a 'step' or 'hmc_step' column
- Simulation log files with optional rank suffixes
- YAML configuration files
- Skipping initial steps from the second directory
"""

import argparse
import os
import shutil
import sys
import yaml
from pathlib import Path


def find_step_column(header_line, delimiter=None):
    """Find the step column name and its index in the header."""
    # Try to detect delimiter
    if delimiter is None:
        if ',' in header_line:
            delimiter = ','
        else:
            delimiter = None  # whitespace
    
    if delimiter:
        columns = [c.strip() for c in header_line.strip().split(delimiter)]
    else:
        columns = header_line.strip().split()
    
    for step_name in ['step', 'hmc_step']:
        if step_name in columns:
            return step_name, columns.index(step_name), delimiter
    return None, None, delimiter


def read_data_file(filepath):
    """Read a data file and return header, data lines, and step column info."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Check if first line is a comment header with column names
    if lines and lines[0].strip().startswith('#'):
        first_line = lines[0].strip()[1:].strip()
        step_name_test, step_idx_test, delimiter_test = find_step_column(first_line)
        if step_name_test is not None:
            # First line is a column header comment
            header = lines[0]
            header_idx = 0
            preamble = lines[:1]
            data_lines = lines[1:]
            step_name, step_idx, delimiter = step_name_test, step_idx_test, delimiter_test
            return preamble, header, data_lines, step_name, step_idx, delimiter
    
    # Otherwise, find header (first non-comment line)
    header_idx = 0
    header = None
    for i, line in enumerate(lines):
        if not line.strip().startswith('#') and line.strip():
            header_idx = i
            header = line
            break
    
    if header is None:
        # All comments or empty
        return lines[:1] if lines else [], "", [], None, None, None
    
    # Normal case: header is not a comment
    step_name, step_idx, delimiter = find_step_column(header)
    preamble = lines[:header_idx]
    data_lines = lines[header_idx + 1:]
    
    return preamble, header, data_lines, step_name, step_idx, delimiter


def get_max_step(data_lines, step_idx, delimiter):
    """Extract the maximum step value from data lines."""
    max_step = 0
    for line in data_lines:
        if line.strip() and not line.strip().startswith('#'):
            if delimiter:
                parts = [p.strip() for p in line.split(delimiter)]
            else:
                parts = line.split()
            if len(parts) > step_idx:
                try:
                    step_val = int(parts[step_idx])
                    max_step = max(max_step, step_val)
                except ValueError:
                    continue
    return max_step


def get_min_step(data_lines, step_idx, delimiter):
    """Extract the minimum step value from data lines."""
    min_step = None
    for line in data_lines:
        if line.strip() and not line.strip().startswith('#'):
            if delimiter:
                parts = [p.strip() for p in line.split(delimiter)]
            else:
                parts = line.split()
            if len(parts) > step_idx:
                try:
                    step_val = int(parts[step_idx])
                    if min_step is None or step_val < min_step:
                        min_step = step_val
                except ValueError:
                    continue
    return min_step if min_step is not None else 0


def adjust_step_column(data_lines, step_idx, offset, skip_steps, delimiter):
    """
    Adjust step column values by adding offset and skipping initial steps.
    
    Args:
        data_lines: List of data line strings
        step_idx: Column index of the step column
        offset: Value to add to each step
        skip_steps: Number of initial steps to skip from dir2
        delimiter: Delimiter used in the file (e.g., ',' or None for whitespace)
    
    Returns:
        List of adjusted data lines
    """
    adjusted_lines = []
    for line in data_lines:
        if line.strip() and not line.strip().startswith('#'):
            if delimiter:
                parts = [p.strip() for p in line.split(delimiter)]
            else:
                parts = line.split()
            
            if len(parts) > step_idx:
                try:
                    step_val = int(parts[step_idx])
                    
                    # Skip if step is below skip_steps threshold
                    if step_val < skip_steps:
                        continue
                    
                    # Adjust step: subtract skip_steps, then add offset
                    adjusted_step = step_val - skip_steps + offset
                    parts[step_idx] = str(adjusted_step)
                    
                    if delimiter:
                        adjusted_lines.append(delimiter.join(parts) + '\n')
                    else:
                        adjusted_lines.append(' '.join(parts) + '\n')
                except ValueError:
                    # Keep non-numeric lines as-is
                    adjusted_lines.append(line)
        else:
            # Keep comment/empty lines
            adjusted_lines.append(line)
    
    return adjusted_lines


def merge_text_file(file1_path, file2_path, output_path, skip_steps):
    """Merge two text data files with continuous step numbering."""
    print(f"  Merging: {os.path.basename(file1_path)} + {os.path.basename(file2_path)}")
    
    # Read both files
    preamble1, header1, data1, step_name1, step_idx1, delim1 = read_data_file(file1_path)
    preamble2, header2, data2, step_name2, step_idx2, delim2 = read_data_file(file2_path)
    
    if step_name1 != step_name2 or step_idx1 != step_idx2:
        print(f"    WARNING: Step column mismatch ({step_name1} vs {step_name2})")
        if step_idx1 is not None:
            print(f"    Using step column from file1: {step_name1}")
    
    if step_idx1 is None:
        print(f"    ERROR: No step column found in {os.path.basename(file1_path)}")
        print(f"    Skipping this file pair")
        return
    
    # Use delimiter from file1
    delimiter = delim1
    
    # Get max step from first file
    max_step1 = get_max_step(data1, step_idx1, delimiter)
    min_step2 = get_min_step(data2, step_idx2, delim2)
    
    print(f"    Dir1: max_step = {max_step1}")
    print(f"    Dir2: min_step = {min_step2}, skip = {skip_steps}")
    
    # Find typical step interval from dir2 (after skipping)
    step_interval = 1
    step_values = []
    for line in data2:
        if line.strip() and not line.strip().startswith('#'):
            if delim2:
                parts = [p.strip() for p in line.split(delim2)]
            else:
                parts = line.split()
            if len(parts) > step_idx2:
                try:
                    step_val = int(parts[step_idx2])
                    if step_val >= skip_steps:
                        step_values.append(step_val)
                        if len(step_values) >= 10:  # Get enough samples
                            break
                except ValueError:
                    continue
    
    if len(step_values) >= 2:
        step_interval = step_values[1] - step_values[0]
    
    # Offset calculation: make first valid step from dir2 = max_step1 + step_interval
    # After skipping, the first step from dir2 will be at position: skip_steps (or first step >= skip_steps)
    # We want: (first_step_after_skip - skip_steps) + offset = max_step1 + step_interval
    # So: offset = max_step1 + step_interval - first_step_after_skip + skip_steps
    first_step_after_skip = min([s for s in step_values if s >= skip_steps]) if step_values else skip_steps
    offset = max_step1 + step_interval - first_step_after_skip + skip_steps
    
    print(f"    Step interval: {step_interval}, first step after skip: {first_step_after_skip}")
    print(f"    Offset for dir2: {offset}")
    
    # Adjust second file's steps
    data2_adjusted = adjust_step_column(data2, step_idx2, offset, skip_steps, delim2)
    
    print(f"    Dir1: {len(data1)} lines, Dir2: {len(data2_adjusted)} lines (after skip)")
    
    if len(data2_adjusted) == 0:
        print(f"    WARNING: No data from dir2 after skipping {skip_steps} steps")
    
    # Write merged file
    with open(output_path, 'w') as out:
        # Use preamble from first file
        out.writelines(preamble1)
        # Only write header if it's not already in preamble (i.e., if it's not a comment header)
        if not (preamble1 and preamble1[-1].strip().startswith('#')):
            if not header1.endswith('\n'):
                out.write(header1 + '\n')
            else:
                out.write(header1)
        out.writelines(data1)
        out.writelines(data2_adjusted)
    
    print(f"    -> Merged: {len(data1) + len(data2_adjusted)} total lines")


def get_yaml_filenames(yaml_path):
    """Extract all data filenames mentioned in the YAML configuration."""
    try:
        with open(yaml_path, 'r') as f:
            config = yaml.safe_load(f)
        
        filenames = set()
        
        # Add filenames from GaugeObservableParams
        if 'GaugeObservableParams' in config:
            params = config['GaugeObservableParams']
            for key, value in params.items():
                if key.endswith('_filename') and isinstance(value, str) and value.endswith('.txt'):
                    filenames.add(value)
        
        # Add filenames from SimulationLoggingParams
        if 'SimulationLoggingParams' in config:
            params = config['SimulationLoggingParams']
            if 'log_filename' in params and params['log_filename'].endswith('.txt'):
                filenames.add(params['log_filename'])
        
        # Add filenames from PTBCSimulationLoggingParams
        if 'PTBCSimulationLoggingParams' in config:
            params = config['PTBCSimulationLoggingParams']
            if 'log_filename' in params and params['log_filename'].endswith('.txt'):
                filenames.add(params['log_filename'])
        
        # Add filenames from FermionObservableParams
        if 'FermionObservableParams' in config:
            params = config['FermionObservableParams']
            for key, value in params.items():
                if key.endswith('_filename') and isinstance(value, str) and value.endswith('.txt'):
                    filenames.add(value)
        
        return filenames
    except Exception as e:
        print(f"    WARNING: Could not parse YAML file: {e}")
        return set()


def get_txt_files(directory, yaml_filenames=None):
    """Get .txt files in a directory, optionally filtered by YAML config."""
    txt_files = {}
    
    # Get all .txt files
    all_files = list(Path(directory).glob("*.txt"))
    
    for filepath in all_files:
        filename = filepath.name
        
        # Skip summary/analysis files (not raw data)
        if any(skip in filename for skip in ['_summary', '_bootstrap', '_autocorr']):
            continue
        
        # Check if it's a ranked simulation log
        # Format: simulation_log.rank<N>.txt or similar
        if '.rank' in filename:
            base = filename.split('.rank')[0] + '.txt'
            # Check if base file is in yaml (if filtering)
            if yaml_filenames and base not in yaml_filenames:
                continue
            
            base_name = filename.split('.rank')[0]
            rank = filename.split('.rank')[1].replace('.txt', '')
            if base_name not in txt_files:
                txt_files[base_name] = {}
            txt_files[base_name][rank] = filepath
        else:
            # Regular file - check if it's in yaml (if filtering)
            if yaml_filenames and filename not in yaml_filenames:
                continue
            
            # Regular file
            base = filename.replace('.txt', '')
            txt_files[base] = filepath
    
    return txt_files


def merge_directories(dir1, dir2, skip_steps):
    """
    Merge two simulation directories into a new merged directory.
    
    Args:
        dir1: Path to first directory (base data)
        dir2: Path to second directory (data to append)
        skip_steps: Number of steps to skip from dir2
    """
    dir1_path = Path(dir1).resolve()
    dir2_path = Path(dir2).resolve()
    
    if not dir1_path.exists():
        print(f"ERROR: Directory 1 does not exist: {dir1_path}")
        sys.exit(1)
    
    if not dir2_path.exists():
        print(f"ERROR: Directory 2 does not exist: {dir2_path}")
        sys.exit(1)
    
    # Create merged directory name
    merged_name = f"{dir1_path.name}_merged_{dir2_path.name}"
    merged_path = dir1_path.parent / merged_name
    
    print(f"Creating merged directory: {merged_path}")
    
    if merged_path.exists():
        response = input(f"Directory {merged_path} already exists. Overwrite? [y/N]: ")
        if response.lower() != 'y':
            print("Aborted.")
            sys.exit(0)
        shutil.rmtree(merged_path)
    
    merged_path.mkdir(parents=True)
    
    # Find and copy YAML file from dir1
    yaml_files1 = list(dir1_path.glob("*.yaml"))
    yaml_filenames = None
    if yaml_files1:
        yaml_file1 = yaml_files1[0]
        shutil.copy2(yaml_file1, merged_path / yaml_file1.name)
        print(f"Copied YAML: {yaml_file1.name}")
        
        # Extract filenames from YAML for filtering
        yaml_filenames = get_yaml_filenames(yaml_file1)
        print(f"Found {len(yaml_filenames)} data files mentioned in YAML")
    else:
        print("WARNING: No YAML file found in dir1")
    
    # Get text files from both directories (filtered by YAML)
    print("\nScanning text files...")
    txt_files1 = get_txt_files(dir1_path, yaml_filenames)
    txt_files2 = get_txt_files(dir2_path, yaml_filenames)
    
    print(f"Dir1: {len(txt_files1)} file groups")
    print(f"Dir2: {len(txt_files2)} file groups")
    
    # Merge files
    print("\nMerging files...")
    
    for base_name in txt_files1:
        if base_name not in txt_files2:
            # Only in dir1, just copy
            file1 = txt_files1[base_name]
            if isinstance(file1, dict):
                # Ranked files
                for rank, filepath in file1.items():
                    output_path = merged_path / filepath.name
                    shutil.copy2(filepath, output_path)
                    print(f"  Copied (dir1 only): {filepath.name}")
            else:
                # Regular file
                output_path = merged_path / file1.name
                shutil.copy2(file1, output_path)
                print(f"  Copied (dir1 only): {file1.name}")
        else:
            # Present in both directories, merge them
            file1 = txt_files1[base_name]
            file2 = txt_files2[base_name]
            
            if isinstance(file1, dict) and isinstance(file2, dict):
                # Ranked files - merge each rank separately
                all_ranks = set(file1.keys()) | set(file2.keys())
                for rank in sorted(all_ranks):
                    if rank in file1 and rank in file2:
                        output_name = f"{base_name}.rank{rank}.txt"
                        output_path = merged_path / output_name
                        merge_text_file(file1[rank], file2[rank], output_path, skip_steps)
                    elif rank in file1:
                        # Only in dir1
                        output_path = merged_path / file1[rank].name
                        shutil.copy2(file1[rank], output_path)
                        print(f"  Copied (dir1 only): {file1[rank].name}")
                    else:
                        # Only in dir2 - skip it (shouldn't append without dir1 data)
                        print(f"  WARNING: Rank {rank} only in dir2, skipping")
            elif not isinstance(file1, dict) and not isinstance(file2, dict):
                # Regular files - merge them
                output_name = f"{base_name}.txt"
                output_path = merged_path / output_name
                merge_text_file(file1, file2, output_path, skip_steps)
            else:
                print(f"  WARNING: Mismatched file types for {base_name}, copying from dir1")
                if isinstance(file1, dict):
                    for rank, filepath in file1.items():
                        output_path = merged_path / filepath.name
                        shutil.copy2(filepath, output_path)
                else:
                    output_path = merged_path / file1.name
                    shutil.copy2(file1, output_path)
    
    # Handle files only in dir2
    for base_name in txt_files2:
        if base_name not in txt_files1:
            print(f"  WARNING: {base_name} only in dir2, skipping")
    
    print(f"\n✓ Merge complete: {merged_path}")
    return merged_path


def main():
    parser = argparse.ArgumentParser(
        description="Merge two simulation data directories with continuous step numbering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Merge two directories, skipping first 1000 steps from dir2
  %(prog)s dir1 dir2 --skip 1000
  
  # Merge without skipping any steps
  %(prog)s dir1 dir2 --skip 0
        """
    )
    
    parser.add_argument('dir1', type=str,
                        help='First simulation directory (base data)')
    parser.add_argument('dir2', type=str,
                        help='Second simulation directory (data to append)')
    parser.add_argument('--skip', type=int, default=0,
                        help='Number of initial steps to skip from dir2 (default: 0)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Simulation Data Merger")
    print("=" * 70)
    print(f"Dir1:  {args.dir1}")
    print(f"Dir2:  {args.dir2}")
    print(f"Skip:  {args.skip} steps from dir2")
    print("=" * 70)
    
    merged_path = merge_directories(args.dir1, args.dir2, args.skip)
    
    print("\nMerged directory ready for analysis!")


if __name__ == "__main__":
    main()
