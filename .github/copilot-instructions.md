# Master Analysis - AI Agent Instructions

## Project Overview
Lattice QCD analysis pipeline combining `klft` (C++/Kokkos HMC simulation) with Python/R analysis. The workflow: generate HMC configurations → analyze observables → bootstrap statistics.

## Key Architecture
- **`extern/`**: Git submodules (`klft` simulation engine, `hadron` R package for analysis)
- **`inputs/`**: YAML configs defining HMC parameters, observables, and output filenames
- **`analysis/`**: Cross-language IO (Python/R) and analysis scripts
- **`scripts/`**: SLURM job submission orchestrating simulation → analysis pipeline

## YAML Configuration Pattern
All analysis flows from `input.yaml` structure:
```yaml
HMCParams: { Ndims: 4, L0: 8, Nc: 2, seed: 1234, coldStart: false }
GaugeObservableParams: { 
  plaquette_filename: "plaquette_output.txt",
  action_density_filename: "action_density.txt",
  measurement_interval: 2
}
SimulationLoggingParams: { 
  log_filename: "simulation_log.txt",
  log_delta_H: true, log_acceptance: true
}
```

## Data File Conventions
- **Headers**: First line contains column names (not comments)
- **HMC Step**: Always first column, rename to `hmc_step` in data_io functions
- **Format**: Whitespace-separated, comments start with `#`
- **Key observables**: `delta_H`, `acceptance`, `plaquette`, `action_density`

## Critical Development Patterns

### Cross-Language IO
Use `data_io.py` and `data_io.R` with matching function signatures:
```python
# Python pattern
def read_data_plaquette_filename(config, base_dir):
    filename = config['GaugeObservableParams']['plaquette_filename']
    return read_data_file(os.path.join(base_dir, filename))
```
```r
# R equivalent  
read_data_plaquette_filename <- function(config, base_dir) {
  filename <- config$GaugeObservableParams$plaquette_filename
  read_data_file(file.path(base_dir, filename))
}
```

### Hadron Package Integration
- **Bootstrap functions**: `bootstrap.cf()`, `bootstrap.analysis()`, `bootstrap.meanerror()`  
- **Installation**: Use `scripts/check_and_install_hadron.sh` for dependencies
- **Documentation**: Generate with `R CMD Rd2pdf .` in `extern/hadron/man/`

### SLURM Workflow
`scripts/run_job.sh` creates unique output dirs, copies input.yaml, runs simulation, then submits dependent analysis job. Always work with output directories containing both `input.yaml` and generated data files.

## Analysis Script Patterns
1. Parse YAML config with `read_yaml_config()`
2. Map YAML params to filenames via `GaugeObservableParams` 
3. Handle thermalization by skipping initial steps
4. Use hadron's bootstrap routines for error analysis
5. Generate histograms for `exp(-delta_H)` acceptance analysis