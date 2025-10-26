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
A `cf` (correlation function) object is a specialized R object for storing time-series data from lattice simulations, along with metadata. The `hadron` package provides several functions for analyzing these objects and other time-series data.

- **`bootstrap.cf(cf, boot.R, boot.l, ...)`**: Performs bootstrap resampling on a `cf` object to prepare it for statistical error analysis.
  - `cf`: The correlation function object.
  - `boot.R`: Number of bootstrap samples to generate (default: 400).
  - `boot.l`: Block length for the resampling to handle autocorrelations (default: 2).
  - `seed`: Seed for random number generation (default: 1234).
  - `sim`: Simulation type - "geom" (default) or "fixed" for block bootstrap.
  - Returns: A `cf` object with bootstrap samples added.

- **`bootstrap.analysis(data, boot.R, boot.l, skip, ...)`**: Conducts a bootstrap analysis on a generic time series (a numeric vector) to determine the mean, error, and integrated autocorrelation time.
  - `data`: A numeric vector containing the time series.
  - `boot.R`: Number of bootstrap samples (default: 100).
  - `boot.l`: Block length for resampling (default: 2).
  - `skip`: Number of initial measurements to discard for thermalization (default: 0).
  - `pl`: Logical, whether to plot the result (default: FALSE).
  - Returns: A data frame with mean, error, error of error, tau int, and bias for all block sizes.

- **`bootstrap.meanerror(data, R, l)`**: Returns the bootstrap error (standard deviation) of the mean for a time series.
  - `data`: The input time series (numeric vector).
  - `R`: Number of bootstrap replicates (default: 400).
  - `l`: Block length (default: 20).
  - **Returns**: A numeric value (or vector) with the standard error only (NOT a named vector with mean and error).

- **`uwerrprimary(data, nrep, S, pl)`**: Analyzes a primary time series using the Gamma method to determine the mean, error, and integrated autocorrelation time (`tauint`). This is an alias for `uwerr` without a function argument.
  - `data`: A numeric vector containing the time series.
  - `nrep`: Vector of replica lengths (default: all data in one replica).
  - `S`: Initial guess for tau/tauint ratio (default: 1.5).
  - `pl`: Logical, whether to plot autocorrelation analysis (default: FALSE).
  - Returns: A `uwerr` object with `value`, `dvalue`, `tauint`, `dtauint`, etc.

- **`uwerrderived(f, data, nrep, S, pl, ...)`**: Calculates a secondary observable and its error from one or more primary time series, correctly propagating statistical errors using the Gamma method.
  - `f`: A user-defined function that computes the secondary observable from a single measurement of the primary observables.
  - `data`: A matrix where each column is a primary observable and each row is a measurement.
  - `nrep`: Vector of replica lengths (default: all data in one replica).
  - `S`: Initial guess for tau/tauint ratio (default: 1.5).
  - `pl`: Logical, whether to plot (default: FALSE).
  - `...`: Additional arguments passed to function `f`.
  - Returns: A `uwerr` object containing the derived observable's statistics.

- **`computeacf(tseries, W.max, Lambda)`**: Computes the autocorrelation function (ACF) and integrated autocorrelation time for a time series.
  - `tseries`: The input time series data.
  - `W.max`: The maximum time lag to use for the calculation.
  - `Lambda`: Cut-off for estimating standard error of ACF (default: 100).
  - Returns: A `hadronacf` object with `tau`, `dtau`, `Gamma`, `dGamma`, etc.

- **Installation**: Use `scripts/check_and_install_hadron.sh` for dependencies.
- **Documentation**: Generate with `R CMD Rd2pdf .` in `extern/hadron/man/`.

### SLURM Workflow
`scripts/run_job.sh` creates unique output dirs, copies input.yaml, runs simulation, then submits dependent analysis job. Always work with output directories containing both `input.yaml` and generated data files.

## Analysis Script Patterns
1. Parse YAML config with `read_yaml_config()`
2. Map YAML params to filenames via `GaugeObservableParams` 
3. Handle thermalization by skipping initial steps
4. Use hadron's bootstrap routines for error analysis
5. Generate histograms for `exp(-delta_H)` acceptance analysis