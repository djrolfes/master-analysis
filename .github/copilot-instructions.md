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

#### Correlation Function (`cf`) Objects

A `cf` (correlation function) object is a specialized R object for storing time-series data from lattice simulations, along with metadata. It uses a **mixin pattern** where you build up the object by adding layers of functionality.

**Core Concept**: The `cf` object stores data in a matrix format where:
- **Rows** = individual measurements (gauge configurations)
- **Columns** = time slices OR different observables

**When to use `cf` objects vs. raw vectors**:
- Use `cf` objects when you have structured data with time/spatial dimensions (e.g., correlation functions, Wilson flow observables at different flow times)
- Use raw vectors with `uwerrprimary()` or `bootstrap.analysis()` for scalar observables per configuration (e.g., plaquette, topological charge, single action density value)

#### Building a `cf` Object (Layer-by-Layer)

```r
# Step 1: Create empty cf object
my_cf <- cf()

# Step 2: Add original data (REQUIRED)
#   data_matrix: N_configs × N_timeslices (or N_observables)
my_cf <- cf_orig(my_cf, cf = data_matrix)

# Step 3: Add metadata (REQUIRED for most operations)
my_cf <- cf_meta(my_cf,
                 nrObs = 1,              # Number of different observables
                 Time = 48,              # Time extent (or number of columns if not time-based)
                 nrStypes = 1,           # Number of smearing types
                 symmetrised = FALSE)    # Whether data is symmetrized

# Step 4: Add bootstrap samples (REQUIRED for error analysis)
my_cf <- bootstrap.cf(my_cf,
                      boot.R = 400,      # Number of bootstrap samples
                      boot.l = 2,        # Block length
                      seed = 1234)       # RNG seed
```

#### Single Observable (Scalar per Configuration)

For observables like **topological charge** or **plaquette** (one value per configuration):

```r
# Read data
data <- read.table("topological_charge.txt", header = TRUE)
topo_charge <- data$topological_charge  # Vector of N values

# Create cf object: N rows × 1 column
cf_data <- matrix(topo_charge, ncol = 1)
topo_cf <- cf() %>%
  cf_orig(cf = cf_data) %>%
  cf_meta(nrObs = 1, Time = 1, nrStypes = 1, symmetrised = FALSE) %>%
  bootstrap.cf(boot.R = 400, boot.l = 2, seed = 1234)

# Use it
plot(topo_cf)
summary(topo_cf)
```

#### Multiple Time Slices (e.g., Wilson Flow Observable)

For observables measured at different **flow times** or **time slices**:

```r
# Data structure: each row = configuration, each column = flow time
# Example: 100 configs × 20 flow times
flow_data <- matrix(your_data, nrow = 100, ncol = 20)

flow_cf <- cf() %>%
  cf_orig(cf = flow_data) %>%
  cf_meta(nrObs = 1, Time = 20, nrStypes = 1, symmetrised = FALSE) %>%
  bootstrap.cf(boot.R = 400, boot.l = 2, seed = 1234)

# Plot will show observable vs. flow time with error bands
plot(flow_cf)
```

#### Multiple Observables (Combining Different Measurements)

To combine **multiple different observables** into one `cf` object (needed for GEVP or simultaneous analysis):

```r
# Create individual cf objects first
obs1_cf <- cf() %>% 
  cf_orig(cf = matrix(obs1_data, ncol = Time)) %>%
  cf_meta(nrObs = 1, Time = Time, nrStypes = 1, symmetrised = FALSE)

obs2_cf <- cf() %>%
  cf_orig(cf = matrix(obs2_data, ncol = Time)) %>%
  cf_meta(nrObs = 1, Time = Time, nrStypes = 1, symmetrised = FALSE)

# Concatenate BEFORE bootstrapping (important!)
combined_cf <- c(obs1_cf, obs2_cf)  # Now nrObs = 2
combined_cf <- cf_meta(combined_cf, nrObs = 2, Time = Time, nrStypes = 1, symmetrised = FALSE)
combined_cf <- bootstrap.cf(combined_cf, boot.R = 400, boot.l = 2, seed = 1234)

# Requirements for concatenation (c.cf):
# - Same Time extent
# - Same number of measurements (rows)
# - Same symmetrisation
# - Same nrStypes
```

#### Key `hadron` Functions for `cf` Objects

- **`cf()`**: Creates an empty `cf` object.

- **`cf_orig(.cf, cf, icf = NULL)`**: Adds original data to a `cf` object.
  - `.cf`: The `cf` object to extend.
  - `cf`: Numeric matrix (N_measurements × N_timeslices) of real data.
  - `icf`: Optional imaginary part (use with caution, many functions ignore it).

- **`cf_meta(.cf, nrObs, Time, nrStypes, symmetrised)`**: Adds metadata.
  - `nrObs`: Number of different observables in the object.
  - `Time`: Time extent (or number of columns for non-time data).
  - `nrStypes`: Number of smearing types (usually 1).
  - `symmetrised`: Logical, whether data has been symmetrized.

- **`bootstrap.cf(cf, boot.R, boot.l, ...)`**: Performs bootstrap resampling on a `cf` object to prepare it for statistical error analysis.
  - `cf`: The correlation function object (must have `cf_orig` mixin).
  - `boot.R`: Number of bootstrap samples to generate (default: 400).
  - `boot.l`: Block length for the resampling to handle autocorrelations (default: 2).
  - `seed`: Seed for random number generation (default: 1234).
  - `sim`: Simulation type - "geom" (default) or "fixed" for block bootstrap.
  - Returns: A `cf` object with bootstrap samples added (adds `cf_boot` mixin).

- **`c.cf(...)`**: Concatenates multiple `cf` objects into one (increases `nrObs`).
  - `...`: Multiple `cf` objects with compatible dimensions.
  - Returns: A single `cf` object containing all observables.

- **`plot.cf(x, ...)`**: Plots the correlation function with error bands.
  - `x`: A `cf` object with bootstrap samples (`cf_boot` mixin required).
  - Returns: Invisibly returns a data frame with `t`, `CF`, and `Err` columns.

#### Non-`cf` Functions for Simple Time Series

For scalar observables (one value per configuration) that don't need the `cf` structure:

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

#### Practical Guidelines for Analysis Scripts

**Use `cf` objects when**:
- You have Wilson flow observables at multiple flow times
- You have correlation functions with time structure
- You need to combine multiple observables (e.g., for GEVP)
- You want to leverage `hadron`'s built-in plotting for time-dependent data

**Use raw vectors with `uwerrprimary()`/`bootstrap.analysis()` when**:
- You have scalar observables (plaquette, single action density, acceptance rate)
- You have topological charge (one value per configuration)
- You want detailed autocorrelation analysis
- You're analyzing thermalization or simple time series

**Example: Rewriting `analysis_topological_charge.R` with `cf` objects**:
```r
# Current approach (raw vector)
topo_data <- read_data_file(topo_path)$topological_charge
uw_result <- uwerrprimary(topo_data, pl = TRUE)

# Alternative with cf object
cf_data <- matrix(topo_data, ncol = 1)
topo_cf <- cf() %>%
  cf_orig(cf = cf_data) %>%
  cf_meta(nrObs = 1, Time = 1, nrStypes = 1, symmetrised = FALSE) %>%
  bootstrap.cf(boot.R = 400, boot.l = 2, seed = 1234)
plot(topo_cf)  # Automatic plotting with error bands
```

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