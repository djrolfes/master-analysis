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
The `hadron` R package provides specialized tools for analyzing lattice QCD time series data, with emphasis on bootstrap error estimation and autocorrelation analysis.

#### Core Concepts

**Correlation Function Objects (`cf`)**
- `cf` objects are S3 class containers designed for correlation functions from QFT/statistical simulations
- Structure: List-based with metadata (`Time`, `nrObs`, `nrStypes`, `symmetrised`)
- Data stored as arrays: `cf$cf` has dimension `c(N, nrObs*nrStypes*(Time/2+1))` where N = number of measurements
- Bootstrap samples added via `cf_boot` mixin: `cf$cf.tsboot$t` contains bootstrap replicas
- Arithmetic operations defined: `+`, `-`, `*`, `/` work element-wise
- Related: `raw_cf` for complex/matrix-valued correlators before averaging

#### Key Functions Used in Analysis

**1. Bootstrap Error Estimation**

`bootstrap.meanerror(data, R = 400, l = 20)`
- **Purpose**: Compute bootstrap error of the mean for 1D time series
- **Usage**: Quick error estimates for single observables
- **Parameters**:
  - `data`: Numeric vector of observations
  - `R`: Number of bootstrap samples (default 400)
  - `l`: Block length for block bootstrap (default 20)
- **Returns**: Single numeric value (standard error of mean)
- **Current use**: Wilson flow observables, PTBC analysis
- **Example**:
```r
error <- hadron::bootstrap.meanerror(topological_charge, R = 200)
mean_with_error <- c(mean(topological_charge), error)
```

**2. Full Bootstrap Analysis**

`bootstrap.analysis(data, skip = 0, boot.R = 100, tsboot.sim = "geom", pl = FALSE, boot.l = 2)`
- **Purpose**: Comprehensive bootstrap analysis with blocking for autocorrelation
- **Usage**: Main tool for analyzing primary observables with autocorrelation
- **Parameters**:
  - `data`: Numeric vector (time series)
  - `skip`: Number of initial measurements to skip (thermalization)
  - `boot.R`: Number of bootstrap samples
  - `tsboot.sim`: Bootstrap type - `"geom"` (geometric block) or `"fixed"` (fixed block)
  - `pl`: If TRUE, generates diagnostic plots (history, blocking analysis)
  - `boot.l`: Block length for bootstrap
- **Returns**: Data frame with columns:
  - `mean`: Mean value of observable
  - `error`: Error estimate
  - `error.error`: Error of the error (uncertainty in error estimate)
  - `tau_int`: Integrated autocorrelation time
  - `bias`: Bootstrap bias estimate
  - (One row per block size tested)
- **Algorithm**: Tests increasing block sizes from 1 until `(length(data)-skip)/blocksize < 20`
- **Current use**: Topological charge mean value analysis
- **Example**:
```r
pdf("bootstrap_output.pdf")
boot_result <- hadron::bootstrap.analysis(data$topo, skip = 100, pl = TRUE)
dev.off()
# boot_result is a data frame with analysis for different block sizes
```

**3. Autocorrelation Analysis (Gamma Method)**

`uwerrprimary(data, pl = FALSE)` / `uwerr(f, data, nrep, S = 1.5, pl = FALSE, ...)`
- **Purpose**: Autocorrelation analysis via Gamma method (Wolff, hep-lat/0306017)
- **Usage**: Gold standard for computing autocorrelation times and properly accounting for correlations
- **Parameters**:
  - `data`: For primary: numeric vector; For derived: matrix (N × Nalpha)
  - `nrep`: Vector or single integer giving replica length(s). For single replica: `nrep = length(data)`
  - `S`: Initial guess for tau/tauint ratio (default 1.5)
  - `pl`: If TRUE, plots ACF, integrated autocorr time, and time history
  - `f`: Function for derived quantities (takes data vector, returns numeric/vector)
- **Returns**: Object of class `uwerr` with:
  - `$value`: Expectation value
  - `$dvalue`: Error estimate
  - `$ddvalue`: Error on the error
  - `$tauint`: Integrated autocorrelation time
  - `$dtauint`: Error of tauint
  - `$Qval`: p-value (for multiple replicas)
  - `$Gamma`: Normalized autocorrelation function
  - `$Wopt`: Optimal cut-off for Gamma integration
- **Advantages**: More sophisticated than simple blocking; provides tauint directly
- **Current use**: Topological charge autocorrelation time estimation, Wilson loop analysis
- **Example (primary observable)**:
```r
# Single replica analysis
uw <- hadron::uwerrprimary(topological_charge, pl = TRUE)
tauint <- uw$tauint
error <- uw$dvalue
plot(uw)  # Generates diagnostic plots
```
- **Example (derived observable)**:
```r
# For derived quantities f(x) with error propagation
uw <- hadron::uwerr(f = function(x) mean(x^2), data = data, nrep = length(data))
```
- **Example (with proper nrep for single replica)**:
```r
# When analyzing W(L,T) for each (L,T) combination
W_clean <- W_temp[!is.na(W_temp)]
uwerr_result <- uwerr(data = W_clean, nrep = length(W_clean), S = 1.5, pl = FALSE)
tau_int <- uwerr_result$tauint  # Integrated autocorrelation time
effective_n <- length(W_clean) / (2 * tau_int)  # Effective sample size
```

**4. Simple ACF Computation**

`computeacf(tseries, W.max, Lambda = 100)`
- **Purpose**: Compute autocorrelation function and integrated autocorr time
- **Usage**: Direct ACF calculation with error estimates
- **Parameters**:
  - `tseries`: Time series vector
  - `W.max`: Maximum time lag (typically `length(data)/5`)
  - `Lambda`: Cut-off for ACF standard error estimation
- **Returns**: Object of class `hadronacf` with:
  - `$Gamma`: Normalized ACF
  - `$dGamma`: Error of ACF
  - `$tau`: Integrated autocorrelation time
  - `$dtau`: Error of tau
  - `$lags`: Time lags
  - `$W`: Cut-off used for integration
- **Algorithm**: Uses Madras-Sokal formula for tau error (via hep-lat/0409106)
- **Current use**: Fallback ACF analysis for topological charge
- **Example**:
```r
acf_result <- hadron::computeacf(data, floor(length(data)/5))
plot(acf_result)
tau <- acf_result$tau
```

**5. Correlation Function Bootstrap**

`**`bootstrap.cf(cf, boot.R = 400, boot.l = 2, seed = 1234, sim = "geom", endcorr = TRUE)`**
- **Purpose**: Bootstrap entire correlation function objects (for multi-time analysis)
- **Usage**: When analyzing correlation functions C(t) across multiple time slices
- **CRITICAL**: `bootstrap.cf()` requires a `cf` object (correlation function), NOT a numeric vector
- **Parameters**:
  - `cf`: Correlation function object (class `cf` with subclass `cf_orig`)
  - `boot.R`: Number of bootstrap samples
  - `boot.l`: Block length
  - `seed`: RNG seed for reproducibility
  - `sim`: `"geom"` (geometric) or `"fixed"` (fixed block)
  - `endcorr`: End corrections for fixed blocks
- **Returns**: Enhanced `cf` object with:
  - `$cf.tsboot`: Bootstrap samples (from `tsboot()`)
  - `$cf0`: Original average
  - `$tsboot.se`: Bootstrap standard errors
  - `$boot.R`, `$boot.l`, `$seed`: Parameters used
- **Common Mistake**: DO NOT use `bootstrap.cf()` on scalar observables (numeric vectors)
  - ❌ WRONG: `boot_E <- bootstrap.cf(E_data, boot.R=1500)` where `E_data` is numeric vector
  - ✅ CORRECT: Use `uwerr()` for scalar observables, which provides autocorrelation-corrected errors directly
  - For error propagation: If `y = f(x)` and `uw <- uwerr(x)`, then `error_y = |df/dx| * uw$dvalue`
- **Example (correct usage for correlation functions)**:
```r
# For time-dependent correlators C(t)
cf_obj <- create_cf_object(correlator_matrix)  # Must be cf object
cf_boot <- bootstrap.cf(cf_obj, boot.R = 1500, boot.l = 1)
# Access bootstrap samples: cf_boot$cf.tsboot$t
```
- **Example (scalar observable error propagation)**:
```r
# For scalar observable like action density
E_data <- data$action_density
uw_E <- uwerrprimary(E_data, pl=FALSE)
mean_E <- uw_E$value
error_E <- uw_E$dvalue  # Already accounts for autocorrelation!

# Error propagation for t0^2 * E
t0_squared_E <- t0^2 * mean_E
error_t0_squared_E <- t0^2 * error_E  # Simple propagation
````
- **Purpose**: Bootstrap entire correlation function objects (for multi-time analysis)
- **Usage**: When analyzing correlation functions C(t) across multiple time slices
- **Parameters**:
  - `cf`: Correlation function object (class `cf`)
  - `boot.R`: Number of bootstrap samples
  - `boot.l`: Block length
  - `seed`: RNG seed for reproducibility
  - `sim`: `"geom"` (geometric) or `"fixed"` (fixed block)
  - `endcorr`: End corrections for fixed blocks
- **Returns**: Enhanced `cf` object with:
  - `$cf.tsboot`: Bootstrap samples (from `tsboot()`)
  - `$cf0`: Original average
  - `$tsboot.se`: Bootstrap standard errors
  - `$boot.R`, `$boot.l`, `$seed`: Parameters used
- **Current use**: Action density analysis (time-dependent energy measurements)
- **Example**:
```r
cf_obj <- bootstrap.cf(E_data, boot.R = 1500, boot.l = 1)
# Access bootstrap samples: cf_obj$cf.tsboot$t
# Original mean: cf_obj$cf0
```

**6. Non-Linear Least-Squares Fitting**

Hadron provides several NLS fitting functions with different bootstrap strategies:

**`parametric.nlsfit(fn, par.guess, boot.R, y, dy, x, dx, ...)`**
- **Purpose**: Parametric bootstrap NLS fit (recommended for most use cases)
- **Usage**: Fit theoretical models to data with Gaussian error propagation
- **Parameters**:
  - `fn`: Fit function `fn(par, x, boot.r, ...)` where `boot.r` is 0 for original data
  - `par.guess`: Initial parameter values (vector)
  - `boot.R`: Number of bootstrap samples (typical: 200-1500)
  - `y`, `dy`: Data and errors on dependent variable
  - `x`, `dx`: Independent variable and errors (dx optional)
  - `lower`, `upper`: Parameter bounds (vectors, default: ±Inf)
  - `bootstrap`: If FALSE, uses Jacobian for errors instead of bootstrap
- **Returns**: Object of class `bootstrapfit` with:
  - `$t0`: Fitted parameters on original data + chisqr (last element)
  - `$t`: Bootstrap samples (matrix: boot.R × (npar+1))
  - `$se`: Bootstrap standard errors for parameters
  - `$Qval`: p-value of fit quality
  - `$chisqr`, `$dof`: Chi-squared and degrees of freedom
- **Current use**: Static potential V(L) fitting in `analysis_wilson_temp.R`
- **Example**:
```r
# Exponential fit: W(T) = a * exp(-V*T)
fn_exp <- function(par, x, boot.r, ...) {
  par[1] * exp(-par[2] * x)
}

fit_result <- parametric.nlsfit(
  fn = fn_exp,
  par.guess = c(a_init, V_init),
  y = W_data,
  dy = W_errors,
  x = T_values,
  boot.R = 200,
  lower = c(0, 0)  # Force positive parameters
)

# Extract results
V_fit <- fit_result$t0[2]
V_error <- fit_result$se[2]
chi2_per_dof <- fit_result$chisqr / fit_result$dof
```

**`bootstrap.nlsfit(fn, par.guess, y, x, bsamples, ...)`**
- **Purpose**: Full bootstrap NLS fit with pre-generated bootstrap samples
- **Usage**: Maximum flexibility, handles x-y correlations
- **Parameters**:
  - `bsamples`: Pre-generated bootstrap samples (array: boot.R × n)
  - `priors`: Constrain parameters from previous fits
  - `CovMatrix`: Full covariance matrix for correlated fits
  - `use.minpack.lm`: Use minpack.lm for faster convergence (default TRUE)
- **When to use**: Correlated data, x-y error correlation, Bayesian priors
- **Example**:
```r
# Generate bootstrap samples first
bsamples <- parametric.bootstrap(boot.R = 500, c(y, x), c(dy, dx))

fit <- bootstrap.nlsfit(
  fn = fn,
  par.guess = c(1, 1),
  y = y,
  x = x,
  bsamples = bsamples,
  CovMatrix = NULL  # Will be computed from bsamples
)
```

**`simple.nlsfit(fn, par.guess, y, x, errormodel, ...)`**
- **Purpose**: Quick fits without bootstrap (faster)
- **Usage**: When bootstrap errors not needed, or for initial parameter estimation
- **Parameters**:
  - `errormodel`: "yerrors" or "xyerrors"
  - `boot.R`: If > 0, generates parametric bootstrap after fit (optional)
- **Returns**: Same structure as `bootstrap.nlsfit`

**`matrixfit(cf, t1, t2, parlist, sym.vec, neg.vec, useCov = FALSE, model = "single", ...)`**
- **Purpose**: Factorizing matrix fits for correlator matrices (GEVP results)
- **Usage**: Fit multiple correlators simultaneously with shared parameters
- **Models**: `"single"`, `"shifted"`, `"pc"` (principal correlator)
- **Returns**: Fitted parameters for each matrix element

**`fit.effectivemass(cf, t1, t2, useCov = FALSE, ...)`**
- **Purpose**: Fit constant to effective mass plateau
- **Usage**: Extract masses from correlation functions
- **Requires**: `cf` must be `effectivemass` object from `bootstrap.effectivemass()`

#### Common Analysis Workflow Patterns

**Pattern 1: Simple Observable Error (No Autocorrelation Needed)**
```r
# For quick error estimates (Wilson flow, distances, etc.)
mean_val <- mean(observable)
error_val <- hadron::bootstrap.meanerror(observable, R = 200)
```

**Pattern 2: Full Analysis with Autocorrelation**
```r
# For critical observables (topological charge, plaquette)
pdf("analysis_output.pdf")

# Bootstrap blocking analysis
boot_result <- hadron::bootstrap.analysis(data, skip = therm_steps, 
                                          boot.R = 500, pl = TRUE)
dev.off()

# Autocorrelation time via Gamma method
uw <- hadron::uwerrprimary(data[(therm_steps+1):length(data)], pl = TRUE)
tauint <- uw$tauint
error <- uw$dvalue

# Fallback: Direct ACF
acf_obj <- hadron::computeacf(data, floor(length(data)/5))
```

**Pattern 3: Correlation Function Analysis**
```r
# For time-dependent observables
cf_data <- create_cf_object(time_series_data)  # Custom function
cf_boot <- bootstrap.cf(cf_data, boot.R = 1500, boot.l = 1)

# Extract effective mass
effmass <- bootstrap.effectivemass(cf_boot, type = "solve")
effmass_fit <- fit.effectivemass(effmass, t1 = 5, t2 = 15)
plot(effmass_fit)
```

#### Plotting with Hadron

Hadron provides specialized plotting functions that handle error bars and bootstrap fits elegantly:

**`plot(fit_object, ...)`** - S3 method for fitted objects
- Works with `bootstrapfit`, `uwerr`, `matrixfit`, `effectivemass` objects
- Automatically plots data points, fit curve, and error bands
- **Key parameters**:
  - `plot.range`: Vector `c(xmin, xmax)` to extend plot beyond data range
  - `col.line`: Fit line color (default: "black")
  - `col.band`: Error band color (default: "gray")
  - `opacity.band`: Error band transparency (default: 0.65)
  - `lwd`: Line width for fit curve
  - `rep`: If TRUE, adds to existing plot (like gnuplot's replot)
- **Example**:
```r
# Plot NLS fit with extended x-range
plot(potential_fit,
     xlab = "Distance L",
     ylab = "Potential V(L)",
     main = "Static Potential Fit",
     plot.range = c(0, 10),
     col.line = "black",
     col.band = "gray",
     opacity.band = 0.5,
     lwd = 2)

# Add vertical line with error band
abline(v = r_0, col = "blue", lwd = 1)
rect(r_0 - r_0_error, par("usr")[3], 
     r_0 + r_0_error, par("usr")[4],
     col = rgb(0, 0, 1, alpha = 0.15), border = NA)
```

**`plotwitherror(x, y, dy, dx, ...)`** - General error bar plotting
- **Purpose**: Create scatter plots with error bars (y-errors, x-errors, or both)
- **Parameters**:
  - `x`, `y`: Coordinate vectors
  - `dy`, `dx`: Error vectors (can be matrix for multiple error types)
  - `mdx`, `mdy`: Non-symmetric errors (negative direction)
  - `errsum.method`: How to sum multiple errors ("linear", "quadrature", "linear.quadrature")
  - `rep`: If TRUE, adds to existing plot
- **Example**:
```r
plotwitherror(L_values, V_values, V_errors,
              xlab = "L", ylab = "V(L)",
              main = "Potential vs Distance",
              col = "red", pch = 19)
```

**Hadron Plotting Style Guidelines**

The hadron package uses a consistent, professional plotting style suitable for publications:

**Visual Design**:
- **Frame**: Complete box around plot (`bty = "o"` or default, showing all four sides)
- **Colors**: Black lines (`col.line = "black"`), gray error bands (`col.band = "gray"`)
- **Transparency**: Error bands use `opacity.band = 0.65` (65% opacity)
- **Line width**: `lwd = 2` for fit curves, standard for data
- **Error bars**: Include end caps (perpendicular lines at bar endpoints)
- **Point markers**: Standard filled circles (`pch = 19`) or as appropriate

**Implementation**:
```r
# Base R plotting with hadron style
par(bty = "o")  # Box type: complete frame around plot (default)
plotwitherror(x, y, dy, 
              xlab = "Variable [units]",
              ylab = "Observable [units]",
              pch = 19, col = "black")

# Add fitted curve with error band
pcol <- rgb(128/255, 128/255, 128/255, alpha = 0.65)
polygon(c(x_fit, rev(x_fit)), 
        c(y_fit + err_fit, rev(y_fit - err_fit)),
        col = pcol, border = NA)
lines(x_fit, y_fit, col = "black", lwd = 2)
```

**Best Practices**:
- Use `pdf()` / `dev.off()` to save plots instead of screen display in scripts
- Set reasonable `xlim` and `ylim` to show error bands clearly
- Use `plot.range` to extend fits beyond data (e.g., to show scale settings)
- Combine hadron plots with base R graphics (`abline`, `rect`, `text`, `legend`)
- Keep frame visible (`bty = "o"`) for professional appearance
- For publications: `lwd = 2`, semi-transparent error bands (`opacity.band = 0.3-0.5`)

#### Bootstrap Parameter Guidelines

**Number of Bootstrap Samples (`boot.R`)**
- **Quick tests**: 99-200 samples (fast, reasonable errors)
- **Standard analysis**: 400-500 samples (good balance of speed/accuracy)
- **Publication quality**: 1000-1500 samples (reliable error estimates)
- **Heavy computations** (GEVP, matrix fits): 200-400 (computational cost)
- **Rule of thumb**: Use more samples when:
  - Data has strong autocorrelations (high τ_int)
  - Fit has many parameters (need stable covariance matrix)
  - Results are sensitive to outliers

**Block Length (`boot.l`)**
- **Purpose**: Account for autocorrelation in time series bootstrap
- **Short correlations** (τ_int < 5): `boot.l = 1-2`
- **Medium correlations** (τ_int ~ 10-20): `boot.l = 2-5`
- **Long correlations** (τ_int > 20): `boot.l = 5-10` or use `uwerr` instead
- **Geometric blocking** (`tsboot.sim = "geom"`): More robust than fixed blocks
- **Rule of thumb**: 
  - Set `boot.l ≈ τ_int / 2` for fixed blocks
  - Use `boot.l = 1-2` with geometric blocking (default)
  - Always check that `(N - skip) / boot.l > 20` (sufficient independent blocks)

**Common Parameter Combinations**
```r
# Fast exploratory analysis
bootstrap.analysis(data, skip = 100, boot.R = 200, boot.l = 1, pl = TRUE)

# Standard production analysis  
bootstrap.cf(cf, boot.R = 400, boot.l = 2, sim = "geom")

# High-precision final results
parametric.nlsfit(fn, par.guess, boot.R = 1500, y, dy, x)

# Heavy correlation functions (action density, GEVP)
bootstrap.cf(cf, boot.R = 1500, boot.l = 1, sim = "geom")
```

**Validation Checks**
- Compare `bootstrap.analysis` with `uwerr` for same data (should agree)
- Check that error plateaus with increasing `boot.R` (convergence)
- Verify `(N - skip) / boot.l > 20` to ensure sufficient independent blocks
- Monitor `τ_int`: if `τ_int > boot.l`, increase `boot.l` or use `uwerr`

#### Installation & Documentation
- **Installation**: Use `scripts/check_and_install_hadron.sh` for dependencies
- **Manual pages**: Located in `extern/hadron/man/`
- **Documentation**: Generate PDF with `R CMD Rd2pdf extern/hadron`
- **Key references**:
  - Gamma method: Wolff, Comput.Phys.Commun. 156 (2004) 143, hep-lat/0306017
  - Covariance matrix inversion: Michael & McKerrell, Phys.Rev. D51 (1995) 3745, hep-lat/9412087

#### Hadron Functions by Analysis Script

**`analysis_topological_charge.R`**
- `bootstrap.analysis()`: Full bootstrap analysis of topological charge with plots
- `uwerrprimary()`: Gamma method autocorrelation analysis for tauint estimation
- `computeacf()`: Fallback ACF computation for tau_int
- Applied to: raw topological charge, rounded charge, Q², rounded Q²

**`analysis_Plaquette.R`**
- `bootstrap.analysis()`: Bootstrap analysis of plaquette values with blocking

**`analysis_SimulationLoggingParams_log_file.R`**
- `bootstrap.analysis()`: Analyze exp(-delta_H) and acceptance rate
- `uwerrprimary()`: Autocorrelation analysis of delta_H

**`analysis_ptbc_log.R`**
- `bootstrap.meanerror()`: Quick error estimates for PTBC log data per rank
- Applied to: acceptance rates per rank

**`analysis_action_density.R`**
- `uwerrprimary()`: Compute mean energy density with autocorrelation correction
- `bootstrap.cf()`: Bootstrap correlation function for action density time series

**`analysis_wilsonflow_tests.R`**
- `bootstrap.meanerror()`: Error estimates for flow observables
- Applied to: flow time data, distance to integer, action density, Wilson W(t)

**`analysis_wilsonflow_improv.R`**
- `bootstrap.meanerror()`: Compare normal vs improved Wilson flow action densities

**`analysis_wilsonflow_details.R`**
- Uses hadron for observables analysis (library loaded but specific functions vary by context)

**`analysis_wilson_temp.R`**
- Uses hadron for static potential analysis (library loaded for future NLS fits)

### SLURM Workflow
`scripts/run_job.sh` creates unique output dirs, copies input.yaml, runs simulation, then submits dependent analysis job. Always work with output directories containing both `input.yaml` and generated data files.

## Analysis Script Patterns
1. Parse YAML config with `read_yaml_config()`
2. Map YAML params to filenames via `GaugeObservableParams` 
3. Handle thermalization by skipping initial steps
4. Use hadron's bootstrap routines for error analysis
5. Generate histograms for `exp(-delta_H)` acceptance analysis