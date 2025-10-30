# Cluster "Illegal Operation" Error - Fix Instructions

## Problem
The R analysis script crashes on the cluster with:
```
*** caught illegal operation ***
address 0x7f47420b5393, cause 'illegal operand'
```

This occurs in the QR decomposition step of `nls()` during static potential fitting, typically caused by numerical instability (very small errors, NaN values, or extreme weight ranges).

## Solution

### 1. Updated Files to Deploy
Make sure these updated files are on the cluster:

- **`analysis/analysis_wilson_temp.R`** - Contains enhanced numerical safeguards:
  - Pre-checks for non-finite values
  - Error floor with multiple safety levels
  - Error normalization to prevent overflow/underflow
  - Weight range validation (rejects if range > 1e6)
  - Fallback to unweighted fitting if weights are problematic
  - Enhanced logging showing error and weight ranges

- **`scripts/check_and_install_hadron.sh`** - Updated to install additional R packages:
  - Added `minpack.lm` (optional but recommended)
  - Added `ggplot2` (required for plotting)
  - Added `yaml` (required for config reading)

### 2. Deploy Steps on Cluster

```bash
# 1. Navigate to your project directory
cd /path/to/master-analysis

# 2. Pull latest changes (if using git)
git pull origin main

# 3. Re-run the dependency installation script
./scripts/check_and_install_hadron.sh

# 4. Verify the updated analysis script is in place
grep "error_scale" analysis/analysis_wilson_temp.R
# Should return a line with "error_scale <- median(safe_errors)"
```

### 3. What Changed in the Numerical Safeguards

#### Before (OLD CODE - causes crash):
```r
weights <- 1 / (error^2)
```

#### After (NEW CODE - safe):
```r
# Multiple safety layers:
min_error <- min(L_data$error[L_data$error > 0], na.rm = TRUE)
error_floor <- max(1e-8, min_error * 0.01, 1e-10)
safe_errors <- pmax(L_data$error, error_floor)

# Normalize to prevent overflow/underflow
error_scale <- median(safe_errors)
normalized_errors <- safe_errors / error_scale

weights <- 1 / (normalized_errors^2)

# Validate weights
if (any(!is.finite(weights)) || max(weights) / min(weights) > 1e6) {
    weights <- rep(1, nrow(L_data))  # Fall back to unweighted
}
```

### 4. Understanding the Logging

The updated script logs error and weight ranges for debugging:

```
fit_static_potential: L=1 error range: [5.4e-05, 0.00012], weight range: [0.39, 1.78]
```

This helps diagnose:
- **Error range**: Shows if errors are extremely small (< 1e-8 triggers extra safeguards)
- **Weight range**: Shows if weights span too wide a range (> 1e6 triggers unweighted fit)

### 5. Test on Cluster

After deploying, test with a small dataset:

```bash
cd analysis
Rscript analysis_wilson_temp.R /path/to/test/data 0
```

Check the log file for:
- ✅ "error range" and "weight range" lines (confirms new code is running)
- ✅ "fit results" lines (confirms fits are completing)
- ❌ Any "ERROR fitting" or "non-finite" messages (indicates data issues)

### 6. If Still Crashes

If the crash persists even with updated code:

1. **Check the log file** to see which L value crashes:
   ```bash
   tail -100 wilson_temp_analysis.log
   ```

2. **Inspect the data** for that L value - look for:
   - NaN or Inf values
   - Extremely small errors (< 1e-10)
   - Identical error values (suggests zero variance)

3. **Try unweighted fit** by forcing this in the code:
   ```r
   weights <- rep(1, nrow(L_data))  # Force unweighted
   ```

4. **Check R version and libraries**:
   ```bash
   R --version
   R -e "sessionInfo()"
   ```
   - Older R versions may have different numerical stability
   - Check if `BLAS/LAPACK` libraries are optimized versions (Intel MKL, OpenBLAS)

## Technical Details

### Why This Happens
The QR decomposition in `nls()` operates on the matrix `.swts * gr` where:
- `.swts` = square root of weights
- `gr` = gradient matrix

When errors are extremely small (< 1e-8) or vary by many orders of magnitude:
- Weights become extremely large (> 1e16)
- Matrix elements can overflow or create denormalized numbers
- CPU raises "illegal operand" signal when encountering these

### Our Solutions
1. **Error floor**: Prevents division by tiny numbers
2. **Normalization**: Scales errors to ~1.0 median, preventing extreme weights
3. **Range check**: Detects when weights span > 1e6 (indicates problems)
4. **Unweighted fallback**: Uses equal weights when numerical issues detected

These safeguards maintain fit quality while preventing crashes.

## Success Indicators

After deploying the fix, you should see:
- ✅ Script completes without crashes
- ✅ PDF plots generated (static_potential_fit.pdf)
- ✅ Scale setting file created (scale_setting_sommer.txt)
- ✅ Log shows "error range" and "weight range" for each L value
- ✅ Reasonable fit parameters (A ~ 0.4, B ~ 0.16, sigma ~ 0.03)

## Contact
If problems persist after deploying these changes, provide:
1. Full error message
2. Last 100 lines of wilson_temp_analysis.log
3. Output of `head -20 wilson_temp_output.txt` (to check data format)
4. R version: `R --version`
