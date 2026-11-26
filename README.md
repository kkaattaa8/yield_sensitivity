# First-Difference Z-Score Bootstrap Regression (Pixel-Based & Country-Level)

Reproducible R workflow accompanying the manuscript:
> **Heat and drought exposure in the European agrarian sector: crop concentration amplifies climate sensitivity in national farming systems**  
>  Zsolt Pinke, Tamás Ács, Bence Decsi, Péter Kalicz, Zoltán Kern, Katalin Takács, Attila Jámbor (2026). *submitted to Scientific Reports*

## Overview
This repository provides a single, modular R script (`yield_sensitivity_bootstrap_regression.R`) to evaluate relationships between climate indicators and crop yield anomalies using:
1. **First-difference transformation** (annual lag for crop yields; 12-month lag for monthly climate variables),
2. **Z-score standardization** on valid paired observations, and
3. **Non-parametric case-resampling bootstrap OLS regression** (5,000 resamples by default) to estimate standardized slopes ($\beta$), $R^2$, empirical two-sided $p$-values, and 95% percentile confidence intervals across consecutive growing-season month windows.

The script supports two execution modes controlled via `cfg$mode`:
- **`pixel_based`**: Gridded spatial analysis using `terra`, exporting multi-band GeoTIFF stacks (`beta`, `r2`, `pval`, `ci_lower`, `ci_upper`) per crop.
- **`country`**: Country-level tabular analysis using parallelized workers (`parallel`), exporting a combined CSV table.

## Requirements
- **R** ($\ge$ 4.1.0 recommended)
- **Packages:** `terra`, `parallel` (base R)

Install `terra` if not already available:
```r
install.packages("terra")# yield_sensitivity
 Yield sensitivity estimation by bootstrap linear regression
