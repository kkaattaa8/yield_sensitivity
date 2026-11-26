# Yield sensitivity estimation by bootstrap linear regression

# This script contains the analyses for the following study:

# Zsolt Pinke, Tamas Acs, Bence Decsi, Peter Kalicz, Zoltan Kern, Katalin Takacs, Attila Jambor: Heat and drought exposure in the European agrarian sector: crop concentration amplifies climate sensitivity in national farming systems (submitted to Scientific Reports)

library(terra)
library(parallel)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------
cfg <- list(
  # --- Workflow & reproducibility settings ---
  # "pixel_based" = gridded spatial rasters (GeoTIFF output)
  # "country" = country-level tables (CSV output)
  mode          = "pixel_based",
  n_boot        = 5000,    # Number of bootstrap resamples per regression
  rng_seed      = 42,      # Random seed for reproducibility
  min_n_obs     = 5,       # Minimum valid observations required to fit OLS
  min_boot_frac = 0.50,    # Minimum fraction of valid bootstrap slopes required
  max_win_len   = 5,       # Maximum consecutive-month window length (1 to max_win_len)
  
  # --- Crop variables & growing-season month windows (1 = Jan, ..., 12 = Dec) ---
  # Each entry defines [start, end] month indices for candidate windows
  crop_windows = list(
    wheat     = c(start = 2, end = 8),   # Feb-Aug
    wheat_w   = c(start = 2, end = 7),   # Feb-Jul (winter wheat)
    wheat_s   = c(start = 3, end = 8),   # Mar-Aug (spring wheat)
    maize     = c(start = 4, end = 8),   # Apr-Aug
    barley    = c(start = 2, end = 7),   # Feb-Jul
    rapeseed  = c(start = 2, end = 7),   # Feb-Jul
    sunflower = c(start = 4, end = 8)    # Apr-Aug
  ),
  
  # --- Pixel-based mode inputs & outputs (used when mode == "pixel_based") ---
  raster = list(
    prec_path      = "PATH/TO/precipitation.nc",
    pet_path       = "PATH/TO/pet.nc",
    tmax_path      = "PATH/TO/tmax.nc",
    tmax90p_path   = "PATH/TO/tmax90p.nc",
    crop_dirs      = list(
      wheat   = "PATH/TO/wheat_nc4",
      wheat_w = "PATH/TO/winter_wheat_nc4",
      wheat_s = "PATH/TO/spring_wheat_nc4",
      maize   = "PATH/TO/maize_nc4"
    ),
    crop_labels    = list(
      wheat   = "wheat",
      wheat_w = "winter wheat",
      wheat_s = "spring wheat",
      maize   = "maize"
    ),
    n_diff_years   = 27,        # 1994-2020 after first differencing
    output_dir     = "output/pixel_based"
  ),
  
  # --- Country CSV mode inputs & outputs (used when mode == "country") ---
  country = list(
    raw_crop_files = list(
      wheat     = "PATH/TO/crop_wheat.csv",
      maize     = "PATH/TO/crop_maize.csv",
      barley    = "PATH/TO/crop_barley.csv",
      rapeseed  = "PATH/TO/crop_rapeseed.csv",
      sunflower = "PATH/TO/crop_sunflower.csv"
    ),
    raw_climate_files = list(
      tmax     = "PATH/TO/climate_tmax.csv",
      tmax90p  = "PATH/TO/climate_tmax90.csv",
      prec     = "PATH/TO/climate_prec.csv",
      pet      = "PATH/TO/climate_pet.csv"
    ),
    diff_dir   = "data/country_first_diff",
    output_csv = "output/country/bootstrap_regression_specific_windows.csv"
  )
)

# --- Basic configuration validation ---
validate_config <- function(cfg) {
  if (!cfg$mode %in% c("pixel_based", "country")) {
    stop("cfg$mode must be 'pixel_based' or 'country', got: ", cfg$mode)
  }
  if (cfg$n_boot < 100L) {
    stop("cfg$n_boot should be at least 100, got: ", cfg$n_boot)
  }
  if (cfg$min_n_obs < 3L) {
    stop("cfg$min_n_obs must be at least 3, got: ", cfg$min_n_obs)
  }
  if (cfg$max_win_len < 1L || cfg$max_win_len > 12L) {
    stop("cfg$max_win_len must be between 1 and 12, got: ", cfg$max_win_len)
  }
  for (cname in names(cfg$crop_windows)) {
    w <- cfg$crop_windows[[cname]]
    s <- w[["start"]]
    e <- w[["end"]]
    if (is.na(s) || is.na(e) || s < 1L || e > 12L || s > e) {
      stop("Invalid month window for crop '", cname, "': start=", s, ", end=", e)
    }
  }
  invisible(TRUE)
}

validate_config(cfg)

# ------------------------------------------------------------------------------
# 2. STATISTICAL CORE: Z-SCORE STANDARDIZATION & BOOTSTRAP OLS
# ------------------------------------------------------------------------------
# Fits OLS on z-score standardized paired vectors (y, x) and estimates
# empirical two-sided p-value and 95% percentile CI via case resampling.
run_zscore_bootstrap <- function(y, x, n_boot = 5000, min_obs = 5, min_boot_frac = 0.5) {
  na_out <- c(beta = NA_real_, r2 = NA_real_, pval = NA_real_,
              ci_lower = NA_real_, ci_upper = NA_real_)
  
  valid <- complete.cases(y, x)
  if (sum(valid) < min_obs) return(na_out)
  
  y_sub <- as.numeric(y[valid])
  x_sub <- as.numeric(x[valid])
  
  sd_y <- sd(y_sub)
  sd_x <- sd(x_sub)
  if (is.na(sd_y) || sd_y == 0 || is.na(sd_x) || sd_x == 0) return(na_out)
  
  # Z-score standardization across valid observations
  y_z <- (y_sub - mean(y_sub)) / sd_y
  x_z <- (x_sub - mean(x_sub)) / sd_x
  
  # OLS fit on standardized variables
  fit_orig <- tryCatch(lm.fit(x = cbind(1, x_z), y = y_z), error = function(e) NULL)
  if (is.null(fit_orig) || length(fit_orig$coefficients) < 2 || is.na(fit_orig$coefficients[2])) {
    return(na_out)
  }
  
  beta_orig <- unname(fit_orig$coefficients[2])
  ss_res    <- sum(fit_orig$residuals^2)
  ss_tot    <- sum((y_z - mean(y_z))^2)
  r2_orig   <- if (!is.na(ss_tot) && ss_tot > 0) 1 - (ss_res / ss_tot) else NA_real_
  
  # Non-parametric bootstrap resampling of paired observations
  n_obs <- length(y_z)
  boot_betas <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n_obs, size = n_obs, replace = TRUE)
    xb  <- x_z[idx]
    yb  <- y_z[idx]
    
    if (sd(xb) == 0) {
      boot_betas[b] <- NA_real_
      next
    }
    
    b_fit <- lm.fit(x = cbind(1, xb), y = yb)
    boot_betas[b] <- unname(b_fit$coefficients[2])
  }
  
  boot_betas <- boot_betas[!is.na(boot_betas)]
  if (length(boot_betas) < (n_boot * min_boot_frac)) {
    return(c(beta = beta_orig, r2 = r2_orig, pval = NA_real_,
             ci_lower = NA_real_, ci_upper = NA_real_))
  }
  
  p_val    <- 2 * min(mean(boot_betas >= 0), mean(boot_betas <= 0))
  ci_bound <- quantile(boot_betas, probs = c(0.025, 0.975), na.rm = TRUE)
  
  c(
    beta     = beta_orig,
    r2       = r2_orig,
    pval     = p_val,
    ci_lower = unname(ci_bound[1]),
    ci_upper = unname(ci_bound[2])
  )
}

# Wrapper for terra::app() expecting a vector of length 2 * n_years
calc_reg_boot_pixel_based <- function(x, n_years = 27, n_boot = 5000,
                                 min_obs = 5, min_boot_frac = 0.5) {
  y_raw <- x[seq_len(n_years)]
  x_raw <- x[n_years + seq_len(n_years)]
  run_zscore_bootstrap(
    y = y_raw, x = x_raw,
    n_boot = n_boot, min_obs = min_obs, min_boot_frac = min_boot_frac
  )
}

# Aggregate monthly climate raster layers into multi-year seasonal means
get_aggregated_climate_raster <- function(climate_rast, start_m, end_m, n_years = 27) {
  yearly_list <- vector("list", n_years)
  for (y in seq_len(n_years)) {
    offset <- (y - 1) * 12
    m_idx  <- (offset + start_m):(offset + end_m)
    yearly_list[[y]] <- if (length(m_idx) == 1) {
      climate_rast[[m_idx]]
    } else {
      mean(climate_rast[[m_idx]])
    }
  }
  rast(yearly_list)
}

# ------------------------------------------------------------------------------
# 3. PIXEL-BASED WORKFLOW (TERRA)
# ------------------------------------------------------------------------------
run_pixel_based_workflow <- function(cfg) {
  set.seed(cfg$rng_seed)
  p_cfg <- cfg$pixel_based
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  prec         <- rast(p_cfg$prec_path)
  pet_daily    <- rast(p_cfg$pet_daily_path)
  tmax         <- rast(p_cfg$tmax_path)
  tmax90p_full <- rast(p_cfg$tmax90p_path)
  waterdef     <- pet - prec
  
  # First differences: 12-month lag for climate, 1-year lag for annual crops
  n_m <- nlyr(tmax)
  climate_vars <- list(
    tmax     = tmax[[13:n_m]]     - tmax[[seq_len(n_m - 12)]],
    tmax90p  = tmax90p[[13:n_m]]  - tmax90p[[seq_len(n_m - 12)]],
    waterdef = waterdef[[13:n_m]] - waterdef[[seq_len(n_m - 12)]]
  )
  
  crop_vars <- list()
  for (cname in names(p_cfg$crop_files)) {
    r <- rast(p_cfg$crop_files[[cname]])
    n_y <- nlyr(r)
    crop_vars[[cname]] <- r[[2:n_y]] - r[[seq_len(n_y - 1)]]
  }
  
  num_cores <- max(1, detectCores() - 1)
  message("Running pixel-based bootstrap regressions on ", num_cores, " cores...")
  
  for (cname in names(crop_vars)) {
    bounds <- windows[[cname]]
    if (is.null(bounds)) next
    
    win_start <- bounds["start"]
    win_end   <- bounds["end"]
    curr_crop <- crop_vars[[cname]]
    
    beta_layers <- list(); r2_layers <- list(); pval_layers <- list()
    ci_lo_layers <- list(); ci_hi_layers <- list()
    layer_names <- character()
    
    for (clim_name in names(climate_vars)) {
      curr_clim <- climate_vars[[clim_name]]
      for (win_len in seq_len(cfg$max_win_len)) {
        max_start <- win_end - win_len + 1
        if (win_start > max_start) next
        
        for (m_start in win_start:max_start) {
          m_end <- m_start + win_len - 1
          win_str <- if (m_start == m_end) {
            month.abb[m_start]
            } else {
              paste0(month.abb[m_start], "-", month.abb[m_end])
          lbl <- paste0(clim_name, "_", win_str)
          
          aggr_clim <- get_aggregated_climate_raster(curr_clim, m_start, m_end, n_years = p_cfg$n_diff_years)
          res_rast  <- app(
            c(curr_crop, aggr_clim),
            fun = calc_reg_boot_pixel_based,
            n_years = p_cfg$n_diff_years,
            n_boot = cfg$n_boot,
            min_obs = cfg$min_n_obs,
            min_boot_frac = cfg$min_boot_frac,
            cores = num_cores
          )
          
          beta_layers[[length(beta_layers) + 1]]   <- res_rast[[1]]
          r2_layers[[length(r2_layers) + 1]]       <- res_rast[[2]]
          pval_layers[[length(pval_layers) + 1]]   <- res_rast[[3]]
          ci_lo_layers[[length(ci_lo_layers) + 1]] <- res_rast[[4]]
          ci_hi_layers[[length(ci_hi_layers) + 1]] <- res_rast[[5]]
          layer_names <- c(layer_names, lbl)
        }
      }
    }
    
    out_stacks <- list(
      beta     = rast(beta_layers),
      r2       = rast(r2_layers),
      pval     = rast(pval_layers),
      ci_lower = rast(ci_lo_layers),
      ci_upper = rast(ci_hi_layers)
    )
    
    for (metric in names(out_stacks)) {
      names(out_stacks[[metric]]) <- layer_names
      out_path <- file.path(cfg$output_dir, paste0(cname, "_", metric, ".tif"))
      writeRaster(out_stacks[[metric]], filename = out_path, overwrite = TRUE)
    }
    message("Saved GeoTIFF outputs for crop: ", cname)
  }
  }
}

# ------------------------------------------------------------------------------
# 4. COUNTRY CSV WORKFLOW
# ------------------------------------------------------------------------------
run_country_workflow <- function(cfg) {
  c_cfg <- cfg$country
  dir.create(cfg$diff_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(cfg$output_csv), showWarnings = FALSE, recursive = TRUE)
  
  climate_tables <- c_cfg$climate_files
  if (!is.null(climate_tables$pet) && !is.null(climate_tables$prec) &&
      file.exists(climate_tables$pet) && file.exists(climate_tables$prec)) {
    df_pet  <- read.csv(climate_tables$pet, check.names = FALSE)
    df_prec <- read.csv(climate_tables$prec, check.names = FALSE)
    common_c <- intersect(df_pet$Country, df_prec$Country)
    pet_sub  <- df_pet[match(common_c, df_pet$Country), , drop = FALSE]
    prec_sub <- df_prec[match(common_c, df_prec$Country), , drop = FALSE]
    wdef_df  <- data.frame(
      Country = common_c,
      as.matrix(pet_sub[, -1]) - as.matrix(prec_sub[, -1]),
      check.names = FALSE
    )
    colnames(wdef_df)[-1] <- colnames(pet_sub)[-1]
    wdef_path <- file.path(c_cfg$diff_dir, "waterdef_input.csv")
    write.csv(wdef_df, wdef_path, row.names = FALSE)
    climate_tables$waterdef <- wdef_path
  }
  climate_tables$prec <- NULL
  climate_tables$pet  <- NULL
  
  diff_crop_files <- list()
  for (crop_code in names(c_cfg$crop_files)) {
    fpath <- c_cfg$crop_files[[crop_code]]
    if (!file.exists(fpath)) {
      warning("Missing crop input file: ", fpath)
      next
    }
    df <- read.csv(fpath, sep = ";", check.names = FALSE)
    mat <- as.matrix(df[, -1])
    diff_mat <- t(apply(mat, 1, diff))
    out_df <- data.frame(Country = df$Country, diff_mat, check.names = FALSE)
    colnames(out_df)[-1] <- paste0(colnames(mat)[-1], "_diff")
    out_path <- file.path(c_cfg$diff_dir, paste0(crop_code, "_first_diff.csv"))
    write.csv(out_df, out_path, row.names = FALSE)
    diff_crop_files[[crop_code]] <- out_path
  }
  
  diff_clim_files <- list()
  for (clim_code in names(climate_tables)) {
    cpath <- climate_tables[[clim_code]]
    if (!file.exists(cpath)) {
      warning("Missing climate input file: ", cpath)
      next
    }
    df_c <- read.csv(cpath, check.names = FALSE)
    mat_c <- as.matrix(df_c[, -1])
    diff_mat_c <- t(apply(mat_c, 1, function(row) diff(row, lag = 12)))
    out_df_c <- data.frame(Country = df_c$Country, diff_mat_c, check.names = FALSE)
    colnames(out_df_c)[-1] <- paste0(colnames(mat_c)[13:ncol(mat_c)], "_diff")
    out_path_c <- file.path(cfg$diff_dir, paste0(clim_code, "_first_diff.csv"))
    write.csv(out_df_c, out_path_c, row.names = FALSE)
    diff_clim_files[[clim_code]] <- out_path_c
  }
  
  tasks <- list()
  for (cname in names(diff_crop_files)) {
    bounds <- windows[[cname]]
    if (is.null(bounds)) next
    df_crop <- read.csv(diff_crop_files[[cname]], check.names = FALSE)
    
    for (clim_name in names(diff_clim_files)) {
      df_clim <- read.csv(diff_clim_files[[clim_name]], check.names = FALSE)
      common_countries <- intersect(df_crop$Country, df_clim$Country)
      
      for (win_len in 1:5) {
        max_start <- bounds["end"] - win_len + 1
        if (bounds["start"] > max_start) next
        
        for (m_start in bounds["start"]:max_start) {
          m_end <- m_start + win_len - 1
          win_str <- if (m_start == m_end) month.abb[m_start] else paste0(month.abb[m_start], "-", month.abb[m_end])
          lbl <- paste0(clim_name, "_", win_str)
          
          for (cntry in common_countries) {
            tasks[[length(tasks) + 1]] <- list(
              crop_name     = cname,
              clim_name     = clim_name,
              cntry         = cntry,
              target_months = m_start:m_end,
              window_label  = lbl,
              crop_row      = df_crop[df_crop$Country == cntry, , drop = FALSE],
              clim_row      = df_clim[df_clim$Country == cntry, , drop = FALSE]
            )
          }
        }
      }
    }
  }
  
  num_cores <- max(1, detectCores() - 1)
  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterSetRNGStream(cl, iseed = seed)
  clusterExport(cl, varlist = c("run_zscore_bootstrap", "n_boot", "min_obs", "min_boot_frac"),
                envir = environment())
  
  res_list <- parLapply(cl, tasks, function(task) {
    crop_cols <- colnames(task$crop_row)[-1]
    clim_cols <- colnames(task$clim_row)[-1]
    y_vals    <- as.numeric(task$crop_row[1, -1])
    x_vals    <- as.numeric(task$clim_row[1, -1])
    
    crop_years   <- sub(".*([0-9]{4}).*", "\\1", crop_cols)
    clim_years   <- sub(".*([0-9]{4}).*", "\\1", clim_cols)
    clim_months  <- rep(1:12, length.out = length(clim_cols))
    unique_years <- unique(clim_years)
    
    annual_x <- numeric()
    annual_y <- numeric()
    
    for (yr in unique_years) {
      m_idx <- which(clim_years == yr & clim_months %in% task$target_months)
      y_idx <- which(crop_years == yr)
      if (length(m_idx) > 0 && length(y_idx) > 0) {
        vx <- mean(x_vals[m_idx], na.rm = TRUE)
        vy <- y_vals[y_idx[1]]
        if (!is.na(vx) && !is.na(vy)) {
          annual_x <- c(annual_x, vx)
          annual_y <- c(annual_y, vy)
        }
      }
    }
    
    stats_out <- run_zscore_bootstrap(
      y = annual_y, x = annual_x,
      n_boot = cfg$n_boot,
      min_obs = cfg$min_n_obs,
      min_boot_frac = cfg$min_boot_frac
    )
    if (is.na(stats_out["beta"])) return(NULL)
    
    data.frame(
      Country    = task$cntry,
      Crop       = task$crop_name,
      TimeWindow = task$window_label,
      Beta       = stats_out["beta"],
      R2         = stats_out["r2"],
      pval       = stats_out["pval"],
      CI_lower   = stats_out["ci_lower"],
      CI_upper   = stats_out["ci_upper"],
      stringsAsFactors = FALSE
    )
  })
  
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) > 0) {
    final_df <- do.call(rbind, res_list)
    write.csv(final_df, cfg$output_csv, row.names = FALSE)
    message("Saved country bootstrap results to: ", cfg$output_csv)
  } else {
    warning("No valid country regression results produced.")
  }
}

# ------------------------------------------------------------------------------
# 5. EXECUTE SELECTED MODE
# ------------------------------------------------------------------------------
if (identical(cfg$mode, "pixel_based")) {
    run_pixel_based_workflow(cfg)
} else if (identical(cfg$mode, "country")) {
    run_country_workflow(cfg)
}