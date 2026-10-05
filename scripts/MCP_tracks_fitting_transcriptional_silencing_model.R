### Libraries required ----------------------------------------------------------
library(readxl) #to read excels
library(writexl) #to write excels with results
library(dplyr) #to handle tables
library(tidyr) #to handle tables
library(minpack.lm) #for linear model
library(ggplot2) #for plots ###may not be needed for publication
library(ggprism) #to make plots beautiful ###may not be needed for publication
library(rstatix) #for dunn's test ###may not be needed for publication
library(numDeriv) #for correlation analysis ###may not be needed for publication

### Load data -------------------------------------------------------------------
#Data is organized as follows:
#one tab per nucleus
#Column A: Time.to.NEBD..min., with the time
#Column B: Whole_sub, subtracted fluorescence intensity
#Column C: Norm, fluorescence normalized on the 5 values around the maximum
#Tabs are named as: Day_embryo number_nucleus letter
#The following script created a list in which every item contain the data about one nucleus, and the embryo to which it belongs

excel_file <- "path_to_data/data_file.xlsx"
sheet_names <- excel_sheets(excel_file)

nuclei_data <- lapply(sheet_names, function(s) {
  df <- read_excel(excel_file, sheet = s)
  colnames(df)[1:3] <- c("Time.to.NEBD..min.","Whole_sub", "Norm")
  df <- df[, c("Time.to.NEBD..min.","Whole_sub", "Norm")]
  df <- df[complete.cases(df), ]
  df$nucleus <- s
  df$embryo  <- sub("_[^_]+$", "", s)
  df
})
names(nuclei_data) <- sheet_names

n_nuclei <- length(nuclei_data)
cat("Loaded", n_nuclei, "nuclei.\n")

embryo_ids     <- sapply(nuclei_data, function(df) df$embryo[1])
unique_embryos <- unique(embryo_ids)
cat("Found", length(unique_embryos), "embryos:", paste(unique_embryos, collapse = ", "), "\n")


### Define model and residual function ------------------------------------------
# Fluorescence model
# All arguments named explicitly to avoid positional ambiguity.
fluorescence_model_2 <- function(t, kf, tm, tau, ts, kappa) {
  switch_term <- ifelse(t <= ts, 1, exp(-kappa * (t-ts)))
  
  switch_term * (1 + 1 / (tau * kf) * log(
    (1 + exp(kf * (t - tau - tm)))/
      (1 + exp(kf* (t - tm)))
  ))
}

# Residual function for a single nucleus.
nucleus_residuals <- function(params, df) {
  F_hat <- fluorescence_model_2(
    t   = df$Time.to.NEBD..min.,
    kf  = params["kf"],
    tm  = params["tm"],
    tau = params["tau"],
    ts = params["ts"],
    kappa = params["kappa"]
    )
  df$Norm - F_hat
}

###Multi start definition
fit_nucleus_multistart <- function(df, lower_n, upper_n, best_start_informed,
                                   n_starts = 50, seed = NULL,
                                   sse_tol = 0.05) {
  
  if (!is.null(seed)) set.seed(seed)
  
  param_names <- names(lower_n)
  
  starts <- matrix(NA_real_, nrow = n_starts, ncol = length(param_names),
                   dimnames = list(NULL, param_names))
  
  # Clamp the informed start into [lower_n, upper_n] before using it as the
  # first row. This matters because bounds (e.g. ts's upper bound) are now
  # nucleus-specific: a short trace can have ts_upper_n below the fixed
  # informed value (ts = 0.5), which would otherwise make nls.lm error out
  # on start #1 with "par out of bounds". A tiny epsilon keeps it strictly
  # inside rather than exactly on the boundary.
  eps <- 1e-6
  informed_clamped <- pmin(pmax(best_start_informed[param_names],
                                lower_n[param_names] + eps),
                           upper_n[param_names] - eps)
  starts[1, ] <- informed_clamped #the first attempt is with the given starting parameters, specified later
  
  #The other starts: the values are drawn from the range of the parameters specified in the loop
  for (p in param_names) {
    starts[-1, p] <- runif(n_starts - 1, min = lower_n[p], max = upper_n[p])
  }
  
  fits <- vector("list", n_starts)
  sse  <- rep(Inf, n_starts)
  
  for (i in seq_len(n_starts)) {
    par_i <- starts[i, ]
    
    fit_i <- tryCatch(
      nls.lm(
        par     = par_i,
        fn      = nucleus_residuals,
        df      = df,
        lower   = lower_n,
        upper   = upper_n,
        control = nls.lm.control(maxiter = 500, ftol = 1e-8, ptol = 1e-8)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_i) && fit_i$info %in% c(1, 2, 3)) {
      fits[[i]] <- fit_i
      sse[i]    <- sum(fit_i$fvec^2)
    }
  }
  
  converged_idx <- which(is.finite(sse))
  n_converged   <- length(converged_idx)
  
  if (n_converged == 0) {
    return(list(best_fit = NULL, best_sse = NA, n_converged = 0, n_starts = n_starts,
                sse_spread = NA, param_range_all = NA, param_range_near_best = NA,
                n_near_best = 0))
  }
  
  best_idx <- converged_idx[which.min(sse[converged_idx])]
  best_sse <- sse[best_idx]
  
  # Parameter matrix across ALL converged fits
  param_matrix <- t(sapply(fits[converged_idx], coef))
  
  # Range (max-min) per parameter across all converged fits
  param_range_all <- apply(param_matrix, 2, function(x) diff(range(x)))
  
  # Same, but restricted to fits "near" the best SSE (within sse_tol, e.g. 5%)
  # This isolates: among solutions that fit roughly equally well, how much
  # do the parameters themselves disagree?
  near_best_idx <- converged_idx[sse[converged_idx] <= best_sse * (1 + sse_tol)]
  n_near_best   <- length(near_best_idx)
  
  if (n_near_best >= 2) {
    param_matrix_near <- t(sapply(fits[near_best_idx], coef))
    param_range_near_best <- apply(param_matrix_near, 2, function(x) diff(range(x)))
  } else {
    param_range_near_best <- setNames(rep(NA, length(param_names)), param_names)
  }
  
  list(
    best_fit              = fits[[best_idx]],
    best_sse              = best_sse,
    n_converged           = n_converged,
    n_starts              = n_starts,
    sse_spread            = diff(range(sse[converged_idx])),
    param_range_all       = param_range_all,
    param_range_near_best = param_range_near_best,
    n_near_best           = n_near_best
  )
}

### Correlation analysis
get_param_correlation <- function(fit, df, resid_fn) {
  params <- coef(fit)
  p <- length(params)
  n <- nrow(df)
  
  resid <- resid_fn(params, df)
  J <- jacobian(function(par) resid_fn(par, df), x = params)
  
  sigma2 <- sum(resid^2) / (n - p)
  JTJ <- t(J) %*% J
  
  covar <- tryCatch(sigma2 * solve(JTJ), error = function(e) NULL)
  if (is.null(covar)) {
    cat("  JTJ singular/near-singular — strong identifiability red flag.\n")
    return(NULL)
  }
  
  colnames(covar) <- rownames(covar) <- names(params)
  cov2cor(covar)
}

corr_list <- vector("list", n_nuclei)
names(corr_list) <- sapply(nuclei_data, function(d) d$nucleus[1])


### Fit each nucleus independently given the best starting parameters just determined----------------------------------------------
output_plot_dir <- "your_path/output_directory"
dir.create(output_plot_dir, showWarnings = FALSE)

results_list <- vector("list", n_nuclei)

for (n in seq_len(n_nuclei)) {

  df          <- nuclei_data[[n]]
  nucleus_name <- df$nucleus[1]
  embryo_name  <- df$embryo[1]

  cat(sprintf("\n--- Fitting nucleus %d/%d: %s ---\n", n, n_nuclei, nucleus_name))

  #Per-nucleus bounds
  #Data-driven bound for ts
  nonzero_times_n <- df$Time.to.NEBD..min.[df$Norm != 0]
  if (length(nonzero_times_n) > 0) {
    ts_upper_n <- max(nonzero_times_n, na.rm = TRUE)
  } else {
    ts_upper_n <- 3  # fallback to fixed bound if a nucleus has no non-zero points
  }
  
  #Theoretical bounds determined using imaging conditions (Figure 1 and Table S1)
  #lower_n <- c(kf = 0.3, tau = 0.5, tm = -4, ts = -4, kappa = 0.01)
  #upper_n <- c(kf = 9, tau = 15, tm = 2, ts = 2, kappa = 9.2)
  
  #Experimentally determined bounds for kf, tau and tm, theoretical for ts and kappa (Figure 4 and Table S2)
  #lower_n <- c(kf = 0.81, tau = 2.43, tm = -2.92, ts = -1, kappa = 0.01) 
  #upper_n <- c(kf = 2.61, tau = 3.53, tm = -1.96, ts = ts_upper_n, kappa = 9.2)
  
  #Experimentally determined bounds for kf, tau, tm and ts, theoretical for kappa; for Lds- and Luciferase nuclei (Figure 5 and Table S3)
  lower_n <- c(kf = 0.81, tau = 2.43, tm = -2.92, ts = -0.07, kappa = 0.01) 
  upper_n <- c(kf = 2.61, tau = 3.53, tm = -1.96, ts = 2, kappa = 9.2) 
  
  #Array of values to start the fitting, chosen as the median of the distributions of each parameter
  informed_start <- c(kf = 1.71, tau = 2.97, tm = -2.44, ts = 0.5, kappa = 4.6)
  
  #Fitting
  ms_result <- fit_nucleus_multistart(
    df, lower_n, upper_n, informed_start,
    n_starts = 50, #n of starts per nucleus
    seed = n   # makes the restart reproducible per nucleus
  )
  
  #Shows the number of restarts that converged for each nucleus
  cat(sprintf("  %d/%d starts converged, SSE spread = %.4f\n",
              ms_result$n_converged, ms_result$n_starts, ms_result$sse_spread))
  
  if (is.null(ms_result$best_fit)) {
    cat(sprintf("  WARNING: no successful fit for %s across %d starts\n",
                nucleus_name, ms_result$n_starts))
    next
  }
  
  fit_n <- ms_result$best_fit

  if (is.null(fit_n)) next

  params_n <- coef(fit_n)
  
  #Multistart: parameter ranges and sse
  best_sse <- sum(fit_n$fvec^2)
  sse_spread_relative <- ms_result$sse_spread / ms_result$best_sse
  #Parameter ranges across ALL converged starts
  pr_all  <- ms_result$param_range_all
  #Parameter ranges among near-best-SSE starts only (the more diagnostic one)
  pr_near <- ms_result$param_range_near_best
  
  #Correlation matrix build
  corr_list[[nucleus_name]] <- get_param_correlation(fit_n, df, nucleus_residuals)
  
  #Convergence and bounds checks
  conv_n <- switch(as.character(fit_n$info),
    "1" = "converged (SSE)",
    "2" = "converged (params)",
    "3" = "converged (both)",
    "4" = "max iterations reached",
    "5" = "WARNING: Jacobian rank-deficient",
    "other"
  )

  at_upper <- names(params_n)[params_n >= upper_n[names(params_n)] - 1e-6]
  at_lower <- names(params_n)[params_n <= lower_n[names(params_n)] + 1e-6]
  bounds_warning <- ""
  if (length(at_upper) > 0)
    bounds_warning <- paste("AT UPPER BOUND:", paste(at_upper, collapse = ", "))
  if (length(at_lower) > 0)
    bounds_warning <- paste(bounds_warning, "AT LOWER BOUND:", paste(at_lower, collapse = ", "))
  bounds_flag <- ifelse(bounds_warning == "", "none", trimws(bounds_warning))

  #R² for the nucleus
  F_hat <- fluorescence_model_2(
    t   = df$Time.to.NEBD..min.,
    kf  = params_n["kf"],
    tm  = params_n["tm"],
    tau = params_n["tau"],
    ts = params_n["ts"],
    kappa = params_n["kappa"]
  )
  ss_r <- sum((df$Norm - F_hat)^2)
  ss_t <- sum((df$Norm - mean(df$Norm))^2)
  r2_n <- 1 - ss_r / ss_t 
  
  #RMSE in the decay region only
  #The decay region is defined as timepoints between tm and tm + tau, where the model predicts the signal should be actively falling
  # A large RMSE here indicates the model cannot follow the observed decay shape, even if the overall R² looks acceptable.
  decay_idx  <- df$Time.to.NEBD..min. >= params_n["tm"] & 
    df$Time.to.NEBD..min. <= params_n["tm"] + params_n["tau"]
  
  n_decay    <- sum(decay_idx)
  
  if (n_decay >= 2) {
    rmse_decay <- sqrt(mean((df$Norm[decay_idx] - F_hat[decay_idx])^2))
  } else {
    rmse_decay <- NA
    cat("  WARNING: fewer than 2 points in decay window — RMSE decay set to NA.\n")
  }
  
  #Overall RMSE for comparison
  rmse_total <- sqrt(mean((df$Norm - F_hat)^2))
  
  #Print results
  cat(sprintf(" kf=%.3f  tau=%.3f  tm=%.3f  R²=%.3f  RMSE_decay=%.3f  [%s]  %s\n",
              params_n["kf"], params_n["tau"], params_n["tm"],
              r2_n, ifelse(is.na(rmse_decay), NaN, rmse_decay), conv_n, bounds_flag))
  
  # --- Residuals ---------------------------------------------------------------
  df$residual <- df$Norm - F_hat
  df$in_decay <- decay_idx   # flag decay region points for coloring in plot
  
  # --- Time of 90 and 10% fluorescence -------------------------------------------------------------------
  t_dense   <- seq(min(df$Time.to.NEBD..min.), max(df$Time.to.NEBD..min.), by = 0.001)
  F_dense   <- fluorescence_model_2(t_dense, kf= params_n["kf"], tm= params_n["tm"], tau= params_n["tau"],
                                    ts = params_n["ts"], kappa = params_n["kappa"]) #pay attention to the model number and parameters
  F_max <- max(F_dense)
  
  F90 <- (F_max/100)*90 ###10% drop
  F10 <- (F_max/100)*10 ###90% drop
  
  t90 <- t_dense[which.min(abs(F_dense-F90))] #it does not work with an absolute number, do not know why
  t10 <- t_dense[which.min(abs(F_dense-F10))]
  

  # --- Store results ----------------------------------------------------------
  results_list[[n]] <- data.frame(
    nucleus         = nucleus_name,
    embryo          = embryo_name,
    kf              = params_n["kf"],
    tau             = params_n["tau"],
    tm              = params_n["tm"],
    elongation_rate_515 = 5.15 / params_n["tau"],
    elongation_rate_58 = 5.8 / params_n["tau"],
    Time_90_F = t90,
    Time_10_F = t10,
    ts         = params_n["ts"],
    kappa            = params_n["kappa"],
    rmse_total          = round(rmse_total, 4),
    rmse_decay          = round(rmse_decay, 4),  # key metric for model failure
    n_decay_points      = n_decay,
    r2              = round(r2_n, 4),
    n_points        = nrow(df),
    convergence     = conv_n,
    bounds_warning  = bounds_flag,
    row.names       = NULL,
    sse_spread_relative = round(sse_spread_relative, 4),
    n_starts_converged = ms_result$n_converged,
    sse_spread_relative  = round(sse_spread_relative, 4),
    n_near_best          = ms_result$n_near_best,
    kf_range_all         = round(pr_all["kf"], 4), #spread across all run, even some that were terrible
    tau_range_all        = round(pr_all["tau"], 4),
    tm_range_all         = round(pr_all["tm"], 4),
    ts_range_all         = round(pr_all["ts"], 4),
    kappa_range_all      = round(pr_all["kappa"], 4),
    kf_range_near_best    = round(pr_near["kf"], 4), #spead across only good runs
    tau_range_near_best   = round(pr_near["tau"], 4),
    tm_range_near_best    = round(pr_near["tm"], 4),
    ts_range_near_best    = round(pr_near["ts"], 4),
    kappa_range_near_best = round(pr_near["kappa"], 4)
  )

  # --- Plot -------------------------------------------------------------------
  df$Fit <- F_hat

  g <- ggplot(df, aes(x = Time.to.NEBD..min., y = Norm)) +
    geom_point(size = 4, alpha = 0.5) +
    geom_line(aes(y = Fit), linewidth = 1.5, color = "dodgerblue3") +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("R² = %.3f\nB = %.3f\nkf = %.3f\ntau = %.3f\ntm = %.3f",
                             r2_n, params_n["B"], params_n["kf"], params_n["tau"], params_n["tm"]),
             hjust = -0.1, vjust = 1.3, size = 5, color = "dodgerblue3") +
    scale_x_continuous(breaks = seq(floor(t_min), ceiling(t_max), by = 1)) +
    geom_hline(yintercept = F90)+
    geom_hline(yintercept = F10)+
    xlab("Time to NEBD (min)") +
    ylab("Subtracted Fluorescence intensity") +
    ggtitle(nucleus_name) +
    theme_prism(base_size = 25)

  ggsave(g,
         filename = paste0(nucleus_name, "_per_nucleus_fit.png"),
         path     = output_plot_dir,
         width    = 8, height = 6)
}


### Summarise results -----------------------------------------------------------
results_df <- bind_rows(results_list)
write_xlsx(results_df, "your_path/results_file.xlsx")


### Results evaluation ----------------------------------------------------------

#Evaluates the spread of the SSE across restarts over the best SSE obtained per nucleus: the smaller the spread, the smaller the more all re-runs agree and end up in similar places
#see_spread_relative is the spead/best
ggplot(results_df, aes(x = sse_spread_relative)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = quantile(results_df$sse_spread_relative, 0.9, na.rm = TRUE),
             linetype = "dashed", color = "firebrick") +
  #scale_x_continuous(limits=c(-0.5, 10), breaks = seq(0, 10, by=1))+
  xlab("SSE spread / best SSE (per nucleus)") +
  ylab("Number of nuclei") +
  ggtitle("Multistart robustness across nuclei") +
  theme_prism(base_size = 20)

#Correlation matrix
pairs <- combn(c("kf", "tau", "tm", "ts", "kappa"), 2, simplify = FALSE)
summary_df <- do.call(rbind, lapply(pairs, function(pp) {
  vals <- sapply(corr_list, function(m) if (!is.null(m)) m[pp[1], pp[2]] else NA)
  data.frame(
    pair             = paste(pp, collapse = "-"),
    median_abs_corr  = median(abs(vals), na.rm = TRUE),
    pct_nuclei_gt_09 = mean(abs(vals) > 0.9, na.rm = TRUE),
    n_valid          = sum(!is.na(vals))
  )
}))
summary_df
