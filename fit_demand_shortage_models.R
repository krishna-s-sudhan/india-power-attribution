setwd("~/code-location") #ensure your working directory is set
# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)
predictor = 'cdd'
# ----- Specify the states in the study area -----
# study_area_states <- c()
study_area_states <- c("PJ","HR","RJ","DL","UP","OD","MP","GJ","CG","TG","JH","UK","BR") # full study area
# ----- Specify the months in the study area -----
study_period_months <- c(1,2,3,4,5,6,7,8,9,10,11,12)
# ----- Pick dependent variable ('electricity_shortage', 'electricity_requirement')
  dependent_variable <- 'electricity_requirement'

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
if (dependent_variable == 'electricity_requirement'){
  data_path  <- "data/prep_panel/electricity_panel_with_cdd.csv"
  output_dir <- "data/fit_statistics/demand_models"
} else if (dependent_variable == 'electricity_shortage'){
  data_path  <- "data/prep_panel/electricity_panel_with_cdd.csv"
  output_dir <- "data/fit_statistics/shortage_models"
} else {
  stop("invalid dependent variable")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# ---------------------------------------------------------------
# 3. Load and prepare data
# ---------------------------------------------------------------
cat("Loading data...\n")
data_panel <- read_csv(data_path, show_col_types=FALSE)
# --- Restrict to study-area states ---
if (!is.null(study_area_states)) {
  data_panel <- data_panel %>% filter(state %in% study_area_states)
}
if (!is.null(study_period_months)) {
  data_panel <- data_panel %>% filter(month %in% study_period_months)
}
# Sanity check: confirm the filter kept what you expect
cat(sprintf("Study area: %d states, %d state-months retained\n",
            dplyr::n_distinct(data_panel$state), nrow(data_panel)))
missing <- setdiff(study_area_states, unique(data_panel$state))
if (length(missing) > 0)
  cat("WARNING: requested codes not found in data:", paste(missing, collapse = ", "), "\n")
spline_df <- 4L    # default df; override per-variable below if you like
# Build an ns basis ONCE on the estimation sample, store the object (knots kept
# for reuse in prediction/attribution), add its columns, return the column names.
ns_bases <- list()                       # registry of fitted basis objects, keyed by variable
add_ns <- function(var, df = spline_df) {
  b <- ns(data_panel[[var]], df = df)
  ns_bases[[var]] <<- b                   # <<- writes to the registry in the enclosing scope
  nm <- paste0(var, "_ns", seq_len(ncol(b)))
  data_panel[nm] <<- as.data.frame(unclass(b))
  nm
}
cdd_ns   <- add_ns("cdd")
dewpt_ns <- add_ns("cdd_dewpt")
swbgt_ns <- add_ns("cdd_swbgt")
tmean_ns <- add_ns("t_mean",  df = 4L)   # temperature response is hump-shaped — keep an interior knot
rh_ns    <- add_ns("rh",      df = 3L)   # humidity covariate — fewer df is usually plenty
# keep the CDD aliases the rest of the script already uses, so nothing downstream changes
cdd_basis   <- ns_bases[["cdd"]]
cdd_knots          <- attr(cdd_basis, "knots")
cdd_boundary_knots <- attr(cdd_basis, "Boundary.knots")

# [EDIT 1] `basis_names` is used in sections 9, 10d and 11 but was never
# assigned — it is just the CDD basis column names add_ns() already returned.
# Defining it here (a character vector, used in no formula) changes no estimate.
basis_names <- cdd_ns                    # cdd_ns1 ... cdd_ns4

cat("Data prepared.\n\n")
# ---------------------------------------------------------------
# 4. Model definitions - keeping it simple for now, then improving
# ---------------------------------------------------------------

models_list <- list(
  'S-Y' = c(
    'state',
    'year'
  ),
  'Y-M' = c(
    'month',
    'year'
  ),
  'S-M' = c(
    'state',
    'month'
  ),
  'Ys-Ms'= c(
    'state_month',
    'state_year'
  ),
  'CDD' = c(
    cdd_ns
  ),
  'CDD-S' = c(
    cdd_ns,
    'state'
  ),
  'CDD-Ms' = c(
    cdd_ns,
    'state_month'
  ),
  'CDD-S-Y' = c(
    cdd_ns,
    'state',
    'year'
  ),
  'CDD-Y-M' = c(
    cdd_ns,
    'year',
    'month'
  ),
  'CDD-S-M' = c(
    cdd_ns,
    'state',
    'month'
  ),
  'CDD-Ys-Ms' = c(
    cdd_ns,
    'state_month',
    'state_year'
  ),
  'CDD-RH-Ys-Ms' = c(
    cdd_ns,
    rh_ns,
    'state_month',
    'state_year'
  ),
  'dCDD-Ys-Ms' = c(
    dewpt_ns,
    'state_month',
    'state_year'
  ),
  'wCDD-Ys-Ms' = c(
    swbgt_ns,              # was 'cdd_swbgt'
    'state_month',
    'state_year'
  ),
  'T-Ys-Ms' = c(
    tmean_ns,              # was 't_mean'
    'state_month',
    'state_year'
  )
)

# ---------------------------------------------------------------
# 5. Helper: build formula
# ---------------------------------------------------------------
build_formula <- function(vars) {
  fixed_effects <- intersect(vars, c("month",
                                     "year",
                                     'state_month',
                                     'state_year',
                                     'state'))
  predictors <- setdiff(vars, fixed_effects)
  
  if (length(predictors) > 0 && length(fixed_effects) > 0) {
    formula_str <- paste(dependent_variable, "~", paste(predictors, collapse = " + "),
                         "|", paste(fixed_effects, collapse = " + "))}
  else if (length(predictors) > 0) {
    formula_str <- paste(dependent_variable, "~", paste(predictors, collapse = " + "))}
  else if (length(fixed_effects) > 0) {
    formula_str <- paste(dependent_variable, "~ 1 |", paste(fixed_effects, collapse = " + "))}
  else {
    formula_str <- paste(dependent_variable, "~ 1")
  }
  as.formula(formula_str)
}
# ---------------------------------------------------------------
# 6. Helper: significance stars from p-value
# ---------------------------------------------------------------
sig_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE       ~ ""
  )
}
# ---------------------------------------------------------------
# 6b. Standard-error clustering
# ---------------------------------------------------------------
# Two-way cluster by state and year: lets shortage/demand be correlated
# within a state across time (persistent grid/capacity issues) AND across
# states within the same year (a shared heatwave or policy shock), instead
# of assuming independence along whichever dimension isn't clustered on.
# fixest computes this via Cameron-Gelbach-Miller (2011): V_state + V_year -
# V_state_year. Defined once here and reused in every feols() call below so
# every fit stat, coefficient table, and curve uses the same SE convention.
cluster_formula <- ~state + year
# ---------------------------------------------------------------
# 7. Fit all models and collect results
# ---------------------------------------------------------------
all_coefs      <- list()
all_fit_stats  <- list()
all_model_fit  <- list()
for (model_idx in seq_along(models_list)) {
  
  model_name <- names(models_list)[model_idx]
  model_vars <- models_list[[model_idx]]
  
  cat(sprintf("[%2d/%2d] Fitting: %s\n",
              model_idx, length(models_list), model_name))
  
  model_formula <- build_formula(model_vars)
  print(model_formula)
  
  fit_result <- tryCatch({
    fixest::feols(
      model_formula,
      data     = data_panel,
      cluster  = cluster_formula,
      nthreads = 1
    )
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    NULL
  })
  
  print(fit_result)
  
  if (is.null(fit_result)) next
  
  # --- Model-level fit statistics ---
  fit_stats <- data.frame(
    model_index      = model_idx,
    model_name       = model_name,
    aic              = AIC(fit_result),
    bic              = BIC(fit_result),
    r_squared        = r2(fit_result, type = "cor2"),
    adj_r_squared    = unname(r2(fit_result, type = "ar2")),
    within_r_squared = unname(r2(fit_result, type = "wr2")),
    n_obs            = fit_result$nobs,
    n_fe = length(fit_result$fixef_vars),   # NULL -> length 0 for no-FE models,
    stringsAsFactors = FALSE
  )
  all_fit_stats[[model_idx]] <- fit_stats
  
  # --- Coefficient table ---
  # coeftable() returns estimate, std error, t/z stat, p-value for non-FE terms
  ct <- coeftable(fit_result)
  
  if (!is.null(ct) && nrow(ct) > 0) {
    coef_df <- as.data.frame(ct)
    coef_df$term        <- rownames(ct)
    rownames(coef_df)   <- NULL
    
    # Standardise column names (fixest uses "Std. Error", "z value", "Pr(>|z|)")
    colnames(coef_df) <- make.names(colnames(coef_df))  # safe names
    # Rename to consistent names regardless of fixest version
    col_map <- c(
      "Estimate"      = "estimate",
      "Std..Error"    = "std_error",
      "z.value"       = "z_value",
      "t.value"       = "z_value",      # alias
      "Pr...z.."      = "p_value",
      "Pr...t.."      = "p_value"       # alias
    )
    for (old in names(col_map)) {
      if (old %in% colnames(coef_df)) {
        colnames(coef_df)[colnames(coef_df) == old] <- col_map[[old]]
      }
    }
    
    coef_df <- coef_df %>%
      mutate(
        model_index  = model_idx,
        model_name   = model_name,
        significance = sig_stars(p_value)
      ) %>%
      select(model_index, model_name, term,
             estimate, std_error, z_value, p_value, significance)
    
    all_coefs[[model_idx]] <- coef_df
  }
  
  # --- observed vs fitted shortage stored ---
  
  fitted_full <- tryCatch(
    predict(fit_result, sample = "original"),
    error = function(e) {
      cat("  Prediction error:", e$message, "\n")
      NULL
    }
  )
  
  if (!is.null(fitted_full)) {
    prediction_df <- data_panel %>%
      mutate(.fitted = as.numeric(fitted_full)) %>%
      filter(!is.na(.fitted)) %>%
      {
        required_columns <- c(
          "state",
          "year",
          "month",
          dependent_variable
        )
        
        if (all(required_columns %in% names(.))) . else NULL
      }
    
    if (!is.null(prediction_df) && nrow(prediction_df) > 0) {
      prediction_agg <- prediction_df %>%
        mutate(
          model_index = model_idx,
          model_name  = model_name
        ) %>%
        select(
          model_index,
          model_name,
          state,
          year,
          month,
          all_of(dependent_variable),
          .fitted
        )
      
      all_model_fit[[model_idx]] <- prediction_agg
    }
  }
  cat(sprintf("  AIC=%.1f  BIC=%.1f  R²=%.4f  n=%d\n",
              fit_stats$aic, fit_stats$bic, fit_stats$r_squared, fit_stats$n_obs))
  
}
cat("Models with stored predictions:", length(all_model_fit), "\n")
# ---------------------------------------------------------------
# 8. Plot Results
# ---------------------------------------------------------------
library(ggplot2)
# --- Temperatures with highlight on Mar-Apr-May 2022 ---
core <- data_panel %>%
  mutate(date = as.Date(sprintf("%d-%02d-01", year, month)))
p <- ggplot(core, aes(date, .data[[dependent_variable]])) +
  annotate("rect", xmin = as.Date("2022-03-01"), xmax = as.Date("2022-05-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#c0392b") +
  geom_line(linewidth = 0.4, colour = "grey30") +
  geom_point(size = 0.5, colour = "grey30") +
  facet_wrap(~ state, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = dependent_variable,
       title = "Northwest core: shortage over time (2022 heatwave shaded)") +
  theme_minimal(base_size = 10)
print(p)
# --- Observed vs fitted, by model ---
fit_all <- bind_rows(all_model_fit)
if (nrow(fit_all) == 0) {
  warning("No fitted predictions were stored.")
}
p_fit <- ggplot(
  fit_all,
  aes(
    x = .data[[dependent_variable]],
    y = .fitted
  )
) +
  geom_point(alpha = 0.3, size = 0.6, colour = "grey30") +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "#c0392b",
    linewidth = 0.5
  ) +
  facet_wrap(~ model_name) +
  coord_equal() +
  labs(
    x = paste("Observed", dependent_variable),
    y = paste("Fitted", dependent_variable),
    title = ""
  ) +
  theme_minimal(base_size = 10)
print(p_fit)
ggsave(file.path(output_dir, "obs_vs_fitted_by_model.pdf"), p_fit,
       width = 10, height = 8, dpi = 150)
# --- Observed vs fitted, Timeseries ---
fit_long <- fit_all %>%
  mutate(date = as.Date(sprintf("%d-%02d-01", year, month))) %>%
  tidyr::pivot_longer(
    cols = c(all_of(dependent_variable), .fitted),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = if_else(series == ".fitted", "Fitted", "Observed")
  )
p_ts <- ggplot(fit_long, aes(date, value, colour = series)) +
  geom_line(linewidth = 0.4) +
  facet_grid(state ~ model_name, scales = "free_y") +
  scale_colour_manual(
    values = c(Observed = "grey30", Fitted = "#c0392b")
  ) +
  labs(
    x = NULL,
    y = dependent_variable,
    colour = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(legend.position = "top")
print(p_ts)
ggsave(file.path(output_dir, "timeseries_obs_vs_fitted_by_model.pdf"), p_ts,
       width = 10, height = 8, dpi = 150)
# ---------------------------------------------------------------
# 9. Combine and save
# ---------------------------------------------------------------
coef_table        <- bind_rows(all_coefs)
fit_stats_table   <- bind_rows(all_fit_stats)
shortage_table   <- bind_rows(all_model_fit)
coef_path         <- file.path(output_dir, "all_models_coefficients.csv")
fit_stats_path    <- file.path(output_dir, "all_models_fit_stats.csv")
shortage_path    <- file.path(output_dir, "all_models_fit.csv")
write_csv(coef_table,      coef_path)
write_csv(fit_stats_table, fit_stats_path)
write_csv(shortage_table, shortage_path)
cat("\n=============================================\n")
cat("Saved coefficient table to:      ", coef_path,      "\n")
cat("Saved fit statistics table to:   ", fit_stats_path, "\n")
cat("Saved shortage table to:", shortage_path, "\n")
cat("Study Area:", study_area_states, "\n")
cat("Study Months:", study_period_months, "\n")
# --- Basis details stored ---
# (this block needs `basis_names`, defined in section 3 — see [EDIT 1])
basis_path_grid <- file.path(output_dir, "cdd_basis_grid.csv")
basis_path_long <- file.path(output_dir, "basis_vcov_long.csv")
# 1) exact R basis on a fine CDD grid (sidesteps any R/Python basis mismatch)
grid  <- seq(min(data_panel$cdd, na.rm = TRUE), max(data_panel$cdd, na.rm = TRUE), length.out = 400)
Bgrid <- predict(cdd_basis, newx = grid)
colnames(Bgrid) <- basis_names                      # cdd_ns1 … cdd_ns4
write.csv(cbind(cdd = grid, Bgrid), basis_path_grid, row.names = FALSE)
spline_models <- names(models_list)[sapply(models_list, \(v) all(basis_names %in% v))]
vcov_long <- purrr::map_dfr(spline_models, function(nm) {
  V  <- vcov(fixest::feols(build_formula(models_list[[nm]]), data = data_panel,
                           cluster = cluster_formula))[basis_names, basis_names]
  df <- as.data.frame(as.table(as.matrix(V)))
  stats::setNames(cbind(df, model_name = nm), c("term_i","term_j","cov","model_name"))
})
write.csv(vcov_long, basis_path_long, row.names = FALSE)
# ---------------------------------------------------------------
# 9. Print a quick summary to console
# ---------------------------------------------------------------
cat("\n--- Model fit statistics summary ---\n")
fit_stats_table %>%
  select(model_index, model_name, aic, bic, r_squared, adj_r_squared, within_r_squared, n_obs, n_fe) %>%
  mutate(across(c(aic, bic, r_squared, adj_r_squared, within_r_squared), ~round(., 4))) %>%
  print()
cat("\nDone!\n")
# ================================================================
# 10. Climate attribution: predict demand under ALL vs NAT scenarios
# ================================================================
# Feeds climate-model CDD (all-forcings "ALL" vs natural-only "NAT")
# through one fitted response function. For a given state-month the only
# thing differing between the two predictions is the CDD term, so
# ALL - NAT isolates the demand attributable to anthropogenic warming.
# The fixed-effect (non-climate) baseline is held at the study-period
# average so it cancels in the difference.
#
# CDD file (all_state_cdd.csv): columns year, month, state, cdd, leg, member.
#   - `leg`    = scenario (ALL / NAT)
#   - `member` = ensemble member (525-member ensemble)
#   - `state`  = FULL names -> remapped to panel codes below
#   - single year (2022); months 3,4,5 (5 dropped by study_period_months)
# ----------------------------------------------------------------
# ----- 10a. Attribution configuration -----
# Response function to use (must be a models_list entry containing `predictor`)
attribution_model_name <- "CDD-Ys-Ms"
# Climate-model CDD file
attribution_data_path <- file.path(
  "data/attribution_runs/all_state_cdd.csv"
)
# Column names in that file
scenario_col <- "leg"        # holds the ALL / NAT labels
member_col   <- "member"     # ensemble-member id
cdd_col      <- "cdd"        # the CDD column (mapped onto the model's predictor)
# How the two scenarios are labelled in scenario_col
all_label <- "ALL"   # factual / all-forcings   (orange "ACT" in the figure)
nat_label <- "NAT"   # counterfactual / natural (blue)
# Full state name (as in the CDD file) -> panel state code.
# Only study-area states need mapping; anything unmapped is dropped below.
name_to_code <- c(
  "Punjab"         = "PJ",
  "Haryana"        = "HR",
  "Rajasthan"      = "RJ",
  "NCT of Delhi"   = "DL",
  "Uttar Pradesh"  = "UP",
  "Odisha"         = "OD",
  "Madhya Pradesh" = "MP",
  "Gujarat"        = "GJ",
  "Chhattisgarh"   = "CG",
  "Telangana"      = "TG",
  "Jharkhand"      = "JH",
  "Uttarakhand"    = "UK",
  "Bihar"          = "BR"
)
# Baseline reference years (state_year FE averaged over these). Default = all fitted years.
attribution_reference_years <- sort(unique(data_panel$year))
attribution_out_path <- file.path(output_dir, "attribution_predictions.csv")
# ----- 10b. Re-fit the chosen response function -----
stopifnot(attribution_model_name %in% names(models_list))
att_vars <- models_list[[attribution_model_name]]
att_formula <- build_formula(att_vars)
cat("Attribution response function:\n"); print(att_formula)
att_fit  <- fixest::feols(att_formula, data = data_panel, cluster = cluster_formula, nthreads = 1)
fe <- fixef(att_fit)
state_month_fe <- data.frame(
  state_month = names(fe$state_month),
  fe_value    = as.numeric(fe$state_month)
)
state_year_fe <- data.frame(
  state_year = names(fe$state_year),
  fe_value   = as.numeric(fe$state_year)
)
write_csv(state_month_fe, file.path(output_dir, "fe_state_month_CDD-Ys-Ms.csv"))
write_csv(state_year_fe,  file.path(output_dir, "fe_state_year_CDD-Ys-Ms.csv"))

# [EDIT 2] `coef(att_fit)[predictor]` returned NA: in a spline model `cdd` is
# not a regressor, its basis columns are. This line is print-only — 10d predicts
# through the basis — so the fix changes the console output, nothing else.
beta_cdd <- coef(att_fit)[basis_names]
cat("CDD spline coefficients used:\n"); print(beta_cdd)

# ----- 10c. Load and prepare counterfactual CDD -----
cat("Loading counterfactual CDD from:\n  ", attribution_data_path, "\n")
cf <- read_csv(attribution_data_path, show_col_types = FALSE)
required_cf <- c("state", "month", scenario_col, member_col, cdd_col)
missing_cf  <- setdiff(required_cf, names(cf))
if (length(missing_cf) > 0)
  stop("Counterfactual file missing columns: ", paste(missing_cf, collapse = ", "))
cf <- cf %>%
  rename(any_of(c(scenario = scenario_col, member = member_col, cdd_cf = cdd_col)))
# Full names -> panel codes (unmapped states keep their name and are filtered out)
cf$state <- dplyr::recode(cf$state, !!!name_to_code)
cf <- cf %>%
  filter(state %in% study_area_states,
         month %in% study_period_months,
         scenario %in% c(all_label, nat_label))
# ---- keep only months present in BOTH scenarios (matched coverage) ----
common_months <- Reduce(intersect, split(cf$month, cf$scenario))
dropped_months <- setdiff(unique(cf$month), common_months)
if (length(dropped_months) > 0)
  cat("Dropping months not in both scenarios:",
      paste(sort(dropped_months), collapse = ", "), "\n")
cf <- cf %>% filter(month %in% common_months)
# Checks
present_scn <- unique(cf$scenario)
if (!all(c(all_label, nat_label) %in% present_scn))
  stop("scenario column must contain both '", all_label, "' and '", nat_label,
       "'. Found: ", paste(present_scn, collapse = ", "))
missing_states <- setdiff(study_area_states, unique(cf$state))
if (length(missing_states) > 0)
  cat("WARNING: study states absent from CDD file:",
      paste(missing_states, collapse = ", "), "\n")
# Collapse to per-member climatology: one CDD per state-month-scenario-member
cf_clim <- cf %>%
  group_by(state, month, scenario, member) %>%
  summarise(cdd_cf = mean(cdd_cf, na.rm = TRUE), .groups = "drop")
cat(sprintf("Counterfactual: %d states, %d months, %d members, scenarios: %s\n",
            n_distinct(cf_clim$state), n_distinct(cf_clim$month),
            n_distinct(cf_clim$member),
            paste(sort(unique(cf_clim$scenario)), collapse = ", ")))
# ----- 10d. Predict demand per scenario / member -----
sm_map <- data_panel %>% distinct(state, month, state_month)
sy_map <- data_panel %>%
  distinct(state, year, state_year) %>%
  filter(year %in% attribution_reference_years)
newdata <- cf_clim %>%
  inner_join(sm_map, by = c("state", "month")) %>%                 # + state_month
  inner_join(sy_map, by = "state", relationship = "many-to-many")  # + state_year (one row / ref year)
# rebuild the identical basis on counterfactual CDD (stored knots → linear tails)
cf_basis <- predict(cdd_basis, newx = newdata$cdd_cf)
newdata[basis_names] <- as.data.frame(cf_basis)
newdata$.pred <- as.numeric(predict(att_fit, newdata = newdata))
n_na <- sum(is.na(newdata$.pred))
if (n_na > 0)
  cat(sprintf("  Note: %d predictions NA (unseen FE level) and dropped\n", n_na))
# Average the baseline over reference years -> one prediction per member
pred_member <- newdata %>%
  filter(!is.na(.pred)) %>%
  group_by(state, month, scenario, member) %>%
  summarise(cdd              = mean(cdd_cf),
            predicted_demand = mean(.pred),
            .groups = "drop") %>%
  arrange(state, month, scenario, member)
# ----- 10e. Save -----
write_csv(pred_member, attribution_out_path)
cat("Saved attribution predictions to:\n  ", attribution_out_path, "\n")
cat(sprintf("Rows: %d (state x month x scenario x member)\n", nrow(pred_member)))
pred_member %>%
  group_by(state, scenario) %>%
  summarise(mean_demand = mean(predicted_demand), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = scenario, values_from = mean_demand) %>%
  mutate(attributable = .data[[all_label]] - .data[[nat_label]]) %>%
  print(n = Inf)

# ----- 10f. Monte Carlo propagation of CDD-coefficient uncertainty -----
# Everything above uses a single point estimate of the CDD response
# (beta_cdd). This propagates the regression's OWN parameter uncertainty --
# the same two-way clustered vcov used throughout, already PD-patched by
# fixest (see the earlier VCOV-not-positive-definite discussion -- that patch
# is exactly what lets MASS::mvrnorm() draw from it below without erroring)
# -- through to the attributable-shortage number itself, by drawing many
# plausible coefficient vectors and recomputing the attribution for each one.
# Currently the ensemble spread (525 members) is the only source of spread in
# the ALL-vs-NAT numbers above; this adds regression uncertainty on top,
# which matters more the noisier the fitted response is (i.e. more for
# shortage than for demand).
#
# Kept fast by exploiting linearity: predicted demand is FE_i + B_i %*% beta,
# and every aggregation step from there (mean over reference years, then mean
# over month & member) is itself a mean -- so "mean prediction" is an AFFINE
# function of beta for a fixed (state, scenario) group. Averaging the basis
# rows once per (state, scenario) and multiplying by each drawn beta is
# mathematically identical to re-running the full ~490k-row prediction
# pipeline once per draw, without the cost of actually doing so.
#
# Requires the MASS package (a base "recommended" package -- normally already
# installed alongside R itself).
MC_N <- 1000
set.seed(20240601)

mc_beta_mat <- MASS::mvrnorm(
  MC_N,
  mu    = as.numeric(beta_cdd),
  Sigma = vcov(att_fit)[basis_names, basis_names]
)
colnames(mc_beta_mat) <- basis_names

# FE_i, recovered exactly from the point-estimate prediction (predict.fixest
# gives FE_i + B_i %*% beta_cdd, so subtracting the basis part back out
# leaves FE_i on its own).
newdata$.fe_contribution <- newdata$.pred -
  as.numeric(as.matrix(newdata[basis_names]) %*% beta_cdd)

mc_cell <- newdata %>%
  filter(!is.na(.pred)) %>%
  group_by(state, month, scenario, member) %>%
  summarise(across(all_of(basis_names), mean),
            fe_contribution = mean(.fe_contribution),
            .groups = "drop")

# per (state, month): mean over members only. This is what a month-resolved MC
# interval needs -- mirrors cell_basis in the outage/district script.
mc_state_month <- mc_cell %>%
  group_by(state, month, scenario) %>%
  summarise(across(all_of(basis_names), mean),
            fe_contribution = mean(fe_contribution),
            .groups = "drop")

# per state, additionally collapsed over month (the existing whole-period
# export below). Averaging the already-per-month rows equally reproduces
# exactly what the single-step version computed before this was split in two.
mc_agg <- mc_state_month %>%
  group_by(state, scenario) %>%
  summarise(across(all_of(basis_names), mean),
            fe_contribution = mean(fe_contribution),
            .groups = "drop")

mc_all <- mc_agg %>% filter(scenario == all_label) %>% arrange(state)
mc_nat <- mc_agg %>% filter(scenario == nat_label) %>% arrange(state)
stopifnot(identical(mc_all$state, mc_nat$state))

B_diff  <- as.matrix(mc_all[basis_names]) - as.matrix(mc_nat[basis_names])   # n_states x 4
fe_diff <- mc_all$fe_contribution - mc_nat$fe_contribution                  # n_states

# attributable_r(state) = B_diff(state) %*% beta_r + fe_diff(state), all draws at once
mc_attributable <- fe_diff + B_diff %*% t(mc_beta_mat)   # n_states x MC_N
rownames(mc_attributable) <- mc_all$state

attribution_uncertainty <- data.frame(
  state     = mc_all$state,
  point_est = fe_diff + as.numeric(B_diff %*% as.numeric(beta_cdd)),   # matches "attributable" above
  mc_mean   = rowMeans(mc_attributable),
  mc_lo     = apply(mc_attributable, 1, quantile, probs = 0.025),
  mc_hi     = apply(mc_attributable, 1, quantile, probs = 0.975)
)
cat("\nAttributable", dependent_variable, "with Monte Carlo (95%) interval over CDD-coefficient uncertainty:\n")
print(attribution_uncertainty, row.names = FALSE)

attribution_uncertainty_path <- file.path(output_dir, "attribution_uncertainty_by_state.csv")
write_csv(attribution_uncertainty, attribution_uncertainty_path)
cat("\nSaved attribution uncertainty to:\n  ", attribution_uncertainty_path, "\n")

# ---- month-resolved version: same mc_beta_mat draws, evaluated per (state,
# month) instead of collapsed over the whole period. Not independent monthly
# noise -- every month shares the same underlying coefficient draws, so
# adjacent months are correlated; this is that one shared uncertainty
# evaluated at each month's own basis vector. ----
mc_all_month <- mc_state_month %>% filter(scenario == all_label) %>% arrange(state, month)
mc_nat_month <- mc_state_month %>% filter(scenario == nat_label) %>% arrange(state, month)
stopifnot(identical(mc_all_month$state, mc_nat_month$state),
          identical(mc_all_month$month, mc_nat_month$month))

B_diff_month  <- as.matrix(mc_all_month[basis_names]) - as.matrix(mc_nat_month[basis_names])
fe_diff_month <- mc_all_month$fe_contribution - mc_nat_month$fe_contribution

mc_attributable_month <- fe_diff_month + B_diff_month %*% t(mc_beta_mat)   # n_state_months x MC_N

attribution_uncertainty_month <- data.frame(
  state     = mc_all_month$state,
  month     = mc_all_month$month,
  point_est = fe_diff_month + as.numeric(B_diff_month %*% as.numeric(beta_cdd)),
  mc_mean   = rowMeans(mc_attributable_month),
  mc_lo     = apply(mc_attributable_month, 1, quantile, probs = 0.025),
  mc_hi     = apply(mc_attributable_month, 1, quantile, probs = 0.975)
)

attribution_uncertainty_month_path <- file.path(output_dir, "attribution_uncertainty_by_state_month.csv")
write_csv(attribution_uncertainty_month, attribution_uncertainty_month_path)
cat("Saved month-resolved attribution uncertainty to:\n  ",
    attribution_uncertainty_month_path, "\n")

# ---- mean relative attribution-CI width -----------------------------------
# The direct analog of mean_rel_pi_pct (section 11), but applied to the
# MC-propagated ATTRIBUTION interval itself (the dashed band in the timeseries
# figures) rather than the model's own prediction interval. One scale-free
# percentage per model, computed the same way for demand and shortage, so the
# two numbers are directly quotable side by side: "the mean relative width of
# the propagated attribution uncertainty was X% for demand vs Y% for
# shortage." Denominator is each state-month's own ACT (factual) predicted
# level -- the same fe_contribution + basis %*% beta_cdd construction used for
# the point estimates above -- not a single global mean, so a state-month with
# a naturally small or large level is not over- or under-weighted.
mc_all_month_level <- mc_all_month$fe_contribution +
  as.numeric(as.matrix(mc_all_month[basis_names]) %*% as.numeric(beta_cdd))
bad_level <- mc_all_month_level <= 0
if (any(bad_level)) {
  cat(sprintf("  NOTE: %d of %d state-months have a non-positive ACT predicted level and are excluded from the ratio below\n",
              sum(bad_level), length(bad_level)))
}
rel_width <- ifelse(bad_level, NA_real_,
                    (attribution_uncertainty_month$mc_hi - attribution_uncertainty_month$mc_lo) /
                      mc_all_month_level)
mean_rel_attrib_ci_pct <- 100 * mean(rel_width, na.rm = TRUE)
cat(sprintf("Mean relative width of the attribution MC 95%% CI: %.1f%% (mean over %d state-months of (mc_hi - mc_lo) / ACT predicted level)\n",
            mean_rel_attrib_ci_pct, sum(!bad_level)))

# =============================================================================
# [EDIT 3 — new] 11. Export the model-implied DEMAND response to CDD
# =============================================================================
# Writes the three CSVs the Python response-curve figure reads. It only reads
# from the objects above and writes new files — no table the other plots use is
# touched.
#   response_curve_<MODEL>.csv        cdd, fit, se, lo, hi, pred_lo, pred_hi,
#                                      resid_sigma, baseline
#   response_curve_knots_<MODEL>.csv  value, kind
#   response_curve_observed_cdd.csv   state, year, month, cdd, weight
#
# TWO DIFFERENT BANDS ARE WRITTEN — DO NOT CONFUSE THEM:
#
# 1) `lo`/`hi` — a CONFIDENCE band on the mean curve: spline-coefficient
#    uncertainty only, se(x) = sqrt(B(x) V B(x)'), V = vcov[ns, ns] (two-way
#    clustered — see cluster_formula). It answers "how precisely do we know
#    the average response function?" It excludes uncertainty in the
#    fixed-effect baseline, which shifts the whole curve vertically without
#    changing its shape.
#
# 2) `pred_lo`/`pred_hi` — a PREDICTION band: adds the model's residual
#    variance to the line's se, pred_se(x) = sqrt(se(x)^2 + resid_sigma^2).
#    It answers "how far can one actual state-month land from the fitted
#    curve?" For a low-R^2 outcome (e.g. shortage) this is the band that
#    shows the reviewer how much variation CDD leaves unexplained, even where
#    the average curve itself (band 1) is estimated fairly tightly. This
#    band is necessarily wider than the confidence band, often by a lot when
#    R^2 is low, and that gap IS the point.
#
# Both bands still exclude fixed-effect baseline uncertainty (see baseline,
# below) and any specification uncertainty (how much the curve would move
# under a different model in models_list) — those are separate exercises.
#
# WHAT THE Y AXIS IS: predicted electricity requirement for a state-month
# sitting at the observation-weighted mean state-month + state-year fixed
# effect -- i.e. an average state-month, not any particular state.
# =============================================================================
RESPONSE_MODEL <- "CDD-Ys-Ms"
CI_Z <- 1.96

stopifnot(RESPONSE_MODEL %in% names(models_list),
          all(basis_names %in% models_list[[RESPONSE_MODEL]]))

rc_formula <- build_formula(models_list[[RESPONSE_MODEL]])
cat("\nResponse function:\n"); print(rc_formula)
rc_fit <- fixest::feols(rc_formula, data = data_panel, cluster = cluster_formula, nthreads = 1)

rc_beta <- coef(rc_fit)[basis_names]
rc_V    <- vcov(rc_fit)[basis_names, basis_names]

# ---- baseline ----
# With FE specified as `| state_month + state_year`, feols absorbs the
# intercept into those group levels instead of reporting one scalar -- there
# is no "(Intercept)" in coef(rc_fit) to read off. "Baseline" below is
# therefore a CONSTRUCTED number, not a model parameter: the
# observation-weighted average of every fitted state_month/state_year level.
# fitted_i = FE_i + B_i %*% beta, so mean(fitted) - colMeans(B) %*% beta gives
# exactly that average FE level, without unpacking fixef() by hand -- and it
# stays correct whatever FE structure the chosen model uses.
rc_fitted <- predict(rc_fit, sample = "original")
rc_keep   <- !is.na(rc_fitted)
rc_B_obs  <- as.matrix(data_panel[rc_keep, basis_names, drop = FALSE])
rc_baseline <- mean(rc_fitted[rc_keep]) - as.numeric(colMeans(rc_B_obs) %*% rc_beta)

cat(sprintf("baseline (obs-weighted mean FE): %.2f  [%s, n = %d]\n",
            rc_baseline, dependent_variable, sum(rc_keep)))

# ---- curve on a fine grid ----
rc_grid <- seq(min(data_panel$cdd, na.rm = TRUE),
               max(data_panel$cdd, na.rm = TRUE), length.out = 400)
rc_B <- predict(cdd_basis, newx = rc_grid)      # stored knots -> identical basis
colnames(rc_B) <- basis_names

rc_fit_vals <- rc_baseline + as.numeric(rc_B %*% rc_beta)
rc_se_vals  <- sqrt(pmax(rowSums((rc_B %*% rc_V) * rc_B), 0))   # sqrt(diag(B V B'))

# Residual sd of the fitted model (unexplained variation left after CDD + FEs).
# This is what turns the confidence band above into a prediction band: a
# single state-month's shortage scatters around the curve with roughly this
# much spread, on top of how uncertain the curve's mean is.
rc_sigma      <- sigma(rc_fit)
rc_pred_se    <- sqrt(rc_se_vals^2 + rc_sigma^2)

# ---- single-number summaries, for comparing THIS run (e.g. shortage) against
# another run of the same script (e.g. demand) where the outcome is on a
# different scale ----
#
# CV(RMSE): residual noise as a % of the mean observed outcome. Standard
# scale-free noise metric (the same idea used e.g. in ASHRAE's building-energy
# calibration guideline) -- dividing by the mean, rather than reporting sigma
# in raw units, is what makes it comparable across two dependent variables
# with different magnitudes (shortage is typically much smaller in absolute
# terms than total demand).
rc_mean_y   <- mean(data_panel[[dependent_variable]][rc_keep])
cv_rmse_pct <- 100 * rc_sigma / rc_mean_y

# Mean relative half-width of the PREDICTION interval (pred_lo/pred_hi),
# evaluated at the panel's actual observed CDD values rather than the 400-point
# plotting grid above. The plotting grid spans min-to-max CDD uniformly, which
# oversamples the extremes where few real observations sit and the band is
# inflated by near-extrapolation; averaging over the real data instead gives
# the number that matches what "the band is wider for shortage" is actually
# claiming about your sample.
obs_B           <- as.matrix(data_panel[rc_keep, basis_names, drop = FALSE])
obs_se_line     <- sqrt(pmax(rowSums((obs_B %*% rc_V) * obs_B), 0))
obs_pred_se     <- sqrt(obs_se_line^2 + rc_sigma^2)
mean_rel_pi_pct <- 100 * mean(2 * CI_Z * obs_pred_se) / rc_mean_y

# ---- cluster (block) bootstrap band ----
# An empirical alternative to the two analytic bands above: resample whole
# STATES with replacement (not individual rows -- a state's observations are
# correlated over time, so resampling rows would understate uncertainty the
# same way un-clustered SEs would), refit on each resample, and take
# percentiles of the resulting curves. This doesn't assume asymptotic
# normality the way the delta-method band does, and it additionally captures
# uncertainty in the fixed-effect baseline (rc_baseline), which both bands
# above hold fixed at its point estimate. It's also unaffected by the
# non-positive-definite CGM issue, since it never builds that matrix at all --
# each replicate's point estimate is all that's used.
#
# A drawn state is relabelled with a unique suffix so that a state drawn
# twice becomes two distinct synthetic units (otherwise feols would silently
# collapse the duplicate rows onto the same fixed-effect levels instead of
# treating them as two independent draws). state_month/state_year are
# rebuilt from the relabelled state so their FE levels stay consistent with
# it -- reusing the original columns here would leave every duplicate copy
# pointing at the SAME state_month/state_year level, silently merging draws
# that are supposed to be independent.
BOOT_REPS <- 500   # drop to ~100 for a quick first pass, then rerun at 500+
set.seed(20240601) # bootstrap draws are random -- fix the seed for reproducibility

boot_states <- unique(data_panel$state)
boot_curves <- matrix(NA_real_, nrow = length(rc_grid), ncol = BOOT_REPS)
boot_ok     <- 0

for (b in seq_len(BOOT_REPS)) {
  drawn <- sample(boot_states, length(boot_states), replace = TRUE)
  
  boot_panel <- purrr::map2_dfr(drawn, seq_along(drawn), function(s, i) {
    data_panel %>%
      filter(state == s) %>%
      mutate(
        state       = paste0(s, "__", i),
        state_month = paste(state, month, sep = "_"),
        state_year  = paste(state, year, sep = "_")
      )
  })
  
  boot_fit <- tryCatch(
    fixest::feols(rc_formula, data = boot_panel, nthreads = 1),
    error = function(e) NULL
  )
  if (is.null(boot_fit)) next
  
  boot_beta <- tryCatch(coef(boot_fit)[basis_names], error = function(e) NULL)
  if (is.null(boot_beta) || anyNA(boot_beta)) next
  
  boot_fitted <- predict(boot_fit, sample = "original")
  boot_keep   <- !is.na(boot_fitted)
  if (sum(boot_keep) == 0) next
  
  boot_B_obs    <- as.matrix(boot_panel[boot_keep, basis_names, drop = FALSE])
  boot_baseline <- mean(boot_fitted[boot_keep]) - as.numeric(colMeans(boot_B_obs) %*% boot_beta)
  
  boot_curves[, b] <- boot_baseline + as.numeric(rc_B %*% boot_beta)  # rc_B: same fixed grid basis
  boot_ok <- boot_ok + 1
}

cat(sprintf("Cluster bootstrap: %d/%d replicates fit successfully\n", boot_ok, BOOT_REPS))

boot_lo <- apply(boot_curves, 1, quantile, probs = 0.025, na.rm = TRUE)
boot_hi <- apply(boot_curves, 1, quantile, probs = 0.975, na.rm = TRUE)

curve_out <- data.frame(
  cdd                = rc_grid,
  fit                = rc_fit_vals,
  se                 = rc_se_vals,
  lo                 = rc_fit_vals - CI_Z * rc_se_vals,
  hi                 = rc_fit_vals + CI_Z * rc_se_vals,
  pred_lo            = rc_fit_vals - CI_Z * rc_pred_se,
  pred_hi            = rc_fit_vals + CI_Z * rc_pred_se,
  boot_lo            = boot_lo,
  boot_hi            = boot_hi,
  boot_reps_ok        = boot_ok,
  resid_sigma        = rc_sigma,
  cv_rmse_pct        = cv_rmse_pct,
  mean_rel_pi_pct    = mean_rel_pi_pct,
  baseline           = rc_baseline,
  model_name         = RESPONSE_MODEL,
  dependent_variable = dependent_variable
)
curve_path <- file.path(output_dir, sprintf("response_curve_%s.csv", RESPONSE_MODEL))
write_csv(curve_out, curve_path)

cat(sprintf("Residual sigma: %.2f | mean CI half-width: %.2f | mean PI half-width: %.2f | mean bootstrap half-width: %.2f\n",
            rc_sigma, CI_Z * mean(rc_se_vals), CI_Z * mean(rc_pred_se), mean((boot_hi - boot_lo) / 2)))
cat(sprintf("CV(RMSE): %.1f%%  |  mean relative PI width (obs-weighted): %.1f%%  [%s]\n",
            cv_rmse_pct, mean_rel_pi_pct, dependent_variable))

# ---- knots ----
knots_out <- rbind(
  data.frame(value = as.numeric(attr(cdd_basis, "knots")),          kind = "interior"),
  data.frame(value = as.numeric(attr(cdd_basis, "Boundary.knots")), kind = "boundary")
)
knots_path <- file.path(output_dir,
                        sprintf("response_curve_knots_%s.csv", RESPONSE_MODEL))
write_csv(knots_out, knots_path)

# ---- the panel's own CDD, for the density panel ----
obs_cdd <- data_panel %>%
  filter(!is.na(cdd), !is.na(.data[[dependent_variable]])) %>%
  transmute(state, year, month, cdd, weight = 1)
obs_path <- file.path(output_dir, "response_curve_observed_cdd.csv")
write_csv(obs_cdd, obs_path)

cat("\nWrote:\n ", curve_path, "\n ", knots_path, "\n ", obs_path, "\n")
cat(sprintf("CDD range %.1f-%.1f | fit range %.1f-%.1f | mean se %.2f\n",
            min(rc_grid), max(rc_grid), min(rc_fit_vals), max(rc_fit_vals),
            mean(rc_se_vals)))

# =============================================================================
# 12. K-fold cross-validated (out-of-sample) R^2
# =============================================================================
# Every R^2 reported so far (r2_squared, adj_r_squared, within_r_squared, and
# the CV(RMSE)/relative-PI numbers above) is IN-SAMPLE: the model is fit and
# evaluated on the same rows. With ~150+ state_month and state_year fixed
# effects soaking up degrees of freedom, in-sample R^2 can look better than
# the model's real predictive power warrants -- more so for a noisier outcome
# like shortage, since there's more residual variance for the FEs to
# opportunistically fit. This reruns the fit K times, holding out a different
# 1/K of rows each time and scoring predictions the fit never saw, which is a
# more honest number to put next to the in-sample R^2 the reviewer questioned.
#
# Folds are assigned by INDIVIDUAL ROW, not by whole state or year (unlike the
# bootstrap above) -- each state_month/state_year fixed-effect level has
# several rows behind it (one per year, or one per month, respectively), so a
# random ~10% holdout leaves nearly every FE level identified in training.
# Held-out rows whose FE level is nonetheless entirely absent from training
# (rare, but possible for a sparsely populated level) get NA from predict()
# and are dropped from the score -- the count below reports how many that was.
# =============================================================================
CV_K <- 10
set.seed(20240601)

cv_ids  <- which(rc_keep)              # same estimation sample as the response curve
cv_fold <- sample(rep_len(seq_len(CV_K), length(cv_ids)))
cv_pred <- rep(NA_real_, nrow(data_panel))

for (k in seq_len(CV_K)) {
  test_idx  <- cv_ids[cv_fold == k]
  
  fit_k <- tryCatch(
    fixest::feols(rc_formula, data = data_panel[-test_idx, ], nthreads = 1),
    error = function(e) NULL
  )
  if (is.null(fit_k)) next
  
  pred_k <- tryCatch(
    as.numeric(predict(fit_k, newdata = data_panel[test_idx, ])),
    error = function(e) rep(NA_real_, length(test_idx))
  )
  cv_pred[test_idx] <- pred_k
}

cv_scored <- !is.na(cv_pred) & rc_keep
y_actual  <- data_panel[[dependent_variable]]

oos_r2 <- 1 - sum((y_actual[cv_scored] - cv_pred[cv_scored])^2) /
  sum((y_actual[cv_scored] - mean(y_actual[cv_scored]))^2)

cat(sprintf("\n%d-fold cross-validated R^2: %.4f  (in-sample R^2 was %.4f)  [%d/%d obs scored, %s]\n",
            CV_K, oos_r2, r2(rc_fit, type = "cor2"),
            sum(cv_scored), sum(rc_keep), dependent_variable))

cv_out <- data.frame(
  model_name         = RESPONSE_MODEL,
  dependent_variable = dependent_variable,
  cv_k               = CV_K,
  in_sample_r2       = r2(rc_fit, type = "cor2"),
  oos_r2             = oos_r2,
  n_scored           = sum(cv_scored),
  n_total            = sum(rc_keep)
)
cv_path <- file.path(output_dir, sprintf("cv_r2_%s.csv", RESPONSE_MODEL))
write_csv(cv_out, cv_path)
cat("Saved cross-validated R^2 to:\n  ", cv_path, "\n")