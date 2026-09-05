setwd("~/code-location") #ensure your working directory is set
# ---------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------
library(tidyverse)
library(fixest)
library(splines)
library(R.utils)

predictor = 'cdd'

# ----- Specify the states in the study area -----
# study_area_states <- c()
study_area_states <- c("Punjab","Haryana","Rajasthan","NCT of Delhi","Uttar Pradesh","Odisha","Madhya Pradesh","Gujarat","Chhattisgarh","Telangana","Jharkhand","Uttarakhand", "Bihar") # full study area

# ----- Specify the months in the study area -----
study_period_months <- c(1,2,3,4,5,6,7,8,9,10,11,12)


dependent_variable <- 'outage_minutes'

# ---------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------
data_path  <- "data/prayas_out/esmi_outage_cdd_panel.csv"
output_dir <- "data/fit_statistics/outage_models"

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

# --- Restrict to rows with a full day's worth of observation ---
observed_minutes_col <- "observed_minutes"   # column recording minutes actually monitored
full_day_minutes     <- 1440L                # 24h * 60min

if (!observed_minutes_col %in% names(data_panel)) {
  stop(sprintf(
    "Column '%s' not found in data_panel. Available columns:\n  %s",
    observed_minutes_col, paste(names(data_panel), collapse = ", ")
  ))
}

n_before <- nrow(data_panel)

completeness_report <- data_panel %>%
  mutate(is_full_day = .data[[observed_minutes_col]] == full_day_minutes) %>%
  group_by(state) %>%
  summarise(
    n_total      = n(),
    n_full_day   = sum(is_full_day, na.rm = TRUE),
    n_dropped    = n_total - n_full_day,
    pct_full_day = round(100 * n_full_day / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(state)

cat(sprintf("\n--- Observation completeness (rows where %s == %d) ---\n",
            observed_minutes_col, full_day_minutes))
print(completeness_report, n = Inf)

data_panel <- data_panel %>% filter(.data[[observed_minutes_col]] == full_day_minutes)

n_after <- nrow(data_panel)
cat(sprintf(
  "\nKept %d of %d rows (%.1f%%) with a full %d-minute observation window; dropped %d rows.\n\n",
  n_after, n_before, 100 * n_after / n_before, full_day_minutes, n_before - n_after
))

# --- Drop LOCATIONS with too few observed days -------------------------------
# Filtering at location level rather than district level. A district is not
# dropped directly: it disappears only if every one of its locations falls
# short, and a district with a mix keeps the locations that clear the bar.
#
# Counts DISTINCT DATES per location. Applied AFTER the full-1440-minute rule,
# so the threshold refers to complete days that enter the estimation sample.
min_location_days <- 365L
location_col      <- "Location name"

for (needed in c("date", location_col))
  if (!needed %in% names(data_panel))
    stop("Column '", needed, "' not found -- cannot count observed days per location.")

data_panel$.loc <- data_panel[[location_col]]

location_days <- data_panel %>%
  group_by(state, district, .loc) %>%
  summarise(n_days   = n_distinct(date),
            n_rows   = n(),
            first_day = min(date, na.rm = TRUE),
            last_day  = max(date, na.rm = TRUE),
            .groups  = "drop")

keep_loc <- location_days %>% filter(n_days >= min_location_days)
drop_loc <- location_days %>% filter(n_days <  min_location_days) %>%
  arrange(state, district, desc(n_days))

# district-level before/after, so the knock-on effect is visible
district_before <- location_days %>% count(state, district, name = "n_loc_before")
district_after  <- keep_loc      %>% count(state, district, name = "n_loc_after")
district_change <- district_before %>%
  left_join(district_after, by = c("state", "district")) %>%
  mutate(n_loc_after = coalesce(n_loc_after, 0L),
         n_loc_lost  = n_loc_before - n_loc_after)

gone    <- district_change %>% filter(n_loc_after == 0)
thinned <- district_change %>% filter(n_loc_after > 0, n_loc_lost > 0)

cat(sprintf("\n--- Location coverage filter (>= %d observed days) ---\n",
            min_location_days))
cat(sprintf("Locations: %d kept, %d dropped, of %d\n",
            nrow(keep_loc), nrow(drop_loc), nrow(location_days)))

if (nrow(drop_loc) > 0) {
  cat("\nDropped locations:\n")
  drop_loc %>%
    mutate(span = sprintf("%s to %s", first_day, last_day)) %>%
    select(state, district, location = .loc, n_days, n_rows, span) %>%
    print(n = Inf)
}

if (nrow(gone) > 0) {
  cat(sprintf("\nDistricts lost entirely (no location reached %d days): %d\n",
              min_location_days, nrow(gone)))
  gone %>% select(state, district, n_loc_before) %>% print(n = Inf)
} else {
  cat("\nNo district lost entirely.\n")
}

if (nrow(thinned) > 0) {
  cat("\nDistricts that keep some locations but lose others:\n")
  thinned %>% select(state, district, n_loc_before, n_loc_after, n_loc_lost) %>%
    print(n = Inf)
}

if (nrow(keep_loc) == 0)
  stop("No location clears ", min_location_days, " observed days -- lower the threshold.")

n_before_loc <- nrow(data_panel)
data_panel <- data_panel %>%
  semi_join(keep_loc, by = c("state", "district", ".loc")) %>%
  select(-.loc)

cat(sprintf("\nRetained: %d locations in %d districts across %d states; %s of %s rows (%.1f%%)\n",
            nrow(keep_loc),
            n_distinct(paste(data_panel$state, data_panel$district)),
            n_distinct(data_panel$state),
            format(nrow(data_panel), big.mark = ","),
            format(n_before_loc, big.mark = ","),
            100 * nrow(data_panel) / n_before_loc))
cat(sprintf("Observed days per kept location: %d-%d (median %d)\n\n",
            min(keep_loc$n_days), max(keep_loc$n_days),
            as.integer(median(keep_loc$n_days))))

# --- Export locations that pass the filter and are within the study area ---
locations_path <- file.path(output_dir, "locations_kept.csv")
keep_loc %>%
  transmute(
    state,
    district,
    location  = .loc,
    n_days,                # distinct observed days (post-filters)
    first_day,
    last_day
  ) %>%
  arrange(state, district, location) %>%
  write_csv(locations_path)
cat("Saved kept-locations table to:", locations_path, "\n")

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

sapply(data_panel[c("cdd_tmean","cdd_dewpt","cdd_swbgt","tmean_c","rh")],
       function(x) sum(!is.na(x)))

cdd_ns   <- add_ns("cdd_tmean")
dewpt_ns <- add_ns("cdd_dewpt")
swbgt_ns <- add_ns("cdd_swbgt")
tmean_ns <- add_ns("tmean_c",  df = 4L)
rh_ns    <- add_ns("rh",      df = 4L)

# keep the CDD aliases the rest of the script already uses, so nothing downstream changes
cdd_basis          <- ns_bases[["cdd_tmean"]]
cdd_knots          <- attr(cdd_basis, "knots")
cdd_boundary_knots <- attr(cdd_basis, "Boundary.knots")

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
  'D-Y' = c(
    'district',
    'year'
  ),
  'D-M' = c(
    'district',
    'month'
  ),
  'Yd-Md'= c(
    'district_month',
    'district_year'
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
    swbgt_ns,
    'state_month',
    'state_year'
  ),
  'T-Ys-Ms' = c(
    tmean_ns,
    'state_month',
    'state_year'
  ),
  'CDD-D-Y' = c(
    cdd_ns,
    'district',
    'year'
  ),
  'CDD-D-M' = c(
    cdd_ns,
    'district',
    'month'
  ),
  'CDD-Yd-Md' = c(
    cdd_ns,
    'district_month',
    'district_year'
  ),
  'CDD-RH-Yd-Md' = c(
    cdd_ns,
    rh_ns,
    'district_month',
    'district_year'
  ),
  'dCDD-Yd-Md' = c(
    dewpt_ns,
    'district_month',
    'district_year'
  ),
  'wCDD-Yd-Md' = c(
    swbgt_ns,
    'district_month',
    'district_year'
  ),
  'T-Yd-Md' = c(
    tmean_ns,
    'district_month',
    'district_year'
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
                                     'state',
                                     'district',
                                     'district_month',
                                     'district_year'))
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
# Two-way cluster by district and year -- the district-level analogue of the
# state x year clustering used in the demand/shortage script, so the two
# pipelines make the same inference assumptions. Lets outage minutes be
# correlated within a district across time (a persistent local grid issue)
# AND across districts within the same year (a shared heatwave or monsoon
# event), instead of assuming independence along whichever dimension isn't
# clustered on. Defined once here and reused in every feols() call below,
# including the one that feeds the attribution Monte Carlo in section 10h.
cluster_formula <- ~district + year

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
  
  # --- observed vs fitted outage stored ---
  
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
       title = "Northwest core: outage over time (2022 heatwave shaded)") +
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

ggsave(file.path(output_dir, "obs_vs_fitted_by_model.png"), p_fit,
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

ggsave(file.path(output_dir, "timeseries_obs_vs_fitted_by_model.png"), p_ts,
       width = 10, height = 8, dpi = 150)

# ---------------------------------------------------------------
# 9. Combine and save
# ---------------------------------------------------------------

coef_table        <- bind_rows(all_coefs)
fit_stats_table   <- bind_rows(all_fit_stats)
outage_table   <- bind_rows(all_model_fit)

coef_path         <- file.path(output_dir, "all_models_coefficients.csv")
fit_stats_path    <- file.path(output_dir, "all_models_fit_stats.csv")
outage_path    <- file.path(output_dir, "all_models_fit.csv")

write_csv(coef_table,      coef_path)
write_csv(fit_stats_table, fit_stats_path)
write_csv(outage_table, outage_path)

cat("\n=============================================\n")
cat("Saved coefficient table to:      ", coef_path,      "\n")
cat("Saved fit statistics table to:   ", fit_stats_path, "\n")
cat("Saved outage table to:", outage_path, "\n")
cat("Study Area:", study_area_states, "\n")
cat("Study Months:", study_period_months, "\n")

# --- Basis details stored ---

basis_path_grid <- file.path(output_dir, "cdd_basis_grid.csv")
basis_path_long <- file.path(output_dir, "basis_vcov_long.csv")

# 1) exact R basis on a fine CDD grid (sidesteps any R/Python basis mismatch)
grid  <- seq(min(data_panel$cdd_tmean, na.rm = TRUE), max(data_panel$cdd_tmean, na.rm = TRUE), length.out = 400)
Bgrid <- predict(cdd_basis, newx = grid)
colnames(Bgrid) <- cdd_ns                      # cdd_ns1 … cdd_ns4
write.csv(cbind(cdd = grid, Bgrid), basis_path_grid, row.names = FALSE)


spline_models <- names(models_list)[sapply(models_list, \(v) all(cdd_ns %in% v))]
vcov_long <- purrr::map_dfr(spline_models, function(nm) {
  V  <- vcov(fixest::feols(build_formula(models_list[[nm]]), data = data_panel,
                           cluster = cluster_formula))[cdd_ns, cdd_ns]
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
# 10. Climate attribution: predict outage under ALL vs NAT scenarios
#     DISTRICT version -- drop-in replacement for the old section 10.
# ================================================================
# Feeds climate-model district CDD (all-forcings "ALL" vs natural-only "NAT")
# through one fitted response function. For a given district-day the only thing
# differing between the two predictions is the CDD term, so ALL - NAT isolates
# the outage attributable to anthropogenic warming.
#
# SCOPE
# -----
# Attribution is produced ONLY for district x year-month cells that are BOTH
#   (a) observed in the ESMI panel (after every filter section 3 applied,
#       including the full-1440-minute-day rule), and
#   (b) present in the ensemble.
# Enforced by one inner join on (gadm_key, year, month). A district monitored
# only in 2014-15 gets nothing, because the ensemble starts in 2016; a district
# monitored throughout gets only its Jun-Dec 2016-2018 cells. Nothing is
# predicted for a year-month that was never observed.
#
# Matching is at MONTH granularity, not day: the ensemble runs on a 360-day
# calendar, so model day 2017-07-15 is not Gregorian 15 July 2017. Within a
# matched month, every model day is used.
#
# Because each cell carries its own year, the fixed-effect baseline is the ACTUAL
# fitted FE for that district-year-month. It is identical in both legs, so it
# cancels in ALL-NAT and only sets the level of the reported predictions.
#
# CDD file (all_district_cdd.csv.gz), one row per district-day-member-leg:
#   year, month, day, gadm_key, tmean_c, cdd_cellmean, cdd_tmean, leg, member
#   - gadm_key = "NAME_1||NAME_2", the same key the district weights were built
#     on, so district-name spelling cannot silently break the merge.
#   - member = 524 matched pairs (r012i1p2 dropped: corrupt NAT input file).
#   - year/month/day are the model's 360-day calendar stamp.
# Build it with:
#   python combine_district_cdd.py --check-only          # QC first
#   python combine_district_cdd.py --csv --gzip --lean   # -> .csv.gz
#
# ENSEMBLE SAMPLE: 14 year-months -- Jun-Aug 2016, Jun-Nov 2017, Jun + Sep-Dec
# 2018. Jul/Aug 2018 were dropped so the p2-p5 physics variants could be kept
# (105 -> 524 matched members). There is NO March-May, so this speaks to monsoon
# and post-monsoon heat, not the pre-monsoon peak.
#
# Two further differences from the old state-level section 10:
#   - Model CDD is DAILY and the response is a nonlinear spline, so predictions
#     are made per district-DAY and averaged after. Averaging CDD first would
#     evaluate f(E[CDD]) instead of E[f(CDD)] and bias the result.
#   - Predictions are formed as X %*% beta + fe_baseline rather than through
#     predict(newdata=), so ~14M daily rows never get expanded further.
# ----------------------------------------------------------------

# ----- 10a. Attribution configuration -----
attribution_model_name <- "CDD-Yd-Md"

# The panel variable the ns() basis was fitted on, and the CDD file column that
# feeds it. Keep them the same variable unless you know why you're crossing them:
#   cdd_tmean    = CDD of the district-mean temperature (panel default)
#   cdd_cellmean = area-weighted mean of per-cell CDD
panel_predictor     <- "cdd_tmean"
attribution_cdd_col <- "cdd_tmean"

attribution_data_path <- "data/attribution_runs/all_district_cdd.csv.gz"

scenario_col <- "leg"
member_col   <- "member"
all_label    <- "ALL"        # factual / all-forcings
nat_label    <- "NAT"        # counterfactual / natural

min_obs_days      <- 1L      # min observed days for a district-year-month to count
member_chunk_size <- 25L     # members per prediction chunk; lower if RAM is tight
clamp_to_support  <- FALSE   # TRUE truncates ensemble CDD to the fitted range

attribution_out_path   <- file.path(output_dir, "attribution_predictions_district.csv")
attribution_sum_path   <- file.path(output_dir, "attribution_summary_district.csv")
attribution_cells_path <- file.path(output_dir, "attribution_cells_used.csv")
attribution_dist_path  <- file.path(output_dir, "attribution_by_district.csv")

# ----- 10b. Checks, and re-fit the chosen response function -----
stopifnot(attribution_model_name %in% names(models_list))
att_vars <- models_list[[attribution_model_name]]

fe_all_names  <- c("month", "year", "state", "state_month", "state_year",
                   "district", "district_month", "district_year")
att_fe_vars   <- intersect(att_vars, fe_all_names)
att_pred_vars <- setdiff(att_vars, att_fe_vars)

att_basis <- ns_bases[[panel_predictor]]
if (is.null(att_basis))
  stop("No stored ns() basis for '", panel_predictor, "' -- add_ns() it in section 3.")
att_ns_cols <- paste0(panel_predictor, "_ns", seq_len(ncol(att_basis)))

# The ensemble carries temperature only, so humidity models cannot be driven by it.
if (!setequal(att_pred_vars, att_ns_cols))
  stop("Attribution model '", attribution_model_name, "' has predictors {",
       paste(att_pred_vars, collapse = ", "), "} but the ensemble can only supply {",
       paste(att_ns_cols, collapse = ", "), "}. The HadGEM3-A run provides tas only, ",
       "so RH / dew-point / sWBGT models can't be used for attribution.")

if (attribution_cdd_col != panel_predictor)
  cat("WARNING: pushing ensemble '", attribution_cdd_col, "' through a basis fitted on '",
      panel_predictor, "' -- these are different CDD definitions.\n", sep = "")

att_formula <- build_formula(att_vars)
cat("Attribution response function:\n"); print(att_formula)
att_fit  <- fixest::feols(att_formula, data = data_panel,
                          cluster = cluster_formula, nthreads = 1)
att_beta <- coef(att_fit)[att_ns_cols]
if (anyNA(att_beta))
  stop("Spline coefficients are NA (collinear basis?) -- check the fitted model.")

# Export district_month/district_year FE, mirroring the demand/shortage script's
# fe_state_month/fe_state_year export -- lets a Python per-district figure build
# a district baseline (avg district_month FE + avg district_year FE over a chosen
# study period) the same way the state-level figure does, without needing the
# ~14M-row ensemble/attribution machinery just to get at a fixed effect.
fe_att <- fixef(att_fit)
district_month_fe <- data.frame(
  district_month = names(fe_att$district_month),
  fe_value       = as.numeric(fe_att$district_month)
)
district_year_fe <- data.frame(
  district_year = names(fe_att$district_year),
  fe_value      = as.numeric(fe_att$district_year)
)
write_csv(district_month_fe, file.path(output_dir, "fe_district_month_CDD-Yd-Md.csv"))
write_csv(district_year_fe,  file.path(output_dir, "fe_district_year_CDD-Yd-Md.csv"))

# ----- 10c. Observed district x year-month cells (the ESMI side of the overlap) -----
panel_keyed <- data_panel %>% mutate(gadm_key = paste0(state, "||", district))

obs_cells <- panel_keyed %>%
  group_by(gadm_key, state, district, year, month) %>%
  summarise(n_obs_days      = n_distinct(date),
            n_locations     = n_distinct(.data[["Location name"]]),
            n_location_days = n(),
            mean_outage     = mean(.data[[dependent_variable]], na.rm = TRUE),
            .groups = "drop") %>%
  filter(n_obs_days >= min_obs_days)

cat(sprintf("\nObserved panel cells: %d district-year-month(s), %d district(s), years %s\n",
            nrow(obs_cells), n_distinct(obs_cells$gadm_key),
            paste(range(obs_cells$year), collapse = "-")))

# ----- 10d. Fixed-effect value for each observed cell -----
if (length(att_fe_vars) == 0) {
  cat("Response function has no fixed effects -- baseline = intercept\n")
  cell_fe <- obs_cells %>%
    transmute(gadm_key, year, month,
              fe_baseline = unname(coef(att_fit)["(Intercept)"]))
} else {
  fe_list <- fixef(att_fit)
  fe_dims <- names(fe_list)
  cat("FE dimensions in the response function:", paste(fe_dims, collapse = " + "), "\n")
  
  fe_rows <- panel_keyed %>%
    select(all_of(unique(c("gadm_key", "year", "month", fe_dims)))) %>%
    distinct()
  for (d in fe_dims)
    fe_rows[[paste0(".fe_", d)]] <- unname(fe_list[[d]][as.character(fe_rows[[d]])])
  fe_rows$.fe_sum <- rowSums(as.matrix(fe_rows[paste0(".fe_", fe_dims)]))
  
  n_fe_na <- sum(is.na(fe_rows$.fe_sum))
  if (n_fe_na > 0)
    cat(sprintf("  Note: %d cell(s) had an unseen FE level and are dropped\n", n_fe_na))
  
  cell_fe <- fe_rows %>%
    filter(!is.na(.fe_sum)) %>%
    group_by(gadm_key, year, month) %>%
    summarise(fe_baseline = mean(.fe_sum), .groups = "drop")
}

# ---- VERIFY the FE reconstruction against fixest itself ----
# pred = X %*% beta + fe_baseline must reproduce predict(att_fit) on the panel.
# If it does not, the fixef() decomposition is being summed wrongly and every
# predicted LEVEL below is meaningless. (The ALL-NAT difference would survive,
# since the baseline cancels, but no percentage of ACT outage would.)
fe_check <- panel_keyed %>%
  select(all_of(c("gadm_key", "year", "month", att_ns_cols))) %>%
  left_join(cell_fe, by = c("gadm_key", "year", "month"))
manual_pred <- as.numeric(as.matrix(fe_check[, att_ns_cols]) %*% att_beta) +
  fe_check$fe_baseline
fixest_pred <- as.numeric(predict(att_fit, sample = "original"))
if (length(fixest_pred) == length(manual_pred)) {
  dif <- abs(manual_pred - fixest_pred)
  cat(sprintf("FE reconstruction check: max |manual - fixest| = %.3e over %d rows\n",
              max(dif, na.rm = TRUE), sum(is.finite(dif))))
  if (max(dif, na.rm = TRUE) > 1e-6)
    warning("FE reconstruction does not match fixest::predict -- predicted LEVELS ",
            "(and every percentage) are wrong. The ALL-NAT difference is still valid.")
} else {
  cat("FE reconstruction check skipped (predict() returned a different length).\n")
}
cat(sprintf("Panel: observed mean %s = %.2f, fitted mean = %.2f, fitted range %.1f to %.1f\n",
            dependent_variable, mean(data_panel[[dependent_variable]], na.rm = TRUE),
            mean(fixest_pred, na.rm = TRUE),
            min(fixest_pred, na.rm = TRUE), max(fixest_pred, na.rm = TRUE)))
neg_fit <- mean(fixest_pred < 0, na.rm = TRUE)
if (neg_fit > 0.01)
  cat(sprintf("  NOTE: %.1f%% of the model's own fitted values are negative outage minutes\n",
              100 * neg_fit),
      "        (OLS on an outcome bounded at zero). Shares are computed as ratios of\n",
      "        SUMS below, which stays well defined as long as the total is positive.\n",
      sep = "")

cells <- obs_cells %>% inner_join(cell_fe, by = c("gadm_key", "year", "month"))
cat(sprintf("Cells with a usable FE baseline: %d of %d\n", nrow(cells), nrow(obs_cells)))

# ----- 10e. Load the ensemble and intersect with the observed cells -----
cat("\nLoading counterfactual district CDD from:\n  ", attribution_data_path, "\n")
cf_cols <- c("year", "month", "day", "gadm_key",
             attribution_cdd_col, scenario_col, member_col)

# fread needs R.utils to read .gz; fall back to the shell's gzip, then to readr.
read_cf <- function(path, cols) {
  is_gz   <- grepl("\\.gz$", path)
  have_dt <- requireNamespace("data.table", quietly = TRUE)
  if (have_dt && (!is_gz || requireNamespace("R.utils", quietly = TRUE)))
    return(as_tibble(data.table::fread(path, select = cols, showProgress = FALSE)))
  if (have_dt && is_gz) {
    unz <- if (nzchar(Sys.which("gzcat"))) "gzcat" else "gunzip -c"
    cat("  (fread via", unz, "-- install R.utils to read .gz directly)\n")
    return(as_tibble(data.table::fread(cmd = paste(unz, shQuote(path)),
                                       select = cols, showProgress = FALSE)))
  }
  cat("  (readr fallback -- slower; install data.table for large files)\n")
  read_csv(path, show_col_types = FALSE, col_select = all_of(cols), progress = FALSE)
}
cf <- read_cf(attribution_data_path, cf_cols)

missing_cf <- setdiff(cf_cols, names(cf))
if (length(missing_cf) > 0)
  stop("Counterfactual file missing columns: ", paste(missing_cf, collapse = ", "),
       ". If you combined WITHOUT --lean, the file has date/state/district instead.")

cf <- cf %>%
  rename(any_of(c(scenario = scenario_col,
                  member   = member_col,
                  cdd_cf   = attribution_cdd_col))) %>%
  filter(scenario %in% c(all_label, nat_label), !is.na(cdd_cf))

present_scn <- unique(cf$scenario)
if (!all(c(all_label, nat_label) %in% present_scn))
  stop("scenario column must contain both '", all_label, "' and '", nat_label,
       "'. Found: ", paste(present_scn, collapse = ", "))

ens_cells <- cf %>% distinct(gadm_key, year, month)
cat(sprintf("Ensemble: %s rows, %d district(s), %d year-month(s), %d member(s)\n",
            format(nrow(cf), big.mark = ","), n_distinct(cf$gadm_key),
            nrow(distinct(ens_cells, year, month)), n_distinct(cf$member)))

# ---- the overlap, reported before anything is dropped ----
obs_ym <- obs_cells %>% distinct(year, month)
ens_ym <- ens_cells %>% distinct(year, month)
cat("\nYear-months -- panel:", nrow(obs_ym), " ensemble:", nrow(ens_ym),
    " overlapping:", nrow(inner_join(obs_ym, ens_ym, by = c("year", "month"))), "\n")
ens_only_ym <- anti_join(ens_ym, obs_ym, by = c("year", "month"))
if (nrow(ens_only_ym))
  cat("  ensemble year-months with NO panel observation anywhere (not attributed):\n   ",
      paste(sprintf("%d-%02d", ens_only_ym$year, ens_only_ym$month), collapse = ", "), "\n")

panel_d <- distinct(obs_cells, gadm_key)
ens_d   <- distinct(ens_cells, gadm_key)
cat("Districts  -- panel:", nrow(panel_d), " ensemble:", nrow(ens_d),
    " both:", nrow(inner_join(panel_d, ens_d, by = "gadm_key")), "\n")
no_ens <- anti_join(panel_d, ens_d, by = "gadm_key")
if (nrow(no_ens)) {
  cat("  panel districts with NO ensemble CDD (add to DISTRICT_CROSSWALK and rebuild",
      "the weights if these matter):\n")
  print(no_ens, n = Inf)
}

# THE overlap join: district AND year AND month must all match.
cf <- cf %>% inner_join(cells, by = c("gadm_key", "year", "month"))
if (nrow(cf) == 0)
  stop("No overlap between the ensemble and the observed panel cells.")

used <- cf %>% distinct(gadm_key, state, district, year, month, n_obs_days, n_locations)
cat(sprintf("\nATTRIBUTING %d district-year-month cell(s), %d district(s), %s daily rows\n",
            nrow(used), n_distinct(used$gadm_key), format(nrow(cf), big.mark = ",")))
cat("Year-months actually used:\n")
used %>% count(year, month, name = "n_districts") %>%
  arrange(year, month) %>% print(n = Inf)
dropped_districts <- setdiff(unique(obs_cells$gadm_key), unique(used$gadm_key))
if (length(dropped_districts))
  cat(sprintf("Panel districts with no overlapping year-month (%d): %s\n",
              length(dropped_districts), paste(head(dropped_districts, 10), collapse = ", ")))
write_csv(used, attribution_cells_path)
cat("Cells used ->", attribution_cells_path, "\n")

# ---- how far outside the fitted CDD support does the ensemble sit? ----
bk <- attr(att_basis, "Boundary.knots")
frac_out <- mean(cf$cdd_cf < bk[1] | cf$cdd_cf > bk[2])
cat(sprintf("Ensemble CDD range %.1f-%.1f vs fitted support %.1f-%.1f; %.2f%% outside%s\n",
            min(cf$cdd_cf), max(cf$cdd_cf), bk[1], bk[2], 100 * frac_out,
            if (clamp_to_support) " (clamped)" else " (linear extrapolation)"))
if (clamp_to_support) cf$cdd_cf <- pmin(pmax(cf$cdd_cf, bk[1]), bk[2])

# ----- 10f. Predict per district-DAY, chunked over members -----
members <- sort(unique(cf$member))
chunks  <- split(members, ceiling(seq_along(members) / member_chunk_size))
cat(sprintf("\nPredicting in %d chunk(s) of <=%d members...\n", length(chunks), member_chunk_size))

pred_chunks <- vector("list", length(chunks))
basis_chunks <- vector("list", length(chunks))   # mean basis vector per cell-member, for section 10h
for (i in seq_along(chunks)) {
  d <- cf %>% filter(member %in% chunks[[i]])
  B <- predict(att_basis, newx = d$cdd_cf)          # stored knots -> identical basis
  d$.pred <- as.numeric(B %*% att_beta) + d$fe_baseline
  colnames(B) <- att_ns_cols
  d[att_ns_cols] <- as.data.frame(B)
  
  pred_chunks[[i]] <- d %>%
    group_by(gadm_key, state, district, year, month, scenario, member) %>%
    summarise(cdd                      = mean(cdd_cf),
              predicted_outage_minutes = mean(.pred),        # mean per model day
              n_model_days             = n(),
              n_obs_days               = first(n_obs_days),
              obs_outage               = first(mean_outage), # OBSERVED mean, ESMI
              .groups = "drop")
  # Mean basis vector per cell-member (same "mean per model day" step as above,
  # just on the 4 basis columns instead of the finished prediction -- this is
  # what lets section 10h redo the attribution cheaply for many draws of beta
  # without re-running predict() over the full ~14M daily rows again).
  basis_chunks[[i]] <- d %>%
    group_by(gadm_key, state, district, year, month, scenario, member) %>%
    summarise(across(all_of(att_ns_cols), mean), .groups = "drop")
  cat(sprintf("  chunk %d/%d: %s daily rows\n", i, length(chunks),
              format(nrow(d), big.mark = ",")))
}
pred_member <- bind_rows(pred_chunks) %>%
  arrange(gadm_key, year, month, scenario, member)
basis_member <- bind_rows(basis_chunks) %>%
  arrange(gadm_key, year, month, scenario, member)

# ----- 10g. Save + summarise -----
write_csv(pred_member, attribution_out_path)
cat("\nSaved attribution predictions to:\n  ", attribution_out_path, "\n")
cat(sprintf("Rows: %d (district x year-month x scenario x member)\n", nrow(pred_member)))

wide <- pred_member %>%
  select(gadm_key, state, district, year, month, scenario, member,
         cdd, predicted_outage_minutes, n_obs_days, obs_outage) %>%
  tidyr::pivot_wider(names_from = scenario,
                     values_from = c(cdd, predicted_outage_minutes)) %>%
  rename(any_of(setNames(
    c(paste0("predicted_outage_minutes_", all_label),
      paste0("predicted_outage_minutes_", nat_label),
      paste0("cdd_", all_label),
      paste0("cdd_", nat_label)),
    c("pred_all", "pred_nat", "cdd_all", "cdd_nat")))) %>%
  mutate(d_cdd    = cdd_all - cdd_nat,
         d_outage = pred_all - pred_nat)

need_wide <- c("pred_all", "pred_nat", "cdd_all", "cdd_nat")
if (!all(need_wide %in% names(wide)))
  stop("pivot_wider did not produce ", paste(setdiff(need_wide, names(wide)), collapse = ", "),
       " -- check that all_label/nat_label match the `leg` values in the CDD file.")

# ---- percentage shares -------------------------------------------------------
# Shares are RATIOS OF SUMS, not means of per-cell ratios:
#     pct = 100 * sum(w * d_outage) / sum(w * pred_all)
# Per-cell ratios are unusable because pred_all is a linear prediction of a
# zero-bounded outcome and can be near zero or negative for low-outage cells,
# which makes individual ratios explode or flip sign. The ratio of sums is also
# what "X% of outage minutes are attributable" should mean.
#
# IMPLEMENTATION NOTE: every aggregate below is built from explicit weighted SUMS
# inside summarise(), then divided in a SEPARATE mutate(). Never compute a share
# in the same summarise() that defines its own numerator -- dplyr evaluates those
# expressions sequentially, so a later one sees the already-collapsed length-1
# value of an earlier one while other columns are still full length, and the
# mismatched lengths silently produce NA.

cat("\n--- predicted outage levels (these set the percentage denominator) ---\n")
print(summary(select(wide, pred_all, pred_nat, d_outage)))
n_bad_denom <- sum(wide$pred_all <= 0, na.rm = TRUE)
cat(sprintf("cell-member predictions with ACT outage <= 0: %d of %d (%.1f%%)\n",
            n_bad_denom, nrow(wide), 100 * n_bad_denom / nrow(wide)))

# --- per district-year-month (ensemble mean over members) ---
cell_summary <- wide %>%
  group_by(gadm_key, state, district, year, month) %>%
  summarise(n_members           = n_distinct(member),
            n_obs_days          = first(n_obs_days),
            obs_outage          = first(obs_outage),
            attributable_cdd    = mean(d_cdd),
            attributable_outage = mean(d_outage),
            pred_outage_all     = mean(pred_all),
            pred_outage_nat     = mean(pred_nat),
            .groups = "drop") %>%
  mutate(pct_of_all      = ifelse(pred_outage_all > 0,
                                  100 * attributable_outage / pred_outage_all, NA_real_),
         pct_of_observed = ifelse(obs_outage > 0,
                                  100 * attributable_outage / obs_outage, NA_real_)) %>%
  arrange(desc(attributable_outage))
write_csv(cell_summary, attribution_sum_path)
cat("Saved cell-level summary to:\n  ", attribution_sum_path, "\n")

# --- per district, weighted by observed days ---
district_summary <- cell_summary %>%
  group_by(state, district) %>%
  summarise(n_cells  = n(),
            obs_days = sum(n_obs_days),
            w_attr   = sum(n_obs_days * attributable_outage),
            w_cdd    = sum(n_obs_days * attributable_cdd),
            w_all    = sum(n_obs_days * pred_outage_all),
            w_nat    = sum(n_obs_days * pred_outage_nat),
            w_obs    = sum(n_obs_days * obs_outage),
            .groups  = "drop") %>%
  mutate(attributable_cdd    = w_cdd  / obs_days,
         attributable_outage = w_attr / obs_days,
         pred_outage_all     = w_all  / obs_days,
         pred_outage_nat     = w_nat  / obs_days,
         obs_outage          = w_obs  / obs_days,
         pct_of_all          = ifelse(w_all > 0, 100 * w_attr / w_all, NA_real_),
         pct_of_observed     = ifelse(w_obs > 0, 100 * w_attr / w_obs, NA_real_)) %>%
  select(state, district, n_cells, obs_days, attributable_cdd, attributable_outage,
         pred_outage_all, pred_outage_nat, obs_outage, pct_of_all, pct_of_observed) %>%
  arrange(desc(pct_of_all))
write_csv(district_summary, attribution_dist_path)
cat("\n--- Attributable outage by district (per observed day) ---\n")
print(district_summary, n = Inf)

cat("\n--- By year-month (observation-weighted across districts) ---\n")
cell_summary %>%
  group_by(year, month) %>%
  summarise(n_districts = n(),
            obs_days    = sum(n_obs_days),
            w_attr      = sum(n_obs_days * attributable_outage),
            w_cdd       = sum(n_obs_days * attributable_cdd),
            w_all       = sum(n_obs_days * pred_outage_all),
            w_obs       = sum(n_obs_days * obs_outage),
            .groups     = "drop") %>%
  mutate(attributable_cdd    = w_cdd  / obs_days,
         attributable_outage = w_attr / obs_days,
         pred_outage_all     = w_all  / obs_days,
         pct_of_all          = ifelse(w_all > 0, 100 * w_attr / w_all, NA_real_),
         pct_of_observed     = ifelse(w_obs > 0, 100 * w_attr / w_obs, NA_real_)) %>%
  select(year, month, n_districts, obs_days, attributable_cdd, attributable_outage,
         pred_outage_all, pct_of_all, pct_of_observed) %>%
  arrange(year, month) %>% print(n = Inf)

# --- ensemble spread: one number per member ---
# Members are paired across legs (same id in both), so the difference is taken
# within member and the CI comes from the spread of those paired differences.
member_delta <- wide %>%
  group_by(member) %>%
  summarise(w       = sum(n_obs_days),
            w_d     = sum(n_obs_days * d_outage),
            w_all   = sum(n_obs_days * pred_all),
            w_obs   = sum(n_obs_days * obs_outage),
            .groups = "drop") %>%
  mutate(d       = w_d / w,
         pct     = ifelse(w_all > 0, 100 * w_d / w_all, NA_real_),
         pct_obs = ifelse(w_obs > 0, 100 * w_d / w_obs, NA_real_))

m  <- mean(member_delta$d)
se <- sd(member_delta$d) / sqrt(nrow(member_delta))
cat(sprintf("\nEnsemble-mean attributable outage: %.3f min per observed district-day\n", m))
cat(sprintf("  95%% CI %.3f to %.3f  (%d paired members, sd %.3f)\n",
            m - 1.96 * se, m + 1.96 * se, nrow(member_delta), sd(member_delta$d)))
cat(sprintf("  member 5th-95th percentile: %.3f to %.3f\n",
            quantile(member_delta$d, 0.05), quantile(member_delta$d, 0.95)))

tot_w    <- sum(cell_summary$n_obs_days)
tot_attr <- sum(cell_summary$n_obs_days * cell_summary$attributable_outage)
tot_all  <- sum(cell_summary$n_obs_days * cell_summary$pred_outage_all)
tot_obs  <- sum(cell_summary$n_obs_days * cell_summary$obs_outage)
overall_pct     <- if (tot_all > 0) 100 * tot_attr / tot_all else NA_real_
overall_pct_obs <- if (tot_obs > 0) 100 * tot_attr / tot_obs else NA_real_

pm  <- mean(member_delta$pct, na.rm = TRUE)
pse <- sd(member_delta$pct, na.rm = TRUE) / sqrt(sum(is.finite(member_delta$pct)))

cat("\n================ HEADLINE ================\n")
if (is.finite(overall_pct)) {
  cat(sprintf("Attributable share of predicted ACT (%s) outage: %.2f%%\n",
              all_label, overall_pct))
  cat(sprintf("  per-member shares: mean %.2f%%, 95%% CI %.2f%% to %.2f%%\n",
              pm, pm - 1.96 * pse, pm + 1.96 * pse))
} else {
  cat(sprintf("Attributable share of predicted ACT (%s) outage: NOT DEFINED\n", all_label))
  cat("  (summed ACT prediction is <= 0 -- see the levels printed above)\n")
}
cat(sprintf("Attributable share of OBSERVED ESMI outage:      %.2f%%\n", overall_pct_obs))
cat(sprintf("  levels: observed %.2f, predicted ACT %.2f, predicted NAT %.2f min/day\n",
            tot_obs / tot_w, tot_all / tot_w, (tot_all - tot_attr) / tot_w))
cat(sprintf("Ensemble-mean attributable CDD: %.3f degree-days/day\n",
            sum(cell_summary$n_obs_days * cell_summary$attributable_cdd) / tot_w))
cat("=========================================\n")
cat("\nAttribution done.\n")

# ================================================================
# 10h. Monte Carlo propagation of CDD-coefficient uncertainty into
#      the district and overall attributable-outage estimates
# ================================================================
# (Numbered as an extension of section 10, not "12" -- it runs before the
# existing section 11 below and depends only on section 10's objects, so this
# keeps the file's top-to-bottom order matching its section numbers.)
# Everything above uses a single point estimate of the CDD spline
# coefficients (att_beta). This is a DIFFERENT, complementary source of
# uncertainty to the ensemble-based CI already computed just above
# (member_delta, pm/pse): that CI reflects spread ACROSS the 524 climate-model
# members, holding att_beta fixed; this section instead holds the ensemble
# fixed (using every member, exactly as already aggregated above) and asks how
# much the attributable-outage estimate would move if the REGRESSION itself
# had come out slightly differently -- i.e. parameter/estimation uncertainty,
# not climate-model uncertainty. Report both, not one instead of the other.
#
# Kept fast and memory-safe (this ensemble has ~14M daily rows) by exploiting
# linearity, as in the state-level demand/shortage script's attribution
# section, extended two steps further here because two things differ there:
#
#   1. Fixed effects are handled differently in this script: fe_baseline was
#      already computed and subtracted out (predictions are formed as
#      X %*% beta + fe_baseline, never through predict(newdata=)), and it is
#      IDENTICAL between ALL and NAT for the same district-year-month cell --
#      this section's own comment above says so ("It is identical in both
#      legs, so it cancels in ALL-NAT"). So attributable_outage needs no FE
#      term at all -- it cancels exactly, not just on average, for every MC
#      draw. FE only has to be recovered for the ACT LEVEL (pred_outage_all),
#      which sets the denominator of the percentage shares, and it's
#      recovered by subtracting the (point-estimate) basis contribution back
#      out of district_summary's already-computed pred_outage_all -- the same
#      trick used for the fixed-effect baseline in the state-level script.
#
#   2. pct_of_all is a RATIO of two beta-dependent quantities (attributable
#      outage over predicted ACT outage), so unlike the state-level script's
#      plain attributable-demand number, it is not itself a linear function of
#      beta. It is still cheap per draw: both the numerator and denominator
#      ARE linear in beta, so each draw needs only two small matrix multiplies
#      and one division -- nothing close to the cost of re-running predict()
#      over the daily data again.
#
# basis_member (built in the chunked loop above) already holds the "mean per
# model day" basis vector for every (district, year, month, scenario, member)
# -- exactly parallel to pred_member's predicted_outage_minutes column, just
# for the 4 spline terms instead of the finished scalar prediction.
MC_N <- 1000
set.seed(20240601)

mc_beta_mat <- MASS::mvrnorm(
  MC_N,
  mu    = as.numeric(att_beta),
  Sigma = vcov(att_fit)[att_ns_cols, att_ns_cols]
)
colnames(mc_beta_mat) <- att_ns_cols

# ---- pivot the mean basis vectors to ALL/NAT columns, one row per cell-member ----
basis_wide <- basis_member %>%
  tidyr::pivot_wider(names_from = scenario, values_from = all_of(att_ns_cols)) %>%
  left_join(
    wide %>% select(gadm_key, state, district, year, month, member, n_obs_days),
    by = c("gadm_key", "state", "district", "year", "month", "member")
  )

all_basis_cols <- paste0(att_ns_cols, "_", all_label)
nat_basis_cols <- paste0(att_ns_cols, "_", nat_label)
if (!all(c(all_basis_cols, nat_basis_cols) %in% names(basis_wide)))
  stop("pivot_wider on basis_member did not produce the expected ALL/NAT columns -- ",
       "check that all_label/nat_label match the `leg` values in the CDD file.")

# ---- mean over members, within each district-year-month cell (mirrors cell_summary) ----
cell_basis <- basis_wide %>%
  group_by(gadm_key, state, district, year, month) %>%
  summarise(across(all_of(c(all_basis_cols, nat_basis_cols)), mean),
            n_obs_days = first(n_obs_days),
            .groups = "drop")

# ---- weighted (by n_obs_days) up to district level (mirrors district_summary) ----
district_basis <- cell_basis %>%
  group_by(state, district) %>%
  summarise(across(all_of(c(all_basis_cols, nat_basis_cols)),
                   ~ sum(.x * n_obs_days) / sum(n_obs_days)),
            obs_days = sum(n_obs_days),
            .groups = "drop") %>%
  arrange(state, district)

# ---- weighted (by n_obs_days) up to YEAR-MONTH, across all districts (mirrors
# B_all_overall/B_diff_overall below, but resolved by calendar time instead of
# collapsed over the whole study period). Same mc_beta_mat draws as everywhere
# else in this section -- a month's interval is that shared coefficient
# uncertainty evaluated at the CDD mixture that month actually saw, not an
# independent per-month estimate, so its width differs across months only
# because of where each month sits on the response curve. ----
month_basis <- cell_basis %>%
  group_by(year, month) %>%
  summarise(across(all_of(c(all_basis_cols, nat_basis_cols)),
                   ~ sum(.x * n_obs_days) / sum(n_obs_days)),
            obs_days = sum(n_obs_days),
            .groups = "drop") %>%
  arrange(year, month)

B_all_month  <- as.matrix(month_basis[all_basis_cols]); colnames(B_all_month) <- att_ns_cols
B_nat_month  <- as.matrix(month_basis[nat_basis_cols]); colnames(B_nat_month) <- att_ns_cols
B_diff_month <- B_all_month - B_nat_month

mc_attr_month <- B_diff_month %*% t(mc_beta_mat)   # n_month x MC_N; FE cancels exactly, as above

month_uncertainty <- data.frame(
  year                 = month_basis$year,
  month                = month_basis$month,
  obs_days             = month_basis$obs_days,
  point_attributable   = as.numeric(B_diff_month %*% as.numeric(att_beta)),
  mc_attributable_mean = rowMeans(mc_attr_month),
  mc_attributable_lo   = apply(mc_attr_month, 1, quantile, probs = 0.025),
  mc_attributable_hi   = apply(mc_attr_month, 1, quantile, probs = 0.975)
)

month_uncertainty_path <- file.path(output_dir, "attribution_uncertainty_by_month.csv")
write_csv(month_uncertainty, month_uncertainty_path)
cat("Saved month-resolved attribution uncertainty (coefficient MC) to:\n  ",
    month_uncertainty_path, "\n")

# ---- mean relative attribution-CI width ------------------------------------
# The direct analog of mean_rel_pi_pct (model-fit section below), but applied
# to the MC-propagated ATTRIBUTION interval itself (the dashed band in the
# outage timeseries figure) rather than the model's own prediction interval.
# One scale-free percentage, directly comparable to the same statistic in the
# demand/shortage scripts. Denominator is each month's own ACT (factual)
# predicted level, observation-weighted across districts the same way as the
# year-month summary printed above -- not a single global mean, so a month
# with a naturally small or large level is not over- or under-weighted.
month_pred_all <- cell_summary %>%
  group_by(year, month) %>%
  summarise(w_all = sum(n_obs_days * pred_outage_all),
            obs_days = sum(n_obs_days),
            .groups = "drop") %>%
  mutate(pred_outage_all = w_all / obs_days) %>%
  select(year, month, pred_outage_all) %>%
  arrange(year, month)

stopifnot(identical(month_pred_all$year, month_uncertainty$year),
          identical(month_pred_all$month, month_uncertainty$month))

bad_level <- month_pred_all$pred_outage_all <= 0
if (any(bad_level)) {
  cat(sprintf("  NOTE: %d of %d months have a non-positive ACT predicted level and are excluded from the ratio below\n", sum(bad_level), length(bad_level)))
}
rel_width <- ifelse(bad_level, NA_real_,
                    (month_uncertainty$mc_attributable_hi - month_uncertainty$mc_attributable_lo) /
                      month_pred_all$pred_outage_all)
mean_rel_attrib_ci_pct <- 100 * mean(rel_width, na.rm = TRUE)
cat(sprintf("Mean relative width of the attribution MC 95%% CI: %.1f%% (mean over %d months of (mc_hi - mc_lo) / ACT predicted level)\n", mean_rel_attrib_ci_pct, sum(!bad_level)))

B_all_district  <- as.matrix(district_basis[all_basis_cols]); colnames(B_all_district) <- att_ns_cols
B_nat_district  <- as.matrix(district_basis[nat_basis_cols]); colnames(B_nat_district) <- att_ns_cols
B_diff_district <- B_all_district - B_nat_district

# Recover the district-level FE contribution to the ACT level, and carry the
# point-estimate numbers along for comparison, by joining onto district_basis's
# row order (NOT district_summary's -- it's sorted differently, by
# desc(pct_of_all), and using it positionally would silently misalign rows).
district_join <- district_basis %>%
  left_join(
    district_summary %>%
      select(state, district, pred_outage_all, obs_outage,
             point_attributable = attributable_outage,
             point_pct_of_all = pct_of_all,
             point_pct_of_observed = pct_of_observed),
    by = c("state", "district")
  )
fe_all_district <- district_join$pred_outage_all -
  as.numeric(B_all_district %*% as.numeric(att_beta))

# ---- overall (headline), weighted the same way as tot_w/tot_attr/tot_all above ----
B_all_overall <- colSums(as.matrix(district_basis[all_basis_cols]) * district_basis$obs_days) / sum(district_basis$obs_days)
B_nat_overall <- colSums(as.matrix(district_basis[nat_basis_cols]) * district_basis$obs_days) / sum(district_basis$obs_days)
names(B_all_overall) <- names(B_nat_overall) <- att_ns_cols
B_diff_overall     <- B_all_overall - B_nat_overall
fe_all_overall     <- (tot_all / tot_w) - as.numeric(B_all_overall %*% as.numeric(att_beta))
obs_outage_overall <- tot_obs / tot_w

# ---- Monte Carlo: attributable outage and pct-of-ACT for every draw ----
# District level: n_district x MC_N matrices. R recycles a length-n_district
# vector down each column of an n_district x MC_N matrix correctly (row-
# aligned), so fe_all_district / district_basis$obs_outage add and divide
# elementwise without needing sweep().
mc_attr_district    <- B_diff_district %*% t(mc_beta_mat)                    # FE cancels exactly
mc_all_district     <- B_all_district %*% t(mc_beta_mat) + fe_all_district   # FE recovered above
mc_pct_all_district <- 100 * mc_attr_district / mc_all_district
# obs_outage doesn't depend on beta at all, so this share is just a linear
# rescaling of mc_attr_district -- no separate "denominator draw" needed.
mc_pct_obs_district <- 100 * mc_attr_district / district_join$obs_outage

district_uncertainty <- data.frame(
  state                   = district_join$state,
  district                = district_join$district,
  point_attributable      = district_join$point_attributable,
  point_pct_of_all        = district_join$point_pct_of_all,
  point_pct_of_observed   = district_join$point_pct_of_observed,
  mc_attributable_mean    = rowMeans(mc_attr_district),
  mc_attributable_lo      = apply(mc_attr_district, 1, quantile, probs = 0.025),
  mc_attributable_hi      = apply(mc_attr_district, 1, quantile, probs = 0.975),
  mc_pct_of_all_mean      = rowMeans(mc_pct_all_district, na.rm = TRUE),
  mc_pct_of_all_lo        = apply(mc_pct_all_district, 1, quantile, probs = 0.025, na.rm = TRUE),
  mc_pct_of_all_hi        = apply(mc_pct_all_district, 1, quantile, probs = 0.975, na.rm = TRUE),
  mc_pct_of_observed_mean = rowMeans(mc_pct_obs_district, na.rm = TRUE),
  mc_pct_of_observed_lo   = apply(mc_pct_obs_district, 1, quantile, probs = 0.025, na.rm = TRUE),
  mc_pct_of_observed_hi   = apply(mc_pct_obs_district, 1, quantile, probs = 0.975, na.rm = TRUE)
)

district_uncertainty_path <- file.path(output_dir, "attribution_uncertainty_by_district.csv")
write_csv(district_uncertainty, district_uncertainty_path)
cat("\nSaved district-level attribution uncertainty (coefficient MC) to:\n  ",
    district_uncertainty_path, "\n")

# Overall headline: MC_N-length vectors
mc_attr_overall    <- as.numeric(B_diff_overall %*% t(mc_beta_mat))
mc_all_overall     <- as.numeric(B_all_overall %*% t(mc_beta_mat)) + fe_all_overall
mc_pct_all_overall <- 100 * mc_attr_overall / mc_all_overall
mc_pct_obs_overall <- 100 * mc_attr_overall / obs_outage_overall

# Nationwide headline as its own CSV row -- so a plotting script can pick up
# the properly-weighted overall MC interval directly, rather than trying to
# reconstruct it by averaging the per-district lo/hi in attribution_uncertainty_
# by_district.csv (which would NOT be the same number: quantiles of a weighted
# sum are not a weighted sum of quantiles).
overall_uncertainty <- data.frame(
  point_attributable      = as.numeric(B_diff_overall %*% as.numeric(att_beta)),
  mc_attributable_mean    = mean(mc_attr_overall),
  mc_attributable_lo      = as.numeric(quantile(mc_attr_overall, 0.025)),
  mc_attributable_hi      = as.numeric(quantile(mc_attr_overall, 0.975)),
  mc_pct_of_all_mean      = mean(mc_pct_all_overall, na.rm = TRUE),
  mc_pct_of_all_lo        = as.numeric(quantile(mc_pct_all_overall, 0.025, na.rm = TRUE)),
  mc_pct_of_all_hi        = as.numeric(quantile(mc_pct_all_overall, 0.975, na.rm = TRUE)),
  mc_pct_of_observed_mean = mean(mc_pct_obs_overall, na.rm = TRUE),
  mc_pct_of_observed_lo   = as.numeric(quantile(mc_pct_obs_overall, 0.025, na.rm = TRUE)),
  mc_pct_of_observed_hi   = as.numeric(quantile(mc_pct_obs_overall, 0.975, na.rm = TRUE))
)
overall_uncertainty_path <- file.path(output_dir, "attribution_uncertainty_overall.csv")
write_csv(overall_uncertainty, overall_uncertainty_path)
cat("Saved nationwide attribution uncertainty (coefficient MC) to:\n  ",
    overall_uncertainty_path, "\n")

cat("\n---- Coefficient-uncertainty (Monte Carlo) headline, alongside the ensemble CI above ----\n")
cat(sprintf("Attributable outage: mean %.3f min/day, 95%% MC interval %.3f to %.3f\n",
            mean(mc_attr_overall), quantile(mc_attr_overall, 0.025), quantile(mc_attr_overall, 0.975)))
cat(sprintf("Pct of predicted ACT outage: mean %.2f%%, 95%% MC interval %.2f%% to %.2f%%\n",
            mean(mc_pct_all_overall, na.rm = TRUE),
            quantile(mc_pct_all_overall, 0.025, na.rm = TRUE),
            quantile(mc_pct_all_overall, 0.975, na.rm = TRUE)))
cat(sprintf("Pct of observed ESMI outage: mean %.2f%%, 95%% MC interval %.2f%% to %.2f%%\n",
            mean(mc_pct_obs_overall), quantile(mc_pct_obs_overall, 0.025), quantile(mc_pct_obs_overall, 0.975)))

# ================================================================
# 11. Export the fitted response function for plotting
# ================================================================
# Run AFTER section 10 -- it reuses att_fit, att_basis, att_beta, att_ns_cols,
# cells and used from there.
#
# Writes the CURVE ITSELF rather than the basis matrix. The old figure exported
# the basis on a grid so Python could rebuild the prediction; exporting the
# finished curve removes the last place an R/Python basis mismatch could hide,
# and keeps the confidence band computed where the vcov actually lives.
#
# The band is spline-coefficient uncertainty only:
#     se(x) = sqrt( B(x) V B(x)' ),  V = vcov(att_fit)[ns, ns]
# It excludes uncertainty in the fixed-effect baseline, exactly as the caption
# on the earlier state-level figure said. That is a deliberate choice -- the
# baseline shifts the whole curve up or down without changing its shape, and the
# shape is what the figure is about.
#
# BASELINE. With FE specified as `| district_month + district_year`, feols
# absorbs the intercept into those group levels instead of reporting one
# scalar -- there is no "(Intercept)" in coef(att_fit) to read off. Like the
# demand/shortage script, "baseline" here is a CONSTRUCTED number, not a model
# parameter: an observation-weighted average across many fixed-effect levels,
# not a single value the model hands you.
#
# It uses a DIFFERENT averaging scheme than the demand/shortage script's,
# though, and deliberately so. The curve is drawn at the NATION-WIDE average
# baseline: the mean of the district-year-month fixed effects (from `cells`,
# reconstructed via fixef(att_fit) in section 10d), restricted to and weighted
# by the cells that were actually attributed (`used`) -- i.e. the same
# district-year-month universe the headline ALL-NAT numbers are computed
# over, rather than a raw average across every row of data_panel. So the y
# axis reads as "predicted outage minutes for a district-day at the average
# non-climate baseline of the attributed cells".
# ----------------------------------------------------------------

curve_out_path <- file.path(output_dir,
                            paste0("response_curve_", attribution_model_name, ".csv"))
knots_out_path <- file.path(output_dir,
                            paste0("response_curve_knots_", attribution_model_name, ".csv"))

stopifnot(exists("att_fit"), exists("att_basis"), exists("att_beta"),
          exists("cells"), exists("used"))

# ---- nation-wide baseline: obs-day-weighted average FE level over the
# attributed cells (constructed, not read from the model -- see note above) ----
base_cells <- cells %>%
  semi_join(distinct(used, gadm_key, year, month), by = c("gadm_key", "year", "month"))
baseline <- weighted.mean(base_cells$fe_baseline, base_cells$n_obs_days)
cat(sprintf("Nation-wide FE baseline: %.2f minutes (over %d attributed cells, %d obs-days)\n",
            baseline, nrow(base_cells), sum(base_cells$n_obs_days)))

# ---- grid across the fitted support ----
x <- data_panel[[panel_predictor]]
grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 400)
B <- predict(att_basis, newx = grid)          # stored knots -> identical basis
V <- vcov(att_fit)[att_ns_cols, att_ns_cols]

CI_Z <- 1.96

fit <- as.numeric(B %*% att_beta) + baseline
se  <- sqrt(rowSums((B %*% V) * B))           # diag(B V B') without forming it

# ---- PREDICTION band ----
# Adds the model's residual variance to the line's se, exactly as in the
# demand/shortage script: pred_se(x) = sqrt(se(x)^2 + resid_sigma^2). `lo`/`hi`
# above answer "how precisely do we know the average curve"; `pred_lo`/
# `pred_hi` answer "how far can one actual district-day land from it" -- the
# band that shows how much variation CDD leaves unexplained even where the
# average curve itself is estimated tightly.
att_sigma <- sigma(att_fit)
pred_se   <- sqrt(se^2 + att_sigma^2)

# ---- single-number summaries, scale-free so they sit next to the demand/
# shortage script's numbers even though outage minutes is on a different
# scale (CV(RMSE), and mean relative PI width evaluated at the panel's actual
# observed CDD values rather than the 400-point plotting grid, for the same
# reason given in the demand/shortage script: the grid oversamples the
# extremes where few real district-days sit) ----
pi_fitted <- predict(att_fit, sample = "original")
pi_keep   <- !is.na(pi_fitted)
mean_y    <- mean(data_panel[[dependent_variable]][pi_keep])

cv_rmse_pct <- 100 * att_sigma / mean_y

obs_B           <- as.matrix(data_panel[pi_keep, att_ns_cols, drop = FALSE])
obs_se_line     <- sqrt(pmax(rowSums((obs_B %*% V) * obs_B), 0))
obs_pred_se     <- sqrt(obs_se_line^2 + att_sigma^2)
mean_rel_pi_pct <- 100 * mean(2 * CI_Z * obs_pred_se) / mean_y

# ---- cluster (block) bootstrap band ----
# Empirical alternative to the two analytic bands above: resample whole
# DISTRICTS with replacement (not individual rows -- a district's observations
# are correlated over time, so resampling rows would understate uncertainty
# the same way un-clustered SEs would, exactly the logic behind cluster_formula
# itself), refit att_formula on each resample, and take percentiles of the
# resulting curve SHAPE. This doesn't assume asymptotic normality the way the
# two bands above do, and it's unaffected by the non-positive-definite CGM
# issue, since it never builds that matrix at all.
#
# ONE DELIBERATE DIFFERENCE from the demand/shortage script's version: there,
# the bootstrap baseline is re-derived per replicate and so the band also
# picks up fixed-effect/baseline resampling uncertainty. Here the band is
# added to the EXISTING nation-wide baseline (held fixed at its point
# estimate) instead, because that baseline comes from the ~14M-row ensemble
# attribution pipeline above (sections 10c-10h) -- refitting THAT per
# bootstrap replicate is infeasible. So this band captures resampled
# spline-coefficient (shape) uncertainty under district-cluster resampling
# only, same as bands 1 and 2 above, just without assuming normality.
#
# A drawn district is relabelled with a unique suffix so a district drawn
# twice becomes two distinct synthetic units (otherwise feols would silently
# collapse the duplicate rows onto the same FE levels instead of treating them
# as two independent draws). district_month/district_year are rebuilt from the
# relabelled district so their FE levels stay consistent with it -- this
# assumes district_month/district_year follow the same paste(district,
# month/year, sep = "_") convention as state_month/state_year did in the
# demand/shortage script; double check that against how they're actually
# built upstream of this script if the bootstrap errors out or looks off.
BOOT_REPS <- 500   # drop to ~100 for a quick first pass, then rerun at 500+
set.seed(20240601)

boot_districts <- unique(data_panel$district)
boot_curves    <- matrix(NA_real_, nrow = length(grid), ncol = BOOT_REPS)
boot_ok        <- 0

for (b in seq_len(BOOT_REPS)) {
  drawn <- sample(boot_districts, length(boot_districts), replace = TRUE)
  
  boot_panel <- purrr::map2_dfr(drawn, seq_along(drawn), function(dd, i) {
    data_panel %>%
      filter(district == dd) %>%
      mutate(
        district       = paste0(dd, "__", i),
        district_month = paste(district, month, sep = "_"),
        district_year  = paste(district, year, sep = "_")
      )
  })
  
  boot_fit <- tryCatch(
    fixest::feols(att_formula, data = boot_panel, nthreads = 1),
    error = function(e) NULL
  )
  if (is.null(boot_fit)) next
  
  boot_beta <- tryCatch(coef(boot_fit)[att_ns_cols], error = function(e) NULL)
  if (is.null(boot_beta) || anyNA(boot_beta)) next
  
  # shape-only deviation from the point estimate, so it can be added to the
  # existing `fit` (nation-wide baseline) below without re-deriving a baseline
  boot_curves[, b] <- as.numeric(B %*% boot_beta) - as.numeric(B %*% att_beta)
  boot_ok <- boot_ok + 1
}

cat(sprintf("Cluster bootstrap: %d/%d replicates fit successfully\n", boot_ok, BOOT_REPS))

boot_shape_lo <- apply(boot_curves, 1, quantile, probs = 0.025, na.rm = TRUE)
boot_shape_hi <- apply(boot_curves, 1, quantile, probs = 0.975, na.rm = TRUE)
boot_lo <- fit + boot_shape_lo
boot_hi <- fit + boot_shape_hi

curve <- tibble::tibble(
  cdd             = grid,
  fit             = fit,
  se              = se,
  lo              = fit - CI_Z * se,
  hi              = fit + CI_Z * se,
  pred_lo         = fit - CI_Z * pred_se,
  pred_hi         = fit + CI_Z * pred_se,
  boot_lo         = boot_lo,
  boot_hi         = boot_hi,
  boot_reps_ok    = boot_ok,
  resid_sigma     = att_sigma,
  cv_rmse_pct     = cv_rmse_pct,
  mean_rel_pi_pct = mean_rel_pi_pct,
  baseline        = baseline,
  model           = attribution_model_name,
  predictor       = panel_predictor
)
readr::write_csv(curve, curve_out_path)

cat(sprintf("Residual sigma: %.2f | mean CI half-width: %.2f | mean PI half-width: %.2f | mean bootstrap half-width: %.2f\n",
            att_sigma, CI_Z * mean(se), CI_Z * mean(pred_se), mean((boot_hi - boot_lo) / 2)))
cat(sprintf("CV(RMSE): %.1f%%  |  mean relative PI width (obs-weighted): %.1f%%  [%s]\n",
            cv_rmse_pct, mean_rel_pi_pct, dependent_variable))

# knots, so the figure can mark where the spline bends and where the linear
# tails take over
kn <- attr(att_basis, "knots")
bk <- attr(att_basis, "Boundary.knots")
readr::write_csv(
  tibble::tibble(kind = c(rep("interior", length(kn)), rep("boundary", length(bk))),
                 value = c(as.numeric(kn), as.numeric(bk))),
  knots_out_path)

cat(sprintf("Curve: %d points over CDD %.2f-%.2f; predicted outage %.1f to %.1f min\n",
            nrow(curve), min(grid), max(grid), min(fit), max(fit)))
cat(sprintf("Interior knots: %s | boundary: %s\n",
            paste(sprintf("%.2f", kn), collapse = ", "),
            paste(sprintf("%.2f", bk), collapse = ", ")))
cat("Saved:\n  ", curve_out_path, "\n  ", knots_out_path, "\n")

# ---- observed CDD distribution, for the density panel under the curve ----
# The panel's own (ERA5-based) CDD over the attributed cells, day-weighted.
obs_cdd <- panel_keyed %>%
  semi_join(distinct(used, gadm_key, year, month), by = c("gadm_key", "year", "month")) %>%
  transmute(cdd = .data[[panel_predictor]]) %>%
  filter(is.finite(cdd))
readr::write_csv(obs_cdd, file.path(output_dir, "response_curve_observed_cdd.csv"))
cat(sprintf("Observed CDD for the density panel: %d district-days, mean %.2f\n",
            nrow(obs_cdd), mean(obs_cdd$cdd)))

# =============================================================================
# 12. K-fold cross-validated (out-of-sample) R^2
# =============================================================================
# Same idea, same caveats, as the demand/shortage script's version. Every R^2
# reported elsewhere for this model (r2(att_fit), the CV(RMSE)/relative-PI
# numbers above) is IN-SAMPLE: fit and evaluated on the same rows. With many
# district_month and district_year fixed effects soaking up degrees of
# freedom, in-sample R^2 can look better than the model's real predictive
# power warrants. This reruns the fit K times, holding out a different 1/K of
# district-days each time and scoring predictions the fit never saw.
#
# Folds are assigned by INDIVIDUAL ROW, not by whole district or year (unlike
# the cluster bootstrap above) -- each district_month/district_year FE level
# has many rows behind it, so a random ~10% holdout leaves nearly every FE
# level identified in training. Held-out rows whose FE level is nonetheless
# entirely absent from training (rare, but possible for a sparsely populated
# level) get NA from predict() and are dropped from the score -- the count
# below reports how many that was.
# =============================================================================
CV_K <- 10
set.seed(20240601)

cv_ids  <- which(pi_keep)              # same estimation sample as the response curve
cv_fold <- sample(rep_len(seq_len(CV_K), length(cv_ids)))
cv_pred <- rep(NA_real_, nrow(data_panel))

for (k in seq_len(CV_K)) {
  test_idx <- cv_ids[cv_fold == k]
  
  fit_k <- tryCatch(
    fixest::feols(att_formula, data = data_panel[-test_idx, ], nthreads = 1),
    error = function(e) NULL
  )
  if (is.null(fit_k)) next
  
  pred_k <- tryCatch(
    as.numeric(predict(fit_k, newdata = data_panel[test_idx, ])),
    error = function(e) rep(NA_real_, length(test_idx))
  )
  cv_pred[test_idx] <- pred_k
}

cv_scored <- !is.na(cv_pred) & pi_keep
y_actual  <- data_panel[[dependent_variable]]

oos_r2 <- 1 - sum((y_actual[cv_scored] - cv_pred[cv_scored])^2) /
  sum((y_actual[cv_scored] - mean(y_actual[cv_scored]))^2)

cat(sprintf("\n%d-fold cross-validated R^2: %.4f  (in-sample R^2 was %.4f)  [%d/%d obs scored, %s]\n",
            CV_K, oos_r2, r2(att_fit, type = "cor2"),
            sum(cv_scored), sum(pi_keep), dependent_variable))

cv_out <- data.frame(
  model_name         = attribution_model_name,
  dependent_variable = dependent_variable,
  cv_k               = CV_K,
  in_sample_r2       = r2(att_fit, type = "cor2"),
  oos_r2             = oos_r2,
  n_scored           = sum(cv_scored),
  n_total            = sum(pi_keep)
)
cv_path <- file.path(output_dir, sprintf("cv_r2_%s.csv", attribution_model_name))
write_csv(cv_out, cv_path)
cat("Saved cross-validated R^2 to:\n  ", cv_path, "\n")