source("setup.r")
source(paste0(redir, "data.r"))

parts <- c("wake_min", "sol_min", "light_min", "stage3_min", "rem_min", "waso_min")

# Filter to adults only
dpsg <- dpsg[age >= 18]

dpsg[, ethncg := NA]
dpsg[, ethncg := ifelse(ethnicity == "Caucasian", "White", ethncg)]
dpsg[, ethncg := ifelse(ethnicity %in% c("Northeast Asian", "South Asian", "Southeast Asian"), "Asian (Northeast / Southeast / South)", ethncg)]
dpsg[, ethncg := ifelse(ethnicity %in% c("Indigenous/Pacific Islander"), "Indigenous / Pacific", ethncg)]
dpsg[, ethncg := ifelse(ethnicity %in% c("Central/South American"), "Latinx / Hispanic", ethncg)]
dpsg[, ethncg := ifelse(ethnicity %in% c("Middle Eastern/North African"), "Middle Eastern / North African", ethncg)]
dpsg[, ethncg := ifelse(ethnicity %in% c("Mixed"), "Mixed / Multiracial", ethncg)]
dpsg[, ethncg := ifelse(ethnicity %in% c("Sub-Saharan African"), "Black / African", ethncg)]

clr <- complr(
  data = dpsg[isig == 1],
  parts = parts,
  # idvar = "record_id",
  total = 1440
)

colnames(clr$dataout)
clr$dataout[, age45 := ifelse(age >= 45, 1, 0)]
clr$dataout[, so_min := 1440 - twake_min]
clr$dataout[, tst_min := tlight_min + tstage3_min + trem_min]
clr$dataout[, se := (tst_min / so_min) * 100]

so_min_levels <- c("≤P10", "P25–P75", "≥P90")
so_min_percentile <- quantile(
  clr$dataout$so_min,
  probs = c(0.1, 0.25, 0.75, 0.9),
  na.rm = TRUE
)

q10 <- so_min_percentile[1]
q25 <- so_min_percentile[2]
q75 <- so_min_percentile[3]
q90 <- so_min_percentile[4]

clr$dataout[, so_min_p := fifelse(
  so_min <= q10,
  so_min_levels[1],
  fifelse(
    so_min >= q90,
    so_min_levels[3],
    fifelse(between(so_min, q25, q75), so_min_levels[2], NA_character_)
  )
)]

clr$dataout[, so_min_p := factor(so_min_p, levels = so_min_levels, ordered = TRUE)]

# descriptive summary
vars <- c(
 "age", 
 
 "sex", "ethncg", "coupled", "educ", "working", "income3cat",  "antidep", "labpsg",
 "bmi", "isi",
 
#  "perHrAHSleep",
 "so_min", "tst_min", "se",
 "twake_min", "tsol_min", "tlight_min", "tstage3_min", "trem_min", 
 "twaso_min"
)

labels <- c(
 "Age (years)",
  "Sex", "Ethnicity", "Married or de facto", "Education", "Currently working", 
  "Income", "Antidepressant use", "In-lab PSG",
  "Body mass index (kg/m²)", "Insomnia Severity Index",

  "Sleep opportunity (min)", "Total sleep time (min)", "Sleep efficiency (%)",
  "Daytime wake (min)", "Sleep onset latency (min)", 
  "Sleep stages N1+N2 (min)", "Sleep stage N3 (min)", "Sleep stage REM (min)", "Wake after sleep onset (min)"
)

setnames(clr$dataout, vars, labels)

# descriptives for overall sample
summary <- rbindlist(lapply(labels, function(var) {
  values <- clr$dataout[[var]]
  n_total <- nrow(clr$dataout)
  n_valid <- sum(!is.na(values))

  # Check if variable is continuous (numeric with many unique values)
  is_continuous <- is.numeric(values) && length(unique(na.omit(values))) > 2
  
  if (is_continuous) {
    # For continuous variables: show mean, sd, range
    data.table(
      variable = var,
      type = "continuous",
      level = "",
      summary = paste0(round(mean(values, na.rm = TRUE), 2), " (", 
                       round(sd(values, na.rm = TRUE), 2), ")"),
      range = paste0(round(min(values, na.rm = TRUE), 0), " - ", 
                     round(max(values, na.rm = TRUE), 0)),
      n = n_valid
    )
  } else {
    # For categorical/binary variables: show n and % for each level
    level_counts <- table(clr$dataout[[var]], useNA = "ifany")
    rbindlist(lapply(names(level_counts)[!is.na(names(level_counts))], function(level) {
      count <- as.numeric(level_counts[level])
      data.table(
        variable = var,
        type = "categorical",
        level = as.character(level),
        summary = paste0(count, " (", round(100 * count / n_total, 1), "%)"),
        range = "",
        n = n_valid

      )
    }))
  }
}), use.names = TRUE, fill = TRUE)

summary[, variable := factor(variable, levels = labels, ordered = TRUE)]
setorder(summary, variable, level)

print(summary[, .(variable, level, summary, range, n)])
saveRDS(summary, file = file.path(out, "desc_all.rds"))

# descriptives by quantiles of so_min

summary_by_so_min_p <- rbindlist(lapply(labels, function(var) {
  values <- clr$dataout[[var]]
  n_total <- nrow(clr$dataout)
  n_valid <- sum(!is.na(values))

  # Check if variable is continuous (numeric with many unique values)
  is_continuous <- is.numeric(values) && length(unique(na.omit(values))) > 2
  
  if (is_continuous) {
    # For continuous variables: show mean, sd, range by so_min_p
    clr$dataout[!is.na(so_min_p), .(
      summary = paste0(round(mean(get(var), na.rm = TRUE), 2), " (", 
                       round(sd(get(var), na.rm = TRUE), 2), ")"),
      range = paste0(round(min(get(var), na.rm = TRUE), 0), " - ", 
                     round(max(get(var), na.rm = TRUE), 0)),
      n = sum(!is.na(get(var)))
    ), by = so_min_p][
      , `:=`(
        variable = var,
        sleep_opportunity_percentile = as.character(so_min_p),
        level = "",
        type = "continuous"
      )
    ]
  } else {
    # For categorical/binary variables: show n and % for each level by so_min_p
    clr$dataout[!is.na(so_min_p), .N, by = .(so_min_p, category = get(var))][
      , `:=`(
        variable = var,
        sleep_opportunity_percentile = as.character(so_min_p),
        level = as.character(category),
        summary = paste0(N, " (", round(100 * N / n_total, 1), "%)"),
        range = "",
        type = "categorical"
      )
    ][
      , n := sum(N), by = so_min_p
    ][
      , `:=`(N = NULL, category = NULL)
    ]
  }
}), use.names = TRUE, fill = TRUE)

summary_by_so_min_p[, `:=`(
  variable = factor(variable, levels = labels, ordered = TRUE),
  sleep_opportunity_percentile = factor(sleep_opportunity_percentile, levels = so_min_levels, ordered = TRUE)
)]
setorder(summary_by_so_min_p, variable, sleep_opportunity_percentile, level)

print(summary_by_so_min_p[, .(variable, sleep_opportunity_percentile, level, summary, range, n)])
saveRDS(summary_by_so_min_p, file = file.path(out, "desc_by_so_min_p.rds"))
