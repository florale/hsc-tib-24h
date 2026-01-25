source("setup.r")
source(paste0(redir, "data.r"))

parts <- c("wake_min", "sol_min", "light_min", "stage3_min", "rem_min", "waso_min")

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

vars <- c(
 "age", 
 "bmi",
 "sex", "working", "white", "coupled", "income3cat", "edu",
 "labpsg",
 "perHrAHSleep",
 "isi",
 "twake_min", "tsol_min", "tlight_min", "tstage3_min", "trem_min", 
 "twaso_min", "tst_min", "so_min"
)

summary <- rbindlist(lapply(vars, function(var) {
  values <- clr$dataout[[var]]
  n_total <- nrow(clr$dataout)
  n_valid <- sum(!is.na(values))

  # Check if variable is continuous (numeric with many unique values)
  is_continuous <- is.numeric(values) && length(unique(na.omit(values))) > 10
  
  if (is_continuous) {
    # For continuous variables: show mean, sd, range
    data.table(
      variable = var,
      type = "continuous",
      level = "",
      summary = paste0(round(mean(values, na.rm = TRUE), 1), " (", 
                       round(sd(values, na.rm = TRUE), 1), ")"),
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
}))

print(summary[, .(variable, level, summary, range, n)])
