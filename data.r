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

# descriptive table
tmp <- clr$dataout
for (v in names(tmp)) {
 vals <- unique(na.omit(tmp[[v]]))
 if (length(vals) == 2 && all(sort(vals) == c(0, 1))) {
  tmp[[v]] <- factor(tmp[[v]], levels = c(0, 1))
 }
}

# ranges for selected continuous variables
vars_for_range <- c(
 "age", 
#  "sex", 
 "bmi",
#  "working", "white", "coupled", "income",
#  "labpsg",
 "perHrAHSleep",
 "isi",
 "twake_min", "tsol_min", "tlight_min", "tstage3_min", "trem_min", 
 "twaso_min", "tst_min", "spt_min", "so_min"
)

range_summary <- rbindlist(lapply(vars_for_range, function(var) {
  values <- clr$dataout[[var]]

  data.table(
    variable = var,
    min = min(values, na.rm = TRUE),
    max = max(values, na.rm = TRUE),
    completeness = nrow(clr$dataout[!is.na(get(var))])
  )
}))

setorder(range_summary, variable)
range_summary[, range := paste0(round(min, digits = 2), " - ", round(max, digits = 2))]
range_summary[, completeness := completeness]
print(range_summary[, .(variable, range, completeness)])
nrow(clr$dataout[!is.na(coupled)])

egltable(c(
 "age", "sex", "bmi",
 "working", "white", "coupled", "labpsg",
 "perHrAHSleep",
 "isi",
 "twake_min", "tsol_min", "tlight_min", "tstage3_min", "trem_min", 
 "twaso_min", "tst_min", "spt_min", "so_min"
), data = tmp)

# plot histogram of so_min by tib_group

# png(file.path(out, "hist-bedrest.png"), width = 6, height = 4, units = "in", res = 300)
# hist(dpsg$so_min, breaks = 100, main = "SleepPeriod (minutes)", xlab = "Minutes", ylab = "Frequency")
# dev.off()

# png(file.path(out, "hist-bedrest-by-labpsg.png"), width = 8, height = 4, units = "in", res = 300)
# labpsgs <- unique(na.omit(dpsg$labpsg))
# par(mfrow = c(1, length(labpsgs)))
# for (stype in labpsgs) {
#   label <- ifelse(stype == 1, "in-lab", "in-home")
#   hist(
#     dpsg[labpsg == stype]$so_min,
#     breaks = 100,
#     main = paste("SleepPeriod:", label),
#     xlab = "Minutes",
#     ylab = "Frequency"
#   )
# }
# dev.off()

# png(file.path(out, "hist-sol-by-labpsg.png"), width = 8, height = 4, units = "in", res = 300)
# labpsgs <- unique(na.omit(dpsg$labpsg))
# par(mfrow = c(1, length(labpsgs)))
# for (stype in labpsgs) {
#   label <- ifelse(stype == 1, "in-lab", "in-home")
#   hist(
#     dpsg[labpsg == stype]$sol_min,
#     breaks = 100,
#     main = paste("SOL:", label),
#     xlab = "Minutes",
#     ylab = "Frequency"
#   )
# }
# dev.off()

# data check export
# d_sub <- d[tp == "T1" & !is.na(studyDate), .(
#   reference, patient_type,
#   redcap_record_id = record_id, 
#   labpsg, studyDate, SleepPeriod, totalSleepPeriod, 
#   lightOut, lightOn, sleepStart, sleepEnd,

#   sleepLatency, WASO, TST, stage1Time, stage2Time, stage3Time, REMTime, awakeTime,
#   timeAvailSleep
# )]
# write.csv(d_sub, file.path(out, "psg_data_check.csv"), row.names = FALSE)
