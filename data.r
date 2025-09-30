source("setup.r")
source(paste0(redir, "data.r"))

parts <- c("wake_min", "sol_min", "light_min", "stage3_min", "rem_min", "waso_min")

clr <- complr(
  data = dpsg,
  parts = parts,
  # idvar = "record_id",
  total = 1440
)
colnames(dpsg)
nrow(dpsg)
nrow(dpsg[!is.na(white)])
nrow(dpsg[!is.na(sex)])
nrow(dpsg[!is.na(working)])
nrow(dpsg[!is.na(coupled)])

dpsg[, age45 := ifelse(age >= 45, 1, 0)]
table(dpsg[isi > 7]$age45)
table(dpsg[isi > 7]$servicetype)

# descriptive table
tmp <- dpsg
for (v in names(tmp)) {
 vals <- unique(na.omit(tmp[[v]]))
 if (length(vals) == 2 && all(sort(vals) == c(0, 1))) {
  tmp[[v]] <- factor(tmp[[v]], levels = c(0, 1))
 }
}

egltable(c(
 "age", "sex", "bmi",
 "working", "white", "coupled",
 "perHrAHSleep",
 "isi",
 "wake_min", "sol_min", "light_min", "stage3_min", "rem_min", 
 "waso_min", "tst_min", "sleep_min", "spt_min", "bed_min", "bedrest_min"
), g = "isig", data = tmp)
