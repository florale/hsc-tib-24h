source("setup.r")
parts <- c("wake_min", "sol_min", "light_min", "stage3_min", "rem_min", "waso_min")

clr <- complr(
  data = dpsg,
  parts = parts,
  # idvar = "record_id",
  total = 1440
)

m_tib_isig <- brmcoda(clr,
  mvbind(z1_1, z2_1, z3_1, z4_1, z5_1) ~
    s(bedrest_min, by = isig) +
    s(age) + female + bmi + white + working + s(perHrAHSleep) + servicetype,
  iter = 4000, chains = 6, cores = 6, seed = 123, warmup = 1000,
  backend = "cmdstanr"
)
summary(m_tib_isig)
saveRDS(m_tib_isig, file.path(out, "m_tib_isig.rds"))

m_tib_isig <- readRDS(file.path(out, "m_tib_isig.rds"))

## prediction -----------------------
quantile(model.frame(m_tib_isig)$bedrest_min, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

d_tibq5_isig <- emmeans::ref_grid(m_tib_isig$model,
  at = list(
    bedrest_min = quantile(model.frame(m_tib_isig)$bedrest_min, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  )
)@grid
d_tibq5_isig <- as.data.table(d_tibq5_isig)
d_tibq5_isig <- d_tibq5_isig[!duplicated(d_tibq5_isig[, .(isig, bedrest_min, age, female, white, perHrAHSleep)]), ]

pred_tibq5_isig <- fitted(m_tib_isig, newdata = d_tibq5_isig, scale = "response", re_formula = NA, summary = FALSE)
pred_tibq5_isig <- apply(pred_tibq5_isig, c(1), function(x) cbind(d_tibq5_isig, x))

pred_tibq5_isig_draws <- lapply(pred_tibq5_isig, function(d) {
  d <- as.data.table(d)

  setnames(d, paste0("t", parts), parts)
  d[, tib_q := factor(bedrest_min, labels = c("q1", "q2", "q3", "q4", "q5"))]
  d[, isi_g := ifelse(isig == 1, "isi1", "isi0")]

  d[, tst_min := light_min + stage3_min + rem_min]
  # cal perc
  # d[, sol_perc := sol_min / (sol_min + waso_min + light_min + stage3_min + rem_min)]
  # d[, waso_perc := (waso_min / (waso_min + light_min + stage3_min + rem_min)) * 100]
  d[, light_perc := (light_min / (light_min + stage3_min + rem_min)) * 100]
  d[, stage3_perc := (stage3_min / (light_min + stage3_min + rem_min)) * 100]
  d[, rem_perc := (rem_min / (light_min + stage3_min + rem_min)) * 100]

  # make wide
  d <- dcast(d, . ~ isi_g + tib_q, fun.aggregate = mean, 
  value.var = c("wake_min", "sol_min", "waso_min", "light_min", "stage3_min", "rem_min", "tst_min",
                "light_perc", "stage3_perc", "rem_perc"))
  d <- d[, -1, with = FALSE] # remove the dot column

  # grand mean min and percentages
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
    d[, paste0(part, "_min_mean") := rowMeans(.SD, na.rm = TRUE),
    .SDcols = paste0(part, "_min_isi", rep(0:1, 5), "_q", rep(1:5, 2))]
  }
  for (part in c("light", "stage3", "rem")) {
    d[, paste0(part, "_perc_mean") := rowMeans(.SD, na.rm = TRUE),
    .SDcols = paste0(part, "_perc_isi", rep(0:1, 5), "_q", rep(1:5, 2))]
  }

  # mean of quantiles across isi
  for (q in 1:5) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
      d[, paste0(part, "_min_tib_q", q) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = c(paste0(part, "_min_isi0_q", q), paste0(part, "_min_isi1_q", q))
      ]
    }
  }
  # percentage of quantiles across isi
  for (q in 1:5) {
    for (part in c("light", "stage3", "rem")) {
      d[, paste0(part, "_perc_tib_q", q) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = c(paste0(part, "_perc_isi0_q", q), paste0(part, "_perc_isi1_q", q))
      ]
    }
  }

  # means across quartiles and ISI group
  for (isi in c("isi1", "isi0")) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
      d[, paste0(part, "_min_", isi) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = paste0(part, "_min_", isi, "_q", 1:5)
      ]
    }
  }
  for (isi in c("isi1", "isi0")) {
    for (part in c("light", "stage3", "rem")) {
      d[, paste0(part, "_perc_", isi) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = paste0(part, "_perc_", isi, "_q", 1:5)
      ]
    }
  }

  # contrasts
  # avg isi1 compared to avg isi0
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
    d[, paste0(part, "_min_isi1_isi0") := get(paste0(part, "_min_isi1")) - get(paste0(part, "_min_isi0"))]
  }
  for (part in c("light", "stage3", "rem")) {
    d[, paste0(part, "_perc_isi1_isi0") := get(paste0(part, "_perc_isi1")) - get(paste0(part, "_perc_isi0"))]
  }

  # isi1 q compared to avg isi0
  for (q in 1:5) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
      d[, paste0(part, "_min_isi1_q", q, "_isi0") := get(paste0(part, "_min_isi1_q", q)) - get(paste0(part, "_min_isi0"))]
    }
  }
  for (q in 1:5) {
    for (part in c("light", "stage3", "rem")) {
      d[, paste0(part, "_perc_isi1_q", q, "_isi0") := get(paste0(part, "_perc_isi1_q", q)) - get(paste0(part, "_perc_isi0"))]
    }
  }

  # contrasts
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
    d[, paste0(part, "_min_isi1_q1_q234") :=
      get(paste0(part, "_min_isi1_q1"))
      - ((get(paste0(part, "_min_isi1_q2")) + get(paste0(part, "_min_isi1_q3")) + get(paste0(part, "_min_isi1_q4"))) / 3)]

    d[, paste0(part, "_min_isi1_q5_q234") :=
      get(paste0(part, "_min_isi1_q5"))
      - ((get(paste0(part, "_min_isi1_q2")) + get(paste0(part, "_min_isi1_q3")) + get(paste0(part, "_min_isi1_q4"))) / 3)]

    d[, paste0(part, "_min_isi1_q1_q5") :=
      get(paste0(part, "_min_isi1_q1"))
      - get(paste0(part, "_min_isi1_q5"))]
  }
  for (part in c("light", "stage3", "rem")) {
    d[, paste0(part, "_perc_isi1_q1_q234") :=
      get(paste0(part, "_perc_isi1_q1"))
      - ((get(paste0(part, "_perc_isi1_q2")) + get(paste0(part, "_perc_isi1_q3")) + get(paste0(part, "_perc_isi1_q4"))) / 3)]

    d[, paste0(part, "_perc_isi1_q5_q234") :=
      get(paste0(part, "_perc_isi1_q5"))
      - ((get(paste0(part, "_perc_isi1_q2")) + get(paste0(part, "_perc_isi1_q3")) + get(paste0(part, "_perc_isi1_q4"))) / 3)]

    d[, paste0(part, "_perc_isi1_q1_q5") :=
      get(paste0(part, "_perc_isi1_q1"))
      - get(paste0(part, "_perc_isi1_q5"))]
  }
  d
})
saveRDS(pred_tibq5_isig_draws, file.path(out, "pred_tibq5_isig_draws.rds"))

# prep final data
pred_tibq5_isig_draws <- readRDS(file.path(out, "pred_tibq5_isig_draws.rds"))

pred_tibq5_isig_draws <- as.data.table(abind(pred_tibq5_isig_draws, along = 1))
# pred_tibq5_isig_draws <- split(pred_tibq5_isig_draws, pred_tibq5_isig_draws$bedrest_min)

pred_tibq5_isig_sum <- apply(pred_tibq5_isig_draws, 2, function(x) {
  describe_posterior(as.numeric(x), centrality = "mean", ci = 0.95)
})
pred_tibq5_isig_sum <- rbindlist(pred_tibq5_isig_sum)

pred_tibq5_isig_sum[, par := colnames(pred_tibq5_isig_draws)]

pred_tibq5_isig_sum[, part_label := NA]
pred_tibq5_isig_sum[, part_label := ifelse(grepl("sol", par) & grepl("min|perc", par), "SOL", part_label)]
pred_tibq5_isig_sum[, part_label := ifelse(grepl("waso", par) & grepl("min|perc", par), "WASO", part_label)]
pred_tibq5_isig_sum[, part_label := ifelse(grepl("light", par) & grepl("min|perc", par), "N1+2", part_label)]
pred_tibq5_isig_sum[, part_label := ifelse(grepl("stage3", par) & grepl("min|perc", par), "N3", part_label)]
pred_tibq5_isig_sum[, part_label := ifelse(grepl("rem", par) & grepl("min|perc", par), "REM", part_label)]
pred_tibq5_isig_sum[, part_label := ifelse(grepl("wake", par) & grepl("min|perc", par), "WD", part_label)]

pred_tibq5_isig_sum[, part_label := factor(part_label, ordered = TRUE, levels = c(
  "WD",
  "SOL", 
  "N1+2",
  "N3", 
  "REM",
  "WASO"
))]
table(pred_tibq5_isig_sum$part_label, useNA = "always")

pred_tibq5_isig_sum[, tib_group := NA]
pred_tibq5_isig_sum[, tib_group := ifelse(grepl("isi[0:1]_q1", par), "P10 TIB", tib_group)]
pred_tibq5_isig_sum[, tib_group := ifelse(grepl("isi[0:1]_q2", par), "P25 TIB", tib_group)]
pred_tibq5_isig_sum[, tib_group := ifelse(grepl("isi[0:1]_q3", par), "P50 TIB", tib_group)]
pred_tibq5_isig_sum[, tib_group := ifelse(grepl("isi[0:1]_q4", par), "P75 TIB", tib_group)]
pred_tibq5_isig_sum[, tib_group := ifelse(grepl("isi[0:1]_q5", par), "P90 TIB", tib_group)]
table(pred_tibq5_isig_sum$tib_group, useNA = "always")

pred_tibq5_isig_sum[, isi_group := NA]
pred_tibq5_isig_sum[, isi_group := ifelse(grepl("isi1", par) & !grepl("isi0", par), "Insomnia+", isi_group)]
pred_tibq5_isig_sum[, isi_group := ifelse(grepl("isi1", par) & grepl("isi0", par), "Insomnia+", isi_group)]
pred_tibq5_isig_sum[, isi_group := ifelse(grepl("isi0", par) & !grepl("isi1", par), "Insomnia-", isi_group)]
table(pred_tibq5_isig_sum$isi_group, useNA = "always")

pred_tibq5_isig_sum[, yintercept_healthy := NA]
pred_tibq5_isig_sum[, yintercept_healthy := ifelse(part_label == "SOL", pred_tibq5_isig_sum[par == "sol_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sum[, yintercept_healthy := ifelse(part_label == "WASO", pred_tibq5_isig_sum[par == "waso_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sum[, yintercept_healthy := ifelse(part_label == "N1+2", pred_tibq5_isig_sum[par == "light_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sum[, yintercept_healthy := ifelse(part_label == "N3", pred_tibq5_isig_sum[par == "stage3_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sum[, yintercept_healthy := ifelse(part_label == "REM", pred_tibq5_isig_sum[par == "rem_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sum[, yintercept_healthy := ifelse(part_label == "WD", pred_tibq5_isig_sum[par == "wake_isi0"]$Mean, yintercept_healthy)]

pred_tibq5_isig_sum[, yintercept_insom := NA]
pred_tibq5_isig_sum[, yintercept_insom := ifelse(part_label == "SOL", pred_tibq5_isig_sum[par == "sol_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sum[, yintercept_insom := ifelse(part_label == "WASO", pred_tibq5_isig_sum[par == "waso_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sum[, yintercept_insom := ifelse(part_label == "N1+2", pred_tibq5_isig_sum[par == "light_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sum[, yintercept_insom := ifelse(part_label == "N3", pred_tibq5_isig_sum[par == "stage3_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sum[, yintercept_insom := ifelse(part_label == "REM", pred_tibq5_isig_sum[par == "rem_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sum[, yintercept_insom := ifelse(part_label == "WD", pred_tibq5_isig_sum[par == "wake_isi1"]$Mean, yintercept_insom)]

# types of estimates
pred_tibq5_isig_sum[, contrast := NA]
pred_tibq5_isig_sum[, contrast := ifelse(grepl("isi1", par) & grepl("isi0", par), 1, contrast)]
pred_tibq5_isig_sum[, contrast := ifelse(grepl("q1_q5", par), 1, contrast)]
pred_tibq5_isig_sum[, contrast := ifelse(grepl("q1_q234", par), 1, contrast)]
pred_tibq5_isig_sum[, contrast := ifelse(grepl("q5_q234", par), 1, contrast)]
pred_tibq5_isig_sum[, contrast := ifelse(is.na(contrast), 0, contrast)]
pred_tibq5_isig_sum[, mean := ifelse(contrast == 0, 1, 0)]

pred_tibq5_isig_sum[, sig := ifelse(!between(0, CI_low, CI_high), "$^a$", "$\\phantom{^a}$")]
table(pred_tibq5_isig_sum$sig, useNA = "always")

part_labels <- c(
  "SOL" = "sol",
  "WASO" = "waso",
  "N1=2" = "light",
  "N3" = "stage3",
  "REM" = "rem",
  "WD" = "wake"
)
quartiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
pred_tibq5_isig_sum[, Contrast := NA]

for (q in quartiles) {
  for (pl in names(part_labels)) {
    var_prefix <- part_labels[[pl]]
    min_par_name <- paste0(var_prefix, "_min_isi1_", tolower(q), "_isi0")
    pred_tibq5_isig_sum[
      isi_group == "Insomnia+" & tib_group == q & part_label == pl,
      `:=`(
        Contrast = pred_tibq5_isig_sum[par == min_par_name]$Mean,
        Contrast_CI_low = pred_tibq5_isig_sum[par == min_par_name]$CI_low,
        Contrast_CI_high = pred_tibq5_isig_sum[par == min_par_name]$CI_high
      )
    ]

    if (var_prefix %in% c("wake", "sol", "waso")) next

    perc_par_name <- paste0(var_prefix, "_perc_isi1_", tolower(q), "_isi0")
    pred_tibq5_isig_sum[
      isi_group == "Insomnia+" & tib_group == q & part_label == pl,
      `:=`(
        Contrast_perc = pred_tibq5_isig_sum[par == perc_par_name]$Mean,
        Contrast_perc_CI_low = pred_tibq5_isig_sum[par == perc_par_name]$CI_low,
        Contrast_perc_CI_high = pred_tibq5_isig_sum[par == perc_par_name]$CI_high
      )
    ]
  }
}

pred_tibq5_isig_sum[, contrast_min_sig := ifelse(!between(0, Contrast_CI_low, Contrast_CI_high), "$*$", "$\\phantom{*}$")]
pred_tibq5_isig_sum[, contrast_perc_sig := ifelse(!between(0, Contrast_perc_CI_low, Contrast_perc_CI_high), "$*$", "$\\phantom{*}$")]

saveRDS(pred_tibq5_isig_sum, file = file.path(out, "pred_tibq5_isig_sum.rds"))
pred_tibq5_isig_sum <- readRDS(file.path(out, "pred_tibq5_isig_sum.rds"))

pred_tibq5_isig_all_sum <- pred_tibq5_isig_sum
pred_tibq5_isig_sum <- pred_tibq5_isig_sum[!is.na(tib_group) & isi_group == "Insomnia+" & contrast == 0]

pred_tibq5_isig_sum[grepl("min", par) & !grepl("perc", par), est_min := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
pred_tibq5_isig_sum[grepl("perc", par), est_perc := paste0(format(round(Mean, 2), nsmall = 2), " [", format(round(CI_low, 2), nsmall = 2), ", ", format(round(CI_high, 2), nsmall = 2), "]")]
pred_tibq5_isig_sum[, est_perc := zoo::na.locf(est_perc), by = .(part_label, tib_group, isi_group)]
pred_tibq5_isig_sum[, est := paste0(est_min, ", ", est_perc)]

### plot -----------------------
## plot_min -----------------------
# make individual then patch
plot_min_params <- list(
    "WD" = list(
    limits = c(700, 1450), breaks = c(900, 1000, 1100), breaks2 = c(.60, .70, .80), name = "WD", y_offset = 1.5
  ),
  "SOL" = list(
    limits = c(-30, 170), breaks = c(30, 60, 90), breaks2 = c(.01, .04, .07), name = "SOL", y_offset = 1.5
  ),
  "WASO" = list(
    limits = c(-10, 190), breaks = c(50, 80, 110), breaks2 = c(.03, .05, .07, .09), name = "WASO", y_offset = 1.5
  ),
  "N1+2" = list(
    limits = c(130, 390), breaks = c(200, 240, 280), breaks2 = c(.15, .20, .25), name = "N1+2", y_offset = 1.5
  ),
  "N3" = list(
    limits = c(-10, 120), breaks = c(30, 50, 70), breaks2 = c(.03, .05, .07, .09), name = "N3", y_offset = 1.5
  ),
  "REM" = list(
    limits = c(-10, 120), breaks = c(30, 50, 70), breaks2 = c(.03, .05, .07), name = "REM", y_offset = 1.5
  )
)

make_min_plot <- function(part_label) {
  params <- plot_min_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isig_sum[!is.na(tib_group) & part_label == part & grepl("min", par)],
    aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
    # geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      size = 0.75,
      linewidth = 1,
      position = position_dodge(width = 0.5)
    ) +
    geom_text(aes(y = Mean + params$y_offset, label = latex2exp::TeX(contrast_min_sig, output = "character")),
      parse = TRUE,
      hjust = 0.5, nudge_x = .2,
      family = "Arial Narrow",
      size = 7,
      show.legend = FALSE
    ) +
    geom_text(aes(y = max(params$limits), label = est_min),
      hjust = 1, nudge_x = 0,
      family = "Arial Narrow",
      fontface = "bold",
      size = 5.5,
      show.legend = FALSE
    ) +
    geom_text(aes(y = min(params$limits), label = tib_group),
      hjust = 0, nudge_x = 0,
      family = "Arial Narrow",
      fontface = "plain",
      size = 5.5,
      colour = "black",
      show.legend = FALSE
    ) +
    scale_y_continuous(
      limits = params$limits,
      breaks = params$breaks,
      labels = paste0(params$breaks),
    ) +
    scale_colour_manual(values = col) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.background    = element_rect(fill = "transparent", colour = NA, linewidth = 1),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.spacing       = unit(0.5, "lines"),
      axis.title.x        = element_blank(),
      axis.text.x         = element_text(size = 15, face = "plain", family = "Arial Narrow"),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 15, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(1, 1.5, 1.5, 0), "lines")
    )
}
plots_min <- lapply(names(plot_min_params), make_min_plot)
names(plots_min) <- names(plot_min_params)

grDevices::cairo_pdf(
  file = file.path(out, "plot_min.pdf"),
  width = 9,
  height = 11,
)
ggarrange(
  plots_min[["WD"]],
  plots_min[["N1+2"]],
  plots_min[["SOL"]],
  plots_min[["N3"]],
  plots_min[["WASO"]],
  plots_min[["REM"]],
  labels = c(
    "  A. Wake During the Day (min)",
    "  D. Light Sleep (min)",
    "  B. Sleep Onset Latency (min)",
    "  E. Slow Wave Sleep (min)",
    "  C. Wake After Sleep Onset (min)",
    "  F. REM Sleep (min)"
  ), 
  hjust = 0,
  ncol = 2, nrow = 3, common.legend = TRUE, legend = "none",
  font.label = list(size = 16, face = "bold", family = "Arial Narrow")
) + theme(plot.margin = unit(c(1, 0, 0, 1), "lines"))
dev.off()

## plot_perc -----------------------
plot_perc_params <- list(
  # "WASO" = list(
  #   limits = c(5, 25), breaks = c(5, 10, 15, 20), name = "WASO", y_offset = 0.005
  # ),
  "N1+2" = list(
    limits = c(40, 120), breaks = c(60, 70, 80), name = "N1+2", y_offset = 0.01
  ),
  "N3" = list(
    limits = c(-5, 50), breaks = c(10, 15, 20), name = "N3", y_offset = 0.005
  ),
  "REM" = list(
    limits = c(-5, 50), breaks = c(10, 15, 20), name = "REM", y_offset = 0.005
  )
)

make_perc_plot <- function(part_label) {
  params <- plot_perc_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isig_sum[!is.na(tib_group) & part_label == part & grepl("perc", par)],
    aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
    # geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
    geom_pointrange(
      aes(
      ymin = CI_low,
      ymax = CI_high
      ),
      size = 0.75,
      linewidth = 1,
      position = position_dodge(width = 0.5),
      shape = shape_perc # square
    ) +
    geom_text(aes(y = Mean + params$y_offset, label = latex2exp::TeX(contrast_perc_sig, output = "character")),
      parse = TRUE,
      hjust = 0.5, nudge_x = .2,
      family = "Arial Narrow",
      size = 7,
      show.legend = FALSE
    ) +
    geom_text(aes(y = max(params$limits), label = est_perc),
      hjust = 1, nudge_x = 0,
      family = "Arial Narrow",
      fontface = "bold",
      size = 5.5,
      show.legend = FALSE
    ) +
    geom_text(aes(y = min(params$limits), label = tib_group),
      hjust = 0, nudge_x = 0,
      family = "Arial Narrow",
      fontface = "plain",
      size = 5.5,
      colour = "black",
      show.legend = FALSE
    ) +
    scale_y_continuous(
      limits = params$limits,
      breaks = params$breaks,
      labels = paste0(params$breaks),
    ) +
    scale_colour_manual(values = col) +
    # scale_shape_manual(values = shape_perc) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.background    = element_rect(fill = "transparent", colour = NA, linewidth = 1),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.spacing       = unit(0.5, "lines"),
      axis.title.x        = element_blank(),
      axis.text.x         = element_text(size = 15, face = "plain", family = "Arial Narrow"),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 15, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(1, 1.5, 1.5, 0), "lines")
    )
}
plots_perc <- lapply(names(plot_perc_params), make_perc_plot)
names(plots_perc) <- names(plot_perc_params)

grDevices::cairo_pdf(
  file = file.path(out, "plot_perc.pdf"),
  width = 4.5,
  height = 11,
)
ggarrange(
  # plots_perc[["WASO"]],
  plots_perc[["N1+2"]],
  plots_perc[["N3"]],
  plots_perc[["REM"]],
  labels = c(
    # "  A. Wake After Sleep Onset (%)",
    "  A. Light Sleep (%)",
    "  B. Slow Wave Sleep (%)",
    "  C. REM Sleep (%)"
  ),
  hjust = 0,
  ncol = 1, nrow = 3, common.legend = TRUE, legend = "none",
  font.label = list(size = 16, face = "bold", family = "Arial Narrow")
) + theme(plot.margin = unit(c(1, 0, 0, 1), "lines"))
dev.off()

# combined plot ------------------------
grDevices::cairo_pdf(
  file = file.path(out, "plot_combined.pdf"),
  width = 13,
  height = 13,
)
ggarrange(
  plots_min[["WD"]],
  plots_min[["N1+2"]],
  plots_perc[["N1+2"]],

  plots_min[["SOL"]],
  plots_min[["N3"]],
  plots_perc[["N3"]],

  plots_min[["WASO"]],
  plots_min[["REM"]],
  plots_perc[["REM"]],

  labels = c(
    "  A. WD (min)",
    "  D. N1+2 (min)",
    "  G. N1+2 (%)",

    "  B. SOL (min)",
    "  E. N3 (min)",
    "  H. N3 (%)",

    "  C. WASO (min)",
    "  F. REM (min)",
    "  I. REM (%)"
  ),
  hjust = 0,
  ncol = 3, nrow = 3, 
  common.legend = TRUE, legend = "none",
  font.label = list(size = 16, face = "bold", family = "Arial Narrow")
) + theme(plot.margin = unit(c(1, 0, 0, 1), "lines"))
dev.off()

# table -----------------------------------------------
pred_tibq5_isig_all_sum[grepl("min", par) & !grepl("perc", par), est_min := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
pred_tibq5_isig_all_sum[grepl("perc", par), est_perc := paste0(format(round(Mean, 2), nsmall = 2), " [", format(round(CI_low, 2), nsmall = 2), ", ", format(round(CI_high, 2), nsmall = 2), "]")]

pred_tibq5_isig_sum_all[grepl("isi1", par) & grepl("isi0", par)][, est := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")][, .(par, est)]

pred_tibq5_isig_sum_all[grepl("q1_q234", par)][, est := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")][, .(par, est)]
pred_tibq5_isig_sum_all[grepl("q5_q234", par)][, est := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")][, .(par, est)]
pred_tibq5_isig_sum_all[grepl("q1_q5", par)][, est := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")][, .(par, est)]

pred_tibq5_isig_sum[grepl("q1_q234", par)][, .(par, est_min)]
pred_tibq5_isig_sum[grepl("q1_q234", par)][, .(par, est_perc)]

pred_tibq5_isig_sum[grepl("q5_q234", par)][, .(par, est_min)]
pred_tibq5_isig_sum[grepl("q5_q234", par)][, .(par, est_perc)]

pred_tibq5_isig_sum[grepl("q1_q5", par)][, .(par, est_min)]
pred_tibq5_isig_sum[grepl("q1_q5", par)][, .(par, est_perc)]
