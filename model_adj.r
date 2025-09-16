source("setup.r")
parts <- c("wake_min", "sol_min", "light_min", "stage3_min", "rem_min", "waso_min")

# clr <- complr(
#   data = dpsg,
#   parts = parts,
#   # idvar = "record_id",
#   total = 1440
# )

# m_tib_isig_sleep <- brmcoda(clr,
#   mvbind(z1_1, z2_1, z3_1, z4_1, z5_1) ~
#     s(bedrest_min, by = isig) +
#     s(age) + female + white + s(perHrAHSleep),
#   iter = 4000, chains = 6, cores = 6, seed = 123, warmup = 1000,
#   backend = "cmdstanr"
# )
# summary(m_tib_isig_sleep)
# saveRDS(m_tib_isig_sleep, file.path(out, "m_tib_isig_sleep.rds"))

m_tib_isig_sleep <- readRDS(file.path(out, "m_tib_isig_sleep.rds"))

## 5 quartiles -----------------------
quantile(model.frame(m_tib_isig_sleep)$bedrest_min, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

d_tibq5_isig_sleep <- emmeans::ref_grid(m_tib_isig_sleep$model,
  at = list(
    bedrest_min = quantile(model.frame(m_tib_isig_sleep)$bedrest_min, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  )
)@grid
d_tibq5_isig_sleep <- as.data.table(d_tibq5_isig_sleep)
d_tibq5_isig_sleep <- d_tibq5_isig_sleep[!duplicated(d_tibq5_isig_sleep[, .(isig, bedrest_min, age, female, white, perHrAHSleep)]), ]

pred_tibq5_isig_sleep <- fitted(m_tib_isig_sleep, newdata = d_tibq5_isig_sleep, scale = "response", re_formula = NA, summary = FALSE)
pred_tibq5_isig_sleep <- apply(pred_tibq5_isig_sleep, c(1), function(x) cbind(d_tibq5_isig_sleep, x))

pred_tibq5_isig_sleep_draws <- lapply(pred_tibq5_isig_sleep, function(d) {
  d <- as.data.table(d)

  setnames(d, paste0("t", parts), parts)
  d[, tib_q := factor(bedrest_min, labels = c("q1", "q2", "q3", "q4", "q5"))]
  d[, isi_g := ifelse(isig == 1, "isi1", "isi0")]

  d[, tst_min := waso_min + light_min + stage3_min + rem_min]
  # cal perc
  # d[, sol_perc := sol_min / (sol_min + waso_min + light_min + stage3_min + rem_min)]
  d[, waso_perc := (waso_min / (waso_min + light_min + stage3_min + rem_min)) * 100]
  d[, light_perc := (light_min / (waso_min + light_min + stage3_min + rem_min)) * 100]
  d[, stage3_perc := (stage3_min / (waso_min + light_min + stage3_min + rem_min)) * 100]
  d[, rem_perc := (rem_min / (waso_min + light_min + stage3_min + rem_min)) * 100]

  # make wide
  d <- dcast(d, . ~ isi_g + tib_q, fun.aggregate = mean, 
  value.var = c("wake_min", "sol_min", "waso_min", "light_min", "stage3_min", "rem_min", "tst_min",
                "waso_perc", "light_perc", "stage3_perc", "rem_perc"))
  d <- d[, -1, with = FALSE] # remove the dot column

  # grand mean min and percentages
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
    d[, paste0(part, "_min_mean") := rowMeans(.SD, na.rm = TRUE),
    .SDcols = paste0(part, "_min_isi", rep(0:1, 5), "_q", rep(1:5, 2))]
  }
  for (part in c("waso", "light", "stage3", "rem")) {
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
    for (part in c("waso", "light", "stage3", "rem")) {
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
    for (part in c("waso", "light", "stage3", "rem")) {
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
  for (part in c("waso", "light", "stage3", "rem")) {
    d[, paste0(part, "_perc_isi1_isi0") := get(paste0(part, "_perc_isi1")) - get(paste0(part, "_perc_isi0"))]
  }

  # isi1 q compared to avg isi0
  for (q in 1:5) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
      d[, paste0(part, "_min_isi1_q", q, "_isi0") := get(paste0(part, "_min_isi1_q", q)) - get(paste0(part, "_min_isi0"))]
    }
  }
  for (q in 1:5) {
    for (part in c("waso", "light", "stage3", "rem")) {
      d[, paste0(part, "_perc_isi1_q", q, "_isi0") := get(paste0(part, "_perc_isi1_q", q)) - get(paste0(part, "_perc_isi0"))]
    }
  }

  # # difference between isi1, q1 vs q234
  # for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
  #   d[, paste0(part, "_isi1_q234_q1") :=
  #     get(paste0(part, "_min_isi1_q1")) -
  #     ((get(paste0(part, "_min_isi1_q2")) + get(paste0(part, "_min_isi1_q3")) + get(paste0(part, "_min_isi1_q4"))) / 3)]
  # }

  # # difference between isi1, q5 vs q234
  # for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
  #   d[, paste0(part, "_isi1_q234_q5") := get(paste0(part, "_min_isi1_q5")) -
  #     ((get(paste0(part, "_min_isi1_q2")) + get(paste0(part, "_min_isi1_q3")) + get(paste0(part, "_min_isi1_q4"))) / 3)]
  # }

  # # difference between isi1, q1 vs q5
  # for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
  #   d[, paste0(part, "_isi1_q1_q5") := get(paste0(part, "_min_isi1_q1")) - get(paste0(part, "_min_isi1_q5"))]
  # }
  d
})
saveRDS(pred_tibq5_isig_sleep_draws, file.path(out, "pred_tibq5_isig_sleep_draws.rds"))

pred_tibq5_isig_sleep_draws <- readRDS(file.path(out, "pred_tibq5_isig_sleep_draws.rds"))

pred_tibq5_isig_sleep_draws <- as.data.table(abind(pred_tibq5_isig_sleep_draws, along = 1))
# pred_tibq5_isig_sleep_draws <- split(pred_tibq5_isig_sleep_draws, pred_tibq5_isig_sleep_draws$bedrest_min)

pred_tibq5_isig_sleep_sum <- apply(pred_tibq5_isig_sleep_draws, 2, function(x) {
  describe_posterior(as.numeric(x), centrality = "mean", ci = 0.95)
})
pred_tibq5_isig_sleep_sum <- rbindlist(pred_tibq5_isig_sleep_sum)

pred_tibq5_isig_sleep_sum[, par := colnames(pred_tibq5_isig_sleep_draws)]

pred_tibq5_isig_sleep_sum[, part_label := NA]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("sol", par) & grepl("min|perc", par), "SOL", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("waso", par) & grepl("min|perc", par), "WASO", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("light", par) & grepl("min|perc", par), "N1+2", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("stage3", par) & grepl("min|perc", par), "N3", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("rem", par) & grepl("min|perc", par), "REM", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("wake", par) & grepl("min|perc", par), "WD", part_label)]

pred_tibq5_isig_sleep_sum[, part_label := factor(part_label, ordered = TRUE, levels = c(
  "WD",
  "SOL", 
  "N1+2",
  "N3", 
  "REM",
  "WASO"
))]
table(pred_tibq5_isig_sleep_sum$part_label, useNA = "always")

pred_tibq5_isig_sleep_sum[, tib_group := NA]
pred_tibq5_isig_sleep_sum[, tib_group := ifelse(grepl("isi[0:1]_q1", par), "Q1", tib_group)]
pred_tibq5_isig_sleep_sum[, tib_group := ifelse(grepl("isi[0:1]_q2", par), "Q2", tib_group)]
pred_tibq5_isig_sleep_sum[, tib_group := ifelse(grepl("isi[0:1]_q3", par), "Q3", tib_group)]
pred_tibq5_isig_sleep_sum[, tib_group := ifelse(grepl("isi[0:1]_q4", par), "Q4", tib_group)]
pred_tibq5_isig_sleep_sum[, tib_group := ifelse(grepl("isi[0:1]_q5", par), "Q5", tib_group)]
table(pred_tibq5_isig_sleep_sum$tib_group, useNA = "always")

pred_tibq5_isig_sleep_sum[, isi_group := NA]
pred_tibq5_isig_sleep_sum[, isi_group := ifelse(grepl("isi1", par) & !grepl("isi0", par), "Insomnia+", isi_group)]
pred_tibq5_isig_sleep_sum[, isi_group := ifelse(grepl("isi1", par) & grepl("isi0", par), "Insomnia+", isi_group)]
pred_tibq5_isig_sleep_sum[, isi_group := ifelse(grepl("isi0", par) & !grepl("isi1", par), "Insomnia-", isi_group)]
table(pred_tibq5_isig_sleep_sum$isi_group, useNA = "always")

pred_tibq5_isig_sleep_sum[, yintercept_healthy := NA]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "SOL", pred_tibq5_isig_sleep_sum[par == "sol_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "WASO", pred_tibq5_isig_sleep_sum[par == "waso_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "N1+2", pred_tibq5_isig_sleep_sum[par == "light_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "N3", pred_tibq5_isig_sleep_sum[par == "stage3_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "REM", pred_tibq5_isig_sleep_sum[par == "rem_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "WD", pred_tibq5_isig_sleep_sum[par == "wake_isi0"]$Mean, yintercept_healthy)]

pred_tibq5_isig_sleep_sum[, yintercept_insom := NA]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "SOL", pred_tibq5_isig_sleep_sum[par == "sol_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "WASO", pred_tibq5_isig_sleep_sum[par == "waso_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "N1+2", pred_tibq5_isig_sleep_sum[par == "light_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "N3", pred_tibq5_isig_sleep_sum[par == "stage3_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "REM", pred_tibq5_isig_sleep_sum[par == "rem_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "WD", pred_tibq5_isig_sleep_sum[par == "wake_isi1"]$Mean, yintercept_insom)]

# types of estimates
pred_tibq5_isig_sleep_sum[, contrast := ifelse(grepl("isi1", par) & grepl("isi0", par), 1, 0)]

# pred_tibq5_isig_sleep_sum[, perc := ifelse(grepl("perc", par), 1, 0)]

table(pred_tibq5_isig_sleep_sum$mean, useNA = "always")
table(pred_tibq5_isig_sleep_sum$contrast, useNA = "always")
# table(pred_tibq5_isig_sleep_sum$perc, useNA = "always")

pred_tibq5_isig_sleep_sum[grepl("_isi1_q[1-5]_isi0", par), nonsig_isi10 := between(0, CI_low, CI_high)]
# pred_tibq5_isig_sleep_sum[grepl("_q1_q5", par), nonsig_q1q5 := between(0, CI_low, CI_high)]
# pred_tibq5_isig_sleep_sum[grepl("_q234_q1", par), nonsig_q234_q1 := between(0, CI_low, CI_high)]
# pred_tibq5_isig_sleep_sum[grepl("_q234_q5", par), nonsig_q234_q5 := between(0, CI_low, CI_high)]

pred_tibq5_isig_sleep_sum[, sig := ifelse(nonsig_isi10 == FALSE, "$^a$", "$\\phantom{^a}$")]

part_labels <- c(
  "SOL" = "sol",
  "WASO" = "waso",
  "N1=2" = "light",
  "N3" = "stage3",
  "REM" = "rem",
  "WD" = "wake"
)
quartiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
pred_tibq5_isig_sleep_sum[, Contrast := NA]

for (q in quartiles) {
  for (pl in names(part_labels)) {
    var_prefix <- part_labels[[pl]]
    min_par_name <- paste0(var_prefix, "_min_isi1_", tolower(q), "_isi0")
    pred_tibq5_isig_sleep_sum[
      isi_group == "Insomnia+" & tib_group == q & part_label == pl,
      `:=`(
        Contrast = pred_tibq5_isig_sleep_sum[par == min_par_name]$Mean,
        Contrast_CI_low = pred_tibq5_isig_sleep_sum[par == min_par_name]$CI_low,
        Contrast_CI_high = pred_tibq5_isig_sleep_sum[par == min_par_name]$CI_high
      )
    ]
    if (var_prefix %in% c("wake", "sol")) next
    perc_par_name <- paste0(var_prefix, "_perc_isi1_", tolower(q), "_isi0")
    pred_tibq5_isig_sleep_sum[
      isi_group == "Insomnia+" & tib_group == q & part_label == pl,
      `:=`(
        Contrast_perc = pred_tibq5_isig_sleep_sum[par == perc_par_name]$Mean,
        Contrast_perc_CI_low = pred_tibq5_isig_sleep_sum[par == perc_par_name]$CI_low,
        Contrast_perc_CI_high = pred_tibq5_isig_sleep_sum[par == perc_par_name]$CI_high
      )
    ]
  }
}

pred_tibq5_isig_sleep_sum[, contrast_min_nonsig := between(0, Contrast_CI_low, Contrast_CI_high)]
pred_tibq5_isig_sleep_sum[, contrast_min_sig := ifelse(contrast_min_nonsig == FALSE, "$*$", "$\\phantom{*}$")]

pred_tibq5_isig_sleep_sum[, contrast_perc_nonsig := between(0, Contrast_perc_CI_low, Contrast_perc_CI_high)]
pred_tibq5_isig_sleep_sum[, contrast_perc_sig := ifelse(contrast_perc_nonsig == FALSE, "$*$", "$\\phantom{*}$")]

pred_tibq5_isig_sleep_sum_all <- pred_tibq5_isig_sleep_sum
pred_tibq5_isig_sleep_sum <- pred_tibq5_isig_sleep_sum[!is.na(tib_group) & isi_group == "Insomnia+" & contrast == 0]

pred_tibq5_isig_sleep_sum[grepl("min", par) & !grepl("perc", par), est_min := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
pred_tibq5_isig_sleep_sum[grepl("perc", par), est_perc := paste0(round(Mean, 2), " [", round(CI_low, 2), ", ", round(CI_high, 2), "]")]
pred_tibq5_isig_sleep_sum[, est_perc := zoo::na.locf(est_perc), by = .(part_label, tib_group, isi_group)]
pred_tibq5_isig_sleep_sum[, est := paste0(est_min, ", ", est_perc)]

### plot -----------------------
# make individual then patch
plot_min_params <- list(
    "WD" = list(
    limits = c(850, 1300), breaks = c(900, 1000, 1100, 1200), breaks2 = c(.60, .70, .80), name = "WD", y_offset = 1.5
  ),
  "SOL" = list(
    limits = c(10, 100), breaks = c(20, 40, 60, 80), breaks2 = c(.01, .04, .07), name = "SOL", y_offset = 1.5
  ),
  "WASO" = list(
    limits = c(30, 120), breaks = c(40, 60, 80, 100), breaks2 = c(.03, .05, .07, .09), name = "WASO", y_offset = 1.5
  ),
  "N1+2" = list(
    limits = c(180, 360), breaks = c(200, 240, 280, 320), breaks2 = c(.15, .20, .25), name = "N1+2", y_offset = 1.5
  ),
  "N3" = list(
    limits = c(30, 120), breaks = c(40, 60, 80, 100), breaks2 = c(.03, .05, .07, .09), name = "N3", y_offset = 1.5
  ),
  "REM" = list(
    limits = c(30, 120), breaks = c(40, 60, 80, 100), breaks2 = c(.03, .05, .07), name = "REM", y_offset = 1.5
  )
)

make_min_plot <- function(part_label) {
  params <- plot_min_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isig_sleep_sum[!is.na(tib_group) & part_label == part & grepl("min", par)],
    aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
    # geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
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
      size = 5,
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
      axis.text.x         = element_text(size = 14, face = "plain", family = "Arial Narrow"),
      axis.text.y         = element_text(size = 14, face = "plain", family = "Arial Narrow", margin = margin(l = 20)),
      strip.text          = element_text(size = 14, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(2, 1.5, 1.5, -2), "lines")
    )
}
plots_min <- lapply(names(plot_min_params), make_min_plot)
names(plots_min) <- names(plot_min_params)

grDevices::cairo_pdf(
  file = file.path(out, "plot_min.pdf"),
  width = 10,
  height = 12,
)
ggarrange(
  plots_min[["WD"]],
  plots_min[["SOL"]],
  plots_min[["WASO"]],
  plots_min[["N1+2"]],
  plots_min[["N3"]],
  plots_min[["REM"]],
  labels = c(
    "  A. Wake During the Day (min)",
    "  B. Sleep Onset Latency (min)",
    "  C. Wake After Sleep Onset (min)",
    "  D. Light Sleep (min)",
    "  E. Slow Wave Sleep (min)",
    "  F. REM Sleep (min)"
  ), 
  hjust = 0,
  ncol = 2, nrow = 3, common.legend = TRUE, legend = "none",
  font.label = list(size = 14, face = "bold", family = "Arial Narrow")
) + theme(plot.margin = unit(c(1, 0, 0, 1), "lines"))
dev.off()

# plot_perc
# make individual then patch
plot_perc_params <- list(
  "WASO" = list(
    limits = c(5, 25), breaks = c(5, 10, 15, 20), name = "WASO", y_offset = 0.005
  ),
  "N1+2" = list(
    limits = c(50, 70), breaks = c(50, 55, 60, 65), name = "N1+2", y_offset = 0.01
  ),
  "N3" = list(
    limits = c(5, 25), breaks = c(5, 10, 15, 20), name = "N3", y_offset = 0.005
  ),
  "REM" = list(
    limits = c(5, 25), breaks = c(5, 10, 15, 20), name = "REM", y_offset = 0.005
  )
)

make_perc_plot <- function(part_label) {
  params <- plot_perc_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isig_sleep_sum[!is.na(tib_group) & part_label == part & grepl("perc", par)],
    aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
    # geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      size = 0.75,
      linewidth = 1,
      position = position_dodge(width = 0.5)
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
      size = 5,
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
      axis.text.x         = element_text(size = 14, face = "plain", family = "Arial Narrow"),
      axis.text.y         = element_text(size = 14, face = "plain", family = "Arial Narrow", margin = margin(l = 20)),
      strip.text          = element_text(size = 14, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(2, 1, 1.5, -2), "lines")
    )
}
plots_perc <- lapply(names(plot_perc_params), make_perc_plot)
names(plots_perc) <- names(plot_perc_params)

grDevices::cairo_pdf(
  file = file.path(out, "plot_perc.pdf"),
  width = 10,
  height = 8,
)
ggarrange(
  plots_perc[["WASO"]],
  plots_perc[["N1+2"]],
  plots_perc[["N3"]],
  plots_perc[["REM"]],
  labels = c(
    "  A. Wake After Sleep Onset (%)",
    "  B. Light Sleep (%)",
    "  C. Slow Wave Sleep (%)",
    "  D. REM Sleep (%)"
  ),
  hjust = 0,
  ncol = 2, nrow = 2, common.legend = TRUE, legend = "none",
  font.label = list(size = 14, face = "bold", family = "Arial Narrow")
) + theme(plot.margin = unit(c(1, 0, 0, 1), "lines"))
dev.off()

## table -----------------------------------------------
# min diff
pred_tibq5_isig_sleep_sum_all[grepl("isi1", par) & grepl("isi0", par)  & grepl("min", par)][, est := paste0(round(Mean, 1), "[", round(CI_low, 1), ", ", round(CI_high, 1), "]")][, .(par, est)]

# perc diff
pred_tibq5_isig_sleep_sum_all[grepl("isi1", par) & grepl("isi0", par)  & grepl("perc", par)][, est := paste0(round(Mean, 2), "[", round(CI_low, 2), ", ", round(CI_high, 2), "]")][, .(par, est)]
