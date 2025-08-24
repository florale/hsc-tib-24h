source("data.r")
parts <- c("wake_min", "sol_min", "waso_min", "light_min", "stage3_min", "rem_min")

clr <- complr(
  data = dpsg,
  parts = parts,
  # idvar = "record_id",
  total = 1440
)

m_tib_isig_sleep <- brmcoda(clr,
  mvbind(ilr1, ilr2, ilr3, ilr4, ilr5) ~
    s(bedrest_min, by = isig) +
    s(age) + female + white + s(perHrAHSleep),
  iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
  backend = "cmdstanr"
)
summary(m_tib_isig_sleep)
saveRDS(m_tib_isig_sleep, file.path(out, "m_tib_isig_sleep.rds"))

m_tib_isig_sleep <- readRDS(file.path(out, "m_tib_isig_sleep.rds"))

## 5 quartiles -----------------------
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

  d[, tib_q := factor(bedrest_min, labels = c("q1", "q2", "q3", "q4", "q5"))]
  d[, isi_g := ifelse(isig == 1, "isi1", "isi0")]

  # make wide
  d <- dcast(d, . ~ isi_g + tib_q, fun.aggregate = mean, value.var = c("wake_min", "sol_min", "waso_min", "light_min", "stage3_min", "rem_min"))
  d <- d[, -1, with = FALSE] # remove the dot column

  # mean
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
    d[, paste0(part, "_mean") := rowMeans(.SD, na.rm = TRUE), .SDcols = paste0("sol_min_isi", rep(0:1, 5), "_q", rep(1:5, 2))]
  }

  # mean of quantiles across isi
  for (q in 1:5) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
      d[, paste0(part, "_tib_q", q) := rowMeans(.SD, na.rm = TRUE), .SDcols = c(paste0(part, "_min_isi0_q", q), paste0(part, "_min_isi1_q", q))]
    }
  }

  # means across quartiles and ISI group
  for (isi in c("isi1", "isi0")) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
      d[, paste0(part, "_", isi) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = paste0(part, "_min_", isi, "_q", 1:5)
      ]
    }
  }

  # contrasts
  # avg isi1 compared to avg isi0
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
    d[, paste0(part, "_isi1_isi0") := get(paste0(part, "_isi1")) - get(paste0(part, "_isi0"))]
  }

  # isi1 q compared to avg isi0
  for (q in 1:5) {
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
      d[, paste0(part, "_isi1_q", q, "_isi0") := get(paste0(part, "_min_isi1_q", q)) - get(paste0(part, "_isi0"))]
    }
  }

  # difference between isi1, q1 vs q234
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
    d[, paste0(part, "_isi1_q234_q1") :=
      get(paste0(part, "_min_isi1_q1")) -
      ((get(paste0(part, "_min_isi1_q2")) + get(paste0(part, "_min_isi1_q3")) + get(paste0(part, "_min_isi1_q4"))) / 3)]
  }

  # difference between isi1, q5 vs q234
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
    d[, paste0(part, "_isi1_q234_q5") := get(paste0(part, "_min_isi1_q5")) -
      ((get(paste0(part, "_min_isi1_q2")) + get(paste0(part, "_min_isi1_q3")) + get(paste0(part, "_min_isi1_q4"))) / 3)]
  }

  # difference between isi1, q1 vs q5
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake")) {
    d[, paste0(part, "_isi1_q1_q5") := get(paste0(part, "_min_isi1_q1")) - get(paste0(part, "_min_isi1_q5"))]
  }
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
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("sol", par) & grepl("min", par), "Sleep onset latency", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("waso", par) & grepl("min", par), "Wake after sleep onset", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("light", par) & grepl("min", par), "Light sleep", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("stage3", par) & grepl("min", par), "Slow wave sleep", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("rem", par) & grepl("min", par), "REM sleep", part_label)]
pred_tibq5_isig_sleep_sum[, part_label := ifelse(grepl("wake", par) & grepl("min", par), "Daytime wake", part_label)]

pred_tibq5_isig_sleep_sum[, part_label := factor(part_label, ordered = TRUE, levels = c(
  "Daytime wake",
  "Sleep onset latency", "Wake after sleep onset", "Light sleep",
  "Slow wave sleep", "REM sleep"
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
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "Sleep onset latency", pred_tibq5_isig_sleep_sum[par == "sol_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "Wake after sleep onset", pred_tibq5_isig_sleep_sum[par == "waso_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "Light sleep", pred_tibq5_isig_sleep_sum[par == "light_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "Slow wave sleep", pred_tibq5_isig_sleep_sum[par == "stage3_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "REM sleep", pred_tibq5_isig_sleep_sum[par == "rem_isi0"]$Mean, yintercept_healthy)]
pred_tibq5_isig_sleep_sum[, yintercept_healthy := ifelse(part_label == "Daytime wake", pred_tibq5_isig_sleep_sum[par == "wake_isi0"]$Mean, yintercept_healthy)]

pred_tibq5_isig_sleep_sum[, yintercept_insom := NA]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "Sleep onset latency", pred_tibq5_isig_sleep_sum[par == "sol_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "Wake after sleep onset", pred_tibq5_isig_sleep_sum[par == "waso_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "Light sleep", pred_tibq5_isig_sleep_sum[par == "light_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "Slow wave sleep", pred_tibq5_isig_sleep_sum[par == "stage3_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "REM sleep", pred_tibq5_isig_sleep_sum[par == "rem_isi1"]$Mean, yintercept_insom)]
pred_tibq5_isig_sleep_sum[, yintercept_insom := ifelse(part_label == "Daytime wake", pred_tibq5_isig_sleep_sum[par == "wake_isi1"]$Mean, yintercept_insom)]

# types of estimates
pred_tibq5_isig_sleep_sum[, mean := ifelse(grepl("min", par) | grepl("mean", par), 1, 0)]
pred_tibq5_isig_sleep_sum[, contrast := ifelse(grepl("isi1", par) & grepl("isi0", par), 1, 0)]

table(pred_tibq5_isig_sleep_sum$mean, useNA = "always")
table(pred_tibq5_isig_sleep_sum$contrast, useNA = "always")

pred_tibq5_isig_sleep_sum[grepl("_isi1_q[1-5]_isi0", par), nonsig_isi10 := between(0, CI_low, CI_high)]
pred_tibq5_isig_sleep_sum[grepl("_q1_q5", par), nonsig_q1q5 := between(0, CI_low, CI_high)]
pred_tibq5_isig_sleep_sum[grepl("_q234_q1", par), nonsig_q234_q1 := between(0, CI_low, CI_high)]
pred_tibq5_isig_sleep_sum[grepl("_q234_q5", par), nonsig_q234_q5 := between(0, CI_low, CI_high)]

pred_tibq5_isig_sleep_sum[, sig := ifelse(nonsig_isi10 == FALSE, "$^a$", "$\\phantom{^a}$")]

part_labels <- c(
  "Sleep onset latency" = "sol",
  "Wake after sleep onset" = "waso",
  "Light sleep" = "light",
  "Slow wave sleep" = "stage3",
  "REM sleep" = "rem",
  "Daytime wake" = "wake"
)
quartiles <- c("Q1", "Q2", "Q3", "Q4", "Q5")
pred_tibq5_isig_sleep_sum[, Contrast := NA]

for (q in quartiles) {
  for (pl in names(part_labels)) {
    var_prefix <- part_labels[[pl]]
    par_name <- paste0(var_prefix, "_isi1_", tolower(q), "_isi0")
    pred_tibq5_isig_sleep_sum[
      isi_group == "Insomnia+" & tib_group == q & part_label == pl,
      `:=`(
        Contrast = pred_tibq5_isig_sleep_sum[par == par_name]$Mean,
        Contrast_CI_low = pred_tibq5_isig_sleep_sum[par == par_name]$CI_low,
        Contrast_CI_high = pred_tibq5_isig_sleep_sum[par == par_name]$CI_high
      )
    ]
  }
}
pred_tibq5_isig_sleep_sum[, nonsig_contrast := between(0, Contrast_CI_low, Contrast_CI_high)]
pred_tibq5_isig_sleep_sum[, sig_contrast := ifelse(nonsig_contrast == FALSE, "$*$", "$\\phantom{*}$")]

pred_tibq5_isig_sleep_sum <- pred_tibq5_isig_sleep_sum[!is.na(tib_group) & isi_group == "Insomnia+" & mean == 1]

### plot -----------------------
(plot_tibq5_isig_sleep <-
  ggplot(
    data = pred_tibq5_isig_sleep_sum,
    aes(x = tib_group, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.5, linetype = "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      y = Mean,
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  #  geom_text(aes(y = 500, label = cancer_time_since_diag_other),
  #            hjust = 0, nudge_x = 0,
  #            family = "Arial Narrow", size = 3,
  #            show.legend = FALSE) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  #  geom_text(aes(y = 650, label = Cases),
  #            hjust = 1, nudge_x = 0,
  #            family = "Arial Narrow", size = 3,
  #            show.legend = FALSE) +
  #  geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
  #  geom_segment(aes(x = 0, yend = 650), col = "black", linewidth = 0.5) +
  # scale_y_continuous(
  #   # limits = c(500, 650),
  #   # breaks = c(500, 650),
  #   name = "Sleep onset latency (mins/day)"
  # ) +
  facet_wrap(~part_label, ncol = 2, scales = "free") +
  scale_colour_manual(values = pal5) +
  labs(x = "", y = "", colour = "") +
  coord_flip() +
  theme_ipsum() +
  theme(
    axis.ticks          = element_blank(),
    # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 1),
    plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
    panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 1),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    # axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
    axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
    axis.text.x         = element_text(size = 10, face = "bold"),
    axis.text.y         = element_text(size = 10, face = "bold"),
    strip.text          = element_text(size = 9, face = "bold", hjust = .5),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)

grDevices::cairo_pdf(
  file = file.path(out, "plot_tibq5_isig_sleep.pdf"),
  width = 7,
  height = 10,
)
plot_tibq5_isig_sleep
dev.off()

# make individual then patch
plot_params <- list(
  "Sleep onset latency" = list(
    limits = c(0, 80), breaks = c(0, 80), name = "Sleep onset latency (min)", y_offset = 1.5
  ),
  "Daytime wake" = list(
    limits = c(800, 1200), breaks = c(800, 1200), name = "Daytime wake (min)", y_offset = 1.5
  ),
  "Wake after sleep onset" = list(
    limits = c(40, 80), breaks = c(40, 80), name = "Wake after sleep onset (min)", y_offset = 1.5
  ),
  "Light sleep" = list(
    limits = c(150, 350), breaks = c(150, 350), name = "Light sleep (min)", y_offset = 1.5
  ),
  "Slow wave sleep" = list(
    limits = c(40, 80), breaks = c(40, 80), name = "Slow wave sleep (min)", y_offset = 1.5
  ),
  "REM sleep" = list(
    limits = c(40, 100), breaks = c(40, 100), name = "REM sleep (min)", y_offset = 1.5
  )
)

make_tibq5_isig_plot <- function(part_label) {
  params <- plot_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isig_sleep_sum[!is.na(tib_group) & part_label == part],
    aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
    geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
    # geom_line(
    #   aes(y = Mean),
    #   size = 0.75,
    #   position = position_dodge(width = 0.5),
    #   colour = "#CBD5D0", linetype = "dashed"
    # ) +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      size = 0.75,
      linewidth = 0.75,
      position = position_dodge(width = 0.5)
    ) +
    geom_text(aes(y = Mean + params$y_offset, label = TeX(sig_contrast, output = "character")),
      parse = TRUE,
      hjust = 0.5, nudge_x = .2,
      family = "Arial Narrow", size = 6,
      show.legend = FALSE
    ) +
    scale_y_continuous(
      limits = params$limits,
      breaks = params$breaks,
      name = params$name,
      position = "right"
    ) +
    scale_colour_manual(values = pal5) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 1),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 13, face = "bold"),
      axis.text.y         = element_text(size = 13, face = "bold"),
      strip.text          = element_text(size = 13, face = "bold", hjust = .5),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5, 0.5, 1, 0.5), "lines")
    )
}
plots_tibq5_isig <- lapply(names(plot_params), make_tibq5_isig_plot)
names(plots_tibq5_isig) <- names(plot_params)

grDevices::cairo_pdf(
  file = file.path(out, "plot_tibq5_isig_all.pdf"),
  width = 12,
  height = 8,
)
ggarrange(
  plots_tibq5_isig[["Daytime wake"]],
  plots_tibq5_isig[["Sleep onset latency"]],
  plots_tibq5_isig[["Wake after sleep onset"]],
  plots_tibq5_isig[["Light sleep"]],
  plots_tibq5_isig[["Slow wave sleep"]],
  plots_tibq5_isig[["REM sleep"]],
  ncol = 3, nrow = 2, common.legend = TRUE, legend = "none"
)
dev.off()

