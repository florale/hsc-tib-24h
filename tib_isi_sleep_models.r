source("data.r")
parts <- c("wake_min", "sol_min", "waso_min", "light_min", "stage3_min", "rem_min")

clr <- complr(
  data = dpsg,
  parts = parts,
  # idvar = "record_id",
  total = 1440
)

m_tib_isig_sleep <- brmcoda(clr,
  mvbind(ilr1, ilr2, ilr3, ilr4, ilr5) ~ s(bedrest_min, by = isig) + age + female,
  iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
  backend = "cmdstanr"
)
summary(m_tib_isig_sleep)
saveRDS(m_tib_isig_sleep, paste0(out, "m_tib_isig_sleep", ".RDS"))

m_tib_isig_sleep <- readRDS(paste0(out, "m_tib_isig_sleep", ".RDS"))

## 5 quartiles -----------------------
d_tibq5_isig_sleep <- emmeans::ref_grid(m_tib_isig_sleep$model,
  at = list(
    bedrest_min = quantile(model.frame(m_tib_isig_sleep)$bedrest_min, probs = seq(0, 1, length.out = 5), na.rm = TRUE)
  )
)@grid
d_tibq5_isig_sleep <- as.data.table(d_tibq5_isig_sleep)
d_tibq5_isig_sleep <- d_tibq5_isig_sleep[!duplicated(d_tibq5_isig_sleep[, .(isig, bedrest_min, age, female)]), ]

pred_tibq5_isig_sleep <- fitted(m_tib_isig_sleep, newdata = d_tibq5_isig_sleep, scale = "response", re_formula = NA, summary = FALSE)
pred_tibq5_isig_sleep <- apply(pred_tibq5_isig_sleep, c(1), function(x) cbind(d_tibq5_isig_sleep, x))

pred_tibq5_isig_sleep_draws <- lapply(pred_tibq5_isig_sleep, function(d) {
  d <- as.data.table(d)

  d[, tib_q := factor(bedrest_min, labels = c("q1", "q2", "q3", "q4", "q5"))] # 5 quartiles
  d[, isi_g := ifelse(isig == 1, "isi1", "isi0")] # 5 quartiles

  # make wide
  d <- dcast(d, . ~ isi_g + tib_q, fun.aggregate = mean, value.var = parts)
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
        .SDcols = paste0(part, "_min_", isi, "_q", 1:5)]
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

  # difference between isi1,q5 vs q234
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
saveRDS(pred_tibq5_isig_sleep_draws, paste0(out, "pred_tibq5_isig_sleep_draws", ".RDS"))
pred_tibq5_isig_sleep_draws <- readRDS(paste0(out, "pred_tibq5_isig_sleep_draws", ".RDS"))

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

# assign tib_group based on both isi and q[1:5] in par
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

pred_tibq5_isig_sleep_sum[, mean := ifelse(grepl("min", par) | grepl("mean", par), 1, 0)]
pred_tibq5_isig_sleep_sum[, contrast := ifelse(grepl("isi1", par) & grepl("isi0", par), 1, 0)]

table(pred_tibq5_isig_sleep_sum$mean, useNA = "always")
table(pred_tibq5_isig_sleep_sum$contrast, useNA = "always")
pred_tibq5_isig_sleep_sum[!is.na(tib_group) & isi_group == "Insomnia+" & contrast == 1]

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
(plot_comp_tibq5_isig_sleep <-
  ggplot(
    data = pred_tibq5_isig_sleep_sum,
    aes(x = tib_group, group = isi_group, colour = interaction(tib_group, isi_group))
  ) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
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
  file = paste0(out, "plot_comp_tibq5_isig_sleep", ".pdf"),
  width = 7,
  height = 10,
)
plot_comp_tibq5_isig_sleep
dev.off()

# Make individual then patch
plot_params <- list(
  "Sleep onset latency" = list(
    limits = c(0, 800), breaks = c(0, 800), name = "Sleep onset latency (min)", y_offset = 1.5
  ),
  "Daytime wake" = list(
    limits = c(350, 1400), breaks = c(400, 1400), name = "Daytime wake (min)", y_offset = 1.5
  ),
  "Wake after sleep onset" = list(
    limits = c(0, 150), breaks = c(0, 150), name = "Wake after sleep onset (min)", y_offset = 1.5
  ),
  "Light sleep" = list(
    limits = c(0, 350), breaks = c(0, 350), name = "Light sleep (min)", y_offset = 1.5
  ),
  "Slow wave sleep" = list(
    limits = c(0, 200), breaks = c(0, 200), name = "Slow wave sleep (min)", y_offset = 1.5
  ),
  "REM sleep" = list(
    limits = c(0, 150), breaks = c(0, 150), name = "REM sleep (min)", y_offset = 1.5
  )
)

make_tibq5_isig_plot <- function(part_label) {
  params <- plot_params[[part_label]]
  part <- part_label
  ggplot(pred_tibq5_isig_sleep_sum[!is.na(tib_group) & part_label == part],
         aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      size = .25,
      linewidth = 0.5,
      position = position_dodge(width = 0.5)
    ) +
    geom_text(aes(y = Mean + params$y_offset, label = TeX(sig_contrast, output = "character")),
      parse = TRUE,
      hjust = 0.5, nudge_x = .15,
      family = "Arial Narrow", size = 4,
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
      axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 10, face = "bold"),
      axis.text.y         = element_text(size = 10, face = "bold"),
      strip.text          = element_text(size = 9, face = "bold", hjust = .5),
      legend.text         = element_text(size = 10, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5, 0.5, 1, 0.5), "lines")
    )
}
plots_tibq5_isig <- lapply(names(plot_params), make_tibq5_isig_plot)
names(plots_tibq5_isig) <- names(plot_params)

grDevices::cairo_pdf(
  file = paste0(out, "plot_comp_tibq5_isig_all", ".pdf"),
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

## 4 quartiles -----------------------
d_tibq4_isig_sleep <- emmeans::ref_grid(m_tib_isig_sleep$model,
  at = list(
    bedrest_min = quantile(model.frame(m_tib_isig_sleep)$bedrest_min, probs = seq(0, 1, length.out = 4), na.rm = TRUE)
  )
)@grid
d_tibq4_isig_sleep <- as.data.table(d_tibq4_isig_sleep)
d_tibq4_isig_sleep <- d_tibq4_isig_sleep[!duplicated(d_tibq4_isig_sleep[, .(isig, bedrest_min, age, female)]), ]

pred_tibq4_isig_sleep <- fitted(m_tib_isig_sleep, newdata = d_tibq4_isig_sleep, scale = "response", re_formula = NA, summary = FALSE)
pred_tibq4_isig_sleep <- apply(pred_tibq4_isig_sleep, c(1), function(x) cbind(d_tibq4_isig_sleep, x))

pred_tibq4_isig_sleep_draws <- lapply(pred_tibq4_isig_sleep, function(d) {
  d <- as.data.table(d)

  d[, tib_q := factor(bedrest_min, labels = c("q1", "q2", "q3", "q4"))] # 4 quartiles

  # mean
  d[, sol_tib := mean(sol_min, na.rm = TRUE), by = isig]
  d[, waso_tib := mean(waso_min, na.rm = TRUE), by = isig]
  d[, light_tib := mean(light_min, na.rm = TRUE), by = isig]
  d[, stage3_tib := mean(stage3_min, na.rm = TRUE), by = isig]
  d[, rem_tib := mean(rem_min, na.rm = TRUE), by = isig]
  d[, wake_tib := mean(wake_min, na.rm = TRUE), by = isig]

  d[, sol_tib1 := mean(sol_min[tib_q == "q1"], na.rm = TRUE), by = isig]
  d[, waso_tib1 := mean(waso_min[tib_q == "q1"], na.rm = TRUE), by = isig]
  d[, light_tib1 := mean(light_min[tib_q == "q1"], na.rm = TRUE), by = isig]
  d[, stage3_tib1 := mean(stage3_min[tib_q == "q1"], na.rm = TRUE), by = isig]
  d[, rem_tib1 := mean(rem_min[tib_q == "q1"], na.rm = TRUE), by = isig]
  d[, wake_tib1 := mean(wake_min[tib_q == "q1"], na.rm = TRUE), by = isig]

  d[, sol_tib2 := mean(sol_min[tib_q == "q2"], na.rm = TRUE), by = isig]
  d[, waso_tib2 := mean(waso_min[tib_q == "q2"], na.rm = TRUE), by = isig]
  d[, light_tib2 := mean(light_min[tib_q == "q2"], na.rm = TRUE), by = isig]
  d[, stage3_tib2 := mean(stage3_min[tib_q == "q2"], na.rm = TRUE), by = isig]
  d[, rem_tib2 := mean(rem_min[tib_q == "q2"], na.rm = TRUE), by = isig]
  d[, wake_tib2 := mean(wake_min[tib_q == "q2"], na.rm = TRUE), by = isig]

  d[, sol_tib3 := mean(sol_min[tib_q == "q3"], na.rm = TRUE), by = isig]
  d[, waso_tib3 := mean(waso_min[tib_q == "q3"], na.rm = TRUE), by = isig]
  d[, light_tib3 := mean(light_min[tib_q == "q3"], na.rm = TRUE), by = isig]
  d[, stage3_tib3 := mean(stage3_min[tib_q == "q3"], na.rm = TRUE), by = isig]
  d[, rem_tib3 := mean(rem_min[tib_q == "q3"], na.rm = TRUE), by = isig]
  d[, wake_tib3 := mean(wake_min[tib_q == "q3"], na.rm = TRUE), by = isig]

  d[, sol_tib4 := mean(sol_min[tib_q == "q4"], na.rm = TRUE), by = isig]
  d[, waso_tib4 := mean(waso_min[tib_q == "q4"], na.rm = TRUE), by = isig]
  d[, light_tib4 := mean(light_min[tib_q == "q4"], na.rm = TRUE), by = isig]
  d[, stage3_tib4 := mean(stage3_min[tib_q == "q4"], na.rm = TRUE), by = isig]
  d[, rem_tib4 := mean(rem_min[tib_q == "q4"], na.rm = TRUE), by = isig]
  d[, wake_tib4 := mean(wake_min[tib_q == "q4"], na.rm = TRUE), by = isig]

  d[, sol_isi := mean(sol_min, na.rm = TRUE), by = tib_q]
  d[, waso_isi := mean(waso_min, na.rm = TRUE), by = tib_q]
  d[, light_isi := mean(light_min, na.rm = TRUE), by = tib_q]
  d[, stage3_isi := mean(stage3_min, na.rm = TRUE), by = tib_q]
  d[, rem_isi := mean(rem_min, na.rm = TRUE), by = tib_q]
  d[, wake_isi := mean(wake_min, na.rm = TRUE), by = tib_q]

  d[, sol_isi1 := mean(sol_min[isig == 1], na.rm = TRUE), by = tib_q]
  d[, waso_isi1 := mean(waso_min[isig == 1], na.rm = TRUE), by = tib_q]
  d[, light_isi1 := mean(light_min[isig == 1], na.rm = TRUE), by = tib_q]
  d[, stage3_isi1 := mean(stage3_min[isig == 1], na.rm = TRUE), by = tib_q]
  d[, rem_isi1 := mean(rem_min[isig == 1], na.rm = TRUE), by = tib_q]
  d[, wake_isi1 := mean(wake_min[isig == 1], na.rm = TRUE), by = tib_q]

  d[, sol_isi0 := mean(sol_min[isig == 0], na.rm = TRUE), by = tib_q]
  d[, waso_isi0 := mean(waso_min[isig == 0], na.rm = TRUE), by = tib_q]
  d[, light_isi0 := mean(light_min[isig == 0], na.rm = TRUE), by = tib_q]
  d[, stage3_isi0 := mean(stage3_min[isig == 0], na.rm = TRUE), by = tib_q]
  d[, rem_isi0 := mean(rem_min[isig == 0], na.rm = TRUE), by = tib_q]
  d[, wake_isi0 := mean(wake_min[isig == 0], na.rm = TRUE), by = tib_q]

  # contrasts
  d[, sol_isi1_isi0 := sol_isi1 - sol_isi0, by = tib_q]
  d[, waso_isi1_isi0 := waso_isi1 - waso_isi0, by = tib_q]
  d[, light_isi1_isi0 := light_isi1 - light_isi0, by = tib_q]
  d[, stage3_isi1_isi0 := stage3_isi1 - stage3_isi0, by = tib_q]
  d[, rem_isi1_isi0 := rem_isi1 - rem_isi0, by = tib_q]
  d[, wake_isi1_isi0 := wake_isi1 - wake_isi0, by = tib_q]

  # difference between q1 vs q23
  d[, sol_q23_q1 := (sol_tib2 + sol_tib3) / 2 - sol_tib1, by = isig]
  d[, waso_q23_q1 := (waso_tib2 + waso_tib3) / 2 - waso_tib1, by = isig]
  d[, light_q23_q1 := (light_tib2 + light_tib3) / 2 - light_tib1, by = isig]
  d[, stage3_q23_q1 := (stage3_tib2 + stage3_tib3) / 2 - stage3_tib1, by = isig]
  d[, rem_q23_q1 := (rem_tib2 + rem_tib3) / 2 - rem_tib1, by = isig]
  d[, wake_q23_q1 := (wake_tib2 + wake_tib3) / 2 - wake_tib1, by = isig]

  # difference between q4 vs q23
  d[, sol_q23_q4 := (sol_tib2 + sol_tib3) / 2 - sol_tib4, by = isig]
  d[, waso_q23_q4 := (waso_tib2 + waso_tib3) / 2 - waso_tib4, by = isig]
  d[, light_q23_q4 := (light_tib2 + light_tib3) / 2 - light_tib4, by = isig]
  d[, stage3_q23_q4 := (stage3_tib2 + stage3_tib3) / 2 - stage3_tib4, by = isig]
  d[, rem_q23_q4 := (rem_tib2 + rem_tib3) / 2 - rem_tib4, by = isig]
  d[, wake_q23_q4 := (wake_tib2 + wake_tib3) / 2 - wake_tib4, by = isig]

  # difference between q1 and q4
  d[, sol_q1_q4 := sol_tib1 - sol_tib4, by = isig]
  d[, waso_q1_q4 := waso_tib1 - waso_tib4, by = isig]
  d[, light_q1_q4 := light_tib1 - light_tib4, by = isig]
  d[, stage3_q1_q4 := stage3_tib1 - stage3_tib4, by = isig]
  d[, rem_q1_q4 := rem_tib1 - rem_tib4, by = isig]
  d[, wake_q1_q4 := wake_tib1 - wake_tib4, by = isig]

  d <- d[female == 1 & tib_q == "q1"] # same across tib_q, so just take one
  d
})
saveRDS(pred_tibq4_isig_sleep_draws, paste0(out, "pred_tibq4_isig_sleep_draws", ".RDS"))

pred_tibq4_isig_sleep_draws <- as.data.table(abind(pred_tibq4_isig_sleep_draws, along = 1))
# pred_tibq4_isig_sleep_draws <- split(pred_tibq4_isig_sleep_draws, pred_tibq4_isig_sleep_draws$bedrest_min)

# mean
parts_mean <- c(
  "sol_tib", "waso_tib", "light_tib", "stage3_tib", "rem_tib", "wake_tib",
  "sol_tib1", "waso_tib1", "light_tib1", "stage3_tib1", "rem_tib1", "wake_tib1",
  "sol_tib2", "waso_tib2", "light_tib2", "stage3_tib2", "rem_tib2", "wake_tib2",
  "sol_tib3", "waso_tib3", "light_tib3", "stage3_tib3", "rem_tib3", "wake_tib3",
  "sol_tib4", "waso_tib4", "light_tib4", "stage3_tib4", "rem_tib4", "wake_tib4"
)

pred_mean_tibq4_isig_sleep <- list()
for (i in 0:1) {
  pred <- pred_tibq4_isig_sleep_draws[isig == i]
  pred <- as.data.frame(pred[, ..parts_mean])
  pred <- apply(pred, 2, as.numeric)
  pred <- apply(pred, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.95))
  pred <- Map(cbind, pred, part = names(pred))
  pred <- rbindlist(pred)
  pred[, isig := i]

  pred_mean_tibq4_isig_sleep[[i + 1]] <- pred
}
pred_mean_tibq4_isig_sleep <- rbindlist(pred_mean_tibq4_isig_sleep)

# constrast
parts_contrast <- c(
  # "sol_isi1_isi0", "waso_isi1_isi0", "light_isi1_isi0", "stage3_isi1_isi0", "rem_isi1_isi0", "wake_isi1_isi0",
  "sol_q23_q1", "waso_q23_q1", "light_q23_q1", "stage3_q23_q1", "rem_q23_q1", "wake_q23_q1",
  "sol_q23_q4", "waso_q23_q4", "light_q23_q4", "stage3_q23_q4", "rem_q23_q4", "wake_q23_q4",
  "sol_q1_q4", "waso_q1_q4", "light_q1_q4", "stage3_q1_q4", "rem_q1_q4", "wake_q1_q4"
)
pred_contrast_tibq4_isig_sleep <- list()
for (i in 0:1) {
  pred <- pred_tibq4_isig_sleep_draws[isig == i]
  pred <- as.data.frame(pred[, ..parts_contrast])
  pred <- apply(pred, 2, as.numeric)
  pred <- apply(pred, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.95))
  pred <- Map(cbind, pred, part = names(pred))
  pred <- rbindlist(pred)
  pred[, isig := i]
  pred_contrast_tibq4_isig_sleep[[i + 1]] <- pred
}
pred_contrast_tibq4_isig_sleep <- rbindlist(pred_contrast_tibq4_isig_sleep)

d_plot_comp_tibq4_isig_sleep <- rbind(
  pred_mean_tibq4_isig_sleep,
  pred_contrast_tibq4_isig_sleep
)
table(d_plot_comp_tibq4_isig_sleep$part)

d_plot_comp_tibq4_isig_sleep[, part_label := NA]
d_plot_comp_tibq4_isig_sleep[, part_label := ifelse(grepl("sol", part), "Sleep onset latency", part_label)]
d_plot_comp_tibq4_isig_sleep[, part_label := ifelse(grepl("waso", part), "Wake after sleep onset", part_label)]
d_plot_comp_tibq4_isig_sleep[, part_label := ifelse(grepl("light", part), "Light sleep", part_label)]
d_plot_comp_tibq4_isig_sleep[, part_label := ifelse(grepl("stage3", part), "Slow wave sleep", part_label)]
d_plot_comp_tibq4_isig_sleep[, part_label := ifelse(grepl("rem", part), "REM sleep", part_label)]
d_plot_comp_tibq4_isig_sleep[, part_label := ifelse(grepl("wake", part), "Daytime wake", part_label)]

d_plot_comp_tibq4_isig_sleep[, part_label := factor(part_label, ordered = TRUE, levels = c(
  "Daytime wake",
  "Sleep onset latency", "Wake after sleep onset", "Light sleep",
  "Slow wave sleep", "REM sleep"
))]
table(d_plot_comp_tibq4_isig_sleep$part_label, useNA = "always")

d_plot_comp_tibq4_isig_sleep[, tib_group := NA]
d_plot_comp_tibq4_isig_sleep[, tib_group := ifelse(grepl("tib1", part), "q1", tib_group)]
d_plot_comp_tibq4_isig_sleep[, tib_group := ifelse(grepl("tib2", part), "q2", tib_group)]
d_plot_comp_tibq4_isig_sleep[, tib_group := ifelse(grepl("tib3", part), "q3", tib_group)]
d_plot_comp_tibq4_isig_sleep[, tib_group := ifelse(grepl("tib4", part), "q4", tib_group)]
d_plot_comp_tibq4_isig_sleep[, tib_group := ifelse(grepl("tib5", part), "q5", tib_group)]
table(d_plot_comp_tibq4_isig_sleep$tib_group, useNA = "always")

d_plot_comp_tibq4_isig_sleep[, isi_group := ifelse(isig == 1, "Insomnia+", "Insomnia-")]

d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := NA]
d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := ifelse(part_label == "Sleep onset latency", d_plot_comp_tibq4_isig_sleep[isig == 0 & part == "sol_tib"]$Mean, yintercept_healthy)]
d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := ifelse(part_label == "Wake after sleep onset", d_plot_comp_tibq4_isig_sleep[isig == 0 & part == "waso_tib"]$Mean, yintercept_healthy)]
d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := ifelse(part_label == "Light sleep", d_plot_comp_tibq4_isig_sleep[isig == 0 & part == "light_tib"]$Mean, yintercept_healthy)]
d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := ifelse(part_label == "Slow wave sleep", d_plot_comp_tibq4_isig_sleep[isig == 0 & part == "stage3_tib"]$Mean, yintercept_healthy)]
d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := ifelse(part_label == "REM sleep", d_plot_comp_tibq4_isig_sleep[isig == 0 & part == "rem_tib"]$Mean, yintercept_healthy)]
d_plot_comp_tibq4_isig_sleep[, yintercept_healthy := ifelse(part_label == "Daytime wake", d_plot_comp_tibq4_isig_sleep[isig == 0 & part == "wake_tib"]$Mean, yintercept_healthy)]

### plot -----------------------
(plot_comp_tibq4_isig_sleep <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group)], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
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
  #  geom_text(aes(y = 627.5, label = TeX(est_sig, output = "character")), parse = TRUE,
  #            hjust = 0.5, nudge_x = 0,
  #            family = "Arial Narrow", size = 3,
  #            show.legend = FALSE) +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)

grDevices::cairo_pdf(
  file = paste0(out, "plot_comp_tibq4_isig_sleep", ".pdf"),
  width = 7,
  height = 10,
)
plot_comp_tibq4_isig_sleep
dev.off()

(plot_comp_tibq4_isig_sol <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "Sleep onset latency"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 800),
    breaks = c(0, 800),
    name = "Sleep onset latency (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)
(plot_comp_tibq4_isig_wake <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "Daytime wake"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(350, 1400),
    breaks = c(400, 1400),
    name = "Daytime wake (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)
(plot_comp_tibq4_isig_waso <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "Wake after sleep onset"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 150),
    breaks = c(0, 150),
    name = "Wake after sleep onset (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)
(plot_comp_tibq4_isig_wake <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "Daytime wake"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(350, 1400),
    breaks = c(400, 1400),
    name = "Daytime wake (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)
(plot_comp_tibq4_isig_light <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "Light sleep"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 350),
    breaks = c(0, 350),
    name = "Light sleep (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)
(plot_comp_tibq4_isig_sws <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "Slow wave sleep"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 200),
    breaks = c(0, 200),
    name = "Slow wave sleep (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)
(plot_comp_tibq4_isig_rem <-
  ggplot(d_plot_comp_tibq4_isig_sleep[!is.na(tib_group) & part_label == "REM sleep"], aes(x = tib_group, y = Mean, group = isi_group, colour = interaction(tib_group, isi_group))) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
  #  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
  #  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
  geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
  geom_pointrange(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    size = .25,
    linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_text(aes(y = Mean + 1.5, label = TeX(sig_contrast, output = "character")),
    parse = TRUE,
    hjust = 0.5, nudge_x = .15,
    family = "Arial Narrow", size = 4,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 150),
    breaks = c(0, 150),
    name = "REM sleep (min)"
  ) +
  # facet_wrap(~part_label, ncol = 2, scales = "free") +
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
    axis.text.x         = element_text(size = 10),
    axis.text.y         = element_text(size = 10),
    strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
    legend.text         = element_text(size = 10, face = "bold", hjust = .5),
    legend.position     = "none",
    plot.margin         = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
  )
)

grDevices::cairo_pdf(
  file = paste0(out, "plot_comp_tibq4_isig_all", ".pdf"),
  width = 7,
  height = 10,
)
ggarrange(plot_comp_tibq4_isig_wake, plot_comp_tibq4_isig_sol,
  plot_comp_tibq4_isig_waso, plot_comp_tibq4_isig_light,
  plot_comp_tibq4_isig_sws, plot_comp_tibq4_isig_rem,
  ncol = 2, nrow = 3, common.legend = TRUE, legend = "none"
)
dev.off()
