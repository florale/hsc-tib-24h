source("data.r")
parts <- c("wake_min", "sol_min", "light_min", "stage3_min", "rem_min", "waso_min")

dpsg[, age45 := ifelse(age >= 45, "45u", "u45")]
dpsg[, age45 := factor(age45, levels = c("u45", "45u"))]

# moderation ---------------------------------------------------------------
clr_isi8 <- complr(
  data = dpsg[isi > 7],
  parts = parts,
  # idvar = "record_id",
  total = 1440
)
clr_isi8$dataout[, so_min := 1440 - twake_min]

m_tib_isi8_age <- brmcoda(clr_isi8,
  mvbind(z1_1, z2_1, z3_1, z4_1, z5_1) ~
    s(so_min, by = age45) +
    female + bmi + white + working + labpsg + s(perHrAHSleep) + antidep + nap,
  iter = 4000, chains = 6, cores = 6, seed = 123, warmup = 1000,
  backend = "cmdstanr"
)
summary(m_tib_isi8_age)
saveRDS(m_tib_isi8_age, file.path(out, "m_tib_isi8_age.rds"))
m_tib_isi8_age <- readRDS(file.path(out, "m_tib_isi8_age.rds"))

d_tib_isi8_age <- emmeans::ref_grid(m_tib_isi8_age$model,
  at = list(
    so_min = quantile(model.frame(m_tib_isi8_age)$so_min, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE),
    age45 = c("u45", "45u")
  )
)@grid
d_tib_isi8_age <- as.data.table(d_tib_isi8_age)
d_tib_isi8_age <- d_tib_isi8_age[!duplicated(d_tib_isi8_age[, .(so_min, age45, female, white, working, antidep, labpsg)]), ]

pred_tib_isi8_age <- fitted(m_tib_isi8_age, newdata = d_tib_isi8_age, scale = "response", re_formula = NA, summary = FALSE)
pred_tib_isi8_age_draws <- apply(pred_tib_isi8_age, 1, function(d) {
  d <- cbind(d_tib_isi8_age, d)
  d <- as.data.table(d)

  setnames(d, paste0("t", parts), parts)
  d[, tib_p := factor(so_min, labels = c("p10", "p25", "p50", "p75", "p90"))]
  d[, isi_g := "isi8"]
  d[, age := age45][, age45 := NULL]

  d[, tst_min := light_min + stage3_min + rem_min]

  d[, light_perc := (light_min / (light_min + stage3_min + rem_min)) * 100]
  d[, stage3_perc := (stage3_min / (light_min + stage3_min + rem_min)) * 100]
  d[, rem_perc := (rem_min / (light_min + stage3_min + rem_min)) * 100]

  # calculated weighted means by .wgt. for min
  d[, sol_min := weighted.mean(sol_min, .wgt.), by = .(isi_g, tib_p, age)]
  d[, waso_min := weighted.mean(waso_min, .wgt.), by = .(isi_g, tib_p, age)]
  d[, light_min := weighted.mean(light_min, .wgt.), by = .(isi_g, tib_p, age)]
  d[, stage3_min := weighted.mean(stage3_min, .wgt.), by = .(isi_g, tib_p, age)]
  d[, rem_min := weighted.mean(rem_min, .wgt.), by = .(isi_g, tib_p, age)]
  d[, wake_min := weighted.mean(wake_min, .wgt.), by = .(isi_g, tib_p, age)]
  d[, tst_min := weighted.mean(tst_min, .wgt.), by = .(isi_g, tib_p, age)]

  # calculated weighted means by .wgt. for all perc
  d[, light_perc := weighted.mean(light_perc, .wgt.), by = .(isi_g, tib_p, age)]
  d[, stage3_perc := weighted.mean(stage3_perc, .wgt.), by = .(isi_g, tib_p, age)]
  d[, rem_perc := weighted.mean(rem_perc, .wgt.), by = .(isi_g, tib_p, age)]

  # make wide
  d <- dcast(d, . ~ isi_g + tib_p + age,
    fun.aggregate = mean,
    value.var = c(
      "wake_min", "sol_min", "waso_min", "light_min", "stage3_min", "rem_min", "tst_min",
      "light_perc", "stage3_perc", "rem_perc"
    )
  )
  d <- d[, -1, with = FALSE] # remove the dot column

  for (age in c("45u", "u45")) {
    # grand mean min and percentages
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
      d[, paste0(part, "_min_mean", "_", age) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = paste0(part, "_min_isi", rep(8, 5), "_p", rep(c(10, 25, 50, 75, 90), 2), "_", age)
      ]
    }
    for (part in c("light", "stage3", "rem")) {
      d[, paste0(part, "_perc_mean", "_", age) := rowMeans(.SD, na.rm = TRUE),
        .SDcols = paste0(part, "_perc_isi", rep(8, 3), "_p", rep(c(10, 25, 50, 75, 90), 2), "_", age)
      ]
    }

    # difference between isi8, p10 vs p2575
    for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
      d[, paste0(part, "_min_isi8_p10_p2575", "_", age) :=
        get(paste0(part, "_min_isi8_p10", "_", age))
        - ((get(paste0(part, "_min_isi8_p25", "_", age)) + get(paste0(part, "_min_isi8_p50", "_", age)) + get(paste0(part, "_min_isi8_p75", "_", age))) / 3)]

      d[, paste0(part, "_min_isi8_p90_p2575", "_", age) :=
        get(paste0(part, "_min_isi8_p90", "_", age))
        - ((get(paste0(part, "_min_isi8_p25", "_", age)) + get(paste0(part, "_min_isi8_p50", "_", age)) + get(paste0(part, "_min_isi8_p75", "_", age))) / 3)]

      d[, paste0(part, "_min_isi8_p10_p90", "_", age) :=
        get(paste0(part, "_min_isi8_p10", "_", age))
        - get(paste0(part, "_min_isi8_p90", "_", age))]
    }
    for (part in c("light", "stage3", "rem")) {
      d[, paste0(part, "_perc_isi8_p10_p2575", "_", age) :=
        get(paste0(part, "_perc_isi8_p10", "_", age))
        - ((get(paste0(part, "_perc_isi8_p25", "_", age)) + get(paste0(part, "_perc_isi8_p50", "_", age)) + get(paste0(part, "_perc_isi8_p75", "_", age))) / 3)]

      d[, paste0(part, "_perc_isi8_p90_p2575", "_", age) :=
        get(paste0(part, "_perc_isi8_p90", "_", age))
        - ((get(paste0(part, "_perc_isi8_p25", "_", age)) + get(paste0(part, "_perc_isi8_p50", "_", age)) + get(paste0(part, "_perc_isi8_p75", "_", age))) / 3)]
      d[, paste0(part, "_perc_isi8_p10_p90", "_", age) :=
        get(paste0(part, "_perc_isi8_p10", "_", age))
        - get(paste0(part, "_perc_isi8_p90", "_", age))]
    }
  }

  # difference in min an perc between sex
  for (part in c("sol", "waso", "light", "stage3", "rem", "wake", "tst")) {
    for (p in c(10, 25, 50, 75, 90)) {
      d[, paste0(part, "_min_isi8_p", p, "_45u_u45") :=
        get(paste0(part, "_min_isi8_p", p, "_u45")) -
        get(paste0(part, "_min_isi8_p", p, "_45u"))]
    }
    d[, paste0(part, "_min_mean", "_45u_u45") :=
      get(paste0(part, "_min_mean", "_u45")) -
      get(paste0(part, "_min_mean", "_45u"))]
  }
  for (part in c("light", "stage3", "rem")) {
    for (p in c(10, 25, 50, 75, 90)) {
      d[, paste0(part, "_perc_isi8_p", p, "_45u_u45") :=
        get(paste0(part, "_perc_isi8_p", p, "_u45")) -
        get(paste0(part, "_perc_isi8_p", p, "_45u"))]
    }
    d[, paste0(part, "_perc_mean", "_45u_u45") :=
      get(paste0(part, "_perc_mean", "_u45")) -
      get(paste0(part, "_perc_mean", "_45u"))]
  }
  d
})
saveRDS(pred_tib_isi8_age_draws, file.path(out, "pred_tib_isi8_age_draws.rds"))

# prep final data
pred_tib_isi8_age_draws <- readRDS(file.path(out, "pred_tib_isi8_age_draws.rds"))

pred_tib_isi8_age_draws <- as.data.table(abind(pred_tib_isi8_age_draws, along = 1))
# pred_tib_isi8_draws <- split(pred_tib_isi8_draws, pred_tib_isi8_draws$so_min)

pred_tib_isi8_age_sum <- apply(pred_tib_isi8_age_draws, 2, function(x) {
  describe_posterior(as.numeric(x), centrality = "mean", ci = 0.95)
})
pred_tib_isi8_age_sum <- rbindlist(pred_tib_isi8_age_sum)

pred_tib_isi8_age_sum[, par := colnames(pred_tib_isi8_age_draws)]

pred_tib_isi8_age_sum[, part_label := NA]
pred_tib_isi8_age_sum[, part_label := ifelse(grepl("sol", par) & grepl("min|perc", par), "SOL", part_label)]
pred_tib_isi8_age_sum[, part_label := ifelse(grepl("waso", par) & grepl("min|perc", par), "WASO", part_label)]
pred_tib_isi8_age_sum[, part_label := ifelse(grepl("light", par) & grepl("min|perc", par), "N1+2", part_label)]
pred_tib_isi8_age_sum[, part_label := ifelse(grepl("stage3", par) & grepl("min|perc", par), "N3", part_label)]
pred_tib_isi8_age_sum[, part_label := ifelse(grepl("rem", par) & grepl("min|perc", par), "REM", part_label)]
pred_tib_isi8_age_sum[, part_label := ifelse(grepl("wake", par) & grepl("min|perc", par), "WD", part_label)]

pred_tib_isi8_age_sum[, part_label := factor(part_label, ordered = TRUE, levels = c(
  "WD",
  "SOL",
  "N1+2",
  "N3",
  "REM",
  "WASO"
))]
table(pred_tib_isi8_age_sum$part_label, useNA = "always")

pred_tib_isi8_age_sum[, tib_group := NA]
pred_tib_isi8_age_sum[, tib_group := ifelse(grepl("isi8_p100", par), "P10 SO", tib_group)]
pred_tib_isi8_age_sum[, tib_group := ifelse(grepl("isi8_p25", par), "P25 SO", tib_group)]
pred_tib_isi8_age_sum[, tib_group := ifelse(grepl("isi8_p50", par), "P50 SO", tib_group)]
pred_tib_isi8_age_sum[, tib_group := ifelse(grepl("isi8_p75", par), "P75 SO", tib_group)]
pred_tib_isi8_age_sum[, tib_group := ifelse(grepl("isi8_p90", par), "P90 SO", tib_group)]

table(pred_tib_isi8_age_sum$tib_group, useNA = "always")

pred_tib_isi8_age_sum[, age := NA]
pred_tib_isi8_age_sum[, age := ifelse(grepl("_u45", par) & !grepl("_45u_u45", par), "u45", age)]
pred_tib_isi8_age_sum[, age := ifelse(grepl("_45u", par) & !grepl("_45u_u45", par), "45u", age)]

# types of estimates
table(pred_tib_isi8_age_sum$par, useNA = "always")
pred_tib_isi8_age_sum[, contrast := NA]
pred_tib_isi8_age_sum[, contrast := ifelse(grepl("p10_p90", par), 1, contrast)]
pred_tib_isi8_age_sum[, contrast := ifelse(grepl("p10_p2575", par), 1, contrast)]
pred_tib_isi8_age_sum[, contrast := ifelse(grepl("p90_p2575", par), 1, contrast)]
pred_tib_isi8_age_sum[, contrast := ifelse(is.na(contrast), 0, contrast)]
pred_tib_isi8_age_sum[, mean := ifelse(contrast == 0, 1, 0)]

table(pred_tib_isi8_age_sum$mean, useNA = "always")
table(pred_tib_isi8_age_sum$contrast, useNA = "always")
# table(pred_tib_isi8_age_sum$perc, useNA = "always")

pred_tib_isi8_age_sum[, sig := ifelse(!between(0, CI_low, CI_high), "$^a$", "$\\phantom{^a}$")]
table(pred_tib_isi8_age_sum$sig, useNA = "always")

pred_tib_isi8_age_sum[grepl("min", par) & !grepl("perc", par), est_min := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
pred_tib_isi8_age_sum[grepl("perc", par), est_perc := paste0(format(round(Mean, 2), nsmall = 2), " [", format(round(CI_low, 2), nsmall = 2), ", ", format(round(CI_high, 2), nsmall = 2), "]")]

saveRDS(pred_tib_isi8_age_sum, file.path(out, "pred_tib_isi8_age_sum.rds"))
pred_tib_isi8_age_sum <- readRDS(file.path(out, "pred_tib_isi8_age_sum.rds"))

### plot -----------------------
## plot_min -----------------------
# make individual then patch
plot_min_params <- list(
  "WD" = list(
    limits = c(800, 1350), breaks = c(1000, 1100, 1200), breaks2 = c(.60, .70, .80), name = "WD", y_offset = 1.5
  ),
  "SOL" = list(
    limits = c(0, 50), breaks = c(10, 30, 50), breaks2 = c(.01, .04, .07), name = "SOL", y_offset = 1.5
  ),
  "WASO" = list(
    limits = c(0, 175), breaks = c(25, 50, 75), breaks2 = c(.03, .05, .07, .09), name = "WASO", y_offset = 1.5
  ),
  "N1+2" = list(
    limits = c(130, 390), breaks = c(200, 250, 300), breaks2 = c(.15, .20, .25), name = "N1+2", y_offset = 1.5
  ),
  "N3" = list(
    limits = c(20, 100), breaks = c(50, 75, 100), breaks2 = c(.03, .05, .07), name = "N3", y_offset = 1.5
  ),
  "REM" = list(
    limits = c(20, 110), breaks = c(50, 75, 100), breaks2 = c(.03, .05, .07), name = "REM", y_offset = 1.5
  )
)

make_min_plot <- function(part_label) {
  params <- plot_min_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tib_isi8_age_sum[!is.na(tib_group) & !is.na(age) & part_label == part & grepl("min", par) & mean == 1],
    aes(x = tib_group, y = Mean, group = age, shape = age, colour = interaction(tib_group, age), fill = interaction(tib_group, age))
  ) +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      size = 1.5,
      linewidth = 0.75,
      position = position_dodge(width = 0.75)
    ) +
    geom_text(aes(y = min(params$limits), label = tib_group),
      hjust = 0, nudge_x = 0,
      family = "Arial Narrow",
      fontface = "plain",
      size = 5,
      colour = "black",
      show.legend = FALSE
    ) +
    geom_text(aes(y = max(params$limits), label = est_min),
      hjust = 1,
      position = position_dodge(width = 0.75),
      family = "Arial Narrow",
      fontface = "bold",
      size = 5,
      show.legend = FALSE
    ) +
    scale_y_continuous(
      limits = params$limits,
      breaks = params$breaks,
      labels = paste0(params$breaks),
      position = "right"
    ) +
    scale_colour_manual(values = col_age) +
    scale_fill_manual(values = col_age) +
    scale_shape_manual(values = shape_min_age) +
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
      axis.text.x         = element_blank(),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 16, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(1.5, 1.5, 0, 0), "lines")
    )
}
plots_min <- lapply(names(plot_min_params), make_min_plot)
names(plots_min) <- names(plot_min_params)

plot_age_min_isi8 <- ggarrange(
  plots_min[["WD"]],
  plots_min[["SOL"]],
  plots_min[["WASO"]],
  plots_min[["N1+2"]],
  plots_min[["N3"]],
  plots_min[["REM"]],
  labels = c(
    "  Daytime Wake (min)",
    "  Sleep Onset Latency (min)",
    "  Wake After Sleep Onset (min)",
    "  N1+2 (min)",
    "  N3 (min)",
    "  REM (min)"
  ),
  hjust = 0,
  ncol = 3, nrow = 2, common.legend = TRUE, legend = "none",
  font.label = list(size = 16, face = "italic", family = "Arial Narrow")
) + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
ggsave(file.path(out, paste0("plot_age_min_isi8", ".pdf")), plot_age_min_isi8, device = cairo_pdf, width = 12, height = 8, dpi = 300)

## plot_perc -----------------------
plot_perc_params <- list(
  # "WASO" = list(
  #   limits = c(5, 25), breaks = c(5, 10, 15, 20), name = "WASO", y_offset = 0.005
  # ),
  "N1+2" = list(
    limits = c(55, 90), breaks = c(55, 65, 75), name = "N1+2", y_offset = 0.01
  ),
  "N3" = list(
    limits = c(5, 40), breaks = c(15, 25, 35), name = "N3", y_offset = 0.005
  ),
  "REM" = list(
    limits = c(5, 40), breaks = c(15, 25, 35), name = "REM", y_offset = 0.005
  )
)
make_perc_plot <- function(part_label) {
  params <- plot_perc_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tib_isi8_age_sum[!is.na(tib_group) & !is.na(age) & part_label == part & grepl("perc", par) & mean == 1],
    aes(x = tib_group, y = Mean, group = age, shape = age, colour = interaction(tib_group, age))
  ) +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      size = 1.5,
      linewidth = 0.75,
      position = position_dodge(width = 0.75)
    ) +
    geom_text(aes(y = min(params$limits), label = tib_group),
      hjust = 0, nudge_x = 0,
      family = "Arial Narrow",
      fontface = "plain",
      size = 5,
      colour = "black",
      show.legend = FALSE
    ) +
    geom_text(aes(y = max(params$limits), label = est_perc),
      hjust = 1,
      position = position_dodge(width = 0.75),
      family = "Arial Narrow",
      fontface = "bold",
      size = 5,
      show.legend = FALSE
    ) +
    scale_y_continuous(
      limits = params$limits,
      breaks = params$breaks,
      labels = paste0(params$breaks),
      position = "right"
    ) +
    scale_colour_manual(values = col_age) +
    # scale_fill_manual(values = col_age) +
    scale_shape_manual(values = shape_age) +
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
      axis.text.x         = element_blank(),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 16, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      # axis.line.x         = element_blank(),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(1.5, 1.5, 0, 0), "lines")
    )
}
plots_perc <- lapply(names(plot_perc_params), make_perc_plot)
names(plots_perc) <- names(plot_perc_params)

plot_age_perc_isi8 <- ggarrange(
  # plots_perc[["WASO"]],
  plots_perc[["N1+2"]],
  plots_perc[["N3"]],
  plots_perc[["REM"]],
  labels = c(
    # "  A. Wake After Sleep Onset (%)",
    "  N1+2 (%)",
    "  N3 (%)",
    "  REM (%)"
  ),
  hjust = 0,
  ncol = 3, nrow = 1, common.legend = TRUE, legend = "none",
  font.label = list(size = 16, face = "italic", family = "Arial Narrow")
) + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
ggsave(file.path(out, paste0("plot_age_perc_isi8", ".pdf")), plot_age_perc_isi8, device = cairo_pdf, width = 12, height = 4, dpi = 300)

# patch into one figure -----------------------------------------------
subtitle_plot <- function(label, margin_cm = c(0, 0, 2, 0)) {
  ggplot() +
    labs(title = label) +
    theme_void() +
    theme(
      plot.title = element_text(
        family = "Arial Narrow",
        face   = "bold",
        size   = 16,
        hjust  = 0 # 0 = left, 0.5 = centre, 1 = right
      )
    )
}

plot_age_min_perc_isi8 <-
  subtitle_plot(" A. 24h Sleep-Wake Architecture Composition") / plot_age_min_isi8 /
  subtitle_plot(" B. Sleep Architecture Composition") / plot_age_perc_isi8 +
  plot_layout(heights = c(0, 2, 0, 1))

ggsave(file.path(out, paste0("plot_age_min_perc_isi8", ".pdf")), plot_age_min_perc_isi8, device = cairo_pdf, width = 12, height = 12, dpi = 300)
saveRDS(plot_age_min_perc_isi8, file.path(out, paste0("plot_age_min_perc_isi8", ".rds")))