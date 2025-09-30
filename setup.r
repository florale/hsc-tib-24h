
base <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/monash/projects/healthy-sleep-clinic/data")
out <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/monash/projects/healthy-sleep-clinic/hsc-tib-24h/output")
redir <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/github/projects/hsc/")
source(paste0(redir, "setup.r"))

scales::show_col(tvthemes:::hilda_palette$Day)
scales::show_col(tvthemes:::hilda_palette$Night)
scales::show_col(tvthemes:::hilda_palette$Dusk)
library(extrafont)
library(colorspace)
col <- c(
  "#4F5D5B", 
  "#708885", 
  # "#3d251e", 
  "#978787",
  "#EAD3BF", 
  "#c48462"
)

col_sex <- c(
  darken("#4F5D5B", 0.3),
  darken("#708885", 0.3), 
  darken("#978787", 0.3),
  darken("#EAD3BF", 0.3),
  darken("#c48462", 0.3), 

  lighten("#4F5D5B", 0.3), 
  lighten("#708885", 0.3),
  lighten("#978787", 0.3), 
  lighten("#EAD3BF", 0.3), 
  lighten("#c48462", 0.3)
)

col_age <- c(
  darken("#4F5D5B", 0.3),
  darken("#708885", 0.3), 
  darken("#978787", 0.3),
  darken("#EAD3BF", 0.3),
  darken("#c48462", 0.3), 

  lighten("#4F5D5B", 0.3), 
  lighten("#708885", 0.3),
  lighten("#978787", 0.3), 
  lighten("#EAD3BF", 0.3), 
  lighten("#c48462", 0.3)
)
col_serv <- col_age

shape_perc <- c(21)
shape_sex <- c("Female" = 24, "Male" = 25)
shape_age <- c("u45" = "\u25c1", "45u" = "\u25b7")

shape_min_sex <- c("Female" = "\u25b3", "Male" = "\u25bd")
shape_min_age <- c("u45" = "\u25C0", "45u" = "\u25B6")
shape_min_serv <- c("home" = "\u25d2", "hospital" = "\u25d3")

shape_perc_sex <- c("Female" = 24, "Male" = 25)
  

# col <- c(
#   # "#F6E0D2",
#   "#c48462",
#   "#9C6755",
#   "#659794",
#   "#586085",
#   "#F5C98E" # B87474
# )
# colf <- c(
#   `REM` = "#DFA398",
#   `N3` = "#9C6755",
#   `N1+2` = "#659794",
#   `WASO` = "#586085", # #586085
#   `SOL` = "#F5C98E",
#   `WD` = "#F6E0D2"
# )

protan <- dichromat::dichromat(col, type = "protan")
deutan <- dichromat::dichromat(col, type = "deutan")
tritan <- dichromat::dichromat(col, type = "tritan")

# plot for comparison
layout(matrix(1:4, nrow = 4)); par(mar = rep(1, 4))
recolorize::plotColorPalette(col, main = "Trichromacy")
recolorize::plotColorPalette(protan, main = "Protanopia")
recolorize::plotColorPalette(deutan, main = "Deutanopia")
recolorize::plotColorPalette(tritan, main = "Tritanopia")

make_min_plot <- function(part_label) {
  params <- plot_min_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isi1_sex_sum[!is.na(tib_group) & part_label == part & grepl("min", par) & mean == 1],
    aes(x = tib_group, y = Mean, group = female, colour = interaction(tib_group, female))
  ) +
    # geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high,
        shape = female
      )
        , size = 1
        , linewidth = 1
        , position = position_dodge(width = 1)
      ) +
    # geom_text(aes(y = Mean + params$y_offset, label = latex2exp::TeX(contrast_min_sig, output = "character")),
    #   parse = TRUE,
    #   hjust = 0.5, nudge_x = .2,
    #   family = "Arial Narrow",
    #   size = 7,
    #   show.legend = FALSE
    # ) +
    geom_text(aes(y = max(params$limits), label = est_min),
      hjust = 1, 
      position = position_dodge(width = 1),
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
    scale_shape_manual(values = c("Female" = "\u25BC", "Male" = "\u25B2")) +
    scale_colour_manual(values = col_sex) +
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
      axis.text.x         = element_text(size = 13, face = "plain", family = "Arial Narrow"),
      axis.text.y         = element_text(size = 13, face = "plain", family = "Arial Narrow", margin = margin(l = 20)),
      strip.text          = element_text(size = 13, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(2, 1.5, 1.5, -2), "lines")
    )
}

make_perc_plot <- function(part_label) {
  params <- plot_perc_params[[part_label]]
  part <- part_label
  ggplot(
    pred_tibq5_isi1_sex_sum[!is.na(tib_group) & part_label == part & grepl("perc", par) & mean == 1],
    aes(x = tib_group, y = Mean, group = female, colour = interaction(tib_group, female))
  ) +
    # geom_hline(aes(yintercept = yintercept_insom), linewidth = 0.75, linetype = "dashed", colour = "#DCD5CE") +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.75, linetype = "dashed", colour = "#A9A9A9") +
    geom_pointrange(
      aes(
        ymin = CI_low,
        ymax = CI_high,
        shape = female
      ),
      size = 1,
      linewidth = 1,
      position = position_dodge(width = 1)
    ) +
    # geom_text(aes(y = Mean + params$y_offset, label = latex2exp::TeX(contrast_perc_sig, output = "character")),
    #   parse = TRUE,
    #   hjust = 0.5, nudge_x = .2,
    #   family = "Arial Narrow",
    #   size = 7,
    #   show.legend = FALSE
    # ) +
    geom_text(aes(y = max(params$limits), label = est_perc),
      hjust = 1,
      position = position_dodge(width = 1),
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
    scale_shape_manual(values = c("Female" = "\u25BC", "Male" = "\u25B2")) +
    scale_colour_manual(values = col_sex) +
    scale_fill_manual(values = col_sex) +
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
      axis.text.x         = element_text(size = 13, face = "plain", family = "Arial Narrow"),
      axis.text.y         = element_text(size = 13, face = "plain", family = "Arial Narrow", margin = margin(l = 20)),
      strip.text          = element_text(size = 13, face = "plain", family = "Arial Narrow", hjust = .5),
      axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      legend.text         = element_blank(),
      legend.position     = "none",
      plot.margin         = unit(c(2, 1, 1.5, -2), "lines")
    )
}