source("setup.R")

col <- c(
  `REM` = "#708885",
  `N3` = "#4F5D5B",
  `N1+2` = "#9DB3A8",
  `WASO` = "#D2BDAB",
  `SOL` = "#c48462", # 7083A4
  `WD` = "#FBE1B1"
)

family <- "Arial Narrow"

day <- c("WD", "SOL", "N1+2", "N3", "REM", "WASO")

composition_a <- data.frame(
 "day" = day,
 "prop" = c(14 / 24, 1.5 / 24, 3.5 / 24, 2 / 24, 2 / 24, 1 / 24),
 "prop_labels" = c("16h", "1h", "4h", "1.5h", "1.5h", "0.5h"),
 label = "A. Reference composition"
)
composition_b <- data.frame(
 "day" = day,
 "prop" = c(14 / 24, 1 / 24, 4 / 24, 2 / 24, 2 / 24, 1 / 24),
 "prop_labels" = c("16h", "0.5h", "3.5h", "1.5h", "1.5h", "0.5h"),
#  label = paste0("B. ", "30-min SOL ", "\u2192", " N1+2")
 label = paste0("B. 30-min from SOL to N1+2")

)
composition_c <- data.frame(
 "day" = day,
 "prop" = c(14.5 / 24, 1 / 24, 3.5 / 24, 2 / 24, 2 / 24, 1 / 24),
 "prop_labels" = c("16.5h", "0.5h", "4h", "1.5h", "1.5h", "0.5h"),
#  label = paste0("C. ", "30-min SOL ", "\u2192", " WD")
  label = paste0("C. 30-min from SOL to WD")
)

composition_d <- data.frame(
 "day" = day,
 "prop" = c(14 / 24, 1.5 / 24, 3.5 / 24, 2 / 24, 2 / 24, 1 / 24),
 "prop_labels" = c("?", "?", "?", "?", "?", "?"),
 label = "D. Redistribution"
)

composition_e <- data.frame(
 "day" = day,
 "prop" = c(14 / 24, 1.5 / 24, 3.5 / 24, 2 / 24, 2 / 24, 1 / 24),
 "prop_labels" = c("", "", "", "", "", ""),
 label = "E."
)

composition <- rbind(
 composition_a,
 composition_b,
 composition_c
)

composition$day <- factor(composition$day,
 ordered = TRUE,
 levels = c(
  "WD",
  "SOL",
  "N1+2",
  "N3",
  "REM",
  "WASO"
 )
)

grDevices::cairo_pdf(
 file = file.path(out, "hsc-tib-24h-media.pdf"),
 width = 8.5,
 height = 4,
)
ggplot(composition, aes(x = 2, y = prop, fill = day)) +
 geom_bar(stat = "identity", width = 4, linewidth = 0.25, colour = "white") +
 coord_polar("y", start = 130, direction = -1) +
 geom_text(
  aes(x = 2.75, label = prop_labels),
  position = position_stack(vjust = 0.5),
  size = 3.5,
  family = family
 ) +
 facet_wrap(~label) +
 scale_fill_manual(values = scales::alpha(col, 0.75), name = "Time in") +
 guides(fill = guide_legend(nrow = 1)) +
 xlim(-3, 4) +
 hrbrthemes::theme_ipsum() +
 theme(
  panel.background = element_blank(),
  panel.border     = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # panel.spacing    = unit(0.5, "lines"),
  plot.background  = element_rect(fill = "transparent", colour = NA),
  strip.background = element_rect(fill = "transparent", colour = NA),
  axis.title.x     = element_blank(),
  axis.title.y     = element_blank(),
  axis.text.x     = element_blank(),
  axis.text.y     = element_blank(),
  strip.text      = element_text(size = 12, hjust = 0.5, face = "bold"),
  legend.title    = element_text(size = 12, face = "bold"),
  legend.text     = element_text(size = 12),
  plot.margin     = margin(0.5, 0, 0.5, 0, "cm"),
  legend.position = "bottom"
 )
dev.off()

composition <- rbind(
 composition_a,
 composition_d,
 composition_e
)

composition$day <- factor(composition$day,
 ordered = TRUE,
 levels = c(
  "WD",
  "SOL",
  "N1+2",
  "N3",
  "REM",
  "WASO"
 )
)

grDevices::cairo_pdf(
 file = file.path(out, "hsc-tib-24h-media-redistribution.pdf"),
 width = 8.5,
 height = 4,
)
ggplot(composition, aes(x = 2, y = prop, fill = day)) +
 geom_bar(stat = "identity", width = 4, linewidth = 0.25, colour = "white") +
 coord_polar("y", start = 130, direction = -1) +
 geom_text(
  aes(x = 2.75, label = prop_labels),
  position = position_stack(vjust = 0.5),
  size = 3.5,
  family = family
 ) +
 facet_wrap(~label) +
 scale_fill_manual(values = scales::alpha(col, 0.75), name = "Time in") +
 guides(fill = guide_legend(nrow = 1)) +
 xlim(-3, 4) +
 hrbrthemes::theme_ipsum() +
 theme(
  panel.background = element_blank(),
  panel.border     = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # panel.spacing    = unit(0.5, "lines"),
  plot.background  = element_rect(fill = "transparent", colour = NA),
  strip.background = element_rect(fill = "transparent", colour = NA),
  axis.title.x     = element_blank(),
  axis.title.y     = element_blank(),
  axis.text.x     = element_blank(),
  axis.text.y     = element_blank(),
  strip.text      = element_text(size = 12, hjust = 0.5, face = "bold"),
  legend.title    = element_text(size = 12, face = "bold"),
  legend.text     = element_text(size = 12),
  plot.margin     = margin(0.5, 0, 0.5, 0, "cm"),
  legend.position = "bottom"
 )
dev.off()
