
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

col_serv <- c(
  "#4F5D5B", 
  "#708885", 
  # "#3d251e", 
  "#978787",
  "#EAD3BF", 
  "#c48462",

  "#4F5D5B", 
  "#708885", 
  # "#3d251e", 
  "#978787",
  "#EAD3BF", 
  "#c48462"
)

shape_perc <- c(21)
shape_sex <- c("Female" = 24, "Male" = 25)
shape_age <- c("u45" = "\u25c1", "45u" = "\u25b7")
shape_serv <- c("home" = 22, "lab" = 23)

shape_min_sex <- c("Female" = "\u25b3", "Male" = "\u25bd")
shape_min_age <- c("u45" = "\u25C0", "45u" = "\u25B6")
shape_min_serv <- c("home" = 22, "lab" = 23)

shape_perc_sex <- c("Female" = 24, "Male" = 25)
  
col_age <- c(
  "#4F5D5B", 
  "#708885", 
  "#978787", 
  "#EAD3BF", 
  "#c48462",

  "#4F5D5B", 
  "#708885", 
  "#978787", 
  "#EAD3BF", 
  "#c48462"
)
col_sex <- col_serv <- col_age

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
