
base <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/monash/projects/healthy-sleep-clinic/data")
out <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/monash/projects/healthy-sleep-clinic/hsc-tib-24h/output")
redir <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/github/projects/hsc/")
source(paste0(redir, "setup.r"))

scales::show_col(tvthemes:::hilda_palette$Day)
scales::show_col(tvthemes:::hilda_palette$Night)
scales::show_col(tvthemes:::hilda_palette$Dusk)
library(extrafont)

col <- c(
  "#708885", 
  "#3d251e", 
  "#978787",
  "#EAD3BF", 
  "#c48462"
)

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

