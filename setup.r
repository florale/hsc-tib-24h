
base <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/monash/projects/healthy-sleep-clinic/data")
out <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/monash/projects/healthy-sleep-clinic/hsc-tib-24h/output")
redir <- paste0("/Users/", Sys.info()[["user"]], "/Library/CloudStorage/OneDrive-Personal/github/projects/hsc/")
source(paste0(redir, "setup.r"))

pal <- c("#708885", "#A9A9A9", "#4E2F26", "#ba6c6e", "#CA8F90", "#E4C7C7")
pal5 <- c("#E4C7C7", "#CA8F90", "#ba6c6e", "#4E2F26", "#708885")
col10 <- c(
  "#2A3E59",
  "#9c8aa4",
  "#ab8b8b",
  "#4F7375",

  "#8CAACB",
  "#ABA2C3",
  "#D1ACA5",
  "#769798"
)

col20 <- c(
  # "#1C1718",
  # "#5A6367",
  "#2A3E59",
  "#456691",
  "#647F9A",
  "#8CAACB",
  "#9c8aa4",
  "#ABA2C3",
  "#9A5C7D",
  "#B98AA3",
  "#cc8a8c",
  "#A54E50",
  "#DCD5CE",
  "#B49797",
  "#C99696",
  "#DAA5AE",
  "#d18d9a",
  "#b6485d",
  "#D1ACA5",
  "#C7AAA5",
  "#4F7375",
  "#769798",
  "#944C4C",
  "#ba6c6e",
  "#bf5b4b",
  "#bb847a",
  "#A69188",
  "#EAD3BF",
  "#FAD899",
  "#353D60",
  "#6171a9",
  "#8DA290",
  "#133A1B",
  "#6d765b",
  "#3b4031",
  "#c48462",
  "#3d251e",
  "#ab8b8b",
  "#D1ACA5"
)
col_brmcoda_d5 <-
  c(
    "#9A5C7D", "#B98AA3",
    "#DCD5CE", "#8DA290",
    "#708885", "#5A6367",
    "#456691", "#2A3E59",
    "#9c8aa4", "#5E4F65",
    "#1C1718"
  )
  col10 <-
  c(
    "#9A5C7D",
    "#5E4F65",
    "#456691",
    "#FAD899",
    "#708885",
    
    "#B98AA3",
    "#9c8aa4",
    "#8CAACB", 
    "#DCD5CE",
    "#8DA290",


    "#133A1B", "#6d765b",

    "#DAA5AE", "#b6485d",
    "#944C4C", "#C99696",
    "#bf5b4b", "#bb847a",
     
    "#FAD899", "#8DA290",
    "#133A1B", "#6d765b",
    "#3b4031", "#3d251e"
  )
