
source("setup.r")

desc_by_so_min_p <- readRDS(file.path(out, "desc_by_so_min_p.rds"))
desc_all <- readRDS(file.path(out, "desc_all.rds"))

m_tib_isi8 <- readRDS(file.path(out, "m_tib_isi8.rds"))

pred_tib_isi8_sum <- readRDS(file.path(out, "pred_tib_isi8_sum.rds"))

pred_tib_isi8_age_sum <- readRDS(file.path(out, "pred_tib_isi8_age_sum.rds"))
pred_tib_isi8_sex_sum <- readRDS(file.path(out, "pred_tib_isi8_sex_sum.rds"))
pred_tib_isi8_serv_sum <- readRDS(file.path(out, "pred_tib_isi8_serv_sum.rds"))
pred_tib_isi11_sum <- readRDS(file.path(out, "pred_tib_isi11_sum.rds"))

plot_min_perc_isi8 <- pdftools::pdf_convert(file.path(out, "plot_min_perc_isi8.pdf"), filenames = file.path(out, "plot_min_perc_isi8.png"), dpi = 300)
plot_age_min_perc_isi8 <- pdftools::pdf_convert(file.path(out, "plot_age_min_perc_isi8.pdf"), filenames = file.path(out, "plot_age_min_perc_isi8.png"), dpi = 300)
plot_sex_min_perc_isi8 <- pdftools::pdf_convert(file.path(out, "plot_sex_min_perc_isi8.pdf"), filenames = file.path(out, "plot_sex_min_perc_isi8.png"), dpi = 300)
plot_serv_min_perc_isi8 <- pdftools::pdf_convert(file.path(out, "plot_serv_min_perc_isi8.pdf"), filenames = file.path(out, "plot_serv_min_perc_isi8.png"), dpi = 300)
plot_min_perc_isi11 <- pdftools::pdf_convert(file.path(out, "plot_min_perc_isi11.pdf"), filenames = file.path(out, "plot_min_perc_isi11.png"), dpi = 300)

## Descriptive statistics --------------------------------------------------
desc_all <- desc_all[, .(variable, level, summary, range)]
desc_by_so_min_p <- desc_by_so_min_p[, .(sleep_opportunity_percentile, variable, level, summary, range)]

# if level = 0, drop the row
desc_all <- desc_all[!(level == "0")]
desc_by_so_min_p <- desc_by_so_min_p[!(level == "0")]

desc_all[, level := fifelse(level == 1, "Yes", level)]
desc_by_so_min_p[, level := fifelse(level == 1, "Yes", level)]

desc_by_so_min_p[, sleep_opportunity_percentile := factor(
  sleep_opportunity_percentile,
  levels = c("≤P10", "P25–P75", "≥P90"),
  ordered = TRUE
)]
setorder(desc_by_so_min_p, variable, sleep_opportunity_percentile, level)

setnames(desc_all, "variable", "Variable")
setnames(desc_all, "level", "Category")
setnames(desc_all, "summary", "M (SD) / n (%)")
setnames(desc_all, "range", "Range")

setnames(desc_by_so_min_p, "sleep_opportunity_percentile", "Sleep Opportunity Percentile")
setnames(desc_by_so_min_p, "variable", "Variable")
setnames(desc_by_so_min_p, "level", "Category")
setnames(desc_by_so_min_p, "summary", "M (SD) / n (%)")
setnames(desc_by_so_min_p, "range", "Range")

## Estimated difference by sleep opportunity
### main -------------------------------------
pred_tib_isi8_sum[, `Sleep Opportunity Contrast` := fifelse(
 grepl("p10_p2575", par), "Short vs Medium",
 fifelse(
  grepl("p90_p2575", par), "Long vs Medium",
  fifelse(grepl("p10_p90", par), "Short vs Long", NA_character_)
 )
)]
pred_tib_isi8_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tib_isi8_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tib_isi8_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tib_isi8_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tib_isi8_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tib_isi8_sum <- pred_tib_isi8_sum[grepl("p10_p2575|p90_p2575|p10_p90", par)][, .(`Sleep Opportunity Contrast`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_sum <- sup_tib_isi8_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`, Sig)]
sup_tib_isi8_sum <- sup_tib_isi8_sum[!is.na(`Sleep-Wake Part`)]

### age -------------------------------------
pred_tib_isi8_age_sum[, `Sleep Opportunity Contrast` := fifelse(
 grepl("p10_p2575", par), "Short vs Medium",
 fifelse(
  grepl("p90_p2575", par), "Long vs Medium",
  fifelse(grepl("p10_p90", par), "Short vs Long", NA_character_)
 )
)]

pred_tib_isi8_age_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tib_isi8_age_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tib_isi8_age_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tib_isi8_age_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tib_isi8_age_sum[, Age := fifelse(age == "u45", "Below 45y", "45y and above")]
pred_tib_isi8_age_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]
pred_tib_isi8_age_sum[, `Age Contrast` := "45y and above vs Below 45y"]
pred_tib_isi8_age_sum[, `Sleep Opportunity Percentile` := tib_group]
pred_tib_isi8_age_sum[, `Sleep Opportunity Percentile` := ifelse(!is.na(tib_group), tib_group, "Overall")]

sup_tib_isi8_age_sum <- pred_tib_isi8_age_sum[grepl("p10_p2575|p90_p2575|p10_p90", par)][, .(Age, `Sleep-Wake Part`, `Sleep Opportunity Contrast`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_age_sum <- sup_tib_isi8_age_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`)]
sup_tib_isi8_age_sum <- sup_tib_isi8_age_sum[!is.na(`Sleep-Wake Part`)]

sup_tib_isi8_age_dif <- pred_tib_isi8_age_sum[grepl("_45u_u45", par)][, .(`Age Contrast`, `Sleep-Wake Part`, `Sleep Opportunity Percentile`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_age_dif <- sup_tib_isi8_age_dif[order( In, `Age Contrast`, `Sleep-Wake Part`)]
sup_tib_isi8_age_dif <- sup_tib_isi8_age_dif[!is.na(`Sleep-Wake Part`)]

### sex -------------------------------------
pred_tib_isi8_sex_sum[, `Sleep Opportunity Contrast` := fifelse(
  grepl("p10_p2575", par), "Short vs Medium",
  fifelse(
    grepl("p90_p2575", par), "Long vs Medium",
    fifelse(grepl("p10_p90", par), "Short vs Long", NA_character_)
  )
)]
pred_tib_isi8_sex_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tib_isi8_sex_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tib_isi8_sex_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tib_isi8_sex_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tib_isi8_sex_sum[, Sex := sex]
pred_tib_isi8_sex_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]
pred_tib_isi8_sex_sum[, `Sex Contrast` := "Female vs Male"]
pred_tib_isi8_sex_sum[, `Sleep Opportunity Percentile` := ifelse(!is.na(tib_group), tib_group, "Overall")]

sup_tib_isi8_sex_sum <- pred_tib_isi8_sex_sum[grepl("p10_p2575|p90_p2575|p10_p90", par)][, .(Sex, `Sleep-Wake Part`, `Sleep Opportunity Contrast`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_sex_sum <- sup_tib_isi8_sex_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`)]
sup_tib_isi8_sex_sum <- sup_tib_isi8_sex_sum[!is.na(`Sleep-Wake Part`)]

sup_tib_isi8_sex_dif <- pred_tib_isi8_sex_sum[grepl("female_male", par)][, .(`Sex Contrast`, `Sleep-Wake Part`, `Sleep Opportunity Percentile`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_sex_dif <- sup_tib_isi8_sex_dif[order( In, `Sex Contrast`, `Sleep-Wake Part`)]
sup_tib_isi8_sex_dif <- sup_tib_isi8_sex_dif[!is.na(`Sleep-Wake Part`)]
sup_tib_isi8_sex_dif[, ]
### service -------------------------------------
pred_tib_isi8_serv_sum[, `Sleep Opportunity Contrast` := fifelse(
  grepl("p10_p2575", par), "Short vs Medium",
  fifelse(
    grepl("p90_p2575", par), "Long vs Medium",
    fifelse(grepl("p10_p90", par), "Short vs Long", NA_character_)
  )
)]
pred_tib_isi8_serv_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tib_isi8_serv_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tib_isi8_serv_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tib_isi8_serv_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tib_isi8_serv_sum[, Service := fifelse(service == "Home", "Home PSG", "Lab PSG")]
pred_tib_isi8_serv_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]
pred_tib_isi8_serv_sum[, `Service Contrast` := "Lab PSG vs Home PSG"]
pred_tib_isi8_serv_sum[, `Sleep Opportunity Percentile` := ifelse(!is.na(tib_group), tib_group, "Overall")]

sup_tib_isi8_serv_sum <- pred_tib_isi8_serv_sum[grepl("p10_p2575|p90_p2575|p10_p90", par)][, .(Service, `Sleep-Wake Part`, `Sleep Opportunity Contrast`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_serv_sum <- sup_tib_isi8_serv_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`)]
sup_tib_isi8_serv_sum <- sup_tib_isi8_serv_sum[!is.na(`Sleep-Wake Part`)]

sup_tib_isi8_serv_dif <- pred_tib_isi8_serv_sum[grepl("lab_home", par)][, .(`Service Contrast`, `Sleep-Wake Part`, `Sleep Opportunity Percentile`, `Estimated Difference`, In, Sig)]
sup_tib_isi8_serv_dif <- sup_tib_isi8_serv_dif[order( In, `Service Contrast`, `Sleep-Wake Part`)]
sup_tib_isi8_serv_dif <- sup_tib_isi8_serv_dif[!is.na(`Sleep-Wake Part`)]

### cut off of 11 -------------------------------------
pred_tib_isi11_sum[, `Sleep Opportunity Contrast` := fifelse(
 grepl("p10_p2575", par), "Short vs Medium",
 fifelse(
  grepl("p90_p2575", par), "Long vs Medium",
  fifelse(grepl("p10_p90", par), "Short vs Long", NA_character_)
 )
)]
pred_tib_isi11_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tib_isi11_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tib_isi11_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tib_isi11_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tib_isi11_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tib_isi11_sum <- pred_tib_isi11_sum[grepl("p10_p2575|p90_p2575|p10_p90", par)][, .(`Sleep Opportunity Contrast`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tib_isi11_sum <- sup_tib_isi11_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`, Sig)]
sup_tib_isi11_sum <- sup_tib_isi11_sum[!is.na(`Sleep-Wake Part`)]