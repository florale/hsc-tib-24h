
source("setup.r")
library(pdftools)

m_tib_isi1 <- readRDS(file.path(out, "m_tib_isi1.rds"))

pred_tibq5_isi1_sum <- readRDS(file.path(out, "pred_tibq5_isi1_sum.rds"))

pred_tibq5_isi1_age_sum <- readRDS(file.path(out, "pred_tibq5_isi1_age_sum.rds"))
pred_tibq5_isi1_sex_sum <- readRDS(file.path(out, "pred_tibq5_isi1_sex_sum.rds"))
pred_tibq5_isi1_serv_sum <- readRDS(file.path(out, "pred_tibq5_isi1_serv_sum.rds"))
pred_tibq5_isi1_11c_sum <- readRDS(file.path(out, "pred_tibq5_isi1_11c_sum.rds"))

plot_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_min_perc_isi1.pdf"), filenames = file.path(out, "plot_min_perc_isi1.png"), dpi = 300)
plot_age_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_age_min_perc_isi1.pdf"), filenames = file.path(out, "plot_age_min_perc_isi1.png"), dpi = 300)
plot_sex_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_sex_min_perc_isi1.pdf"), filenames = file.path(out, "plot_sex_min_perc_isi1.png"), dpi = 300)
plot_serv_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_serv_min_perc_isi1.pdf"), filenames = file.path(out, "plot_serv_min_perc_isi1.png"), dpi = 300)
plot_min_perc_isi1_11c <- pdftools::pdf_convert(file.path(out, "plot_min_perc_isi1_11c.pdf"), filenames = file.path(out, "plot_min_perc_isi1_11c.png"), dpi = 300)

## Estimated difference by sleep opportunity
### main -------------------------------------
pred_tibq5_isi1_sum[, `Sleep Opportunity Contrast` := fifelse(
 grepl("q1_q234", par), "Short vs Medium",
 fifelse(
  grepl("q5_q234", par), "Long vs Medium",
  fifelse(grepl("q1_q5", par), "Short vs Long", NA_character_)
 )
)]
pred_tibq5_isi1_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tibq5_isi1_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tibq5_isi1_sum <- pred_tibq5_isi1_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(`Sleep Opportunity Contrast`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_sum <- sup_tibq5_isi1_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`, Sig)]
sup_tibq5_isi1_sum <- sup_tibq5_isi1_sum[!is.na(`Sleep-Wake Part`)]

### age -------------------------------------
pred_tibq5_isi1_age_sum[, `Sleep Opportunity Contrast` := fifelse(
 grepl("q1_q234", par), "Short vs Medium",
 fifelse(
  grepl("q5_q234", par), "Long vs Medium",
  fifelse(grepl("q1_q5", par), "Short vs Long", NA_character_)
 )
)]

pred_tibq5_isi1_age_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tibq5_isi1_age_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_age_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_age_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_age_sum[, Age := fifelse(age == "u45", "Below 45y", "45y and above")]
pred_tibq5_isi1_age_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]
pred_tibq5_isi1_age_sum[, `Age Contrast` := "45y and above vs Below 45y"]
pred_tibq5_isi1_age_sum[, `Sleep Opportunity Percentile` := tib_group]
pred_tibq5_isi1_age_sum[, `Sleep Opportunity Percentile` := ifelse(!is.na(tib_group), tib_group, "Overall")]

sup_tibq5_isi1_age_sum <- pred_tibq5_isi1_age_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(Age, `Sleep-Wake Part`, `Sleep Opportunity Contrast`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_age_sum <- sup_tibq5_isi1_age_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`)]
sup_tibq5_isi1_age_sum <- sup_tibq5_isi1_age_sum[!is.na(`Sleep-Wake Part`)]

sup_tibq5_isi1_age_dif <- pred_tibq5_isi1_age_sum[grepl("_45u_u45", par)][, .(`Age Contrast`, `Sleep-Wake Part`, `Sleep Opportunity Percentile`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_age_dif <- sup_tibq5_isi1_age_dif[order( In, `Age Contrast`, `Sleep-Wake Part`)]
sup_tibq5_isi1_age_dif <- sup_tibq5_isi1_age_dif[!is.na(`Sleep-Wake Part`)]

### sex -------------------------------------
pred_tibq5_isi1_sex_sum[, `Sleep Opportunity Contrast` := fifelse(
  grepl("q1_q234", par), "Short vs Medium",
  fifelse(
    grepl("q5_q234", par), "Long vs Medium",
    fifelse(grepl("q1_q5", par), "Short vs Long", NA_character_)
  )
)]
pred_tibq5_isi1_sex_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tibq5_isi1_sex_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_sex_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_sex_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_sex_sum[, Sex := sex]
pred_tibq5_isi1_sex_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]
pred_tibq5_isi1_sex_sum[, `Sex Contrast` := "Female vs Male"]
pred_tibq5_isi1_sex_sum[, `Sleep Opportunity Percentile` := ifelse(!is.na(tib_group), tib_group, "Overall")]

sup_tibq5_isi1_sex_sum <- pred_tibq5_isi1_sex_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(Sex, `Sleep-Wake Part`, `Sleep Opportunity Contrast`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_sex_sum <- sup_tibq5_isi1_sex_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`)]
sup_tibq5_isi1_sex_sum <- sup_tibq5_isi1_sex_sum[!is.na(`Sleep-Wake Part`)]

sup_tibq5_isi1_sex_dif <- pred_tibq5_isi1_sex_sum[grepl("female_male", par)][, .(`Sex Contrast`, `Sleep-Wake Part`, `Sleep Opportunity Percentile`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_sex_dif <- sup_tibq5_isi1_sex_dif[order( In, `Sex Contrast`, `Sleep-Wake Part`)]
sup_tibq5_isi1_sex_dif <- sup_tibq5_isi1_sex_dif[!is.na(`Sleep-Wake Part`)]
sup_tibq5_isi1_sex_dif[, ]
### service -------------------------------------
pred_tibq5_isi1_serv_sum[, `Sleep Opportunity Contrast` := fifelse(
  grepl("q1_q234", par), "Short vs Medium",
  fifelse(
    grepl("q5_q234", par), "Long vs Medium",
    fifelse(grepl("q1_q5", par), "Short vs Long", NA_character_)
  )
)]
pred_tibq5_isi1_serv_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tibq5_isi1_serv_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_serv_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_serv_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_serv_sum[, Service := fifelse(service == "Home", "Home PSG", "Lab PSG")]
pred_tibq5_isi1_serv_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]
pred_tibq5_isi1_serv_sum[, `Service Contrast` := "Lab PSG vs Home PSG"]
pred_tibq5_isi1_serv_sum[, `Sleep Opportunity Percentile` := ifelse(!is.na(tib_group), tib_group, "Overall")]

sup_tibq5_isi1_serv_sum <- pred_tibq5_isi1_serv_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(Service, `Sleep-Wake Part`, `Sleep Opportunity Contrast`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_serv_sum <- sup_tibq5_isi1_serv_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`)]
sup_tibq5_isi1_serv_sum <- sup_tibq5_isi1_serv_sum[!is.na(`Sleep-Wake Part`)]

sup_tibq5_isi1_serv_dif <- pred_tibq5_isi1_serv_sum[grepl("lab_home", par)][, .(`Service Contrast`, `Sleep-Wake Part`, `Sleep Opportunity Percentile`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_serv_dif <- sup_tibq5_isi1_serv_dif[order( In, `Service Contrast`, `Sleep-Wake Part`)]
sup_tibq5_isi1_serv_dif <- sup_tibq5_isi1_serv_dif[!is.na(`Sleep-Wake Part`)]

### cut off of 11 -------------------------------------
pred_tibq5_isi1_11c_sum[, `Sleep Opportunity Contrast` := fifelse(
 grepl("q1_q234", par), "Short vs Medium",
 fifelse(
  grepl("q5_q234", par), "Long vs Medium",
  fifelse(grepl("q1_q5", par), "Short vs Long", NA_character_)
 )
)]
pred_tibq5_isi1_11c_sum[, `Sleep Opportunity Contrast` := factor(`Sleep Opportunity Contrast`, levels = c("Short vs Medium", "Long vs Medium", "Short vs Long"))]
pred_tibq5_isi1_11c_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_11c_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_11c_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_11c_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tibq5_isi1_11c_sum <- pred_tibq5_isi1_11c_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(`Sleep Opportunity Contrast`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_11c_sum <- sup_tibq5_isi1_11c_sum[order( In, `Sleep Opportunity Contrast`, `Sleep-Wake Part`, Sig)]
sup_tibq5_isi1_11c_sum <- sup_tibq5_isi1_11c_sum[!is.na(`Sleep-Wake Part`)]