
source("setup.r")
library(pdftools)

m_tib_isi1 <- readRDS(file.path(out, "m_tib_isi1.rds"))

pred_tibq5_isi1_sum <- readRDS(file.path(out, "pred_tibq5_isi1_sum.rds"))

pred_tibq5_isi1_age_u45_sum <- readRDS(file.path(out, "pred_tibq5_isi1_age_u45_sum.rds"))
pred_tibq5_isi1_age_45u_sum <- readRDS(file.path(out, "pred_tibq5_isi1_age_45u_sum.rds"))

pred_tibq5_isi1_female_sum <- readRDS(file.path(out, "pred_tibq5_isi1_female_sum.rds"))
pred_tibq5_isi1_male_sum <- readRDS(file.path(out, "pred_tibq5_isi1_male_sum.rds"))

pred_tibq5_isi1_serv_home_sum <- readRDS(file.path(out, "pred_tibq5_isi1_serv_home_sum.rds"))
pred_tibq5_isi1_serv_lab_sum <- readRDS(file.path(out, "pred_tibq5_isi1_serv_lab_sum.rds"))

plot_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_min_perc_isi1.pdf"), filenames = file.path(out, "plot_min_perc_isi1.png"), dpi = 300)
plot_age_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_age_min_perc_isi1.pdf"), filenames = file.path(out, "plot_age_min_perc_isi1.png"), dpi = 300)
plot_sex_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_sex_min_perc_isi1.pdf"), filenames = file.path(out, "plot_sex_min_perc_isi1.png"), dpi = 300)
plot_serv_min_perc_isi1 <- pdftools::pdf_convert(file.path(out, "plot_serv_min_perc_isi1.pdf"), filenames = file.path(out, "plot_serv_min_perc_isi1.png"), dpi = 300)

## Estimated difference by sleep opportunity
### adj -------------------------------------
pred_tibq5_isi1_sum[, `Sleep Opportunity` := fifelse(
 grepl("q1_q234", par), "Short - Medium",
 fifelse(
  grepl("q5_q234", par), "Long - Medium",
  fifelse(grepl("q1_q5", par), "Short - Long", NA_character_)
 )
)]
pred_tibq5_isi1_sum[, `Sleep Opportunity` := factor(`Sleep Opportunity`, levels = c("Short - Medium", "Long - Medium", "Short - Long"))]
pred_tibq5_isi1_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tibq5_isi1_sum <- pred_tibq5_isi1_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(`Sleep Opportunity`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_sum <- sup_tibq5_isi1_sum[order( In, `Sleep Opportunity`, `Sleep-Wake Part`, Sig)]
sup_tibq5_isi1_sum <- sup_tibq5_isi1_sum[!is.na(`Sleep-Wake Part`)]

### age -------------------------------------
pred_tibq5_isi1_age_u45_sum[, age := "u45"]
pred_tibq5_isi1_age_45u_sum[, age := "45u"]
colnames(pred_tibq5_isi1_age_u45_sum)
colnames(pred_tibq5_isi1_age_45u_sum)

pred_tibq5_isi1_age_sum <- rbind(
  pred_tibq5_isi1_age_u45_sum,
  pred_tibq5_isi1_age_45u_sum
)
pred_tibq5_isi1_age_sum[, `Sleep Opportunity` := fifelse(
 grepl("q1_q234", par), "Short - Medium",
 fifelse(
  grepl("q5_q234", par), "Long - Medium",
  fifelse(grepl("q1_q5", par), "Short - Long", NA_character_)
 )
)]
pred_tibq5_isi1_age_sum[, `Sleep Opportunity` := factor(`Sleep Opportunity`, levels = c("Short - Medium", "Long - Medium", "Short - Long"))]
pred_tibq5_isi1_age_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_age_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_age_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_age_sum[, Age := fifelse(age == "u45", "Below 45y", "45y and above")]
pred_tibq5_isi1_age_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tibq5_isi1_age_sum <- pred_tibq5_isi1_age_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(Age, `Sleep Opportunity`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_age_sum <- sup_tibq5_isi1_age_sum[order( In, `Sleep Opportunity`, `Sleep-Wake Part`, Age, Sig)]
sup_tibq5_isi1_age_sum <- sup_tibq5_isi1_age_sum[!is.na(`Sleep-Wake Part`)]

### sex -------------------------------------
pred_tibq5_isi1_female_sum[, sex := "Female"]
pred_tibq5_isi1_male_sum[, sex := "Male"]

pred_tibq5_isi1_sex_sum <- rbind(
  pred_tibq5_isi1_female_sum,
  pred_tibq5_isi1_male_sum
)
pred_tibq5_isi1_sex_sum[, `Sleep Opportunity` := fifelse(
  grepl("q1_q234", par), "Short - Medium",
  fifelse(
    grepl("q5_q234", par), "Long - Medium",
    fifelse(grepl("q1_q5", par), "Short - Long", NA_character_)
  )
)]
pred_tibq5_isi1_sex_sum[, `Sleep Opportunity` := factor(`Sleep Opportunity`, levels = c("Short - Medium", "Long - Medium", "Short - Long"))]
pred_tibq5_isi1_sex_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_sex_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_sex_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_sex_sum[, Sex := sex]
pred_tibq5_isi1_sex_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tibq5_isi1_sex_sum <- pred_tibq5_isi1_sex_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(Sex, `Sleep Opportunity`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_sex_sum <- sup_tibq5_isi1_sex_sum[order(In, `Sleep Opportunity`, `Sleep-Wake Part`, Sex, Sig)]
sup_tibq5_isi1_sex_sum <- sup_tibq5_isi1_sex_sum[!is.na(`Sleep-Wake Part`)]

### service -------------------------------------
pred_tibq5_isi1_serv_home_sum[, service := "Home"]
pred_tibq5_isi1_serv_lab_sum[, service := "Lab"]

pred_tibq5_isi1_serv_sum <- rbind(
  pred_tibq5_isi1_serv_home_sum,
  pred_tibq5_isi1_serv_lab_sum
)
pred_tibq5_isi1_serv_sum[, `Sleep Opportunity` := fifelse(
  grepl("q1_q234", par), "Short - Medium",
  fifelse(
    grepl("q5_q234", par), "Long - Medium",
    fifelse(grepl("q1_q5", par), "Short - Long", NA_character_)
  )
)]
pred_tibq5_isi1_serv_sum[, `Sleep Opportunity` := factor(`Sleep Opportunity`, levels = c("Short - Medium", "Long - Medium", "Short - Long"))]
pred_tibq5_isi1_serv_sum[, `Sleep-Wake Part` := factor(part_label, levels = c("WD", "SOL", "WASO", "N1+2", "N3", "REM"))]
pred_tibq5_isi1_serv_sum[, `Estimated Difference` := fifelse(!is.na(est_min), est_min, est_perc)]
pred_tibq5_isi1_serv_sum[, In := fifelse(!is.na(est_min), "Minutes", "Percentage")]
pred_tibq5_isi1_serv_sum[, Service := fifelse(service == "Home", "Home PSG", "Lab PSG")]
pred_tibq5_isi1_serv_sum[, Sig := ifelse(!between(0, CI_low, CI_high), "$\\ast$", "")]

sup_tibq5_isi1_serv_sum <- pred_tibq5_isi1_serv_sum[grepl("q1_q234|q5_q234|q1_q5", par)][, .(Service, `Sleep Opportunity`, `Sleep-Wake Part`, `Estimated Difference`, In, Sig)]
sup_tibq5_isi1_serv_sum <- sup_tibq5_isi1_serv_sum[order(In, `Sleep Opportunity`, Service, `Sleep-Wake Part`, Sig)]
sup_tibq5_isi1_serv_sum <- sup_tibq5_isi1_serv_sum[!is.na(`Sleep-Wake Part`)]
