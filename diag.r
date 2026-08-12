source("setup.r")

library(brms)
library(data.table)

m_tib_isi8 <- readRDS(file.path(out, "m_tib_isi8.rds"))
summary(m_tib_isi8)

# Summarise R-hat values to diagnose convergence
rhats <- brms::rhat(m_tib_isi8)
rhat_summary <- data.table(parameter = names(rhats), rhat = as.numeric(rhats))
setorder(rhat_summary, -rhat)

print(rhat_summary[1:20])
cat(sprintf("Max R-hat: %.3f\n", max(rhat_summary$rhat)))
cat(sprintf("Parameters with R-hat > 1.01: %d\n", rhat_summary[rhat > 1.01, .N]))

pp_check(m_tib_isi8, resp = "z11", ndraws = 100) +
  ggtitle("Posterior Predictive Check - ILR coordinate 1/5") +
  theme_minimal()
