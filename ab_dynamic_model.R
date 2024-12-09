rm(list = ls())
# script to run viral dynamic models
library(ggplot2)
library(nlmixr2)
library(gridExtra)
#
# Read data
ab.dat <-  read.csv("paxlovid_nlmixr_data.csv",
                    head = TRUE, stringsAsFactors = FALSE)
ab.dat <- ab.dat[ab.dat$dvid == "s_antibody",]

#### Add covariates for the final model: ####
# 1. Vaccine effect as number of vaccines
ab.dat[ab.dat$n_vacccine == -99,]$n_vacccine <- median(ab.dat$n_vacccine)
ab.dat$log_vac <- log(ab.dat$n_vacccine/median(ab.dat$n_vacccine))
# 2. Dichotomous drug effect, drug_eff gives as-treated analysis
# ab.dat$drug_eff
si.mod.full.ab <- function() {
  ini({
    tdelta.ab <- -2.11  # Increase in s antibody
    tab0  <- 3.112      # Viral load at treatment initiation
    eta.delta.ab ~ 0.15
    eta.ab.ab0  ~ 0.13
    
    add.err.ab <- 0.22     # residual variability
    
    beta_pax <- -0.2
    beta_1_vac <- -0.11
    beta_2_vac <- 0.3
  })
  model({
    delta.ab <- exp(tdelta.ab + eta.delta.ab +
                      beta_pax * drug_eff+
                      beta_1_vac * log_vac)  # individual value of delta
    a0    <- 10^(tab0 + eta.ab.ab0+
                   beta_2_vac * log_vac)        # individual value of v0
    A_a(0) = a0
    
    d/dt(A_a) =  delta.ab * A_a 
    
    s_antibody = log10(A_a) 
    s_antibody ~ add(add.err.ab)       # define error model
    
    treat_group = treat_group
    sample_group = sample_group
  })
}
# si.fit.11.ab <- nlmixr2(si.mod.full.ab , ab.dat,
#                         table = tableControl(cwres = TRUE, npde = TRUE),
#                         est = "foce")
# saveRDS(si.fit.11.ab, "si.fit.11.ab.rds")

# final model saved as rds
si.fit.11.ab <- readRDS("si.fit.11.ab.rds")
###### make parameter estimate table #######
pe <- si.fit.11.ab$parFixedDf
# take exp of all covariates to get fractional change
pe$`Back-transformed`[4] <- exp(pe$`Back-transformed`[4])
pe$`CI Lower`[4] <- exp(pe$`CI Lower`[4])
pe$`CI Upper`[4] <- exp(pe$`CI Upper`[4])
pe$`Back-transformed` <- round(pe$`Back-transformed`, 3)
pe$`CI Lower` <- round(pe$`CI Lower`, 3)
pe$`CI Upper` <- round(pe$`CI Upper`, 3)
pe$est <- paste0(pe$`Back-transformed`, " (", pe$`CI Lower`, ", ", pe$`CI Upper`, ")")
pe$est[3] <- pe$`Back-transformed`[3]
pe$`BSV(CV% or SD)`[1] <- signif(pe$`BSV(CV% or SD)`[1], 3)
pe$`BSV(CV% or SD)`[2] <- signif(pe$`BSV(CV% or SD)`[2]/pe$Estimate[2] * 100, 3)
pe$`BSV(CV% or SD)`[is.na(pe$`BSV(CV% or SD)`)] <- "-"

ab.parest <- data.frame(matrix(ncol = 3, nrow = 6))
colnames(ab.parest) <- c("Parameter", "Estimate (95%CI)", "IIV (%CV)")
ab.parest$Parameter <- c("Delta (/d)", "A0 (log10 U/mL)", "Additive error","Beta_pax",
                         "Beta_1_vac", "Beta_2_vac")
ab.parest$`Estimate (95%CI)` <- pe$est
ab.parest$`IIV (%CV)` <- pe$`BSV(CV% or SD)`
write.csv(ab.parest, "antibody_dynamic_model_parameter_table.csv", row.names = FALSE)
