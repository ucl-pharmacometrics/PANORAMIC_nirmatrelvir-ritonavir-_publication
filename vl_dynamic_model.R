rm(list = ls())
# script to run viral dynamic models
library(ggplot2)
library(nlmixr2)
library(gridExtra)
#
# Read data
vl.dat <-  read.csv("paxlovid_nlmixr_data.csv",
                    head = TRUE, stringsAsFactors = FALSE)
#
# set dv to limit not LOQ/2 for cens = 1
vl.dat$dv[vl.dat$cens == 1 & vl.dat$dvid == "virus"] <- log10(112)

#### Add covariates for the final model: ####
# 1. Time varying treatment effect need on and off treatment as decline rate slows after treatment stops
vl.dat$pax_dich_on <- 0
vl.dat$pax_dich_on[vl.dat$treat_group2 == "Paxlovid" & vl.dat$t_start_therapy < 5.9] <- 1
# multiply by drug_eff to get as-treated analysis
vl.dat$mol_dich_on <- vl.dat$mol_dich_on * vl.dat$drug_eff
# 2. Allometric effect for age
vl.dat$log_age_cov <- log(vl.dat$age_y / median(vl.dat$age_y))
# 3. Allometric effect for antibody
vl.dat$log_ab_cov <- log(vl.dat$bl_antibody / median(vl.dat$bl_antibody,na.rm = T))
# 4. Allometric effect for time since symptom onset
vl.dat$t_symp_day <- vl.dat$t_symp_enrol + 1 
vl.dat$log_symp_cov <- log(vl.dat$t_symp_day  / median(vl.dat$t_symp_day))

#### run 11 full on baseline, only mol on slope - final model #####
si.mod.pax.dich.full  <- function() {
  ini({
    tdelta <-  -0.055   # Death rate of infected cells
    tv0    <- 6        # Viral load at treatment initiation
    eta.delta ~ 0.1
    eta.v0 ~ 0.1
    add.err <- 1.4     # residual variability
    
    beta_pax_on <- 0.3
    # beta_pax_off <- -0.42
    beta_age <- 0.9
    beta_ab <- -0.36 
    beta_sym <- -1.2
    
  })
  model({
    delta <- exp(tdelta + eta.delta + 
                   beta_pax_on * pax_dich_on)  # individual value of delta
    v0    <- 10^(tv0 + eta.v0 +
                   beta_sym * log_symp_cov+
                   beta_ab * log_ab_cov +
                   beta_age * log_age_cov)    # individual value of v0
    A_v(0) = v0
    
    d/dt(A_v) =  - delta * A_v 
    
    virus = log10(A_v) 
    virus ~ add(add.err)       # define error model
    
    treat_group = treat_group
    sample_group = sample_group
  })
}
# si.fit2.11 <- nlmixr2(si.mod.pax.dich.full, vl.dat[vl.dat$dvid == "virus", ],
#                       table = tableControl(cwres = TRUE, npde = TRUE),
#                       est = "foce")
# saveRDS(si.fit2.11, "si.fit2.11.rds")

# paramter estimation table
si.fit2.11 <- readRDS("si.fit2.11.rds")
pe <- si.fit2.11$parFixedDf
# take exp of beta pax on to get fractional change
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

vl.parest <- data.frame(matrix(ncol = 3, nrow = 7))
colnames(vl.parest) <- c("Parameter", "Estimate (95%CI)", "IIV (%CV)")
vl.parest$Parameter <- c("Delta (/d)", "V0 (log10 cp/mL)", "Additive error",
                         "Beta_pax_on", "Beta_age", "Beta_ab", "Beta_sym")
vl.parest$`Estimate (95%CI)` <- pe$est
vl.parest$`IIV (%CV)` <- pe$`BSV(CV% or SD)`
write.csv(vl.parest, "viral_dynamic_model_parameter_table.csv", row.names = FALSE)