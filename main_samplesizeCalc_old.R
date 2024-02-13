rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')
source("main_infections.R")

# Parameters --------------------------------------------------------------

font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation for L/L ----------------------------------------------------

cl.perarm.inf_ll <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_ll <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

# Do power calculation for L/H ----------------------------------------------------

cl.perarm.inf_lh <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_lh <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

# Do power calculation for H/D ----------------------------------------------------

cl.perarm.inf_hd <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_hd <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

# Do power calculation for L/D ----------------------------------------------------

cl.perarm.inf_ld <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_ld <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)


## calc power with baseline sample size for each scenario -----

# Do power calculation for L/L ----------------------------------------------------

cl.perarm.inf_ll <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_ll <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

# Do power calculation for L/H ----------------------------------------------------

cl.perarm.inf_lh <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_lh <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

# Do power calculation for H/D ----------------------------------------------------

cl.perarm.inf_hd <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_hd <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

# Do power calculation for L/D ----------------------------------------------------

cl.perarm.inf_ld <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases_ld <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)











