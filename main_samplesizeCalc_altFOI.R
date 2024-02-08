rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')
source("main_infections_altFOI.R")

# Parameters --------------------------------------------------------------

font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation for default ----------------------------------------------------

cl.perarm.inf_d <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_d, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.cases_d <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_d, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# Do power calculation for L/L ----------------------------------------------------

cl.perarm.inf_ll <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_ll, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.cases_ll <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_ll, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# Do power calculation for L/H ----------------------------------------------------

cl.perarm.inf_lh <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_lh, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.cases_lh <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_lh, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# Do power calculation for H/D ----------------------------------------------------

cl.perarm.inf_hd <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_hd, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.cases_hd <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_hd, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# Do power calculation for L/D ----------------------------------------------------

cl.perarm.inf_ld <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_ld, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.cases_ld <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_ld, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

##plot as bar chart

cl_df <- data.frame(row.names=c("active", "passive"),
                    d=c(cl.perarm.inf_d, cl.perarm.cases_d),
                    ll=c(cl.perarm.inf_ll, cl.perarm.cases_ll),
                    lh=c(cl.perarm.inf_lh, cl.perarm.cases_lh),
                    hd=c(cl.perarm.inf_hd, cl.perarm.cases_hd),
                    ld=c(cl.perarm.inf_ld, cl.perarm.cases_ld))
#cl_df <- do.call(rbind, cl_df)
barplot(as.matrix(cl_df), beside = TRUE, legend=rownames(cl_df),
        args.legend = list(x = "topright", bty="n"), ylim=c(0,2000),
        names.arg=c("baseline", "low/\nlow", "low/\nhigh",
                    "high/\nbaseline", "low/\nbaseline"),
        xlab = "scenario", ylab = "no. of pt per cluster")

##plot relative difference in sample size by FOI scenario
cl_df$ll_fold <- cl_df$ll/cl_df$d
cl_df$lh_fold <- cl_df$lh/cl_df$d
cl_df$hd_fold <- cl_df$hd/cl_df$d
cl_df$ld_fold <- cl_df$ld/cl_df$d

cl_df$ll_perc <- (cl_df$ll - cl_df$d)/cl_df$d
cl_df$lh_perc <- (cl_df$lh - cl_df$d)/cl_df$d
cl_df$hd_perc <- (cl_df$hd - cl_df$d)/cl_df$d
cl_df$ld_perc <- (cl_df$ld - cl_df$d)/cl_df$d

barplot(as.matrix(cl_df[c("ll_perc", "lh_perc", "hd_perc", "ld_perc")]),
        ylim=c(-0.5, 1), 
        beside = TRUE, legend=rownames(cl_df),
        args.legend = list(x = "topright", bty="n"),
        names.arg=c("low/\nlow", "low/\nhigh",
                    "high/\nbaseline", "low/\nbaseline"),
        xlab = "scenario", ylab = "% change relative to baseline")
abline(h=1,lty=2)





