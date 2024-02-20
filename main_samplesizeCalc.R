rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')
source("main_infections.R")

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
        names.arg=c("Baseline", "Low/\nlow", "Low/\nhigh",
                    "High/\nbaseline", "Low/\nbaseline"),
        xlab = "Scenario", ylab = "Participants per cluster")

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
        ylim=c(-0.5, 1.1), 
        beside = TRUE, legend=rownames(cl_df),
        args.legend = list(x = "topright", bty="n"),
        names.arg=c("low/\nlow", "low/\nhigh",
                    "high/\nbaseline", "low/\nbaseline"),
        xlab = "scenario", ylab = "% difference relative to baseline")

# calculate relative power for each scenario given baseline sample size --------

# baseline -------
power.inf_d <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_d, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_d, clusters_perarm = 15)
power.d_d <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_d, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.cases_d, clusters_perarm = 15)

# L/L --------
power.inf_ll <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_ll, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_d, clusters_perarm = 15)
power.d_ll <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_ll, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.cases_d, clusters_perarm = 15)

# L/H --------
power.inf_lh <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_lh, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_d, clusters_perarm = 15)
power.d_lh <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_lh, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.cases_d, clusters_perarm = 15)

# H/D --------
power.inf_hd <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_hd, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_d, clusters_perarm = 15)
power.d_hd <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_hd, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.cases_d, clusters_perarm = 15)

# L/D --------
power.inf_ld <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_ld, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_d, clusters_perarm = 15)
power.d_ld <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_ld, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.cases_d, clusters_perarm = 15)

cl_df$baseline_power <- c(power.inf_d, power.d_d)
cl_df$ll_power <- c(power.inf_ll, power.d_ll)
cl_df$lh_power <- c(power.inf_lh, power.d_lh)
cl_df$hd_power <- c(power.inf_hd, power.d_hd)
cl_df$ld_power <- c(power.inf_ld, power.d_ld)

barplot(as.matrix(cl_df[c("baseline_power", "ll_power", "lh_power", "hd_power", "ld_power")]), 
        beside = TRUE, legend=rownames(cl_df),
        args.legend = list(x = 16.5, y = 1, bty="n"), ylim=c(0,1),
        names.arg=c("Baseline", "Low/\nlow", "Low/\nhigh",
                    "High/\nbaseline", "Low/\nbaseline"),
        xlab = "Scenario", ylab = "Power")

cl_df$ll_power_perc <- (cl_df$ll_power-cl_df$baseline_power)/cl_df$baseline_power
cl_df$lh_power_perc <- (cl_df$lh_power-cl_df$baseline_power)/cl_df$baseline_power
cl_df$hd_power_perc <- (cl_df$hd_power-cl_df$baseline_power)/cl_df$baseline_power
cl_df$ld_power_perc <- (cl_df$ld_power-cl_df$baseline_power)/cl_df$baseline_power

barplot(as.matrix(cl_df[c("ll_power_perc", "lh_power_perc", "hd_power_perc", "ld_power_perc")]),
        ylim=c(-0.5,0.5), 
        beside = TRUE, legend=rownames(cl_df),
        args.legend = list(x = "topright", bty="n"),
        names.arg=c("Low/\nlow", "Low/\nhigh",
                    "High/\nbaseline", "Low/\nbaseline"),
        xlab = "Scenario", ylab = "Power relative to baseline",
        axes=F)
labs <- seq(-0.5,1,0.25)
axis(side=2, at=labs, labels=paste0(labs*100, "%"))

{
  ##need to fix legend plotting
  layout(matrix(1:2, ncol=2, byrow=T))
  
  par(omi = c(0, 0, 0, 0.6),
      xpd=T)
  
  barplot(as.matrix(cl_df[c("ll_perc", "lh_perc", "hd_perc", "ld_perc")]),
          ylim=c(-0.5, 1.1), 
          beside = TRUE,
          names.arg=c("Low/\nlow", "Low/\nhigh",
                      "High/\nbaseline", "Low/\nbaseline"),
          xlab = "Scenario", ylab = "Participants per cluster relative to baseline",
          axes=F)
  labs <- seq(-0.5,1,0.25)
  axis(side=2, at=labs, labels=paste0(labs*100, "%"))
  
  barplot(as.matrix(cl_df[c("ll_power_perc", "lh_power_perc", "hd_power_perc", "ld_power_perc")]),
          ylim=c(-0.5,0.5),
          beside = TRUE,
          names.arg=c("Low/\nlow", "Low/\nhigh",
                      "High/\nbaseline", "Low/\nbaseline"),
          xlab = "Scenario", ylab = "Power relative to baseline",
          axes=F, xpd=F)
  labs <- seq(-0.5,1,0.25)
  axis(side=2, at=labs, labels=paste0(labs*100, "%"), xpd=F)
  
  legend(par('usr')[2], par('usr')[4],
         c("Active", "Passive"), title="Surveillance type",
         col=grey.colors(2), pch=15,
         bty="n", xpd=NA)
  
}

























