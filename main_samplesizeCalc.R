rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')
source('main_infections.R')

# Parameters --------------------------------------------------------------

font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation ----------------------------------------------------
# Pi_0 equals inf_c (proportion of population experiencing the outcome) or case_c_2yr

cl.perarm.inf <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.094, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

effectsizes <- c(0.25,0.3,0.35)
k_vec <- c(0.02,0.15,0.25)

# plots for presentations -------------------------------------------------
# infection
par(mar = c(5.1, 5.1, 4.1, 2.1))

cl.perarm.inf_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.094, 
                                  treatment_effect = effectsizes[2], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.case_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[2], 
                                  nr.percluster = 20:2000)

plot(x=20:200, y=cl.perarm.inf_k2_e2, type="l", col="black", ylim=c(0,20), lwd=lwd,
     xlab="Number of children per cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.3", cex.axis=font.size, cex.lab=font.size, cex.main=font.size)
abline(h=15, lwd=lwd, lty=2)
abline(v=c(20:200)[124], lwd=lwd, lty=2)

# disease
plot(x=20:2000, y=cl.perarm.case_k2_e2, type="l", col="black", lwd=lwd,
     xlab="Number of people per cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.3", cex.axis=font.size, cex.lab=font.size, cex.main=font.size,
     log="y", ylim=c(10,500), xlim=c(0,1000))
abline(h=15, lwd=lwd, lty=2)
abline(v=c(20:2000)[948], lwd=lwd, lty=2)


# plots for report --------------------------------------------------------

plot_sample_infections()
plot_sample_cases()

# text for report:
#nr of participants per cluster 

  # comparing impact of different effect sizes:
## infection endpoint, k=0.15, effect size=0.3
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
           treatment_effect = effectsizes[2], k = k_vec[2], 
           nr.percluster = 20:200) # 149 per arm

## infection endpoint, k=0.15, effect size=0.35
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
           treatment_effect = effectsizes[3], k = k_vec[2], 
           nr.percluster = 20:200) #98 per arm (34% less)

## infection endpoint, k=0.15, effect size=0.25
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
           treatment_effect = effectsizes[1], k = k_vec[2], 
           nr.percluster = 20:300) # 259 per arm


## case endpoint, k=0.15, effect size=0.3
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[2], k = k_vec[2], 
           nr.percluster = 20:2000) # 967 per arm

## case endpoint, k=0.15, effect size=0.35
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[3], k = k_vec[2], 
           nr.percluster = 20:2000) # 635 per arm

## case endpoint, k=0.15, effect size=0.25
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[1], k = k_vec[2], 
           nr.percluster = 1000:3000) # 1681 per arm


#* comparing impact of different k: ----

## infection endpoint, k=0.15, effect size=0.3
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
           treatment_effect = effectsizes[2], k = k_vec[2], 
           nr.percluster = 20:200) # 149 per arm

## infection endpoint, k=0.2, effect size=0.3
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
           treatment_effect = effectsizes[2], k = k_vec[1], 
           nr.percluster = 20:200) # 116 per arm (22% less)

## infection endpoint, k=0.25, effect size=0.3
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
           treatment_effect = effectsizes[2], k = k_vec[3], 
           nr.percluster = 20:310) # 307 per arm (106% more)


## case endpoint, k=0.15, effect size=0.3
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[2], k = k_vec[2], 
           nr.percluster = 20:2000) # 967 per arm

## case endpoint, k=0.2, effect size=0.35
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[2], k = k_vec[1], 
           nr.percluster = 20:2000) # 753 per arm (22% less)

## case endpoint, k=0.25, effect size=0.25
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[2], k = k_vec[3], 
           nr.percluster = 1000:3000) # 1994 per arm (106% more)

# sens and spec analysis ----------
library(lattice)

pi_0 <- inf_c_2yr

#pi_0_adj <- pi_0*sens + (1-pi_0)*(1-sp)

s_vals <- seq(0, 1, 0.05)
s_df <- expand.grid(s_vals, s_vals)

s_df$pi_0 <- pi_0

names(s_df)[names(s_df) == "Var1"] <- "sens"
names(s_df)[names(s_df) == "Var2"] <- "sp"

s_df$app_pi_0 <- numeric(nrow(s_df))

s_df$app_pi_0 <- (pi_0 * s_df$sens) + ((1-pi_0) * (1-s_df$sp))

#calc sample size based on app_pi (ie observed pi)
s_df$app_pi_nrpercl <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0, treatment_effect = 0.3, k=0.15, clusters_perarm=15)

#calc sample size based on true pi
truepi_nrpercl <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=pi_0, treatment_effect = 0.3, k=0.15, clusters_perarm=15)

#relative difference between the two sample sizes
s_df$nrpercl_fold <- (s_df$app_pi_nrpercl)/truepi_nrpercl

#log10
s_df$log10_nrpercl_fold <- log10(s_df$nrpercl_fold)

levelplot(log10_nrpercl_fold ~ sens*sp, s_df, col.regions=terrain.colors(100))

# levelplot(nrpercl_diff ~ sens*sp, s_df,
#           at=c(-1,-0.5,0,1,2,5,10,15,25),
#           colorkey=list((at=c(-1,-0.5,0,1,2,5,10,15,25)),
#           labels=list(at=c(-1,-0.5,0,1,2,5,10,15,25))))

pt.percl.inf <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.094, 
                           treatment_effect = effectsizes[2], k = k_vec[2], 
                           clusters_perarm=15)
pt.percl.case <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                            treatment_effect = effectsizes[2], k = k_vec[2], 
                            clusters_perarm=15)










