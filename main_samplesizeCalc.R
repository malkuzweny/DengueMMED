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
library(latticeExtra)

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

# sens and spec analysis v2 --------------
# allowing sn/sp to vary by inf valency --------

pi_0_inf_pr <- inf_pr
pi_0_inf_sec <- inf_sec

s_vals <- seq(0.5, 1, 0.05)
s_df <- expand.grid(s_vals, s_vals)

s_df$pi_0_inf_pr <- pi_0_inf_pr
s_df$pi_0_inf_sec <- pi_0_inf_sec

names(s_df)[names(s_df) == "Var1"] <- "sens_sec"
names(s_df)[names(s_df) == "Var2"] <- "sp_sec"

s_df$sens_pr_50 <- 0.5
s_df$sp_pr_50 <- 0.5

s_df$sens_pr_75 <- 0.75
s_df$sp_pr_75 <- 0.75

s_df$sens_pr_100 <- 1
s_df$sp_pr_100 <- 1

#sn = 50, sp = 50
s_df$app_pi_0_pr_1 <- (pi_0_inf_pr * s_df$sens_pr_50) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_50))
#sn = 50, sp = 75
s_df$app_pi_0_pr_2 <- (pi_0_inf_pr * s_df$sens_pr_50) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_75))
#sn = 50, sp = 100
s_df$app_pi_0_pr_3 <- (pi_0_inf_pr * s_df$sens_pr_50) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_100))

#sn = 75, sp = 50
s_df$app_pi_0_pr_4 <- (pi_0_inf_pr * s_df$sens_pr_75) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_50))
#sn = 75, sp = 75
s_df$app_pi_0_pr_5 <- (pi_0_inf_pr * s_df$sens_pr_75) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_75))
#sn = 75, sp = 100
s_df$app_pi_0_pr_6 <- (pi_0_inf_pr * s_df$sens_pr_75) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_100))

#sn = 100, sp = 50
s_df$app_pi_0_pr_7 <- (pi_0_inf_pr * s_df$sens_pr_100) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_50))
#sn = 100, sp = 75
s_df$app_pi_0_pr_8 <- (pi_0_inf_pr * s_df$sens_pr_100) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_75))
#sn = 100, sp = 100
s_df$app_pi_0_pr_9 <- (pi_0_inf_pr * s_df$sens_pr_100) + ((1-pi_0_inf_pr) * (1-s_df$sp_pr_100))

s_df$app_pi_0_sec <- (pi_0_inf_sec * s_df$sens_sec) + ((1-pi_0_inf_sec) * (1-s_df$sp_sec))

s_df$app_pi_0_1_total <- s_df$app_pi_0_pr_1 + s_df$app_pi_0_sec
s_df$app_pi_0_2_total <- s_df$app_pi_0_pr_2 + s_df$app_pi_0_sec
s_df$app_pi_0_3_total <- s_df$app_pi_0_pr_3 + s_df$app_pi_0_sec
s_df$app_pi_0_4_total <- s_df$app_pi_0_pr_4 + s_df$app_pi_0_sec
s_df$app_pi_0_5_total <- s_df$app_pi_0_pr_5 + s_df$app_pi_0_sec
s_df$app_pi_0_6_total <- s_df$app_pi_0_pr_6 + s_df$app_pi_0_sec
s_df$app_pi_0_7_total <- s_df$app_pi_0_pr_7 + s_df$app_pi_0_sec
s_df$app_pi_0_8_total <- s_df$app_pi_0_pr_8 + s_df$app_pi_0_sec
s_df$app_pi_0_9_total <- s_df$app_pi_0_pr_9 + s_df$app_pi_0_sec

s_df$app_pi_0_1_total_2yr <- 1-(1-s_df$app_pi_0_1_total)^2
s_df$app_pi_0_2_total_2yr <- 1-(1-s_df$app_pi_0_2_total)^2
s_df$app_pi_0_3_total_2yr <- 1-(1-s_df$app_pi_0_3_total)^2
s_df$app_pi_0_4_total_2yr <- 1-(1-s_df$app_pi_0_4_total)^2
s_df$app_pi_0_5_total_2yr <- 1-(1-s_df$app_pi_0_5_total)^2
s_df$app_pi_0_6_total_2yr <- 1-(1-s_df$app_pi_0_6_total)^2
s_df$app_pi_0_7_total_2yr <- 1-(1-s_df$app_pi_0_7_total)^2
s_df$app_pi_0_8_total_2yr <- 1-(1-s_df$app_pi_0_8_total)^2
s_df$app_pi_0_9_total_2yr <- 1-(1-s_df$app_pi_0_9_total)^2

#calc sample size based on app_pi (ie observed pi)
s_df$app_pi_nrpercl_1 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_1_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_2 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_2_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_3 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_3_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_4 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_4_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_5 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_5_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_6 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_6_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_7 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_7_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_8 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_8_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)
s_df$app_pi_nrpercl_9 <- 
  run.sscalc(z_a2=1.96, z_b=0.84, pi_0=s_df$app_pi_0_9_total_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)

#calc sample size based on true pi
truepi_nrpercl <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr, treatment_effect = 0.3, k=0.15, clusters_perarm=15)

#relative difference between the two sample sizes
s_df$nrpercl_fold <- (s_df$app_pi_nrpercl)/truepi_nrpercl

#log10
s_df$log10_nrpercl_fold <- log10(s_df$nrpercl_fold)

levelplot(log10_nrpercl_fold ~ sens*sp, s_df, col.regions=terrain.colors(100))

#layout(matrix(c(3,6,9,2,5,8,1,4,7), nrow=3, byrow=T))

#x axis = sn, y axis = sp
plot1 <- levelplot(log10(app_pi_nrpercl_1/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("50%"), ylab=c("50%"))

plot2 <- levelplot(log10(app_pi_nrpercl_2/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("50%"), ylab=c("75%"))

plot3 <- levelplot(log10(app_pi_nrpercl_3/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("50%"), ylab=c("100%"))

plot4 <- levelplot(log10(app_pi_nrpercl_4/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("75%"), ylab=c("50%"))

plot5 <- levelplot(log10(app_pi_nrpercl_5/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("75%"), ylab=c("75%"))

plot6 <- levelplot(log10(app_pi_nrpercl_6/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("75%"), ylab=c("100%"))

plot7 <- levelplot(log10(app_pi_nrpercl_7/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("100%"), ylab=c("50%"))

plot8 <- levelplot(log10(app_pi_nrpercl_8/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("100%"), ylab=c("75%"))

plot9 <- levelplot(log10(app_pi_nrpercl_9/truepi_nrpercl) ~ sens_sec*sp_sec, s_df, 
                   cuts=10,
                   at=seq(-2,0.5,by=0.25),
                   colorkey=list((at=seq(-2,0.5,by=0.25)),
                                 labels=list(at=seq(-2,0.5,by=0.25))),
                   col.regions=c("#9E0142","#D53E4F","#F46D43","#FDAE61",
                                 "#FEE08B","#FFFFBF","#E6F598","#ABDDA4",
                                 "#66C2A5","#48a36c","#5E4FA2"),
                   xlab=c("100%"), ylab=c("100%"))

full_plot <- c(plot1, plot4, plot7, plot2, plot5, plot8, plot3, plot6, plot9,
               merge.legends=F)

update(full_plot, 
       xlab="Sensitivity (secondary)\n
       50%                        75%                        100%\n
       Sensitivity (primary)", 
       ylab="Specificity (primary)\n
       50%                        75%                        100%\n
       Specificity (secondary)")
#save as 650x600










