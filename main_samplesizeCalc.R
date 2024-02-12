rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')
source("main_infections.R")

# Parameters --------------------------------------------------------------

font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation for 3 serotype scenario -----

cl.perarm.inf_3st <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_3st, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.d_3st <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_3st, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# Do power calculation for 3 to 4 serotype scenario -----

cl.perarm.inf_3to4st <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_3to4st, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.d_3to4st <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_3to4st, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# Do power calculation for 4 serotype scenario -----

cl.perarm.inf_4st <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=inf_c_2yr_4st, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)
cl.perarm.d_4st <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=case_c_2yr_4st, treatment_effect = 0.3, k=0.15, clusters_perarm = 15)

# plot -----

cl_df <- data.frame(row.names=c("active", "passive"),
                    first=c(cl.perarm.inf_3st, cl.perarm.d_3st),
                    second=c(cl.perarm.inf_3to4st, cl.perarm.d_3to4st),
                    third=c(cl.perarm.inf_4st, cl.perarm.d_4st))
#cl_df <- do.call(rbind, cl_df)
barplot(as.matrix(cl_df), beside = TRUE, ylim=c(0,2000),
        names.arg=c("3 serotypes", "3 -> 4 serotypes", "4 serotypes"),
        xlab = "scenario", ylab = "no. of pt per cluster")

##plot relative difference in sample size by FOI scenario
cl_df$second_fold <- cl_df$second/cl_df$first
cl_df$third_fold <- cl_df$third/cl_df$first

cl_df$second_perc <- (cl_df$second-cl_df$first)/cl_df$first
cl_df$third_perc <- (cl_df$third-cl_df$first)/cl_df$first

barplot(as.matrix(cl_df[c("second_perc", "third_perc")]), 
        beside = TRUE,
        ylim=c(-0.2, 0.2),
        names.arg=c("3 -> 4 serotypes", "4 serotypes"),
        xlab = "scenario", ylab = "% change relative to baseline",
        legend=rownames(cl_df),
        args.legend = list(x = "topright", bty="n"))






