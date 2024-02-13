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

# plot -----

cl_df <- data.frame(row.names=c("active", "passive"),
                    first=c(cl.perarm.inf_3st, cl.perarm.d_3st),
                    second=c(cl.perarm.inf_3to4st, cl.perarm.d_3to4st),
                    third=c(cl.perarm.inf_4st, cl.perarm.d_4st))
#cl_df <- do.call(rbind, cl_df)
barplot(as.matrix(cl_df[c("first", "second")]), beside = TRUE, ylim=c(0,1000),
        names.arg=c("3 serotypes\n(baseline)", "4th serotype emerges\nduring trial"),
        xlab = "Scenario", ylab = "Participants per cluster",
        legend=rownames(cl_df),
        args.legend = list(x = "topright", bty="n"))

##plot relative difference in sample size by scenario
cl_df$second_fold <- cl_df$second/cl_df$first
cl_df$third_fold <- cl_df$third/cl_df$first

cl_df$second_perc <- (cl_df$second-cl_df$first)/cl_df$first
cl_df$third_perc <- (cl_df$third-cl_df$first)/cl_df$first

barplot(as.matrix(cl_df[c("second_perc")]), 
        beside = TRUE,
        ylim=c(-0.2, 0.2),
        names.arg=c("3 -> 4 serotypes"),
        xlab = "scenario", ylab = "% change relative to baseline",
        legend=rownames(cl_df),
        axes=F,
        args.legend = list(x = "topright", bty="n"))
axis(side = 2, at = seq(-0.2, 0.2, 0.1), labels = paste0(seq(-0.2, 0.2, 0.1)*100, "%"))

# calculate relative power by scenario ------
#given baseline sample size
#3 serotype scenario
power.inf_3st <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_3st, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_3st, clusters_perarm = 15)
power.d_3st <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_3st, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.d_3st, clusters_perarm = 15)

#3->4 serotype scenario
power.inf_3to4st <- run.sscalc(z_a2=1.96, pi_0=inf_c_2yr_3to4st, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.inf_3st, clusters_perarm = 15)
power.d_3to4st <- run.sscalc(z_a2=1.96, pi_0=case_c_2yr_3to4st, treatment_effect = 0.3, k=0.15, nr.percluster=cl.perarm.d_3st, clusters_perarm = 15)

#for active surveillance (inf) study would be slightly overpowered if sample size according to baseline scenario
#for passive surveillance (dis) study would be underpowered if sample size according to baseline scenario

cl_df$baseline_power <- c(power.inf_3st, power.d_3st)
cl_df$second_power <- c(power.inf_3to4st, power.d_3to4st)

cl_df$second_power_perc <- (cl_df$second_power-cl_df$baseline_power)/cl_df$baseline_power

barplot(as.matrix(cl_df[c("baseline_power", "second_power")]), 
        beside = TRUE,
        ylim=c(0,1),
        names.arg=c("3 serotypes\n (baseline)", "4th serotype emerges\nduring trial"),
        xlab = "Scenario", ylab = "Power with baseline sample size",
        legend=rownames(cl_df),
        args.legend = list(bty="n", x=6, y=1.2))

barplot(as.matrix(cl_df[c("second_power_perc")]), 
        beside = TRUE,
        names.arg=c("3 -> 4 serotypes"),
        xlab = "scenario", ylab = "% diff in power relative to baseline",
        legend=rownames(cl_df),
        axes=F,
        args.legend = list(x = "topright", bty="n"))
axis(side = 2, at = seq(0, 0.08, 0.01), labels = paste0(seq(0, 0.08, 0.01)*100, "%"))

#plot percentage difference
tiff(filename = "~/Downloads/PerkinsLab/DengueTrial/st_emerge.tiff", 
    width = 600, height = 300)

{
  layout(matrix(1:2, ncol=2, byrow=T))
  
  par(omi = c(0, 0, 0, 1.1),
      xpd=TRUE)
  
  barplot(as.matrix(cl_df[c("second_perc")]), 
          beside = TRUE,
          ylim=c(-0.2, 0.2),
          names.arg=c(""),
          xlab = "", ylab = "% difference relative to baseline",
          main="Participants per cluster",
          axes=F)
  axis(side = 2, at = seq(-0.2, 0.2, 0.1), labels = paste0(seq(-0.2, 0.2, 0.1)*100, "%"))
  
  barplot(as.matrix(cl_df[c("second_power_perc")]), 
          beside = TRUE,
          names.arg=c(""),
          xlab = "", ylab = "% difference relative to baseline",
          main="Power",
          axes=F)
  axis(side = 2, at = seq(0, 0.08, 0.01), labels = paste0(seq(0, 0.08, 0.01)*100, "%"))
  
  legend("topright", inset = c(-0.9, 0), 
         legend=c("Active", "Passive"),
         pch=c(15,15), title="Surveillance type",
         col=grey.colors(2),
         bty="n", xpd=NA)
  
  mtext("Fourth serotype emerges during trial", side=3, line=-1.5, outer=T, cex=1.5, font=2)
}

dev.off()



