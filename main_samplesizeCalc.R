rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')

# Parameters --------------------------------------------------------------

font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation ----------------------------------------------------

cl.perarm.inf <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

effectsizes <- c(0.25,0.3,0.35)
k_vec <- c(0.02,0.15,0.25)



# plots for presentations -------------------------------------------------
# infection
par(mar = c(5.1, 5.1, 4.1, 2.1))

plot(x=20:200, y=cl.perarm.inf_k2_e2, type="l", col="black", ylim=c(0,20), lwd=lwd,
     xlab="Number of children per cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.3", cex.axis=font.size, cex.lab=font.size, cex.main=font.size)
abline(h=15, lwd=lwd, lty=2)
abline(v=148, lwd=lwd, lty=2)

# disease
plot(x=20:2000, y=cl.perarm.case_k2_e2, type="l", col="black", lwd=lwd,
     xlab="Number of people per cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.3", cex.axis=font.size, cex.lab=font.size, cex.main=font.size,
     log="y", ylim=c(10,500), xlim=c(0,1000))
abline(h=15, lwd=lwd, lty=2)
abline(v=980, lwd=lwd, lty=2)


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

## infection endpoint, k=0.02, effect size=0.3
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

## case endpoint, k=0.02, effect size=0.35
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[2], k = k_vec[1], 
           nr.percluster = 20:2000) # 753 per arm (22% less)

## case endpoint, k=0.25, effect size=0.25
run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
           treatment_effect = effectsizes[2], k = k_vec[3], 
           nr.percluster = 1000:3000) # 1994 per arm (106% more)
