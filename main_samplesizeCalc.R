rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('functions.R')

# Parameters --------------------------------------------------------------

font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation ----------------------------------------------------
# Pi_0 equals inf_c (proportion of population experiencing the outcome) or case_c_2yr

cl.perarm.inf <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

effectsizes <- c(0.25,0.3,0.35)
k_vec <- c(0.02,0.15,0.25)



# Different age groups ----------------------------------------------------

# * Infection (active surv) -----------------------------------------------

cl.perarm.inf <- vector()

for (ii in prop.infection$prop.infection) {
  cl.perarm.inf.ii <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=ii, treatment_effect = 0.3, k=0.15, nr.percluster = 200)
  cl.perarm.inf <- rbind(cl.perarm.inf, cl.perarm.inf.ii)
}

prop.infection$cl.perarm.inf <- cl.perarm.inf
prop.infection$sampleSize.inf <- prop.infection$cl.perarm.inf * 200
prop.infection$prop.population <- prop.infection$sampleSize.inf / prop.infection$population * 2 # *2 bc two arms
names(prop.infection) <- c("prop.infection", "minAge", "population", "maxAge", "clusters.perarm", "people.perarm", "prop.population")


plot(x=prop.infection$minAge, y=prop.infection$clusters.perarm, type="l", col="black", 
     xlab="Minimum age", ylab="Number of clusters per arm",
     main="Effect size: 0.3, k:0.15, nr.percluster:200")

plot(x=prop.infection$minAge, y=prop.infection$people.perarm, type="l", col="black", 
     xlab="Minimum age", ylab="Number of participants per arm",
     main="Effect size: 0.3, k:0.15, nr.percluster:200")

plot(x=prop.infection$minAge, y=prop.infection$prop.population, type="l", col="black", 
     xlab="Minimum age", ylab="Proportion of total population included in study",
     main="Effect size: 0.3, k:0.15, nr.percluster:200")


# only people aged 20+ is about 4500 per arm -> 9000 in total. we had 4500 in report (with 4 serotypes, and fewer per cluster)

# using pre-specified age groups

cl.perarm.inf <- vector()

for (ii in prop.infection$prop.infection) {
  cl.perarm.inf.ii <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=ii, treatment_effect = 0.3, k=0.15, nr.percluster = 200)
  cl.perarm.inf <- rbind(cl.perarm.inf, cl.perarm.inf.ii)
}

prop.infection$cl.perarm.inf <- cl.perarm.inf
prop.infection$sampleSize.inf <- prop.infection$cl.perarm.inf * 200
prop.infection$prop.population <- prop.infection$sampleSize.inf / prop.infection$population * 2 # *2 bc two arms
# names(prop.infection) <- c("prop.infection", "minAge", "population", "maxAge", "clusters.perarm", "people.perarm", "prop.population")

ggplot(data=prop.infection) +
  geom_tile(aes(x=minAge, y=maxAge, fill=sampleSize.inf)) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)) +
  scale_y_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank())   +
  scale_fill_gradient2(midpoint=mean(prop.infection$sampleSize.inf), low="blue", high="red", mid="white") +
  labs(x='Minimum age', y="Maximum age", fill="Required no. of participants")


ggplot(data=prop.infection) +
  geom_tile(aes(x=minAge, y=maxAge, fill=prop.population)) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)) +
  scale_y_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #color = "gray", linetype = "dashed"),
        panel.grid.minor = element_blank())

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
