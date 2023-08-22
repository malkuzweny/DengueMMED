
font.size <- 1.5 #presentation used 2
lwd <- 3 #presentation used 4

# Do power calculation ----------------------------------------------------

run.sscalc <- function(z_a2, z_b, pi_0, treatment_effect, k, nr.percluster) {

  # pi_0 = proportion of population experiencing primary & second infection during trial at baseline
  # pi_a = proportion of population experiencing primary & second infection during trial in intervention group
  # k is between cluster variation coefficient
  
  pi_a <- pi_0 * (1-treatment_effect)
  
  clusters_perarm <- 2 + 
                    (z_a2 + z_b)^2 * 
                    ((pi_0 * (1-pi_0)/nr.percluster) + (pi_a * (1-pi_a)/nr.percluster) + k^2*(pi_0^2 + pi_a^2)) / (pi_0 - pi_a)^2 
  
  #clusters_perarm <- ceiling(clusters_perarm)
  # print(clusters_perarm)
  return(clusters_perarm)
}

cl.perarm.inf <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)

effectsizes <- c(0.25,0.3,0.35)
k_vec <- c(0.02,0.15,0.25)

plot_sample_infections <- function(){
  
cl.perarm.inf_k1_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[1], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k2_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[1], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k3_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[1], k = k_vec[3], 
                                  nr.percluster = 20:200)

layout(matrix(1:3, ncol=3, byrow=T))
par(oma=c(0,0,2,0))

plot(x=20:200, y=cl.perarm.inf_k1_e1, type="l", col="black", ylim=c(0,40), lwd=lwd, 
     xlab="Number of children per cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.25", cex.axis=font.size, cex.lab=font.size)
lines(x=20:200, y=cl.perarm.inf_k2_e1, col="red", lwd=lwd)
lines(x=20:200, y=cl.perarm.inf_k3_e1, col="green", lwd=lwd)

cl.perarm.inf_k1_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[2], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[2], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k3_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[2], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.inf_k1_e2, type="l", col="black", ylim=c(0,40), lwd=lwd,
     xlab="Number of children per cluster", ylab="",
     main="Effect size: 0.3", cex.axis=font.size, cex.lab=font.size)
lines(x=20:200, y=cl.perarm.inf_k2_e2, col="red", lwd=lwd)
lines(x=20:200, y=cl.perarm.inf_k3_e2, col="green", lwd=lwd)

cl.perarm.inf_k1_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[3], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k2_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[3], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k3_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[3], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.inf_k1_e3, type="l", col="black", ylim=c(0,40), lwd=lwd, 
     xlab="Number of children per cluster", ylab="",
     main="Effect size: 0.35",  cex.axis=font.size, cex.lab=font.size)
lines(x=20:200, y=cl.perarm.inf_k2_e3, col="red", lwd=lwd)
lines(x=20:200, y=cl.perarm.inf_k3_e3, col="green", lwd=lwd)
legend("topright", legend=c("k=0.02", 
                            "k=0.15",
                            "k=0.25"), col=c("black", "red", "green"), lty=1)

mtext("Subject enrollment using infections as events", line=0, side=3, outer=TRUE, cex=1.2)

}

# using cases as events
plot_sample_cases <- function(){
cl.perarm.case_k1_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[1], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.case_k2_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[1], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.case_k3_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[1], k = k_vec[3], 
                                  nr.percluster = 20:200)

layout(matrix(1:3, ncol=3, byrow=T))
par(oma=c(0,0,2,0))

plot(x=20:200, y=cl.perarm.case_k1_e1, type="l", col="black", lwd=lwd,
     xlab="Number of children per cluster", ylab="Number of clusters per arm", ylim=c(0,200), 
     main="Effect size: 0.25", cex.axis=font.size, cex.lab=font.size)
lines(x=20:200, y=cl.perarm.case_k2_e1, col="red", lwd=lwd)
lines(x=20:200, y=cl.perarm.case_k3_e1, col="green", lwd=lwd)

cl.perarm.case_k1_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.case_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.case_k3_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.case_k1_e2, type="l", col="black", lwd=lwd,
     xlab="Number of children per cluster", ylab="", ylim=c(0,200), 
     main="Effect size: 0.3", cex.axis=font.size, cex.lab=font.size)
lines(x=20:200, y=cl.perarm.case_k2_e2, col="red", lwd=lwd)
lines(x=20:200, y=cl.perarm.case_k3_e2, col="green", lwd=lwd)

cl.perarm.case_k1_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[3], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.case_k2_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[3], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.case_k3_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[3], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.case_k1_e3, type="l", col="black", lwd=lwd,
     xlab="Number of children per cluster", ylab="", ylim=c(0,200), 
     main="Effect size: 0.35", cex.axis=font.size, cex.lab=font.size)
lines(x=20:200, y=cl.perarm.case_k2_e3, col="red", lwd=lwd)
lines(x=20:200, y=cl.perarm.case_k3_e3, col="green", lwd=lwd)
legend("topright", legend=c("k=0.02", 
                            "k=0.15",
                            "k=0.25"), col=c("black", "red", "green"), lty=1)
 
 mtext("Subject enrollment using cases as events", line=0, side=3, outer=TRUE, cex=1.2)


}

plot_sample_cases_log10 <- function(){
  cl.perarm.case_k1_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[1], k = k_vec[1], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k2_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[1], k = k_vec[2], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k3_e1 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[1], k = k_vec[3], 
                                     nr.percluster = 20:200)
  
  layout(matrix(1:3, ncol=3, byrow=T))
  par(oma=c(0,0,2,0))
  
  plot(x=20:200, y=cl.perarm.case_k1_e1, type="l", col="black", 
       xlab="Number of children per cluster", ylab="Number of clusters per arm", 
       main="Effect size: 0.25",
       log="y")
  lines(x=20:200, y=cl.perarm.case_k2_e1, col="red")
  lines(x=20:200, y=cl.perarm.case_k3_e1, col="green")
  
  cl.perarm.case_k1_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[2], k = k_vec[1], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[2], k = k_vec[2], 
                                     nr.percluster = 20:2000)
  cl.perarm.case_k3_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[2], k = k_vec[3], 
                                     nr.percluster = 20:200)
  
  plot(x=20:200, y=cl.perarm.case_k1_e2, type="l", col="black", 
       xlab="Number of children per cluster", 
       main="Effect size: 0.3",
       log="y")
  lines(x=20:200, y=cl.perarm.case_k2_e2, col="red")
  lines(x=20:200, y=cl.perarm.case_k3_e2, col="green")
  
  cl.perarm.case_k1_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[3], k = k_vec[1], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k2_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[3], k = k_vec[2], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k3_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[3], k = k_vec[3], 
                                     nr.percluster = 20:200)
  
  plot(x=20:200, y=cl.perarm.case_k1_e3, type="l", col="black", 
       xlab="Number of children per cluster", 
       main="Effect size: 0.35",
       log="y")
  lines(x=20:200, y=cl.perarm.case_k2_e3, col="red")
  lines(x=20:200, y=cl.perarm.case_k3_e3, col="green")
  legend("topright", legend=c("k=0.02", 
                              "k=0.15",
                              "k=0.25"), col=c("black", "red", "green"), lty=1)
  
  mtext("Subject enrollment using cases as events", line=0, side=3, outer=TRUE, cex=1.2)
  
  
}


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

###############
  # comparing impact of different k:
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
