
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

plot(x=20:200, y=cl.perarm.inf_k1_e1, type="l", col="black", ylim=c(0,20), 
     xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.25")
lines(x=20:200, y=cl.perarm.inf_k2_e1, col="red")
lines(x=20:200, y=cl.perarm.inf_k3_e1, col="green")

cl.perarm.inf_k1_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[2], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[2], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k3_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[2], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.inf_k1_e2, type="l", col="black", ylim=c(0,20), 
     xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.3")
lines(x=20:200, y=cl.perarm.inf_k2_e2, col="red")
lines(x=20:200, y=cl.perarm.inf_k3_e2, col="green")

cl.perarm.inf_k1_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[3], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k2_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[3], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.inf_k3_e3 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, 
                                  treatment_effect = effectsizes[3], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.inf_k1_e3, type="l", col="black", ylim=c(0,20), 
     xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.35")
lines(x=20:200, y=cl.perarm.inf_k2_e3, col="red")
lines(x=20:200, y=cl.perarm.inf_k3_e3, col="green")
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

plot(x=20:200, y=cl.perarm.case_k1_e1, type="l", col="black", 
     xlab="Number of children\nper cluster", ylab="Number of clusters per arm", ylim=c(0,20), 
     main="Effect size: 0.25")
lines(x=20:200, y=cl.perarm.case_k2_e1, col="red")
lines(x=20:200, y=cl.perarm.case_k3_e1, col="green")

cl.perarm.case_k1_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[1], 
                                  nr.percluster = 20:200)
cl.perarm.case_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[2], 
                                  nr.percluster = 20:200)
cl.perarm.case_k3_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                  treatment_effect = effectsizes[2], k = k_vec[3], 
                                  nr.percluster = 20:200)

plot(x=20:200, y=cl.perarm.case_k1_e2, type="l", col="black", 
     xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.3")
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
     xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
     main="Effect size: 0.35")
lines(x=20:200, y=cl.perarm.case_k2_e3, col="red")
lines(x=20:200, y=cl.perarm.case_k3_e3, col="green")
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
       xlab="Number of children\nper cluster", ylab="Number of clusters per arm", 
       main="Effect size: 0.25",
       log="y")
  lines(x=20:200, y=cl.perarm.case_k2_e1, col="red")
  lines(x=20:200, y=cl.perarm.case_k3_e1, col="green")
  
  cl.perarm.case_k1_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[2], k = k_vec[1], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k2_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[2], k = k_vec[2], 
                                     nr.percluster = 20:200)
  cl.perarm.case_k3_e2 <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, 
                                     treatment_effect = effectsizes[2], k = k_vec[3], 
                                     nr.percluster = 20:200)
  
  plot(x=20:200, y=cl.perarm.case_k1_e2, type="l", col="black", 
       xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
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
       xlab="Number of children\nper cluster", ylab="Number of clusters per arm",
       main="Effect size: 0.35",
       log="y")
  lines(x=20:200, y=cl.perarm.case_k2_e3, col="red")
  lines(x=20:200, y=cl.perarm.case_k3_e3, col="green")
  legend("topright", legend=c("k=0.02", 
                              "k=0.15",
                              "k=0.25"), col=c("black", "red", "green"), lty=1)
  
  mtext("Subject enrollment using cases as events", line=0, side=3, outer=TRUE, cex=1.2)
  
  
}

