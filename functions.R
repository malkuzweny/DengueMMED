
# Calculate no of infections by age and year ------------------------------------------------------

infection_probs <- function( lambda, max_age, P_D_mat, output.df, verbose=1 ){
  if(verbose==1){browser()}
  
  for (y in 1:(length(years)-1)) {

    # infection probabilities
    P_I0_I1 <- 1-(1 - lambda[y])^3
    P_I1_I2 <- 1-(1 - lambda[y])^2
    P_I2_I3 <- 1-(1 - lambda[y])^1
    P_I3_I4 <- 0L
    P_I4_I5 <- 0L
    
    # S4
    output.df[y+1, 1:max_age + 1 , "S3"] <- output.df[y, 1:max_age, "S3"] * (1-P_I0_I1) # probability of individual staying susceptible
    
    # S3
    output.df[y+1, 1:max_age + 1, "S2"] <- output.df[y, 1:max_age, "S3"] * (P_I0_I1) + # prob of people with first infection 
      (output.df[y, 1:max_age, "S2"] * (1-P_I1_I2) ) # probability of not getting a second infection
    
    # I1 (reminder: first infections occur during this year) - that's why there's no +1's
    output.df[y, 1:max_age, "I1"] <- output.df[y,1:max_age, "S3"] * (P_I0_I1)  # Number of new first infections
    
    # S2
    output.df[y+1, 1:max_age+1, "S1"] <- output.df[y, 1:max_age, "S2"] * (P_I1_I2) + # prob of people with second infection 
      output.df[y, 1:max_age, "S1"] * (1-P_I2_I3) # probability of not getting a third infection
    
    # I2
    output.df[y, 1:max_age, "I2"] <- output.df[y, 1:max_age, "S2"] * (P_I1_I2)   # Number of new second infections
    
    # S1
    output.df[y+1, 1:max_age+1, "S0"] <- output.df[y, 1:max_age, "S1"] * (P_I2_I3) + # prob of people with third infection 
      (output.df[y, 1:max_age, "S0"] * (1-P_I3_I4) ) # probability of not getting a fourth infection
    
    #I3
    output.df[y, 1:max_age, "I3"] <- output.df[y, 1:max_age, "S1"] * (P_I2_I3) # number of new third infections

    #calculate the total population across susceptible statuses
    output.df[y+1,1:max_age+1,"total"] <- rowSums(output.df[y+1,1:max_age+1,which(col_names %like% "S")])
    
    output.df[y,1:max_age,which(col_names %like% "D")] <- output.df[y,1:max_age,which(col_names %like% "I")]*
      P_D_mat
    
  } # ends year loop
  
  #removing last year and last age group
  output.df <- output.df[-length(years),,]
  output.df <- output.df[,-(max_age+1),]
  
  return(output.df)
}



# Sample size calculations ------------------------------------------------

run.sscalc <- function(z_a2, z_b, pi_0, treatment_effect, k, nr.percluster=NULL, clusters_perarm=NULL) {
  
  # pi_0 = proportion of population experiencing primary & second infection during trial at baseline
  # pi_a = proportion of population experiencing primary & second infection during trial in intervention group
  # k is between cluster variation coefficient
  
  pi_a <- pi_0 * (1-treatment_effect)
  
  if(is.null(clusters_perarm)){
    clusters_perarm <- 2 + 
      (z_a2 + z_b)^2 * 
      ((pi_0 * (1-pi_0)/nr.percluster) + (pi_a * (1-pi_a)/nr.percluster) + k^2*(pi_0^2 + pi_a^2)) / (pi_0 - pi_a)^2 
    return(clusters_perarm)
  }
  
  if(is.null(nr.percluster)){
    nr.percluster <- (pi_0*(1-pi_0) + pi_a*(1-pi_a))/
      ((((clusters_perarm-2)*(pi_0-pi_a)^2)/(z_a2 + z_b)^2) - k^2*(pi_0^2 + pi_a^2))
    return(nr.percluster)
  }
  
}


# * plots -----------------------------------------------------------------

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




