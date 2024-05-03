
rm(list = ls())

#change # of weeks to 104
#change how # of households are calculated for passive surveillance
#need to account for age dist and proportion S3/S2


#functions
passive_surveillance_cost <- function(sample_size, intervention_cost, 
                                      length_of_trial, duration_of_intervention, pt_per_hh){
  
  #number of households in the trial based on sample size
  no_hh <- sample_size/pt_per_hh
  
  #number of interventions needed for the entire duration of the trial per household
  #"intervention" = cost of installing all required shields depending on house size 
  #does not explicitly account for # of shields
  #accounts for # of times intervention needs to be replaced
  no_intervention_needed_per_hh <- length_of_trial/duration_of_intervention
  
  cost <- intervention_cost*no_hh*no_intervention_needed_per_hh
  return(cost)
}


###need to edit pt_per_hh - not 3.7
##age dist*(S3+S2)*3.7
active_surveillance_cost <- function(sample_size, intervention_cost, 
                                     length_of_trial, duration_of_intervention, pt_per_hh,
                                     test_cost_multiplier, freq_of_test){
  
  #number of households in the trial based on sample size
  no_hh <- sample_size/pt_per_hh
  
  #number of interventions needed for the entire duration of the trial per household
  no_intervention_needed_per_hh <- length_of_trial/duration_of_intervention
  
  #total cost of interventions needed
  total_int_cost <- intervention_cost*no_hh*no_intervention_needed_per_hh
  
  #total cost of each test needed
  #in terms of cost of intervention
  #freq_of_test = number of times test performed during entire length of trial
  total_test_cost <- (test_cost_multiplier*intervention_cost)*freq_of_test*sample_size
  
  cost <- total_int_cost + total_test_cost
  return(cost)
}

test_cost_multiplier_vec <- seq(1,50,1)
pt_ratio_vec <- seq(1,20,1)

mat <- matrix(data=0, nrow=length(test_cost_multiplier_vec), 
              ncol=length(pt_ratio_vec))
#mat[1,] <- pt_ratio_vec
#mat[,1] <- test_cost_multiplier_vec

for(i in 1:(length(test_cost_multiplier_vec))){
  for(j in 1:(length(pt_ratio_vec))){
    
    active_cost <- active_surveillance_cost(1000, 1,
                                            24*4, 4, 3.7,
                                            test_cost_multiplier_vec[i], 3)
    passive_cost <- passive_surveillance_cost(1000*pt_ratio_vec[j],1,
                                              24*4, 4, 3.7)
    
    mat[i,j] <- passive_cost/active_cost
      
  }
}

library(colorRamps)

filled.contour(x = pt_ratio_vec, y = test_cost_multiplier_vec, z = t(mat), 
               nlevels=18,
               col = c(colorRampPalette(c("brown1", "orange", "gold", "yellow", "yellowgreen", "green3", "green4"))(14)),
               key.title = {
                 title(main = c("Fold\ndifference\nin cost"), line=0.7, cex.main=1)
               },
               plot.axes = {
                 axis(1, cex.axis = 1.2)
                 axis(2, cex.axis = 1.2)
                 contour(x = pt_ratio_vec, y = test_cost_multiplier_vec, z = t(mat),
                         levels=1, labels = "1.0",
                         vfont=c("sans serif", "bold"),
                         labcex=1.3,
                         add = TRUE, lwd=1, lty=2)
               },
               plot.title = {
                 title(xlab="Ratio of participants between passive and active",
                       cex.lab=1.2)
                 title(ylab="Ratio of cost of test to intervention",
                       cex.lab=1.2)
               })







