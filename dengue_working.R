# we want to recreate the dengue model

# we first create the four probabilities of infection 
getwd()
library(tidyverse)
library(ggplot2)

#============ infection P

lambda <- 0.01
max_age <- 100

N0 <- 100 

no_infection <- N0

col_length <- rep(0, max_age+1)

col_names <- list("age",
               "no_infected",
               "first_inf",
               "sec_inf",
               "third_inf",
               "fourth_inf",
               "new_first",
               "new_second",
               "new_third",
               "new_fourth",
               "total")

dengue <- array(data = 0,
                dim = c(length(col_length),length(col_names)))
colnames(dengue) <- col_names
dengue[,"age"] <- seq( 0, max_age, by = 1)
dengue[1,"no_infected"] <- 100

infection_probs <- function( lambda, max_age, df, verbose=1 ){
  if(verbose==1){browser()}
  
  # infection probabilities
  infect1_0 <- 1-(1- lambda)^4
  infect2_1 <- 1-(1- lambda)^3
  infect3_2 <- 1-(1- lambda)^2
  infect4_3 <- 1-(1- lambda)^1
  infect5_4 <- 0L
  
  
  # first column(Number of people who are susceptible)
  for( ii in 2:(max_age+1)) {
    
    df[ii, 2] <- df[ii -1, 2] * (1-infect1_0) # probability 
    
  }
  
  # second column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 3] <- df[ii -1, 2] * (infect1_0) + # prob of people with first infection 
      (df[ii-1, 3] * (1-infect2_1) ) # probability of not getting a second infection
    
    
    df[ii, 7] <- df[ii -1, 2] * (infect1_0)  # Number of new cases

    
  }
  
  # third column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 4] <- df[ii -1, 3] * (infect2_1) + # prob of people with second infection 
      (df[ii-1, 4] * (1-infect3_2) ) # probability of not getting a third infection
    
    df[ii, 8] <- df[ii -1, 3] * (infect2_1)   # Number of new cases
  }
  
  
  # fourth column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 5] <- df[ii -1, 4] * (infect3_2) + # prob of people with third infection 
      (df[ii-1, 5] * (1-infect4_3) ) # probability of not getting a fourth infection
    
    df[ii, 9] <- df[ii -1, 4] * (infect3_2) # number of new cases 
  }
  
  
  # fifth column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 6] <- df[ii -1, 5] * (infect4_3) + # prob of people with third infection 
      (df[ii-1, 6] * (1-infect5_4) ) # probability of not getting a fourth infection
    
    df[ii, 10] <- df[ii -1, 5] * (infect4_3) # number of new cases 
  }
  
  df[,"total"] <- rowSums(df[, 2:6])
  
  return(df)
}


dengue_lambda_Lo <- infection_probs( lambda = 0.01, max_age = max_age, df = dengue, verbose = 0 )

dengue_lambda_Hi <- infection_probs( lambda = 0.05, max_age = max_age, df = dengue, verbose = 0 )


plot( x= dengue_lambda_Lo[,"age"], y = dengue_lambda_Lo[,"new_second"], type = "l", ylim = c(0, 5))
lines( x = dengue_lambda_Hi[,"age"], y = dengue_lambda_Hi[,"new_second"])

 

