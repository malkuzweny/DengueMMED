# we want to recreate the dengue model

# we first create the four probabilities of infection 
getwd()
library(tidyverse)

# Set parameters ----------------------------------------------------------

max_age <- 100
N0 <- 100 # population size
no_infection <- N0


# Initialise dataframe ----------------------------------------------------

col_length <- rep(0, max_age+1)

col_names <- list("age",
               "S4",
               "S3",
               "S2",
               "S1",
               "S0",
               "I1",
               "I2",
               "I3",
               "I4",
               "total")

dengue <- array(data = 0,
                dim = c(length(col_length),length(col_names)))
colnames(dengue) <- col_names
dengue[,"age"] <- seq( 0, max_age, by = 1)
dengue[1,"no_infected"] <- 100


# Populate dataframe ------------------------------------------------------

infection_probs <- function( lambda, max_age, df, verbose=1 ){
  if(verbose==1){browser()}
  
  # infection probabilities
  infect1_0 <- 1-(1- lambda)^4
  infect2_1 <- 1-(1- lambda)^3
  infect3_2 <- 1-(1- lambda)^2
  infect4_3 <- 1-(1- lambda)^1
  infect5_4 <- 0L
  

  for(age in seq_len(max_age)) {
    # S4
    df[y+1, age + 1 , "S4"] <- df[y, age, "S4"] * (1-infect1_0) # probability of individual staying susceptible
    
    # S3
    df[y+1, age + 1, "S3"] <- df[y, age, "S4"] * (infect1_0) + # prob of people with first infection 
      (df[y, age, "S3"] * (1-infect2_1) ) # probability of not getting a second infection

    # I1 (reminder: infections occur during this age and year)
    df[y, age, "I1"] <- df[y,age, "S4"] * (infect1_0)  # Number of new first infections
    
    # S2
    df[age, "S2"] <- df[age -1, "S3"] * (infect2_1) + # prob of people with second infection 
      (df[age-1, "S2"] * (1-infect3_2) ) # probability of not getting a third infection
    
    # I2
    df[age, "I2"] <- df[age -1, "S3"] * (infect2_1)   # Number of new second infections
    
    # S1
    df[age, "S1"] <- df[age -1, "S2"] * (infect3_2) + # prob of people with third infection 
      (df[age-1, "S1"] * (1-infect4_3) ) # probability of not getting a fourth infection
    
    #I3
    df[age, "I3"] <- df[age -1, "S2"] * (infect3_2) # number of new third infections
    
    # S0
    df[age, "S0"] <- df[age -1, "S1"] * (infect4_3) + # prob of people with fourth infection 
      (df[age-1, "S0"] * (1-infect5_4) ) # probability of not getting a fifth infection
    
    # I4
    df[age, "I4"] <- df[age -1, "S1"] * (infect4_3) # number of new fourth infections
    
  }
  
  df[,"total"] <- rowSums(df[, "S4":"S0"])
  
  return(df)
}


# Plotting ----------------------------------------------------------------

dengue_lambda_Lo <- infection_probs( lambda = 0.01, max_age = max_age, df = dengue, verbose = 0 )

dengue_lambda_Hi <- infection_probs( lambda = 0.05, max_age = max_age, df = dengue, verbose = 0 )


plot( x= dengue_lambda_Lo[,"age"], y = dengue_lambda_Lo[,"new_second"], type = "l", ylim = c(0, 5))
lines( x = dengue_lambda_Hi[,"age"], y = dengue_lambda_Hi[,"new_second"])

 

