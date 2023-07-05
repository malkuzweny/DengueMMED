#rm(list=ls())
# we want to recreate the dengue model

# we first create the four probabilities of infection 
getwd()
library(tidyverse)
library(data.table)

# Set parameters ----------------------------------------------------------

#will remove last age
max_age <- 100
#will remove last year
years <- 2010:2030
N0 <- 100 # population size
no_infection <- N0
lambda <- rep(0.01, length=length(years))
P_D <- c(0.005, 0.05, 0.00001, 0.00001)

P_D_mat <- t(replicate(max_age,P_D))

# Initialise dataframe ----------------------------------------------------

col_names <- list(
               "S4",
               "S3",
               "S2",
               "S1",
               "S0",
               "I1",
               "I2",
               "I3",
               "I4",
               "D1",
               "D2",
               "D3",
               "D4",
               "total")

array_names <- vector("list", 3)
array_names[[1]] <- paste("year", years, sep="_")
array_names[[2]] <- paste("age", (1:(max_age+1) -1), sep="_")
array_names[[3]] <- col_names

# dims are year - age - outcomes
dengue <- array(data = 0,
                dim = c(length(years), max_age+1,length(col_names)),
                dimnames=array_names)



dengue[1,,"S4"] <- N0 # add 100 to S4 to each age group in first year
dengue[,1,"S4"] <- N0 # add 100 to S4 of 0 age group in all years (births)

dengue[1,,"total"] <- dengue[1,,"S4"] # add total to first year
dengue[,1,"total"] <- dengue[,1,"S4"] # add total to first age

# Populate dataframe ------------------------------------------------------

infection_probs <- function( lambda, max_age, P_D_mat, df, verbose=1 ){
  if(verbose==1){browser()}
  
  for (y in 1:(length(years)-1)) {

    # infection probabilities
    P_I0_I1 <- 1-(1 - lambda[y])^4
    P_I1_I2 <- 1-(1 - lambda[y])^3
    P_I2_I3 <- 1-(1 - lambda[y])^2
    P_I3_I4 <- 1-(1 - lambda[y])^1
    P_I4_I5 <- 0L
    
    # S4
    df[y+1, 1:max_age + 1 , "S4"] <- df[y, 1:max_age, "S4"] * (1-P_I0_I1) # probability of individual staying susceptible
    
    # S3
    df[y+1, 1:max_age + 1, "S3"] <- df[y, 1:max_age, "S4"] * (P_I0_I1) + # prob of people with first infection 
      (df[y, 1:max_age, "S3"] * (1-P_I1_I2) ) # probability of not getting a second infection
    
    # I1 (reminder: first infections occur during this year) - that's why there's no +1's
    df[y, 1:max_age, "I1"] <- df[y,1:max_age, "S4"] * (P_I0_I1)  # Number of new first infections
    
    # S2
    df[y+1, 1:max_age+1, "S2"] <- df[y, 1:max_age, "S3"] * (P_I1_I2) + # prob of people with second infection 
      df[y, 1:max_age, "S2"] * (1-P_I2_I3) # probability of not getting a third infection
    
    # I2
    df[y, 1:max_age, "I2"] <- df[y, 1:max_age, "S3"] * (P_I1_I2)   # Number of new second infections
    
    # S1
    df[y+1, 1:max_age+1, "S1"] <- df[y, 1:max_age, "S2"] * (P_I2_I3) + # prob of people with third infection 
      (df[y, 1:max_age, "S1"] * (1-P_I3_I4) ) # probability of not getting a fourth infection
    
    #I3
    df[y, 1:max_age, "I3"] <- df[y, 1:max_age, "S2"] * (P_I2_I3) # number of new third infections
    
    # S0
    df[y+1, 1:max_age+1, "S0"] <- df[y, 1:max_age, "S1"] * (P_I3_I4) + # prob of people with fourth infection 
      (df[y, 1:max_age, "S0"] * (1-P_I4_I5) ) # probability of not getting a fifth infection
    
    # I4
    df[y, 1:max_age, "I4"] <- df[y, 1:max_age, "S1"] * (P_I3_I4) # number of new fourth infections
    
    #calculate the total population across susceptible statuses
    df[y+1,1:max_age+1,"total"] <- rowSums(df[y+1,1:max_age+1,which(col_names %like% "S")])
    
    df[y,1:max_age,which(col_names %like% "D")] <- df[y,1:max_age,which(col_names %like% "I")]*
      P_D_mat
    
  } # ends year loop
  
  #removing last year and last age group
  df <- df[-length(years),,]
  df <- df[,-(max_age+1),]
  
  return(df)
}

# Plotting ----------------------------------------------------------------

dengue_lambda_Lo <- infection_probs( lambda = lambda, max_age = max_age, P_D_mat = P_D_mat, df = dengue, verbose = 0 )

#old_array <- dengue_lambda_Lo

dengue_lambda_Lo[, , "total"]

dengue_lambda_Hi <- infection_probs( lambda = 0.05, max_age = max_age, df = dengue, verbose = 0 )

dengue_lambda_Lo["2010",,"I2"]

plot( x=1:(length(years)-1), y = dengue_lambda_Lo[,3,"D2" ], type = "b", 
      ylim = c(0, 0.02), xlab="Years", ylab="Number of new clinical secondary infections")

for(i in 2:(length(years)-1)){
  lines( x = (1:(max_age) -1), y = dengue_lambda_Lo[i,,"I1" ])
  
}


