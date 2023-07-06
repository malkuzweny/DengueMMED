rm(list=ls())
# we want to recreate the dengue model

# we first create the four probabilities of infection 
getwd()
library(tidyverse)
library(data.table)

# Set parameters ----------------------------------------------------------

# Calculate lambda in Colombo

# probability first infection = 14.1%
# 0.141 <- 1-(1 - lambda)^4
# lambda = 0.037


#will remove last age
max_age <- 61
#will remove last year
years <- 2010:2110
N0 <- 1 # population size
no_infection <- N0
lambda <- rep(0.037, length=length(years))
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

###set baseline
dengue_baseline <- infection_probs( lambda = lambda, max_age = max_age, P_D_mat = P_D_mat, df = dengue, verbose = 0 )

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

dengue[1,1:max_age,"S4"] <- dengue_baseline[100,,"S4"] # add equilibrium susc to S4 to each age group in first year
dengue[1,1:max_age,"S3"] <- dengue_baseline[100,,"S3"] # add equilibrium susc to S3 to each age group in first year
dengue[1,1:max_age,"S2"] <- dengue_baseline[100,,"S2"] # add equilibrium susc to S2 to each age group in first year
dengue[1,1:max_age,"S1"] <- dengue_baseline[100,,"S1"] # add equilibrium susc to S1 to each age group in first year
dengue[1,1:max_age,"S0"] <- dengue_baseline[100,,"S0"] # add equilibrium susc to S0 to each age group in first year

dengue[,1,"S4"] <- N0 # add 100 to S4 of 0 age group in all years (births)

dengue[1,,"total"] <- dengue[1,,"S4"] + dengue[1,,"S3"] +
  dengue[1,,"S2"] + dengue[1,,"S1"] + 
  dengue[1,,"S0"] # add total to first year
dengue[,1,"total"] <- dengue[,1,"S4"] # add total to first age

#-------------------------

lambda_vals <- seq(0, 0.1, by=0.001)
I2_vec <- numeric(length(lambda_vals))

AgeDist <- c(0.01518000, 0.01518000, 0.01518000, 0.01518000, 0.01518000, 0.01560000, 0.01560000, 0.01560000, 0.01560000, 0.01560000,
             0.01478000, 0.01478000, 0.01478000, 0.01478000, 0.01478000, 0.01617775, 0.01617775, 0.01617775, 0.01617775, 0.01617775,
             0.01506456, 0.01506456, 0.01506456, 0.01506456, 0.01506456, 0.01525617, 0.01525617, 0.01525617, 0.01525617, 0.01525617,
             0.01610475, 0.01610475, 0.01610475, 0.01610475, 0.01610475, 0.01383275, 0.01383275, 0.01383275, 0.01383275, 0.01383275)

AgeDist <- c(AgeDist, rep(AgeDist[length(AgeDist)],20))

for(i in 1:length(lambda_vals)){
  dengue_df <- infection_probs( lambda = rep(lambda_vals[i], length(years)), 
                                max_age = max_age, 
                                P_D_mat = P_D_mat, 
                                df = dengue, 
                                verbose = 0 )
  
  I2_vec[i] <- sum(AgeDist[1:30]*dengue_df[length(years)-1,1:30,"I2"])
  
}

#set 1 for each age group (not 100)
#multiply I2 by the number of children in each age group
#max age is 60

plot(x = I2_vec/I2_vec[38], y = lambda_vals, type = "l")
abline(v = 1.0)
abline(h = lambda_vals[38])
abline(v = 0.77)
abline(h = lambda_vals[24])



#divide I2_vec by the I2_vec val when lambda = 0.037

lambda_vals[38]

#passive surveillance (cases = I2*rho, solve for rho) -- entire pop
#we know lambda --> I2 --> find rho to go to cases (known)
#active surveillance (D1 + D2 + D3 + D4) -- entire pop
#tracking infections (I1 + I2, pick age range and number of years)

#goal - get to:
#in 1 year -- among 4-16 y/o, 0.053 experience first or second infection (on avg across age)
#^weighted average by age
#agedist*sum in line 240/(agedist)
#in 2 years -- among 4-16 y/o, 0.103 experience first or second infection (on avg across age)
# ^^^ p (probability of an event in the control group)
# ^^^ effect size will be 30% (1-0.3)*(0.103)
#1-(1-0.053)^(# years)

#normal-distributed
#null_mu +/- 1.96*(sd/(sqrt(n))) = 95% CI
#P(X_a > cutoff) = power (known)
#n = (z_alpha/2 + z_beta)^2/(delta^2)
#z_beta = z for power beta
#delta = mu_null - mu_alt (aka effect size)

#for cRCTS
#c = 1+(z_alpha/2 + z_beta)^2*((lambda_0 + lambda_A)/y + k^2(lambda_0 + lambda_A)^2) / (lambda_0 - lambda_A)^2
#c = clusters per arm
#y = number of people per cluster
#k = between-cluster coefficient of variation
# more clusters w/ fewer people is preferable if outcomes correlated within a cluster

#cRCTs for binary outcomes
#c = 1 + (z_alpha/2 + z_beta)^2((pi_0*(1-pi_0)))


#things to vary -- effect size, k, n, cluster
#k = 0.15

##other ideas -- imbalanced, by events





