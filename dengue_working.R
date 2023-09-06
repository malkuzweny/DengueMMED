rm(list=ls())
# we want to recreate the dengue model

# we first create the four probabilities of infection
library(data.table)

# Set parameters ----------------------------------------------------------

# Calculate lambda in Colombo

# probability first infection = 14.1%
# 0.141 = 1-(1 - lambda)^3
# solve the above equation for lambda manually (eg Wolfram Alpha)
# lambda = 0.049


#will remove last age
max_age <- 61
#will remove last year
years <- 2010:2110
N0 <- 1 # population size
no_infection <- N0
lambda <- rep(0.049, length=length(years))
P_D <- c(0.005, 0.05, 0.00001, 0.00001)
#probability of disease given inf for 3 serotype model (vector needs to be length 3)
P_D <- P_D[1:3]

P_D_mat <- t(replicate(max_age,P_D))

# Initialise dataframe ----------------------------------------------------

col_names <- list(
               "S3",
               "S2",
               "S1",
               "S0",
               "I1",
               "I2",
               "I3",
               "D1",
               "D2",
               "D3",
               "total")

array_names <- vector("list", 3)
array_names[[1]] <- paste("year", years, sep="_")
array_names[[2]] <- paste("age", (1:(max_age+1) -1), sep="_")
array_names[[3]] <- col_names

# dims are year - age - outcomes
dengue <- array(data = 0,
                dim = c(length(years), max_age+1,length(col_names)),
                dimnames=array_names)


dengue[1,,"S3"] <- N0 # add 100 to S3 to each age group in first year
dengue[,1,"S3"] <- N0 # add 100 to S3 of 0 age group in all years (births)

dengue[1,,"total"] <- dengue[1,,"S3"] # add total to first year
dengue[,1,"total"] <- dengue[,1,"S3"] # add total to first age

# Populate dataframe ------------------------------------------------------

infection_probs <- function( lambda, max_age, P_D_mat, df, verbose=1 ){
  if(verbose==1){browser()}
  
  for (y in 1:(length(years)-1)) {

    # infection probabilities
    P_I0_I1 <- 1-(1 - lambda[y])^3
    P_I1_I2 <- 1-(1 - lambda[y])^2
    P_I2_I3 <- 1-(1 - lambda[y])^1
    P_I3_I4 <- 0L
    P_I4_I5 <- 0L
    
    # S4
    df[y+1, 1:max_age + 1 , "S3"] <- df[y, 1:max_age, "S3"] * (1-P_I0_I1) # probability of individual staying susceptible
    
    # S3
    df[y+1, 1:max_age + 1, "S2"] <- df[y, 1:max_age, "S3"] * (P_I0_I1) + # prob of people with first infection 
      (df[y, 1:max_age, "S2"] * (1-P_I1_I2) ) # probability of not getting a second infection
    
    # I1 (reminder: first infections occur during this year) - that's why there's no +1's
    df[y, 1:max_age, "I1"] <- df[y,1:max_age, "S3"] * (P_I0_I1)  # Number of new first infections
    
    # S2
    df[y+1, 1:max_age+1, "S1"] <- df[y, 1:max_age, "S2"] * (P_I1_I2) + # prob of people with second infection 
      df[y, 1:max_age, "S1"] * (1-P_I2_I3) # probability of not getting a third infection
    
    # I2
    df[y, 1:max_age, "I2"] <- df[y, 1:max_age, "S2"] * (P_I1_I2)   # Number of new second infections
    
    # S1
    df[y+1, 1:max_age+1, "S0"] <- df[y, 1:max_age, "S1"] * (P_I2_I3) + # prob of people with third infection 
      (df[y, 1:max_age, "S0"] * (1-P_I3_I4) ) # probability of not getting a fourth infection
    
    #I3
    df[y, 1:max_age, "I3"] <- df[y, 1:max_age, "S1"] * (P_I2_I3) # number of new third infections

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


# I2 vs. lambda plot ------------------------------------------------------

lambda_vals <- seq(0, 0.1, by=0.001)
I2_vec <- numeric(length(lambda_vals))
I1_vec <- numeric(length(lambda_vals))

AgeDist <- c(0.01518000, 0.01518000, 0.01518000, 0.01518000, 0.01518000, 0.01560000, 0.01560000, 0.01560000, 0.01560000, 0.01560000,
             0.01478000, 0.01478000, 0.01478000, 0.01478000, 0.01478000, 0.01617775, 0.01617775, 0.01617775, 0.01617775, 0.01617775,
             0.01506456, 0.01506456, 0.01506456, 0.01506456, 0.01506456, 0.01525617, 0.01525617, 0.01525617, 0.01525617, 0.01525617,
             0.01610475, 0.01610475, 0.01610475, 0.01610475, 0.01610475, 0.01383275, 0.01383275, 0.01383275, 0.01383275, 0.01383275)

AgeDist <- c(AgeDist, rep(AgeDist[length(AgeDist)],20)) # fill age vector to age 60

# obtain prop of people with second infection in a year at equilibrium
for(i in 1:length(lambda_vals)){
  
  dengue_df <- infection_probs(lambda = rep(lambda_vals[i], length(years)), 
                               max_age = max_age, 
                               P_D_mat = P_D_mat, 
                               df = dengue,
                               verbose = 0)
  
  I2_vec[i] <- sum(AgeDist[1:60]*dengue_df[length(years)-1,1:60,"I2"])/sum(AgeDist[1:60])
  I1_vec[i] <- sum(AgeDist[1:60]*dengue_df[length(years)-1,1:60,"I1"])/sum(AgeDist[1:60])
  
}

#if using 4-serotype model, use I2_vec[38] and lambda_vals[38]
plot(x = I2_vec/I2_vec[50], y = lambda_vals, type = "l", 
     xlab = "Modelled cases proportional to Colombo" , 
     ylab = "Force of Infection")
abline(v = 1.0)
abline(h = lambda_vals[50])
abline(v = 0.77)
#value of lambda where I2_vec/I2_vec[38] is closest to 0.77
abline(h = lambda_vals[21])


# plot proportion of people with second infection at equilibrium relative to prop of people
# with second infection in Colombo for a range of FOIs. this is used to find out the FOI
# in the region where prop of people with second infection is 77% of that in Colombo

font.size <- 2
lwd <- 4

plot(x = I2_vec/I2_vec[50], y = lambda_vals, type = "l", lwd=lwd,
     xlab="Secondary infections relative to Colombo", ylab="Force of infection", cex.axis=font.size, cex.lab=font.size)
abline(v = 1.0, lwd=lwd, lty=2)
abline(h = lambda_vals[50], lwd=lwd, lty=2)
abline(v = 0.77, lwd=lwd, lty=2)
#value of lambda where I2_vec/I2_vec[38] is closest to 0.77
abline(h = lambda_vals[21], lwd=lwd, lty=2)

#* Calculating rho ---------------------------------------------------------

rho <- (100/10000)/I2_vec[50]

# Calculating I1 and I2 in Gampaha ----------------------------------------

dengue_g_df <- infection_probs(lambda = rep(lambda_vals[21], length(years)), 
                               max_age = max_age, 
                               P_D_mat = P_D_mat, 
                               df = dengue,
                               verbose = 0)

#proportion of population who experience primary or secondary infections per year
ageRange <- 4:16
inf_c <- sum(AgeDist[ageRange]*
               rowSums(dengue_g_df[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_2yr <- 1-(1-inf_c)^2

case_c_2yr <- 1-(1-(77/10000))^2





