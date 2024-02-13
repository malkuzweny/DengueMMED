rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(RColorBrewer)

source('functions.R')

# Set parameters ----------------------------------------------------------

# Calculate lambda in Colombo

# probability first infection = 14.1%
# 0.141 = 1-(1 - lambda)^3
# solve the above equation for lambda manually (eg Wolfram Alpha)
# lambda = 0.049


# will remove last age
max_age <- 61
# will remove last year
years <- 2010:2110
# population size
N0 <- 1 
no_infection <- N0
#FOI
lambda <- rep(0.049, length=length(years))
#probability of disease given infection
P_D <- c(0.005, 0.05, 0.00001, 0.00001)
#probability of disease given inf for 3 serotype model (vector needs to be length 3)
P_D <- P_D[1:3]

#create matrix of probability of disease with one line per age 
P_D_mat <- t(replicate(max_age,P_D))

lambda_vals <- seq(0, 0.1, by=0.001)
I2_vec <- numeric(length(lambda_vals))
I1_vec <- numeric(length(lambda_vals))

AgeDist <- c(0.01518000, 0.01518000, 0.01518000, 0.01518000, 0.01518000, 0.01560000, 0.01560000, 0.01560000, 0.01560000, 0.01560000,
             0.01478000, 0.01478000, 0.01478000, 0.01478000, 0.01478000, 0.01617775, 0.01617775, 0.01617775, 0.01617775, 0.01617775,
             0.01506456, 0.01506456, 0.01506456, 0.01506456, 0.01506456, 0.01525617, 0.01525617, 0.01525617, 0.01525617, 0.01525617,
             0.01610475, 0.01610475, 0.01610475, 0.01610475, 0.01610475, 0.01383275, 0.01383275, 0.01383275, 0.01383275, 0.01383275)

AgeDist <- c(AgeDist, rep(AgeDist[length(AgeDist)],20)) # fill age vector to age 60


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



# Estimate FOI in Gampaha -------------------------------------------------

# * I2 vs. lambda plot ------------------------------------------------------

# obtain prop of people with second infection in a year at equilibrium
for(i in 1:length(lambda_vals)){
  
  dengue_df <- infection_probs_3st(lambda = rep(lambda_vals[i], length(years)), 
                                   max_age = max_age, 
                                   P_D_mat = P_D_mat, 
                                   output.df = dengue,
                                   verbose = 0)
  
  I2_vec[i] <- sum(AgeDist[1:60]*dengue_df[length(years)-1,1:60,"I2"])/sum(AgeDist[1:60])
  I1_vec[i] <- sum(AgeDist[1:60]*dengue_df[length(years)-1,1:60,"I1"])/sum(AgeDist[1:60])
  
}

# plot proportion of people with second infection at equilibrium relative to prop of people
# with second infection in Colombo for a range of FOIs. this is used to find out the FOI
# in the region where prop of people with second infection is 77% of that in Colombo

#if using 4-serotype model, use I2_vec[38] and lambda_vals[38]
plot(x = I2_vec/I2_vec[50], y = lambda_vals, type = "l", 
     xlab = "Modelled cases proportional to Colombo" , 
     ylab = "Force of Infection")
abline(v = 1.0)
abline(h = lambda_vals[50])
abline(v = 0.77)
#value of lambda where I2_vec/I2_vec[38] is closest to 0.77
abline(h = lambda_vals[21])

# * pretty version of plot  -----------------------------------------------
# FOI values used for labels are still for 4 serotype model
font.size <- 2
lwd <- 4
par(mar = c(5.1, 5.1, 4.1, 2.1))

plot(x = I2_vec/I2_vec[50], y = lambda_vals, type = "l", lwd=lwd,
     xlab="Secondary infections relative to Colombo", ylab="Force of infection", cex.axis=font.size, cex.lab=font.size)

abline(v = 1.0, lwd=lwd, lty=2)
abline(h = lambda_vals[50], lwd=lwd, lty=2)
abline(v = 0.77, lwd=lwd, lty=2)
#value of lambda where I2_vec/I2_vec[38] is closest to 0.77
abline(h = lambda_vals[21], lwd=lwd, lty=2)

abline(v = 1.0, lwd=lwd, lty=2, col="blue")
abline(h = lambda_vals[50], lwd=lwd, lty=2, col="blue")
abline(v = 0.77, lwd=lwd, lty=2, col="green")
#value of lambda where I2_vec/I2_vec[38] is closest to 0.77
abline(h = lambda_vals[21], lwd=lwd, lty=2, col="green")
text(x=0.1, y=0.018, "FOI 0.014", col="green", cex=font.size)
text(x=0.1, y=0.041, "FOI 0.037", col="blue", cex=font.size)

# Calculating rho ---------------------------------------------------------

# nr cases = nr infections * detection prob
# detection prob (rho) = nr cases/ nr infections
# det prob = (100/10000) / (I2/1)
rho <- (100/10000)/I2_vec[50]

# when using both 1st and 2nd infections can lead to cases:
# nr cases = (nr 1st inf * detection prob 1) + (nr 2nd inf * det prop 2)
# should set a relationship between det prob 1 and 2 such that
# nr cases = (nr 1st inf * detection prob 2 * relative severity) + (nr 2nd inf * det prop 2)

# Calculate no of infections in Gampaha ----------------------------------------

# Only focus in I1 and I2

#scenario - only 3 serotypes circulating (before and during trial)

dengue_df.Gamp_3st <- infection_probs_3st(lambda = rep(lambda_vals[21], length(years)), 
                                          max_age = max_age, 
                                          P_D_mat = P_D_mat, 
                                          output.df = dengue,
                                          verbose = 0)

#proportion of population who experience primary or secondary infections per year
ageRange <- 4:16
inf_c_3st <- sum(AgeDist[ageRange]*
                   rowSums(dengue_df.Gamp_3st[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

# proportion of population who experience primary or secondary infections over a period of 2 years
inf_c_2yr_3st <- 1-(1-inf_c_3st)^2

# proportion of population who shows up as case over a period of 2 years
case_c_2yr_3st <- 1-(1-(77/10000))^2

# 4th serotype emerges, same FoI ----------
#re-initialize new dataframe for 4th serotype introduction
#will remove last age
max_age <- 61
#will remove last year
years <- 2010:2110
N0 <- 1 # population size
no_infection <- N0
#same per-serotype FOI
lambda <- rep(lambda_vals[21], length=length(years))
P_D <- c(0.005, 0.05, 0.00001, 0.00001)

P_D_mat <- t(replicate(max_age,P_D))

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
dengue_df.Gamp_3to4st <- array(data = 0,
                               dim = c(length(years), max_age+1,length(col_names)),
                               dimnames=array_names)

#when 4th serotype is introduced
#S3 --> S4 (susceptible to 3 become susceptible to 4)
#S2 --> S3 (susceptible to 2 become susceptible to 3)
#S1 --> S2 (susceptible to 1 become susceptible to 2)
#S0 --> S1 (susceptible to 0 become susceptible to 1)
#nobody is S0 (everyone susceptible to at least 1 serotype)

#not everyone is susceptible to all 4 ST at start of simulation
#initialized from 3-serotype simulation
#dengue_df.Gamp_yr2[1,,"S4"] <- N0 # add 100 to S4 to each age group in first year
dengue_df.Gamp_3to4st[,1,"S4"] <- N0 # add 100 to S4 of 0 age group in all years (births)

#dengue_df.Gamp_yr2[1,,"total"] <- dengue_df.Gamp_yr2[1,,"S4"] # add total to first year
#dengue_df.Gamp_yr2[,1,"total"] <- dengue_df.Gamp_yr2[,1,"S4"] # add total to first age

#fill in matrix with data from dengue_df.Gamp
#susceptibles
dengue_df.Gamp_3to4st[1,-(max_age+1),"S4"] <- dengue_df.Gamp_3st[length(years)-1,,"S3"]
dengue_df.Gamp_3to4st[1,-(max_age+1),"S3"] <- dengue_df.Gamp_3st[length(years)-1,,"S2"]
dengue_df.Gamp_3to4st[1,-(max_age+1),"S2"] <- dengue_df.Gamp_3st[length(years)-1,,"S1"]
dengue_df.Gamp_3to4st[1,-(max_age+1),"S1"] <- dengue_df.Gamp_3st[length(years)-1,,"S0"]

dengue_df.Gamp_3to4st[1,1:max_age+1,"total"] <- rowSums(dengue_df.Gamp_3to4st[1,1:max_age+1,which(col_names %like% "S")])

#infections
#dengue_df.Gamp_yr2[1,-(max_age+1),"I1"] <- dengue_df.Gamp[length(years)-1,,"I1"]
#dengue_df.Gamp_yr2[1,-(max_age+1),"I2"] <- dengue_df.Gamp[length(years)-1,,"I2"]
#dengue_df.Gamp_yr2[1,-(max_age+1),"I3"] <- dengue_df.Gamp[length(years)-1,,"I3"]
#disease
#dengue_df.Gamp_yr2[1,-(max_age+1),"D1"] <- dengue_df.Gamp[length(years)-1,,"D1"]
#dengue_df.Gamp_yr2[1,-(max_age+1),"D2"] <- dengue_df.Gamp[length(years)-1,,"D2"]
#dengue_df.Gamp_yr2[1,-(max_age+1),"D3"] <- dengue_df.Gamp[length(years)-1,,"D3"]

dengue_df.Gamp_3to4st <- infection_probs_4st(lambda = rep(lambda_vals[21], length(years)), 
                                             max_age = max_age, 
                                             P_D_mat = P_D_mat, 
                                             output.df = dengue_df.Gamp_3to4st,
                                             verbose = 0)

#did I do this right?
#summing equilibrium infected under 3 serotypes (assuming 3 serotypes circulating) + 
#equilibrium under 4 serotypes (initialized with proportion susceptible after 3 serotypes circulating,
#per-serotype FoI stays the same)
ageRange <- 4:16
inf_c_yr1_3st <- sum(AgeDist[ageRange]*
                       rowSums(dengue_df.Gamp_3st[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_yr2_4st <- sum(AgeDist[ageRange]*
                       rowSums(dengue_df.Gamp_3to4st[1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_2yr_3to4st <- inf_c_yr1_3st + inf_c_yr2_4st

#when 4 serotypes are co-circulating ------
###4 st scenario
# probability first infection = 14.1%
# 0.141 <- 1-(1 - lambda)^4
# lambda = 0.037 in Colombo
# in Gampaha - lambda = 0.014
#will remove last age
max_age <- 61
#will remove last year
years <- 2010:2110
N0 <- 1 # population size
no_infection <- N0
lambda <- rep(lambda_vals[15], length=length(years))
P_D <- c(0.005, 0.05, 0.00001, 0.00001)

P_D_mat <- t(replicate(max_age,P_D))

# Initialise dataframe

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

dengue_df.Gamp_2yr_4st <- infection_probs_4st(lambda = rep(lambda_vals[15], length(years)), 
                                              max_age = max_age, 
                                              P_D_mat = P_D_mat, 
                                              output.df = dengue,
                                              verbose = 0)

#proportion of population who experience primary or secondary infections per year
ageRange <- 4:16
inf_c_4st <- sum(AgeDist[ageRange]*
                   rowSums(dengue_df.Gamp_2yr_4st[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

# proportion of population who experience primary or secondary infections over a period of 2 years
inf_c_2yr_4st <- 1-(1-inf_c_4st)^2

# calculate number of cases based on secondary infections --------------
#assuming same rho for scenario where 4th serotype emerges
#3 serotype scenario cases - calc above

#3 to 4 serotype scenario cases
case_c_3to4st_y1 <- rho*sum(AgeDist[1:60]*dengue_df.Gamp_3st[length(years)-1,1:60,"I2"])/sum(AgeDist[1:60])
case_c_3to4st_y2 <- rho*sum(AgeDist[1:60]*dengue_df.Gamp_3to4st[1,1:60,"I2"])/sum(AgeDist[1:60])

case_c_2yr_3to4st <- case_c_3to4st_y1 + case_c_3to4st_y2

#4 serotype scenario cases
case_c_2yr_4st <- rho*sum(AgeDist[1:60]*dengue_df.Gamp_2yr_4st[length(years)-1,1:60,"I2"])/sum(AgeDist[1:60])
#not for 2 years ^^^




