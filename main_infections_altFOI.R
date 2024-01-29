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
  
  dengue_df <- infection_probs(lambda = rep(lambda_vals[i], length(years)), 
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

#default FoI
dengue_df.Gamp <- infection_probs(lambda = rep(lambda_vals[21], length(years)), 
                                  max_age = max_age, 
                                  P_D_mat = P_D_mat, 
                                  output.df = dengue,
                                  verbose = 0)

dengue_df.Gamp_ll <- infection_probs(lambda = c(rep(lambda_vals[21], length(years)-3),
                                                rep(lambda_vals[21]*0.5, 3)), 
                                  max_age = max_age, 
                                  P_D_mat = P_D_mat, 
                                  output.df = dengue,
                                  verbose = 0)

dengue_df.Gamp_lh <- infection_probs(lambda = c(rep(lambda_vals[21], length(years)-3),
                                                lambda_vals[21]*0.5,
                                                lambda_vals[21]*2,
                                                lambda_vals[21]), 
                                     max_age = max_age, 
                                     P_D_mat = P_D_mat, 
                                     output.df = dengue,
                                     verbose = 0)

dengue_df.Gamp_hd <- infection_probs(lambda = c(rep(lambda_vals[21], length(years)-3),
                                                lambda_vals[21]*2,
                                                lambda_vals[21],
                                                lambda_vals[21]), 
                                     max_age = max_age, 
                                     P_D_mat = P_D_mat, 
                                     output.df = dengue,
                                     verbose = 0)

dengue_df.Gamp_ld <- infection_probs(lambda = c(rep(lambda_vals[21], length(years)-3),
                                                lambda_vals[21]*0.5,
                                                lambda_vals[21],
                                                lambda_vals[21]), 
                                     max_age = max_age, 
                                     P_D_mat = P_D_mat, 
                                     output.df = dengue,
                                     verbose = 0)

#proportion of population who experience primary or secondary infections per year
ageRange <- 4:16

#calculate proportion of population who experience primary or secondary infection
#for each FOI scenario
inf_c_d <- sum(AgeDist[ageRange]*
                 rowSums(dengue_df.Gamp[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_ll_y1 <- sum(AgeDist[ageRange]*
                     rowSums(dengue_df.Gamp_ll[length(years)-2,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])
inf_c_ll_y2 <- sum(AgeDist[ageRange]*
                     rowSums(dengue_df.Gamp_ll[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_lh_y1 <- sum(AgeDist[ageRange]*
                     rowSums(dengue_df.Gamp_lh[length(years)-2,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])
inf_c_lh_y2 <- sum(AgeDist[ageRange]*
                     rowSums(dengue_df.Gamp_lh[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_hd_y1 <- sum(AgeDist[ageRange]*
                     rowSums(dengue_df.Gamp_hd[length(years)-2,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])
inf_c_hd_y2 <- sum(AgeDist[ageRange]*
                     rowSums(dengue_df.Gamp_hd[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

inf_c_ld_y1 <- sum(AgeDist[ageRange]*
                  rowSums(dengue_df.Gamp_ld[length(years)-2,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])
inf_c_ld_y2 <- sum(AgeDist[ageRange]*
                  rowSums(dengue_df.Gamp_ld[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

# proportion of population who experience primary or secondary infections over a period of 2 years
inf_c_2yr_d <- 1-(1-inf_c_d)^2
inf_c_2yr_ll <- inf_c_ll_y1 + inf_c_ll_y2
inf_c_2yr_lh <- inf_c_lh_y1 + inf_c_lh_y2
inf_c_2yr_hd <- inf_c_hd_y1 + inf_c_hd_y2
inf_c_2yr_ld <- inf_c_ld_y1 + inf_c_ld_y2

# proportion of population who shows up as case over a period of 2 years
#first calc secondary infections for each FoI scenario
inf_c_ll_sec_y1 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_ll[length(years)-2,ageRange,c("I2")]))/sum(AgeDist[ageRange])
inf_c_ll_sec_y2 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_ll[length(years)-1,ageRange,c("I2")]))/sum(AgeDist[ageRange])

inf_c_lh_sec_y1 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_lh[length(years)-2,ageRange,c("I2")]))/sum(AgeDist[ageRange])
inf_c_lh_sec_y2 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_lh[length(years)-1,ageRange,c("I2")]))/sum(AgeDist[ageRange])

inf_c_hd_sec_y1 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_hd[length(years)-2,ageRange,c("I2")]))/sum(AgeDist[ageRange])
inf_c_hd_sec_y2 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_hd[length(years)-1,ageRange,c("I2")]))/sum(AgeDist[ageRange])

inf_c_ld_sec_y1 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_ld[length(years)-2,ageRange,c("I2")]))/sum(AgeDist[ageRange])
inf_c_ld_sec_y2 <- sum(AgeDist[ageRange]*
                         (dengue_df.Gamp_ld[length(years)-1,ageRange,c("I2")]))/sum(AgeDist[ageRange])

#then calc cases per year based on rho
#for default scenario - just use data
case_c_2yr_d <- 1-(1-(77/10000))^2
case_c_2yr_ll <- rho*(inf_c_ll_sec_y1 + inf_c_ll_sec_y2)
case_c_2yr_lh <- rho*(inf_c_lh_sec_y1 + inf_c_lh_sec_y2)
case_c_2yr_hd <- rho*(inf_c_hd_sec_y1 + inf_c_hd_sec_y2)
case_c_2yr_ld <- rho*(inf_c_ld_sec_y1 + inf_c_ld_sec_y2)







