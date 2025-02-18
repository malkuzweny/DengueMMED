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

dengue_df.Gamp <- infection_probs(lambda = rep(lambda_vals[21], length(years)), 
                                  max_age = max_age, 
                                  P_D_mat = P_D_mat, 
                                  output.df = dengue,
                                  verbose = 0)

#proportion of population who experience primary or secondary infections per year
ageRange <- 4:16
inf_c <- sum(AgeDist[ageRange]*
               rowSums(dengue_df.Gamp[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])

# proportion of population who experience primary or secondary infections over a period of 2 years
inf_c_2yr <- 1-(1-inf_c)^2

# proportion of population who shows up as case over a period of 2 years
case_c_2yr <- 1-(1-(77/10000))^2


# * for different age groups ----------------------------------------------

# create dataframe of proportion of population experiencing 1st or 2nd infection for multiple age groups
# age groups: varying minimum age to age 60
prop.infection <- data.frame(prop.infection=as.numeric(), minAge=as.numeric(), population=as.numeric())
minAge <- 1
maxAge <- 60

while (minAge < 40){
  ageRange <- minAge:maxAge
  
  # proportion of people experiencing infection
  inf_c <- sum(AgeDist[ageRange]*
                 rowSums(dengue_df.Gamp[length(years)-1,ageRange,c("I1","I2")]))/sum(AgeDist[ageRange])
  inf_c_2yr <- 1-(1-inf_c)^2
  
  # number of people in this age group (total population size 2.5m)
  population <- 2500000 * sum(AgeDist[minAge:60])  
  
  output <- c(inf_c_2yr, minAge, population)
  output <- t(as.data.frame(output))
  prop.infection <- rbind(prop.infection, output)
  minAge <- minAge+1
}
names(prop.infection) <- c("prop.infection", "minAge", "population")
prop.infection$maxAge <- maxAge

plot(prop.infection$minAge, prop.infection$prop.infection)

# * plot % infected by age at equilibrium ---------------------------------------------
# use all infections (I1 to I3 or I4)
# needs adapting to 3 serotype model

plotDF <- dengue_df.Gamp
par(mar = c(5.1, 5.1, 4.1, 2.1))

font.size <- 1.3
lwd <- 4

# par(mar = c(bottom, left, top, right)) , alue for mar is c(5.1, 4.1, 4.1, 2.1)  


{
  cols <- brewer.pal(4, "Set1")
  
  plot( x=1:(max_age), y = plotDF[100,,"I1" ], type = "l", col=cols[1], lwd=lwd,
        ylim = c(0, 0.06), xlab="Age", ylab="% of population infected", cex.axis=font.size, cex.lab=font.size)
  
  lines( x = (1:max_age), y = plotDF[100,,"I2" ], col=cols[2], lwd=lwd)
  lines( x = (1:max_age), y = plotDF[100,,"I3" ], col=cols[3], lwd=lwd)
  # lines( x = (1:max_age), y = plotDF[100,,"I4" ], col=cols[4], lwd=lwd)
  abline(v=4)
  abline(v=16)
  rect(xleft = 4, xright = 16, ybottom = par("usr")[3], ytop = par("usr")[4],
       border = NA, col = adjustcolor("black", alpha = 0.3))
  legend(x = c(20, 55), y = c(0.03, 0.06), legend=c("First infection", "Second infection", "Third infection"),
         col=c(cols[1], cols[2], cols[3]), lty=1, cex=font.size, bg="white", bty="n", lwd=lwd)
  
}

# Adjusted 
png("Figure_1_update.png", width = 3000, height = 2000, res = 300  )  # Increase width to 1200 pixels

{
  # Your plotting code here:
     #cols <- colorRampPalette(c("lightcoral", "darkred"))(3)
     cols <- colorRampPalette(c("mistyrose", "darkred"))(3)
  
  plot(x = 1:max_age, y = (plotDF[100, , "I1"]*100), type = "l", col = cols[1], lwd = lwd,
       ylim = c(0, 0.06*100), xlab = "Age", ylab = "% of population infected",
       cex.axis = font.size, cex.lab = font.size)
  
  lines(x = 1:max_age, y = plotDF[100, , "I2"]*100, col = cols[2], lwd = lwd)
  lines(x = 1:max_age, y = plotDF[100, , "I3"]*100, col = cols[3], lwd = lwd)
  
  abline(v = 4)
  abline(v = 16)
  
  rect(xleft = 4, xright = 16, ybottom = par("usr")[3], ytop = par("usr")[4],
       border = NA, col = adjustcolor("black", alpha = 0.3))
  
  legend("topright", 
         legend = c("First infection", "Second infection", "Third infection"),
         col = cols, lty = 1, lwd = lwd, cex = font.size, bg = "white", bty = "o")
}
dev.off()






