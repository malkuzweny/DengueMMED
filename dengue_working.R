# we want to recreate the Dengue model. 

# we frist create the four probabiilties of infection 
getwd()
library(tidyverse)
library(ggplot)

#============ infection P

lambda <- 0.01
max_age <- 100

N0 <- 100 

no_infection <- N0

col_length <- rep(0, max_age +1)

dengue <- tibble( 
  age = seq( 0, max_age, by = 1), 
  no_infected = c(100, rep(0, max_age )),
  first_inf = col_length,
  sec_inf = col_length, 
  third_inf = col_length,
  fourth_inf = col_length,
  new_first = col_length, 
  new_second = col_length, 
  new_third = col_length, 
  new_fourth = col_length)



infect1_0
infection_probs <- function( lambda, max_age, df, verbose=1 ){
  if(verbose==1){browser()}
  
  # infection probabilities
  infect1_0 <- 1-(1- lambda)^4
  infect2_1 <- 1-(1- lambda)^3
  infect3_2 <- 1-(1- lambda)^2
  infect4_3 <- 1-(1- lambda)^1
  infect5_4 <- 0L
  
  
  # first column(Numbe of people who are susceuptible)
  for( ii in 2:(max_age+1)) {
    
    df[ii, 2] <- df[ii -1, 2] * (1-infect1_0) # probability 
    
  }
  
  # second column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 3] <- df[ii -1, 2] * (infect1_0) + # prob of people with first infection 
      (df[ii-1, 3] * (1-infect2_1) ) # probability of not gettign a second infection
    
    
    df[ii, 7] <- df[ii -1, 2] * (infect1_0)  # Number of new cases

    
  }
  
  # third column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 4] <- df[ii -1, 3] * (infect2_1) + # prob of people with second infection 
      (df[ii-1, 4] * (1-infect3_2) ) # probability of not gettign a third infection
    
    df[ii, 8] <- df[ii -1, 3] * (infect2_1)   # Number of new cases
  }
  
  
  # fourth column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 5] <- df[ii -1, 4] * (infect3_2) + # prob of people with third infection 
      (df[ii-1, 5] * (1-infect4_3) ) # probability of not gettign a fourth infection
    
    df[ii, 9] <- df[ii -1, 4] * (infect3_2) # number of new cases 
  }
  
  
  # fifth column 
  for( ii in 2:(max_age+1)) {
    
    df[ii, 6] <- df[ii -1, 5] * (infect4_3) + # prob of people with third infection 
      (df[ii-1, 6] * (1-infect5_4) ) # probability of not gettign a fourth infection
    
    df[ii, 10] <- df[ii -1, 5] * (infect4_3) # number of new cases 
  }
  
  df<- df %>% mutate( total = rowSums(df[, 2:6])) 
  
  return(df)
}


dengue_lambda_Lo <- infection_probs( lambda = 0.01, max_age = max_age, df = dengue, verbose = 0 )


dengue_lambda_Hi <- infection_probs( lambda = 0.05, max_age = max_age, df = dengue, verbose = 0 )

sum(dengue_lambda_Lo$new_first)
sum(dengue_lambda_Lo$new_second)
sum(dengue_lambda_Lo$new_third)
sum(dengue_lambda_Lo$new_fourth)

infect1_0 <- (1-(1- lambda)^400) * (1-(1- lambda)^400)

dengue_lambda_Lo

dengue_lambda_Lo%>%view()

 plot( x= dengue_lambda_Lo$age, y = dengue_lambda_Lo$new_second, type = "l", ylim = c(0, 5))
 lines( x = dengue_lambda_Hi$age, y = dengue_lambda_Hi$new_second)

 

