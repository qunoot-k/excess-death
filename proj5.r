# Qunoot Khaleeq (s2314595)

# In this project I will implement a demographic model for England and Wales, 
# to predict the expected number of deaths per week from the beginning of 2020 
# assuming death rates had stayed at the levels seen in 2017-19.
# Data for annual probability of death for each one year age group, 
# and the populations of each age group at the start of a period
# in UK are available from the Office for National Statistics (ONS)
# The model will include an adjustment to allow for seasonal variation 
# in mortality rates. I will then compare this to the actual deaths per week 
# that occurred over this period, to obtain an excess deaths time series, 
# and then model this series using a simple Bayesian model in JAGS.

library(rjags)

model <- function(Ni, mi, dj, Ni_1star) {
  # computes death, reduces the population, and applies aging
  
  # Inputs
  # Ni - current population of an age group
  # mi - annual death rate for an age band
  # dj - mortality rate modifier for a week
  # Ni_1star - aging parameter
  
  # returns list of deaths, updated population, and aging parameter

  qi <- 1 - exp(-mi/52)
  Di <- 0.9885*dj*qi*Ni
  Ni_star <- Ni - Di
  Ni_plus <- (Ni_star*51/52) + (Ni_1star/52)
  return(list(Ni_star=Ni_star, Ni_plus=Ni_plus, Di=Di))
}

predictDeath <- function(pop_male, pop_female, mf, mm, dj) {
  # predicts deaths for every age class over certain weeks
  
  # Inputs
  # pop_male - male population at the start of a period in each 1 year age class
  # pop_female - female population at the start of a period in each 1 year age class
  # mf <- female annual death rates for each 1-year age band in any year
  # mm <- male annual death rates for each 1-year age band in any year
  # dj <-  mortality rate modifier for weeks to be predicted
  
  # return the predicted total number of deaths each week as a vector 
  
  weekly_deaths <- c()
  Nf_0star <- pop_female[1] # birth constant for female new born
  Nm_0star <- pop_male[1] # birth constant for male new born
  deaths <- c()
  
  # for each week calculates total death in every age class
  # by running the model separately for both males and females
  for (j in 1:length(dj)){
    Nfi_1star <- Nf_0star # growth rate for age class 0 each week is constant
    Nmi_1star <- Nm_0star
    for (i in 1:length(pop_female)) {
      modf <- model(pop_female[i], mf[i], dj[j], Nfi_1star)
      modm <- model(pop_male[i], mm[i], dj[j], Nmi_1star)
      Nfi_1star <- modf$Ni_star
      Nmi_1star <- modm$Ni_star
      pop_female[i] <- modf$Ni_plus
      pop_male[i] <- modm$Ni_plus
      deaths[i] <- modf$Di + modm$Di # total weekly death for each age class
    }
    weekly_deaths[j] <- sum(deaths) # total weekly deaths
  }
  return(weekly_deaths)
}

deathData <- read.table("death1722uk.dat", header=TRUE)
population <- read.table("lt1720uk.dat", header=TRUE)

# data for 2020 onward starts from week 157
d <- deathData$d[157:nrow(deathData)]
# time over which data is to be predicted
weeks <- length(d) 

# predict deaths for 2020 onward each week
predicted <- predictDeath(population$mpop20, population$fpop20,
                          population$mf, population$mm, d)

# actual deaths from 2020 onward for each week
observedDeath <- deathData$deaths[157:nrow(deathData)]

# excess death 
xi <- observedDeath - predicted

# total excess death in 2020 and overall
totalExcessDeaths <- sum(xi)
totalExcessDeaths2020 <- sum(xi[1:52])

# plot observed and expected deaths agains weeks
plot(observedDeath, type="l", col="red",
     ylim=c(0,max(observedDeath)), xlab="week", ylab="",
     main=paste("excess death 2020 ",round(totalExcessDeaths2020,digits=3), 
       " overall ",round(totalExcessDeaths, digits=3)))
lines(predicted, col="blue")
legend("topright", c("observed","predicted"), fill=c("red","blue"))

# plot of the cumulative excess deaths by week
plot(cumsum(xi), xlab="week", type="l",
     ylab=" ",main="cumulative excess deaths")

# set Christmas/New Year weeks to NA due to delay in reported 
# deaths which can cause various recording problems
index_na <- c(51, 52, 53, 105, 106)
temp <- xi
xi[index_na] = NA

# compile jags model
mod <- jags.model("model.jags",
                  data=list(x=xi, N=weeks))

# implement the model and use it to draw 10000 samples
# from the posterior densities of the mu vector, rho and k 
sam.coda <- coda.samples(mod,c("mu","rho","k"),n.iter=10000)
rho <- sam.coda[[1]][,"rho"][1:10000]
muh <- unname(as.data.frame(sam.coda[[1]][,2:(weeks+1)]))
xi <- temp

# trace plots and histograms for rho
plot(rho, type="l", xlab="iterations", ylab="", main="Trace of rho")
hist(rho, xlab="tau", main="Histogram for rho", probability=TRUE)

# posterior expected value vector for mu
E_mu <- colMeans(muh)

# plot showing every 50th sampled mu vector (against week) in grey, 
# with the estimated expectation for mu overlaid in blue,
# and observed excess deaths as black points (red if set to NA in JAGS call), 
plot(as.numeric(muh[50,]), type="l", col="grey", ylab=""
     , xlab="week",main="model comparison")
fifty <- seq(from = 100, to = 10000, by = 50)
sampled_lines <- sapply(fifty, function(x){
  lines(as.numeric(muh[x,]), col="grey")
})
lines(E_mu, col="blue")
points((1:weeks)[-index_na],xi[-index_na], col="black")
points((1:weeks)[index_na],xi[index_na], col="red")
legend("topright", c("sampled mu","estimated mu", "excess deaths", "NA")
       ,pch=c("-","-","o","o"), col=c("grey","blue", "black", "red"))

# plot for residuals from the model against time
plot(xi-E_mu, ylab="",type="l",xlab="week", main="residuals")