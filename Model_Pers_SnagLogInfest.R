#library(R2WinBUGS)
library(R2jags) #or...
library(abind)
library(R.utils)
setwd("~/Documents/Rdata/LoggingData_Julia")
load("Data_compiled.RData")

# Covariate inputs #
Logging <- (Plot.data$Logging - mean(Plot.data$Logging))/sd(Plot.data$Logging)
snag <- (Plot.data$snag - mean(Plot.data$snag))/sd(Plot.data$snag)

Infest <- as.matrix(Plot.data[,c("EarlInf_2014","EarlInf_2015","EarlInf_2016")])
Infest <- (Infest - mean(Infest))/sd(Infest)
Infest <- cbind(Infest[,1],Infest)

nplot <- length(plot)
nyear <- length(year)
nspec <- length(spp)
Y <- Y.arry

Z.init <- apply(Y,c(1,2,3),function(x) sum(x>=1))

# Assemble the data names list for JAGS.
data <- list("Y","nplot","nyear","n.visits","nspec","Logging","snag","Infest")

# Assemble the initial values for JAGS.  Here z is the actual detection 

inits <- function(){list(Z=Z.init,rho=runif(1,-1,1))}

# Assemble the parameters vector for WinBUGS (What we want to track).
#For model parameters - this one without the PC's
parameters <- c("pspp","beta0","beta1","beta.Logging","beta.snag","beta.Infest",
                "rho","SR","test")

sink("model.txt")
cat("

# The model.
model{

for(i in 1:nspec){
  #Species-specific parameters for detection
  # Species-specific parameters for baseline occupancy (intercept) and related parameters for each location.
  beta0[i] ~ dnorm(mu.beta0,tau.beta0)T(-10,10) # mean probability of occupancy for sites previously unoccupied (= colonization or gamma)
  mu.eta[i] <- mu.alpha0 + (rho*sd.alpha0/sd.beta0)*(beta0[i] - mu.beta0)      #sd.alpha0 = sd(detection), sd.beta0 = sd(occupancy)
  eta[i] ~ dnorm(mu.eta[i], var.eta)T(-10,10) #incorporates correlation within hyperprior mean and variance for detection
  logit(pspp[i]) <- eta[i]      	#detection at intercept with correlation structure

  #Species-specific parameters for occupancy covariates in the northwest
  beta.Logging[i] ~ dnorm(mu.beta.Logging, tau.beta.Logging)T(-10,10)
  beta.snag[i] ~ dnorm(mu.beta.snag, tau.beta.snag)T(-10,10)
  beta.Infest[i] ~ dnorm(mu.beta.Infest, tau.beta.Infest)T(-10,10)

  beta1[i] ~ dnorm(mu.beta1, tau.beta1)T(-10,10) # Persistence offset (added to beta0 if occupied in previous year) 
  # colonization (or gamma) = expit(beta0)
  # persistence (or phi) = expit(beta0 + beta1)

  #Model
  psi0[i] ~ dunif(0,1)    
  for(j in 1:nplot){   
    z0[j,i] ~ dbern(psi0[i]) # Occupancy state in year prior to sampling
    ###Model for occupancy probability ###
    logit(psi[j,i,1]) <- beta0[i] + beta1[i]*z0[j,i] + beta.Logging[i]*Logging[j] + beta.snag[i]*snag[j] +
      beta.Infest[i]*Infest[j,1]
    ###Model for detection probability ###
    Z[j,i,1] ~ dbern(psi[j,i,1])
    pZ[j,i,1] <- pspp[i]*Z[j,i,1]
    Y[j,i,1] ~ dbin(pZ[j,i,1],n.visits[1])
    #____________Bayesian GOF_________________________________
    ynew[j,i,1] ~ dbin(pZ[j,i,1],n.visits[1])  #simulated data y under model
    
    LLsim[j,i,1] <- (ynew[j,i,1]*log(pspp[i])+
      (n.visits[1]-ynew[j,i,1])*log(1-pspp[i]))*Z[j,i,1]  #log-likelihood simulated data
    LLdata[j,i,1]<- (Y[j,i,1]*log(pspp[i])+
      (n.visits[1]-Y[j,i,1])*log(1-pspp[i]))*Z[j,i,1]    #log-likelihood observed data
    #_________________________________________________________
    for(t in 2:nyear){
      ###Model for occupancy probability ###
      logit(psi[j,i,t]) <- beta0[i] + beta1[i]*Z[j,i,(t-1)] + beta.Logging[i]*Logging[j] + beta.snag[i]*snag[j] +
        beta.Infest[i]*Infest[j,t]
      ###Model for detection probability ###
      Z[j,i,t] ~ dbern(psi[j,i,t]) # Occupancy state
      pZ[j,i,t] <- pspp[i]*Z[j,i,t] # Unconditional probability of detection
      Y[j,i,t] ~ dbin(pZ[j,i,t],n.visits[t]) # Data model
      #____________Bayesian GOF_________________________________
      ynew[j,i,t] ~ dbin(pZ[j,i,t],n.visits[t])  #simulated data y under model
    
      LLsim[j,i,t] <- (ynew[j,i,t]*log(pspp[i])+
        (n.visits[t]-ynew[j,i,t])*log(1-pspp[i]))*Z[j,i,t]  #log-likelihood simulated data
      LLdata[j,i,t]<- (Y[j,i,t]*log(pspp[i])+
        (n.visits[t]-Y[j,i,t])*log(1-pspp[i]))*Z[j,i,t]    #log-likelihood observed data
        #_________________________________________________________
        }
      }
    }

for(j in 1:nplot){
  for(t in 1:nyear){
    SR[j,t] <- sum(Z[j,,t])
  }
}

### Hyper-parameters for detection model ###

mu.alpha0 ~ dnorm(0,0.1)
tau.alpha0 <- 1/(sd.alpha0*sd.alpha0)
sd.alpha0 ~ dunif(0,10)

### Hyper-parameters for occupancy model ###

# Hyper-parameters for occupancy intercept and persistence parameters at each location #

mu.beta0 ~ dnorm(0,0.1)         #Hyper-parameter for intercept
sd.beta0 ~ dunif(0,10)
tau.beta0 <- 1/(sd.beta0*sd.beta0)

mu.beta1 ~ dnorm(0,0.1)         #Hyper-parameter for intercept
sd.beta1 ~ dunif(0,10)
tau.beta1 <- 1/(sd.beta1*sd.beta1)

mu.beta.Logging ~ dnorm(0,0.1)
tau.beta.Logging <- 1/(sd.beta.Logging*sd.beta.Logging)
sd.beta.Logging ~ dunif(0,10)

mu.beta.snag ~ dnorm(0,0.1)
tau.beta.snag <- 1/(sd.beta.snag*sd.beta.snag)
sd.beta.snag ~ dunif(0,10)

mu.beta.Infest ~ dnorm(0,0.1)
tau.beta.Infest <- 1/(sd.beta.Infest*sd.beta.Infest)
sd.beta.Infest ~ dunif(0,10)

rho ~ dunif(-1,1)				#correlation parameter
var.eta <- tau.alpha0/(1.-pow(rho,2))  #variance in detection

#________Bayesian GOF___________
#deviance
dev_sim <- (-2)*sum(LLsim[,,])
dev_data <- (-2)*sum(LLdata[,,])

#test statistics should be ~0.5 if model fits
test<-step(dev_data-dev_sim)
#________________________________

}
",fill=TRUE)
sink()

# MCMC values.
nc <- 4
nb <- 5000
ni <- 70000
nt <- 1

# Send it all to JAGS and hope for the best!

# To help track time.
starttime <- Sys.time()

bugout <- jags(data, inits, parameters, model.file="model.txt", n.chains=nc, n.iter=ni,
n.burnin=nb, n.thin=nt)
#run this if it didn't converge and you don't want the model to run the burn in this time
#ni <- 50000
#bugout <- update(bugout, n.iter=ni)

endtime <- Sys.time()
runtime <- endtime - starttime
runtime

#Check n.effectives and R.hats for parameters
length(which(bugout$BUGSoutput$summary[,"n.eff"]<100))/length(bugout$BUGSoutput$summary[,"n.eff"])
min(bugout$BUGSoutput$summary[,"n.eff"])
sort(bugout$BUGSoutput$summary[,"n.eff"])[1:50]
max(bugout$BUGSoutput$summary[,"Rhat"])
sort(bugout$BUGSoutput$summary[,"Rhat"],decreasing=T)[1:50]

#use bugout object to manipulate results in R environment

#save results
library(R.utils)
saveObject(bugout,"Mod_Pers_SnagLogInfest")
