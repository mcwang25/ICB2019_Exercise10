# Code to model tumor growth with drug-resistent mutant subpopulation
# Tumor starts as a single cell (N=1)
# At 100 cells, one cell becomes drug-resistent (M=1)
# After tumor reaches equilibrium, drug is introduced, only killing N-cells and
# reducing growth rate of M-cells
library(ggplot2)
library(cowplot)
library(reshape2)

K <- 1000000
rN <- 0.1
rM <- 0.1
N0 <- 1
M0 <- 0

timesteps <- 500
count <- 0
Npop <- numeric(length=1000)
Npop[1] <- N0
Mpop <- numeric(length=1000)
Mpop[1] <- M0

for (t in 1:timesteps){
  if (count<100) {
    Mpop[t+1] <- 0
    Npop[t+1] <- (Npop[t]+rN*Npop[t]*(1-(Npop[t]+Mpop[t])/K))
    count <- Npop[t+1]
    timer <- t
  }else{
    Mpop[timer+1] <- 1
    Npop[t+1] <- (Npop[t]+rN*Npop[t]*(1-(Npop[t]+Mpop[t])/K))
    Mpop[t+1] <- (Mpop[t]+rM*Mpop[t]*(1-(Npop[t]+Mpop[t])/K))
  }
}

# Drug is introduced at t=501
# The daily rate at which N-cells decrease after drug treatment was changed from
# -0.1 to -1.0 because (Npop[t]+Mpop[t])/K was too small to show significant
# decline, even after 1000000 time units.
rMd <- 0.05
rNd <- 1.0
for (t in 501:999){
  Npop[t+1] <- (Npop[t]-rNd*Npop[t]*(1-(Npop[t]+Mpop[t])/K))
  Mpop[t+1] <- (Mpop[t]+rMd*Mpop[t]*(1-(Npop[t]+Mpop[t])/K))
}

timeaxis=as.numeric(1:1000)

plot(timeaxis, Npop, type="l", xlab="time", ylab="population", main="Tumor Cell Populations vs. Time")+lines(timeaxis,Mpop, type="l", col="blue")+legend(175,4e+05,legend=c("Susceptible cells", "Resistent cells"), col=c("black","blue"), lty=1:1, cex=0.8)
