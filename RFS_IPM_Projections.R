##################################################################################################################
# Modeling Reticulated Flatwoods Salamander populations: Assessing effects of future climate change scenarios
# on long-term population viability. This script uses the functions and parameters in the other script to
# run the IPM forward and calculate population viability and other relevant metrics.
##################################################################################################################

#-----------------------------------------------------------------------------------------------------------------
# 1.0 Setup Workspace
#-----------------------------------------------------------------------------------------------------------------

# Load Packages
library(tidyverse)
library(lubridate)
library(patchwork)


# functions
std <- function(x){
  (x-mean(x)/sd(x))
}

# Function to view a matrix (Ellner et al.)
matrix.image <- function(A, x=NULL, y=NULL, col=rainbow(100,start=0.67,end=0), bw=FALSE, do.contour=FALSE, do.legend=TRUE,...) {
  if(do.legend) layout(mat=cbind(matrix(1,5,5),rep(2,5)));
  par(mar=c(6,5,3,2)); 
  if(is.null(x)) x=1:ncol(A);
  if(is.null(y)) y=1:nrow(A); 
  nx=length(x); ny=length(y); 
  x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
  y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
  if(bw) col=grey( (200:50)/200 ); 
  image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,cex.axis=1.5,cex.lab=1.5,bty="u",...);
  abline(v=range(x1)); abline(h=range(y1)); 
  if(do.contour) contour(x,y,t(A),nlevels=5,labcex=1.2,add=TRUE);   
  
  if(do.legend) {
    l.y=seq(min(A),max(A),length=100);  
    par(mar=c(6,2,3,1))
    image(list(x=1:2,y=l.y,z=rbind(l.y,l.y)),col=col,bty="o",xaxt="n",yaxt="n"); 
    axis(side=2,cex.axis=1.5,at=pretty(seq(min(A),max(A),length=10))); 
  } 
}

# Outline
# 1.0 Setup Workspace - run every time!
# 2.0 Prepare to run IPMs
# 3.0 
# 4.0 
# 5.0
# 6.0
# 7.0
# 8.0
# 9.0 Sensitivity and elasticity analyses


#-----------------------------------------------------------------------------------------------------------------
# 2.0 Prepare to run population simulations 
#-----------------------------------------------------------------------------------------------------------------

# 2.2 Load and format potential starting densities ---------------------------------------------------------------

Starting.Densities <- readRDS("Data/Starting_Densities.rds")

# Smoothed
Starting.Densities.SM <- split(Starting.Densities$Density.Sm, f = Starting.Densities$Dummy)

# Measured
Starting.Densities.ME <- split(Starting.Densities$Density.Act, f = Starting.Densities$Dummy)


# 2.3 Define integration parameters ------------------------------------------------------------------------------

# Define integration parameters
L <- 25 # lower size limit
U <- 80 # upper size limit
m <- 110 # Size of the iteration matrix (number of bins)
h <- (U - L) / m # bin width
meshpoints <- L + (1:m) * h - h/2 # meshpoitns

# Integration Parameters
i.par <- list(meshpoints = meshpoints, m = m, h = h)


# 2.4 Function to calculate lists of iteration matrices ----------------------------------------------------------

# Fully stochastic IPM
FUN.IPM.stoch <- function(i.par, m.par){
  within(i.par, {
    Rec.short <- list()
    Rec.long <- list()
    Rec.short.flood <- list()
    Rec.long.flood <- list()
    Surv <- list()
    K.short <- list()
    K.long <- list()
    K.short.flood <- list()
    K.long.flood <- list()
    for (i in 1:1000){
      Rec.short[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.short, m.par = m.par)
      Rec.long[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.long, m.par = m.par)
      Rec.short.flood[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.short.flood, m.par = m.par)
      Rec.long.flood[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.long.flood, m.par = m.par)
      Surv[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z, m.par = m.par)
      K.short[[i]] <- Surv[[i]] + Rec.short[[i]]
      K.long[[i]] <- Surv[[i]] + Rec.long[[i]]
      K.short.flood[[i]] <- Surv[[i]] + Rec.short.flood[[i]]
      K.long.flood[[i]] <- Surv[[i]] + Rec.long.flood[[i]]
    }
    rm(i)
  })
}

# IPM to examine effects of survival, carrying capacity, and frequency of reproduction
FUN.IPM.surv <- function(i.par, m.par){
  within(i.par, {
    # Create lists
    Rec.short <- list()
    Rec.long <- list()
    Surv10 <- list()
    Surv30 <- list()
    Surv50 <- list()
    Surv70 <- list()
    Surv90 <- list()
    Surv <- list()
    K.short.surv10 <- list()
    K.short.surv30 <- list()
    K.short.surv50 <- list()
    K.short.surv70 <- list()
    K.short.surv90 <- list()
    K.short <- list()
    K.long.surv10 <- list()
    K.long.surv30 <- list()
    K.long.surv50 <- list()
    K.long.surv70 <- list()
    K.long.surv90 <- list()
    K.long <- list()
    # Define kernels
    for (i in 1:1000){
      Rec.short[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.short, m.par = m.par)
      Rec.long[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.long, m.par = m.par)
      Surv10[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z.10, m.par = m.par)
      Surv30[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z.30, m.par = m.par)
      Surv50[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z.50, m.par = m.par)
      Surv70[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z.70, m.par = m.par)
      Surv90[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z.90, m.par = m.par)
      Surv[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z, m.par = m.par)
      K.short.surv10[[i]] <- Surv10[[i]] + Rec.short[[i]]
      K.short.surv30[[i]] <- Surv30[[i]] + Rec.short[[i]]
      K.short.surv50[[i]] <- Surv50[[i]] + Rec.short[[i]]
      K.short.surv70[[i]] <- Surv70[[i]] + Rec.short[[i]]
      K.short.surv90[[i]] <- Surv90[[i]] + Rec.short[[i]]
      K.short[[i]] <- Surv[[i]] + Rec.short[[i]]
      K.long.surv10[[i]] <- Surv10[[i]] + Rec.long[[i]]
      K.long.surv30[[i]] <- Surv30[[i]] + Rec.long[[i]]
      K.long.surv50[[i]] <- Surv50[[i]] + Rec.long[[i]]
      K.long.surv70[[i]] <- Surv70[[i]] + Rec.long[[i]]
      K.long.surv90[[i]] <- Surv90[[i]] + Rec.long[[i]]
      K.long[[i]] <- Surv[[i]] + Rec.long[[i]]
    }
    rm(i)
  })
}

# IPM to look at effects of different parameter values in the reproduction kernel
FUN.IPM.repro <- function(i.par, m.par){
  within(i.par, {
    # Create lists
    Rec.100.01.l <- list()
    Rec.100.02.l <- list()
    Rec.100.03.l <- list()
    Rec.100.04.l <- list()
    Rec.100.05.l <- list()
    Rec.75.01.l <- list()
    Rec.75.02.l <- list()
    Rec.75.03.l <- list()
    Rec.75.04.l <- list()
    Rec.75.05.l <- list()
    Rec.50.01.l <- list()
    Rec.50.02.l <- list()
    Rec.50.03.l <- list()
    Rec.50.04.l <- list()
    Rec.50.05.l <- list()
    Rec.100.01.s <- list()
    Rec.100.02.s <- list()
    Rec.100.03.s <- list()
    Rec.100.04.s <- list()
    Rec.100.05.s <- list()
    Rec.75.01.s <- list()
    Rec.75.02.s <- list()
    Rec.75.03.s <- list()
    Rec.75.04.s <- list()
    Rec.75.05.s <- list()
    Rec.50.01.s <- list()
    Rec.50.02.s <- list()
    Rec.50.03.s <- list()
    Rec.50.04.s <- list()
    Rec.50.05.s <- list()
    Surv <- list()
    K.100.01.l <- list()
    K.100.02.l <- list()
    K.100.03.l <- list()
    K.100.04.l <- list()
    K.100.05.l <- list()
    K.75.01.l <- list()
    K.75.02.l <- list()
    K.75.03.l <- list()
    K.75.04.l <- list()
    K.75.05.l <- list()
    K.50.01.l <- list()
    K.50.02.l <- list()
    K.50.03.l <- list()
    K.50.04.l <- list()
    K.50.05.l <- list()
    K.100.01.s <- list()
    K.100.02.s <- list()
    K.100.03.s <- list()
    K.100.04.s <- list()
    K.100.05.s <- list()
    K.75.01.s <- list()
    K.75.02.s <- list()
    K.75.03.s <- list()
    K.75.04.s <- list()
    K.75.05.s <- list()
    K.50.01.s <- list()
    K.50.02.s <- list()
    K.50.03.s <- list()
    K.50.04.s <- list()
    K.50.05.s <- list()
    # Define kernels
    for (i in 1:1000){
      Surv[[i]] <- h * outer(meshpoints, meshpoints, FUN.surv.z.50, m.par = m.par)
      Rec.100.01.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.01.l, m.par = m.par)
      Rec.100.02.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.02.l, m.par = m.par)
      Rec.100.03.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.03.l, m.par = m.par)
      Rec.100.04.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.04.l, m.par = m.par)
      Rec.100.05.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.05.l, m.par = m.par)
      Rec.75.01.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.01.l, m.par = m.par)
      Rec.75.02.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.02.l, m.par = m.par)
      Rec.75.03.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.03.l, m.par = m.par)
      Rec.75.04.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.04.l, m.par = m.par)
      Rec.75.05.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.05.l, m.par = m.par)
      Rec.50.01.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.01.l, m.par = m.par)
      Rec.50.02.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.02.l, m.par = m.par)
      Rec.50.03.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.03.l, m.par = m.par)
      Rec.50.04.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.04.l, m.par = m.par)
      Rec.50.05.l[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.05.l, m.par = m.par)
      Rec.100.01.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.01.s, m.par = m.par)
      Rec.100.02.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.02.s, m.par = m.par)
      Rec.100.03.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.03.s, m.par = m.par)
      Rec.100.04.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.04.s, m.par = m.par)
      Rec.100.05.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.100.05.s, m.par = m.par)
      Rec.75.01.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.01.s, m.par = m.par)
      Rec.75.02.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.02.s, m.par = m.par)
      Rec.75.03.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.03.s, m.par = m.par)
      Rec.75.04.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.04.s, m.par = m.par)
      Rec.75.05.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.75.05.s, m.par = m.par)
      Rec.50.01.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.01.s, m.par = m.par)
      Rec.50.02.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.02.s, m.par = m.par)
      Rec.50.03.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.03.s, m.par = m.par)
      Rec.50.04.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.04.s, m.par = m.par)
      Rec.50.05.s[[i]] <- h * outer(meshpoints, meshpoints, FUN.fec.z.50.05.s, m.par = m.par)
      K.100.01.l[[i]] <- Surv[[i]] + Rec.100.01.l[[i]]
      K.100.02.l[[i]] <- Surv[[i]] + Rec.100.02.l[[i]]
      K.100.03.l[[i]] <- Surv[[i]] + Rec.100.03.l[[i]]
      K.100.04.l[[i]] <- Surv[[i]] + Rec.100.04.l[[i]]
      K.100.05.l[[i]] <- Surv[[i]] + Rec.100.05.l[[i]]
      K.75.01.l[[i]] <- Surv[[i]] + Rec.75.01.l[[i]]
      K.75.02.l[[i]] <- Surv[[i]] + Rec.75.02.l[[i]]
      K.75.03.l[[i]] <- Surv[[i]] + Rec.75.03.l[[i]]
      K.75.04.l[[i]] <- Surv[[i]] + Rec.75.04.l[[i]]
      K.75.05.l[[i]] <- Surv[[i]] + Rec.75.05.l[[i]]
      K.50.01.l[[i]] <- Surv[[i]] + Rec.50.01.l[[i]]
      K.50.02.l[[i]] <- Surv[[i]] + Rec.50.02.l[[i]]
      K.50.03.l[[i]] <- Surv[[i]] + Rec.50.03.l[[i]]
      K.50.04.l[[i]] <- Surv[[i]] + Rec.50.04.l[[i]]
      K.50.05.l[[i]] <- Surv[[i]] + Rec.50.05.l[[i]]
      K.100.01.s[[i]] <- Surv[[i]] + Rec.100.01.s[[i]]
      K.100.02.s[[i]] <- Surv[[i]] + Rec.100.02.s[[i]]
      K.100.03.s[[i]] <- Surv[[i]] + Rec.100.03.s[[i]]
      K.100.04.s[[i]] <- Surv[[i]] + Rec.100.04.s[[i]]
      K.100.05.s[[i]] <- Surv[[i]] + Rec.100.05.s[[i]]
      K.75.01.s[[i]] <- Surv[[i]] + Rec.75.01.s[[i]]
      K.75.02.s[[i]] <- Surv[[i]] + Rec.75.02.s[[i]]
      K.75.03.s[[i]] <- Surv[[i]] + Rec.75.03.s[[i]]
      K.75.04.s[[i]] <- Surv[[i]] + Rec.75.04.s[[i]]
      K.75.05.s[[i]] <- Surv[[i]] + Rec.75.05.s[[i]]
      K.50.01.s[[i]] <- Surv[[i]] + Rec.50.01.s[[i]]
      K.50.02.s[[i]] <- Surv[[i]] + Rec.50.02.s[[i]]
      K.50.03.s[[i]] <- Surv[[i]] + Rec.50.03.s[[i]]
      K.50.04.s[[i]] <- Surv[[i]] + Rec.50.04.s[[i]]
      K.50.05.s[[i]] <- Surv[[i]] + Rec.50.05.s[[i]]
    }
    rm(i)
  })
}


# 2.5 Compile the IPM iteration matrix ---------------------------------------------------------------------------

# Reproduction IPM
set.seed(6754545)
IPM.iter.repro <- FUN.IPM.repro(i.par, m.par)

# Survival IPM
set.seed(354825)
IPM.iter.surv <- FUN.IPM.surv(i.par, m.par)

# Stochastic IPM
set.seed(944554)
IPM.iter.stoch <- FUN.IPM.stoch(i.par, m.par)


# 2.7 Save iteration matrix --------------------------------------------------------------------------------------

# Save Reproduction Iteration Matrix
saveRDS(IPM.iter.repro, file = "Results/IMP.iter.repro.RData")

# Save Survival Iteration Matrix
saveRDS(IPM.iter.surv, file = "Results/IMP.iter.surv.RData")

# Save Stochastic Iteration Matrix
saveRDS(IPM.iter.stoch, file = "Results/IMP.iter.stoch.RData")


# 2.8 Check to make sure eviction is not occurring in high numbers -----------------------------------------------

# Check to see if small individuals are being evicted
pnorm(80, mean = 33.0, sd = 2.7) - pnorm(25, mean = 33.0, sd = 2.7) # mean and SD based on shorter hydroperiod

# Check to see if large individuals are being evicted

# 1st using one of the constant survival functions
Surv.constant <- FUN.surv.50(meshpoints, m.par)
Surv.kernel.constant <- outer(meshpoints, meshpoints, FUN.surv.z.50, m.par = m.par) * h

# Plot results
plot(meshpoints, Surv.constant, xlab = "Size", type = "l", ylab = "Survival Probability", lwd = 12, ylim = c(0,0.8))
points(meshpoints, apply(Surv.kernel.constant, 2, sum), col = "red", lwd = 3, cex = 0.1, pch = 19)


# 2nd using one of the stochastic survival function
WrongPlace <- function(z, U) {
  fac1 <- FUN.surv.50(z, m.par)
  fac2 <- integrate(function(x) FUN.growth(x, z, m.par), U, Inf)$value
  return(fac1 * fac2)
}

zvals <- Wvals <- seq(25, 80, length = 100)
for(j in seq_along(zvals)) Wvals[j] <- WrongPlace(zvals[j], 80)
plot(zvals, Wvals, type = "l", lty = 1, lwd = 2, col = "black", xlab = "Initial size z", 
     ylab = "Fraction wrongfully evicted", ylim = c(0, 1.1 * max(Wvals)))

# Looks good!



#-----------------------------------------------------------------------------------------------------------------
# 3.0 Evaluate the effects of changes to the reproduction kernel on model output
#-----------------------------------------------------------------------------------------------------------------

# 3.1 Load IPM iteration matrix ----------------------------------------------------------------------------------

#IPM.iter.repro <- readRDS("Results/IMP.iter.repro.RData")

# 3.2 Calculate lambda values for all parameterizations of the kernel --------------------------------------------

# Calculate all lambda values

Repro.lam <- matrix(nrow = 1000, ncol = 30)

# Loop through to calculate lambda using eigen values
for(i in 1:1000){
  Repro.lam[i, 1] <- Re(eigen(IPM.iter.repro$K.100.01.l[[i]])$values[1])
  Repro.lam[i, 2] <- Re(eigen(IPM.iter.repro$K.100.02.l[[i]])$values[1])
  Repro.lam[i, 3] <- Re(eigen(IPM.iter.repro$K.100.03.l[[i]])$values[1])
  Repro.lam[i, 4] <- Re(eigen(IPM.iter.repro$K.100.04.l[[i]])$values[1])
  Repro.lam[i, 5] <- Re(eigen(IPM.iter.repro$K.100.05.l[[i]])$values[1])
  Repro.lam[i, 6] <- Re(eigen(IPM.iter.repro$K.75.01.l[[i]])$values[1])
  Repro.lam[i, 7] <- Re(eigen(IPM.iter.repro$K.75.02.l[[i]])$values[1])
  Repro.lam[i, 8] <- Re(eigen(IPM.iter.repro$K.75.03.l[[i]])$values[1])
  Repro.lam[i, 9] <- Re(eigen(IPM.iter.repro$K.75.04.l[[i]])$values[1])
  Repro.lam[i, 10] <- Re(eigen(IPM.iter.repro$K.75.05.l[[i]])$values[1])
  Repro.lam[i, 11] <- Re(eigen(IPM.iter.repro$K.50.01.l[[i]])$values[1])
  Repro.lam[i, 12] <- Re(eigen(IPM.iter.repro$K.50.02.l[[i]])$values[1])
  Repro.lam[i, 13] <- Re(eigen(IPM.iter.repro$K.50.03.l[[i]])$values[1])
  Repro.lam[i, 14] <- Re(eigen(IPM.iter.repro$K.50.04.l[[i]])$values[1])
  Repro.lam[i, 15] <- Re(eigen(IPM.iter.repro$K.50.05.l[[i]])$values[1])
  Repro.lam[i, 16] <- Re(eigen(IPM.iter.repro$K.100.01.s[[i]])$values[1])
  Repro.lam[i, 17] <- Re(eigen(IPM.iter.repro$K.100.02.s[[i]])$values[1])
  Repro.lam[i, 18] <- Re(eigen(IPM.iter.repro$K.100.03.s[[i]])$values[1])
  Repro.lam[i, 19] <- Re(eigen(IPM.iter.repro$K.100.04.s[[i]])$values[1])
  Repro.lam[i, 20] <- Re(eigen(IPM.iter.repro$K.100.05.s[[i]])$values[1])
  Repro.lam[i, 21] <- Re(eigen(IPM.iter.repro$K.75.01.s[[i]])$values[1])
  Repro.lam[i, 22] <- Re(eigen(IPM.iter.repro$K.75.02.s[[i]])$values[1])
  Repro.lam[i, 23] <- Re(eigen(IPM.iter.repro$K.75.03.s[[i]])$values[1])
  Repro.lam[i, 24] <- Re(eigen(IPM.iter.repro$K.75.04.s[[i]])$values[1])
  Repro.lam[i, 25] <- Re(eigen(IPM.iter.repro$K.75.05.s[[i]])$values[1])
  Repro.lam[i, 26] <- Re(eigen(IPM.iter.repro$K.50.01.s[[i]])$values[1])
  Repro.lam[i, 27] <- Re(eigen(IPM.iter.repro$K.50.02.s[[i]])$values[1])
  Repro.lam[i, 28] <- Re(eigen(IPM.iter.repro$K.50.03.s[[i]])$values[1])
  Repro.lam[i, 29] <- Re(eigen(IPM.iter.repro$K.50.04.s[[i]])$values[1])
  Repro.lam[i, 30] <- Re(eigen(IPM.iter.repro$K.50.05.s[[i]])$values[1])
}

# Format data for plotting
Repro.lam.data <- as.data.frame(Repro.lam)

colnames(Repro.lam.data) <- c("K.100.01.l", "K.100.02.l", "K.100.03.l", "K.100.04.l", "K.100.05.l",
                              "K.75.01.l", "K.75.02.l", "K.75.03.l", "K.75.04.l", "K.75.05.l",
                              "K.50.01.l", "K.50.02.l", "K.50.03.l", "K.50.04.l", "K.50.05.l",
                              "K.100.01.s", "K.100.02.s", "K.100.03.s", "K.100.04.s", "K.100.05.s",
                              "K.75.01.s", "K.75.02.s", "K.75.03.s", "K.75.04.s", "K.75.05.s",
                              "K.50.01.s", "K.50.02.s", "K.50.03.s", "K.50.04.s", "K.50.05.s")

Repro.lam.data <- Repro.lam.data %>% pivot_longer(cols = K.100.01.l:K.50.05.s, names_to = "Model", values_to = "Lambda") %>% 
  mutate(Label = rep(c("1%", "2%", "3%", "4%", "5%"), 6000))

# Save Lambda Data
save(Repro.lam.data, file = "Results/Repro_IPM_Lambda_RAW.rda")

# 3.3 Plot density of lambda values for each parameterization ----------------------------------------------------

# Load Data
load("Results/Repro_IPM_Lambda_RAW.rda")

# 100% of females return to breed each season
Repro.lam.100.plot.l <- ggplot(data = subset(Repro.lam.data, Model == c("K.100.01.l", "K.100.02.l", "K.100.03.l", 
                                                                        "K.100.04.l", "K.100.05.l")), 
                             aes(x = Lambda, fill = Label)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(1.0, 1.5, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs(tag = "(a)")

Repro.lam.100.plot.l


Repro.lam.100.plot.s <- ggplot(data = subset(Repro.lam.data, Model == c("K.100.01.s", "K.100.02.s", "K.100.03.s", 
                                                                        "K.100.04.s", "K.100.05.s")), 
                               aes(x = Lambda, fill = Label)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(1.0, 1.5, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs(tag = "(d)")

Repro.lam.100.plot.s


# 75% of females return to breed each season
Repro.lam.75.plot.l <- ggplot(data = subset(Repro.lam.data, Model == c("K.75.01.l", "K.75.02.l", "K.75.03.l",
                                                                       "K.75.04.l", "K.75.05.l")), 
                            aes(x = Lambda, fill = Label)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(1.0, 1.5, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs(tag = "(b)")

Repro.lam.75.plot.l


Repro.lam.75.plot.s <- ggplot(data = subset(Repro.lam.data, Model == c("K.75.01.s", "K.75.02.s", "K.75.03.s",
                                                                       "K.75.04.s", "K.75.05.s")), 
                            aes(x = Lambda, fill = Label)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(1.0, 1.5, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs(tag = "(e)")

Repro.lam.75.plot.s


# 50% of females return to breed each season
Repro.lam.50.plot.l <- ggplot(data = subset(Repro.lam.data, Model == c("K.50.01.l", "K.50.02.l", "K.50.03.l", 
                                                                       "K.50.04.l", "K.50.05.l")), 
                            aes(x = Lambda, fill = Label)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(1.0, 1.5, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs(tag = "(c)")

Repro.lam.50.plot.l


Repro.lam.50.plot.s <- ggplot(data = subset(Repro.lam.data, Model == c("K.50.01.s", "K.50.02.s", "K.50.03.s",
                                                                     "K.50.04.s", "K.50.05.s")), 
                            aes(x = Lambda, fill = Label)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(1.0, 1.5, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs(tag = "(f)")

Repro.lam.50.plot.s


# Save plot
tiff("Reproduction Lambdas Final.tiff", type = "cairo", width = 6.5, height = 4, units = "in", res = 600, 
     compression = "lzw")
Repro.lam.100.plot.l + Repro.lam.75.plot.l + Repro.lam.50.plot.l + Repro.lam.100.plot.s + Repro.lam.75.plot.s + 
  Repro.lam.50.plot.s + plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "bottom")
dev.off()


# 3.4 Calculate summary statistics for each parameterization -----------------------------------------------------

# Load raw data
load(file = "Results/Repro_IPM_Lambda_RAW.rda")

# Summarize data
Repro.lam.data.summary <- Repro.lam.data %>% 
  group_by(Model) %>% 
  summarize(Mean = mean(Lambda),
            Median = median(Lambda),
            Max = max(Lambda),
            Min = min(Lambda))
  
# Save summary of lambda values for each scenario
save(Repro.lam.data.summary, file = "Results/Repro_IPM_Lambda_summary.rda")


#-----------------------------------------------------------------------------------------------------------------
# 4.0 Examine the effects of survival and reproductive frequency on extinction risk
#-----------------------------------------------------------------------------------------------------------------

# 4.1 Load IPM iteration matrix ----------------------------------------------------------------------------------

IPM.iter.surv <- readRDS("Results/IMP.iter.surv.RData")

# 4.2 Calculate lambda values for each kernel --------------------------------------------------------------------

# Calculate all lambda values

Surv.lam <- matrix(nrow = 1000, ncol = 12)

# Loop through to calculate lambda using eigen values
for(i in 1:1000){
  Surv.lam[i, 1] <- Re(eigen(IPM.iter.surv$K.short.surv10[[i]])$values[1])
  Surv.lam[i, 2] <- Re(eigen(IPM.iter.surv$K.short.surv30[[i]])$values[1])
  Surv.lam[i, 3] <- Re(eigen(IPM.iter.surv$K.short.surv50[[i]])$values[1])
  Surv.lam[i, 4] <- Re(eigen(IPM.iter.surv$K.short.surv70[[i]])$values[1])
  Surv.lam[i, 5] <- Re(eigen(IPM.iter.surv$K.short.surv90[[i]])$values[1])
  Surv.lam[i, 6] <- Re(eigen(IPM.iter.surv$K.short[[i]])$values[1])
  Surv.lam[i, 7] <- Re(eigen(IPM.iter.surv$K.long.surv10[[i]])$values[1])
  Surv.lam[i, 8] <- Re(eigen(IPM.iter.surv$K.long.surv30[[i]])$values[1])
  Surv.lam[i, 9] <- Re(eigen(IPM.iter.surv$K.long.surv50[[i]])$values[1])
  Surv.lam[i, 10] <- Re(eigen(IPM.iter.surv$K.long.surv70[[i]])$values[1])
  Surv.lam[i, 11] <- Re(eigen(IPM.iter.surv$K.long.surv90[[i]])$values[1])
  Surv.lam[i, 12] <- Re(eigen(IPM.iter.surv$K.long[[i]])$values[1])
}

# Format data for plotting
Surv.lam.data <- as.data.frame(Surv.lam)

colnames(Surv.lam.data) <- c("Short.Surv10", "Short.Surv30", "Short.Surv50", "Short.Surv70", "Short.Surv90", "Short.Surv",
                             "Long.Surv10", "Long.Surv30", "Long.Surv50", "Long.Surv70", "Long.Surv90",  "Long.Surv")

Surv.lam.data <- Surv.lam.data %>% pivot_longer(cols = Short.Surv10:Long.Surv, names_to = "Model", values_to = "Lambda")

# Save Lambda Data
save(Surv.lam.data, file = "Results/Surv_IPM_Lambda_RAW.rda")

# 4.3 Plot density of lambda values for each parameterization ----------------------------------------------------

# Plot density of lambda values

# Short breeding seasons
Surv.lam.short.plot <- ggplot(data = subset(Surv.lam.data, Model == c("Short.Surv10", "Short.Surv30", "Short.Surv50", 
                                                                      "Short.Surv70", "Short.Surv90", "Short.Surv")), 
                                      aes(x = Lambda, fill = Model)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(0.8, 1.4, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme(legend.text = element_text(size = 8, color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(tag = "A)")


Surv.lam.short.plot

# Long breeding seasons
Surv.lam.long.plot <- ggplot(data = subset(Surv.lam.data, Model == c("Long.Surv10", "Long.Surv30", "Long.Surv50",
                                                                     "Long.Surv70", "Long.Surv90", "Long.Surv")), 
                              aes(x = Lambda, fill = Model)) + 
  geom_density(alpha = 0.5) + theme_classic() +
  xlab(expression(lambda)) + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 2.1), expand = c(0,0), breaks = c(0.8, 1.4, 2.0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")) + 
  theme(legend.text = element_text(size = 8, color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(tag = "B)")


Surv.lam.long.plot


# Save plot
tiff("Survival Lambdas Final.tiff", type = "cairo", width = 6.5, height = 5, units = "in", res = 600, 
     compression = "lzw")
Surv.lam.short.plot + Surv.lam.long.plot + plot_layout(ncol = 1) & theme(legend.position = "right")
dev.off()


# 4.4 Calculate summary statistics for each parameterization -----------------------------------------------------

# Load raw data
load(file = "Results/Surv_IPM_Lambda_RAW.rda")

# Summarize data
Surv.lam.data.summary <- Surv.lam.data %>% 
  group_by(Model) %>% 
  summarize(Mean = mean(Lambda),
            Median = median(Lambda),
            Max = max(Lambda),
            Min = min(Lambda))

# Save summary of lambda values for each scenario
save(Surv.lam.data.summary, file = "Results/Surv_IPM_Lambda_summary.rda")


# 4.5 Estimate extinction probability under different scenarios - long hydroperiod -------------------------------

# 4.5.1 50% survival and reproduction every other year -----------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv50.50 <- matrix(nrow = 70, ncol = 1000)

for(z in 1:1000){

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

# Create blank lists to store projection results
nt <- mods.nt <- list() 

# Create large loop to run the simulation X number of times
for(j in 1:n.sims){
  
  # Define vector for reproductive conditions in each year
  hydro <- sample(c(1,0), size = 70, replace = TRUE)
  
  # Randomly draw initial density and multiply by starting population size. Can do this one of two ways.
  nt0 <- Starting.Densities.SM[[1]] * start.pop
  
  # Assign starting population values
  nt[[1]] <- nt0
  
  # Run population model forward
  for(i in 1:70){
    
    # Does recruitment happen? If so, what is size distribution of metamorphs?
    if(hydro[i] == 0){
      
      nt[[i + 1]] <- (IPM.iter.surv$Surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    } else {
      
      nt[[i + 1]] <- (IPM.iter.surv$K.long.surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    }
    
    # Calculate pop size and reduce if above some carrying capacity
    n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
    
    if(n > k){
      
      nt[[i + 1]] <- nt[[i + 1]]/(n/k)
      
    }
  }
  
  # Compile results of each model run into a nested list object
  mods.nt[[j]] <- nt
  
}

# Calculate population size in each year
mods.n <- n.AUC <- list()
for(j in 1:n.sims){
  for(i in 1:71){
    
    # Calculate area under the curve for each year
    n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
  }
  
  # Combine n for each year
  mods.n[[j]] <- as.numeric(unlist(n.AUC))
}

# Calculate extinction probability for each year of the simulation
mods.ex.prob <- ex.prob <- list()
for(j in 1:n.sims){
  for(i in 2:71){
    
    # Does population go extinct in each year?
    ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
  }
  
  mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
}

# Create data frame and calculate annual values
mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)

# Final formatting
mods.ex.prob.df <- mods.ex.prob.df %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(Ex.Prob = (1 - value/n.sims))

mods.ex.prob.surv50.50[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))

}

# Create Summary of Results
mods.ex.prob.surv50.50.summary <- as.tibble(mods.ex.prob.surv50.50) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv50.50, file = "Results/Surv_IPM_exProb_50-50.rda")
save(mods.ex.prob.surv50.50.summary, file = "Results/Surv_IPM_exProb_50-50_summary.rda")


# 4.2.2 10% survival and reproduction every other year -----------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv10.50 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){

# Create blank lists to store projection results
nt <- kernels <- mods.nt <- mods.kernels <- list() 

# Create large loop to run the simulation X number of times
for(j in 1:n.sims){
  
  # Define vector for reproductive conditions in each year
  hydro <- sample(c(1, 0), size = 70, replace = TRUE)
  
  # Randomly draw initial density and multiply by starting population size.
  nt0 <- Starting.Densities.SM[[1]] * start.pop

  # Assign starting population values
  nt[[1]] <- nt0
  
  # Run population model forward
  for(i in 1:70){
    
    # Does recruitment happen? If so, what is size distribution of metamorphs?
    if(hydro[i] == 0){
      
      nt[[i + 1]] <- (IPM.iter.surv$Surv10[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    } else {
      
      nt[[i + 1]] <- (IPM.iter.surv$K.long.surv10[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    }
    
    # Calculate pop size and reduce if above some carrying capacity
    n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
    
    if(n > k){
      
      nt[[i + 1]] <- nt[[i + 1]]/(n/k)
      
    }
  }
  
  # Compile results of each model run into a nested list object
  mods.nt[[j]] <- nt
  
}

# Calculate population size in each year
mods.n <- n.AUC <- list()
for(j in 1:n.sims){
  for(i in 1:71){
    
    # Calculate area under the curve for each year
    n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
  }
  
  # Combine n for each year
  mods.n[[j]] <- as.numeric(unlist(n.AUC))
}

# Calculate extinction probability for each year of the simulation
mods.ex.prob <- ex.prob <- list()
for(j in 1:n.sims){
  for(i in 2:71){
    
    # Does population go extinct in each year?
    ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
  }
  
  mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
}

# Create data frame and calculate annual values
mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)

# Final formatting
mods.ex.prob.df <- mods.ex.prob.df %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(Ex.Prob = (1 - value/n.sims))

mods.ex.prob.surv10.50[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))

}

# Create Summary of Results
mods.ex.prob.surv10.50.summary <- as.tibble(mods.ex.prob.surv10.50) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv10.50, file = "Results/Surv_IPM_exProb_10-50.rda")
save(mods.ex.prob.surv10.50.summary, file = "Results/Surv_IPM_exProb_10-50_summary.rda")

# 4.2.3 90% survival and reproduction every other year -----------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv90.50 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){

# Create blank lists to store projection results
nt <- kernels <- mods.nt <- mods.kernels <- list() 

# Create large loop to run the simulation X number of times
for(j in 1:n.sims){
  
  # Define vector for reproductive conditions in each year
  hydro <- sample(c(1, 0), size = 70, replace = TRUE)
  
  # Randomly draw initial density and multiply by starting population size. Can do this one of two ways.
  nt0 <- Starting.Densities.SM[[1]] * start.pop
  #nt0 <-abs(Re(eigen(IPM.iter$K.long[[sample(1:1000, 1)]])$vector[,1])*100)
  
  # Assign starting population values
  nt[[1]] <- nt0
  
  # Run population model forward
  for(i in 1:70){
    
    # Does recruitment happen? If so, what is size distribution of metamorphs?
    if(hydro[i] == 0){
      
      nt[[i + 1]] <- (IPM.iter.surv$Surv90[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    } else {
      
      nt[[i + 1]] <- (IPM.iter.surv$K.long.surv90[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    }
    
    # Calculate pop size and reduce if above some carrying capacity
    n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
    
    if(n > k){
      
      nt[[i + 1]] <- nt[[i + 1]]/(n/k)
      
    }
  }
  
  # Compile results of each model run into a nested list object
  mods.nt[[j]] <- nt
  
}

# Calculate population size in each year
mods.n <- n.AUC <- list()
for(j in 1:n.sims){
  for(i in 1:71){
    
    # Calculate area under the curve for each year
    n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
  }
  
  # Combine n for each year
  mods.n[[j]] <- as.numeric(unlist(n.AUC))
}

# Calculate extinction probability for each year of the simulation
mods.ex.prob <- ex.prob <- list()
for(j in 1:n.sims){
  for(i in 2:71){
    
    # Does population go extinct in each year?
    ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
  }
  
  mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
}

# Create data frame and calculate annual values
mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)

# Final formatting
mods.ex.prob.df <- mods.ex.prob.df %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(Ex.Prob = (1 - value/n.sims))

mods.ex.prob.surv90.50[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))

}

# Create Summary of Results
mods.ex.prob.surv90.50.summary <- as.tibble(mods.ex.prob.surv90.50) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv90.50, file = "Results/Surv_IPM_exProb_90-50.rda")
save(mods.ex.prob.surv90.50.summary, file = "Results/Surv_IPM_exProb_90-50_summary.rda")


# 4.2.4 Variable survival and reproduction every other year ------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.survVAR.50 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){
  
  # Create blank lists to store projection results
  nt <- kernels <- mods.nt <- mods.kernels <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Define vector for reproductive conditions in each year
    hydro <- sample(c(1, 0), size = 70, replace = TRUE)
    
    # Randomly draw initial density and multiply by starting population size
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:70){
      
      # Does recruitment happen? If so, what is size distribution of metamorphs?
      if(hydro[i] == 0){
        
        nt[[i + 1]] <- (IPM.iter.surv$Surv[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      } else {
        
        nt[[i + 1]] <- (IPM.iter.surv$K.long[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      }
      
      # Calculate pop size and reduce if above some carrying capacity
      n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
      
      if(n > k){
        
        nt[[i + 1]] <- nt[[i + 1]]/(n/k)
        
      }
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:71){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  # Calculate extinction probability for each year of the simulation
  mods.ex.prob <- ex.prob <- list()
  for(j in 1:n.sims){
    for(i in 2:71){
      
      # Does population go extinct in each year?
      ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
    }
    
    mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
  }
  
  # Create data frame and calculate annual values
  mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)
  
  # Final formatting
  mods.ex.prob.df <- mods.ex.prob.df %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(where(is.numeric)) %>% 
    mutate(Ex.Prob = (1 - value/n.sims))
  
  mods.ex.prob.survVAR.50[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))
  
}

# Create Summary of Results
mods.ex.prob.survVAR.50.summary <- as.tibble(mods.ex.prob.survVAR.50) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.survVAR.50, file = "Results/Surv_IPM_exProb_VAR-50.rda")
save(mods.ex.prob.survVAR.50.summary, file = "Results/Surv_IPM_exProb_VAR-50_summary.rda")

# 4.2.5 50% survival and reproduction every 3rd year -------------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv50.33 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){

# Create blank lists to store projection results
nt <- mods.nt <- list() 

# Create large loop to run the simulation X number of times
for(j in 1:n.sims){
  
  # Define vector for reproductive conditions in each year
  hydro <- sample(c(1, 0, 0), size = 70, replace = TRUE)
  
  # Randomly draw initial density and multiply by starting population size.
  nt0 <- Starting.Densities.SM[[1]] * start.pop
  
  # Assign starting population values
  nt[[1]] <- nt0
  
  # Run population model forward
  for(i in 1:70){
    
    # Does recruitment happen? If so, what is size distribution of metamorphs?
    if(hydro[i] == 0){
      
      nt[[i + 1]] <- (IPM.iter.surv$Surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    } else {
      
      nt[[i + 1]] <- (IPM.iter.surv$K.long.surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    }
    
    # Calculate pop size and reduce if above some carrying capacity
    n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
    
    if(n > k){
      
      nt[[i + 1]] <- nt[[i + 1]]/(n/k)
      
    }
  }
  
  # Compile results of each model run into a nested list object
  mods.nt[[j]] <- nt
  
}

# Calculate population size in each year
mods.n <- n.AUC <- list()
for(j in 1:n.sims){
  for(i in 1:71){
    
    # Calculate area under the curve for each year
    n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
  }
  
  # Combine n for each year
  mods.n[[j]] <- as.numeric(unlist(n.AUC))
}

# Calculate extinction probability for each year of the simulation
mods.ex.prob <- ex.prob <- list()
for(j in 1:n.sims){
  for(i in 2:71){
    
    # Does population go extinct in each year?
    ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
  }
  
  mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
}

# Create data frame and calculate annual values
mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)

# Final formatting
mods.ex.prob.df <- mods.ex.prob.df %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(Ex.Prob = (1 - value/n.sims))

mods.ex.prob.surv50.33[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))

}

# Create Summary of Results
mods.ex.prob.surv50.33.summary <- as.tibble(mods.ex.prob.surv50.33) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv50.33, file = "Results/Surv_IPM_exProb_50-33.rda")
save(mods.ex.prob.surv50.33.summary, file = "Results/Surv_IPM_exProb_50-33_summary.rda")


# 4.2.6 10% survival and reproduction every 3rd year -------------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv10.33 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){

# Create blank lists to store projection results
nt <- kernels <- mods.nt <- mods.kernels <- list() 

# Create large loop to run the simulation X number of times
for(j in 1:n.sims){
  
  # Define vector for reproductive conditions in each year
  hydro <- sample(c(1, 0, 0), size = 70, replace = TRUE)
  
  # Randomly draw initial density and multiply by starting population size.
  nt0 <- Starting.Densities.SM[[1]] * start.pop

  # Assign starting population values
  nt[[1]] <- nt0
  
  # Run population model forward
  for(i in 1:70){
    
    # Does recruitment happen? If so, what is size distribution of metamorphs?
    if(hydro[i] == 0){
      
      nt[[i + 1]] <- (IPM.iter.surv$Surv10[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    } else {
      
      nt[[i + 1]] <- (IPM.iter.surv$K.long.surv10[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    }
    
    # Calculate pop size and reduce if above some carrying capacity
    n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
    
    if(n > k){
      
      nt[[i + 1]] <- nt[[i + 1]]/(n/k)
      
    }
  }
  
  # Compile results of each model run into a nested list object
  mods.nt[[j]] <- nt
  
}

# Calculate population size in each year
mods.n <- n.AUC <- list()
for(j in 1:n.sims){
  for(i in 1:71){
    
    # Calculate area under the curve for each year
    n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
  }
  
  # Combine n for each year
  mods.n[[j]] <- as.numeric(unlist(n.AUC))
}

# Calculate extinction probability for each year of the simulation
mods.ex.prob <- ex.prob <- list()
for(j in 1:n.sims){
  for(i in 2:71){
    
    # Does population go extinct in each year?
    ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
  }
  
  mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
}

# Create data frame and calculate annual values
mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)

# Final formatting
mods.ex.prob.df <- mods.ex.prob.df %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(Ex.Prob = (1 - value/n.sims))

mods.ex.prob.surv10.33[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))

}

# Create Summary of Results
mods.ex.prob.surv10.33.summary <- as.tibble(mods.ex.prob.surv10.33) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv10.33, file = "Results/Surv_IPM_exProb_10-33.rda")
save(mods.ex.prob.surv10.33.summary, file = "Results/Surv_IPM_exProb_10-33_summary.rda")

# 4.2.7 90% survival and reproduction every 3rd year -----------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv90.33 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 100; quasi.ext = 5; k = 100

for(z in 1:1000){
  
# Create blank lists to store projection results
nt <- kernels <- mods.nt <- mods.kernels <- list() 

# Create large loop to run the simulation X number of times
for(j in 1:n.sims){
  
  # Define vector for reproductive conditions in each year
  hydro <- sample(c(1, 0, 0), size = 70, replace = TRUE)
  
  # Randomly draw initial density and multiply by starting population size.
  nt0 <- Starting.Densities.SM[[1]] * start.pop
  
  # Assign starting population values
  nt[[1]] <- nt0
  
  # Run population model forward
  for(i in 1:70){
    
    # Does recruitment happen? If so, what is size distribution of metamorphs?
    if(hydro[i] == 0){
      
      nt[[i + 1]] <- (IPM.iter.surv$Surv90[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    } else {
      
      nt[[i + 1]] <- (IPM.iter.surv$K.long.surv90[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
      
    }
    
    # Calculate pop size and reduce if above some carrying capacity
    n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
    
    if(n > k){
      
      nt[[i + 1]] <- nt[[i + 1]]/(n/k)
      
    }
  }
  
  # Compile results of each model run into a nested list object
  mods.nt[[j]] <- nt
  
}

# Calculate population size in each year
mods.n <- n.AUC <- list()
for(j in 1:n.sims){
  for(i in 1:71){
    
    # Calculate area under the curve for each year
    n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
  }
  
  # Combine n for each year
  mods.n[[j]] <- as.numeric(unlist(n.AUC))
}

# Calculate extinction probability for each year of the simulation
mods.ex.prob <- ex.prob <- list()
for(j in 1:n.sims){
  for(i in 2:71){
    
    # Does population go extinct in each year?
    ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
  }
  
  mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
}

# Create data frame and calculate annual values
mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)

# Final formatting
mods.ex.prob.df <- mods.ex.prob.df %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(Ex.Prob = (1 - value/n.sims))

mods.ex.prob.surv90.33[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))

}

# Create Summary of Results
mods.ex.prob.surv90.33.summary <- as.tibble(mods.ex.prob.surv90.33) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv90.33, file = "Results/Surv_IPM_exProb_90-33.rda")
save(mods.ex.prob.surv90.33.summary, file = "Results/Surv_IPM_exProb_90-33_summary.rda")

# 4.2.8 Variable survival and reproduction every 3rd year --------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.survVAR.33 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){
  
  # Create blank lists to store projection results
  nt <- kernels <- mods.nt <- mods.kernels <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Define vector for reproductive conditions in each year
    hydro <- sample(c(1, 0, 0), size = 70, replace = TRUE)
    
    # Randomly draw initial density and multiply by starting population size.
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:70){
      
      # Does recruitment happen? If so, what is size distribution of metamorphs?
      if(hydro[i] == 0){
        
        nt[[i + 1]] <- (IPM.iter.surv$Surv[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      } else {
        
        nt[[i + 1]] <- (IPM.iter.surv$K.long[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      }
      
      # Calculate pop size and reduce if above some carrying capacity
      n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
      
      if(n > k){
        
        nt[[i + 1]] <- nt[[i + 1]]/(n/k)
        
      }
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:71){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  # Calculate extinction probability for each year of the simulation
  mods.ex.prob <- ex.prob <- list()
  for(j in 1:n.sims){
    for(i in 2:71){
      
      # Does population go extinct in each year?
      ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
    }
    
    mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
  }
  
  # Create data frame and calculate annual values
  mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)
  
  # Final formatting
  mods.ex.prob.df <- mods.ex.prob.df %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(where(is.numeric)) %>% 
    mutate(Ex.Prob = (1 - value/n.sims))
  
  mods.ex.prob.survVAR.33[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))
  
}

# Create Summary of Results
mods.ex.prob.survVAR.33.summary <- as.tibble(mods.ex.prob.survVAR.33) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.survVAR.33, file = "Results/Surv_IPM_exProb_VAR-33.rda")
save(mods.ex.prob.survVAR.33.summary, file = "Results/Surv_IPM_exProb_VAR-33_summary.rda")



# 4.2.9 50% survival and reproduction every 4th year -------------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv50.25 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){
  
  # Create blank lists to store projection results
  nt <- mods.nt <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Define vector for reproductive conditions in each year
    hydro <- sample(c(1, 0, 0, 0), size = 70, replace = TRUE)
    
    # Randomly draw initial density and multiply by starting population size.
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:70){
      
      # Does recruitment happen? If so, what is size distribution of metamorphs?
      if(hydro[i] == 0){
        
        nt[[i + 1]] <- (IPM.iter.surv$Surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      } else {
        
        nt[[i + 1]] <- (IPM.iter.surv$K.long.surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      }
      
      # Calculate pop size and reduce if above some carrying capacity
      n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
      
      if(n > k){
        
        nt[[i + 1]] <- nt[[i + 1]]/(n/k)
        
      }
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:71){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  # Calculate extinction probability for each year of the simulation
  mods.ex.prob <- ex.prob <- list()
  for(j in 1:n.sims){
    for(i in 2:71){
      
      # Does population go extinct in each year?
      ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
    }
    
    mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
  }
  
  # Create data frame and calculate annual values
  mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)
  
  # Final formatting
  mods.ex.prob.df <- mods.ex.prob.df %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(where(is.numeric)) %>% 
    mutate(Ex.Prob = (1 - value/n.sims))
  
  mods.ex.prob.surv50.25[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))
  
}

# Create Summary of Results
mods.ex.prob.surv50.25.summary <- as.tibble(mods.ex.prob.surv50.25) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv50.25, file = "Results/Surv_IPM_exProb_50-25.rda")
save(mods.ex.prob.surv50.25.summary, file = "Results/Surv_IPM_exProb_50-25_summary.rda")


# 4.2.6 10% survival and reproduction every 4th year -------------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv10.25 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){
  
  # Create blank lists to store projection results
  nt <- kernels <- mods.nt <- mods.kernels <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Define vector for reproductive conditions in each year
    hydro <- sample(c(1, 0, 0, 0), size = 70, replace = TRUE)
    
    # Randomly draw initial density and multiply by starting population size.
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:70){
      
      # Does recruitment happen? If so, what is size distribution of metamorphs?
      if(hydro[i] == 0){
        
        nt[[i + 1]] <- (IPM.iter.surv$Surv10[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      } else {
        
        nt[[i + 1]] <- (IPM.iter.surv$K.long.surv10[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      }
      
      # Calculate pop size and reduce if above some carrying capacity
      n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
      
      if(n > k){
        
        nt[[i + 1]] <- nt[[i + 1]]/(n/k)
        
      }
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:71){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  # Calculate extinction probability for each year of the simulation
  mods.ex.prob <- ex.prob <- list()
  for(j in 1:n.sims){
    for(i in 2:71){
      
      # Does population go extinct in each year?
      ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
    }
    
    mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
  }
  
  # Create data frame and calculate annual values
  mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)
  
  # Final formatting
  mods.ex.prob.df <- mods.ex.prob.df %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(where(is.numeric)) %>% 
    mutate(Ex.Prob = (1 - value/n.sims))
  
  mods.ex.prob.surv10.25[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))
  
}

# Create Summary of Results
mods.ex.prob.surv10.25.summary <- as.tibble(mods.ex.prob.surv10.25) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv10.25, file = "Results/Surv_IPM_exProb_10-25.rda")
save(mods.ex.prob.surv10.25.summary, file = "Results/Surv_IPM_exProb_10-25_summary.rda")

# 4.2.7 90% survival and reproduction every 4th year -----------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.surv90.25 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 100; quasi.ext = 5; k = 100

for(z in 1:1000){
  
  # Create blank lists to store projection results
  nt <- kernels <- mods.nt <- mods.kernels <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Define vector for reproductive conditions in each year
    hydro <- sample(c(1, 0, 0, 0), size = 70, replace = TRUE)
    
    # Randomly draw initial density and multiply by starting population size.
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:70){
      
      # Does recruitment happen? If so, what is size distribution of metamorphs?
      if(hydro[i] == 0){
        
        nt[[i + 1]] <- (IPM.iter.surv$Surv90[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      } else {
        
        nt[[i + 1]] <- (IPM.iter.surv$K.long.surv90[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      }
      
      # Calculate pop size and reduce if above some carrying capacity
      n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
      
      if(n > k){
        
        nt[[i + 1]] <- nt[[i + 1]]/(n/k)
        
      }
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:71){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  # Calculate extinction probability for each year of the simulation
  mods.ex.prob <- ex.prob <- list()
  for(j in 1:n.sims){
    for(i in 2:71){
      
      # Does population go extinct in each year?
      ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
    }
    
    mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
  }
  
  # Create data frame and calculate annual values
  mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)
  
  # Final formatting
  mods.ex.prob.df <- mods.ex.prob.df %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(where(is.numeric)) %>% 
    mutate(Ex.Prob = (1 - value/n.sims))
  
  mods.ex.prob.surv90.25[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))
  
}

# Create Summary of Results
mods.ex.prob.surv90.25.summary <- as.tibble(mods.ex.prob.surv90.25) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.surv90.25, file = "Results/Surv_IPM_exProb_90-25.rda")
save(mods.ex.prob.surv90.25.summary, file = "Results/Surv_IPM_exProb_90-25_summary.rda")

# 4.2.8 Variable survival and reproduction every 4th year --------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
mods.ex.prob.survVAR.25 <- matrix(nrow = 70, ncol = 1000)

# Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
n.sims = 100; start.pop = 50; quasi.ext = 5; k = 100

for(z in 1:1000){
  
  # Create blank lists to store projection results
  nt <- kernels <- mods.nt <- mods.kernels <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Define vector for reproductive conditions in each year
    hydro <- sample(c(1, 0, 0, 0), size = 70, replace = TRUE)
    
    # Randomly draw initial density and multiply by starting population size.
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:70){
      
      # Does recruitment happen? If so, what is size distribution of metamorphs?
      if(hydro[i] == 0){
        
        nt[[i + 1]] <- (IPM.iter.surv$Surv[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      } else {
        
        nt[[i + 1]] <- (IPM.iter.surv$K.long[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
      }
      
      # Calculate pop size and reduce if above some carrying capacity
      n <- DescTools::AUC(x = meshpoints, y = as.numeric(nt[[i + 1]]))
      
      if(n > k){
        
        nt[[i + 1]] <- nt[[i + 1]]/(n/k)
        
      }
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:71){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  # Calculate extinction probability for each year of the simulation
  mods.ex.prob <- ex.prob <- list()
  for(j in 1:n.sims){
    for(i in 2:71){
      
      # Does population go extinct in each year?
      ex.prob[[i]] <- ifelse(mods.n[[j]][[i]] >= quasi.ext, 1, 0)
    }
    
    mods.ex.prob[[j]] <- ifelse(cumany(as.numeric(unlist(ex.prob)) < 1), 0, 1)
  }
  
  # Create data frame and calculate annual values
  mods.ex.prob.df <- do.call(rbind.data.frame, mods.ex.prob)
  
  # Final formatting
  mods.ex.prob.df <- mods.ex.prob.df %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(where(is.numeric)) %>% 
    mutate(Ex.Prob = (1 - value/n.sims))
  
  mods.ex.prob.survVAR.25[,z] <- as.numeric(unlist(mods.ex.prob.df[,3]))
  
}

# Create Summary of Results
mods.ex.prob.survVAR.25.summary <- as.tibble(mods.ex.prob.survVAR.25) %>% 
  rowwise() %>% 
  summarise(Mean = mean(c_across(V1:V1000)),
            UC = quantile(c_across(V1:V1000), probs = 0.95),
            LC = quantile(c_across(V1:V1000), probs = 0.05))

# Save raw predictions and summary
save(mods.ex.prob.survVAR.25, file = "Results/Surv_IPM_exProb_VAR-25.rda")
save(mods.ex.prob.survVAR.25.summary, file = "Results/Surv_IPM_exProb_VAR-25_summary.rda")


# 4.2.9 Plot extinction probability across scenarios -------------------------------------------------------------

# Read in predictions
load("Results/Surv_IPM_exProb_50-50_summary.rda")
load("Results/Surv_IPM_exProb_10-50_summary.rda")
load("Results/Surv_IPM_exProb_90-50_summary.rda")
load("Results/Surv_IPM_exProb_VAR-50_summary.rda")
load("Results/Surv_IPM_exProb_50-33_summary.rda")
load("Results/Surv_IPM_exProb_10-33_summary.rda")
load("Results/Surv_IPM_exProb_90-33_summary.rda")
load("Results/Surv_IPM_exProb_VAR-33_summary.rda")

load("Results/Surv_IPM_exProb_50-25_summary.rda")
load("Results/Surv_IPM_exProb_10-25_summary.rda")
load("Results/Surv_IPM_exProb_90-25_summary.rda")
load("Results/Surv_IPM_exProb_VAR-25_summary.rda")

# Format data for plotting
Surv.ext.stack <- rbind(mods.ex.prob.surv50.50.summary, mods.ex.prob.surv10.50.summary, mods.ex.prob.surv90.50.summary, mods.ex.prob.survVAR.50.summary, 
                        mods.ex.prob.surv50.33.summary, mods.ex.prob.surv10.33.summary, mods.ex.prob.surv90.33.summary, mods.ex.prob.survVAR.33.summary,
                        mods.ex.prob.surv50.25.summary, mods.ex.prob.surv10.25.summary, mods.ex.prob.surv90.25.summary, mods.ex.prob.survVAR.25.summary) %>% 
  mutate(Survival = factor(rep(c("Median", "Low", "High", "Variable", "Median", "Low", "High", "Variable", "Median", "Low", "High", "Variable"), each = 70),
                           levels = c("High", "Median", "Low", "Variable")),
         Reproduction = factor(rep(c("50% Chance of Reproduction", "33% Chance of Reproduction", "25% Chance of Reproduction"), each = 280), 
                               levels = c("50% Chance of Reproduction", "33% Chance of Reproduction", "25% Chance of Reproduction")),
         Year = rep(1:70, 12))

# Plot results
Surv.IPM.graph <- ggplot(data = Surv.ext.stack, aes(x = Year, y = Mean, color = Survival)) + 
  geom_ribbon(aes(ymin = LC, ymax = UC), alpha = 0.2, color = NA) + 
  geom_line(size = 0.5) + theme_classic() +
  scale_color_manual(values = c("#440154FF", "#414487FF", "#2A788EFF", "#7AD151FF")) +
  scale_x_continuous(breaks = c(seq(10, 70, 10)), expand = c(0.02,0)) + 
  scale_y_continuous(limits = c(0, 1), breaks = c(seq(0, 1, 0.25))) + 
  ylab("Extinction Probability") + 
  theme_Publication() +
  theme(legend.position = "none") +
  facet_grid(Survival ~ Reproduction)


Surv.IPM.graph

# Save plot
tiff("Survival Extinction Probabilities k100.tiff", type = "cairo", width = 8, height = 5, units = "in", res = 600, 
     compression = "lzw")
Surv.IPM.graph
dev.off()





#-----------------------------------------------------------------------------------------------------------------
# 11.0 Conduct elasticity and sensitivity analyses
#-----------------------------------------------------------------------------------------------------------------

# 11.1 Calculate mean kernels for analyses -----------------------------------------------------------------------

# Load IPM iterations
IPM.iter.stoch <- readRDS("Results/IMP.iter.stoch.RData")

# Calculate average kernel using plyr functions
Mean.kernel.long <- plyr::aaply(plyr::laply(IPM.iter.stoch$K.long, as.matrix), c(2, 3), mean)
Mean.kernel.short <- plyr::aaply(plyr::laply(IPM.iter.stoch$K.short, as.matrix), c(2, 3), mean)
Mean.kernel.surv <- plyr::aaply(plyr::laply(IPM.iter.stoch$Surv, as.matrix), c(2, 3), mean)
Mean.rec.long <- plyr::aaply(plyr::laply(IPM.iter.stoch$Rec.long, as.matrix), c(2, 3), mean)
Mean.rec.short <- plyr::aaply(plyr::laply(IPM.iter.stoch$Rec.short, as.matrix), c(2, 3), mean)

# 11.2 Calculate Sensitivity and Elasticity of the kernels (long hydroperiod) ------------------------------------

# Calculate lambda
lambda <- Re(eigen(Mean.kernel.long)$values[1])

# Calculate right eigen vector
w.z <- Re(eigen(Mean.kernel.long)$vectors[, 1])

# Calculate left eigen vector
v.z1 <- Re(eigen(t(Mean.kernel.long))$vectors[, 1])

# Sensitivity
k.sens <- outer(v.z1, w.z, "*")/sum(v.z1 * w.z * h)

# Elasticity
k.elas <- k.sens * (Mean.kernel.long / h) / lambda

# Check to make sure it sums to 1
sum(k.elas) * h^2

# Calculate values to look at just the survival portion of the kernel
Pvals <- Mean.kernel.surv / h
P.elas <- Pvals * k.sens / lambda

# Calculate values to look at just the reproduction portion of the kernel
Fvals <- Mean.rec.long / h
F.elas <- Fvals * k.sens / lambda

# Calculate relative contributions of survival and fecundity
sum(P.elas) * h^2
sum(F.elas) * h^2

# Combine data
Sens.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(k.sens)))
Elas.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(k.elas)))
Elas.surv.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(P.elas)))
Elas.rec.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(F.elas)))

# Make contour plots
Sens.long.plot <- ggplot(data = Sens.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(a)")

Sens.long.plot

Elas.long.plot <- ggplot(data = Elas.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2, 
                          label.placer = metR::label_placer_random()) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(b)")

Elas.long.plot

Elas.long.surv.plot <- ggplot(data = Elas.surv.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2, 
                          label.placer = metR::label_placer_random()) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(c)")

Elas.long.surv.plot

Elas.long.rec.plot <- ggplot(data = Elas.rec.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(d)")

Elas.long.rec.plot

# Save plots
tiff("Sensitivity and elasticity long hydroperiod.tiff", type = "cairo", width = 6.5, height = 6.5, units = "in", res = 600, compression = "lzw")
Sens.long.plot + Elas.long.plot + Elas.long.surv.plot + Elas.long.rec.plot + plot_layout(ncol = 2)
dev.off()


# 11.3 Calculate Sensitivity and Elasticity of the kernels (short hydroperiod) ------------------------------------

# Calculate lambda
lambda <- Re(eigen(Mean.kernel.short)$values[1])

# Calculate right eigen vector
w.z <- Re(eigen(Mean.kernel.short)$vectors[, 1])

# Calculate left eigen vector
v.z1 <- Re(eigen(t(Mean.kernel.short))$vectors[, 1])

# Sensitivity
k.sens <- outer(v.z1, w.z, "*")/sum(v.z1 * w.z * h)

# Elasticity
k.elas <- k.sens * (Mean.kernel.short / h) / lambda

# Check to make sure it sums to 1
sum(k.elas) * h^2

# Calculate values to look at just the survival portion of the kernel
Pvals <- Mean.kernel.surv / h
P.elas <- Pvals * k.sens / lambda

# Calculate values to look at just the reproduction portion of the kernel
Fvals <- Mean.rec.short / h
F.elas <- Fvals * k.sens / lambda

# Calculate relative contributions of survival and fecundity
sum(P.elas) * h^2
sum(F.elas) * h^2

# Combine data
Sens.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(k.sens)))
Elas.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(k.elas)))
Elas.surv.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(P.elas)))
Elas.rec.data <- data.frame("X" = rep(meshpoints, 110), "Y" = rep(meshpoints, each = 110), "Value" = as.numeric(t(F.elas)))

# Make contour plots
Sens.short.plot <- ggplot(data = Sens.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(a)")

Sens.short.plot

Elas.short.plot <- ggplot(data = Elas.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(b)")

Elas.short.plot

Elas.short.surv.plot <- ggplot(data = Elas.surv.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(c)")

Elas.short.surv.plot

Elas.short.rec.plot <- ggplot(data = Elas.rec.data, aes(x = X, y = Y)) + theme_classic() +
  geom_tile() + geom_contour_filled(aes(z = Value)) + 
  metR::geom_text_contour(aes(z = Value), stroke = 0.2, stroke.colour = "white", size = 2) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("SVL (t), z") +
  ylab("SVL (t+1), z'") + 
  theme_Publication() +
  theme(legend.position = "none") + 
  labs(tag = "(d)")

Elas.short.rec.plot

# Save plots
tiff("Sensitivity and elasticity short hydroperiod TEST.tiff", type = "cairo", width = 6.5, height = 6.5, units = "in", res = 1200, compression = "lzw")
Sens.short.plot + Elas.short.plot + Elas.short.surv.plot + Elas.short.rec.plot + plot_layout(ncol = 2)
dev.off()


# 11.4 Calculate sensitivity and elasticity for individual parts of the kernel (UNUSED) --------------------------

#* ** THIS IS CURRENTLY CONFOUNDED BY THE STOCHASTIC FUNCTIONS. THE DRAWS ARE NOT THE SAME AS USED TO CALCULATE **
#* ** OVERALL KERNEL SENSITIVITY / ELASTICITY **


# Survival function

Surv.sens <- outer(meshpoints, meshpoints,
                   function(z1, z, m.par){(FUN.growth(z1, z, m.par) + ((1/2) * FUN.prob.breed(z, m.par) * 
                       FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.long(z1, m.par)))/2}, m.par)

Surv.sens.z <- apply(k.sens * Surv.sens, 2, sum) * h
Surv.elas.z <- Surv.sens.z * FUN.surv(meshpoints, m.par)/lambda


# Probability of recruitment

Prob.rec.sens <- outer(meshpoints, meshpoints,
                   function(z1, z, m.par){(FUN.surv(z, m.par) * (1/2) * FUN.eggs(z, m.par) * FUN.prob.breed(z, m.par) * 
                       FUN.metas.long(z1, m.par))/2}, m.par)

Prob.rec.sens.z <- apply(k.sens * Prob.rec.sens, 2, sum) * h
Prob.rec.elas.z <- Prob.rec.sens.z * FUN.prob.rec(100)/lambda


# 11.5 Plot function-specific perturbations for sensitivity and elasticity (UNUSED) ----------------------------------------

par(mfrow = c(1, 1))

# Survival function
plot(meshpoints, Surv.sens.z, type = "l", xlab = "SVL, z", 
     ylab = expression(paste("Sensitivity / Elasticity of  ", italic(s),"(",italic(z),")")))
lines(meshpoints, Surv.elas.z, lty = 2)
legend("topleft", legend = c("Sensitivity","Elasticity"), lty = c(1,2), bty = "n")


# Probability of recruitment
plot(meshpoints, Prob.rec.sens.z, type = "l", xlab = "SVL, z", 
     ylab = expression(paste("Sensitivity / Elasticity of ", italic(p)[italic(b)],"(",italic(z),")")))
lines(meshpoints, Prob.rec.elas.z, lty = 2)
legend("topleft", legend = c("Sensitivity","Elasticity"), lty = c(1,2), bty = "n")




#-----------------------------------------------------------------------------------------------------------------
# 12.0 Examine time to extinction from various carrying capacities
#-----------------------------------------------------------------------------------------------------------------

# 12.1 Load IPM iteration matrix ---------------------------------------------------------------------------------
IPM.iter.surv <- readRDS("Results/IMP.iter.surv.RData")

# 12.2 Run simulations for different values of K -----------------------------------------------------------------

# Create loop to generate confidence intervals and mean prediction
Pop.decline <- matrix(nrow = 20, ncol = 23)

# Create carrying capacity values
sizes <- seq(10, 230, 10)

for(z in 1:length(sizes)){
  
  # Define number of simulations, starting population size, quasi-extinction value, and carrying capacity
  n.sims = 100; start.pop = sizes[z]
  
  # Create blank lists to store projection results
  nt <- mods.nt <- list() 
  
  # Create large loop to run the simulation X number of times
  for(j in 1:n.sims){
    
    # Randomly draw initial density and multiply by starting population size. Can do this one of two ways.
    nt0 <- Starting.Densities.SM[[1]] * start.pop
    
    # Assign starting population values
    nt[[1]] <- nt0
    
    # Run population model forward
    for(i in 1:20){

        nt[[i + 1]] <- (IPM.iter.surv$Surv50[[sample(1:1000, 1)]]%*%nt[[i]])[,,drop=TRUE]
        
    }
    
    # Compile results of each model run into a nested list object
    mods.nt[[j]] <- nt
    
  }
  
  # Calculate population size in each year
  mods.n <- n.AUC <- list()
  for(j in 1:n.sims){
    for(i in 1:20){
      
      # Calculate area under the curve for each year
      n.AUC[[i]] <- DescTools::AUC(x = meshpoints, y = as.numeric(mods.nt[[j]][[i]]))
    }
    
    # Combine n for each year
    mods.n[[j]] <- as.numeric(unlist(n.AUC))
  }
  
  Pop.decline[,z] <- apply(X = as.data.frame(mods.n), MARGIN = 1, FUN = mean)
  
}

# Format data
Pop.decline.data <- ifelse(Pop.decline >= 5, 1, 0)
Pop.decline.data <- data.frame("Years" = colSums(Pop.decline.data),
                               "K" = sizes)

# Save data
save(Pop.decline.data, file = "Results/Pop_decline_data.rda")

# Plot results
Decline.graph <- ggplot(data = Pop.decline.data, aes(x = Years, y = K)) + theme_classic() + 
  #geom_line() + geom_point() + 
  geom_smooth(method = lm, formula = y ~ I(x^3), se = FALSE) +
  scale_x_continuous(limits = c(2, 13), breaks = c(3, 6, 9, 12)) +
  scale_y_continuous(limits = c(0, 230), breaks = c(0, 50, 100, 150, 200)) +
  xlab("Years before Extinction") +
  ylab("Population Size") +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        axis.ticks = element_line(color = "black"))

Decline.graph

# Save plots
tiff("Population Size and Extinction.tiff", type = "cairo", width = 3, height = 2, units = "in", res = 1200, compression = "lzw")
Decline.graph
dev.off()




##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################


matrix.image(x = meshpoints, y = meshpoints, IPM.iter$Rec.short[[sample(1:1000, 1)]])
matrix.image(x = meshpoints, y = meshpoints, IPM.iter$Surv[[sample(1:1000, 1)]])
matrix.image(x = meshpoints, y = meshpoints, IPM.iter$K.short[[sample(1:1000, 1)]])



stable.dist <- w.eigen/sum(w.eigen)

repro.val <- v.eigen/v.eigen[1]

v.dot.w <- sum(stable.dist * repro.val * h)



plot(x = meshpoints, 
     y = Starting.Densities.SM[[sample(1:20, 1)]] * 100)

# Could you just do this in one step and say K * pop density for each year?



sum(Starting.Densities.ME[[sample(1:20, 1)]])

Starting.Densities.SM[[sample(1:20, 1)]]%*%IPM.iter$Rec.short[[sample(1:1000, 1)]]

test <- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)

test%*%c(1,2,2)

  for (i in seq_len(n)) {
    xnew[[1]] <- (xnew[[1]] + Rec[[sample(1:1000, 1)]])}

  for (i in seq_len(n - 1)) {
    xnew[[i + 1]] <- (IPM.iter$Surv[[sample(1:1000, 1)]]%*%xnew[[i]])[,,drop = TRUE]}
  xnew[[n]] <- xnew[[n]] + (Surv[[sample(1:1000, 1)]]%*%xnew[[n - 1]])[,,drop = TRUE]
  return(xnew)
  


IPM.sim <- with(IPM.sys, {
  x <- nt0
  for(i in seq_len(100)){
    x1 <- r_iter(x, n, Rec, Surv)
    lam <- sum(unlist(x1))
    x <- lapply(x1, function(x){ x / lam})
  }
  list(lambda = lam, x = x)
})


plot(density(IPM.sim$x[[5]]))


test <- IPM.sim$x[[1]]

matrix.image(IPM.sim$x[[1]])

plot(density(IPM.sim$x[[1]]))


# starting density
nt0 <- with(i.par, lapply(1, function(i) rep(0, m)))
nt0[[1]] <- with(i.par, rep(1 / m, m))
x = nt0[[1]]
n = 1000

xnew <- list(0)
xnew[[1]] <- x

for (i in seq_len(n-1)) {
  xnew[[i+1]] = (IPM.sys$Surv[[sample(1:1000, 1)]]%*%xnew[[i]])[,,drop=TRUE] 
  xnew[[i+1]] = xnew[[i+1]] + (xnew[[i+1]]%*%IPM.sys$Rec[[sample(1:1000, 1)]])[,,drop=TRUE]
}


# Cohort distributions
## compute the size distribution for each age class with repetitious code...
lam.est  <- Re(eigen(IPM.sys$K[[1]])$values[1])

a0.z.dist.est <- IPM.sys$Rec[[1]] %*% stable.z.dist.est / lam.est
a1.z.dist.est <- IPM.sys$Surv[[1]] %*% a0.z.dist.est / lam.est
a2.z.dist.est <- IPM.sys$Surv[[1]] %*% a1.z.dist.est / lam.est
a3.z.dist.est <- IPM.sys$Surv[[1]] %*% a2.z.dist.est / lam.est

## build a little helper function to compute the means & variances
mk_moments <- function(z.dist, meshpoints) {
  z.dist <- z.dist/sum(z.dist)
  mean.z <- sum(z.dist * meshpoints)
  var.z  <- sum(z.dist * meshpoints^2) - mean.z^2
  return(c(mean = mean.z, sd = sqrt(var.z)))
}

##
set_graph_pars(ptype = "panel1")
## age = 0 (new recruits)
z.dist <- z.dist.by.age[[1]]
plot(meshpoints, z.dist, type="n", xlab="Mass, z", ylab="Density",xlim=c(1.5,4.2))
## age = 1, 2, 3 and 4
for (A in 0:4) {
  z.dist <-  z.dist.by.age[[A+1]]
  lines(meshpoints, z.dist)
  moments <- round(mk_moments(z.dist, meshpoints), 2)
  text(x=moments["mean"], y=max(z.dist)+5e-5, pos=4, cex=0.75,
       labels=paste("A = ", A, " (mean = ", moments["mean"], ", s.d. = ", moments["sd"],")", sep=""))
}


plot(meshpoints, a0.z.dist.est)
lines(meshpoints, a1.z.dist.est)
lines(meshpoints, a2.z.dist.est)
lines(meshpoints, a3.z.dist.est)







lam = vector()
for(i in 1:(n-1)){
  lam[i] = (sum(unlist(xnew[i+1]))/sum(unlist(xnew[i])))*2
}
lam = as.data.frame(lam)
lam = as.data.frame(lam[!is.na(lam)])
names(lam) = "lam"




