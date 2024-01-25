##################################################################################################################
# Modeling Reticulated Flatwoods Salamander populations: Assessing effects of future climate change scenarios
# on long-term population viability. This script uses the functions, parameters, and results from the two other 
# scripts to graph results from the demographic functions and IPM.
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
# 2.0 
# 3.0 
# 4.0 
# 5.0


#-----------------------------------------------------------------------------------------------------------------
# 2.0 Fecundity kernel - Eggs produced, survival until metamorphosis, size at metamorphosis
#-----------------------------------------------------------------------------------------------------------------

# 2.1 Plot relationship between female body size and number of eggs produced -------------------------------------

# Run function a bunch of times
Eggs.reps <- replicate(1000, FUN.eggs(meshpoints, m.par))

# Calculate quantiles and means and format data
Eggs.reps.plot.data <- cbind(as.data.frame(t(apply(Eggs.reps, 1, quantile, c(0.05, 0.95)))), 
                              "Mean" = rowMeans(Eggs.reps),
                              "SVL" = meshpoints) %>% 
  filter(Mean > 0)

# Plot results
Eggs.plot <- ggplot(data = Eggs.reps.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  #scale_y_continuous(limits = c(0, 0.16)) +
  scale_x_continuous(limits = c(45, 80)) +
  xlab("Snout-vent Length") + ylab("Number of Eggs") + 
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

Eggs.plot


# 2.2 Plot how many eggs survive until metamorphosis -------------------------------------------------------------

Prob.rec.plot <- FUN.prob.rec(100)

plot(Prob.rec.plot, ylab = "Probability of Survival until Metamorphosis", ylim = c(0, 0.06))

# 2.3 Plot of salamander breeding potential ----------------------------------------------------------------------

# Run function a bunch of times
prob.breed.reps <- replicate(1000, FUN.prob.breed.100(meshpoints, m.par))

# Calculate quantiles and means and format data
prob.breed.plot.data <- cbind(as.data.frame(t(apply(prob.breed.reps, 1, quantile, c(0.05, 0.95)))), 
                             "Mean" = rowMeans(prob.breed.reps),
                             "SVL" = meshpoints) 

Prob.breed.plot <- ggplot(data = prob.breed.plot.data, aes(x = SVL, y = Mean)) + 
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  theme_classic() +
  xlab("Female Snout-vent Length (mm)") + ylab("Probability of Breeding") + 
  scale_x_continuous(breaks = seq(30, 80, 10)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

Prob.breed.plot


# 2.4 Plot size at metamorphosis based on hydroperiod ------------------------------------------------------------

# Run function a bunch of times
Long.meta.reps <- replicate(1000, FUN.metas.long(meshpoints, m.par))

# Calculate quantiles and means and format data
Long.meta.plot.data <- cbind(as.data.frame(t(apply(Long.meta.reps, 1, quantile, c(0.05, 0.95)))), 
                             "Mean" = rowMeans(Long.meta.reps),
                             "SVL" = meshpoints) 

# Plot results
Long.meta.plot <- ggplot(data = Long.meta.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(0, 0.12)) +
  scale_x_continuous(limits = c(25, 60)) +
  xlab("Snout-vent Length (mm)") + ylab("Density") +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))


Long.meta.plot


# Run function a bunch of times
Short.meta.reps <- replicate(1000, FUN.metas.short(meshpoints, m.par))

# Calculate quantiles and means and format data
Short.meta.plot.data <- cbind(as.data.frame(t(apply(Short.meta.reps, 1, quantile, c(0.05, 0.95)))), 
                             "Mean" = rowMeans(Short.meta.reps),
                             "SVL" = meshpoints) 

# Plot results
Short.meta.plot <- ggplot(data = Short.meta.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(0, 0.16)) +
  scale_x_continuous(limits = c(25, 50)) +
  xlab("Snout-vent Length (mm)") + ylab("Density") +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))


Short.meta.plot



# 2.5 Plot the 'average' fecundity kernel ------------------------------------------------------------------------

# Load IPM
IPM.iter.stoch <- readRDS("Results/IMP.iter.stoch.RData")

# Calculate average kernel using plyr functions
Mean.fec.kernel = plyr::aaply(plyr::laply(IPM.iter.stoch$Rec.long, as.matrix), c(2, 3), mean)

# Plot the kernel using matrix image function
# matrix.image(Mean.fec.kernel, meshpoints, meshpoints)

# Plot matrix using ggplot
Fec.kernel.data.plot <- data.frame("SVL.z1" = rep(meshpoints, 110), 
                                   "SVL.z" = rep(meshpoints, each = 110),
                                   "Value" = c(Mean.fec.kernel))

Fec.kernel.plot <- ggplot(data = Fec.kernel.data.plot, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) + 
  xlab("Snout-vent length at time t (mm)") + ylab("Snout-vent length at time t+1 (mm)") +
  scale_fill_viridis_c(direction = 1) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title.align = 0.5)


Fec.kernel.plot


#-----------------------------------------------------------------------------------------------------------------
# 3.0 Growth kernel - How do salamanders grow after metamorphosis?
#-----------------------------------------------------------------------------------------------------------------

# 3.1 Plot growth function ---------------------------------------------------------------------------------------

# Define values from posterior distributions
Linf.05 <- as.numeric(quantile(m.par$Linf, 0.05))
Linf.50 <- as.numeric(mean(m.par$Linf))
Linf.95 <- as.numeric(quantile(m.par$Linf, 0.95))

K.05 <- as.numeric(quantile(m.par$k.mu, 0.05))
K.50 <- as.numeric(mean(m.par$k.mu))
K.95 <- as.numeric(quantile(m.par$k.mu, 0.95))

# Create data
Growth.plot.data <- matrix(nrow = 15, ncol = 3); Growth.plot.data[1,] <- 30

for(i in 2:15){Growth.plot.data[i,1] <- Growth.plot.data[i-1,1] + (Linf.05 - Growth.plot.data[i-1,1]) * (1 - exp(-K.05))}
for(i in 2:15){Growth.plot.data[i,2] <- Growth.plot.data[i-1,2] + (Linf.50 - Growth.plot.data[i-1,2]) * (1 - exp(-K.50))}
for(i in 2:15){Growth.plot.data[i,3] <- Growth.plot.data[i-1,3] + (Linf.95 - Growth.plot.data[i-1,3]) * (1 - exp(-K.95))}

Growth.plot.data <- as.data.frame(Growth.plot.data) %>% mutate(Age = 1:15)
colnames(Growth.plot.data) <- c("LCL", "Mean", "UCL", "Age")

# Plot graph
Growth.plot <- ggplot(data = Growth.plot.data, aes(x = Age, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(25, 75), expand = c(0,0)) +
  scale_x_continuous(limits = c(0.5, 10.5), expand = c(0,0), breaks = seq(0, 10, 2)) +
  xlab("Age") + ylab("Snout-vent Length (mm)") +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

Growth.plot


# Can also plot by size to show transition probabilities
plot(x = meshpoints, y = FUN.growth(meshpoints, 75, m.par), xlab = "SVL", ylab = "Density", type = "l", xlim = c(45, 80))
lines(x = meshpoints, y = FUN.growth(meshpoints, 75, m.par))
lines(x = meshpoints, y = FUN.growth(meshpoints, 75, m.par))
lines(x = meshpoints, y = FUN.growth(meshpoints, 75, m.par))
lines(x = meshpoints, y = FUN.growth(meshpoints, 75, m.par))

# Plot an example of the growth kernel (this is a stochastic draw so changes every time)
Growth.matrix <- outer(meshpoints, meshpoints, FUN.growth, m.par = m.par) * h

Growth.matrix.plot.data <- data.frame("SVL.z1" = rep(meshpoints, 110), 
                                   "SVL.z" = rep(meshpoints, each = 110),
                                   "Value" = c(Growth.matrix))

Growth.matrix.plot <- ggplot(data = Growth.matrix.plot.data, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_reverse(limits = c(80, 25), breaks = seq(80, 30, -10), expand = c(0,0)) + 
  xlab("Snout-vent length at time t (mm)") + ylab("Snout-vent length at time t+1 (mm)") +
  scale_fill_viridis_c(direction = 1) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title.align = 0.5)

Growth.matrix.plot

#-----------------------------------------------------------------------------------------------------------------
# 4.0 Adult survival kernel - Brooks 2020
#-----------------------------------------------------------------------------------------------------------------

# 4.1 Plot survival function -------------------------------------------------------------------------------------

plot(x = meshpoints, y = FUN.surv(meshpoints, m.par), xlab = "SVL", ylab = "Probability of Survival")

# Run function a bunch of times
Surv.reps <- replicate(1000, FUN.surv(meshpoints, m.par))

# Calculate quantiles and means and format data
Surv.plot.data <- cbind(as.data.frame(t(apply(Surv.reps, 1, quantile, c(0.05, 0.95)))), 
                              "Mean" = rowMeans(Surv.reps),
                              "SVL" = meshpoints) 

# Plot results
Surv.plot <- ggplot(data = Surv.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0.2, 0.8, 0.2)) +
  scale_x_continuous(limits = c(25, 80)) +
  xlab("Snout-vent Length (mm)") + ylab("Survival Probility") +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

Surv.plot


# 4.2 Plot combined survival and growth kernel (mean kernel across all simulations)

# Load IPM
IPM.iter.stoch <- readRDS("Results/IMP.iter.stoch.RData")

# Calculate average kernel using plyr functions
Mean.surv.kernel = plyr::aaply(plyr::laply(IPM.iter.stoch$Surv, as.matrix), c(2, 3), mean)

# Plot matrix using ggplot
Surv.kernel.data.plot <- data.frame("SVL.z1" = rep(meshpoints, 110), 
                                   "SVL.z" = rep(meshpoints, each = 110),
                                   "Value" = c(Mean.surv.kernel))

Surv.growth.matrix.plot <- ggplot(data = Surv.kernel.data.plot, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_reverse(limits = c(80, 25), breaks = seq(80, 30, -10), expand = c(0,0)) + 
  xlab("Snout-vent length at time t (mm)") + ylab("Snout-vent length at time t+1 (mm)") +
  scale_fill_viridis_c(direction = 1) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title.align = 0.5)

Surv.growth.matrix.plot


#-----------------------------------------------------------------------------------------------------------------
# 5.0 Combined kernel
#-----------------------------------------------------------------------------------------------------------------

# 5.1 Plot the kernel with long hydroperiod reproduction ---------------------------------------------------------

# Load IPM
IPM.iter.stoch <- readRDS("Results/IMP.iter.stoch.RData")

# Calculate average kernel using plyr functions
Mean.kernel.long = plyr::aaply(plyr::laply(IPM.iter.stoch$K.long, as.matrix), c(2, 3), mean)

# Plot matrix using ggplot
Mean.kernel.long.data.plot <- data.frame("SVL.z1" = rep(meshpoints, 110), 
                                    "SVL.z" = rep(meshpoints, each = 110),
                                    "Value" = c(Mean.kernel.long))

Mean.kernel.long.matrix.plot <- ggplot(data = Mean.kernel.long.data.plot, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_reverse(limits = c(80, 25), breaks = seq(80, 30, -10), expand = c(0,0)) + 
  xlab("Snout-vent length at time t (mm)") + ylab("Snout-vent length at time t+1 (mm)") +
  scale_fill_viridis_c(direction = 1) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title.align = 0.5)

Mean.kernel.long.matrix.plot


# 5.2 Plot the kernel with short hydroperiod reproduction --------------------------------------------------------

# Load IPM
IPM.iter.stoch <- readRDS("Results/IMP.iter.stoch.RData")

# Calculate average kernel using plyr functions
Mean.kernel.short = plyr::aaply(plyr::laply(IPM.iter.stoch$K.short, as.matrix), c(2, 3), mean)

# Plot matrix using ggplot
Mean.kernel.short.data.plot <- data.frame("SVL.z1" = rep(meshpoints, 110), 
                                         "SVL.z" = rep(meshpoints, each = 110),
                                         "Value" = c(Mean.kernel.short))

Mean.kernel.short.matrix.plot <- ggplot(data = Mean.kernel.short.data.plot, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_reverse(limits = c(80, 25), breaks = seq(80, 30, -10), expand = c(0,0)) + 
  xlab("Snout-vent length at time t (mm)") + ylab("Snout-vent length at time t+1 (mm)") +
  scale_fill_viridis_c(direction = 1) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title.align = 0.5)

Mean.kernel.short.matrix.plot


#-----------------------------------------------------------------------------------------------------------------
# 6.0 Arrange and format about plots for Figure(s)
#-----------------------------------------------------------------------------------------------------------------

# 6.1 Fecundity kernel -------------------------------------------------------------------------------------------

# Eggs per female
Eggs.plot <- ggplot(data = Eggs.reps.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  #scale_y_continuous(limits = c(0, 0.16)) +
  scale_x_continuous(limits = c(45, 80)) +
  xlab("Size (mm)") + ylab("Number \nof eggs") + 
  theme_Publication() + 
  labs(tag = "(c)")

Eggs.plot  

# Probability of breeding
Prob.breed.plot <- ggplot(data = prob.breed.plot.data, aes(x = SVL, y = Mean)) + 
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + 
  geom_line() + 
  theme_classic() +
  xlab("Size (mm)") + ylab("Breeding \nprobability") + 
  scale_x_continuous(breaks = seq(30, 80, 10)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  theme_Publication()+ 
  labs(tag = "(d)")

Prob.breed.plot

# Metamorphs size distributions
Long.meta.plot <- ggplot(data = Long.meta.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(0, 0.12), breaks = c(0, 0.04, 0.08, 0.12)) +
  scale_x_continuous(limits = c(25, 60)) +
  xlab("Size (mm)") + ylab("Density") +
  theme_Publication()

Short.meta.plot <- ggplot(data = Short.meta.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(0, 0.16), breaks = c(0, 0.05, 0.10, 0.15)) +
  scale_x_continuous(limits = c(25, 50)) +
  xlab("Size (mm)") + ylab("Density") +
  theme_Publication()

# Fecundity kernel
Fec.kernel.plot <- ggplot(data = Fec.kernel.data.plot, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) + 
  xlab("Size (t)") + ylab("Size (t+1)") +
  scale_fill_viridis_c(direction = 1) +
  theme_Publication() + 
  theme(legend.title = element_blank()) +
  labs(tag = "(f)")

# Save plots
tiff("Fecundity kernel.tiff", type = "cairo", width = 6.5, height = 4, units = "in", res = 1200, compression = "lzw")
((Eggs.plot | Prob.breed.plot) / (Short.meta.plot + Long.meta.plot)) + Fec.kernel.plot + plot_layout(heights = c(1,1,2.5))
dev.off()


# 6.2 Survival and Growth Kernels --------------------------------------------------------------------------------


# von-Bertalanffy Growth Curve
Growth.plot <- ggplot(data = Growth.plot.data, aes(x = Age, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(25, 75), expand = c(0,0)) +
  scale_x_continuous(limits = c(0.5, 10.5), expand = c(0,0), breaks = seq(0, 10, 2)) +
  xlab("Age (years)") + ylab("Size (mm)") +
  theme_Publication()+ 
  labs(tag = "(a)")

# Survival plot
Surv.plot <- ggplot(data = Surv.plot.data, aes(x = SVL, y = Mean)) + theme_classic() +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_line() +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0.2, 0.8, 0.2)) +
  scale_x_continuous(limits = c(25, 80)) +
  xlab("Size (mm)") + ylab("Survival \nprobability") +
  theme_Publication()+ 
  labs(tag = "(b)")

# Plot the kernel combined 
Surv.growth.matrix.plot <- ggplot(data = Surv.kernel.data.plot, aes(x = SVL.z, y = SVL.z1, fill = Value)) + 
  geom_tile() + theme_classic() + 
  scale_x_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 80), breaks = seq(30, 80, 10), expand = c(0,0)) + 
  xlab("Size (t)") + ylab("Size (t+1)") +
  scale_fill_viridis_c(direction = 1) +
  theme_Publication() + 
  theme(legend.title = element_blank()) +
  labs(tag = "(e)")

# Save plots
tiff("Survival kernel.tiff", type = "cairo", width = 6.5, height = 4, units = "in", res = 1200, compression = "lzw")
(Growth.plot | Surv.plot) / Surv.growth.matrix.plot + plot_layout(heights = c(1,2.5))
dev.off()


# all combined plot
tiff("both kernel.tiff", type = "cairo", width = 6, height = 6.75, units = "in", res = 600, compression = "lzw")
(Growth.plot | Surv.plot) / (Eggs.plot | Prob.breed.plot) / (Surv.growth.matrix.plot + Fec.kernel.plot) + plot_layout(heights = c(1,1,1))
dev.off()


#-----------------------------------------------------------------------------------------------------------------
# 8.0 Starting Densities
#-----------------------------------------------------------------------------------------------------------------

# Load data
Starting.Densities <- readRDS("Data/Starting_Densities.rds")

# Format smooth data
Starting.Densities.SM <- split(Starting.Densities$Density.Sm, f = Starting.Densities$Dummy)
Starting.Densities.SM.plot <- data.frame("Values" = as.numeric(unlist(Starting.Densities.SM)),
                                         "Number" = rep(1:20, each = 110),
                                         "Meshpoints" = rep(meshpoints, 20))



# Create plot

Start.Density.plot <- ggplot(data = Starting.Densities.SM.plot, aes(x = Meshpoints, y = Values)) + theme_classic() + 
  geom_line(size = 1) + facet_wrap(vars(Number)) + 
  xlab("SVL") + 
  ylab("Density") + 
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"))


Start.Density.plot

# Save plots
tiff("Starting Densities.tiff", type = "cairo", width = 6.5, height = 4, units = "in", res = 600, compression = "lzw")
Start.Density.plot
dev.off()


#-----------------------------------------------------------------------------------------------------------------
# 9.0 Stable size distribution
#-----------------------------------------------------------------------------------------------------------------

# 9.1 Load IPM and calculate distribution ------------------------------------------------------------------------
IPM.iter.stoch <- readRDS("Results/IMP.iter.stoch.RData")

# Calculate average kernel using plyr functions
Mean.kernel.long <- plyr::aaply(plyr::laply(IPM.iter.stoch$K.long, as.matrix), c(2, 3), mean)
Mean.kernel.short <- plyr::aaply(plyr::laply(IPM.iter.stoch$K.short, as.matrix), c(2, 3), mean)
Mean.kernel.surv <- plyr::aaply(plyr::laply(IPM.iter.stoch$Surv, as.matrix), c(2, 3), mean)
Mean.kernel.combined <- plyr::aaply(plyr::laply(append(IPM.iter.stoch$Surv, IPM.iter.stoch$K.short),as.matrix), c(2, 3), mean)

# Calculate stable size distribution for long hydrperiods
w.eigen.long <- Re(eigen(Mean.kernel.long)$vectors[,1])
stable.dist.long <- w.eigen.long/sum(w.eigen.long) 
Stable.dist.long.plot <- data.frame("Meshpoints" = meshpoints,
                                    "Values" = stable.dist.long)

# Calculate stable size distribution for long hydrperiods
w.eigen.short <- Re(eigen(Mean.kernel.short)$vectors[,1])
stable.dist.short <- w.eigen.short/sum(w.eigen.short) 
Stable.dist.short.plot <- data.frame("Meshpoints" = meshpoints,
                                    "Values" = stable.dist.short)

# Plot distributions
Stable.size.long.plot <- ggplot(data = Stable.dist.long.plot, aes(x = Meshpoints, y = Values)) + theme_classic() +
  geom_line(size = 1) + 
  xlab("SVL") + 
  ylab("Density") + 
  scale_x_continuous(breaks = c(30, 45, 60, 75)) +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(tag = "A)")


Stable.size.long.plot


Stable.size.short.plot <- ggplot(data = Stable.dist.short.plot, aes(x = Meshpoints, y = Values)) + theme_classic() +
  geom_line(size = 1) + 
  xlab("SVL") + 
  ylab("Density") + 
  scale_x_continuous(breaks = c(30, 45, 60, 75)) +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(tag = "B)")


Stable.size.short.plot


# Save plots
tiff("Stable Size Distributions.tiff", type = "cairo", width = 6.5, height = 3, units = "in", res = 600, compression = "lzw")
Stable.size.long.plot + Stable.size.short.plot
dev.off()
