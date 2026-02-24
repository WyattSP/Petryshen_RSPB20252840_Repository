# Library Imports
library(divDyn)
library(reshape2)
library(chronosphere)
library(dplyr)
library(scales)

##### #### #####
# Data Import  #
##### #### #####
dynMat <- read.csv("./Blind_Review_Code_Repository_2025/data/ForDivDyn/Sim_1000_10.csv",  check.names = FALSE)
colnames(dynMat) <- c("taxon", "bin", "mid")
dynMat$taxon <- dynMat$taxon + 1 # Add 1 to all values of taxon
dynMat$bin <- dynMat$bin + 1 # Add 1 to all values of taxon

times = unique(dynMat$mid) # Time Vector

##### #### #####
# Correct Bins #
##### #### #####
tsdat <- data.frame(
  bin    = seq_len(length(times) - 1),  # Bin ID
  bottom = times[-length(times)],         # Start of bin
  top    = times[-1]                    # End of bin
)

##### #### ##### ###
# Fad-Lad matrix   #
##### #### ##### ###
fl <- fadlad(dynMat, bin="mid", tax="taxon")

########
# Plot #
########
tsplot(tsdat, bottom = "bottom", top = "top", xlim = (range(tsdat$bottom, tsdat$top)))
ranges(fl, tax="taxon", labs=F, labels.args=list(cex=0.2), occs=T)

# Get some stats
samp <-binstat(dynMat, tax="taxon", bin="bin", indices=TRUE)
colnames(samp)

# basic plot of occurrences per bin
oldPar <- par(mar=c(4,4,2,4))
var_p = samp$occs
tsplot(tsdat, bottom = "bottom", top = "top", xlim = (range(tsdat$bottom, tsdat$top)), ylim = c(0,max(var_p)))
lines(times, var_p, col="blue")
par(oldPar)

#################
# Raw Diversity #
#################
# Raw Diversity
dd <- divDyn(dynMat, tax = "taxon", bin = "bin")

# Plot of Diversity Metrics
tsplot(tsdat, bottom = "bottom", top = "top", xlim = (range(tsdat$bottom, tsdat$top)),
       ylab = "Richness", xlab = "Time", ylim = c(0,30))
# lines
lines(times, dd$divRT, col= "red", lwd=2)
lines(times, dd$divBC, col="blue", lwd=2)
lines(times, dd$divSIB, col="green", lwd=2)
# legend
legend("topleft", legend=c("RT", "BC", "SIB"),
       col=c("red", "blue", "green"), lwd=c(2,2,2), bg="white")

# Turnover Metrics
tsplot(tsdat, bottom = "bottom", top = "top", xlim = (range(tsdat$bottom, tsdat$top)),
       ylab = "Richness", xlab = "Time", ylim = c(0,2))
lines(times, dd$extPC, col="black", lwd=2)
lines(times, dd$extGF, col="blue", lwd=2)
lines(times, dd$ext2f3, col="green", lwd=2)
# legend
legend("topright", legend=c("per capita", "gap-filler", "second-for-third"),
       col=c("black", "blue", "green"), lwd=c(2,2,2), bg="white")


# Test for sampling completeness
cor.test(dd$divRT, samp$occs, method="spearman")
# With perfect sampling our diversity and occurrences are identical
# This should be expected because we observe all taxa!
#'''	
#Spearman's rank correlation rho
#
#data:  dd$divRT and samp$occs
#S = 0, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#1 
#'''

##################################
# Classific Rarefaction Sampling #
##################################
subRes <- subsample(
  dynMat, 
  tax = "taxon",   # column for taxa
  bin = "bin",     # column for bins
  q   = 10,       # quota: number of occurrences per sample
  iter = 10,      # number of iterations
  type = "cr"    # Classical Rarefaction
)

subRes <- merge(subRes, unique(dynMat[, c("bin", "mid")]), by = "bin", all.x = TRUE)

plot(subRes$mid, subRes$divRT, type = "b", 
     xlab = "Time (Generations)", ylab = "Subsampled richness", ylim = c(0,35))
