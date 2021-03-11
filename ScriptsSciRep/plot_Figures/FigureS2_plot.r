# Simulation of the mixed distribution, Figure S2

library(viridis)

COL1 = viridis(10, alpha=0.7)
COL2 = viridis(10, alpha=0.7)


x <- seq(5,15,0.01)
subs1 <-  seq(9.5,10.5, length.out=5) # rnorm(5,10,.25)
subs2 <-  seq(8.66359,11.33641, length.out=5)


msubs1 <- mean(subs1)
msubs2 <- mean(subs2)
vsubs1 <- var(subs1)
vsubs2 <- var(subs2)

dsubs1 <- lapply(subs1, function(z) dnorm(x, z, .2 ) )
dsubs1_2 <- lapply(subs1, function(z) dnorm(x, z, 1 ) )
dsubs2 <- lapply(subs2, function(z) dnorm(x, z, .2 ) )

ind1 <- dnorm(x, msubs1, sqrt(vsubs1+.2^2) )
ind1_2 <- dnorm(x, msubs1, sqrt(vsubs1+1^2) )
ind2 <- dnorm(x, msubs2, sqrt(vsubs2+.2^2) )


setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("FigureS2.pdf", paper="a4r", width=9, height=6)

par(mfrow=c(2,3))
plot(x, dsubs1[[1]], type="l", col=COL2[3], lwd=2, las=1, xlab="log2 expression (gene ensemble)", ylab="density", las=1)
for(i in 2:length(subs1)) { lines(x, dsubs1[[i]], col=COL2[i+2], lwd=2) }
segments(x0 = subs1, y0=1/(.2*sqrt(2*pi)), x1=subs1, y1= 1/(.2*sqrt(2*pi)) + .07,  col=COL2[3:10], lwd=1)
plot(x, dsubs1_2[[1]], type="l", col=COL2[3], lwd=2, las=1, xlab="log2 expression (gene ensemble)", ylab="density", las=1)
for(i in 2:length(subs1)) { lines(x, dsubs1_2[[i]], col=COL2[i+2], lwd=2) }
segments(x0 = subs1, y0=1/(1*sqrt(2*pi)), x1=subs1, y1=1/(1*sqrt(2*pi))+.03,  col=COL2[3:10], lwd=1)
plot(x, dsubs2[[1]], type="l", col=COL2[3], lwd=2, las=1, xlab="log2 expression (gene ensemble)", ylab="density", las=1)
for(i in 2:length(subs2)) { lines(x, dsubs2[[i]], col=COL2[i+2], lwd=2) }
segments(x0 = subs2, y0=1/(.2*sqrt(2*pi)), x1=subs2, y1= 1/(.2*sqrt(2*pi))+.07,  col=COL2[3:10], lwd=1)

plot(x, ind1, col=COL2[1], lwd=2, type="l", ylim=c(0,0.9), xlab="log2 expression (compound)", ylab="density", las=1)
plot(x, ind1_2, col=COL2[1], lwd=2, type="l", ylim=c(0,0.9), xlab="log2 expression", ylab="density", las=1)
plot(x, ind2, col=COL2[1], lwd=2, type="l", ylim=c(0,0.9), xlab="log2 expression", ylab="density", las=1)

dev.off()
