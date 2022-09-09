#setwd('/home/laura/Projects/ovrf/temp/oct_hiv')
setwd('/home/laura/Projects/ovrf/')
library(tidyr)
ov <- read.csv('significant_HBV_july28.fa_1375_1840.fa.csv')
new_ov<- ov %>% separate(Selection.detected., c("neg-pos", "p-value"), sep='. p = ',
                         convert = TRUE, remove = FALSE)
# Create new column filled with default colour
new_ov$Colour="black"
# Set new column values to appropriate colours
new_ov$Colour[new_ov$`neg-pos`=="Neg"]="#143D59"
new_ov$Colour[new_ov$`neg-pos`=="Pos"]="#F4B41A"
new_ov$dnds=new_ov$alpha/new_ov$beta

par(mfrow=c(3,1))
par(mar=c(4, 4, 2, 1))
my_pch <-19
my_cex <- 1.2

#alpha
plot(ov$Codon, ov$alpha, ylim = c(0,4), cex=my_cex, pch=my_pch,
     main=paste0("Alpha values. Mean=", round(mean(ov$alpha),digits=2)), 
     xlab = "", col="darkorange", las=1)

#beta
plot(ov$Codon, ov$beta, ylim= c(0,4), cex=my_cex, pch=my_pch,
     main=paste0("Beta values. Mean=", round(mean(ov$beta),digits=2))
     , xlab = "", col="#a95aa1", las=1)

# p values
plot(new_ov$Codon, new_ov$`p-value`,
     col = new_ov$Colour, main = "p-values",
     xlab = "Codon position", pch = 19, las=1)

#alpha/beta = dN/dS
plot(new_ov$Codon, new_ov$dnds,
     col = new_ov$Colour, main = "dN/dS",
     xlab = "Codon position", pch = 19, las=1)

legend("topright", legend =c("Neg", "Pos"), col = c("#143D59", "#F4B41A"), pch = 19, cex=1.2)
summary(new_ov)

###############################################################
# Plot FEL from JSON 
##############################################################
setwd("/home/laura/Projects/detectOVRF")
require(jsonlite)
json <- read_json("cleaned_KY382410_hbvA.fa.FEL.json",
                  simplifyVector = TRUE)
summary(json$MLE)
#json$MLE$content[1]

# row indices correspond to codons
site.stats <- as.data.frame(json$MLE$content[[1]][,1:6])
names(site.stats) <- json$MLE$headers[,1]
site.stats$codon.pos<-row.names(site.stats)
site.stats$dnds<- site.stats$beta/site.stats$alpha # Calculate dN/dS ratio

par(mfrow=c(2,2))

## Paint with limits
# plot(site.stats$codon.pos, site.stats$dnds, main = "dN/dS",
#      xlab = "Codon position", las=1, col="#1e8449", pch = 19, cex=1, xlim=c(0,210), ylim=c(0,5))

# plot(site.stats$codon.pos, site.stats$beta, col="darkorange", pch = 19, cex=1, main="Individual measures", xlim=c(0,210), ylim=c(0,5))
# points(site.stats$codon.pos, site.stats$alpha, col="dodgerblue4", pch = 19, cex=1)
# legend("topleft", legend=c("Beta", "Alpha"), col=c("darkorange", "dodgerblue4"),
#        lty=1:1, cex=1.5, box.lty = 0, lwd=1)



######################################################
# Plot HBV Subtype A for simulated vs real data
######################################################
add.alpha <- function(col, alpha=1){
        if(missing(col))
                stop("Please provide a vector of colours.")
        apply(sapply(col, col2rgb)/255, 2, 
              function(x) 
                      rgb(x[1], x[2], x[3], alpha=alpha))  
}
pal <- c("#ae2012", "#bb3e03", "#ca6702", "#468faf", "#2c7da0", "#2a6f97", "#014f86", "#013a63", "#4f772d")
setwd("/home/laura/Projects/detectOVRF")
require(jsonlite)

# -------------------- Read and process SIMULATED data
sim_json <- read_json("cleaned_KY382410_hbvA.fa.FEL.json",
                  simplifyVector = TRUE)

# row indices correspond to codons
sim.stats <- as.data.frame(sim_json$MLE$content[[1]][,1:6])
names(sim.stats) <- sim_json$MLE$headers[,1]
sim.stats$codon.pos<-row.names(sim.stats)
sim.stats$dnds<- sim.stats$beta/sim.stats$alpha # Calculate dN/dS ratio


par(mfrow=c(2,2))
# Plot without limits
plot(sim.stats$codon.pos, sim.stats$dnds, main = "dN/dS Simulated data",
     xlab = "Codon position", las=1, col="#4f772d", pch = 19)


plot(sim.stats$codon.pos, sim.stats$beta, col=add.alpha("#ca6702",0.7), pch = 19, cex=1, main="simulated HBV")
points(sim.stats$codon.pos, sim.stats$alpha, col=add.alpha("#014f86",0.5), pch = 19, cex=1)
legend("topleft", legend=c("Beta", "Alpha"), col=c("darkorange", "dodgerblue4"),
       lty=1:1, cex=1.5, box.lty = 0, lwd=1)


# -------------------- Read and process REAL data
real_json <- read_json("HBV/A_HBV_hpclean.fasta.FEL.json",
                      simplifyVector = TRUE)

# row indices correspond to codons
real.stats <- as.data.frame(real_json$MLE$content[[1]][,1:6])
names(real.stats) <- real_json$MLE$headers[,1]
real.stats$codon.pos<-row.names(real.stats)
real.stats$dnds<- real.stats$beta/real.stats$alpha # Calculate dN/dS ratio

temp<-real.stats[which(real.stats$beta > 10),]
temp2<-sim.stats[which(sim.stats$beta > 4),]

real.stats <- real.stats[-which(real.stats$beta > 10),]


# Plot without limits
plot(real.stats$codon.pos, real.stats$dnds, main = "dN/dS Real data",
     xlab = "Codon position", las=1, col="#4f772d", pch = 19)

plot(real.stats$codon.pos, real.stats$beta, col=add.alpha("#ca6702",0.7), pch = 19, cex=1, main="real HBV subtype A")
points(real.stats$codon.pos, real.stats$alpha, col=add.alpha("#014f86",0.5), pch = 19, cex=1)
# legend("topleft", legend=c("Beta", "Alpha"), col=c("darkorange", "dodgerblue4"),
#        lty=1:1, cex=1.5, box.lty = 0, lwd=1)





######################################################
# Other
######################################################



## Single 
json <- read_json("cleaned_cds_one_1815-2454.fa.FEL.json",
                  simplifyVector = TRUE)

# row indices correspond to codons
site.stats <- as.data.frame(json$MLE$content[[1]][,1:6])
names(site.stats) <- json$MLE$headers[,1]
site.stats$codon.pos<-row.names(site.stats)
site.stats$dnds<- site.stats$beta/site.stats$alpha # Calculate dN/dS ratio

plot(site.stats$codon.pos, site.stats$dnds, main = "dN/dS on 1815-2454 region (Single ORFs)",
     xlab = "Codon position", las=1, col="#1e8449", pch = 19, cex=1, xlim=c(0,210), ylim=c(0,5))

plot(site.stats$codon.pos, site.stats$beta, col="darkorange", pch = 19, cex=1, main="Single ORFs", xlim=c(0,210), ylim=c(0,5))
points(site.stats$codon.pos, site.stats$alpha, col="dodgerblue4", pch = 19, cex=1)
legend("topleft", legend=c("Beta", "Alpha"), col=c("darkorange", "dodgerblue4"),
       lty=1:1, cex=1.5, box.lty = 0, lwd=1)

