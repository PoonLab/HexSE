setwd('/home/laura/Projects/ovrf/temp/oct_hiv')
#setwd('/home/laura/Downloads')
library(tidyr)
ov <- read.csv('env_oct27_1.csv')
new_ov<- ov %>% separate(Selection.detected., c("neg-pos", "p-value"), sep='. p = ',
                         convert = TRUE, remove = FALSE)
# Create new column filled with default colour
new_ov$Colour="black"
# Set new column values to appropriate colours
new_ov$Colour[new_ov$`neg-pos`=="Neg"]="#143D59"
new_ov$Colour[new_ov$`neg-pos`=="Pos"]="#F4B41A"

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

legend("topright", legend =c("Neg", "Pos"), col = c("#143D59", "#F4B41A"), pch = 19, cex=1.2)
summary(new_ov)

