setwd('/home/laura/Projects/ovrf/temp/oct_hiv')
#setwd('/home/laura/Downloads')
library(tidyr)
no_ov <- read.csv('just_env_hyphy.csv')

ov <- read.csv('env_oct27_1.csv')
no_ov<-read.csv('hyphysim_8mu_iav_24_0.2_2295.csv')
new_ov<- ov %>% separate(Selection.detected., c("neg-pos", "p-value"), sep='. p = ',
                convert = TRUE, remove = FALSE)

new_no_ov<- no_ov %>% separate(Selection.detected., c("neg-pos", "p-value"), sep='. p = ',convert = TRUE, remove = FALSE)
# Create new column filled with default colour
new_ov$Colour="black"
# Set new column values to appropriate colours
new_ov$Colour[new_ov$`neg-pos`=="Neg"]="#143D59"
new_ov$Colour[new_ov$`neg-pos`=="Pos"]="#F4B41A"

# Create new column filled with default colour
new_no_ov$Colour="black"
# Set new column values to appropriate colours
new_no_ov$Colour[new_no_ov$`neg-pos`=="Neg"]="#143D59"
new_no_ov$Colour[new_no_ov$`neg-pos`=="Pos"]="#F4B41A"


#Create plots
par(mfrow=c(3,2))
par(mar=c(4, 4, 2, 1))
plot(no_ov$alpha, main = "alpha values, global rate = 0.2", xlab = "")
plot(ov$Codon, ov$alpha, main = "alpha values, global rate = 1", xlab = "")
plot(no_ov$beta, main = "beta values", xlab = "")
plot(ov$Codon, ov$beta, main = "beta values", xlab = "")
plot(new_no_ov$`p-value`, col = new_no_ov$Colour, main = "p_values", xlab = "codon position", pch = 19)
legend("topright", legend =c("Neg", "Pos"), col = c("#143D59", "#F4B41A"), pch = 19, cex=0.8)
plot(new_ov$Codon, new_ov$`p-value`, col = new_ov$Colour, main = "p_values", xlab = "codon position", pch = 19)
legend("topright", legend =c("Neg", "Pos"), col = c("#143D59", "#F4B41A"), pch = 19, cex=0.8)
summary(new_ov)
table(new_ov$`neg-pos`)
