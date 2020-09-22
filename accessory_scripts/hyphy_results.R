setwd('/home/laura/Downloads')

no_ov <- read.csv('hyphy_no_overlaps.csv')
ov <- read.csv('with_overlaps_hiv.csv')

library(tidyr)
new_ov<- ov %>% separate(Selection.detected., c("neg-pos", "p-value"), sep='. p = ',
                convert = TRUE, remove = FALSE)

new_no_ov<- no_ov %>% separate(Selection.detected., c("neg-pos", "p-value"), sep='. p = ',
                               convert = TRUE, remove = FALSE)
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
plot(no_ov$alpha, main = "alpha values, CDS with no overlap", xlab = "")
plot(ov$alpha, main = "alpha values, CDS with overlap", xlab = "")
plot(no_ov$beta, main = "beta values, CDS with no overlap", xlab = "")
plot(ov$beta, main = "beta values, CDS with no overlap", xlab = "")
plot(new_no_ov$`p-value`, col = new_no_ov$Colour, main = "p_values", xlab = "codon position", pch = 19)
legend("topright", legend =c("Neg", "Pos"), col = c("#143D59", "#F4B41A"), pch = 19, cex=0.8)
plot(new_ov$`p-value`, col = new_ov$Colour, main = "p_values", xlab = "codon position", pch = 19)
legend("topright", legend =c("Neg", "Pos"), col = c("#143D59", "#F4B41A"), pch = 19, cex=0.8)
