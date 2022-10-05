setwd('/home/laura/Projects/ovrf/paper_example/')

library(scales)

one <- read.csv("1gen.HBV.m20.shannon.txt", sep="\t", header = FALSE)
colnames(one) <-  c("pos", "value")

comp <- read.csv("compl.HBV.m20.shannon.txt", sep="\t", header = FALSE)
colnames(comp) <-  c("pos", "value")


# Fit local polynomial regression
lo_one <- loess(one$value~one$pos, span=0.1)
lo_comp <- loess(comp$value~comp$pos, span=0.1)


# pdf(file = "shannon_HBV_compare.pdf",
#     width = 12, # The width of the plot in inches
#     height = 5) # The height of the plot in inches

# Points
plot(comp, xaxt='n', col=alpha("#117fa3", 0.3), pch=16, cex=0.8, ylim=c(0,1),
     ylab = "Shannon Entropy",
     xlab= "Genome position"
     )
points(one$pos, one$value, col=alpha("#f1a257", 0.3), pch=16, cex=0.8)
# Add trend lines
lines(predict(lo_comp), col="#005e7d", lwd=3)
lines(predict(lo_one), col="#f08723", lwd=3)
axis(side=1, at=seq(0, 3500, by=200))
legend(400, 0.2, legend=c("Multiple ORFs", "Single ORF"), col=c("#005e7d", "#f08723"),
       lty=1:1, cex=0.8, box.lty = 0, lwd=4)

# 
# dev.off()

# Only one orf in [1815:2454] (coding vs non-coding)
non_coding_one <- one[-one$pos[1815:2454],]
coding_one <- one[one$pos[1815:2454],]


# Complete genome information (lowest are vs the rest)



# Non coding one 

# Valid ORFs: [1375:1840, 2308:3182, 0: 1625, 2849:3182, 0:837, 1815:2454]
non_ov_c <- comp[comp$pos[1840:2308],] # Non overlapping on gene C
all_non_ov <- c(837:1365, 1625:1815, 1840:2308, 2545:2849) 
comp_non_ov <- comp[comp$pos[t_ov],]  # Total non_overlapping regions
comp_ov <- comp[-comp$pos[t_ov],]  # Total non_overlapping regions

# Comparing differences
comp$region <- ifelse(comp$pos %in% all_non_ov, "no-ov", "ov")
boxplot(value~region, data=comp)
res <- wilcox.test(value~region, data=comp)

# complete, overlap gene C
ov_c <- comp[comp$pos[2308:2454],]
