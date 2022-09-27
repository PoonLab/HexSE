setwd('/home/laura/Projects/ovrf/paper_example/')

library(scales)

one <- read.csv("1gene_HBV.m20.shannon.txt", sep="\t", header = FALSE)
colnames(one) <-  c("pos", "value")

comp <- read.csv("correct_splits_HBV.0.shannon.txt", sep="\t", header = FALSE)
colnames(comp) <-  c("pos", "value")


# Fit local polynomial regression
lo_one <- loess(one$value~one$pos, span=0.1)
lo_comp <- loess(comp$value~comp$pos, span=0.1)


# pdf(file = "NEWshannon_equal_ShapeAndClases.pdf",
#     width = 12, # The width of the plot in inches
#     height = 5) # The height of the plot in inches

# Points
plot(comp, xaxt='n', col=alpha("#117fa3", 0.3), pch=16, cex=0.8, 
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

non_coding_one <- one[-one$pos[1800:2600],]
non_coding_one <- one[-one$pos[1815:2454],]
non_coding_one[which.min(non_coding_one$value),]

# Complete genome information (lowest are vs the rest)

lowest_region <- comp[comp$pos[1374:2454],]
rest_region <- comp[-comp$pos[1374:2454],]

# Only one orf (coding vs non-coding)
non_coding_one <- one[-one$pos[1815:2454],]
coding_one <- one[one$pos[1815:2454],]
