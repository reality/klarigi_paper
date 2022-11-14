HP_0001635 <- read.delim("~/Projects/klarigi_paper/mimic_eval/data/scores_examples/HP_0001635.txt", header=FALSE)
hist(HP_0001635$V1, breaks=50, freq=F)
hist(HP_0001635$V2, breaks=20)
hist(HP_0001635$V2)

HP_0012735 <- read.delim("~/Projects/klarigi_paper/mimic_eval/data/scores_examples/HP_0012735.txt", header=FALSE)
hist(log(HP_0012735$V1))

HP_0031652 <- read.delim("~/Projects/klarigi_paper/mimic_eval/data/scores_examples/HP_0031652.txt", header=FALSE)
hist(HP_0031652$V3)

library(ggpubr)
ggqqplot(HP_0001635$V1)
