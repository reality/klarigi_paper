
multivariate_sim_class <- read.delim("data/sim/multivariate_sim_class.tsv")
enrichment_sim_class <- read.delim("data/sim/enrichment_sim_class.tsv")
univariate_sim_class <- read.delim("data/sim/univariate_sim_class.tsv")

multivariate_sim_class$match1 <- as.factor(multivariate_sim_class$match1)
enrichment_sim_class$match1 <- as.factor(enrichment_sim_class$match1)
univariate_sim_class$match1 <- as.factor(univariate_sim_class$match1)

multivariate_sim_class$match2 <- as.factor(multivariate_sim_class$match2)
enrichment_sim_class$match2 <- as.factor(enrichment_sim_class$match2)
univariate_sim_class$match2 <- as.factor(univariate_sim_class$match2)

print("enrichment")
print(paste("  ", "pulmonary embolism: ", signif(as.numeric(AUC::auc(AUC::roc(enrichment_sim_class$pesim, enrichment_sim_class$match1))), 3)))
pROC::ci(enrichment_sim_class$match1, enrichment_sim_class$pesim)
print(paste("  ", "pneumonia: ", signif(as.numeric(AUC::auc(AUC::roc(enrichment_sim_class$pnsim, enrichment_sim_class$match2))), 3)))
pROC::ci(enrichment_sim_class$match1, enrichment_sim_class$pnsim)

print("multivariate")
print(paste("  ", "pulmonary embolism: ", signif(as.numeric(AUC::auc(AUC::roc(multivariate_sim_class$pesim, multivariate_sim_class$match1))), 3)))
pROC::ci(multivariate_sim_class$match1, multivariate_sim_class$pesim)
print(paste("  ", "pneumonia: ", signif(as.numeric(AUC::auc(AUC::roc(multivariate_sim_class$pnsim, multivariate_sim_class$match2))), 3)))
pROC::ci(multivariate_sim_class$match2, multivariate_sim_class$pnsim)

print("univariate")
print(paste("  ", "pulmonary embolism: ", signif(as.numeric(AUC::auc(AUC::roc(univariate_sim_class$pesim, univariate_sim_class$match1))), 3)))
pROC::ci(univariate_sim_class$match1, univariate_sim_class$pesim)
print(paste("  ", "pneumonia: ", signif(as.numeric(AUC::auc(AUC::roc(univariate_sim_class$pnsim, univariate_sim_class$match2))), 3)))
pROC::ci(univariate_sim_class$match2, univariate_sim_class$pnsim)


enrichment_sim_class <- read.delim("reclassify-pulmonary embolism-scores.lst",header=F)
pROC::auc(enrichment_sim_class$V3, enrichment_sim_class$V2)
pROC::ci(enrichment_sim_class$V3, enrichment_sim_class$V2)
