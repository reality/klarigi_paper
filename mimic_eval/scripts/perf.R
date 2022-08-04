
multivariate_sim_class <- read.delim("~/Projects/klarigi_paper/mimic_eval/data/sim/multivariate_sim_class.tsv")
enrichment_sim_class <- read.delim("~/Projects/klarigi_paper/mimic_eval/data/sim/enrichment_sim_class.tsv")
univariate_sim_class <- read.delim("~/Projects/klarigi_paper/mimic_eval/data/sim/univariate_sim_class.tsv")

multivariate_sim_class$match1 <- as.factor(multivariate_sim_class$match1)
enrichment_sim_class$match1 <- as.factor(enrichment_sim_class$match1)
univariate_sim_class$match1 <- as.factor(univariate_sim_class$match1)

multivariate_sim_class$match2 <- as.factor(multivariate_sim_class$match2)
enrichment_sim_class$match2 <- as.factor(enrichment_sim_class$match2)
univariate_sim_class$match2 <- as.factor(univariate_sim_class$match2)

print("enrichment")
print(AUC::auc(AUC::roc(enrichment_sim_class$pesim, enrichment_sim_class$match1)))
print(AUC::auc(AUC::roc(enrichment_sim_class$pnsim, enrichment_sim_class$match2)))

print("multivariate")
print(AUC::auc(AUC::roc(multivariate_sim_class$pesim, multivariate_sim_class$match1)))
print(AUC::auc(AUC::roc(multivariate_sim_class$pnsim, multivariate_sim_class$match2)))

print("univariate")
print(AUC::auc(AUC::roc(univariate_sim_class$pesim, univariate_sim_class$match1)))
print(AUC::auc(AUC::roc(univariate_sim_class$pnsim, univariate_sim_class$match2)))
