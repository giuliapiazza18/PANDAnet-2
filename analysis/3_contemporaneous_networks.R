# PANDAnet 3 - Estimating contemporaneous networks ====

# Load packages ====
library(bootnet)
library(dplyr)
library(networktools)
library(qgraph)
library(corrplot)
library(NetworkComparisonTest)

# Regress out sex, age, stratification variables (severity, duration, site), baseline values and variables associated with missingness at baseline ====
res.symptom = function(x) {
  base.value = paste0(substr(x, 1, 3), "_base")
  lm.formula = paste0(x, "~  factor(cistotal_cat) + factor(depr_dur) + factor(site) + factor(sex) + factor(age)  + factor(ethnic2) +
                      factor(marstat3) + factor(fin3) + factor(antidepressantsinpast) + factor(lifeevents) +", base.value)
  lm1 = as.vector(rstandard(lm(lm.formula,
                               data = data.node,
                               na.action = na.exclude)))
  return(lm1)
}

lm.models = lapply(colnames(data.node[51:113]), res.symptom)
names(lm.models) = colnames(data.node[51:113])  

data.res = cbind("TREAT" = data.node$TREAT, as.data.frame(lm.models))

# Recode treatment to be 1 and 0
data.res$TREAT = data.res$TREAT - 1

# Aggregate group networks ====
## TWO WEEKS ====
set.seed(1)
twoweeks.mgm.res.data = data.res %>% 
  select(TREAT, BAD_2wk:SLE_2wk)


twoweeks.mgm.res <- estimateNetwork(twoweeks.mgm.res.data, default = "mgm",  
                                    type = c("c", rep("g", 21)), 
                                    level = c(2, rep(1, 21)),    
                                    criterion = "CV",             
                                    nFolds = 10,            
                                    order = 2,                      
                                    binarySign = TRUE)

two.res.graph = qgraph(twoweeks.mgm.res$graph, layout = 'spring', labels = colnames(twoweeks.mgm.res.data))
flow(two.res.graph,1, theme = "colorblind")


## SIX WEEKS ====
sixweeks.mgm.res.data = data.res %>% 
  select(TREAT, BAD_6wk:SLE_6wk)


sixweeks.mgm.res <- estimateNetwork(sixweeks.mgm.res.data, default = "mgm",  
                                    type = c("c", rep("g", 21)), 
                                    level = c(2, rep(1, 21)),    
                                    criterion = "CV",             
                                    nFolds = 10,            
                                    order = 2,                      
                                    binarySign = TRUE)

six.res.graph = qgraph(sixweeks.mgm.res$graph, layout = 'spring', labels = colnames(sixweeks.mgm.res.data))
flow(six.res.graph,1, theme = "colorblind")

## TWELVE WEEKS ====
twelveweeks.mgm.res.data = data.res %>% 
  select(TREAT, BAD_12wk:SLE_12wk)

twelveweeks.mgm.res <-  estimateNetwork(twelveweeks.mgm.res.data, default = "mgm",  
                                        type = c("c", rep("g", 21)), 
                                        level = c(2, rep(1, 21)),    
                                        criterion = "CV",             
                                        nFolds = 10,            
                                        order = 2,                      
                                        binarySign = TRUE)

twelve.res.graph = qgraph(twelveweeks.mgm.res$graph, layout = 'spring', labels = colnames(twelveweeks.mgm.res.data))
flow(twelve.res.graph,1, theme = "colorblind")

# Save aggregate group networks ====

flattenCorrMatrix <- function(cormat) {
  colnames(cormat) = rownames(cormat) = node.names
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

write.csv((flattenCorrMatrix(twoweeks.mgm.res$graph)), file = "results/tables/two.mgm.csv", quote=F, row.names =F)
write.csv((flattenCorrMatrix(sixweeks.mgm.res$graph)), file = "results/tables/six.mgm.csv", quote=F, row.names =F)
write.csv((flattenCorrMatrix(twelveweeks.mgm.res$graph)), file = "results/tables/twelve.mgm.csv", quote=F, row.names =F)

# Sample size in each network
sample.mgm = data.frame(
  network = c("2wk", "6wk", "12wk"),
  n = c(dim(na.omit(twoweeks.mgm.res.data))[1],
        dim(na.omit(sixweeks.mgm.res.data))[1],
        dim(na.omit(twelveweeks.mgm.res.data))[1])
)

# Separate group networks ====
# Comparing placebo and drug at each time point

## 2 WEEKS

twoweeks.ebic.treat.res = twoweeks.mgm.res.data %>%  filter(TREAT == 1) %>% select(-TREAT) %>%
  estimateNetwork(., default = "EBICglasso")

twoweeks.ebic.pla.res = twoweeks.mgm.res.data %>%  filter(TREAT == 0) %>% select(-TREAT) %>%
  estimateNetwork(., default = "EBICglasso")


twoweeks.ebic.nct = NCT(twoweeks.ebic.treat.res, twoweeks.ebic.pla.res, it = 100, 
                        paired = FALSE, test.edges = TRUE, 
                        edges = "all", 
                        p.adjust.methods = "fdr",
                        test.centrality = TRUE, centrality = "strength")

twoweeks.ebic.nct$nwinv.pval # check p value for invariance test

## 6 WEEKS

sixweeks.ebic.treat.res = sixweeks.mgm.res.data %>%  filter(TREAT == 1) %>% select(-TREAT) %>%
  estimateNetwork(., default = "EBICglasso")

sixweeks.ebic.pla.res = sixweeks.mgm.res.data %>%  filter(TREAT == 0) %>% select(-TREAT) %>%
  estimateNetwork(., default = "EBICglasso")

sixweeks.ebic.nct = NCT(sixweeks.ebic.treat.res, sixweeks.ebic.pla.res, it = 100, 
                        paired = FALSE, test.edges = TRUE, 
                        edges = "all", 
                        p.adjust.methods = "fdr",
                        test.centrality = TRUE, centrality = "strength")

sixweeks.ebic.nct$nwinv.pval # check p value for invariance test

## 12 WEEKS

twelveweeks.ebic.treat.res = twelveweeks.mgm.res.data %>%  filter(TREAT == 1) %>% select(-TREAT) %>%
  estimateNetwork(., default = "EBICglasso")

twelveweeks.ebic.pla.res = twelveweeks.mgm.res.data %>%  filter(TREAT == 0) %>% select(-TREAT) %>%
  estimateNetwork(., default = "EBICglasso")


twelveweeks.ebic.nct = NCT(twelveweeks.ebic.treat.res, twelveweeks.ebic.pla.res, it = 100, 
                           paired = FALSE, test.edges = TRUE, 
                           edges = "all", 
                           p.adjust.methods = "fdr",
                           test.centrality = TRUE, centrality = "strength")

twelveweeks.ebic.nct$nwinv.pval # check p value for invariance test


# Save separate group networks ====
flattenCorrMatrixEbic <- function(cormat) {
  colnames(cormat) = rownames(cormat) = node.names[-1]
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
ebic.network = function(placebo, treatment) {
  df = flattenCorrMatrixEbic(placebo)
  df$edge.treat = flattenCorrMatrixEbic(treatment)$cor
  colnames(df) = c("Node 1", "Node 2", "Placebo Edge", "Sertraline Edge")
  return(df)
}


write.csv(ebic.network(twoweeks.ebic.pla.res$graph, twoweeks.ebic.treat.res$graph), 
          file = "results/tables/ebic.two.csv", quote = F, row.names = F)

write.csv(ebic.network(sixweeks.ebic.pla.res$graph, sixweeks.ebic.treat.res$graph), 
          file = "results/tables/ebic.six.csv", quote = F, row.names = F)

write.csv(ebic.network(twelveweeks.ebic.pla.res$graph, twelveweeks.ebic.treat.res$graph), 
          file = "results/tables/ebic.twelve.csv", quote = F, row.names = F)




