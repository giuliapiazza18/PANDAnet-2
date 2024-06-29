# PANDAnet 4 - Estimating temporally lagged networks ====

# Load packages ====
library(dplyr)
library(lavaan)
library(tidyr)
library(qgraph)

# Aggregate group networks ====

# Specify model with time invariant variable ====
ti.model.full = "
                                                                                   
BAD_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
ANX_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
AFR_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
FAI_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
GUI_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
PUN_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
CRY_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
IND_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
LIB_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
PHY_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
IMP_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
TIR_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
WOR_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
ANH_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
DIS_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
RES_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
SUI_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
SAD_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
APP_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
CON_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT
SLE_6wk  ~  BAD_2wk+ANX_2wk+AFR_2wk+FAI_2wk+GUI_2wk+PUN_2wk+CRY_2wk+IND_2wk+LIB_2wk+PHY_2wk+IMP_2wk+TIR_2wk+WOR_2wk+ANH_2wk+DIS_2wk+RES_2wk+SUI_2wk+SAD_2wk+APP_2wk+CON_2wk+SLE_2wk + TREAT

BAD_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
ANX_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
AFR_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
FAI_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
GUI_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
PUN_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
CRY_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
IND_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
LIB_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
PHY_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
IMP_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
TIR_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
WOR_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
ANH_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
DIS_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
RES_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
SUI_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
SAD_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
APP_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
CON_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
SLE_12wk  ~  BAD_6wk+ANX_6wk+AFR_6wk+FAI_6wk+GUI_6wk+PUN_6wk+CRY_6wk+IND_6wk+LIB_6wk+PHY_6wk+IMP_6wk+TIR_6wk+WOR_6wk+ANH_6wk+DIS_6wk+RES_6wk+SUI_6wk+SAD_6wk+APP_6wk+CON_6wk+SLE_6wk + TREAT
"
# Run model ====
ti.model.full.sem = sem(ti.model.full, data = data.res, estimator = "ML", missing = "FIML")



# Save parameters ====
ti.estimates.full = subset( # save all estimated parameters for t2 nodes
  parameterEstimates(ti.model.full.sem, standardized = TRUE),
  op == "~"
)

# Save parameters in each lag ====
ti.estimates.lag1 = ti.estimates.full %>% filter(., grepl("_6wk", lhs))
ti.estimates.lag2 = ti.estimates.full %>% filter(., grepl("_12wk", lhs))

# Create list of edge weights for networks and plot
ti.net.lag1 = ti.estimates.lag1 %>% 
  mutate(est.net = ifelse(pvalue <0.05,est,0),
         t2.node = gsub("_6wk", "", lhs),
         t1.node = gsub("_2wk", "", rhs)) %>%
  select(t1.node,t2.node,est.net)

ti.net.lag1 %>% filter(t1.node=="TREAT" & abs(est.net) >0)

ti.net.lag1.plot = qgraph(ti.net.lag1)


ti.net.lag1.st = ti.estimates.lag1 %>% 
  mutate(est.net = ifelse(pvalue <0.05,std.all,0),
         t2.node = gsub("_6wk", "", lhs),
         t1.node = gsub("_2wk", "", rhs)) %>%
  select(t1.node,t2.node,est.net)

ti.net.lag1.st %>% filter(t1.node=="TREAT" & abs(est.net) >0)

ti.net.lag2 = ti.estimates.lag2 %>% 
  mutate(est.net = ifelse(pvalue <0.05,est,0),
         t2.node = gsub("_12wk", "", lhs),
         t1.node = gsub("_6wk", "", rhs)) %>%
  select(t1.node,t2.node,est.net)

ti.net.lag2 %>% filter(t1.node=="TREAT" & abs(est.net) >0)

ti.net.lag2.plot = qgraph(ti.net.lag2)

ti.net.lag2.st = ti.estimates.lag2 %>% 
  mutate(est.net = ifelse(pvalue <0.05,std.all,0),
         t2.node = gsub("_12wk", "", lhs),
         t1.node = gsub("_6wk", "", rhs)) %>%
  select(t1.node,t2.node,est.net)

# Separate group networks ====
# M1: Free groups, free time points ====
# two lags: 6wk predicted by 2wk, 12wk predicted by 6wk
n.nodes = 21

model_free_lag1<- rep(0, n.nodes) # create empty variables for each lag
model_free_lag2<- rep(0, n.nodes)


for (i in 1:n.nodes){  # for each network node
  name.t2 <- names(data.res[-1])[n.nodes+i] # save the name of the node at t
  name.t3 <- names(data.res[-1])[(2*n.nodes)+i] 
  include.t1 <- names(data.res[-1])[1:n.nodes] # save the (t - 1) nodes to include
  include.t2 <- names(data.res[-1])[(n.nodes+1):(n.nodes+21)]
  model_free_lag1[i] <- paste(name.t2, " ~ ", paste(paste(include.t1, sep = ""), collapse = "+")) # create model per lag
  model_free_lag2[i] <- paste(name.t3, " ~ ", paste(paste(include.t2, sep = ""), collapse = "+"))
}  

model_free_full <- c(model_free_lag1, model_free_lag2) # create full model with all lags

freegroups_freetime = sem(model_free_full, data = data.res, group = "TREAT", estimator = "ML", missing = "FIML")
summary_freegroups = summary(freegroups_freetime)

# M2: Equal groups, free time points ====

a = c(1:n.nodes)
b = letters[1:21]
c = paste0(letters,letters)[1:21]

model_equal_lag1<- rep(0, n.nodes)
model_equal_lag2<- rep(0, n.nodes)

for (i in 1:n.nodes){  
  name.t2 <- names(data.res[-1])[n.nodes+i] 
  name.t3 <- names(data.res[-1])[(2*n.nodes)+i] 
  include.t1 <- names(data.res[-1])[1:n.nodes]
  include.t2 <- names(data.res[-1])[(n.nodes+1):(n.nodes+21)]
  betas.t1 <- paste0(b, a[i]) # save name of coefficients
  betas.t2 <- paste0(c, a[i])
  model_equal_lag1[i] <- paste(name.t2, " ~ ", paste(paste(betas.t1,'*',include.t1, sep = ""), collapse = "+"))
  model_equal_lag2[i] <- paste(name.t3, " ~ ", paste(paste(betas.t2,'*',include.t2, sep = ""), collapse = "+"))
}  

model_equal_full <- c(model_equal_lag1, model_equal_lag2)

equalgroups_freetime = sem(model_equal_full, data = data.res, group = "TREAT", group.equal = c("regressions"), estimator = "ML", missing = "FIML")

summary(equalgroups_freetime)


# Compare models ====

model.comparison = anova(freegroups_freetime,
                         equalgroups_freetime)

write.csv(model.comparison, paste0(results, "tables/model.comparison.csv"), quote = F)

# Save parameters =====
m1_par = subset( # save all estimated regressions
  parameterEstimates(freegroups_freetime),
  op == "~" & lhs %in% c(colnames(data.res))
)

m2_par = subset(
  parameterEstimates(equalgroups_freetime),
  op == "~" & lhs %in% c(colnames(data.res))
)


