# PANDAnet 2 - Linear Mixed Models ====

# Load packages 
library(dplyr)
library(nlme)
library(lmerTest)
library(effectsize)

# Rename data frame  ====
data.node = data.ref.full %>%
  rename(., "TREAT" = "group",
         "BAD_base" = "phq6",
         "ANX_base" = "gad1",
         "AFR_base" = "gad7",
         "FAI_base" = "beck3",
         "GUI_base" = "beck5",
         "PUN_base" = "beck6",
         "CRY_base" = "beck10",
         "IND_base" = "beck13",
         "LIB_base" = "beck21",
         "PHY_base" = "sf_physical",
         "IMP_base" = "global",
         "TIR_base" = "r_tired",
         "WOR_base" = "r_worry",
         "ANH_base" = "r_anh",
         "DIS_base" ="r_dislike",
         "RES_base" = "r_restless",
         "SUI_base" = "r_suicidal",
         "SAD_base" = "r_sad",
         "APP_base" = "r_appetite",
         "CON_base" = "r_concentration",
         "SLE_base" = "r_sleep",
         "BAD_2wk" = "phq6_2wk",
         "ANX_2wk" = "gad1_2wk",
         "AFR_2wk" = "gad7_2wk",
         "FAI_2wk" = "beck3_2wk",
         "GUI_2wk" = "beck5_2wk",
         "PUN_2wk" = "beck6_2wk",
         "CRY_2wk" = "beck10_2wk",
         "IND_2wk" = "beck13_2wk",
         "LIB_2wk" = "beck21_2wk",
         "PHY_2wk" = "sf_physical_2wk",
         "IMP_2wk" = "global_2wk",
         "TIR_2wk" = "r_tired_2wk",
         "WOR_2wk" = "r_worry_2wk",
         "ANH_2wk" = "r_anh_2wk",
         "DIS_2wk" ="r_dislike_2wk",
         "RES_2wk" = "r_restless_2wk",
         "SUI_2wk" = "r_suicidal_2wk",
         "SAD_2wk" = "r_sad_2wk",
         "APP_2wk" = "r_appetite_2wk",
         "CON_2wk" = "r_concentration_2wk",
         "SLE_2wk" = "r_sleep_2wk",
         "BAD_6wk" = "phq6_6wk",
         "ANX_6wk" = "gad1_6wk",
         "AFR_6wk" = "gad7_6wk",
         "FAI_6wk" = "beck3_6wk",
         "GUI_6wk" = "beck5_6wk",
         "PUN_6wk" = "beck6_6wk",
         "CRY_6wk" = "beck10_6wk",
         "IND_6wk" = "beck13_6wk",
         "LIB_6wk" = "beck21_6wk",
         "PHY_6wk" = "sf_physical_6wk",
         "IMP_6wk" = "global_6wk",
         "TIR_6wk" = "r_tired_6wk",
         "WOR_6wk" = "r_worry_6wk",
         "ANH_6wk" = "r_anh_6wk",
         "DIS_6wk" ="r_dislike_6wk",
         "RES_6wk" = "r_restless_6wk",
         "SUI_6wk" = "r_suicidal_6wk",
         "SAD_6wk" = "r_sad_6wk",
         "APP_6wk" = "r_appetite_6wk",
         "CON_6wk" = "r_concentration_6wk",
         "SLE_6wk" = "r_sleep_6wk",
         "BAD_12wk" = "phq6_12wk",
         "ANX_12wk" = "gad1_12wk",
         "AFR_12wk" = "gad7_12wk",
         "FAI_12wk" = "beck3_12wk",
         "GUI_12wk" = "beck5_12wk",
         "PUN_12wk" = "beck6_12wk",
         "CRY_12wk" = "beck10_12wk",
         "IND_12wk" = "beck13_12wk",
         "LIB_12wk" = "beck21_12wk",
         "PHY_12wk" = "sf_physical_12wk",
         "IMP_12wk" = "global_12wk",
         "TIR_12wk" = "r_tired_12wk",
         "WOR_12wk" = "r_worry_12wk",
         "ANH_12wk" = "r_anh_12wk",
         "DIS_12wk" ="r_dislike_12wk",
         "RES_12wk" = "r_restless_12wk",
         "SUI_12wk" = "r_suicidal_12wk",
         "SAD_12wk" = "r_sad_12wk",
         "APP_12wk" = "r_appetite_12wk",
         "CON_12wk" = "r_concentration_12wk",
         "SLE_12wk" = "r_sleep_12wk",
         "site" = `_site_n`) %>%
  apply(., 2, as.numeric) %>% 
  as.data.frame() 



# Wide to long ====
data.node$TREAT_B = data.node$TREAT - 1
data.mlm = reshape(data.node, 
                   direction = "long", 
                   varying = c(colnames(data.node[-c(1:50,114)])),
                   idvar = "identifier_n",
                   sep = "_")

# Recode time to be -10,-6,0 for unit to be the week
data.mlm.rec = data.mlm %>% 
  mutate_at(vars(time), funs(recode(., `2wk`=-10, `6wk`=-6, `12wk`=0)))

# LMMs ====

# Linear mixed models 
## Create function with lmer model for lmm 
mlm.symptom = function(x) {
  base.value = paste0(x, "_base")
  mlm.formula = paste0(x, "~  severity_rand + 
                       depr_dur + 
                       factor(site) + 
                       factor(TREAT_B) + 
                       time + 
                       time:factor(TREAT_B) + 
                       (0 + time|identifier_n) +", 
                       base.value)
  mlm1 = lmer(mlm.formula,
              data = data.mlm.rec,
              na.action = na.exclude,
              REML = T)
  return(mlm1)
}



## Apply function to all symptoms
mlm.models = lapply(node.names[2:22], mlm.symptom)
names(mlm.models) = node.names[2:22]

## Calculate p-values with lmerTest
mlm.anova = lapply(mlm.models, anova)

## Calculate effect sizes (eta sq) and CIs
mlm.eta = lapply(mlm.anova, eta_squared, alternative = "two.sided")


# Create table with all effects

get_summary_anova = function(x) {
  dataframe = data.frame(
    Symptom = names(mlm.anova[x]),
    Effect = rownames(mlm.anova[[x]]),
    Fvalue = round(mlm.anova[[x]]$`F value`,2),
    Df = paste0(mlm.anova[[x]]$`NumDF`, ", ", round(mlm.anova[[x]]$`DenDF`,1)),
    pvalue = round(mlm.anova[[x]]$`Pr(>F)`,5)
  )
  return(dataframe)
}

summary_table_anova = bind_rows(lapply(node.names[2:22], get_summary_anova)) %>%
  mutate_at(vars(Effect), 
            funs(recode(., `time`= "Time", 
                        `factor(TREAT_B)`= "Group", 
                        `factor(TREAT_B):time`= "Group x Time",
                        `severity_rand` = "Depression Severity",
                        `depr_dur` = "Depression Duration",
                        `factor(site)` = "Site")))

# Sample size in each model ====

get_n = function(x) {
  count.n = data.mlm %>% 
    filter_at(vars(x, time, TREAT_B, severity_rand, depr_dur, site), all_vars(!is.na(.))) %>% 
    group_by(time, TREAT_B) %>% 
    summarise(x = n())
  names(count.n) = c("time", "group", paste0(x))
  return(count.n)
}

sample.mlm = lapply(node.names[-1], get_n)

sample.mlm.grouped =merge(sample.mlm[[1]], sample.mlm[[2]], by = c("time", "group")) %>%
  merge(., sample.mlm[[3]], by = c("time", "group")) %>%
  merge(., sample.mlm[[4]], by = c("time", "group")) %>%
  merge(., sample.mlm[[5]], by = c("time", "group")) %>%
  merge(., sample.mlm[[6]], by = c("time", "group")) %>%
  merge(., sample.mlm[[7]], by = c("time", "group")) %>%
  merge(., sample.mlm[[8]], by = c("time", "group")) %>%
  merge(., sample.mlm[[9]], by = c("time", "group")) %>%
  merge(., sample.mlm[[10]], by = c("time", "group")) %>%
  merge(., sample.mlm[[11]], by = c("time", "group")) %>%
  merge(., sample.mlm[[12]], by = c("time", "group")) %>%
  merge(., sample.mlm[[13]], by = c("time", "group")) %>%
  merge(., sample.mlm[[14]], by = c("time", "group")) %>%
  merge(., sample.mlm[[15]], by = c("time", "group")) %>%
  merge(., sample.mlm[[16]], by = c("time", "group")) %>%
  merge(., sample.mlm[[17]], by = c("time", "group")) %>%
  merge(., sample.mlm[[18]], by = c("time", "group")) %>%
  merge(., sample.mlm[[19]], by = c("time", "group")) %>%
  merge(., sample.mlm[[20]], by = c("time", "group")) %>%
  merge(., sample.mlm[[21]], by = c("time", "group"))
  

write.csv(sample.mlm.grouped, file = paste0(results, "/tables/sample.mlm.grouped.csv"), row.names = F, quote = F)

