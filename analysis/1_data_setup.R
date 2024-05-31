# PANDAnet 1 - Data setup ====

# Load packages ====
library(haven)
library(dplyr)
library(naniar)
library(mice)
library(corrplot)
library(readxl)
library(ggplot2)
library(tidyr)
library(readr)
library(networktools)

# Directory Organisation ====
raw.data.url = "data/PANDA trial data merged_clean_13.dta"

output = "output/"

results = "results/"

# Inspect data: save relevant vars, check missingness etc. ====

## Load raw data
raw.data = read_dta(raw.data.url)


## Save variable names and labels
labels = lapply(raw.data, function(x) attributes(x)$label)   

unlist.labels = data.frame("labels" = unlist(labels)) %>% tibble::rownames_to_column(., "variable name")

write.csv(unlist.labels, paste0(output, 'stata.labels.csv'), row.names = F, quote = F)

## Save relevant variables for networks (referring to PANDA documentation on data management & analysis in Stata, notes in stata.labels.inspected.xlsx)

included_variables = read_lines(paste0(output, 'included_variables'))
net.data = raw.data %>% select(all_of(included_variables))

## Exclusions

# Non randomised: check the IDs of those not randomised are not in file (based on analysis script in PANDA doc)
non.rand.ids = c(30007,40139,30022,20001,20155,20227,40071,40049)
non.rand.ids %in% net.data$identifier_n # 40049 is in there
net.data = net.data[-(which(net.data$identifier_n == 40049)),] # removing ID row
rand.n = length(net.data$identifier_n) # save n of randomised

# Missing data
miss.data = net.data %>%  miss_var_summary()
mean(miss.data$pct_miss)

# Mean age and sex

net.data %>% group_by(group) %>%
  summarise(mean = mean(age, na.rm = T),
            sd = sd(age, na.rm = T),
            n = n())

# Reduce number of nodes/ collinearity with goldbricker function ====

baseline = net.data %>% select(identifier_n,
                               group,
                               phq1:phq9,
                               gad1,
                               gad2:gad7,
                               beck1:beck21,
                               sf_physical,
                               global)

gold <- goldbricker(baseline[-c(1,2)], threshold = 0.4, corMin = 0.5, p = 0.05)
reductions = as.data.frame(gold$suggested_reductions)


write.csv(reductions, paste0(output, 'reductions_gold.csv'))

# Recombine items/nodes ====
# Function to recombine data based on goldbricker results in reductions_gold_inspected.csv

recomb = function(data, fu) {
  # tired: beck20,beck15
  data$r_tired = data %>% select(paste(c("beck20", "beck15", "phq4"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # worry: gad3, gad2
  data$r_worry = data %>% select(paste(c("gad2", "gad3"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # anhedonia: beck12, beck4, phq1
  data$r_anh = data %>% select(paste(c("beck12", "beck4", "phq1"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # dislike: beck14, beck7, beck8, beck2
  data$r_dislike = data %>% select(paste(c("beck14", "beck7", "beck8", "beck2"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # restless: gad5, gad6, beck11, gad4, phq8, beck17
  data$r_restless = data %>% select(paste(c("gad5", "gad6", "beck11", "gad4", "phq8", "beck17"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # suicidal: beck9, phq9
  data$r_suicidal = data %>% select(paste(c("beck9", "phq9"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # sad: beck1, phq2
  data$r_sad = data %>% select(paste(c("beck1", "phq2"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # appetite: beck18, phq5
  data$r_appetite = data %>% select(paste(c("beck18", "phq5"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # concentration: beck19, phq7
  data$r_concentration = data %>% select(paste(c("beck19", "phq7"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  # sleep: phq3, beck16
  data$r_sleep = data %>% select(paste(c("beck16", "phq3"), fu, sep = "")) %>% rowMeans(., na.rm = T) %>% ceiling()
  
  to.drop = paste(c("beck20", "beck15", "gad3", "gad2", "beck12", "beck4", "phq1", 
                    "beck14", "beck7", "beck8", "beck2", "gad5", "gad6", "beck11", "gad4",
                    "phq8", "beck17", "beck9", "phq9", "beck1", "phq2", "beck18", "phq5", 
                    "beck19", "phq3", "beck16", "phq4", "phq7"),
                  fu, sep = "")
  
  data.ref= data %>% select(-all_of(to.drop))
  data.ref = data.ref %>% rename_with(~ paste(., fu, sep = ""), r_tired:r_sleep)
  
  return(data.ref)
}

# Node descriptions and abbreviations ====
node.names = c("group" = "TREAT",
               "phq6" = "BAD",
               "gad1" = "ANX",
               "gad7" = "AFR",
               "beck3" = "FAI",
               "beck5" = "GUI",
               "beck6" = "PUN",
               "beck10" = "CRY",
               "beck13" = "IND",
               "beck21" = "LIB",
               "sf_physical" = "PHY",
               "global" = "IMP",
               "r_tired" = "TIR",
               "r_worry" = "WOR",
               "r_anh" = "ANH",
               "r_dislike" = "DIS",
               "r_restless" = "RES",
               "r_suicidal" = "SUI",
               "r_sad" = "SAD",
               "r_appetite" = "APP",
               "r_concentration" = "CON",
               "r_sleep" = "SLE")
node.desc = c("TREAT" = "Treatment",
              "BAD" = "Feeling bad \nabout oneself",
              "ANX" = "Feeling nervous \nor anxious",
              "AFR" = "Feeling afraid",
              "FAI" = "Past failure",
              "GUI" = "Guilty feelings",
              "PUN" = "Punishment \nfeelings",
              "CRY" = "Crying",
              "IND" = "Indecisiveness",
              "LIB" = "Loss of interest \nin sex",
              "PHY" = "General physical \nhealth",
              "IMP" = "Self-reported \nimprovement",
              "TIR" = "Feeling tired",
              "WOR" = "Feeling worried",
              "ANH" = "Anhedonia",
              "DIS" = "Disliking oneself",
              "RES" = "Being restless \nor slow",
              "SUI" = "Suicidal thoughts",
              "SAD" = "Feeling sad \nor depressed",
              "APP" = "Lack of appetite or \neating too much",
              "CON" = "Concentration \nproblems",
              "SLE" = "Sleep problems")

# Recombine data ====

base.ref = recomb(baseline, "") 

twoweeks = net.data %>% select(identifier_n,
                               group,
                               phq1_2wk:phq9_2wk,
                               gad1_2wk:gad7_2wk,
                               beck1_2wk:beck21_2wk,
                               sf_physical_2wk,
                               global_2wk)

twoweeks.ref = recomb(twoweeks, "_2wk")


sixweeks = net.data %>% select(identifier_n,
                               group,
                               phq1_6wk:phq9_6wk,
                               gad1_6wk:gad7_6wk,
                               beck1_6wk:beck21_6wk,
                               sf_physical_6wk,
                               global_6wk)

sixweeks.ref = recomb(sixweeks, "_6wk")

twelveweeks = net.data %>% select(identifier_n,
                                  group,
                                  phq1_12wk:phq9_12wk,
                                  gad1_12wk:gad7_12wk,
                                  beck1_12wk:beck21_12wk,
                                  sf_physical_12wk,
                                  global_12wk)

twelveweeks.ref = recomb(twelveweeks, "_12wk")

# Prepare df for longitudinal nets  ====

data.ref = merge(base.ref, twoweeks.ref, by = c("identifier_n", "group")) %>%
  merge(., sixweeks.ref, by =  c("identifier_n", "group")) %>%
  merge(., twelveweeks.ref, by =  c("identifier_n", "group"))

write.csv(data.ref, paste0(output, 'data.ref.csv'), row.names = F, quote = F)


# adding other relevant variables
data.ref.full = net.data %>% select(identifier_n,
                                    group,
                                    eligible_rand,
                                    age,
                                    agecat,
                                    empstat2,
                                    fin3,
                                    edu3,
                                    lifeevents,
                                    ethnic,
                                    ethnic2,
                                    `_site_n`,
                                    severity_rand,
                                    antidepressantsinpast,
                                    depdiag,
                                    depr,
                                    depr_dur,
                                    depressedinpast,
                                    highestqualifications,
                                    married,
                                    marstat3,
                                    sex,
                                    withdrew,
                                    dateofwithdrawal,
                                    phqtot,
                                    gadtot,
                                    cisdepscore,
                                    cistotal_cat,
                                    norm_sf_physical) %>%
  merge(., data.ref, by =  c("identifier_n", "group"))

write.csv(data.ref.full, paste0(output, 'data.ref.full.csv'), row.names = F, quote = F)


