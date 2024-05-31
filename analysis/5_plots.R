# PANDAnet 5 - Plotting graphs
library(dplyr)
library(qgraph)
library(ggplot2)
library(gridExtra)

# Means over time ====
data.wide = sapply(data.ref, as.numeric) 
data.wide = as.data.frame(data.wide)

colnames(data.wide) = c("ID", "group",
                        rep(paste0(node.names[-1], "_base")),
                        rep(paste0(node.names[-1], "_2wk")),
                        rep(paste0(node.names[-1], "_6wk")),
                        rep(paste0(node.names[-1], "_12wk")))

data.long = reshape(data.wide,
                    direction = "long",
                    varying = c(colnames(data.wide)[-c(1,2)]),
                    idvar = "ID",
                    sep = "_") %>% as.data.frame()

data.means = data.long %>%
  group_by(time, group) %>%
  summarise_at(vars(BAD:SLE), list(mean = mean, sd = sd), na.rm = TRUE)

data.means.se =  data.long %>%
  group_by(time, group) %>%
  summarise_at(vars(BAD:SLE), list(se = function(x) {sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))}))

data.means.full = merge(data.means, data.means.se, by = c("time", "group"))

data.means.long = data.means.full %>%
  gather(key = "stat", value ='value', BAD_mean:SLE_se) %>%
  separate(stat, c("node", "stat")) %>%
  spread(stat,value)


node.desc2 = c("TREAT" = "Treatment",
               "BAD" = "Feeling bad \nabout oneself",
               "ANX" = "Anxious",
               "AFR" = "Afraid",
               "FAI" = "Past \nfailure",
               "GUI" = "Guilt feelings",
               "PUN" = "Punishment \nfeelings",
               "CRY" = "Crying",
               "IND" = "Indecisiveness",
               "LIB" = "Libido \nproblems",
               "PHY" = "Physical \nhealth",
               "IMP" = "Self-reported \nimprovement",
               "TIR" = "Tired",
               "WOR" = "Worried",
               "ANH" = "Anhedonia",
               "DIS" = "Disliking \noneself",
               "RES" = "Restless \nor slow",
               "SUI" = "Suicidal \nthoughts",
               "SAD" = "Feeling \ndepressed",
               "APP" = "Appetite \nproblems",
               "CON" = "Concentration \nproblems",
               "SLE" = "Sleep \nproblems")
label.node = function(x){paste0(x, ": ", "\n", node.desc2[names(node.desc2) == x])}
data.means.long$label = sapply(data.means.long$node, FUN = label.node)

p = data.means.long %>% 
  ggplot(., aes(x = factor(time, levels = c("base", "2wk", "6wk", "12wk")), y = mean, group = as.factor(group), color = as.factor(group)))+ 
  geom_point(size = 0.5) + 
  geom_line(lwd = 0.3) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3) +
  facet_wrap(~label, ncol =7) +
  xlab("Time") +
  ylab("Mean")+
  scale_color_manual(values =  c("red", "black"),
                     name = "Group", labels = c("Placebo", "Sertraline")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 8),
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.position = "bottom")
p


# Plot IMP on a different scale
noimpplot = function(x, ymin, ymax, position) {plot1 = data.means.long %>% 
  filter(node == x) %>%
  ggplot(., aes(x = factor(time, levels = c("base", "2wk", "6wk", "12wk")), 
                y = mean, group = as.factor(group), color = as.factor(group)))+ 
  geom_point(size = 0.5) + 
  geom_line(lwd = 0.3) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3) +
  facet_wrap(~label) +
  scale_color_manual(values =  c("red", "black"),
                     name = "Group", labels = c("Placebo", "Sertraline")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 8),
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.position = position,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_cartesian(ylim = c(ymin,ymax))
return(plot1)
}

singleplots = c(lapply(node.names[2:11], FUN = noimpplot, ymin = 0, ymax = 2, position = "none"),
                lapply(node.names[12], FUN = noimpplot, ymin = 2, ymax = 4, position = "none"),
                lapply(node.names[13:21], FUN = noimpplot, ymin = 0, ymax = 2, position = "none"),
                lapply(node.names[22], FUN = noimpplot, ymin = 0, ymax = 2, position = "none"))

pdf(paste0(results, "plots/", "mean_change_sympt_imp.pdf"), width = 10, height = 7)
do.call("grid.arrange", args = list(grobs = singleplots, heights = c(1,1,1), ncol = 7, nrow = 3, bottom = "Time", left = "Mean"))
dev.off()

## Plot effect sizes and CIs ====

get_eta = function(x) {
  dataframe = data.frame(
    Symptom = names(mlm.eta[x]),
    Eta = mlm.eta[[x]]$Eta2_partial,
    CI_low = mlm.eta[[x]]$CI_low,
    CI_high = mlm.eta[[x]]$CI_high,
    Effect= mlm.eta[[x]]$Parameter
  )
  return(dataframe)
}

eta_list = lapply(node.names[2:22], get_eta) 

eta_df = bind_rows(eta_list)

eta_df2 = eta_df %>% filter(Effect == "time" |
                              Effect == "factor(TREAT_B)" |
                              Effect == "factor(TREAT_B):time") %>%
  mutate_at(vars(Effect), funs(recode(., `time`= "Time", `factor(TREAT_B)`= "Group", `factor(TREAT_B):time`= "Group x Time")))


ggplot(eta_df2,aes(y = Symptom, x = Eta, xmin = CI_low, xmax = CI_high, color = Effect)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbarh(height = 0.2, position = position_dodge(width = 1)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  facet_wrap(~factor(Effect, levels = c("Time", "Group", "Group x Time"))) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#009E73", "#E69F00", "#0072B2")) +
  labs(x = c(expression(eta^2)))


# Plot with significance values

get_pval = function(x) {
  dataframe = data.frame(
    Symptom = names(mlm.anova[x]),
    Pval = mlm.anova[[x]]$`Pr(>F)`,
    Effect = rownames(mlm.anova[[x]])
  )
  return(dataframe)
}

pval_list = lapply(node.names[2:22], get_pval) 
pval_df = bind_rows(pval_list)

pval_df2 = pval_df %>% filter(Effect == "time" |
                                Effect == "factor(TREAT_B)" |
                                Effect == "factor(TREAT_B):time") %>%
  mutate_at(vars(Effect), funs(recode(., `time`= "Time", `factor(TREAT_B)`= "Group", `factor(TREAT_B):time`= "Group x Time")))

# Correct p-values for multiple comparisons (21 tests)
pval_df2$Pval.Adj = rep(NA, 63)

pval_group = pval_df2 %>% 
  filter(Effect == "Group") %>% 
  mutate(Pval.Adj = p.adjust(Pval, n = 21, method = "fdr"))
pval_time = pval_df2 %>% 
  filter(Effect == "Time") %>% 
  mutate(Pval.Adj = p.adjust(Pval, n = 21, method = "fdr"))
pval_timegroup = pval_df2 %>% 
  filter(Effect == "Group x Time") %>% 
  mutate(Pval.Adj = p.adjust(Pval, n = 21, method = "fdr"))

pval_df3 = rbind(pval_group, pval_time, pval_timegroup)
  
eta_pval_df = merge(eta_df2,pval_df3, by = c("Symptom", "Effect"))


pdf(paste0(results, "plots/", "eta_sq_CIs_sig_adj.pdf"), width = 5, height = 5)
print(ggplot(eta_pval_df,aes(y = Symptom, x = Eta, xmin = CI_low, xmax = CI_high, color = Effect)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbarh(height = 0.2, position = position_dodge(width = 1)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~factor(Effect, levels = c("Time", "Group", "Group x Time"))) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#009E73", "#E69F00", "#0072B2")) +
  labs(x = c(expression(eta^2))) +
  geom_point(data = eta_pval_df[eta_pval_df$Pval.Adj < 0.05, ], 
             aes(y = Symptom, x = CI_high + 0.004), shape = "*", size=4.3, color="black")) 
dev.off()

# Get ranges of significant effect sizes
eta_pval_df %>% filter(Pval.Adj < 0.05) %>% filter(Effect == "Group") %>% select(Eta) %>% range() %>% round(., 3)

# Save full results of linear models
summary_table_anova_full = summary_table_anova %>%
  filter(Effect == "Group" | Effect == "Time" | Effect == "Group x Time") %>% 
  select(Symptom, Effect, Fvalue, Df) %>%
  merge(., eta_pval_df, by = c("Symptom", "Effect")) %>%
  mutate_at(vars(Eta:Pval.Adj), round, 3)

write.csv(summary_table_anova_full, paste0(results, "tables/", "summary_table_anova_full.csv"))
# Contemporaneous networks ====

# Treatment node in the centre
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  colnames(cormat) = node.names
  rownames(cormat) = node.names
  data.frame(
    node1 = rownames(cormat)[row(cormat)[ut]],
    node2 = rownames(cormat)[col(cormat)[ut]],
    est  =(cormat)[ut]
  )
}

# Sourcing qgraph support functions and a modified version of qgraph that plots highlighted relevant edges, created for a previous project
ti.layout = ti.net.lag1.plot$layout
avg.layout = averageLayout(twoweeks.mgm.res, sixweeks.mgm.res, twelveweeks.mgm.res,
                           qgraph(ti.net.lag1, layout = "spring"),
                           qgraph(ti.net.lag2, layout = "spring"))
source(paste0(results, "/plots/qgraph.alt.R"))
source(paste0(results, "/plots/support.functions.R"))

unfade.single = function(network, title, curve.val) {
  # flatten correlation matrix to list of edges
  network.sort = flattenCorrMatrix(network) 
  # reorder edges so those connected to TREAT are first (to be plotted last)
  network.sort$edgeSort = c(ifelse(network.sort$node1 == "TREAT",1,2))
  network.sort = network.sort %>% arrange(factor(edgeSort)) %>% select(node1:est)
  
  # find indices of edges connected to TREAT
  group.index = which(network.sort$node1 == "TREAT")
  
  # change color of edges connected to TREAT
  color = c(ifelse(network.sort$est > 0, "lightblue1", "lightgoldenrod1"))
  color[group.index] = c(ifelse(network.sort$est[group.index] > 0, "darkblue", "red"))
  
  n_edges = length(network.sort$est[abs(network.sort$est) > 0]) # find n of edges
  
  # curve edges connected to TREAT
  curve = rep(0, length(network.sort$est))
  curve[group.index] = c(ifelse(network.sort$est[group.index] != 0, 2, 0))
  curve[curve == 2] = curve.val
  
  # custom layout
  layout.centre = avg.layout
  layout.centre[1,] = c(0,0) 
  layout.centre[15,] = c(-0.3,0.3)
  plot2 = qgraph.alt(network.sort, 
                     theme = "colorblind",
                     directed = FALSE,
                     edgeSort = n_edges:1,
                     title = title, 
                     edge.color = color,
                     layout = layout.centre,
                     curveAll = TRUE,
                     curve = curve,
                     shape = c("square", rep("circle", 21)),
                     color = c("grey87", rep("white", 21)),
                     cut = 0.1
  )
}


pdf(file = "results/plots/2wk_centre.pdf", width = 5, height = 5)
unfade.single(twoweeks.mgm.res$graph, title = "", curve.val = c(0.6, 0.3,0.5, -0.5, 0, 0, 0, 0, 0, 0))
dev.off()

pdf(file = "results/plots/6wk_centre.pdf", width = 5, height = 5)
unfade.single(sixweeks.mgm.res$graph, title = "", curve.val = c(0.6, -0.7, 0.5, 0, 0, 0, 0, -0.5, 0, 0))
dev.off()

pdf(file = "results/plots/12wk_centre.pdf", width = 5, height = 5)
unfade.single(twelveweeks.mgm.res$graph, title = "", curve.val = c(0, 1, -1, 0, 0, 0, -0.5))
dev.off()

# Temporally lagged networks ====


# Create specific function to unfade + hightlight treatment edges
unfade = function(network, title, curve.val) {
  
  # reorder edges so those connected to TREAT are first (to be plotted last)
  network.sort = network 
  network.sort$edgeSort = c(ifelse(network.sort$t1.node == "TREAT",1,2))
  network.sort = network.sort %>% arrange(factor(edgeSort)) %>% select(t1.node:est.net)
  
  # find indices of edges connected to TREAT
  group.index = which(network.sort$t1.node == "TREAT")
  
  # change color of edges connected to TREAT (optional: line type)
  lty = c(ifelse(network.sort$est.net > 0, "solid", "solid"))
  color = c(ifelse(network.sort$est.net > 0, "lightblue1", "lightgoldenrod1"))
  
  lty[group.index] = "solid"
  color[group.index] = c(ifelse(network.sort$est.net[group.index] > 0, "darkblue", "red"))
  
  # curve edges connected to TREAT
  curve = rep(0.5, length(network.sort$est))
  curve[group.index] = c(ifelse(network.sort$est[group.index] != 0, 2, 0))
  curve[curve == 2] = curve.val
  
  n_edges = length(network.sort$t1.node[abs(network.sort$est.net) > 0])
  layout.centre = avg.layout
  layout.centre[1,] = c(0,0) 
  layout.centre[15,] = c(-0.3,0.3) 
  plot2 = qgraph.alt(network.sort, 
                     theme = "colorblind", 
                     title = title, 
                     edge.color = color,
                     lty = lty,
                     edgeSort = n_edges:1,
                     layout = layout.centre,
                     shape = c("square", rep("circle", 21)),
                     color = c("gold", rep("white", 21)),
                     cut = 0.1,
                     curveAll = TRUE,
                     curve = curve)
}


write.csv(ti.net.lag1, file = "results/tables/ti_lag1.csv", quote = FALSE, row.names = FALSE)
write.csv(ti.net.lag2, file = "results/tables/ti_lag2.csv", quote = FALSE, row.names = FALSE)


pdf("results/plots/unfaded_ti_lag1_2lags.pdf", width = 5, height = 5)
unfade(ti.net.lag1, title = "", curve.val = c(0,0,-0.8,0.4,1,0,0,0,0,-0.6))
dev.off()

pdf("results/plots/unfaded_ti_lag2_2lags.pdf", width = 5, height = 5)
unfade(ti.net.lag2, title = "", curve.val = c(0,0.5,0,0,0,0))
dev.off()

