rm(list=ls())
library(dplyr); library(ggplot2); library(egg); library(gridExtra); library(grid)
library(batchmeans)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

al = 0.01

n = 250000

load(paste0('ACDAIKS/Liang/n', n, '/simAIKSLiang', 1, '.RData'))
aiksCutoff = quantile(aiksBootLiang, 1-al)


# ==============================================================================
# DMH 
# ==============================================================================
diagLiang = data.frame() 
m = 1:5
# m = 1:4

for(i in m){
  ### ACD
  load(paste0('ACDAIKS/Liang/n', n, '/simACDLiang', i, '.RData'))
  df = data.frame(Method = 'ACD', Algorithm = 'DMH', m = as.factor(i), Value = acdLiang, 
                  Below = ifelse(acdLiang < qchisq(1-al, 3), TRUE, FALSE))
  diagLiang = rbind(diagLiang, df)
  
  ### AIKS
  load(paste0('ACDAIKS/Liang/n', n, '/simAIKSLiang', i, '.RData'))
  df = data.frame(Method = 'AIKS', Algorithm = 'DMH', m = as.factor(i), Value = aiksLiang,
                  Below = ifelse(aiksLiang < aiksCutoff, TRUE, FALSE))
  diagLiang = rbind(diagLiang, df)
  
  ### ESS
  load(paste0('ACDAIKS/Liang/simLiang', i, '.RData'))
  df = data.frame(Method = 'ESS', Algorithm = 'DMH', m = as.factor(i),
                  Value = round(min(apply(Liang[burnin+1:n,], 2, ess))),
                  Below = NA)
  diagLiang = rbind(diagLiang, df)
}
diagLiang




# ==============================================================================
# plot
# ==============================================================================

# ------------------------------------------------------------------------------
# ESS
# ------------------------------------------------------------------------------
plotESSLiang = diagLiang %>% 
  filter(Method == 'ESS') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(size = 2, shape = 16) +
  geom_line(color = 'gray45') +
  coord_cartesian(ylim = c(NA, 17000)) +
  labs(x = 'm', y = 'Minimum ESS')



# ------------------------------------------------------------------------------
# ACD
# ------------------------------------------------------------------------------
plotACDLiang = diagLiang %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = qchisq(1-al, 3), linetype = 2) +
  labs(x = 'm', y = 'ACD') +
  theme(legend.position = 'none')




# ------------------------------------------------------------------------------
# AIKS
# ------------------------------------------------------------------------------
plotAIKSLiang = diagLiang %>% 
  filter(Method == 'AIKS') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = aiksCutoff, linetype = 2) +
  labs(x = 'm', y = 'AIKS') +
  # labs(x = 'm', y = 'AIKS', title = 'DMH') +
  theme(legend.position = 'none')





# ------------------------------------------------------------------------------
# Combine plots
# ------------------------------------------------------------------------------

# plotDiag = ggarrange(plotESSLiang, plotACDLiang, plotAIKSLiang, nrow = 1, top = grid::textGrob('DMH: Diagnostics', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
plotDiag = ggarrange(plotESSLiang, plotACDLiang, plotAIKSLiang, nrow = 1)

ggsave(plot = plotDiag, 
       filename = paste0('figures_paper/ergmDiag', 100*(1-al), '.eps'), 
       width = 6.2, height = 2)


# # ==============================================================================
# # Posterior samples
# # ==============================================================================
# Names = c('ESS-selected', 'ACD-selected', 'AIKS-selected')
# DMH_selected = c(4, 3, 2)
# 
# sampleAll = data.frame()
# 
# ### DMH
# for(i in 1:length(Names)){
#   load(paste0('ACDAIKS/Liang/simLiang', DMH_selected[i], '.RData'))
#   if(i == 'AIKS-selected'){
#     df = data.frame(Sample = Names[i], parameter = Liang[burnin+seq(25, n, by = 25),], time = timeLiang / nrow(Liang) *n)
#   } else {
#     df = data.frame(Sample = Names[i], parameter = Liang[burnin+1:n,], time = timeLiang / nrow(Liang) *n)
#   }
#   sampleAll = rbind(sampleAll, df)
# }
# 
# ### Gold standard
# load(paste0('ACDAIKS/Liang/simLiang', 20, '.RData'))
# df = data.frame(Sample = 'Gold standard', parameter = Liang[burnin+1:n,], time = timeLiang / nrow(Liang) *n)
# sampleAll = rbind(sampleAll, df)
# 
# 
# 
# # ==============================================================================
# # Tail probabilities
# # ==============================================================================
# 
# cut.1 = sampleAll %>% filter(Sample == 'Gold standard') %>% 
#   summarise(Cut = quantile(parameter.1, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
#   # summarise(Cut = quantile(parameter.1, c(0.05, 0.49, 0.5, 0.51, 0.95))) %>% unlist()
# 
# cut.2 = sampleAll %>% filter(Sample == 'Gold standard') %>% 
#   summarise(Cut = quantile(parameter.2, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
#   # summarise(Cut = quantile(parameter.2, c(0.05, 0.49, 0.5, 0.51, 0.95))) %>% unlist()
# 
# sumstatsAll.1 = sampleAll %>%
#   group_by(Sample) %>%
#   summarise(Mean = round(mean(parameter.1), 2),
#             Median = round(median(parameter.1), 2),
#             SD = round(sd(parameter.1), 2),
#             Prob1 = round(mean(parameter.1 <= cut.1[1]), 2),
#             Prob2 = round(mean(parameter.1 >= cut.1[2] & parameter.1 <= cut.1[3]), 2),
#             Prob3 = round(mean(parameter.1 >= cut.1[3] & parameter.1 <= cut.1[4]), 2),
#             Prob4 = round(mean(parameter.1 >= cut.1[5]), 2),
#             Time = format(round(mean(time)/60/60, 2), nsmall = 2))
# 
# sumstatsAll.2 = sampleAll %>%
#   # group_by(Sample, m) %>%
#   group_by(Sample) %>%
#   summarise(Mean = round(mean(parameter.2), 2),
#             Median = round(median(parameter.2), 2),
#             SD = round(sd(parameter.2), 2),
#             Prob1 = round(mean(parameter.2 <= cut.2[1]), 2),
#             Prob2 = round(mean(parameter.2 >= cut.2[2] & parameter.2 <= cut.2[3]), 2),
#             Prob3 = round(mean(parameter.2 >= cut.2[3] & parameter.2 <= cut.2[4]), 2),
#             Prob4 = round(mean(parameter.2 >= cut.2[5]), 2),
#             Time = format(round(mean(time)/60/60, 2), nsmall = 2))
# 
# 
# sumstatsAll.1 %>% select(Sample, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.2 %>% select(Sample, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# 
# 
