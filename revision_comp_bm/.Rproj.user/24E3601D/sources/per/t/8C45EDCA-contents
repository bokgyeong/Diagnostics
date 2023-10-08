rm(list=ls())
library(dplyr); library(ggplot2); library(egg); library(gridExtra); library(grid)
library(batchmeans); library(mcmcse)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

al = 0.01
# al = 0.05

n = 300000
# n = 400000

load(paste0('ACDAIKS/Trunc_knu/n', n, '/bids2AIKSTrunc', 8, '.RData'))
aiksCutoff = quantile(aiksBootTrunc, 1-al)


# ==============================================================================
# Diagnostics
# ==============================================================================
# Truncation algorithm
# ------------------------------------------------------------------------------
diagTrunc = data.frame()
# m = c(3:4, 8:10, 20, 30)
m = c(3, 8:10, 30)
# m = c(8:10, 20, 30)


for(i in m){
  ### ACD
  load(paste0('ACDAIKS/Trunc_knu/n', n, '/bids2ACDTrunc', i, '.RData'))
  df = data.frame(Method = 'ACD', Algorithm = 'Truncation', m = as.factor(i), 
                  Value = acdTrunc, rtime = timeacdTrunc, 
                  Below = ifelse(acdTrunc < qchisq(1-al, 55), TRUE, FALSE))
  diagTrunc = rbind(diagTrunc, df)
  
  ### AIKS
  load(paste0('ACDAIKS/Trunc_knu/n', n, '/bids2AIKSTrunc', i, '.RData'))
  df = data.frame(Method = 'AIKS', Algorithm = 'Truncation', m = as.factor(i),
                  Value = aiksTrunc, rtime = timeaiksTrunc, 
                  Below = ifelse(aiksTrunc < aiksCutoff, TRUE, FALSE))
  diagTrunc = rbind(diagTrunc, df)
  
  ### ESS
  load(paste0('ACDAIKS/Trunc_knu/bids2Trunc', i, '.RData'))
  df = data.frame(Method = 'ESS', Algorithm = 'Truncation', m = as.factor(i),
                  Value = round(min(apply(Trunc[burn+1:n,], 2, ess))),
                  # Value = round(multiESS(Trunc[burn+1:n,])), 
                  rtime = NA, Below = NA)
  diagTrunc = rbind(diagTrunc, df)
}
diagTrunc




# ------------------------------------------------------------------------------
# ESS
# ------------------------------------------------------------------------------

plotESSTrunc = diagTrunc %>% 
  filter(Method == 'ESS') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(size = 2, shape = 16) +
  geom_line(color = 'gray45') +
  # coord_cartesian(ylim = c(NA, 11000)) +
  coord_cartesian(ylim = c(NA, 8000)) +
  labs(x = 'k', y = 'Minimum ESS')



# ------------------------------------------------------------------------------
# ACD
# ------------------------------------------------------------------------------
plotACDTrunc = diagTrunc %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = qchisq(1-al, 55), linetype = 2) +
  # coord_cartesian(ylim = c(NA, 200)) +
  scale_y_log10() +
  labs(x = 'k', y = 'ACD') +
  theme(legend.position = 'none')


# ------------------------------------------------------------------------------
# AIKS
# ------------------------------------------------------------------------------
plotAIKSTrunc = diagTrunc %>% 
  filter(Method == 'AIKS') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = aiksCutoff, linetype = 2) +
  # coord_cartesian(ylim = c(NA, 10)) +
  scale_y_log10() +
  labs(x = 'k', y = 'AIKS') +
  theme(legend.position = 'none')



# ------------------------------------------------------------------------------
# Combine plots
# ------------------------------------------------------------------------------
# plotDiag = ggarrange(plotESSTrunc, plotACDTrunc, plotAIKSTrunc, nrow = 1, top = grid::textGrob('Truncation method: Diagnostics', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
plotDiag = ggarrange(plotESSTrunc, plotACDTrunc, plotAIKSTrunc, nrow = 1)

ggsave(plot = plotDiag, 
       filename = paste0('figures_paper/n', n, '/compDiag', 100*(1-al),'.eps'), 
       # width = 6.2, height = 2)
       width = 6.8, height = 2)






# ==============================================================================
# Posterior samples
# ==============================================================================
sampleAll = data.frame()

### Truncation
for(i in c(3, 9)){
  load(paste0('ACDAIKS/Trunc_knu/bids2Trunc', i, '.RData'))
  df = data.frame(Sample = 'NormTrunc', m = as.factor(i), parameter = Trunc[burn+1:n,], time = rtime / nrow(Trunc) * n)
  sampleAll = rbind(sampleAll, df)
}

### Gold standard
load('ACDAIKS/Murray_knu/bids2Murray.RData')
df = data.frame(Sample = 'Gold standard', m = NA, parameter = Murray, time = NA)
sampleAll = rbind(sampleAll, df)




# ------------------------------------------------------------------------------
# Tail probabilities
# ------------------------------------------------------------------------------

cut.1 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.1, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.2 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.2, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.3 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.3, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.4 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.4, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.5 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.5, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.6 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.6, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.7 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.7, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.8 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.8, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.9 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.9, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()
cut.10 = sampleAll %>% filter(Sample == 'Gold standard') %>% summarise(Cut = quantile(parameter.10, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()

sumstatsAll.1 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.1), 2),
            Median = round(median(parameter.1), 2),
            SD = round(sd(parameter.1), 2),
            Prob1 = round(mean(parameter.1 <= cut.1[1]), 2),
            Prob2 = round(mean(parameter.1 >= cut.1[2] & parameter.1 <= cut.1[3]), 2),
            Prob3 = round(mean(parameter.1 >= cut.1[3] & parameter.1 <= cut.1[4]), 2),
            Prob4 = round(mean(parameter.1 >= cut.1[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.2 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.2), 2),
            Median = round(median(parameter.2), 2),
            SD = round(sd(parameter.2), 2),
            Prob1 = round(mean(parameter.2 <= cut.2[1]), 2),
            Prob2 = round(mean(parameter.2 >= cut.2[2] & parameter.2 <= cut.2[3]), 2),
            Prob3 = round(mean(parameter.2 >= cut.2[3] & parameter.2 <= cut.2[4]), 2),
            Prob4 = round(mean(parameter.2 >= cut.2[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.3 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.3), 2),
            Median = round(median(parameter.3), 2),
            SD = round(sd(parameter.3), 2),
            Prob1 = round(mean(parameter.3 <= cut.3[1]), 2),
            Prob2 = round(mean(parameter.3 >= cut.3[2] & parameter.3 <= cut.3[3]), 2),
            Prob3 = round(mean(parameter.3 >= cut.3[3] & parameter.3 <= cut.3[4]), 2),
            Prob4 = round(mean(parameter.3 >= cut.3[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.4 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.4), 2),
            Median = round(median(parameter.4), 2),
            SD = round(sd(parameter.4), 2),
            Prob1 = round(mean(parameter.4 <= cut.4[1]), 2),
            Prob2 = round(mean(parameter.4 >= cut.4[2] & parameter.4 <= cut.4[3]), 2),
            Prob3 = round(mean(parameter.4 >= cut.4[3] & parameter.4 <= cut.4[4]), 2),
            Prob4 = round(mean(parameter.4 >= cut.4[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.5 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.5), 2),
            Median = round(median(parameter.5), 2),
            SD = round(sd(parameter.5), 2),
            Prob1 = round(mean(parameter.5 <= cut.5[1]), 2),
            Prob2 = round(mean(parameter.5 >= cut.5[2] & parameter.5 <= cut.5[3]), 2),
            Prob3 = round(mean(parameter.5 >= cut.5[3] & parameter.5 <= cut.5[4]), 2),
            Prob4 = round(mean(parameter.5 >= cut.5[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.6 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.6), 2),
            Median = round(median(parameter.6), 2),
            SD = round(sd(parameter.6), 2),
            Prob1 = round(mean(parameter.6 <= cut.6[1]), 2),
            Prob2 = round(mean(parameter.6 >= cut.6[2] & parameter.6 <= cut.6[3]), 2),
            Prob3 = round(mean(parameter.6 >= cut.6[3] & parameter.6 <= cut.6[4]), 2),
            Prob4 = round(mean(parameter.6 >= cut.6[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.7 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.7), 2),
            Median = round(median(parameter.7), 2),
            SD = round(sd(parameter.7), 2),
            Prob1 = round(mean(parameter.7 <= cut.7[1]), 2),
            Prob2 = round(mean(parameter.7 >= cut.7[2] & parameter.7 <= cut.7[3]), 2),
            Prob3 = round(mean(parameter.7 >= cut.7[3] & parameter.7 <= cut.7[4]), 2),
            Prob4 = round(mean(parameter.7 >= cut.7[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.8 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.8), 2),
            # Median = round(median(parameter.8), 2),
            Median = format(round(median(parameter.8), 2), nsmall = 2),
            SD = round(sd(parameter.8), 2),
            Prob1 = round(mean(parameter.8 <= cut.8[1]), 2),
            Prob2 = round(mean(parameter.8 >= cut.8[2] & parameter.8 <= cut.8[3]), 2),
            Prob3 = round(mean(parameter.8 >= cut.8[3] & parameter.8 <= cut.8[4]), 2),
            Prob4 = round(mean(parameter.8 >= cut.8[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.9 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.9), 2),
            Median = round(median(parameter.9), 2),
            SD = round(sd(parameter.9), 2),
            Prob1 = round(mean(parameter.9 <= cut.9[1]), 2),
            Prob2 = round(mean(parameter.9 >= cut.9[2] & parameter.9 <= cut.9[3]), 2),
            Prob3 = round(mean(parameter.9 >= cut.9[3] & parameter.9 <= cut.9[4]), 2),
            Prob4 = round(mean(parameter.9 >= cut.9[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))

sumstatsAll.10 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.10), 2),
            Median = round(median(parameter.10), 2),
            SD = round(sd(parameter.10), 2),
            Prob1 = round(mean(parameter.10 <= cut.10[1]), 2),
            Prob2 = round(mean(parameter.10 >= cut.10[2] & parameter.10 <= cut.10[3]), 2),
            Prob3 = round(mean(parameter.10 >= cut.10[3] & parameter.10 <= cut.10[4]), 2),
            Prob4 = round(mean(parameter.10 >= cut.10[5]), 2),
            Time = format(round(mean(time)/60, 2), nsmall = 2))


sumstatsAll.1 %>% dplyr::select(Sample, m, Median, Prob1, Prob4, Time)
# sumstatsAll.2 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.3 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.4 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.5 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.6 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.7 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.8 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time) # difference!
# sumstatsAll.9 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)
# sumstatsAll.10 %>% dplyr::select(Sample, m, Median, SD, Prob1, Prob2, Prob3, Prob4, Time)




# ==============================================================================
# Variability in Diagnostics
# ==============================================================================
# ACD
# ------------------------------------------------------------------------------
varDiagTrunc = data.frame()
m = c(8:10, 30, 100)


for(i in m){
  ### ACD
  load(paste0('rep100/Trunc', i, '/bids2ACDTrunc', i, '.RData'))
  ACDs = ACD_rep100[ACD_rep100 != 0]
  df = data.frame(Method = 'ACD', Algorithm = 'Truncation', m = as.factor(i), 
                  Mean = round(mean(ACDs), 2), SD = round(sd(ACDs), 2), 
                  Min = round(min(ACDs), 2), Max = round(max(ACDs), 2),
                  Below = ifelse(max(ACDs) < qchisq(1-al, 55), TRUE, FALSE))
  varDiagTrunc = rbind(varDiagTrunc, df)

  ### AIKS
  load(paste0('rep100/Trunc', i, '/bids2AIKSTrunc', i, '.RData'))
  AIKSs = AIKS_rep100[AIKS_rep100 != 0]
  df = data.frame(Method = 'AIKS', Algorithm = 'Truncation', m = as.factor(i), 
                  Mean = round(mean(AIKSs), 3), SD = round(sd(AIKSs), 5), 
                  Min = round(min(AIKSs), 4), Max = round(max(AIKSs), 4),
                  Below = ifelse(max(AIKSs) < aiksCutoff, TRUE, FALSE))
  varDiagTrunc = rbind(varDiagTrunc, df)
}

varDiagTrunc %>% filter(Method == 'ACD')
varDiagTrunc %>% filter(Method == 'AIKS')
