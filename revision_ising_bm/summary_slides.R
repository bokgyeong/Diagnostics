rm(list=ls())
library(dplyr); library(ggplot2); library(egg); library(gridExtra); library(grid); library(scales)
library(batchmeans); library(coda)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  inv <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  # return the transformation
  return(trans_new("squished", trans, inv))
}


al = 0.01

n = 120000
# n = 100000
# n = 50000

load(paste0('ACDAIKS/Liang/n', n, '/simAIKSLiang', 1, '.RData'))
aiksCutoff = quantile(aiksBootLiang, 1-al)


# ==============================================================================
# DMH 
# ==============================================================================
diagLiang = data.frame()
m = 1:5

for(i in m){
  ### ACD
  load(paste0('ACDAIKS/Liang/n', n, '/simACDLiang', i, '.RData'))
  df = data.frame(Method = 'ACD', Algorithm = 'DMH', m = as.factor(i), 
                  Value = acdLiang, rtime = timeacdLiang, 
                  Below = ifelse(acdLiang < qchisq(1-al, 1), TRUE, FALSE))
  diagLiang = rbind(diagLiang, df)
  
  ### AIKS
  load(paste0('ACDAIKS/Liang/n', n, '/simAIKSLiang', i, '.RData'))
  df = data.frame(Method = 'AIKS', Algorithm = 'DMH', m = as.factor(i),
                  Value = aiksLiang, rtime = timeaiksLiang, 
                  Below = ifelse(aiksLiang < aiksCutoff, TRUE, FALSE))
  diagLiang = rbind(diagLiang, df)
  
  ### ESS
  load(paste0('ACDAIKS/Liang/simLiang', i, '.RData'))
  df = data.frame(Method = 'ESS', Algorithm = 'DMH', m = as.factor(i),
                  Value = round(ess(as.vector(Liang[burnin + 1:n,]))),
                  rtime = NA, Below = NA)
  diagLiang = rbind(diagLiang, df)
}




# ==============================================================================
# Atchade
# ==============================================================================
diagAtchade = data.frame()
# d = c(10, 20, 50, 100, 200, 400)
d = c(10, 20, 100, 200, 400)

for(i in d){
  ### ACD
  load(paste0('ACDAIKS/Atchade/n', n, '/simACDAtchade', i, '.RData'))
  df = data.frame(Method = 'ACD', Algorithm = 'ALR', d = as.factor(i),
                  Value = acdAtchade, rtime = timeacdAtchade,
                  Below = ifelse(acdAtchade < qchisq(1-al, 1), TRUE, FALSE))
  diagAtchade = rbind(diagAtchade, df)
  
  ### AIKS
  load(paste0('ACDAIKS/Atchade/n', n, '/simAIKSAtchade', i, '.RData'))
  df = data.frame(Method = 'AIKS', Algorithm = 'ALR', d = as.factor(i),
                  Value = aiksAtchade, rtime = timeaiksAtchade, 
                  Below = ifelse(aiksAtchade < aiksCutoff, TRUE, FALSE))
  diagAtchade = rbind(diagAtchade, df)
  
  ### ESS
  load(paste0('ACDAIKS/Atchade/simAtchade', i, '.RData'))
  df = data.frame(Method = 'ESS', Algorithm = 'ALR', d = as.factor(i),
                  Value = round(ess(as.vector(Atchade[burnin + 1:n,]))),
                  rtime = NA, Below = NA)
  diagAtchade = rbind(diagAtchade, df)
}


# ==============================================================================
# AEX
# ==============================================================================
diagAEX = data.frame()
# d = c(100, 400, 450, 490, 500, 510)
d = c(100, 400, 490, 500, 510)


for(i in d){
  ### ACD
  load(paste0('ACDAIKS/AEX/n', n, '/simACDAEX', i, '.RData'))
  df = data.frame(Method = 'ACD', Algorithm = 'AEX', d = as.factor(i),
                  Value = acdAEX, rtime = timeacdAEX,
                  Below = ifelse(acdAEX < qchisq(1-al, 1), TRUE, FALSE))
  diagAEX = rbind(diagAEX, df)
  
  ### AIKS
  load(paste0('ACDAIKS/AEX/n', n, '/simAIKSAEX', i, '.RData'))
  df = data.frame(Method = 'AIKS', Algorithm = 'AEX', d = as.factor(i),
                  Value = aiksAEX, rtime = timeaiksAEX, 
                  Below = ifelse(aiksAEX < aiksCutoff, TRUE, FALSE))
  diagAEX = rbind(diagAEX, df)
  
  ### ESS
  load(paste0('ACDAIKS/AEX/simAEX', i, '.RData'))
  df = data.frame(Method = 'ESS', Algorithm = 'AEX', d = as.factor(i),
                  Value = round(ess(as.vector(AEX[burnin + 1:n,]))),
                  rtime = NA, Below = NA)
  diagAEX = rbind(diagAEX, df)
}
diagAEX




# ==============================================================================
# plot
# ==============================================================================
rangeDiag = rbind(diagLiang %>% group_by(Method) %>% summarise(range = range(Value)),
                  # diagAEX %>% group_by(Method) %>% summarise(range = range(Value)),
                  diagAtchade %>% group_by(Method) %>% summarise(range = range(Value)))
rangeDiag = rangeDiag %>% group_by(Method) %>% summarise(low = min(range), up = max(range))
rangeDiag[2, 3] = 0.5

# ------------------------------------------------------------------------------
# ESS 
# ------------------------------------------------------------------------------

plotESSLiang = diagLiang %>% 
  filter(Method == 'ESS') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(size = 2, shape = 16) +
  geom_line(color = 'gray45') +
  scale_y_log10() +
  # scale_y_sqrt() +
  # scale_y_continuous(trans = 'log2') +
  # coord_cartesian(ylim = rangeDiag %>% filter(Method == 'ESS') %>% select(low, up) %>% unlist()) +
  coord_cartesian(ylim = c(NA, 10000)) +
  labs(x = 'm', y = 'ESS')

plotESSAtchade = diagAtchade %>%
  filter(Method == 'ESS') %>%
  ggplot(aes(d, Value, group = 1)) +
  geom_point(size = 2, shape = 16) +
  geom_line(color = 'gray45') +
  scale_y_log10() +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'ESS') %>% select(low, up) %>% unlist()) +
  labs(x = 'd', y = 'ESS')

plotESSAEX = diagAEX %>%
  filter(Method == 'ESS') %>%
  ggplot(aes(d, Value, group = 1)) +
  geom_point(size = 2, shape = 16) +
  geom_line(color = 'gray45') +
  scale_y_log10() +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'ESS') %>% select(low, up) %>% unlist()) +
  labs(x = 'd', y = 'ESS')

# ggarrange(plotESSAEX, plotESSAtchade, plotESSLiang, nrow = 1)


# ------------------------------------------------------------------------------
# ACD
# ------------------------------------------------------------------------------

plotACDLiang = diagLiang %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = qchisq(1-al, 1), linetype = 2) +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'ACD') %>% select(low, up) %>% unlist()) +
  labs(x = 'm', y = 'ACD') +
  theme(legend.position = 'none')

plotACDAtchade = diagAtchade %>%
  filter(Method == 'ACD') %>%
  ggplot(aes(d, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = qchisq(1-al, 1), linetype = 2) +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'ACD') %>% select(low, up) %>% unlist()) +
  labs(x = 'd', y = 'ACD') +
  theme(legend.position = 'none')

plotACDAEX = diagAEX %>%
  filter(Method == 'ACD') %>%
  ggplot(aes(d, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = qchisq(1-al, 1), linetype = 2) +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'ACD') %>% select(low, up) %>% unlist()) +
  labs(x = 'd', y = 'ACD') +
  theme(legend.position = 'none')


# ggarrange(plotACDAEX, plotACDAtchade, plotACDLiang, nrow = 1)


# ------------------------------------------------------------------------------
# AIKS
# ------------------------------------------------------------------------------

plotAIKSLiang = diagLiang %>% 
  filter(Method == 'AIKS') %>% 
  ggplot(aes(m, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = aiksCutoff, linetype = 2) +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'AIKS') %>% select(low, up) %>% unlist()) +
  # scale_y_continuous(trans = squish_trans(0.1, 6, 200), breaks = c(0, 0.025, 0.05, 0.075, 0.2, 6)) +
  labs(x = 'm', y = 'AIKS') +
  # labs(x = 'm', y = 'AIKS', title = 'DMH') +
  theme(legend.position = 'none')

plotAIKSAtchade = diagAtchade %>%
  filter(Method == 'AIKS') %>%
  ggplot(aes(d, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = aiksCutoff, linetype = 2) +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'AIKS') %>% select(low, up) %>% unlist()) +
  # scale_y_continuous(trans = squish_trans(0.1, 6, 200), breaks = c(0, 0.025, 0.05, 0.075, 0.2, 6)) +
  labs(x = 'd', y = 'AIKS') +
  # labs(x = 'd', y = 'AIKS', title = 'ALR') +
  theme(legend.position = 'none')

plotAIKSAEX = diagAEX %>%
  filter(Method == 'AIKS') %>%
  ggplot(aes(d, Value, group = 1)) +
  geom_point(aes(shape = Below, color = Below), size = 2) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_line(color = 'gray45') +
  geom_hline(yintercept = aiksCutoff, linetype = 2) +
  coord_cartesian(ylim = rangeDiag %>% filter(Method == 'AIKS') %>% select(low, up) %>% unlist()) +
  # scale_y_continuous(trans = squish_trans(0.1, 6, 200), breaks = c(0, 0.025, 0.05, 0.075, 0.2, 6)) +
  labs(x = 'd', y = 'AIKS') +
  # labs(x = 'd', y = 'AIKS', title = 'AEX') +
  theme(legend.position = 'none')

# ggarrange(plotAIKSAEX, plotAIKSAtchade, plotAIKSLiang, nrow = 1)







# ==============================================================================
# Posterior samples
# ==============================================================================

### DMH
sampleLiang = data.frame()
for(i in unique(diagLiang$m)){
  load(paste0('ACDAIKS/Liang/simLiang', i, '.RData'))
  df = data.frame(m = as.factor(i), parameter = as.vector(Liang[burnin + 1:n]))
  sampleLiang = rbind(sampleLiang, df)
}

### ALR
sampleAtchade = data.frame()
for(i in unique(diagAtchade$d)){
  load(paste0('ACDAIKS/Atchade/simAtchade', i, '.RData'))
  df = data.frame(d = as.factor(i), parameter = as.vector(Atchade[burnin + 1:n]))
  sampleAtchade = rbind(sampleAtchade, df)
}

### AEX
sampleAEX = data.frame()
for(i in unique(diagAEX$d)){
  load(paste0('ACDAIKS/AEX/simAEX', i, '.RData'))
  df = data.frame(d = as.factor(i), parameter = as.vector(AEX[burnin + 1:n]))
  sampleAEX = rbind(sampleAEX, df)
}



# ------------------------------------------------------------------------------
# Kernel density
# ------------------------------------------------------------------------------
rangeSamp = range(c(sampleLiang$parameter,
                    # sampleAEX$parameter,
                    sampleAtchade$parameter
                    ))
rangeDens = c(0, 18)

plotSampLiang = ggplot(sampleLiang, aes(x = parameter)) +
  geom_density(aes(color = m, linetype = m)) +
  coord_cartesian(xlim = rangeSamp, ylim = rangeDens) +
  labs(x = expression(theta), y = 'Density') +
  theme(legend.position = 'right', legend.key.size = unit(0.5, 'cm')) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))

plotSampAtchade = ggplot(sampleAtchade, aes(x = parameter)) +
  geom_density(aes(color = d, linetype = d)) +
  coord_cartesian(xlim = rangeSamp, ylim = rangeDens) +
  labs(x = expression(theta), y = 'Density') +
  theme(legend.position = 'right', legend.key.size = unit(0.5, 'cm')) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))

plotSampAEX = ggplot(sampleAEX, aes(x = parameter)) +
  geom_density(aes(color = d, linetype = d)) +
  coord_cartesian(xlim = rangeSamp, ylim = rangeDens) +
  labs(x = expression(theta), y = 'Density') +
  theme(legend.position = 'right', legend.key.size = unit(0.5, 'cm')) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))



# ------------------------------------------------------------------------------
# Combine plots for ACD and AIKS
# ------------------------------------------------------------------------------

plotDiagAtchade = ggarrange(plotESSAtchade, plotACDAtchade, plotAIKSAtchade, plotSampAtchade, nrow = 1, top = grid::textGrob('(a) ALR: Diagnostics and estimated posterior densities', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
# plotDiagAEX = ggarrange(plotESSAEX, plotACDAEX, plotAIKSAEX, plotSampAEX, nrow = 1, top = grid::textGrob('(b) AEX: Diagnostics and estimated posterior densities', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
# plotDiagLiang = ggarrange(plotESSLiang, plotACDLiang, plotAIKSLiang, plotSampLiang, nrow = 1, top = grid::textGrob('(b) DMH: Diagnostics and estimated posterior densities', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
plotDiagLiang = ggarrange(plotESSLiang, plotACDLiang, plotAIKSLiang, nrow = 1)

# plotDiag = grid.arrange(plotDiagAtchade, plotDiagAEX, plotDiagLiang, ncol = 1)
plotDiag = grid.arrange(plotDiagAtchade, plotDiagLiang, ncol = 1)

ggsave(plot = plotDiag, 
       filename = paste0('figures_slides/isingDiag', 100*(1-al),'.eps'), 
       width = 8.4, height = 4)

# ggsave(width = 8.4, height = 2,
#        # plot = plotDiagAtchade,
#        # filename = paste0('figures_slides/isingDiagALR', 100*(1-al),'.eps'))
#        plot = plotDiagLiang,
#        filename = paste0('figures_slides/isingDiagDMH', 100*(1-al),'.eps'))

ggsave(width = 6.2, height = 1.9,
       plot = plotDiagLiang,
       filename = paste0('figures_slides/isingDiagDMH', 100*(1-al),'_2.eps'))




# ------------------------------------------------------------------------------
# Combine plots for ACD only
# ------------------------------------------------------------------------------

plotDiagAtchade = ggarrange(plotESSAtchade, plotACDAtchade, plotSampAtchade, nrow = 1, top = grid::textGrob('(a) ALR: Diagnostics and estimated posterior densities', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
# plotDiagAEX = ggarrange(plotESSAEX, plotACDAEX, plotSampAEX, nrow = 1, top = grid::textGrob('(b) AEX: Diagnostics and estimated posterior densities', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))
plotDiagLiang = ggarrange(plotESSLiang, plotACDLiang, plotSampLiang, nrow = 1, top = grid::textGrob('(b) DMH: Diagnostics and estimated posterior densities', gp = gpar(fontsize = 13), x = 0.05, hjust = 0))

# plotDiag = grid.arrange(plotDiagAtchade, plotDiagAEX, plotDiagLiang, ncol = 1)
plotDiag = grid.arrange(plotDiagAtchade, plotDiagLiang, ncol = 1)

ggsave(plot = plotDiag, 
       filename = paste0('figures_slides/isingACD', 100*(1-al),'.eps'), 
       width = 6.5, height = 4)

ggsave(width = 6.5, height = 2,
       plot = plotDiagAtchade,
       filename = paste0('figures_slides/isingACDALR', 100*(1-al),'.eps'))
       # plot = plotDiagLiang,
       # filename = paste0('figures_slides/isingACDDMH', 100*(1-al),'.eps'))




# ==============================================================================
# standard MCMC diagnostics
# ==============================================================================
sampleLiang = data.frame()
for(i in 1:3){
  load(paste0('ACDAIKS/Liang/simLiang', i, '.RData'))
  df = data.frame(m = as.factor(i), parameter = as.vector(Liang[burnin + 1:n]))
  sampleLiang = rbind(sampleLiang, df)
}

sampleLiang %>% 
  group_by(m) %>% 
  summarise(ess = round(ess(parameter)), 
            mcse = round(bm(parameter)[[2]], 5))


load('mcmcDiag/simLiang.RData')
gelman.diag(mcmc.list(as.mcmc(Liang[[1]]), as.mcmc(Liang[[2]]), as.mcmc(Liang[[3]])))[1]
gelman.diag(mcmc.list(as.mcmc(Liang[[4]]), as.mcmc(Liang[[5]]), as.mcmc(Liang[[6]])))[1]
gelman.diag(mcmc.list(as.mcmc(Liang[[7]]), as.mcmc(Liang[[8]]), as.mcmc(Liang[[9]])))[1]

