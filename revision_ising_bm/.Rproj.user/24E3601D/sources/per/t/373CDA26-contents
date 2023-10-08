rm(list=ls())
library(dplyr); library(ggplot2); library(egg); library(gridExtra); library(grid); library(scales)
library(batchmeans)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

cite = read.csv('cite/sumCite.csv')
cite$year = as.factor(cite$year)

cite.plot = cite %>% 
  filter(year %in% as.factor(2005:2022)) %>% 
  ggplot(aes(x = year, y = cite)) +
  geom_bar(stat="identity") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = 'Year', y = 'Number of citations')

ggsave(plot = cite.plot, width = 5, height = 3, file = 'cite/cite.eps')
