# loading library
## data handling ##
library(tidyverse)
library(dplyr)

## generating graph ##
library(ggplot2)
library(reshape2)
library(maptools)
library(lattice)
library(rasterVis)
library(gridExtra)
library(ggpubr)


#
data1 <- longterm.stats %>%
  filter(Scenario != 'History', Model !='ACCESS1-0')

# The palette with grey:
cbp1 <- c("#E69F00", "#56B4E9", "#009E73", #"#999999", 
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## mean pr and tasmax
prChange.g <- ggplot(data=data1, aes(x = Country, y = pr_change, fill=Model)) +
  geom_bar(stat='identity',position=position_dodge()) +
  scale_fill_manual(values = cbp1) +
  labs(y = "Annual Mean Precipitation Change (mm)", x = "") +
  facet_grid(Scenario ~ Period) +
  theme_bw()+
  theme(legend.position = 'none',
           axis.text.x = element_text(angle = 45, hjust=1))

tasmaxChange.g <- ggplot(data=data1, aes(x = Country, y = tasmax_change, fill=Model)) +
  geom_bar(stat='identity',position=position_dodge()) +
  scale_fill_manual(values = cbp1) +
  labs(y = "Temperature Change (°C)", x = "") +
  facet_grid(Scenario ~ Period) +
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust=1))

## stddev pr and tasmax
prStdChange.g <- ggplot(data=data1, aes(x = Country, y = prStddev_change, fill=Model)) +
  geom_bar(stat='identity',position=position_dodge()) +
  scale_fill_manual(values = cbp1) +
  labs(y = expression(atop("Change in Interannual","Standard Deviation Precipitation (mm)")), x = "") +
  facet_grid(Scenario ~ Period) +
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust=1))

tasmaxStdChange.g <- ggplot(data=data1, aes(x = Country, y = tasmaxStddev_change, fill=Model)) +
  geom_bar(stat='identity',position=position_dodge()) +
  scale_fill_manual(values = cbp1) +
  labs(y = expression(atop("Change in Temerpature","Standard Deviation (°C)")), x = "") +
  facet_grid(Scenario ~ Period) +
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust=1))

## historical pr and tasmax
data2 <- longterm.stats %>%
  filter(Scenario == 'History')

pr.g <- ggplot(data=data2, aes(x = Country, y = pr_mean), fill='lightgrey') +
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(aes(ymin=pr_mean-pr_stddev, ymax=pr_mean+pr_stddev), width=0.3, size=1.5,
                position=position_dodge()) +
  labs(y = "Annual Precipitation (mm)", x = "") +
  facet_grid(~Period) +
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust=1))

tasmax.g <- ggplot(data=data2, aes(x = Country, y = tasmax_mean), fill='lightgrey') +
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(aes(ymin=tasmax_mean-tasmax_stddev, ymax=tasmax_mean+tasmax_stddev), width=0.3, size=1.5,
                position=position_dodge()) +
  labs(y = "Temerature (°C)", x = "") +
  facet_grid(~Period) +
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust=1))

ggarrange(
  ggarrange(pr.g, ncol=1,nrow=3,labels=c('A'),heights = c(1,1.2,1.1,1.1)), # First column with pr
  ggarrange(prChange.g, prStdChange.g, ncol = 1, nrow=2, align = "v", labels = c("B", "C")), # Second column with change mean and std
  ncol=2,nrow=1, widths = c(1,2.5)
  )
ggarrange(
  ggarrange(tasmax.g, ncol=1,nrow=3,labels=c('A'),heights = c(1,1.2,1.1,1.1)), # First column with tasmax
  ggarrange(tasmaxChange.g, tasmaxStdChange.g, ncol = 1, nrow=2, align = "v", labels = c("B", "C")), # Second column with change mean and std
  ncol=2,nrow=1, widths = c(1,2.5)
)
