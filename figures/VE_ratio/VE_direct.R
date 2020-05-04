library(tidyverse)
library(viridis)
library(cowplot)

textSize = 12
plot_themes  =	theme_classic() +
  theme(axis.ticks = element_line(size=.5, color = 'black'),
        axis.ticks.length = unit(-4,'pt'),
        axis.text.x=element_text(size = textSize, color='black', margin=margin(7,7,0,7,'pt')), 
        axis.text.y=element_text(size = textSize, color='black', margin=margin(7,7,7,0,'pt')),
        axis.title=element_text(size= textSize, color='black'),
        plot.title=element_text(size=textSize+2, color='black'),
        plot.margin=unit(c(9,9,9,9),'pt'),
        legend.margin=margin(l = 8, unit='pt'),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()
  )

#for reference
WRITESAMPLERATEI = 1
WRITESAMPLERATES = 1

get_data = function(filename, simtype){
  raw_data = read.csv(filename) 
  data = raw_data %>% mutate(I = ifelse(is.na(I), 0, I),
                             V = ifelse(is.na(V), 0, V))
  data = data %>% mutate(runId = paste0(runId,simtype)) %>%
    #group_by(runId, time) %>%
    group_by(runId) %>%
    mutate(VI = (V==1 & I==1),
           VU = (V==1 & I==0),
           NI = (V==0 & I==1),
           NU = (V==0 & I==0)) %>% 
    #filter(extinctNow == 0) %>%
    summarise(VI = sum(VI)*1/WRITESAMPLERATEI,
              VU = sum(VU)*1/WRITESAMPLERATES,
              NI = sum(NI)*1/WRITESAMPLERATEI,
              NU = sum(NU)*1/WRITESAMPLERATES,
              rate = mean(VR.init)) %>%
    mutate(VE = 1-(VI/(VU))/(NI/(NU))) %>%
    mutate(VE = ifelse(is.nan(VE), 0, VE)) %>%
    mutate(ARR = NI/(NI+NU) - VI/(VI+VU)) %>%
    ungroup() %>% 
    mutate(type = simtype)
  
  
  return(data)
}

summarize_data = function(data){
  data_summ = data %>%
    group_by(runId) %>%
    summarise(VE_arithmetic = mean(VE, na.rm=T), 
              VE = exp(mean(log(VE), na.rm=T)), 
              ARR = mean(ARR),
              rate = mean(rate),
              type = unique(type),
              n=n())
  return(data_summ)
}

read_data = function(){
  data_orig = rbind(get_data('subsampled_10percent/dynamic_10percent_breadth_1.csv','dynamic'), 
                    get_data('subsampled_10percent/static_10percent_breadth_1.csv', 'static'))
  
  return(data_orig)
}

data_orig = read_data()


data = data_orig %>% filter(rate < .2) %>%
  mutate(rate = factor(rate, levels= c(0, .01, .03, .05, .07, .1))) 
# data_summ = summarize_data(data) %>% 
#   mutate(rate = factor(rate)) %>% 
#   mutate(VE_arithmetic = ifelse(is.infinite(VE_arithmetic), 0, VE_arithmetic))

data_nozero = data %>% filter(rate!=0)

plotdata=data_nozero %>%
  group_by(rate, type) %>%
  do(data.frame(rbind(smean.cl.boot(.$VE)))) %>%
  rename(VE = Mean,
         lowerCI = Lower,
         upperCI = Upper)
pd = position_dodge(width = 0.5)
private = ggplot(plotdata, 
                 aes(x=rate, y = VE, color = type)) +
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.2, position=pd)+
  coord_flip() + ylim(-1,1) +
  plot_themes + scale_color_brewer(palette = 'Dark2') + xlab("Vaccination rate") + ylab('Private benefit (1 - Odds ratio)')

# private = ggplot(data_nozero, 
#                  aes(x=rate, y = VE, color = type)) +
#   stat_summary(fun.data=mean_cl_boot, position=pd) +
#   coord_flip() + ylim(0,1) +
#   plot_themes + scale_color_brewer(palette = 'Dark2') + xlab("Vaccination rate") + ylab('Private benefit (1 - Odds ratio)')
# 

private_ARR = ggplot(data %>% filter(as.numeric(as.character(rate))> 0), aes(x=rate, y = ARR, color = type)) + 
  stat_summary(fun.data = 'mean_cl_boot', aes(color = type), position=position_dodge(width=.5)) +
  coord_flip()+
  plot_themes + scale_color_brewer(palette = 'Dark2') + xlab("Vaccination rate") + ylab('Private benefit (ARR)')


socialdata = data %>% 
  group_by(runId) %>%
  summarise(I = mean(NI),#sum(VI + NI),
            U = mean(NU),#sum(VU + NU),
            rate = unique(rate),
            type = unique(type)) %>%
  mutate(OR = I/U) 

zerovac = socialdata %>% filter(rate==0) %>%
  summarise(meanOR = mean(I/U),
            meanrisk = mean(I/(I+U)))
socialdata = socialdata %>%
  mutate(benefit = 1- OR/zerovac$meanOR,
         ARR = zerovac$meanrisk - I/(I+U))
social_ARR = ggplot(socialdata %>% filter(rate != 0) %>% mutate(rate = factor(rate)), aes(x=rate, y = ARR, color = type)) + 
  stat_summary(fun.data = 'mean_cl_boot', aes(color = type), position=position_dodge(width=.5)) +
  coord_flip()+
  plot_themes + scale_color_brewer(palette = 'Dark2') + xlab("Vaccination rate") + ylab('Social benefit (ARR)')



plotdata=socialdata %>%
  group_by(rate, type) %>%
  do(data.frame(rbind(smean.cl.boot(.$benefit)))) %>%
  rename(benefit = Mean,
         lowerCI = Lower,
         upperCI = Upper)
  
# social = ggplot(socialdata %>% filter(rate != 0) %>% mutate(rate = factor(rate)), aes(x=rate, y = benefit, color = type)) + 
#   stat_summary(fun.data = 'mean_cl_boot', aes(color = type), position=position_dodge(width=.5)) +
#   coord_flip()+ ylim(0,1.5) +
#   plot_themes + scale_color_brewer(palette = 'Dark2') + xlab("Vaccination rate") + ylab('Social benefit (1 - Odds ratio)')
# 

social = ggplot(plotdata %>% filter(rate != 0) , aes(x=rate, y = benefit, color = type)) + 
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.2, position=pd)+
  coord_flip()+ ylim(-.1,1) +
  plot_themes + scale_color_brewer(palette = 'Dark2') + xlab("Vaccination rate") + ylab('Social benefit (1 - Odds ratio)')



legend = get_legend(social + theme(legend.position ='bottom'))
outrow = plot_grid(social + theme(legend.position = 'none'), 
                   private +  theme(legend.position = 'none'),
                   labels=c('A', 'B',''), nrow = 2)
outplot_OR = plot_grid(outrow, legend, ncol=1, rel_heights = c(1,.04))

#outplot_OR = plot_grid(social, private, labels='AUTO', nrow = 2)
outplot_ARR = plot_grid(social_ARR, private_ARR, labels='AUTO', nrow = 2)

save_plot("VE_OR_fromincidence_021520.pdf", outplot_OR, nrow=2, base_aspect_ratio = 1.5)
#save_plot("VE_ARR_fromincidence_021420.pdf", outplot_ARR, nrow=2, base_aspect_ratio = 1.5)
