library(tidyverse)
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

format_VE = function(data){
  VEs = data %>% filter(str_detect(X1, 't-0') | str_detect(X1, 'Social'))
  VEs = VEs %>%
    mutate(rate = str_extract(X1, '..\\%'),
           benefit = substr(X1, 1,7)) %>%
    gather(key = simtype, value = VE, Static_5:Dynamic_1) %>%
    separate(simtype, into = c('simtype','breadth'), sep = '_') %>%
    mutate(benefit = str_trim(benefit),
           rate = str_trim(rate),
           VE = str_replace_all(VE,'\\*',''))
  return(VEs)
}

dat = read_csv('linreg_or_032520.csv')

coefs_raw = dat[2*(1:(nrow(dat)/2))-1,]
confint_raw = dat[2*(1:(nrow(dat)/2)),]
confint_raw$X1 = coefs_raw$X1

VEs = format_VE(coefs_raw) %>% mutate(VE = as.numeric(VE))
CIs = format_VE(confint_raw) %>%
  separate(VE, into = c('lowerCI','upperCI'), sep = ' - ') %>%
  mutate(lowerCI = str_replace(lowerCI,'\\(','') %>% as.numeric(),
         upperCI = str_replace(upperCI,'\\)','') %>% as.numeric())

lowerbound = -1

VE_forplot = merge(VEs,CIs) %>%
  mutate(rate = as.character(rate)) %>%
  mutate(VE = 1-VE, upperCI_tmp = 1-lowerCI, lowerCI_tmp = 1-upperCI) %>%
  select(-upperCI, -lowerCI) %>% 
  rename(upperCI = upperCI_tmp, lowerCI = lowerCI_tmp) %>%
  filter(breadth==1) %>%
  mutate(rate = factor(rate, levels= c( '1%', '3%', '5%', '7%', '10%'))) %>%
  mutate(lower = ifelse(lowerCI < lowerbound, lowerbound, NA),
         label_location = ifelse(lowerCI<lowerbound, lowerbound+.1, NA),
         low_label = ifelse(lowerCI < lowerbound, round(lowerCI,2), NA),
         lowerCI = ifelse(is.na(lower), lowerCI, lower))

pd = position_dodge(width=.5)

privateplot = ggplot(VE_forplot %>% filter(benefit == 'Private'), aes(x=rate,y=VE, color = simtype)) + 
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.2, position=pd) +
  #geom_text(aes(y=label_location, label = low_label), color = 'black', position = position_dodge(width = 1), size = 3) +
  geom_point(position = pd, size = .5)  + 
  ylim(0,1) +
  scale_color_brewer(palette = 'Dark2',name='Type') +
  # geom_segment(aes(xend=rate,y=0, yend=lowerCI, color=simtype, group = interaction(simtype, rate)),
  #   arrow = arrow(type='closed', angle=8, length = unit(0.2,'inches')), 
  #   position = pd, size=1, show.legend= FALSE) +
  coord_flip() +
  ylab("Private benefit (1 - Odds ratio)") + 
  xlab("Annual vaccination rate") + plot_themes

arrowdata = VE_forplot %>% 
  mutate(lower = lower-.05)%>%
  filter(benefit=='Private') %>%
  mutate(arrowstart = ifelse(!is.na(lower), 0, NA)) %>%
  gather(key = type, value=value, arrowstart, lower) %>%
  select(rate, value, simtype, lowerCI) %>%
  arrange(rate, simtype)

arrows = geom_line(data = arrowdata, 
                   aes(x=rate, y = value, group = interaction(rate, simtype)), 
                   arrow=arrow(type='closed', angle=15, length = unit(0.1,'inches')), position=pd,
                   size=0, show.legend= FALSE)

privateplot = privateplot #+ arrows

socialplot = ggplot(VE_forplot %>% filter(benefit == 'Social'), aes(x=rate,y=VE, color = simtype)) + 
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.2, position=pd) +
  geom_point(position = pd, size=.5) + 
  ylim(0,1) +
  scale_color_brewer(palette = 'Dark2', name='Type') +
  coord_flip() +
  ylab("Social benefit (1 - Odds ratio)") + 
  xlab("Annual vaccination rate") + plot_themes

legend = get_legend(socialplot + theme(legend.position ='bottom'))
outrow = plot_grid(socialplot + theme(legend.position = 'none'), 
                   privateplot +  theme(legend.position = 'none'),
                   labels=c('A', 'B',''), nrow = 2)
outplot = plot_grid(outrow, legend, ncol=1, rel_heights = c(1,.04))

save_plot("VE.pdf", outplot, nrow=2, base_aspect_ratio = 1.8, base_height=2.9)

