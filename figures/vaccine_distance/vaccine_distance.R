library(tidyverse)
library(jsonlite)
library(interp)
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

get_data = function(topdir){
  dirs = paste0(topdir,'/vaccine/results/',0:139,'/')
  
  out = data.frame(matrix(ncol=4, nrow = length(dirs)))
  names(out) = c('runId', 'distance','sd', 'rate')
  
  for(i in 1:length(dirs)){
    print(dir)
    dir = dirs[i]
    
    vaccines = read_tsv(paste0(dir, 'vaccine.timeseries'))
    params = fromJSON(paste0(dir,'parameters.json'), flatten=T)
    
    vaccineLag = params$vaccineLag
    vaccineWindow = params$vaccineWindow
    rate=params$vaccinationRate
    breadth=params$vaccineImmuneBreadth
    
    vaccines = vaccines %>% 
      mutate(vaccineA = lag(idealA),
             vaccineB = lag(idealB)) %>%
      mutate(distance = ((vaccineA-idealA)^2 + (vaccineB-idealB)^2)^.5)
    
    mean_distance = mean(vaccines$distance, na.rm=T)
    sd_distance = sd(vaccines$distance, na.rm=T)
    
    out[i,] = c(dir, mean_distance, sd_distance, rate)
  }
  return(out)
}

compute_distance = function(data){
  out = data %>%
    mutate(distance = as.numeric(distance),
           sd = as.numeric(sd),
           rate = as.numeric(rate)) %>%
    group_by(rate) %>%
    summarise(d = mean(distance),
              n = n(),
              sd_avg = sum(sd^2, na.rm=T)^.5 / sum(!is.na(sd)))
  return(out)
}

dynamic_distances = get_data('../../analysis/panelsim_dynamic') %>%
  compute_distance() %>%
  mutate(type = 'Dynamic')

static_distances = get_data('../../analysis/panelsim_static') %>%
  compute_distance() %>%
  mutate(type = 'Static')

combined_distances = rbind(dynamic_distances,static_distances) %>%
  mutate(rate = (round(rate*365,2))) %>%
  filter(as.numeric(rate) < .15) %>%
  mutate(rate = as.factor(rate)) 
  

pos = position_dodge(width = .5)
outplot = ggplot(combined_distances, aes(x=rate, y = d, color = type)) + 
  geom_point(position = pos) +
  geom_errorbar(position = pos, aes(ymin = d-sd_avg, ymax = d+sd_avg), width = .2) + 
  scale_color_brewer(palette= 'Dark2') +
  xlab('Annual vaccination rate') +
  ylab('Mean vaccine distance\n(Antigenic units)') +
  plot_themes

save_plot('~/Downloads/vaccine_distance.pdf', outplot,base_aspect_ratio = 1.5)
