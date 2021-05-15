library(RSQLite)
library(plyr)
library(doParallel)
library(foreach)
library(ggplot2)
library(viridis)
library(cowplot)
library(tidyverse)
library(patchwork)

textSize = 12
pointSize = .5
lineSize = 1
plotDirectory='./plots'
plot_themes  = 	theme_classic() +
  theme(axis.line = element_line(size=1)) +
  theme(axis.ticks = element_line(size=0.5)) +
  theme(axis.ticks.length = unit(-0.1,'cm')) +
  theme(axis.title.x=element_text(size=textSize)) + 
  theme(axis.text.x=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) + 
  theme(axis.title.y=element_text(size= textSize)) +
  theme(axis.text.y=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) +
  theme(plot.title=element_text(size=textSize+2)) +
  theme(plot.margin=unit(c(5,5,5,5),'mm')) +
  theme(legend.title=element_text(size=textSize)) +
  theme(legend.text=element_text(size=textSize)) +
  theme(legend.position ='bottom') +
  theme(legend.direction='horizontal') +
  theme(legend.margin = unit(0,'cm')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.line = element_blank())

sampleRate = 1e-5


compute_coverage = function(vacrate,  p_revac, n, fraction_never_vac){
  print(vacrate)
  vac_ts = data.frame(year=NA,coverage=NA, rate=NA)
  pop = 1:round((1-fraction_never_vac)*n)
  
  vaccinated_total = c()
  
  dat = data.frame(id = pop)
  df2 = data.frame(matrix(ncol=20))
  names(df2) = paste0('year',1:20)
  dat = cbind(dat,df2)
  
  for(year in 1:20){
    if(year==1){
      vac_this_season = sample(pop, vacrate*length(pop))
    }
    
    else{
      n_vac = vacrate*length(pop)
      repeat_vac = sample(vac_last_season, p_revac*length(vac_last_season))
      nonrepeat_vac = sample(not_vac_last_season, n_vac - length(repeat_vac))
      vac_this_season = c(repeat_vac, nonrepeat_vac)
      
      #print(paste0(length(repeat_vac)/length(pop),':',length(nonrepeat_vac)/length(pop)))
      
    }
    #print(length(vac_this_season))
    vaccinated_total = c(vac_this_season, vaccinated_total)
    vac_last_season = vac_this_season
    not_vac_last_season = c(setdiff(vac_last_season, pop), (setdiff(pop, vac_last_season)))
    
    dat[dat$id %in% vac_this_season, paste0('year',year)] = 1
  }
  
  datlong = dat %>% pivot_longer(cols = year1:year20, values_to = 'vaccinated', names_to='year') %>%
    replace_na(list(vaccinated = 0)) %>%
    mutate(year = str_extract(year,'\\d+') %>% as.numeric())
  
  
  #datlong %>% group_by(year) %>%
  # summarise(vaccinated = sum(vaccinated))
  
  for(i in 1:length(unique(datlong$year))){
    tmp = datlong %>%
      filter(year <= datlong$year[i])
    coverage = tmp %>% group_by(id) %>%
      summarise(vac = max(vaccinated)) %>%
      pull(vac) %>%
      sum()
    coverage= coverage/n
    vac_ts[i,] = c(i, coverage, vacrate)
  }
  
  return(vac_ts)
}

get.distance = function(x1,y1, x2,y2){
  return(sqrt((x1-x2)^2 + (y1-y2)^2))
}

get.cumulative.vac = function(runId, rate, sampleRate, resultsDir){
  fn = paste(resultsDir,'results/',runId,'/samples.sqlite',sep='')
  df = dbConnect(SQLite(), dbname = fn)
  initExtension(df)
  
  hosts = dbGetQuery(df, 'SELECT * FROM hosts WHERE date == cast(date as int)')
  vaccines = dbGetQuery(df, 'SELECT date, idealA, idealB FROM vaccines')
  names(vaccines) = c('date','vac1','vac2')
  hosts$sampleHosts = sampleRate*hosts$S
  hosts = merge(hosts, vaccines,by='date')
  hosts$immunity = pmax(0,1 - get.distance(hosts$ag1,hosts$ag2, hosts$vac1, hosts$vac2) *.07)/hosts$sampleHosts
  hosts$vaccinated = 1/hosts$sampleHosts
  
  if(nrow(hosts)  == 0){
    tmp = data.frame(matrix(ncol=5, nrow=0))
    names(tmp) = c('date','nvac','nimm', 'rate', 'runId')
  }else{
    h2 = aggregate(immunity~hostId + date +sampleHosts + vaccinated, hosts, max)
    
    
    tmp = ddply(h2, .(date), summarise, nvac = sum(vaccinated), nimm = sum(immunity))
    
    tmp$rate = rate
    tmp$runId = runId
    
  }

  return(tmp)
  
}

resultsDir = '../../analysis/vaccine_coverage/vaccine/'
resultsDir = './'

resultsdf = dbConnect(SQLite(), dbname = paste(resultsDir,'results.sqlite',sep=''))
initExtension(resultsdf)

vaccineDF = dbGetQuery(resultsdf, 'SELECT * FROM pooled_results')
vaccineDF = vaccineDF[vaccineDF$tmrcaLimit == 0,]
RUNIDS = vaccineDF$runId
RATES = unique(vaccineDF$vaccinationRate)

cumulative.vac = data.frame(matrix(ncol=5, nrow=0))
names(cumulative.vac) = c('date','nvac','nimm','rate','runId')

cl=makeCluster(8)
registerDoParallel(cl)

for(rate in RATES){
  RUNIDS = vaccineDF$runId[vaccineDF$vaccinationRate == rate]

  tmp = foreach(runId = RUNIDS, .combine = rbind, .packages = c('RSQLite','plyr'))%dopar%{
  
    get.cumulative.vac(runId, rate, sampleRate,resultsDir)
  }
  cumulative.vac = rbind(cumulative.vac,tmp)
}



expected_coverage = (c(1,5, 10, 20, 30, 50)/100) %>% 
  map_dfr(compute_coverage, p_revac=.8, n=100000, fraction_never_vac= 0.33)

summaryDF = ddply(cumulative.vac, .(date, rate), summarise,
                  coverage = mean(nvac, na.rm=T),
                  immunity = mean(nimm, na.rm=T))

expected_coverage = expected_coverage %>%
  dplyr::rename(expected_coverage = coverage, date= year)

summaryDF2 = full_join(summaryDF, expected_coverage)

test = expected_coverage %>% group_by(rate) %>%
  summarise(expected_coverage = max(expected_coverage))

coverage = ggplot(data=cumulative.vac, aes(x=date, y=nvac)) + 
  geom_point(aes(color = factor(rate), group=factor(rate)), size=.5)  +
  geom_segment(data = test,aes(color = factor(rate), group = factor(rate), y=expected_coverage, yend = expected_coverage, x = 20, xend =22)) +
#  geom_hline(data = test,aes(color = factor(rate), group = factor(rate), yintercept=expected_coverage), lty=2, size=.5) +
  guides(colour= guide_legend('Annual vaccination coverage',title.position='top')) +
  xlab('Year') + 
  ylab('Fraction vaccinated\nat least once') + ylim(c(0,1))+ plot_themes + 
  scale_color_brewer(palette = 'Dark2')

immunity = ggplot(data=cumulative.vac, aes(x=date, y=nimm)) + 
  geom_point(aes(color = factor(rate), group=factor(rate)), size=.5) + 
  guides(colour= guide_legend('Annual vaccination coverage',title.position='top')) +
  geom_hline(aes(yintercept = 1 - 1 / 1.8)) + 
  xlab('Year') + 
  ylab('Effective vaccine immunity') + ylim(c(0,1))+ plot_themes + 
  scale_color_brewer(palette = 'Dark2')

plot = coverage +immunity + plot_annotation(tag_levels = 'A') +plot_layout(guides = "collect") & theme(legend.position = 'bottom')
plotName = 'vaccine_coverage'
save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 1, base_aspect_ratio = 0.85)
