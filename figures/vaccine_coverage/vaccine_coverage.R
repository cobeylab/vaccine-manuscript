library(RSQLite)
library(plyr)
library(doParallel)
library(foreach)
library(ggplot2)
library(viridis)
library(cowplot)

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

sampleRate = 1e-4

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
resultsdf = dbConnect(SQLite(), dbname = paste(resultsDir,'results.sqlite',sep=''))
initExtension(resultsdf)

vaccineDF = dbGetQuery(resultsdf, 'SELECT * FROM pooled_results')
vaccineDF = vaccineDF[vaccineDF$tmrcaLimit == 0,]
vaccineDF$vaccinationRate = round(vaccineDF$vaccinationRate * 365,2)
RUNIDS = vaccineDF$runId
RATES = unique(vaccineDF$vaccinationRate)

cumulative.vac = data.frame(matrix(ncol=5, nrow=0))
names(cumulative.vac) = c('date','nvac','nimm','rate','runId')

cl=makeCluster(4)
registerDoParallel(cl)

for(rate in RATES){
  RUNIDS = vaccineDF$runId[vaccineDF$vaccinationRate == rate]

  tmp = foreach(runId = RUNIDS, .combine = rbind, .packages = c('RSQLite','plyr'))%dopar%{
  
    get.cumulative.vac(runId, rate, sampleRate,resultsDir)
  }
  cumulative.vac = rbind(cumulative.vac,tmp)
}

summaryDF = ddply(cumulative.vac, .(date, rate), summarise,
                  coverage = mean(nvac, na.rm=T),
                  immunity = mean(nimm, na.rm=T))
coverage = ggplot(data=cumulative.vac, aes(x=date, y=nvac)) + 
  geom_point(aes(color = factor(rate), group=factor(rate)))  +
  guides(colour= guide_legend('Annual vaccination rate',title.position='top')) +
  xlab('Year') + 
  ylab('Fraction vaccinated at least once') + ylim(c(0,1))+ plot_themes

immunity = ggplot(data=cumulative.vac, aes(x=date, y=nimm)) + 
  geom_point(aes(color = factor(rate), group=factor(rate))) + 
  guides(colour= guide_legend('Annual vaccination rate',title.position='top')) +
  geom_hline(aes(yintercept = 1 - 1 / 1.8)) + 
  xlab('Year') + 
  ylab('Effective vaccine immunity') + ylim(c(0,1))+ plot_themes

plot = plot_grid(coverage,immunity,labels = c('A','B'), align='h', ncol=2)
plotName = 'vaccine_coverage'
save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 1, base_aspect_ratio = 0.85)
