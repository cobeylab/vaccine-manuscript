library(tidyverse)
library(RSQLite)
library(cowplot)

get.timeseries = function(runId){
  filename = paste('results/', runId,'/out.hosts', sep='')
  hosts = read.table(filename, sep = ',', header=T) %>%
    mutate(total = S+R+I+V,
           S = S/total,
           #I = I/total,
           R = R/total,
           V = V/total) %>%
    select(date,S,R,V)
  hosts = hosts %>%
    gather(key='class', value='frequency', 2:3 ) %>%
    mutate(runId = runId) %>%
    select(date, frequency, runId, class)
  
  return(hosts)
}

get.its = function(runId){
  filename = paste('results/', runId,'/output.sqlite', sep='')
  tsdb = dbConnect(SQLite(), dbname = filename)
  initExtension(tsdb)
  
  its = dbGetQuery(tsdb, 'SELECT date, totalI FROM timeseries')
  its$runId = runId
  its = its %>%
    rename('frequency' = 'totalI') %>%
    mutate(class = 'I')
  dbDisconnect(tsdb)
  return(its)
}

get.results = function(filename){
  comboDb = dbConnect(SQLite(), dbname = filename)
  initExtension(comboDb)
  
  statusDF = dbGetQuery(comboDb, 'SELECT * FROM status')
  paramsDF = dbGetQuery(comboDb, 'SELECT beta, vaccinationRate, runId FROM parameters')
  vaccineDF = merge(statusDF, paramsDF)
  
  vaccineDF$vaccinationRate = vaccineDF$vaccinationRate*365
  dbDisconnect(comboDb)
  return(vaccineDF)
}

plot.timeseries = function(ts, threshold=0, rate=0, R0=0){
  #transform I
  scalefactorI = 1/max(ts$frequency[ts$class=='I'])
  ts$frequency[ts$class=='I'] = ts$frequency[ts$class=='I'] *scalefactorI
  
  end.date = max(ts$date)
  plot = ggplot(ts, aes(x=date, y=frequency, color=class, group= interaction(runId, class))) + 
    geom_line(alpha=.1) + 
    geom_hline(yintercept = threshold, lty=2, alpha=1) + 
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*1/scalefactorI, name = 'Prevalence (counts)')) +
    xlim(c(0,20)) +
    #ylim(c(0,1)) + 
    ylab('Frequency') +
    xlab('Date (year)') +
    guides(colour = guide_legend(override.aes = list(alpha=1))) +
    ggtitle(bquote(paste(italic('R'[0]),' = ', .(round(R0,2)), '; rate = ',.(round(rate,3)))))
    
    #geom_vline(xintercept = end.date) + 
  return(plot)
}

vaccineDF = get.results('results.sqlite')

mu = 1/(30*365)
gamma = 0.2

vaccineDF  = vaccineDF %>%
  mutate(R0 = round(beta/(gamma + 1/(30*365)),2),
         expected.threshold = (gamma+mu+vaccinationRate/365)/beta) 

#RATES = unique(vaccineDF$vaccinationRate)

plotindex=1
plots = vector('list', length(unique(vaccineDF$R0))*4)
counter = 1
for(thisR0 in c(1.8)){
  subDF = vaccineDF %>% filter(R0==thisR0)
  RATES = unique(subDF$vaccinationRate)
  for(rate in RATES){
    subSubDF = subDF %>% filter(vaccinationRate == rate)
    
    ts = vector('list', nrow(subSubDF))
    i=1
    for(runId in subSubDF$runId){
      print(counter/nrow(vaccineDF))
      ts[[i]] = rbind(get.timeseries(runId), get.its(runId))
      i=i+1
      counter=counter+1
    }
    ts = do.call(rbind, ts)
    beta = unique(subDF$beta)
    susceptible.threshold = (gamma+mu+rate/365)/beta
    plot = plot.timeseries(ts, susceptible.threshold, rate, R0 = beta/(gamma+mu))
    plots[[plotindex]] = plot
    plotindex = plotindex+1
  }
}


grid = plot_grid(plotlist=plots, nrow=1, ncol=3)
save_plot('out.pdf', grid, nrow=1, ncol=3, base_aspect_ratio=1.5, base_height=4)

#grid = plot_grid(plotlist = plots)
