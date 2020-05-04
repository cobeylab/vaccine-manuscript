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
    gather(key='class', value='frequency', 2:4 ) %>%
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

plot.timeseries = function(ts, threshold=0, rate=0, R0=0, extinctDF, equilibria){
  #transform I
  scalefactorI = 1/max(ts$frequency[ts$class=='I'])
  ts$frequency[ts$class=='I'] = ts$frequency[ts$class=='I'] *scalefactorI
  equilibria$frequency[equilibria$class =='I'] = equilibria$frequency[equilibria$class =='I'] * 5e7 *scalefactorI
  equilibria$date = ifelse(equilibria$class=='I', 20.4, -.4)
  end.date = max(ts$date)
  plot = ggplot(ts, aes(x=date, y=frequency)) + 
    geom_line(alpha=.18, aes(color=class, group= interaction(runId, class))) + 
    #geom_point(data = equilibria, alpha=1, aes(color = class), size=1.2) +
    geom_hline(data = equilibria, lty=2, alpha=1, aes(yintercept=frequency, color = class)) + 
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*1/scalefactorI, name = 'Number of infecteds')) +
    xlim(c(0,20)) +
    coord_cartesian(xlim = c(-1, 21), ylim=c(-.025, 1), expand=FALSE) +
    geom_point(data = extinctDF, color='red', size=1, alpha=.6) +
    #ylim(c(0,1)) + 
    ylab('Frequency') +
    xlab('Date (year)') +
    guides(colour = guide_legend(override.aes = list(alpha=1))) +
    ggtitle(paste('Annual vaccination rate = ',(round(rate,3))))
#    ggtitle(bquote(paste(italic('R'[0]),' = ', .(round(R0,2)), '; rate = ',.(round(rate,3)))))
    
    #geom_vline(xintercept = end.date) + 
  return(plot)
}

theoretical.threshold = function(beta,gamma, mu){
  left = -(gamma+2*mu)/2
  right = (gamma^2 + 4*mu*beta)^.5/2
  
  p = left+right
  return(p)
  
}

vaccineDF = get.results('results.sqlite')

mu = 1/(30*365)
gamma = 0.2

vaccineDF  = vaccineDF %>%
  mutate(R0 = round(beta/(gamma + 1/(30*365)),2),
         expected.threshold = (gamma+mu+vaccinationRate/365)/beta) 

#RATES = unique(vaccineDF$vaccinationRate)

get.equilibria = function(beta, gamma, mu, p){
  R0 = beta/(gamma+mu+p)
  S_eq = 1/R0
  I_eq = mu/beta*(R0-1) - p/beta
  R_eq = gamma/(mu+p) * I_eq
  
  if(I_eq < 0){
    I_eq = 0
    S_eq = mu/(mu+p)
    R_eq = 0
  }
  V_eq = 1- S_eq - R_eq - I_eq #p/(mu+p)
  out = data.frame(S=S_eq, I=I_eq,R=R_eq, V=V_eq) %>%
    gather(key='class', value='frequency')
  
  return(out)
}

plotindex=1
plots = vector('list', length(unique(vaccineDF$R0))*4)
counter = 1
for(thisR0 in c(1.8)){
  subDF = vaccineDF %>% filter(R0==thisR0)
  RATES = unique(subDF$vaccinationRate)
  for(rate in RATES){
    subSubDF = subDF %>% filter(vaccinationRate == rate)
    
    tslist = vector('list', nrow(subSubDF))
    i=1
    for(runId in subSubDF$runId){
      if(runId != 99){
        print(counter/nrow(vaccineDF))
        SRVtimeseries = get.timeseries(runId)
        Itimeseries = get.its(runId)
        last.date = max(Itimeseries$date)
        subDF$lastDate[subDF$runId == runId] = last.date
        tslist[[i]] = rbind(get.timeseries(runId), get.its(runId))
        i=i+1
        counter=counter+1
      }
    }
    ts = do.call(rbind, tslist)
    beta = unique(subDF$beta)
    susceptible.threshold = (gamma+mu+rate/365)/beta
    extinct.sims = subDF %>% filter(extinct==1)
    extinctDF = data.frame(frequency = rep(0,nrow(extinct.sims)), date= extinct.sims$lastDate, runId = extinct.sims$runId)
    
    equilibria = get.equilibria(beta, gamma, mu, rate/365)
    print('plotting')
    plot = plot.timeseries(ts, susceptible.threshold, rate, R0 = beta/(gamma+mu), extinctDF, equilibria)
    plots[[plotindex]] = plot
    plotindex = plotindex+1
  }
}


grid = plot_grid(plotlist=plots, nrow=1, ncol=4)
save_plot('out.pdf', grid, nrow=1, ncol=3, base_aspect_ratio=2, base_height=4)

#grid = plot_grid(plotlist = plots)
