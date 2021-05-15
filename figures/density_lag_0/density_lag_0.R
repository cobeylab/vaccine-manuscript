library(RSQLite)
library(ggplot2)
library(sensitivity)
library(gridExtra)
library(grid)
library(scales)
library(viridis)
library(plyr)
library(cowplot)
library(tidyverse)

textSize = 11
pointSize = 1.0
lineSize = 1
plot_themes  = 	theme_classic() +
  theme(axis.line = element_line(size=1)) +
  theme(axis.ticks = element_line(size=0.5)) +
  theme(axis.ticks.length = unit(-0.1,'cm')) +
  theme(axis.title.x=element_text(size=textSize)) + 
  theme(axis.text.x=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) + 
  theme(axis.title.y=element_text(size= textSize)) +
  theme(axis.text.y=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) +
  theme(plot.title=element_text(size=textSize+2)) +
  theme(plot.margin=unit(c(4,4,0,4),'mm')) +
  theme(legend.title=element_text(size=textSize)) +
  theme(legend.text=element_text(size=textSize)) +
  theme(legend.position ='bottom') +
  theme(legend.direction='horizontal') +
  theme(legend.box.margin = margin(0,0,0,0,'mm')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.line = element_blank())

percentile = function(percentile){
  function(x) quantile(x, percentile, na.rm=T)
}

makeDriftDensity = function(vaccineDF,nbreaks=41){
  rates = sort(unique(vaccineDF$vaccinationRate))
  dens = data.frame(matrix(nrow = length(rates) * (nbreaks-1), ncol=4))
  names(dens) = c('vaccinationRate','cumulativeDrift','density','mean')
  for(i in 1:length(rates)){
    rate = rates[i]
    breaks = seq(0,35,length= nbreaks)
    subDF = vaccineDF[vaccineDF$vaccinationRate == rate,]
    drift.hist = hist(subDF$cumulativeDrift, breaks=breaks,xlim=c(0,c(max(vaccineDF$cumulativeDrift))), plot=FALSE)
    drift.dens = drift.hist$counts/(sum(drift.hist$counts))
    tmp = cbind(rep(rate,length(breaks)-1),drift.hist$breaks[1:length(breaks)-1],drift.dens, mean(subDF$cumulativeDrift[subDF$fluLike==1]))
    start.index = ((i-1)*(length(breaks)-1))+1
    dens[start.index:(start.index+length(breaks)-2),] = tmp
  }
  return(dens)
}

makeIncDensity = function(vaccineDF, nbreaks=41){
  rates = sort(unique(vaccineDF$vaccinationRate))
  dens = data.frame(matrix(nrow = length(rates) * (nbreaks-1), ncol=4))
  names(dens) = c('vaccinationRate','cumulativeIncidence','density','mean')
  for(i in 1:length(rates)){
    rate = rates[i]
    breaks = seq(0,3.2,length=nbreaks)
    subDF = vaccineDF[vaccineDF$vaccinationRate == rate,]
    drift.hist = hist(subDF$cumulativeIncidence, breaks=breaks,xlim=c(0,c(max(vaccineDF$cumulativeIncidence))), plot=FALSE)
    drift.dens = drift.hist$counts/(sum(drift.hist$counts))
    tmp = cbind(rep(rate,length(breaks)-1),drift.hist$breaks[1:length(breaks)-1],drift.dens,mean(subDF$cumulativeIncidence[subDF$fluLike==1]))
    start.index = ((i-1)*(length(breaks)-1))+1
    dens[start.index:(start.index+length(breaks)-2),] = tmp
  }
  return(dens)
}

makeDensityPlot = function(driftDens, incDens, plotName, meanDF, vaccineDF){
  
  plot1 = ggplot(data = driftDens, aes(x=vaccinationRate,y=cumulativeDrift)) + 
    xlab('Vaccination coverage') +
    ylab('Cumulative antigenic evolution') +
    geom_tile(aes(fill=density)) + 
    stat_summary(data = vaccineDF, aes(y = cumulativeDrift),
                 fun.data = "mean_cl_boot",
                 geom = 'errorbar',
                 size = .5,
                 color = 'white') + 
    geom_point(data = vaccineDF, stat = 'summary',
               fun.y = 'mean',
               color = 'white', size = .5) + 
    # geom_errorbar(data = vaccineDF, stat = 'summary',
    #               fun.ymin = percentile(.05),
    #               fun.ymax = percentile(.95),
    #               width = 0,
    #               color = 'white', size = .3) +
    scale_fill_viridis(option='plasma', limits = c(0,1)) +
    guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    #geom_smooth(data=vaccineDF, aes(x=vaccinationRate,y=cumulativeDrift), color='white', size=0.6, se=TRUE, alpha=0.7) +
    plot_themes
  
  plot2 = ggplot(data = incDens, aes(x=vaccinationRate,y=cumulativeIncidence)) + 
    xlab('Vaccination coverage') +
    ylab('Cumulative incidence') +
    geom_tile(aes(fill=density)) + 
    stat_summary(data = vaccineDF, aes(y = cumulativeIncidence),
                 fun.data = "mean_cl_boot",
                 geom = 'errorbar',
                 size = .5,
                 color = 'white') + 
    geom_point(data = vaccineDF, stat = 'summary',
               fun.y = 'mean',
               color = 'white', size = .5) + 
    # geom_errorbar(data = vaccineDF, stat = 'summary',
    #               fun.ymin = percentile(.05),
    #               fun.ymax = percentile(.95),
    #               width = 0,
    #               color = 'white', size = .3) +
    scale_fill_viridis(option='plasma', limits = c(0,1)) +
    guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    #geom_smooth(data=vaccineDF, aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.6, se=TRUE, alpha=0.7) +
    #geom_smooth(data=vaccineDF[vaccineDF$fluLike==1,], aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.6, se=TRUE, alpha=0.7, linetype='dotted') +
    #geom_smooth(data=vaccineDF[vaccineDF$extinct==1,], aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.6, se=TRUE, alpha=0.7, linetype='dashed') +
    plot_themes
  
  legend_b <- get_legend(
    plot1 + 
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  
  plotrow = plot_grid(plot1 + theme(legend.position = 'none'), 
                      plot2 + theme(legend.position = 'none'),
                      labels = c('A','B'), align = 'h', nrow = 1, ncol = 2, hjust = -1)
  
  outplot = plot_grid(plotrow, legend_b, ncol=1, rel_heights = c(1,.3))
  return(outplot)
}


resultsDb = '../../analysis/lag_0_b1/vaccine/results.sqlite'
comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)

vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')

baseLineFlux = vaccineDF$meanFluxRate[vaccineDF$vaccinationRate==0 & vaccineDF$fluLike==1]
vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000

vaccineDF$acceleration = vaccineDF$meanFluxRate - mean(baseLineFlux)

driftDens = makeDriftDensity(vaccineDF[vaccineDF$tmrcaLimit==0,])
incDens = makeIncDensity(vaccineDF[vaccineDF$tmrcaLimit==0,])
#extinctDens = makeExtinctDensity(vaccineDF[vaccineDF$extinct==1,])

summaryDF = ddply(vaccineDF[vaccineDF$tmrcaLimit==0,], .(vaccinationRate), summarise, 
                  mean.drift = mean(cumulativeDrift), 
                  mean.incidence = mean(cumulativeIncidence),
                  var.drift = var(cumulativeDrift),
                  var.incidence = var(cumulativeIncidence),
                  extinct = sum(extinct),
                  tmrcaLimit = sum(tmrcaLimit),
                  fluLike = sum(fluLike))

plotName = 'density_lag_0'
outplot = makeDensityPlot(driftDens,incDens,plotName,summaryDF, vaccineDF[vaccineDF$tmrcaLimit == 0,])
save_plot(paste(plotName,'.pdf',sep=''), outplot, base_aspect_ratio = 1.8)


cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeDrift, method = 'spearman')
cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeIncidence, method = 'spearman')

head(vaccineDF)

vaccineDF %>% filter(tmrcaLimit == 0) %>%
  group_by(vaccineLag, vaccinationRate) %>%
  summarise(drift = mean(cumulativeDrift), sddrift = sd(cumulativeDrift),
            inc = mean(cumulativeIncidence), sdinc = sd(cumulativeIncidence)) %>%
  filter(vaccinationRate == "0.05")