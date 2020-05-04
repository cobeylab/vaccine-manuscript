library(RSQLite)
library(ggplot2)
library(sensitivity)
library(gridExtra)
library(grid)
library(scales)
library(viridis)
library(plyr)
library(cowplot)

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

makeExtinctDensity = function(vaccineDF, nbreaks=41){
  rates = sort(unique(vaccineDF$vaccinationRate))
  dens = data.frame(matrix(nrow = length(rates) * (nbreaks-1), ncol=4))
  names(dens) = c('vaccinationRate','lastDate','density','mean')
  for(i in 1:length(rates)){
    rate = rates[i]
    breaks = seq(0,21,length=nbreaks)
    subDF = vaccineDF[vaccineDF$vaccinationRate == rate,]
    drift.hist = hist(subDF$lastDate, breaks=breaks,xlim=c(0,c(max(vaccineDF$lastDate))), plot=FALSE)
    drift.dens = drift.hist$counts/(sum(drift.hist$counts))
    tmp = cbind(rep(rate,length(breaks)-1),drift.hist$breaks[1:length(breaks)-1],drift.dens,mean(subDF$lastDate))
    start.index = ((i-1)*(length(breaks)-1))+1
    dens[start.index:(start.index+length(breaks)-2),] = tmp
  }
  return(dens)
}

make.extinct.plot = function(extinctDens, summaryDF,plotName, vaccineDF){
  extinct.plot = ggplot(summaryDF, aes(x=vaccinationRate, y=extinct/500)) +
    xlab('Vaccination rate') +
    ylab('Fraction extinct') +
    geom_point(size=0.5) +
    plot_themes
  
  time.plot = ggplot(data = extinctDens, aes(x=vaccinationRate,y=lastDate)) + 
    xlab('Vaccination rate') +
    ylab('Time to extinction') +
    geom_tile(aes(fill=density)) + 
    scale_fill_viridis(option='plasma') +
    guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    stat_summary(data = vaccineDF, 
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
    #geom_line(aes(y=mean),color = 'white',size=.5) +
    plot_themes
  
  plot = plot_grid(extinct.plot, time.plot, labels = c('A','B'), align = 'h', axis = 'bt', ncol = 2)
  save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 1, base_aspect_ratio = .9)
}

resultsDir = '../../analysis/breadth_1_density/vaccine/'

resultsDb = paste(resultsDir,'results.sqlite',sep='')
comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)

vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')

baseLineFlux = vaccineDF$meanFluxRate[vaccineDF$vaccinationRate==0 & vaccineDF$fluLike==1]
vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000

vaccineDF$acceleration = vaccineDF$meanFluxRate - mean(baseLineFlux)

vaccineDF$vaccinationRate = vaccineDF$vaccinationRate*365

extinctDens = makeExtinctDensity(vaccineDF[vaccineDF$extinct==1,])

summaryDF = ddply(vaccineDF[vaccineDF$tmrcaLimit==0,], .(vaccinationRate), summarise, 
                  mean.drift = mean(cumulativeDrift), 
                  mean.incidence = mean(cumulativeIncidence),
                  var.drift = var(cumulativeDrift),
                  var.incidence = var(cumulativeIncidence),
                  extinct = sum(extinct),
                  tmrcaLimit = sum(tmrcaLimit),
                  fluLike = sum(fluLike))

make.extinct.plot(extinctDens, summaryDF, 'extinct_breadth_1', vaccineDF)


cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeDrift, method = 'spearman')
cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeIncidence, method = 'spearman')

