library(RSQLite)
library(ggplot2)
library(sensitivity)
library(gridExtra)
library(grid)
library(scales)
library(viridis)
library(plyr)
library(cowplot)

textSize = 12
pointSize = 1.0
lineSize = 1
plotDirectory='./'
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

makeTMRCAdensity = function(vaccineDF, nbreaks = 41){
  rates = sort(unique(vaccineDF$vaccinationRate))
  dens = data.frame(matrix(nrow = length(rates) * (nbreaks-1), ncol=4))
  names(dens) = c('vaccinationRate','meanTMRCA','density','mean')
  for(i in 1:length(rates)){
    rate = rates[i]
    breaks = seq(0,11,length=nbreaks)
    subDF = vaccineDF[vaccineDF$vaccinationRate == rate,]
    drift.hist = hist(subDF$meanTMRCA, breaks=breaks,xlim=c(0,c(max(vaccineDF$meanTMRCA))), plot=FALSE)
    drift.dens = drift.hist$counts/(sum(drift.hist$counts))
    tmp = cbind(rep(rate,length(breaks)-1),drift.hist$breaks[1:length(breaks)-1],drift.dens,mean(subDF$meanTMRCA))
    start.index = ((i-1)*(length(breaks)-1))+1
    dens[start.index:(start.index+length(breaks)-2),] = tmp
  }
  return(dens)
}

make.diversity.plot = function(summaryDF,vaccineDF,plotName){
  div.plot.TE = ggplot(summaryDF, aes(x=vaccinationRate, y=tmrcaLimit/extinct)) +
    xlab('Vaccination rate') +
    ylab('TMRCA > 10 : extinct ratio') +
    geom_point(size=0.5) +
    plot_themes
  
  div.plot.TF = ggplot(summaryDF, aes(x=vaccinationRate, y=tmrcaLimit/fluLike)) +
    xlab('Vaccination rate') +
    ylab('TMRCA > 10 : flu-like ratio') +
    geom_point(size=0.5) +
    plot_themes
  
  div.plot.meanTMRCA = ggplot(vaccineDF, aes(x=vaccinationRate, y=meanTMRCA)) +
    xlab('Vaccination rate') +
    ylab('Mean TMRCA') +
    geom_tile(aes(fill=density)) + 
    scale_fill_viridis(option='plasma') +
    guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    geom_line(aes(y=mean),color = 'white',size=.5) +
    plot_themes
  
  plot = plot_grid(div.plot.TE, div.plot.TF, div.plot.meanTMRCA, labels = c('A','B','C'), align = 'h', ncol = 3)
  save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=3, nrow = 1, base_aspect_ratio = 1)
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

summaryDF.excess.div = ddply(vaccineDF, .(vaccinationRate), summarise, 
                  mean.drift = mean(cumulativeDrift), 
                  mean.incidence = mean(cumulativeIncidence),
                  var.drift = var(cumulativeDrift),
                  var.incidence = var(cumulativeIncidence),
                  extinct = sum(extinct),
                  tmrcaLimit = sum(tmrcaLimit),
                  fluLike = sum(fluLike))
tmrcaDens = makeTMRCAdensity(vaccineDF[vaccineDF$fluLike ==1, ])

make.diversity.plot(summaryDF.excess.div, tmrcaDens, 'diversity_breadth_1')

cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeDrift, method = 'spearman')
cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeIncidence, method = 'spearman')

