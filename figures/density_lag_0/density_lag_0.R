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

plot_themes  = 	theme_classic() +
				theme(axis.line = element_line(size=1)) +
				theme(axis.ticks = element_line(size=0.5)) +
				theme(axis.ticks.length = unit(-0.1,'cm')) +
	#			theme(axis.ticks.margin = unit(0.5,'cm')) +
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
				xlab('Vaccination rate') +
				ylab('Cumulative antigenic evolution') +
				geom_tile(aes(fill=density)) + 
				scale_fill_viridis(option='plasma') +
				guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
				scale_x_continuous(expand = c(0,0)) +
				scale_y_continuous(expand = c(0,0)) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
				theme(axis.ticks = element_line(size=.1))  +
# 				geom_line(aes(y=mean),color = 'white',size=.5, linetype='dotted') +
# 				geom_line(data = meanDF, aes(y=mean.drift),color = 'white',size=.5) +
	      geom_smooth(data=vaccineDF, aes(x=vaccinationRate,y=cumulativeDrift), color='white', size=0.3, se=TRUE, alpha=0.7,) +
	      geom_smooth(data=vaccineDF[vaccineDF$fluLike==1,], aes(x=vaccinationRate,y=cumulativeDrift), color='white', size=0.3, se=TRUE, alpha=0.7, linetype='dotted') +
				plot_themes

	plot2 = ggplot(data = incDens, aes(x=vaccinationRate,y=cumulativeIncidence)) + 
				xlab('Vaccination rate') +
				ylab('Cumulative incidence') +
				geom_tile(aes(fill=density)) + 
				scale_fill_viridis(option='plasma') +
	      guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
	      scale_x_continuous(expand = c(0,0)) +
				scale_y_continuous(expand = c(0,0)) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
				theme(axis.ticks = element_line(size=.1))  +
     	  geom_smooth(data=vaccineDF, aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.3, se=TRUE, alpha=0.7) +
    	  geom_smooth(data=vaccineDF[vaccineDF$fluLike==1,], aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.3, se=TRUE, alpha=0.7, linetype='dotted') +
	      plot_themes

	plot = plot_grid(plot1, plot2, labels = c('A','B'), align = 'h', ncol = 2)
	save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 1, base_aspect_ratio = 0.85)
}

resultsDb = '../../analysis/lag_0_breadth_1_density/vaccine/results.sqlite'
comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)

vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')

baseLineFlux = vaccineDF$meanFluxRate[vaccineDF$vaccinationRate==0 & vaccineDF$fluLike==1]
vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000

vaccineDF$acceleration = vaccineDF$meanFluxRate - mean(baseLineFlux)

vaccineDF$vaccinationRate = vaccineDF$vaccinationRate*365

driftDens = makeDriftDensity(vaccineDF[vaccineDF$tmrcaLimit==0,])
incDens = makeIncDensity(vaccineDF[vaccineDF$tmrcaLimit==0,])
extinctDens = makeExtinctDensity(vaccineDF[vaccineDF$extinct==1,])

summaryDF = ddply(vaccineDF[vaccineDF$tmrcaLimit==0,], .(vaccinationRate), summarise, 
                  mean.drift = mean(cumulativeDrift), 
                  mean.incidence = mean(cumulativeIncidence),
                  var.drift = var(cumulativeDrift),
                  var.incidence = var(cumulativeIncidence),
                  extinct = sum(extinct),
                  tmrcaLimit = sum(tmrcaLimit),
                  fluLike = sum(fluLike))

makeDensityPlot(driftDens,incDens,'density_lag_0',summaryDF, vaccineDF[vaccineDF$tmrcaLimit == 0,])

cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeDrift, method = 'spearman')
cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeIncidence, method = 'spearman')

