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
dynamicDir = '../../analysis/breadth_1_density/vaccine/'
staticDir = '../../analysis/static_density/vaccine/'
resultsDirs = c(dynamicDir,staticDir)

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

makePlot = function(summaryDF, fluDF, plotName){
	plot1 = ggplot(data = summaryDF, aes(x=vaccinationRate,y=meanDrift)) + 
				xlab('Vaccination rate') +
				ylab('Cumulative antigenic evolution') +
	      geom_smooth(aes(colour = sim), size=0.5) +
	      scale_color_viridis(discrete=TRUE) +
	      guides(colour= guide_legend('',title.position='top')) +
	      geom_smooth(data = fluDF, aes(colour=factor(sim), y=fluDrift), size=0.5, linetype = 'dashed', show_guide=FALSE) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
				theme(axis.ticks = element_line(size=.1)) +
	      ylim(c(0,30)) +
				plot_themes

	plot2 = ggplot(data = summaryDF, aes(x=vaccinationRate,y=meanInc)) + 
    	  xlab('Vaccination rate') +
    	  ylab('Cumulative incidence') +
    	  geom_smooth(aes(colour = sim), size=0.5) +
	      scale_color_viridis(discrete=TRUE, labels = c('All simulations', 'Surviving only')) +
	      guides(color = guide_legend('',title.position='top', override.aes = list(colour='black', linetype=c('solid','dashed')))) +
	      geom_smooth(data = fluDF, aes(colour=factor(sim), y=fluInc), size=0.5, linetype = 'dashed',  show_guide=FALSE) +
    	  theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    	  theme(axis.ticks = element_line(size=.1))  +
	      ylim(c(0,3)) +
    	  plot_themes 
	
	plot = plot_grid(plot1, plot2, labels = c('A','B'), align = 'h', ncol = 2)
	print(plot)
	save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 1, base_aspect_ratio = 0.85)
}

for(dir in resultsDirs){
  resultsDb = paste(dir,'results.sqlite',sep='')
  comboDb = dbConnect(SQLite(), dbname = resultsDb)
  initExtension(comboDb)
  if(dir == resultsDirs[1]){
    vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results WHERE tmrcaLimit==0')
    vaccineDF$sim = 'Dynamic'
  }
  else{
    df2 = dbGetQuery(comboDb, 'SELECT * FROM pooled_results WHERE tmrcaLimit == 0')
    df2$sim = 'Static'
  }
}

vaccineDF = rbind(vaccineDF, df2)

baseLineFlux = vaccineDF$meanFluxRate[vaccineDF$vaccinationRate==0 & vaccineDF$fluLike==1]
vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000

vaccineDF$acceleration = vaccineDF$meanFluxRate - mean(baseLineFlux)

vaccineDF$vaccinationRate = vaccineDF$vaccinationRate*365

summaryDF = ddply(vaccineDF, .(vaccinationRate, vaccineImmuneBreadth, sim), summarise,
                  meanDrift = mean(cumulativeDrift),
                  meanInc = mean(cumulativeIncidence))

fluDF = ddply(vaccineDF[vaccineDF$fluLike==1,], .(vaccinationRate, vaccineImmuneBreadth,sim), summarise,
              fluDrift = mean(cumulativeDrift),
              fluInc = mean(cumulativeIncidence))


makePlot(summaryDF, fluDF, 'stat_dyn')

cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeDrift, method = 'spearman')
cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeIncidence, method = 'spearman')
wilcox.test(vaccineDF$cumulativeDrift[vaccineDF$sim=='Dynamic'], vaccineDF$cumulativeDrift[vaccineDF$sim=='Static'], alternative = 'less')
wilcox.test(vaccineDF$cumulativeIncidence[vaccineDF$sim=='Dynamic'], vaccineDF$cumulativeIncidence[vaccineDF$sim=='Static'], alternative = 'less')