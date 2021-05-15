library(RSQLite)
library(ggplot2)
library(sensitivity)
library(gridExtra)
library(grid)
library(scales)
library(viridis)
library(plyr)
library(cowplot)
library(patchwork)

textSize = 12
pointSize = 1.0
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

single.plot = function(df, x, y, ylab, ylim, p.val){
  plot = ggplot(data = df, aes_string(x = x, y = y)) + 
    xlab('Vaccination coverage') +
    ylab(ylab) +
    geom_line(aes(group = vaccineImmuneBreadth, colour = vaccineImmuneBreadth), size=0.4) +
    geom_point(aes_string(alpha = p.val, colour = 'vaccineImmuneBreadth'), size=0.3, show.legend = FALSE) +
    geom_hline(yintercept = mean(df[,y][df$vaccinationRate==0]), alpha = 0.5, linetype = 'longdash', size=0.2) +
    scale_color_viridis(begin = 0, end = 0.9) +
    guides(colour= guide_colourbar('Vaccine immune breadth',title.position='bottom',barwidth=11,barheight=0.5)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    ylim(ylim) +
    plot_themes
  return(plot)
}

makePlot = function(summaryDF, fluDF, plotName){
  plot1 = single.plot(df = summaryDF, x='vaccinationRate', y='meanDrift', ylab = 'Cumulative antigenic evolution',ylim=c(0,25), p.val = 'drift.p') 
  plot2 = single.plot(df = summaryDF, x='vaccinationRate', y='meanInc', ylab = 'Cumulative incidence',ylim=c(0,2), p.val = 'inc.p') 
  plot3 = single.plot(df = fluDF, x='vaccinationRate', y='meanDrift', ylab = 'Cumulative antigenic evolution',ylim=c(17,29.5), p.val = 'drift.p') 
  plot4 = single.plot(df = fluDF, x='vaccinationRate', y='meanInc', ylab = 'Cumulative incidence',ylim=c(0.65,2.2), p.val = 'inc.p') 
  
  plot = (plot1+plot2)/(plot3+plot4) + plot_annotation(tag_levels = 'A')+ plot_layout(guides = "collect") & theme(legend.position = 'bottom')
    #plot_grid(plot1, plot2, plot3, plot4, labels = c('A','B','C','D'), align = 'h', ncol = 2)
  print(plot)
  save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 2, base_aspect_ratio = 0.85)
}

summarize.df = function(df){
  zeroDF = df[df$vaccinationRate==0,]
  out.df = ddply(df, .(vaccinationRate, vaccineImmuneBreadth), here(summarise),
                 meanDrift = mean(cumulativeDrift),
                 meanInc = mean(cumulativeIncidence),
                 count = sum(cumulativeDrift)/mean(cumulativeDrift),
                 drift.p.greater = wilcox.test(cumulativeDrift, zeroDF$cumulativeDrift, alternative='greater')$p.value,
                 drift.p.less = 1-drift.p.greater,
                 inc.p.greater = wilcox.test(cumulativeIncidence, zeroDF$cumulativeIncidence, alternative='greater')$p.value,
                 inc.p.less = 1-inc.p.greater,
                 drift.p = as.numeric(min(drift.p.greater, drift.p.less) < 0.05 & count >= 5),
                 inc.p = as.numeric(min(drift.p.greater, drift.p.less) < 0.05 & count >= 5)
  )
  return(out.df)
}

format.data = function(vaccineDF){
  vaccineDF = vaccineDF[vaccineDF$vaccineImmuneBreadth %in% c(1, 0.3,0.2, 0.1, 0.05,.5,.7),]
  
  vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000
  return(vaccineDF)
}

resultsDirs = c('../../analysis/breadth_1_density/vaccine/','../../analysis/breadth_low_density/vaccine/')

for(dir in resultsDirs){
  resultsDb = paste(dir,'results.sqlite',sep='')
  comboDb = dbConnect(SQLite(), dbname = resultsDb)
  initExtension(comboDb)
  if(dir == resultsDirs[1]){
    vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
  }
  else{
    df2 = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
  }
}

completeDF = rbind(vaccineDF, df2)
completeDF = format.data(completeDF)
vaccineDF = completeDF[completeDF$tmrcaLimit == 0,]

summaryDF = summarize.df(vaccineDF)
fluDF = summarize.df(vaccineDF[vaccineDF$fluLike==1,])

makePlot(summaryDF, fluDF, 'breadthSummary')

