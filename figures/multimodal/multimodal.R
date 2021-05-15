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

multimodal.plots = list()
i=1
for(breadth in sort(unique(vaccineDF$vaccineImmuneBreadth))){
  xdf = vaccineDF[vaccineDF$vaccineImmuneBreadth==breadth,]
  x=vaccineDF[vaccineDF$vaccineImmuneBreadth==breadth,]$cumulativeDrift
  plot1 = ggplot(vaccineDF[vaccineDF$vaccineImmuneBreadth==breadth,]) + 
    geom_freqpoly(aes(x=cumulativeDrift, color = vaccinationRate, group=vaccinationRate), binwidth=1.5) + 
    scale_color_viridis() +
    guides(color=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top', title = 'Annual vaccination coverage')) +
    plot_themes +
    ggtitle(paste("Breadth = ",breadth,sep='')) +
    theme(plot.title = element_text(hjust=0)) + xlab('Cumulative antigenic evolution') + ylim(0,400)
  
  multimodal.plots[[i]] = plot1
  i=i+1
}
plot = wrap_plots(multimodal.plots, ncol = 2) +
plotName = 'multimodal'
save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 4, base_aspect_ratio = 1)


rates = sort(unique(xdf$vaccinationRate))
for(i in 1:length(rates)){
  rate = rates[i]
  test=  wilcox.test(vaccineDF$cumulativeDrift[vaccineDF$cumulativeDrift>12 & vaccineDF$vaccinationRate==0], xdf$cumulativeDrift[xdf$cumulativeDrift>12 & xdf$vaccinationRate==rate])
  print(paste(rate,':',test$p.value))
}
